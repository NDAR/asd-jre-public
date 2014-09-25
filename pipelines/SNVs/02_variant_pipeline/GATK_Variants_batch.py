import sys
import os
import subprocess
import tempfile
from ruffus import *
import pandas as pd
from starpipe.Pipeline import Pipeline, PipelineException
import argparse
import NdarManifestReader as ndar
import multiprocessing

class GATK_Variants_Pipeline(Pipeline):
    """
    GATK_Variants_Pipeline

    Arguments:
        base_dir (required) Path to a local directory in which temporary files,
            reference files and log files will be stored.

        n_threads (default=1) # of threads to use.

        remove_intermediate (default=False) Remove intermediate files if True.
    """
    def __init__(self, base_dir, n_threads, remove_intermediate, **kwargs):
        super(GATK_Variants_Pipeline, self).__init__(base_dir=base_dir,
                                            n_threads=n_threads,
                                            remove_intermediate=remove_intermediate,
                                            **kwargs)
        self.name = "GATK_Variants"
        self.default_stage=self.gatk_HaplotypeCaller

        self.files = {
            "reference_genome": self.as_ref("human_g1k_v37.fasta"),
            "reference_genome_dict": self.as_ref("human_g1k_v37.dict"),
            "interval_list": self.as_ref("nimblegen_solution_V2refseq2010.HG19.list"),
            "hapmap_3.3.b37.sites.vcf": self.as_ref("hapmap_3.3.b37.vcf"),
            "1000G_omni2.5.b37.sites.vcf": self.as_ref("1000G_omni2.5.b37.vcf"),
            "1000G_phase1.snps.high_confidence.vcf": self.as_ref("1000G_phase1.snps.high_confidence.b37.vcf"),
            "mills_indels": self.as_ref("Mills_and_1000G_gold_standard.indels.b37.vcf"),
            "known_dbSNP": self.as_ref("dbsnp_137.b37.vcf"),
        }

        self.cmds = {
            "gatk": "java -d64 -Xmx10g -Djava.io.tmpdir=/mnt/temp -jar /usr/local/bin/gatk/GenomeAnalysisTK.jar",
            "s3cmd": "s3cmd --no-progress",
            "vcf-concat": "vcf-concat",
        }

        self.logs = {
            "stdout": self.create_logger("stdout", self.as_out("gatk_variants.stdout.log")),
            "stderr": self.create_logger("stderr", self.as_out("gatk_variants.stderr.log")),
            "pipeline": self.create_logger("pipeline", self.as_out("gatk_variants.pipeline.log")),
            "sims": self.create_sims_logger(remote_host='master'),
        }
        self.HaplotypeCallerParameters = {
            "--standard_min_confidence_threshold_for_calling": 30.0,
            "--standard_min_confidence_threshold_for_emitting": 25.0,
        }
        self.num_bad = {
            "SNP": 3000,
            "INDEL": 2000
        }
        self.RecalibrationParameters = {
            "SNP": "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap_3_sites}\
                   -resource:omni,known=false,training=true,truth=true,prior=12.0 {Omni_1000G_sites} \
                   -resource:1000G,known=false,training=true,truth=false,prior=10.0 {Phase1_1000G_HC_SNPs} \
                   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {db_snp}\
                   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP \
                   -mode SNP"
            .format(
                hapmap_3_sites=self.files["hapmap_3.3.b37.sites.vcf"],
                Omni_1000G_sites=self.files["1000G_omni2.5.b37.sites.vcf"],
                Phase1_1000G_HC_SNPs=self.files["1000G_phase1.snps.high_confidence.vcf"],
                db_snp=self.files["known_dbSNP"]
            ),
            "INDEL": "--maxGaussians 4 \
                     -resource:mills,known=false,training=true,truth=true,prior=12.0 {mills_indels} \
                     -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {db_snp} \
                     -an DP -an FS -an ReadPosRankSum -an MQRankSum \
                     -mode INDEL"
            .format(
                mills_indels=self.files["mills_indels"],
                db_snp=self.files["known_dbSNP"],
            ),
        }

    ######## BEGIN PIPELINE METHODS

    def gatk_HaplotypeCaller(self, bams_in, vcf_out, interval_list=None, sample_names=None, nct_threads=2):
        """
        STEP 1.
            Input: List of De-duped, realigned, recalibrated bams (optionally ReducedReads bam)
            Output: VCF file

        Optionally takes a list of sample_names, which should be in the exact order of the bam files.
            The output file will use these sample names (instead of the sample names specified 
            in the bam file)
        """
        interval_list = interval_list or self.files["interval_list"] 
        if sample_names is not None:
            tf = tempfile.NamedTemporaryFile(delete=False)
            for bf, sn in zip(bams_in, sample_names):
                tf.write("%s %s\n" % (bf, sn))
            sample_name_file = "--sample_rename_mapping_file %s" % tf.name
        else:
            sample_name_file = ""
        tf.close()
        self.cmd("{gatk_cmd} -T HaplotypeCaller \
                    {bam_list_in} {sample_name_file} \
                    --dbsnp {dbsnp} \
                    -R {reference_genome} \
                    -L {interval_list} \
                    {calling_parameters} \
                    --max_alternate_alleles 10 \
                    -nct {nct_threads} \
                    -o {vcf_out}"
            .format(
                gatk_cmd=self.cmds["gatk"],
                bam_list_in=" ".join(["-I %s" % f for f in bams_in]),
                sample_name_file=sample_name_file,
                reference_genome=self.files["reference_genome"],
                dbsnp=self.files["known_dbSNP"],
                interval_list=interval_list,
                calling_parameters=" ".join(["%s %s" % k for k in self.HaplotypeCallerParameters.iteritems()]),
                nct_threads=nct_threads,
                vcf_out=vcf_out,
            ),
            shell=True,
            n_retry=2)
        self.checkpoint(vcf_out)

    def gatk_BuildVariantRecalibrationModel(self, vcf_in, out_files, mode):
        if mode not in ["SNP", "INDEL"]:
            raise PipelineException("Must specify either SNP or INDEL for mode"
                                    " in gatk_BuildVariantRecalibrationModel!")

        return_code = self.cmd("{gatk_cmd} -T VariantRecalibrator \
                    --input {vcf_in} \
                    -R {reference_genome} \
                    -recalFile {recal_out} \
                    -tranchesFile {tranches_out} \
                    -nt {n_threads} \
                    --rscript_file {rscript_file} \
                    -numBad {num_bad} \
                    {recal_parameters}"
            .format(
                gatk_cmd=self.cmds["gatk"],
                vcf_in=vcf_in,
                reference_genome=self.files["reference_genome"],
                recal_out=out_files[0],
                tranches_out=out_files[1],
                n_threads=8,
                rscript_file=out_files[3],
                num_bad=self.num_bad[mode],
                recal_parameters=self.RecalibrationParameters[mode]
            ),
            shell=True,
            stop_on_error=False)
        self.checkpoint(out_files)
        return return_code

    def gatk_ApplyVariantRecalibrationModel(self, vcf_in, ts_filter_level, ts_file, recal_file, vcf_out, mode):
        if mode not in ["SNP", "INDEL", "BOTH"]:
            raise PipelineException("Must specify either SNP or INDEL for mode"
                                    " in gatk_ApplyVariantRecalibrationModel!")

        self.cmd("{gatk_cmd} -T ApplyRecalibration \
                    -R {reference_genome} \
                    -input {vcf_in} \
                    --ts_filter_level {ts_filter_level} \
                    -tranchesFile {ts_file} \
                    -recalFile {recal_file} \
                    -mode {mode} \
                    -o {vcf_out}"
                    .format(
                        gatk_cmd=self.cmds["gatk"],
                        reference_genome=self.files["reference_genome"],
                        vcf_in=vcf_in,
                        ts_filter_level=ts_filter_level,
                        ts_file=ts_file,
                        recal_file=recal_file,
                        mode=mode,
                        vcf_out=vcf_out,
                        ),
                    shell=True)
        self.checkpoint(vcf_out)

    def concat_vcf(self, vcfs_in, vcf_out):
        """
        Concatenate multiple VCF files using vcf-concat from vcftools
        """
        # self.cmd("{vcf_concat} {in_vcfs} > {out_vcf}"
        #     .format(
        #         vcf_concat=self.cmds["vcf-concat"],
        #         in_vcfs=" ".join(vcfs_in),
        #         out_vcf=vcf_out
        #         ),
        #     shell=True)
        self.cmd("{gatk_cmd} -T CombineVariants \
               -R {reference_genome} \
               {input_vcfs} \
               -o {vcf_out} \
               -genotypeMergeOptions UNSORTED"
               .format(
                    gatk_cmd=self.cmds["gatk"],
                    reference_genome=self.files["reference_genome"],
                    input_vcfs=" ".join(["-V %s" % f for f in vcfs_in]),
                    vcf_out=vcf_out
                ),
               shell=True)

        self.checkpoint(vcf_out)


    def sort_vcf(self, vcf_in, vcf_out):
        """
        Sort a VCF in correct chromosome order. Preserves header
        """
        self.cmd("(cat {in_vcf} | head -5000 | grep ^#; \
                   grep -v ^# {in_vcf} | sort -k1,1V -k2,2n) \
                   > {out_vcf}"
                  .format(
                    in_vcf=vcf_in,
                    out_vcf=vcf_out 
                    ),
                shell=True)

    def split_vcf_into_families(self, vcf_in, fID, vcf_out):
        """
        Split a large VCF into a family vcf based on familyID 
        """
        self.cmd("{gatk_cmd} -T SelectVariants \
                -R {reference_genome} \
                --variant {vcf_in} \
                -o {vcf_out} \
                -se \"{fID}.+\"\
                --excludeNonVariants"
            .format(
                gatk_cmd=self.cmds["gatk"],
                reference_genome=self.files["reference_genome"],
                vcf_in=vcf_in,
                fID=fID,
                vcf_out=vcf_out
            ),
            shell=True)

    def download_bam_file_from_s3(self, remote_bam, local_bam):
        """
        Download input bam files to self.local_temp_dir
        """
        self.cmd("{s3cmd} get {remote} {local}"
            .format(
                s3cmd=self.cmds["s3cmd"],
                remote=remote_bam,
                local=local_bam,
            ),
            shell=True)

    def check_s3_download_finished(self, remote_bam, local_bam):
        """
        Check if downloaded bam size is the same as remote bam size
        TODO: lets move this all to using boto...
        """
        remote_size = subprocess.check_output("{s3cmd} ls {remote} | awk '{{ if ($4 == \"{remote}\") print $3}}'"
            .format(
                s3cmd=self.cmds["s3cmd"],
                remote=remote_bam,
            ),
            shell=True)
        try:
            local_size = os.path.getsize(local_bam)
        except OSError:
            print "File %s does not exist" % local_bam
            return False
        try:
            if int(remote_size) == local_size:
                return True
            else:
                if remote_bam[-3:] != local_bam[-3:]: # check if fileendings are the same (see if file was decompressed)
                    print "File endings of local and remote paths are different! Skipping check for %s" % remote_bam
                    return True
                return False
        except ValueError:
            print "Could not get remote size!", remote_bam, remote_size
            return False

    def compare_s3_filesize(self, remote_bam, local_bam):
        """
        Check if downloaded bam size is the same as remote bam size
        TODO: lets move this all to using boto...
        """
        remote_size = subprocess.check_output("{s3cmd} ls {remote} | awk '{{ if ($4 == \"{remote}\") print $3}}'"
            .format(
                s3cmd=self.cmds["s3cmd"],
                remote=remote_bam,
            ),
            shell=True)
        try:
            local_size = os.path.getsize(local_bam)
        except OSError:
            local_size = None
        try:
            if int(remote_size) == local_size:
                return True
            else:
                return False
        except ValueError:
            print "Could not get remote size!", remote_bam, remote_size
            return False

    def upload_to_s3(self, local_file, remote_path):
        """
        Uploads files to S3 bucket using s3cmd
        """
        self.cmd("{s3cmd} put {local_file} {remote_path}"
            .format(
                s3cmd=self.cmds["s3cmd"],
                local_file=local_file,
                remote_path=remote_path,
            ),
            log_command=False,
            shell=True,
            on_error=lambda: os.rename(local_file, "/data/failed_uploads/%s" % os.path.basename(local_file))
        )

        if not self.compare_s3_filesize(remote_path, local_file):
            raise PipelineException(
                "S3 upload of local file %s failed to destination %s" % (local_file, remote_path))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry", action="store_true",
                        help="Print out task list only. Does not run commands.")
    parser.add_argument("--recent-tasks-only", "--recent", action="store_false", default=True,
                        help="Only run most recent tasks necesseary. Use with caution!"
                             "Corresponds to the gnu_make_maximal_rebuild_mode mode in ruffus.")    
    parser.add_argument("sampletable", type=str,
                        help="path to sample file. Must have 'path', 'sampleID', 'familyID', and 'batchID'")
    parser.add_argument("tempdir", type=str, 
                        help="path to temporary directory")
    parser.add_argument("--checkpoint-path", type=str, default=None,
                         help="path to checkpoint directory or S3 bucket")
    parser.add_argument("--no-vqsr", action="store_true", default=False,
                         help="Do not run VQSR annotation/filtering on variants")
    parser.add_argument("batch_id", type=str, metavar="default_batch or batch file basename",
                                  help="Used to identify this batch of samples, should line up with batchID column in sample-table")
    parser.add_argument("--nct-threads", type=int, default=2)
    args = parser.parse_args()

    #manifest = ndar.NdarManifestReader(args.manifest)
    s3_bucket = "s3://asdjre"
    familyID_list = []

    samples = pd.read_csv(args.sampletable)

    batchsamples = samples[samples.batchID == args.batch_id]

    familyID_list = list(set(batchsamples.familyID.values))
    
    print args
    print familyID_list

    n_threads = multiprocessing.cpu_count()

    p = GATK_Variants_Pipeline(args.tempdir, n_threads=n_threads, remove_intermediate=False, checkpoint_path=args.checkpoint_path)

    p.name = "%s_GATK_Variants" % args.batch_id
    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.batch_id
    }
    p.ts_filter_level = {
        "SNP": 99.0,
        "INDEL": 99.0
    }
    
    batchsamples["local_path"] = map(lambda x: p.as_temp(os.path.basename(x[1]["path"])), batchsamples.iterrows())
    input_files = dict(batchsamples.set_index("path")["local_path"])
    
    bai_files = {}
    for k,v in input_files.iteritems():
        k = k[0:-3] + "bai"
        v = v[0:-3] + "bai"
        bai_files[k] = v
    input_files.update(bai_files)

    sample_names = dict(batchsamples.set_index("sampleID")["local_path"])

    for s,f in sample_names.iteritems():
        print s, f

    files_to_upload = {}

    CONTIGS = ["part%0.3d" % i for i in range(50)]
    
    def do_check_S3_download_finished(remote, local):
        # note this reverses True and False due to ruffus implementation
        if p.check_s3_download_finished(remote, local) == False:
            return True, "S3 download incomplete"
        else:
            return False, "S3 file already fully downloaded"

    #@follows(run_download_reference_file_from_s3)
    @check_if_uptodate(do_check_S3_download_finished)
    @parallel(zip(input_files.keys(), input_files.values()))
    def run_download_bam_file_from_s3(s3_file, local_file):
        p.download_bam_file_from_s3(s3_file, local_file)

    @jobs_limit(n_threads-1)
    @follows(run_download_bam_file_from_s3)
    @check_if_uptodate(p.check_output_with_checkpoint)
    @files([[sample_names.values(), p.as_temp("%s.vcf" % contig), "%s.%s.20bp_slop.list" % (p.files["interval_list"], contig)] for contig in CONTIGS])
    def run_gatk_HaplotypeCaller(bams_in, vcf_out, interval_list):
        p.gatk_HaplotypeCaller(bams_in, vcf_out, interval_list, sample_names.keys(), nct_threads=args.nct_threads)
    
    @merge(run_gatk_HaplotypeCaller, p.as_temp("%s.all_chr.vcf" % args.batch_id))
    def run_concat_vcf(vcfs_in, vcf_out):
        p.concat_vcf(vcfs_in, vcf_out)

    @transform(run_concat_vcf, suffix(".all_chr.vcf"), r"\1.all_chr.sorted.vcf")
    def run_sort_vcf(vcf_in, vcf_out):
        p.sort_vcf(vcf_in, vcf_out)

    if not args.no_vqsr:
        @check_if_uptodate(p.check_output_with_checkpoint)
        @transform(run_sort_vcf, suffix(".all_chr.sorted.vcf"), [r"\1.snps.vqsr.table", r"\1.snps.vqsr.tranches", r"\1.snps.vqsr.tranches.pdf", r"\1.snps.rscript", r"\1.snps.rscript.pdf"])
        def run_gatk_BuildVQSRModel_snps(vcf_in, vqsr_files_out):
            success = False
            n_try = 0
            while n_try <= 3:
                n_try += 1
                success = p.gatk_BuildVariantRecalibrationModel(vcf_in, vqsr_files_out, "SNP")
                if success:
                    break
                else:
                    p.num_bad["SNP"] = p.num_bad["SNP"] + 1000
            if not success:
                p.log(action="%s.error" % p.name, message="Could not find SNP VQSR parameters!")
                sys.exit(1)
            else:
                for f in vqsr_files_out:
                    p.upload_to_s3(f, "%s/GATK_HC_batches/%s" % (s3_bucket, os.path.basename(f)))

        @check_if_uptodate(p.check_output_with_checkpoint)
        @transform(run_sort_vcf, suffix(".all_chr.sorted.vcf"), [r"\1.indel.vqsr.table", r"\1.indel.vqsr.tranches", r"\1.indel.vqsr.tranches.pdf", r"\1.indel.rscript", r"\1.indel.rscript.pdf"])
        def run_gatk_BuildVQSRModel_indels(vcf_in, vqsr_files_out):
            success = False
            n_try = 0
            while n_try <= 3:
                n_try += 1
                success = p.gatk_BuildVariantRecalibrationModel(vcf_in, vqsr_files_out, "INDEL")
                if success:
                    break
                else:
                    p.num_bad["INDEL"] = p.num_bad["INDEL"] + 1000
            if not success:
                p.log(action="%s.error" % p.name, message="Could not find INDEL VQSR parameters!")
                sys.exit(1)
            else:
                for f in vqsr_files_out:
                    p.upload_to_s3(f, "%s/GATK_HC_batches/%s" % (s3_bucket, os.path.basename(f)))

        @check_if_uptodate(p.check_output_with_checkpoint)
        @transform(run_gatk_BuildVQSRModel_snps, # use inputs from these tasks
                   suffix(".snps.vqsr.table"), # break up input here (group1 = batch_id)
                   add_inputs([r"\1.all_chr.sorted.vcf"]), # need to add in the unfiltered VCF
                   r"\1.snps.filtered.vcf") # this is the output
        def run_gatk_ApplyVQSR_snps(files_in, vcf_out):
            recal_file = files_in[0][0]
            ts_file = files_in[0][1]
            vcf_in = files_in[1][0]
            ts_filter_level = p.ts_filter_level["SNP"]
            p.gatk_ApplyVariantRecalibrationModel(vcf_in, ts_filter_level, ts_file,
                recal_file, vcf_out, "SNP")

        @check_if_uptodate(p.check_output_with_checkpoint)
        @transform(run_gatk_BuildVQSRModel_indels, # use inputs from these tasks
                   suffix(".indel.vqsr.table"), # break up input here (group1 = batch_id)
                   add_inputs([r"\1.all_chr.sorted.vcf"]), # need to add in the unfiltered VCF
                   r"\1.indels.filtered.vcf") # this is the output
        def run_gatk_ApplyVQSR_indels(files_in, vcf_out):
            recal_file = files_in[0][0]
            ts_file = files_in[0][1]
            vcf_in = files_in[1][0]
            ts_filter_level = p.ts_filter_level["INDEL"]
            p.gatk_ApplyVariantRecalibrationModel(vcf_in, ts_filter_level, ts_file,
                recal_file, vcf_out, "INDEL")

        @merge([run_gatk_ApplyVQSR_indels, run_gatk_ApplyVQSR_snps], p.as_temp("%s.all_chr.vqsr.vcf" % args.batch_id))
        def run_concat_vqsr(vcfs_in, vcf_out):
            p.concat_vcf(vcfs_in, vcf_out)

        @transform(run_concat_vqsr, suffix(".all_chr.vqsr.vcf"), r"\1.all_chr.vqsr.sorted.vcf")
        def run_sort_final_vcf(vcf_in, vcf_out):
            p.sort_vcf(vcf_in, vcf_out)
            p.upload_to_s3(vcf_out, "%s/GATK_HC_batches/%s" % (s3_bucket, os.path.basename(vcf_out)))

        #@split(run_sort_final_vcf, regex(r"(.+).all_chr.vqsr.sorted.vcf"), r"\1.*.family.sorted.vcf")
        @files(run_sort_final_vcf, "dummy")
        @posttask(touch_file("split_vcf_done.flag"))
        def run_split_vcf_into_families(vcf_in, dummy):
            for fID in familyID_list:
                vcf_out = os.path.join(os.path.dirname(vcf_in), "%s.family.vqsr.sorted.vcf" % fID)
                print vcf_out
                p.split_vcf_into_families(vcf_in, fID, vcf_out)
                p.upload_to_s3(vcf_out, "%s/complete/%s/variants/%s" % (s3_bucket, fID, os.path.basename(vcf_out)))
    
    else:
        # do not run VQSR
        @files(run_sort_vcf, "dummy")
        @posttask(touch_file("split_vcf_done.flag"))
        def run_split_vcf_into_families(vcf_in, dummy):
            p.upload_to_s3(vcf_in, "%s/GATK_HC_batches/%s" % (s3_bucket, os.path.basename(vcf_in)))
            for fID in familyID_list:
                vcf_out = os.path.join(os.path.dirname(vcf_in), "%s.family.sorted.vcf" % fID)
                print vcf_out
                p.split_vcf_into_families(vcf_in, fID, vcf_out)
                p.upload_to_s3(vcf_out, "%s/complete/%s/variants/%s" % (s3_bucket, fID, os.path.basename(vcf_out)))
        

    @follows(run_split_vcf_into_families)
    def finalize():
        #p.cleanup_intermediate_files()
        p.log(action="%s.done" % p.name)

    if args.dry:
        pipeline_printout(sys.stdout, [finalize], gnu_make_maximal_rebuild_mode=args.recent_tasks_only, verbose=9)
    else:
        pipeline_run([finalize], verbose=9, gnu_make_maximal_rebuild_mode=args.recent_tasks_only, multiprocess=n_threads)
