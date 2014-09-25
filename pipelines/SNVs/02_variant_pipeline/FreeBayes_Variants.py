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

class FreeBayes_Variants_Pipeline(Pipeline):
    """
    FreeBayes_Pipeline

    Arguments:
        base_dir (required) Path to a local directory in which temporary files,
            reference files and log files will be stored.

        n_threads (default=1) # of threads to use.

        remove_intermediate (default=False) Remove intermediate files if True.
    """
    def __init__(self, base_dir, n_threads, remove_intermediate, **kwargs):
        super(FreeBayes_Variants_Pipeline, self).__init__(base_dir=base_dir,
                                            n_threads=n_threads,
                                            remove_intermediate=remove_intermediate,
                                            **kwargs)
        self.name = "FreeBayes_Variants"
        
        self.files = {
            "reference_genome": self.as_ref("human_g1k_v37.fasta"),
            "reference_genome_dict": self.as_ref("human_g1k_v37.dict"),
            "interval_list_prefix": self.as_ref("nimblegen_solution_V2refseq2010.HG19"),
            "interval_list": self.as_ref("nimblegen_solution_V2refseq2010.HG19.bed"),
        }

        self.cmds = {
            "freebayes": "/data/bin/freebayes/bin/freebayes",
            "gatk": "java -d64 -Xmx10g -Djava.io.tmpdir=/mnt/temp -jar /usr/local/bin/gatk/GenomeAnalysisTK.jar",
            "s3cmd": "s3cmd --no-progress",
            "vcf-concat": "vcf-concat",
        }

        self.logs = {
            "stdout": self.create_logger("stdout", self.as_out("freebayes_variants.stdout.log")),
            "stderr": self.create_logger("stderr", self.as_out("freebayes_variants.stderr.log")),
            "pipeline": self.create_logger("pipeline", self.as_out("freebayes_variants.pipeline.log")),
            "sims": self.create_sims_logger(remote_host='master'),
        }
        self.FreebayesCallerParameters = {
            "--genotype-qualities": "",
            "--report-genotype-likelihood-max": "",
            "-P": 0.01,
        }


    ######## BEGIN PIPELINE METHODS

    def freebayes_Caller(self, bams_in, vcf_out, interval_list=None, sample_names=None):
        """
        STEP 1.
            Input: List of De-duped, realigned, recalibrated bams
            Output: VCF file

        Optionally takes a list of sample_names, which should be in the exact order of the bam files.
            The output file will use these sample names (instead of the sample names specified 
            in the bam file)
        """
        interval_list = interval_list or self.files["interval_list"] 
        self.cmd("{freebayes_cmd} \
                    --fasta-reference {reference_genome} \
                    --targets {interval_list} \
                    {parameters} \
                    {bam_list_in} \
                    > {vcf_out}"
            .format(
                freebayes_cmd=self.cmds["freebayes"],
                reference_genome=self.files["reference_genome"],
                interval_list=interval_list,
                parameters=" ".join(["%s %s" % k for k in self.FreebayesCallerParameters.iteritems()]),
                bam_list_in=" ".join(bams_in),
                vcf_out=vcf_out,
            ),
            shell=True)
        self.checkpoint(vcf_out)


    def concat_vcf(self, vcfs_in, vcf_out):
        """
        Concatenate multiple VCF files using vcf-concat from vcftools
        """

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
    parser.add_argument("batch_id", type=str, metavar="default_batch or batch file basename",
                                  help="Used to identify this batch of samples, should line up with batchID column in sample-table")
    args = parser.parse_args()

    
    s3_bucket = "s3://asdjre"
    familyID_list = []
    
    samples = pd.read_csv(args.sampletable)

    batchsamples = samples[samples.batchID == args.batch_id]

    familyID_list = list(set(batchsamples.familyID.values))

    print args

    n_threads = multiprocessing.cpu_count()

    p = FreeBayes_Variants_Pipeline(args.tempdir, n_threads=n_threads, remove_intermediate=False, checkpoint_path=args.checkpoint_path)

    p.name = "%s_FB_Variants" % args.batch_id
    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.batch_id
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

    # Ordered by the (descending) time it takes to run each chromosome
    # this maximizes (roughly) the utilization of the multiple cores over time.
    CONTIGS =["1", "2", "11", "6", "17", "3", "14", "4", "7", "19", "5", "12", "10", "8", "9", "20", "16", "15", "X", "18", "22", "13", "21", "Y"]
    print CONTIGS

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
    @files([[sample_names.values(), p.as_temp("%s.vcf" % contig), "%s.%s.20bp_slop.bed" % (p.files["interval_list_prefix"], contig)] for contig in CONTIGS])
    def run_freebayes_Caller(bams_in, vcf_out, interval_list):
        p.freebayes_Caller(bams_in, vcf_out, interval_list)
    
    @merge(run_freebayes_Caller, p.as_temp("%s.all_chr.freebayes.vcf" % args.batch_id))
    def run_concat_vcf(vcfs_in, vcf_out):
        p.concat_vcf(vcfs_in, vcf_out)

    @transform(run_concat_vcf, suffix(".all_chr.freebayes.vcf"), r"\1.all_chr.freebayes.sorted.vcf")
    def run_sort_vcf(vcf_in, vcf_out):
        p.sort_vcf(vcf_in, vcf_out)

    @files(run_sort_vcf, "dummy")
    @posttask(touch_file("split_vcf_done.flag"))
    def run_split_vcf_into_families(vcf_in, dummy):
        for fID in familyID_list:
            vcf_out = os.path.join(os.path.dirname(vcf_in), "%s.family.freebayes.sorted.vcf" % fID)
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
