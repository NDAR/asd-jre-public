import sys
from ruffus import *
import pandas as pd
from starpipe.Pipeline import Pipeline
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
            "reference_genome_contigs": self.as_ref("human_g1k_v37.contigs"),
            "interval_list": self.as_ref("nimblegen_solution_V2refseq2010.HG19.list"),
            "hapmap_3.3.b37.sites.vcf": self.as_ref("hapmap_3.3.b37.sites.vcf"),
            "1000G_omni2.5.b37.sites.vcf": self.as_ref("1000G_omni2.5.b37.sites.vcf"),
            "1000G_phase1.snps.high_confidence.vcf": self.as_ref("1000G_phase1.snps.high_confidence.vcf"),
            "mills_indels": self.as_ref("Mills_and_1000G_gold_standard.indels.b37.vcf"),
            "known_dbSNP": self.as_ref("dbsnp_137.b37.vcf"),
        }

        self.cmds = {
            "gatk": "java -d64 -Xmx10g -Djava.io.tmpdir=/mnt/temp -jar /usr/local/bin/gatk/GenomeAnalysisTK.jar",
            "s3cmd": "s3cmd --no-progress",            
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
            "-minPruning": 3,
        }
        self.RecalibrationParameters = {
            "SNP": "-percentBad 0.01 -minNumBad 1000 \
                   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap_3_sites}\
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
            "INDEL": "--maxGaussians 4 -percentBad 0.01 -minNumBad 1000 \
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

    def gatk_HaplotypeCaller(self, bams_in, vcf_out):
        """
        STEP 1.
            Input: List of De-duped, realigned, recalibrated bams (optionally ReducedReads bam)
            Output: VCF file
        """
        self.cmd("{gatk_cmd} -T HaplotypeCaller \
                    {bam_list_in} \
                    --dbsnp {dbsnp} \
                    -R {reference_genome} \
                    -L {interval_list} \
                    {calling_parameters} \
                    -nct {n_threads} \
                    -o {vcf_out}"
            .format(
                gatk_cmd=self.cmds["gatk"],
                bam_list_in=" ".join(["-I %s" % f for f in bams_in]),
                reference_genome=self.files["reference_genome"],
                dbsnp=self.files["known_dbSNP"],
                interval_list=self.files["interval_list"],
                calling_parameters=" ".join(["%s %s" % k for k in self.HaplotypeCallerParameters.iteritems()]),
                n_threads=n_threads,
                vcf_out=vcf_out,
            ),
            shell=True)
        self.checkpoint(vcf_out)

    def gatk_BuildVariantRecalibrationModel(self, vcf_in, out_files, mode):
        if mode not in ["SNP", "INDEL"]:
            raise PipelineException("Must specify either SNP or INDEL for mode"
                                    " in gatk_BuildVariantRecalibrationModel!")

        self.cmd("{gatk_cmd} -T VariantRecalibrator \
                    --input {vcf_in} \
                    -R {reference_genome} \
                    -recalFile {recal_out} \
                    -tranchesFile {tranches_out} \
                    -nt {n_threads} \
                    {recal_parameters}"
            .format(
                gatk_cmd=self.cmds["gatk"],
                vcf_in=vcf_in,
                reference_genome=self.files["reference_genome"],
                recal_out=out_files[0],
                tranches_out=out_files[1],
                n_threads=8,
                recal_parameters=self.RecalibrationParameters[mode]
            ),
            shell=True)
        self.checkpoint(out_files)

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry", action="store_true",
                        help="Print out task list only. Does not run commands.")
    parser.add_argument("--recent-tasks-only", "--recent", action="store_false", default=True,
                        help="Only run most recent tasks necessary. Use with caution!"
                             "Corresponds to the gnu_make_maximal_rebuild_mode mode in ruffus.")
    parser.add_argument("manifest", type=str,
                        help="path to NDAR manifest file")
    parser.add_argument("tempdir", type=str, 
                        help="path to temporary directory")
    parser.add_argument("--familyIDs", type=str, nargs="+", required=True,
                        help="SSC FamilyIDs to discover variants")
    parser.add_argument("--checkpoint-path", type=str, default=None,
                        help="path to checkpoint directory or S3 bucket")
    args = parser.parse_args()

    n_threads = multiprocessing.cpu_count()

    p = GATK_Variants_Pipeline(args.tempdir, n_threads=n_threads, remove_intermediate=False, checkpoint_path=args.checkpoint_path)

    p.name = "%s_GATK_Variants" % args.familyID
    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.familyID
    }

    manifest = ndar.NdarManifestReader(args.manifest)

    samples = manifest.get_family(args.familyID)

    input_files = []
    for s in samples:
        in_file = p.as_temp("%s.realigned.recal.reduced.bam" % s.sample_id)
        input_files.append(in_file)

    @check_if_uptodate(p.check_output_with_checkpoint)
    @merge(input_files, p.as_temp("%s.vcf" % args.familyID))
    def run_gatk_HaplotypeCaller(bams_in, vcf_out):
        p.gatk_HaplotypeCaller(bams_in, vcf_out)
    
    @check_if_uptodate(p.check_output_with_checkpoint)
    @transform(run_gatk_HaplotypeCaller, suffix(".vcf"), [r"\1.snps.vqsr.table", r"\1.snps.vqsr.tranches"])
    def run_gatk_BuildVQSRModel_snps(vcf_in, vqsr_files_out):
        p.gatk_BuildVariantRecalibrationModel(vcf_in, vqsr_files_out, "SNP")

    @check_if_uptodate(p.check_output_with_checkpoint)
    @transform(run_gatk_HaplotypeCaller, suffix(".vcf"), [r"\1.indel.vqsr.table", r"\1.indel.vqsr.tranches"])
    def run_gatk_BuildVQSRModel_indels(vcf_in, vqsr_files_out):
        p.gatk_BuildVariantRecalibrationModel(vcf_in, vqsr_files_out, "INDEL")
    
    @check_if_uptodate(p.check_output_with_checkpoint)    
    @transform([run_gatk_BuildVQSRModel_snps, run_gatk_BuildVQSRModel_indels], # use inputs from these tasks
               regex(r"(.*)\.(.*)vqsr.table"), # break up input here (group1 = familyID)
               add_inputs([r"\1.vcf"]), # need to add in the unfiltered VCF
               r"\1.snps.filtered.vcf") # this is the output
    def run_gatk_ApplyVQSR_snps(files_in, vcf_out):
        recal_file = files_in[0][0]
        ts_file = files_in[0][1]
        vcf_in = files_in[1][0]
        ts_filter_level = p.ts_filter_level["SNP"]
        p.gatk_ApplyVariantRecalibrationModel(vcf_in, ts_filter_level, ts_file,
            recal_file, vcf_out, "SNP")

    @check_if_uptodate(p.check_output_with_checkpoint)    
    @transform([run_gatk_BuildVQSRModel_snps, run_gatk_BuildVQSRModel_indels], # use inputs from these tasks
               regex(r"(.*)\.(.*)vqsr.table"), # break up input here (group1 = familyID)
               add_inputs([r"\1.vcf"]), # need to add in the unfiltered VCF
               r"\1.indels.filtered.vcf") # this is the output
    def run_gatk_ApplyVQSR_indels(files_in, vcf_out):
        recal_file = files_in[0][0]
        ts_file = files_in[0][1]
        vcf_in = files_in[1][0]
        ts_filter_level = p.ts_filter_level["INDEL"]
        p.gatk_ApplyVariantRecalibrationModel(vcf_in, ts_filter_level, ts_file,
            recal_file, vcf_out, "INDEL")

    @follows(run_gatk_ApplyVQSR_snps, run_gatk_ApplyVQSR_indels)
    def finalize():
        #p.cleanup_intermediate_files()
        p.log(action="%s.done" % p.name)

    if args.dry:
        pipeline_printout(sys.stdout, [finalize], verbose=9)
    else:
        pipeline_run([finalize], verbose=9, multiprocess=n_threads)
