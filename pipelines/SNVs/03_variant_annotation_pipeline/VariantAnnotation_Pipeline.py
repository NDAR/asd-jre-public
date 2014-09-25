import sys
import os
import glob
import subprocess
from ruffus import *
from starpipe.Pipeline import Pipeline
import argparse
import multiprocessing

import NdarManifestReader as ndar


class VariantAnnotation_PipelineException(Exception):
    pass

class VariantAnnotation_Pipeline(Pipeline):
    """
    VariantAnnotation Pipeline

    Arguments:
        base_dir (required) Path to a local directory in which temporary files,
            reference files and log files will be stored.

        n_threads (default=1) # of threads to use.
    """
    def __init__(self, base_dir, n_threads, remove_intermediate, **kwargs):
        super(VariantAnnotation_Pipeline, self).__init__(base_dir=base_dir,
                                            n_threads=n_threads,
                                            remove_intermediate=False,
                                            **kwargs)

        self.cmds = {
            "rsync": "rsync --recursive --perms --group --verbose --compress --update",
            "s3cmd": "s3cmd -c /root/.s3cfg --no-progress --force",
            "snpeff": "java -Xmx4g -jar /data/bin/snpEff/snpEff.jar",
            "snpsift": "java -Xmx6g -jar /data/bin/snpEff/SnpSift.jar",
            "famseq": "/data/bin/FamSeq/FamSeq",
            "vcf-annotate": "/usr/bin/vcftools/vcftools_0.1.11/bin/vcf-annotate",
            "vcftools_lib": "export PERL5LIB=/usr/bin/vcftools/vcftools_0.1.11/perl"
        }

        self.logs = {
            "stdout": None,#self.create_logger("stdout", self.as_out("init.stdout.log")),
            "stderr": None, #self.create_logger("stderr", self.as_out("init.stderr.log")),
            #"pipeline": self.create_logger("pipeline", self.as_out("bwa.pipeline.log")),
            #"sims": self.create_sims_logger(remote_host='master'),
        }
        self.files = {
            "trf_track": self.as_ref("simple_repeats.sorted.header.bed.gz"),
            "segdup_track": self.as_ref("hg19genomicSuperDups.sorted.merged.header.bed.gz"),
            "esp_track": self.as_ref("ESP6500.coverage.all_chrs.sorted.txt.gz"),
            "dbsnp_track": self.as_ref("dbsnp_137.b37.vcf"),
            "dbnsfp_track": self.as_ref("dbNSFP2.1.txt"),
            "cadd_track": self.as_ref("cadd_v1.0.nimblegen_solution_V2refseq2010.HG19.300bp_slop.vcf"),
        }
        self.params = {
            "tstv": 2.0,
        }
        self.snpeff_parameters = "-lof -v hg19 -no-downstream -no-intergenic -no-upstream"
        self.dbnsfp_fields =  ("Uniprot_id,Uniprot_aapos,Interpro_domain,Ancestral_allele,SIFT_score_converted,Polyphen2_HDIV_score,"
                               "Ensembl_geneid,Ensembl_transcriptid,"
                              "Polyphen2_HVAR_score,LRT_score_converted,MutationTaster_score_converted,MutationAssessor_score_converted,"
                              "FATHMM_score_converted,GERP++_NR,GERP++_RS,phyloP,29way_pi,29way_logOdds,LRT_Omega,"
                              "1000Gp1_AC,1000Gp1_AFR_AC,1000Gp1_EUR_AC,1000Gp1_AMR_AC,1000Gp1_ASN_AC,ESP6500_AA_AF,ESP6500_EA_AF,")

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
            raise VariantAnnotation_PipelineException(
                "S3 upload of local file %s failed to destination %s" % (local_file, remote_path))

    def download_from_s3(self, remote_bam, local_bam):
        self.cmd("{s3cmd} get {remote} {local}"
            .format(
                s3cmd=self.cmds["s3cmd"],
                remote=remote_bam,
                local=local_bam,
            ),
            shell=True)

    def cleanup(self):
        """
        Remove all files in the local_temp_dir and local_out_dir
        """
        self.cmd("rm -f {local_temp_dir}/* \
                        {local_out_dir}/*".
            format(
                local_temp_dir=self.local_temp_dir,
                local_out_dir=self.local_out_dir
            ),
            shell=True)
        print "Cleaning up!"

    def family_postprocessor(self, in_vcf, in_ped, out_vcf, postprocessor=None):
        if postprocessor == "polymutt":
            self.cmd("{polymutt} -p {in_ped} \
                        -d <(touch blank.dat) \
                        --nthreads=4 \
                        --poly_tstv={tstv} \
                        --all_sites --chrX --chrY \
                        --in_vcf {in_vcf} \
                        --out_vcf {out_vcf}"
                .format(
                    polymutt=self.cmds["polymutt"],
                    tstv=self.params["tstv"],
                    in_vcf = in_vcf,
                    out_vcf = out_vcf
                    ),
                shell=True
            )
        elif postprocessor == "famseq":
            self.cmd("{famseq} vcf -a \
                        -vcfFile {in_vcf} \
                        -pedFile {in_ped} \
                        -output {out_vcf}"
                .format(
                    famseq=self.cmds["famseq"],
                    in_vcf=in_vcf,
                    in_ped=in_ped,
                    out_vcf=out_vcf,
                    ),
                shell=True)
            # this will clean up the improper formatted VCF file from FamSeq
            self.cmd("sed 's/^##FSP /##FSP=/g' -i %s" % out_vcf,
                shell=True)

        elif postprocessor == None:
            self.cmd("cp {in_vcf} {out_vcf}"
                .format(
                    in_vcf=in_vcf,
                    out_vcf=out_vcf
                ), shell=True
            )
        else:
            raise VariantAnnotation_PipelineException("Unknown postprocessor (%s) selected!" % self.family_postprocessor)

    def ensemble_caller(self, in_vcfs, out_vcf):
        #self.cmd("")
        pass

    def mendelian_annotation(self, in_vcf, out_vcf):
        self.cmd("python mendelian_annotator.py \
                    {in_vcf} {out_vcf}"
                .format(
                    in_vcf=in_vcf,
                    out_vcf=out_vcf
                    ),
                shell=True)
    
    def snpeff_annotation(self, in_vcf, out_vcf, annotation="snpeff", header=None):
        """Use snpEff to annotate the variants in the VCF"""
        
        if annotation == "snpeff":
            self.cmd("{snpeff} eff \
                     {snpeff_parameters} \
                     {in_vcf} > {out_vcf}"
                .format(
                    snpeff=self.cmds["snpeff"],
                    snpeff_parameters=self.snpeff_parameters,
                    in_vcf=in_vcf,
                    out_vcf=out_vcf,
                ), shell=True
            )
        elif annotation == "dbnsfp":
            self.cmd("{snpsift} dbnsfp \
                      -f {dbnsfp_fields} \
                      -v {dbnsfp_file} \
                      {in_vcf} > {out_vcf}"
                    .format(
                        snpsift = self.cmds["snpsift"],
                        dbnsfp_file=self.files["dbnsfp_track"],
                        dbnsfp_fields=self.dbnsfp_fields,
                        in_vcf=in_vcf,
                        out_vcf=out_vcf,
                    ), shell=True)
        else:
            self.cmd("{snpsift} annotate \
                      -v {annotation_vcf} \
                      {in_vcf} > {out_vcf}"
                    .format(
                        snpsift = self.cmds["snpsift"],
                        annotation_vcf=annotation,
                        in_vcf=in_vcf,
                        out_vcf=out_vcf,
                    ), shell=True)

        if header:
            header_str = "\n".join(header)
            self.cmd(r"sed -i 's/^#CHROM/{header}\n#CHROM/' {out_vcf}"
                .format(
                    header=header_str.encode('string_escape'),
                    out_vcf=out_vcf
                ),
                shell=True)

    def vcf_annotation(self, in_vcf, out_vcf, annotation, 
            header_id, header_number, header_type, header_description, header_columns=None, header_key="INFO"):
        """use vcf-annotate to annotate a vcf with another vcf and header"""
        
        header="key=%s,ID=%s,Number=%s,Type=%s,Description='%s'" % (header_key, header_id, header_number, header_type, header_description)
        
        if not header_columns:
            columns="CHROM,FROM,TO,%s/%s" % (header_key, header_id)
        else:
            columns=header_columns

        self.cmd("{vcftools_lib} && \
                cat {in_vcf} | {vcf_annotate} \
                -a {annotation_vcf} \
                -c {columns} \
                -d {header} > {out_vcf}"
            .format(
                in_vcf = in_vcf,
                vcftools_lib = self.cmds["vcftools_lib"],
                vcf_annotate = self.cmds["vcf-annotate"],
                columns=columns,
                header=header,
                annotation_vcf=annotation,
                out_vcf=out_vcf),
            shell=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--dry", action="store_true",
                        help="Print out task list only. Does not run commands.")
    parser.add_argument("--recent-tasks-only", "--recent", action="store_false", default=True,
                        help="Only run most recent tasks necessary. Use with caution!"
                             "Corresponds to the gnu_make_maximal_rebuild_mode mode in ruffus.")        
    parser.add_argument("manifest", type=str,
                        help="path to Ndar manifest file")
    parser.add_argument("familyID", type=str,
                        help="SSC FamilyID to process")
    parser.add_argument("tempdir", type=str,
                        help="path to temporary directory")
    parser.add_argument("vcf", type=str,
                        help="S3 location of VCF file to annotate")
    parser.add_argument("--enable-ensemble-calls", default=False)
    parser.add_argument("--disable-famseq", default=False, action="store_true")

    args = parser.parse_args()
    n_threads = multiprocessing.cpu_count()
    p = VariantAnnotation_Pipeline(args.tempdir, n_threads=n_threads, remove_intermediate=False)
    
    p.name = "%s_VariantAnnotation" % args.familyID

    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.familyID
    }
    
    # manifest = ndar.NdarManifestReader(args.manifest)
    # samples = manifest.get_family(args.familyID)
    
    #s3_bucket_path = "s3://asdjre/complete/%s" % args.familyID

    # vcf_files = {#"s3://asdjre/complete/{familyID}/variants/{familyID}.family.vqsr.sorted.vcf".format(familyID=args.familyID):
    #              #   p.as_temp("%s.family.vqsr.sorted.vcf" % args.familyID)}
    #             "s3://asdjre/complete/{familyID}/variants/{familyID}.family.freebayes.sorted.vcf".format(familyID=args.familyID):
    #                 p.as_temp("%s.family.freebayes.sorted.vcf" % args.familyID)} # get both freebayes and GATK if possible} # get both freebayes and GATK if possible

    vcf_files = {
        args.vcf: p.as_temp("%d.input.vcf" % args.familyID)
    }
    remote_ped_file = "s3://asdjre/complete/{familyID}/info/{familyID}.ped".format(familyID=args.familyID)

    def do_check_S3_upload_finished(local, remote):
        # note this reverses True and False due to ruffus implementation
        if p.compare_s3_filesize(remote, local) == False:
            return True, "S3 upload incomplete"
        else:
            return False, "S3 file fully uploaded"

    def start():
        p.log(action="%s.start" % p.name)
        for remote_vcf, local_vcf in vcf_files.iteritems():    
            p.download_from_s3(remote_vcf, local_vcf)
        p.download_from_s3(remote_ped_file, p.as_temp("%s.ped" % args.familyID))

    @follows(start)
    @transform(vcf_files.values(), suffix(".input.vcf"), add_inputs(p.as_temp("%s.ped" % args.familyID)), r"\1.clean.vcf" )
    def run_family_postprocessor(inputs, out_vcf):
        vcf_file = inputs[0]
        ped_file = inputs[1]
        if args.disable_famseq:
            postprocessor=None
        else:
            postprocessor="famseq"
        p.family_postprocessor(vcf_file, ped_file, out_vcf, postprocessor=postprocessor)

    if args.enable_ensemble_calls:
        @merge(run_family_postprocessor, p.as_temp("%s.family.clean.ensemble.vcf" % args.familyID))
        def run_ensemble(vcf_files, vcf_out):
            p.ensemble_caller(vcf_files, vcf_out)
    else:
        # passthrough
        # @transform(run_family_postprocessor, suffix("*.vcf") p.as_temp("%s.family.clean.vcf" % args.familyID))
        # def run_ensemble(input_file, output):
        #     for f in input_file:
        #         p.cmd("cp {input} {output}"
        #             .format(
        #                 input=f,
        #                 output=output,
        #                 ), shell=True
        #             )
        run_ensemble = run_family_postprocessor

    # @transform(run_ensemble, suffix(".vcf"), r"\1.ann.vcf")
    # def run_mendelian_annotation(vcf_in, vcf_out):
    #     p.mendelian_annotation(vcf_in, vcf_out)

    @transform(run_ensemble, suffix(".vcf"), r"\1.eff.vcf")
    def run_snpeff_annotation(vcf_in, vcf_out):
        p.snpeff_annotation(vcf_in, vcf_out, annotation="snpeff")
    
    @transform(run_snpeff_annotation, suffix(".vcf"), r"\1.dbnsfp.vcf")
    def run_dbnsfp_annotation(vcf_in, vcf_out):
        p.snpeff_annotation(vcf_in, vcf_out, annotation="dbnsfp")

    # @transform(run_dbnsfp_annotation, suffix(".vcf"), r"\1.dbsnp.vcf")
    # def run_dbnsnp_annotation(vcf_in, vcf_out):
    #     p.snpeff_annotation(vcf_in, vcf_out, annotation=self.files["dbsnp_track"])

    @transform(run_dbnsfp_annotation, suffix(".vcf"), r"\1.trf.vcf")
    def run_trf_annotation(vcf_in, vcf_out):
        #key=INFO,ID=TRF_score,Number=.,Type=Integer,Description='Tandem Repeat Finder score[s], if annotated. Higher scores are more repetitive. Multiple scores can be reported'
        params = {"annotation": p.files["trf_track"],
                  "header_id": "TRF_score",
                  "header_number": ".",
                  "header_type": "Integer",
                  "header_description": "Tandem Repeat Finder score[s], if annotated. Higher scores are more repetitive. Multiple scores can be reported",
                  "header_columns": "CHROM,FROM,TO,-,INFO/TRF_score"}
        p.vcf_annotation(vcf_in, vcf_out, **params)
    
    @transform(run_trf_annotation, suffix(".vcf"), r"\1.sd.vcf")
    def run_sd_annotation(vcf_in, vcf_out):
        #key=INFO,ID=SegDup,Number=0,Type=Flag,Description='Variant in a segmental duplication' 
        params = {"annotation": p.files["segdup_track"],
                  "header_id": "SegDup",
                  "header_number": "0",
                  "header_type": "Flag",
                  "header_description": "Variant in a segmental duplication.",
                  "header_columns": "CHROM,FROM,TO,INFO/SegDup"}
        p.vcf_annotation(vcf_in, vcf_out, **params)
    
    @transform(run_sd_annotation, suffix(".vcf"), r"\1.cadd.vcf")
    def run_cadd_annotation(vcf_in, vcf_out):
        header = [r'##INFO=<ID=CADD_p,Number=1,Type=Float,Description="CADD score PHRED-scaled (Kircher et al.)">',
                  r'##INFO=<ID=CADD_score,Number=1,Type=Float,Description="CADD raw score; v1.0; (Kircher et al.)">']
        p.snpeff_annotation(vcf_in, vcf_out, annotation=p.files["cadd_track"], header=header)
    
        
    @files(run_cadd_annotation, "dummy")
    @posttask(touch_file("%d.annotation_done.flag"))
    def run_upload_files(vcf_in, dummy):
        if "vqsr" in vcf_out:
            p.upload_to_s3(vcf_out, "s3://asdjre/complete/{familyID}/variants/{familyID}.family.vqsr.famseq.annotated.vcf"
                .format(familyID=args.familyID))
        elif "freebayes" in vcf_out:
            p.upload_to_s3(vcf_out, "s3://asdjre/complete/{familyID}/variants/{familyID}.family.freebayes.famseq.annotated.vcf"
                .format(familyID=args.familyID))
        else:
            p.upload_to_s3(vcf_out, "s3://asdjre/complete/{familyID}/variants/{familyID}.family.annotated.vcf"
                .format(familyID=args.familyID))

    @follows(run_upload_files)
    def run_cleanup():
        #p.cleanup()
        #p.log(action="%s.done" % p.name)
        pass

    if args.dry:
        pipeline_printout(sys.stdout, [run_cleanup],  gnu_make_maximal_rebuild_mode=args.recent_tasks_only, verbose=9)
    else:
        pipeline_run([run_cleanup], verbose=9,  gnu_make_maximal_rebuild_mode=args.recent_tasks_only, multiprocess=8)
