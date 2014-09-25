import sys
import os
import pysam
from collections import defaultdict
import subprocess
from ruffus import *
import pandas as pd
from starpipe.Pipeline import Pipeline
import argparse
import NdarManifestReader as ndar

class QC_Pipeline(Pipeline):
    """
    QC_Pipeline

    Arguments:
        base_dir (required) Path to a local directory in which temporary files,
            reference files and log files will be stored.

        n_threads (default=1) # of threads to use.
    """
    def __init__(self, base_dir, n_threads, remove_intermediate):
        super(QC_Pipeline, self).__init__(base_dir=base_dir,
                                           n_threads=n_threads,
                                           remove_intermediate=remove_intermediate)

        self.name = "QC"

        self.files = {
            "qplot_reference": self.as_ref("human.g1k.v37.umfa"),
            "qplot_dbsnp": self.as_ref("dbSNP130.UCSC.coordinates.tbl"),
            "qplot_gc_data": self.as_ref("human.g1k.v37.umfa.winsize100.gc"),
            "qplot_target_regions": self.as_ref("nimblegen_solution_V2refseq2010.HG19.bed"),
        }
        
        self.cmds = {
            "qplot": "/usr/local/bin/qplot",
            "s3cmd": "s3cmd --no-progress --force",
        }

        self.logs = {
            "stdout": None,#self.create_logger("stdout", self.as_out("init.stdout.log")),
            "stderr": None, #self.create_logger("stderr", self.as_out("init.stderr.log")),
            #"pipeline": self.create_logger("pipeline", self.as_out("bwa.pipeline.log")),
            "sims": self.create_sims_logger(remote_host='master'),
        }


    def qplot(self, bam_files, out_prefix, labels):
        self.cmd("{qplot} \
            --reference {qplot_reference} \
            --dbsnp {qplot_dbsnp} \
            --region {qplot_regions} \
            --plot {plot_file} \
            --stats {stats_file} \
            --Rcode {Rcode_file} \
            --xml {xml_file} \
            --bamLabel {labels} \
            {files}".
            format(
                qplot=self.cmds["qplot"],
                qplot_reference=self.files["qplot_reference"],
                qplot_dbsnp=self.files["qplot_dbsnp"],
                qplot_regions=self.files["qplot_target_regions"],
                plot_file=self.as_temp("%s.qplot.pdf" % out_prefix),
                stats_file=self.as_temp("%s.qplot.txt" % out_prefix),
                Rcode_file=self.as_temp("%s.qplot.R" % out_prefix),
                xml_file=self.as_temp("%s.qplot.xml" % out_prefix),
                labels=",".join(labels),
                files=" ".join(bam_files)
            ),
            shell=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--dry", action="store_true",
                        help="Print out task list only. Does not run commands.")
    parser.add_argument("manifest", type=str,
                        help="path to Ndar manifest file")
    parser.add_argument("familyID", type=str,
                        help="SSC FamilyID to download")    
    parser.add_argument("tempdir", type=str,
                        help="path to temporary directory")
    args = parser.parse_args()

    p = QC_Pipeline(args.tempdir, n_threads=8, remove_intermediate=False)
    p.name = "%s_QC" % args.familyID

    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.familyID
    }

    input_files = {}
    file_to_sample = {}
    manifest = ndar.NdarManifestReader(args.manifest)

    input_files = []
    samples = manifest.get_family(args.familyID)
    for s in samples:
        in_file = p.as_temp("%s.realigned.recal.bam" % s.sample_id)
        input_files.append(in_file)

    @merge(input_files, "%s.qplot.pdf" % args.familyID)
    def run_qplot(bam_files, dummy):
        p.qplot(bam_files, out_prefix=args.familyID, labels=[s.sample_id for s in samples])

    if args.dry:
        pipeline_printout(sys.stdout, [run_qplot], verbose=9)
    else:
        pipeline_run([run_qplot], verbose=9, multiprocess=8)



