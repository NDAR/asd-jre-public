import sys
import os
import re
import subprocess
from ruffus import *
import pandas as pd
from starpipe.Pipeline import Pipeline
import argparse
from collections import defaultdict
import glob
import pysam
import NdarManifestReader as ndar
import multiprocessing

class BWA_Pipeline(Pipeline):
    """
    BWA_pipeline

    Arguments:
        input_files (required) List of local paths to bam files to process

        base_dir (required) Path to a local directory in which temporary files,
            reference files and log files will be stored.

        n_threads (default=1) # of threads to use.

        remove_intermediate (default=False) Remove intermediate files if True.
    """
    def __init__(self, base_dir, n_threads, remove_intermediate, **kwargs):
        super(BWA_Pipeline, self).__init__(base_dir=base_dir,
                                           n_threads=n_threads,
                                           remove_intermediate=remove_intermediate,
                                           **kwargs)

        self.name = "BWA"
        self.default_stage = self.picard_mark_duplicates

        self.files = {
            "reference_genome": self.as_ref("human_g1k_v37.fasta"),
            "reference_genome_contigs": self.as_ref("human_g1k_v37.contigs"),
        }

        self.cmds = {
            "picard": "java -d64 -Xmx10g -XX:ParallelGCThreads=4 -jar /usr/local/bin/picard",
            "picard_custom_limits": lambda max_mem, threads: "java -d64 -Xmx%dg -XX:ParallelGCThreads=%d -jar /usr/local/bin/picard" % (max_mem, threads),
            "readgroup_mover": "/data/asd-jre/Experiments/ruffus/readgroup_mover.py",
            "samtools": "/usr/bin/samtools",
            "samtools_sort_params": "-@8 -m5200000000",
            "bedtools": "bedtools",
            "s3cmd": "s3cmd --no-progress",
            "bamtools": "/data/bamtools/bamtools"
        }

        self.logs = {
            "stdout": self.create_logger("stdout", self.as_out("bwa.stdout.log")),
            "stderr": self.create_logger("stderr", self.as_out("bwa.stderr.log")),
            #"pipeline": self.create_logger("pipeline", self.as_out("bwa.pipeline.log")),
            "sims": self.create_sims_logger(remote_host='master'),
        }


    def bwa_align(self, input_bam, out_sais):
        """
        BWA aln command for reads 1 and 2.
        """
        for read_num in [1, 2]:
            self.cmd("bwa aln -t{threads} -b -{read_num} {reference_genome} {input_bam} > {out_sai}"
                .format(
                    threads=self.n_threads,
                    read_num=read_num,
                    reference_genome=self.files["reference_genome"],
                    input_bam=input_bam,
                    out_sai=out_sais[read_num-1]  # account for indexing
                ),
                on_error=lambda: self.create_error_file(out_sais[read_num-1]),
                shell=True)

    def create_rg_dict(self, bam_in, rg_dict_out_files, rg_prefix=False):
        """
        Create readgroup dictionary files. These files are used downstream
        in order to re-assign original RG tags to reads in re-mapped files
        """
        if rg_prefix:
            sampleID = os.path.basename(bam_in).rstrip(".input.bam")
            prefix_string = "--prefix %s" % (sampleID)
        else:
            prefix_string = ""

        self.cmd("python {readgroup_mover} create\
                {prefix_string} \
                --input {bam_in}\
                --output {dict_out}"
            .format(
                readgroup_mover=self.cmds["readgroup_mover"],
                prefix_string=prefix_string,
                bam_in=bam_in,
                dict_out=rg_dict_out_files[0]
            ),
            on_error=lambda: self.create_error_file(rg_dict_out_files[0]),
            shell=True)

        self.checkpoint(rg_dict_out_files[0])
        self.checkpoint(rg_dict_out_files[1])
        self.checkpoint(rg_dict_out_files[2])

    def bwa_sampe(self, files_in, bam_out):
        """
        BWA sampe step. Merges *.sai files into a bam file. Uses RG dict to
        re-append correct RGs based on readnames.
        """
        self.cmd("bwa sampe -n 0 -N 0 -P \
                {fasta} {sai_in_files} '{input_bam}' '{input_bam}'\
                | python {readgroup_mover} translate --dictfile {rg_dict} \
                | {samtools} view -b -S > {bam_out}"
            .format(
                readgroup_mover=self.cmds["readgroup_mover"],
                fasta=self.files["reference_genome"],
                sai_in_files="%s %s" % (files_in[0], files_in[1][0]),
                input_bam=files_in[1][2],
                rg_dict=files_in[1][1],
                samtools=self.cmds["samtools"],
                bam_out=bam_out),
            on_error=lambda: self.create_error_file(bam_out),
            shell=True)
        if self.remove_intermediate:
            self.rm(sai_in)

    def bwa_mem_align(self, files_in, bam_out, bwa_mem_seed=19):
        """
        BWA MEM alignment. Sorts, then streams a bam file into bwa-mem and converts
        it back to a bam file at then end.
        """
        self.cmd("{samtools} sort -n -o {samtools_sort_params} {bam_in} {bam_in}.sorted.\
                    | {bedtools} bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout \
                    | bwa mem -p -M -t {threads} -k {bwa_mem_seed} -v 1 {ref_file} - \
                    | python {readgroup_mover} translate --dictfile {rg_dict} \
                    | {samtools} view -@4 -b -S - > {bam_out}".
                format(
                    samtools = self.cmds["samtools"],
                    samtools_sort_params=self.cmds["samtools_sort_params"],
                    bedtools = self.cmds["bedtools"],
                    readgroup_mover=self.cmds["readgroup_mover"],
                    bam_in=files_in[0],
                    threads=self.n_threads,
                    bwa_mem_seed=bwa_mem_seed,
                    ref_file=self.files["reference_genome"],
                    rg_dict=files_in[1],
                    bam_out=bam_out,
                ),
                on_error=lambda: self.create_error_file(bam_out),
                checkpoint_file=bam_out,
                shell=True)

    def split_bam_by_rg(self, bam_in):
        """
        Helper function to split a bam file into individual readgroup bams.
        Not normally neccessary given the readgroup_mover.py approach, but
        needed in cases where PE and SE libraries are mixed in the same bam.
        Run when --map-se is set in command line.
        Requires @RG header line(s).
        If only a single RG is present, command saves time by moving file instead.
        """
        bam_header = pysam.Samfile(bam_in,'rb',check_header=False, check_sq=False).header
        split_bams = []
        read_group_ids = []

        if len(bam_header["RG"]) == 0:
            raise Exception("Bam file %s does not contain any RG header lines!" % bam_in)
        elif len(bam_header["RG"]) == 1:
            #  no need to split anything, only one RG in this bam file
            rg = bam_header["RG"][0]["ID"]
            bam_split = re.sub("bam$","%s.bam" % rg, bam_in)
            self.cmd("mv {bam_in} {bam_split}"
                .format(
                    bam_in=bam_in,
                    bam_split=bam_split,
                ), 
                shell=True)
            split_bams.append(bam_split)
            read_group_ids.append(rg)
        else:
            for rg_line in bam_header["RG"]:
                rg = rg_line["ID"]
                bam_split = re.sub("bam$","%s.bam" % rg, bam_in)
                self.cmd("samtools view -1 -r {rg} {bam_in} > {bam_split}".
                    format(
                        bam_in=bam_in,
                        rg=rg,
                        bam_split=bam_split,
                    ),shell=True)
                split_bams.append(bam_split)
                read_group_ids.append(rg)

        return split_bams, read_group_ids

    def bam_preflight_check(self, bam_file, n_reads=1e6):
        """
        Do a series of quick checks on first n_reads of bam_file.
        Currently only checks flags for paired/unpaired status.
        Returns:
            - <dictionary> counts "flags" and "readlengths"
            - <int> total # of reads examined
        """
        bam_iter = pysam.Samfile(bam_file, 'rb').fetch(until_eof=True)
        flags = defaultdict(int)
        readlengths = defaultdict(int)
        ix = 0
        while ix < n_reads:
            try:
                read = bam_iter.next()
                ix += 1
            except StopIteration:
                return False
            flags[read.flag] += 1
            readlengths[len(read.seq)] += 1

        return {"flags": flags, "readlengths": readlengths}, ix

    def bwa_mem_align_by_rg(self, bam_in, bam_out, rg, pairing="auto", bwa_mem_seed=19):
        if pairing == "auto":
            preflight_check, total_reads_checked = self.bam_preflight_check(bam_in)
            #  TODO we need to check the actual binary flag here, not the int representation
            #  Also this should probably ensure that there is no "mixed" PE/SE situation.
            #  (i.e., unpaired_cnt == total_cnt || paired_cnt == total_cnt)
            unpaired_cnt = sum([preflight_check["flags"][f] for f in [0,4,16]])
            if unpaired_cnt > 0:
                _pairing = "se"
            else:
                _pairing = "pe"
        else:
            _pairing = pairing
        if _pairing == "pe":
            cmd_str = "{samtools} sort -n -o {samtools_sort_params} {bam_in} {bam_in}.sorted.\
                        | {bedtools} bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout \
                        | bwa mem -p -M -t {threads} -R\"{readgroup_str}\" -k {bwa_mem_seed} -v 1 {ref_file} - \
                        | {samtools} view -b -S - > {bam_out}"
        elif _pairing == "se":
            cmd_str = "{bedtools} bamtofastq -i {bam_in} -fq /dev/stdout \
                        | bwa mem -M -t {threads} -R\"{readgroup_str}\" -k {bwa_mem_seed} -v 1 {ref_file} - \
                        | {samtools} view -b -S - > {bam_out}"
        else:
            raise Exception("incorrect pairing argument to bwa_mem_align_by_rg (pe|se|auto)")

        # get the readgroup header line based on the RG id from arguments
        readgroup_str = filter(lambda x: "ID:%s" % rg in x, pysam.Samfile(bam_in, 'rb').text.split("\n"))[0]

        self.cmd(cmd_str.format(
                samtools=self.cmds["samtools"],
                samtools_sort_params=self.cmds["samtools_sort_params"],
                bedtools=self.cmds["bedtools"],
                bam_in=bam_in,
                threads=self.n_threads,
                bwa_mem_seed=bwa_mem_seed,
                readgroup_str=readgroup_str.encode('string-escape'),
                ref_file=self.files["reference_genome"],
                bam_out=bam_out,
            ),
            on_error=lambda: self.create_error_file(bam_out),
            checkpoint_file=bam_out,
            shell=True)
    
    def concat_bam(self, bams_in, bam_out, header_file=None):
        """
        Concatenate a set of bams, or just rename if only one file
        TODO: improve this function (header_file inconsistency)
        """
        if len(bams_in) == 1:
            self.cmd("mv {bam_in} {bam_out}"
                .format(
                    bam_in=bams_in[0],
                    bam_out=bam_out,
                    ), shell=True)
        else:
            self.cmd("{samtools} cat \
                    -h {header_file} \
                    -o {bam_out} \
                    {input_bams}".
                    format(
                        samtools=self.cmds["samtools"],
                        header_file=header_file,
                        bam_out=bam_out,
                        input_bams=" ".join(bams_in),
                        ),
                shell=True)

    def make_clean_bam(self, bam_in, bam_out):
        """
        Convert SAM to BAM and clean using Picard
        """
        self.cmd("{picard_cmd}/CleanSam.jar\
                INPUT='{bam_in}'\
                OUTPUT='/dev/stdout'\
                QUIET=false\
                VALIDATION_STRINGENCY=LENIENT\
                COMPRESSION_LEVEL=5 \
                > {bam_out}"
            .format(
                bam_in=bam_in,
                picard_cmd=self.cmds["picard"],
                bam_out=bam_out,
            ), 
            on_error=lambda: self.create_error_file(bam_out),
            checkpoint_file=bam_out,
            shell=True)
        if self.remove_intermediate:
            self.rm(bam_in)

    def sort_bam(self, bam_in, bam_sorted_out):
        """
        Sort bam using Picard
        """
        self.cmd("{picard_cmd}/SortSam.jar\
                INPUT='{bam_in}'\
                OUTPUT='{bam_sorted_out}'\
                SORT_ORDER=coordinate CREATE_MD5_FILE=false\
                CREATE_INDEX=false MAX_RECORDS_IN_RAM=5000000\
                VALIDATION_STRINGENCY=LENIENT\
                QUIET=false COMPRESSION_LEVEL=5 TMP_DIR='{local_temp_dir}'"
            .format(
                picard_cmd=self.cmds["picard"],
                bam_in=bam_in,
                bam_sorted_out=bam_sorted_out,
                local_temp_dir=self.local_temp_dir,
                ),
            on_error=lambda: self.create_error_file(bam_sorted_out),
            shell=True)
        if self.remove_intermediate:
            self.rm(bam_in)

    def fix_mates(self, bam_in, fixed_bam_out):
        """
        Fix mate information using Picard
        """
        self.cmd("{picard_cmd}/FixMateInformation.jar \
                INPUT='{bam_in}' \
                OUTPUT='/dev/stdout'\
                VALIDATION_STRINGENCY=LENIENT \
                SORT_ORDER=coordinate \
                MAX_RECORDS_IN_RAM=10000000 \
                TMP_DIR={local_temp_dir} \
                > {fixed_bam_out}"
            .format(
                picard_cmd=self.cmds["picard_custom_limits"](60, 8),
                bam_in=bam_in,
                fixed_bam_out=fixed_bam_out,
                local_temp_dir=self.local_temp_dir
                ),
            on_error=lambda: self.create_error_file(fixed_bam_out),
            checkpoint_file=fixed_bam_out,
            shell=True)
        if self.remove_intermediate:
            self.rm(bam_in)

    def picard_mark_duplicates(self, bam_in, files_out):
        """
        Mark duplicates and re-index bam
        """
        self.cmd("{picard_cmd}/MarkDuplicates.jar \
                INPUT='{bam_in}' \
                OUTPUT='/dev/stdout' \
                METRICS_FILE={metricfile_out} \
                REMOVE_DUPLICATES=false ASSUME_SORTED=true COMPRESSION_LEVEL=5 \
                VALIDATION_STRINGENCY=LENIENT \
                MAX_RECORDS_IN_RAM=5000000 \
                CREATE_INDEX=false \
                TMP_DIR={local_temp_dir} \
                > {bam_out}"
            .format(
                picard_cmd=self.cmds["picard"],
                bam_in=bam_in,
                bam_out=files_out[0],
                metricfile_out=files_out[1],
                local_temp_dir=self.local_temp_dir
                ),
            on_error=lambda: self.create_error_file(files_out),
            checkpoint_file=files_out[0],
            shell=True)

        self.checkpoint(files_out[1])

        if self.remove_intermediate:
            self.rm(bam_in)

    def qc_metrics(self, files_in, qc_files):
        """
        Calculate QC metrics on bam file
        """
        self.cmd("{samtools} index {bam_in}"
            .format(
                samtools=self.cmds["samtools"],
                bam_in=files_in[0],
            ),
            shell=True)
        self.cmd("{samtools} idxstats {bam_in} | tee {qc_file}"
            .format(
                samtools=self.cmds["samtools"],
                bam_in = files_in[0],
                qc_file = qc_files[0],
            ),
            shell=True,
            log_output=True)
        self.cmd("{samtools} flagstat {bam_in} | tee {qc_file}"
            .format(
                samtools=self.cmds["samtools"],
                bam_in = files_in[0],
                qc_file = qc_files[1],
            ),
            shell=True,
            log_output=True)
        
        self.checkpoint(qc_files[0])
        self.checkpoint(qc_files[1])
        self.checkpoint(qc_files[2])

    def cleanup_intermediate_files(self):
        """
        Cleans up all intermediate files
        """
        self.cmd("rm -f {local_temp_dir}/*rg_dict* \
                        {local_temp_dir}/*aln* \
                        {local_temp_dir}/snappy*".
            format(
                local_temp_dir=self.local_temp_dir
            ),
            shell=True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--dry", action="store_true",
                        help="Print out task list only. Does not run commands.")
    parser.add_argument("--recent-tasks-only", "--recent", action="store_false", default=True,
                        help="Only run most recent tasks necessary. Use with caution!"
                             "Corresponds to the gnu_make_maximal_rebuild_mode mode in ruffus.")    
    parser.add_argument("--bwa-mem", action="store_true", default=False,
                        help="Use bwa-mem")
    parser.add_argument("--map-se", action="store_true", default=False,
                        help="Map single-end reads")
    parser.add_argument("--prefix-rg", action="store_true", default=False, help="prefix the readgroup IDs with the sample name (from filename)")
    parser.add_argument("--bwa-seed", type=int, default=19,
                        help="Set the bwa initial seed (bwa -k <int>). Default 19 per BWA mem.")
    parser.add_argument("manifest", type=str,
                        help="path to NDAR manifest file")
    parser.add_argument("familyID", type=str,
                        help="SSC FamilyID to map")    
    parser.add_argument("tempdir", type=str,
                        help="path to temporary directory")
    parser.add_argument("--checkpoint-path", type=str, default=None,
                        help="path to checkpoint directory or S3 bucket")
    args = parser.parse_args()
    print args

    n_threads = multiprocessing.cpu_count()
    p = BWA_Pipeline(args.tempdir, n_threads=n_threads, remove_intermediate=False, checkpoint_path=args.checkpoint_path)
    if args.bwa_mem:
        p.name = "%s_BWA_mem" % args.familyID
    else:
        p.name = "%s_BWA_aln" % args.familyID

    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.familyID
    }

    manifest = ndar.NdarManifestReader(args.manifest)
    
    samples = manifest.get_family(args.familyID)

    input_files = []
    for s in samples:
        in_file = p.as_temp("%s.input.bam" % s.sample_id)
        input_files.append(in_file)

    def start():
        p.log(action="%s.start" % p.name)

    if args.bwa_mem:
        if args.map_se:
            @follows(start)
            @jobs_limit(1)
            @transform(input_files, suffix(".input.bam"), r"\1.aln.bam")
            def run_bwa(input_file, output_file):
                in_split_files,readgroups=p.split_bam_by_rg(input_file)
                out_split_files = []
                for f,rg in zip(in_split_files, readgroups):
                    bam_out = f + ".split.aln.bam"
                    p.bwa_mem_align_by_rg(f, bam_out, rg=rg, bwa_mem_seed=args.bwa_seed)
                    out_split_files.append(bam_out)

                # todo, write a header merge function
                # prepare a new header, merging @RG lines
                header_filename=input_file + ".header"
                # get the header from the first split file
                header = pysam.Samfile(out_split_files[0], 'rb').text.split("\n")
                # then add additional RG's from the other files
                for f in out_split_files[1:]:
                    header.append(filter(lambda x: x.startswith("@RG"), pysam.Samfile(f, 'rb').text.split("\n"))[0])
                # filter out empty lines which causes samtools to stop reading the header file
                header = filter(lambda x: (x != "\n") & (x != ""), header)
                # write to a file
                with open(header_filename, "w") as header_file:
                    header_file.write('\n'.join(header) + '\n')
                # concatenate the bam, supplying merged header
                p.concat_bam(out_split_files, output_file, header_file=header_filename)
        else:
            @follows(start)
            @transform(input_files, suffix(".input.bam"), [r"\1.rg_dict",r"\1.rg_dict.header",r"\1.rg_dict.lookups"])
            @check_if_uptodate(p.check_output_with_checkpoint)
            @jobs_limit(2)
            def run_create_rg_dict(bam_in, rg_dict_out):
                p.create_rg_dict(bam_in, rg_dict_out, prefix_rg=args.prefix_rg)
            
            @jobs_limit(1)
            @follows(run_create_rg_dict)
            @transform(input_files, suffix(".input.bam"), add_inputs(r"\1.rg_dict"), r"\1.aln.bam")
            @check_if_uptodate(p.check_output_with_checkpoint)
            def run_bwa(files_in, bam_out):
                p.bwa_mem_align(files_in, bam_out, bwa_mem_seed=args.bwa_seed)

    else:
        if args.map_se:
            raise Exception("Not Implemented")
        @split(input_files, regex(r"(.+)\.input.bam"), [r'\1.aln.1.sai', r'\1.aln.2.sai'])
        @jobs_limit(1)
        def run_bwa_align(input_bam, out_sais):
            p.bwa_align(input_bam, out_sais)

        @follows(run_create_rg_dict)
        @transform(run_bwa_align, suffix(".aln.1.sai"), add_inputs([r"\1.aln.2.sai", r"\1.rg_dict", r"\1.input.bam"]), r'\1.aln.bam')
        @jobs_limit(2)
        def run_bwa(files_in, bam_out):
            p.bwa_sampe(files_in, bam_out)

    @jobs_limit(4)
    @transform(run_bwa, suffix(".aln.bam"), r"\1.aln.clean.bam")
    @check_if_uptodate(p.check_output_with_checkpoint)
    def run_make_clean_bam(bam_in, bam_out):
        p.make_clean_bam(bam_in, bam_out)

    @jobs_limit(1)
    @transform(run_make_clean_bam, suffix(".bam"), r"\1.sorted.mates_fixed.bam")
    @check_if_uptodate(p.check_output_with_checkpoint)
    def run_fix_mates(bam_in, fixed_bam_out):
        p.fix_mates(bam_in, fixed_bam_out)

    @jobs_limit(4)
    @transform(run_fix_mates, suffix(".aln.clean.sorted.mates_fixed.bam"), [r"\1.out.bam", r"\1.duplicate_metrics.txt"])
    @check_if_uptodate(p.check_output_with_checkpoint)
    def run_picard_mark_duplicates(bam_in, files_out):
        p.picard_mark_duplicates(bam_in, files_out)

    @jobs_limit(4)
    @transform(run_picard_mark_duplicates, suffix(".out.bam"), [r"\1.idxstats.txt",r"\1.flagstat.txt",r"\1.out.bam.bai"] )
    @check_if_uptodate(p.check_output_with_checkpoint)
    def run_qc_metrics(files_in, qc_files):
        p.qc_metrics(files_in, qc_files)

    @follows(run_qc_metrics)
    def run_cleanup_files():
        p.cleanup_intermediate_files()
        p.log(action="%s.done" % p.name)
    
    if args.dry:
        pipeline_printout(sys.stdout, [run_cleanup_files], gnu_make_maximal_rebuild_mode=args.recent_tasks_only, verbose=9)
    else:
        pipeline_run([run_cleanup_files], verbose=9, gnu_make_maximal_rebuild_mode=args.recent_tasks_only, multiprocess=8)
