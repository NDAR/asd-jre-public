import sys
import os
import pysam
import copy
from collections import defaultdict
import subprocess
import re
from ruffus import *
import pandas as pd
from starpipe.Pipeline import Pipeline, PipelineException
import argparse
import NdarManifestReader as ndar

class Init_Pipeline(Pipeline):
    """
    Init_Pipeline

    Arguments:
        base_dir (required) Path to a local directory in which temporary files,
            reference files and log files will be stored.

        n_threads (default=1) # of threads to use.
    """
    def __init__(self, base_dir, n_threads, remove_intermediate):
        super(Init_Pipeline, self).__init__(base_dir=base_dir,
                                           n_threads=n_threads,
                                           remove_intermediate=remove_intermediate)

        self.name = "Init"

        self.default_stage = self.download_bam_file_from_s3

        self.files = {
            "remote_reference_dir": "/data/reference/",
            "s3_reference_dir": "s3://asdjre/REFERENCE/"
        }
        
        self.cmds = {
            "rsync": "rsync --recursive --perms --group --verbose --compress --update",
            "s3cmd": "s3cmd --no-progress --force",
            "samtools": "samtools",
            "picard": "java -d64 -Xmx10g -XX:ParallelGCThreads=4 -jar /usr/local/bin/picard",
        }

        self.logs = {
            "stdout": None,#self.create_logger("stdout", self.as_out("init.stdout.log")),
            "stderr": None, #self.create_logger("stderr", self.as_out("init.stderr.log")),
            #"pipeline": self.create_logger("pipeline", self.as_out("bwa.pipeline.log")),
            "sims": self.create_sims_logger(remote_host='master'),
        }

    def setup_directories(self):
        """
        Create local temporary directories:
            <base_dir>/temp/
            <base_dir>/reference/
            <base_dir>/out/
        """
        for directory in [self.local_temp_dir, self.local_out_dir, self.local_reference_dir]:
            if not os.path.exists(directory):
                self.cmd("mkdir -p %s" % directory, shell=True)
        
        self.logs["stdout"] = self.create_logger("stdout", self.as_out("init.stdout.log"))
        self.logs["stderr"] = self.create_logger("stderr", self.as_out("init.stderr.log"))

    def rsync_reference_files(self):
        """
        Download reference files if needed to self.local_reference_dir
        """
        self.cmd("{rsync_cmd} {remote} {local} && gunzip {local}/*.gz"
            .format(
                rsync_cmd=self.cmds["rsync"],
                remote=self.files["remote_reference_dir"],
                local=self.local_reference_dir,
            ),
            shell=True)

    def download_reference_file_from_s3(self, s3_file, local_file, unzip=False):
        """
        Use s3cmd and streaming gunzip to speed up download from s3
        """
        if unzip:
            rc = subprocess.call("{s3cmd} get {s3_file} - | gunzip > {local_file}"
                .format(
                    s3cmd=self.cmds["s3cmd"],
                    s3_file=s3_file,
                    local_file=local_file, 
                ), shell=True)
        else:
            rc = subprocess.call("{s3cmd} get {s3_file} {local_file}"
                .format(
                    s3cmd=self.cmds["s3cmd"],
                    s3_file=s3_file,
                    local_file=local_file,
                ), shell=True)
        return rc

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

    def rename_bam_files(self, old_name, new_name):
        self.cmd("mv {old} {new}"
            .format(
                old=old_name,
                new=new_name
                ),
            shell=True)

    def bam_preflight_check(self, bam_file, n_reads=1e6):
        """
        Do a series of quick checks on first n_reads of bam_file.
        Currently only checks flags for paired/unpaired status.
        Returns:
            - <dictionary> counts "flags" and "readlengths"
            - <int> total # of reads examined
        """
        bam_iter = pysam.Samfile(bam_file, 'rb', check_header=False, check_sq=False).fetch(until_eof=True)
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

    def modify_bam_header(self, in_bam, out_bam):
        """
        This (for now) replaces the SM tag of the @RG line in the header
        with a new sample_id derived from the filename. 
        It moves the existing SM tag into a new @CO line at the end of the header
        """
        #bam_header = pysam.Samfile(in_bam,'rb',check_header=False, check_sq=False).header
        bam_header_raw = pysam.Samfile(in_bam,'rb',check_header=False, check_sq=False).text.replace("\t\n","\n")
        temp_header = in_bam + ".tempheader"
        with open(temp_header ,"w") as f:
            f.write(bam_header_raw)

        bam_header = pysam.Samfile(temp_header,'r', check_header=False, check_sq=False).header
        sample_id = os.path.basename(in_bam).replace(".pre.bam", "")
        try:
            original_SM = list(set([x["SM"] for x in bam_header["RG"]]))[0]
        except:
            raise PipelineException("@RG header line not found in %s!" % bam_in)

        # make sure SM tags in RG line are consistent with sample_id
        rgs = copy.copy(bam_header["RG"])
        bam_header["RG"] = []
        for rg in rgs:
            rg["SM"] = sample_id
            bam_header["RG"].append(rg)

        # save original SM tage
        if "CO" not in bam_header:
            bam_header["CO"] = ["Original RG/SM tag: %s" % original_SM]
        else:
            bam_header["CO"].append("Original RG/SM tag: %s" % original_SM)

        # write out header
        header_filename = self.as_temp("%s.header" % in_bam)
        header_file = pysam.Samfile(header_filename, 'wh', header=bam_header)
        header_file.close()

        self.cmd("{samtools} reheader \
                  {header_file} \
                  {in_bam} > {out_bam}"
                .format(
                    samtools = self.cmds["samtools"],
                    in_bam=in_bam,
                    out_bam=out_bam,
                    header_file=header_filename,
                ),
            shell=True)

        self.rm(in_bam)

    # TODO, rename and rework this (merge+RG tag adding function)
    def merge_bam_files(self, inputs, output, sample_id, rg_id=None,
            platform='illumina', library='A', sort_order="readname"):
        """
        Add @RG tabs and merge multiple SORTED bam files into one output file. 
        If only providing a single file, file will be re-headered with @RG information and RG tags.
        Arguments:
            inputs : list of paths/bam files to merge (or just one file)
            output : path of output merged file
            sort_order : specify "readname" or "coordinate" for samtools merge order. 
                Note: sort_order ignored when only providing a single BAM file to merge
        """
        if len(inputs) > 1:
            if sort_order == "readname":
                sort_options = "-n"
            else:
                sort_options = ""
    
            header_file = p.as_temp("%s.header" % output)

            with open(header_file, "w") as header:
                for ix, input_file in enumerate(inputs):
                    # TODO use pysam here
                    in_header = pysam.Samfile(input_file,'rb',check_header=False, check_sq=False).text
                    RG_lines = filter(lambda x: x.startswith("@RG"), in_header.split("\n"))
                    if len(RG_lines) == 1:
                        rg_id = re.findall("ID:([a-zA-Z0-9_\-\.]*)", RG_lines[0])[0]
                    else:
                        rg_id = re.sub("\.bam$", "", os.path.basename(input_file))
                    header.write("@RG\tID:%s\tPU:%s\tDS:%s\tLB:%s\tPL:%s\tSM:%s\n" % (rg_id, rg_id, input_file, library, platform, sample_id))
            merge_options = "-h %s" % (header_file)

            self.cmd("{samtools} merge \
                      {sort_options} \
                      {merge_options} \
                      {output_bam} {input_bam_list}"
                    .format(
                        samtools=self.cmds["samtools"],
                        sort_options=sort_options,
                        merge_options=merge_options,
                        output_bam=output,
                        input_bam_list=" ".join(inputs),
                        ),
                shell=True)
        else:
            # TODO use pysam here
            input_file = inputs[0]
            in_header = pysam.Samfile(input_file,'rb',check_header=False, check_sq=False).text
            RG_lines = filter(lambda x: x.startswith("@RG"), in_header.split("\n"))
            if len(RG_lines) == 1:
                rg_id = re.findall("ID:([a-zA-Z0-9_\-\.]*)", RG_lines[0])[0]
            else:
                rg_id = re.sub("\.bam$", "", os.path.basename(input_file))
            with open(p.as_temp("%s.header" % output), "w") as header:
                header.write("@RG\tID:%s\tPU:%s\tDS:%s\tLB:%s\tPL:%s\tSM:%s\n" % (rg_id, rg_id, input_file, library, platform, sample_id))
            
            self.cmd("{picard}/AddOrReplaceReadGroups.jar \
                        INPUT={in_bam} \
                        OUTPUT={out_bam} \
                        QUIET=false \
                        VALIDATION_STRINGENCY=LENIENT\
                        COMPRESSION_LEVEL=5 \
                        RGID={rg_id} \
                        RGSM={sample_id} \
                        RGPU={rg_id} \
                        RGLB=A \
                        RGPL=illumina \
                        RGDS={in_bam}"
                .format(
                    picard=self.cmds["picard"],
                    in_bam=inputs[0],
                    out_bam=output,
                    sample_id=sample_id,
                    rg_id=rg_id,
                    ),
                shell=True)

    def bam_check_rg(self, inputs):
        """
        Check if all files in inputs have at least one @RG header line
        """
        for f in inputs:
            bam_header = subprocess.check_output("{samtools} view -H {in_bam}"
                .format(samtools="samtools", in_bam=f),
                shell=True)
            if len(filter(lambda x: x.startswith("@RG\t"), bam_header.split("\n"))) == 0:
                return False
        return True

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--dry", action="store_true",
                        help="Print out task list only. Does not run commands.")
    parser.add_argument("--recent-tasks-only", "--recent", action="store_false", default=True,
                        help="Only run most recent tasks necessary. Use with caution!"
                             "Corresponds to the gnu_make_maximal_rebuild_mode mode in ruffus.")        
    parser.add_argument("--map-se", action="store_true", default=False,
                        help="Skips preflight check for read pairedness" 
                             "(useful when enabling --map-se for BWA_Pipeline")
    parser.add_argument("manifest", type=str,
                        help="path to Ndar manifest file")
    parser.add_argument("familyID", type=str,
                        help="SSC FamilyID to download")    
    parser.add_argument("tempdir", type=str,
                        help="path to temporary directory")
    args = parser.parse_args()

    p = Init_Pipeline(args.tempdir, n_threads=8, remove_intermediate=False)
    p.name = "%s_Init" % args.familyID

    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.familyID
    }

    input_files = defaultdict(list)
    manifest = ndar.NdarManifestReader(args.manifest)

    for s in manifest.get_family(args.familyID):
        for i in s.files:
            f = p.as_temp("%s.pre.bam" % s.sample_id)
            input_files[f].append(i)

    print input_files
    # TODO, it would be best to order these largest to smallest
    ref_file_list = subprocess.check_output( \
            "s3cmd ls {s3_reference_dir} | awk '{{print $4}}'" \
            .format(s3_reference_dir=p.files["s3_reference_dir"]), shell=True) \
            .split()
    
    ref_files = {}
    for f in ref_file_list:
        ref_files[f] = p.as_ref(os.path.basename(f.rstrip(".gz")))


    def do_check_S3_download_finished(remote, local):
        # note this reverses True and False due to ruffus implementation
        if p.check_s3_download_finished(remote, local) == False:
            return True, "S3 download incomplete"
        else:
            return False, "S3 file already fully downloaded"

    def run_setup_directories():
        p.log(action="%s.start" % p.name)
        p.setup_directories()
    
    @follows(run_setup_directories)
    @check_if_uptodate(do_check_S3_download_finished)
    @parallel(zip(ref_files.keys(), ref_files.values()))
    @jobs_limit(4)
    def run_download_reference_file_from_s3(s3_file, local_file):
        # todo, use BOTO to get list of files here...
        unzip = s3_file[-3:] == ".gz"
        p.download_reference_file_from_s3(s3_file, local_file, unzip=unzip)

    @follows(run_download_reference_file_from_s3)
    @check_if_uptodate(do_check_S3_download_finished)
    @parallel(zip(input_files.values(), input_files.keys()))
    def run_download_bam_file_from_s3(s3_files, local_file):        
        local_file_parts = []
        for s3_file in s3_files:
            # download each part separately
            l = p.as_temp(os.path.basename(s3_file))
            p.download_bam_file_from_s3(s3_file, l)
            local_file_parts.append(l)

        if (len(local_file_parts) > 1) or (p.bam_check_rg(local_file_parts) == False):
            # now merge into the parts into the final file
            sample_id = re.sub("\.pre\.bam$", "", os.path.basename(local_file))
            p.merge_bam_files(local_file_parts, local_file, sample_id, sort_order="readname")
        else:
            # rename single bam file into local_file:
            p.rename_bam_files(local_file_parts[0], local_file)
    
        preflight_check, total_reads_checked = p.bam_preflight_check(local_file)

        # TODO there's a real potential for failure here, if the sample has 
        # multiple merged files, and some are SE-- this means that just reading the first 
        # million reads would /not/ discover the SE reads.
        #
        # unpaired reads are flags 0, 4 and 16
        unpaired_cnt = sum([preflight_check["flags"][f] for f in [0,4,16]])
        if (unpaired_cnt > (0.1 * total_reads_checked)) \
            and not args.map_se:
            p.log(action="%s.qc_failure" % p.name,\
                message="Preflight check failed for %s! (%d/%d reads unpaired)" % (local_file, unpaired_cnt, total_reads_checked),\
                flags=preflight_check["flags"], readlengths=preflight_check["readlengths"])
            sys.exit(2) # exit on status 2 to indicate known exit condition
        else:
            p.log(action="%s.qc_OK" % p.name, message="Preflight check OK",\
                flags=preflight_check["flags"], readlengths=preflight_check["readlengths"])


    
    @transform(run_download_bam_file_from_s3, suffix(".pre.bam"), r"\1.input.bam")
    def run_modify_bam_header(local_file, modified_bam):
        p.modify_bam_header(local_file, modified_bam)

    @follows(run_modify_bam_header)
    def finalize():
        p.log(action="%s.done" % p.name)

    if args.dry:
        pipeline_printout(sys.stdout, [finalize], gnu_make_maximal_rebuild_mode=args.recent_tasks_only, verbose=9)
    else:
        pipeline_run([finalize], verbose=9, gnu_make_maximal_rebuild_mode=args.recent_tasks_only, multiprocess=8)
