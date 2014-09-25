import sys
import os
import glob
import subprocess
from ruffus import *
from starpipe.Pipeline import Pipeline
import argparse
import NdarManifestReader as ndar


class Cleanup_PipelineException(Exception):
    pass

class Cleanup_Pipeline(Pipeline):
    """
    Cleanup_Pipeline

    Arguments:
        base_dir (required) Path to a local directory in which temporary files,
            reference files and log files will be stored.

        n_threads (default=1) # of threads to use.
    """
    def __init__(self, base_dir, n_threads):
        super(Cleanup_Pipeline, self).__init__(base_dir=base_dir, n_threads=n_threads)

        self.cmds = {
            "rsync": "rsync --recursive --perms --group --verbose --compress --update",
            "s3cmd": "s3cmd -c /root/.s3cfg --no-progress --force",
        }

        self.logs = {
            "stdout": None,#self.create_logger("stdout", self.as_out("init.stdout.log")),
            "stderr": None, #self.create_logger("stderr", self.as_out("init.stderr.log")),
            #"pipeline": self.create_logger("pipeline", self.as_out("bwa.pipeline.log")),
            "sims": self.create_sims_logger(remote_host='master'),
        }
    
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
            raise Cleanup_PipelineException(
                "S3 upload of local file %s failed to destination %s" % (local_file, remote_path))

    def list_bucket_files(self, familyID):
        """
        get listing of all files in s3 bucket
        """
        self.cmd("{s3cmd} ls -r s3://asdjre/complete/{family_id}"
            .format(
                s3cmd=self.cmds["s3cmd"],
                family_id=familyID
            ),
            log_output=True,
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
    args = parser.parse_args()

    p = Cleanup_Pipeline(args.tempdir, n_threads=4)
    
    p.name = "%s_Cleanup" % args.familyID

    p.metadata = {
        "project_id": "default_project",
        "sample_id": args.familyID
    }
    
    manifest = ndar.NdarManifestReader(args.manifest)
    samples = manifest.get_family(args.familyID)
    
    s3_bucket_path = "s3://asdjre/complete/%s" % args.familyID

    files = {#p.as_temp("{sample_id}.realigned.recal.reduced.bam"):
             #   "{s3_bucket_path}/reduced_bams/{sample_id}.realigned.recal.reduced.bam",
             #p.as_temp("{sample_id}.realigned.recal.reduced.bai"):
             #   "{s3_bucket_path}/reduced_bams/{sample_id}.realigned.recal.reduced.bai",
             p.as_temp("{sample_id}.realigned.recal.bam"):
                "{s3_bucket_path}/complete_bams/{sample_id}.realigned.recal.bam",
             p.as_temp("{sample_id}.realigned.recal.bai"):
                "{s3_bucket_path}/complete_bams/{sample_id}.realigned.recal.bai",
             p.as_temp("{sample_id}.flagstat.txt"):
                "{s3_bucket_path}/qc/{sample_id}.flagstat.txt",
             p.as_temp("{sample_id}.idxstats.txt"):
                "{s3_bucket_path}/qc/{sample_id}.idxstats.txt",
             p.as_temp("{sample_id}.duplicate_metrics.txt"):
                "{s3_bucket_path}/qc/{sample_id}.duplicate_metrics.txt",
             p.as_temp("*.qplot.*"):
                "{s3_bucket_path}/qc/{basename}",
             p.as_out("*.log"):
                "{s3_bucket_path}/log/{basename}",
             }

    files_to_upload = {}    
    for s in samples:
        for local_path_template, remote_path_template in files.iteritems():
            for local_path in glob.glob(local_path_template.format(sample_id=s.sample_id)):
                remote_path = remote_path_template.format(
                        sample_id=s.sample_id,
                        s3_bucket_path=s3_bucket_path,
                        basename=os.path.basename(local_path))
                files_to_upload[local_path] = remote_path
    
    print files_to_upload

    def do_check_S3_upload_finished(local, remote):
        # note this reverses True and False due to ruffus implementation
        if p.compare_s3_filesize(remote, local) == False:
            return True, "S3 upload incomplete"
        else:
            return False, "S3 file fully uploaded"

    def start():
        p.log(action="%s.start" % p.name)

    @follows(start)
    @parallel(zip(files_to_upload.keys(), files_to_upload.values()))
    @check_if_uptodate(do_check_S3_upload_finished)
    def run_upload_to_s3(local_file, remote_path):
        p.upload_to_s3(local_file, remote_path)

    @follows(run_upload_to_s3)
    def run_list_bucket_files():
        p.list_bucket_files(args.familyID)

    @follows(run_list_bucket_files)
    def run_cleanup():
        p.cleanup()
        p.log(action="%s.done" % p.name)

    if args.dry:
        pipeline_printout(sys.stdout, [run_cleanup],  gnu_make_maximal_rebuild_mode=args.recent_tasks_only, verbose=9)
    else:
        pipeline_run([run_cleanup], verbose=9,  gnu_make_maximal_rebuild_mode=args.recent_tasks_only, multiprocess=8)
