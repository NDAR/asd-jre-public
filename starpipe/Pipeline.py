import subprocess
import sys
import signal
import time
import os
import collections
import re
import ruffus
from ruffus.proxy_logger import *

try:
    import sims.SimsConnection
    sims_available = True
except ImportError:
    print "Could not import sims.SimsConnection package. Disabling!"
    sims_available = False


class PipelineException(Exception):
    pass


class Pipeline(object):
    """A Ruffus pipeline with extra pizzaz"""
    def __init__(self, base_dir, n_threads=1, remove_intermediate=False,
                 default_stage=None, checkpoint_path=None, dry_run=False):
        super(Pipeline, self).__init__()
        # Set up signal termination handles
        signal.signal(signal.SIGINT, self.termination_handler)
        signal.signal(signal.SIGTERM, self.termination_handler)

        self.name = "DefaultPipeline"
        self.metadata = {
            "project_id": "default_project",
            "sample_id": "default_sample"
        }
        # set up directories
        self.base_dir = base_dir
        self.local_out_dir = os.path.join(self.base_dir, "out")
        self.local_temp_dir = os.path.join(self.base_dir, "temp")
        self.local_reference_dir = os.path.join(self.base_dir, "reference")

        # runtime options
        self.n_threads = n_threads
        self.remove_intermediate = remove_intermediate
        self.dry_run = dry_run

        # set up checkpointing if desired
        self.checkpoint_path = checkpoint_path
        if self.checkpoint_path is not None:
            self.checkpoint_path = checkpoint_path.rstrip("/")
            if checkpoint_path[0:5] == "s3://":
                # checkpoint to amazon bucket using s3cmd
                self.stream_checkpoint = \
                    lambda filename: ">(tee %s) | s3cmd --no-progress --force put - %s/%s" % \
                    (filename, checkpoint_path, os.path.basename(filename))
            elif os.path.exists(checkpoint_path):
                self.stream_checkpoint = \
                    lambda filename: os.path.join(checkpoint_path, os.path.basename(filename))
            else:
                raise PipelineException("Could not find checkpoint path (%s)" % checkpoint_path)

        # set up default loggers
        self.logs = {
            "stdout": sys.stdout,
            "stderr": sys.stderr,
        }

        self._excluded_input_files = []
        self.log()

    def as_temp(self, file):
        return os.path.join(self.local_temp_dir, file)

    def as_out(self, file):
        return os.path.join(self.local_out_dir, file)

    def as_ref(self, file):
        return os.path.join(self.local_reference_dir, file)

    def cmd(self, cmd_string, shell=False, verbose=True, 
        on_error=None, on_success=None, stop_on_error=True, 
        log_output=False, log_command=True, checkpoint_file=False, n_retry=0):
        t1 = time.time()
        if checkpoint_file and self.checkpoint_path:
            if log_output:
                raise PipelineException("Cannot checkpoint_file and log_output in same command!")
            if cmd_string.count(checkpoint_file) != 1:
                raise PipelineException("Cannot find checkpoint_file (%s)"
                                        " in cmd_string!" % checkpoint_file)
            cmd_string = cmd_string.replace(checkpoint_file, self.stream_checkpoint(checkpoint_file))
        if not shell:
            arglist = cmd_string.split(" ")
        else:
            arglist = cmd_string
        cmd_string = re.sub(r"\s{5,}", '\n\t', cmd_string)  # this replaces ugly line breaks in commands!
        if verbose:
            print "    Command:", cmd_string
        if "pipeline" in self.logs:
            with self.logs["pipeline"]["mutex"]:
                self.logs["pipeline"]["log"].info("\tCommand:", cmd_string)
        n_tries = 0
        if log_output:
            while n_tries <= n_retry:
                n_tries += 1
                p = subprocess.Popen(arglist,
                                 shell=shell,
                                 stdout=subprocess.PIPE)
                output = p.communicate()[0]
                return_code = p.returncode
                if return_code == 0:
                    break
        else:
            while n_tries <= n_retry:
                n_tries += 1
                return_code = subprocess.call(arglist,
                                          shell=shell,
                                          stdout=self.logs["stdout"],
                                          stderr=self.logs["stderr"],
                                          executable='/bin/bash') # required for checkpointing magic)
                output = None
                if return_code == 0:
                    break

        if return_code != 0:
            if on_error is not None:
                on_error()

            error_msg = "Fatal exception:\n\tCommand: %s\n\tReturn code: %d" % (cmd_string, return_code)

            self.log(action="%s.error" % self.name,
                     command=cmd_string,
                     runtime=time.time() - t1,
                     output=output,
                     return_code=return_code)

            if stop_on_error:
                raise ruffus.JobSignalledBreak(error_msg)
            else:
                # do not stop the pipeline
                return False
        else:
            if on_success is not None:
                on_success()
            if log_command:
                self.log(action="%s.success" % self.name,
                         command=cmd_string,
                         runtime=time.time() - t1,
                         output=output)
                return True

    def create_sims_logger(self, remote_host, remote_port=8080):
        return sims.SimsConnection.SimsConnection(remote_host, remote_port)

    def create_logger(self, logger_type, log_filename):
        if logger_type == 'pipeline':
            logger_proxy, logging_mutex = make_shared_logger_and_proxy(setup_std_shared_logger,
                                                                       "pipeline_log",
                                                                       {"file_name": log_filename})
            return {"log": logger_proxy, "mutex": logging_mutex}
        elif logger_type in ["stderr", "stdout"]:
            return open(log_filename, 'w')


    def rm(self, infile):
        if isinstance(infile, collections.Iterable) and not isinstance(infile, basestring):
            for f in infile:
                os.unlink(f)
        else:
            os.unlink(infile)

    def create_error_file(self, infile, suffix='.error'):
        if not isinstance(infile, collections.Iterable) and not isinstance(infile, basestring):
            infile = [infile]

        for f in infile:
            error_filename = "%s.%s" % (f, suffix)
            if not os.path.exists(error_filename):
                try:
                    os.rename(f, error_filename)
                except OSError:
                    pass
            else:
                i = 1
                error_filename = "%s%s.%d" % (f, suffix, i)
                while os.path.exists(error_filename):
                    error_filename = "%s%s.%d" % (f, suffix, i)
                try:
                    os.rename(f, error_filename)
                except OSError:
                    pass

    def termination_handler(self, signal, frame):
        """
        Runs if pipeline receives a SIGTERM
        """
        self.log(action="%s.Terminated" % self.name,
                 message="Pipeline received SIGTERM (%d) and aborted." % signal,
                 line_number=frame.f_lineno,
                 local_vars=frame.f_locals)
                 #traceback=frame.f_trace)
        sys.exit(0)


    def log(self, action=None, project_id=None, sample_id=None, **kwargs):
        if action is None:
            action = "%s.info" % self.name
        if project_id is None:
            project_id = self.metadata["project_id"]
        if sample_id is None:
            sample_id = self.metadata["sample_id"]

        if ("sims" in self.logs) and sims_available:
            d = {}
            for key, val in kwargs.iteritems():
                #{"command": cmd_string, "runtime": runtime}
                if val is not None:
                    d[key] = val
            self.logs["sims"].record_action(
                project_id=project_id,
                sample_id=sample_id,
                action=action, #"%s.success" % self.name,
                data=d)

    def checkpoint(self, local_files):
        # todo, use BOTO here
        if self.checkpoint_path is None:
            # no checkpoint_path specified
            return False
        else:
            if type(local_files) not in [list, tuple]:
                local_files = [local_files]
            for local_file in local_files:
                checkpoint_file = os.path.join(self.checkpoint_path, os.path.basename(local_file)) # todo, this could be a self.as_checkpoint() func
                if self.checkpoint_path[0:5] == "s3://":
                    print "CHECKPOINTING %s --> %s" % (local_file, checkpoint_file)
                    subprocess.call("s3cmd --no-progress --force put %s %s" % (local_file, checkpoint_file), shell=True)
                elif os.path.exists(self.checkpoint_path):
                    subprocess.call("cp %s %s" % (local_file, checkpoint_file), shell=True)
                else:
                    raise PipelineException("Could not find checkpoint path (%s)" % self.checkpoint_path)


    # def check_output(self, input_files, output_files):
    #     if type(output_files) is not list:
    #         output_files = [output_files]
        
    #     missing = filter(lambda x: os.path.exists(x) is False, output_files)
    #     if len(missing) > 0:
    #         return True, "Missing files: %s" % ", ".join(missing)
    #     else:    
    #         return False, "All output files exists"

    def check_s3_file_exists(self, s3_path):
        # TODO, again use boto here
        file_count = subprocess.check_output("{s3cmd} ls {s3_path} | wc -l"
                                .format(
                                    s3cmd=self.cmds["s3cmd"],
                                    s3_path=s3_path,
                                ), shell=True)
        if int(file_count) >= 1: # s3cmd returns all files with prefix and longer
            return True
        else:
            return False

    def download_s3_file(self, s3_path, local_path):
        cmd_string = "{s3cmd} get {s3_path} {local_path}"\
                    .format(
                        s3cmd=self.cmds["s3cmd"],
                        s3_path=s3_path,
                        local_path=local_path)

        return_code = subprocess.call(cmd_string,
                    shell=True)
        print cmd_string, return_code        
        return return_code

    def check_file_in_checkpoint(self, cp_path):
        if cp_path[0:5] == "s3://":
            return self.check_s3_file_exists(cp_path)
        else:
            return os.path.exists(cp_path)

    def get_file_in_checkpoint(self, cp_path, local_path):
        if cp_path[0:5] == "s3://":
            return self.download_s3_file(cp_path, local_path)
        else:
            return subprocess.call("cp %s %s" % (cp_path, local_path), shell=True)
    
    def check_output_and_input_with_checkpoint(self, input_files, output_files):
        return self.check_output_with_checkpoint(input_files, output_files, check_inputs=True)

    def _flatten(self, foo):
        for x in foo:
            if hasattr(x, '__iter__') and not isinstance(x, str):
                for y in self._flatten(x):
                    yield y
            else:
                yield x

    def check_output_with_checkpoint(self, input_files, output_files, check_inputs=False):
        if type(input_files) not in [list, tuple]:
            input_files = [input_files]
        if type(output_files) not in [list, tuple]:
            output_files = [output_files]
        if check_inputs:
            output_files.extend(input_files)
        
        missing = []
        for local_path in self._flatten(output_files):
            if not os.path.exists(local_path):  # file does not exist locally
                if self.checkpoint_path: # 
                    cp_path = os.path.join(self.checkpoint_path, os.path.basename(local_path))
                    if not self.check_file_in_checkpoint(cp_path):  # and not in checkpoint either
                        print "%s Not found in checkpoint as %s" % (local_path, cp_path)
                        missing.append(local_path)
                    else:
                        rc = self.get_file_in_checkpoint(cp_path, local_path)  # does exist in checkpoint, attempt to download
                        if rc != 0:  # failure to download-- consider file as missing
                            missing.append(local_path)
                else:
                    missing.append(local_path) # no local file, no local path

        if len(missing) > 0:
            return True, "Missing files: %s" % ", ".join(missing)
        else:
            return False, "All output files exist"

    # def pipeline_run(self, stage=None):
    #     if stage is None and self.default_stage is not None:
    #         stage = self.default_stage
    #     else:
    #         raise PipelineException("Please specify a stage or set the default_stage!")
    #     ruffus.pipeline_run([stage], verbose=9, logger=self.logger, multiprocess=self.n_threads)
