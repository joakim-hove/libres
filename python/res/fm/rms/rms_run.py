from __future__ import print_function
import json
import os
import sys
import os.path
import time
import random

from contextlib import contextmanager
from .rms_config import RMSConfig


@contextmanager
def pushd(path):
    cwd0 = os.getcwd()
    os.chdir(path)

    yield

    os.chdir(cwd0)


def _remove_python_from_path(path):
    python_exec = 'python'
    path_elems = []
    for elem in path.split(os.pathsep):
        if os.path.isfile(os.path.join(elem, python_exec)):
            print("**Warning: the path: {} is removed from $PATH to avoid conflict with RMS internal python interprete".format(elem))
        else:
            path_elems.append(elem)

    return os.pathsep.join(path_elems)


class RMSRun(object):
    _single_seed_file = "RMS_SEED"
    _multi_seed_file = "random.seeds"
    _max_seed = 2146483648
    _seed_factor = 7907

    def __init__(self, iens, project, workflow, run_path="rms",target_file=None, export_path="rmsEXPORT", import_path="rmsIMPORT", version=None, readonly=True):
        if not os.path.isdir(project):
            raise OSError("The project:{} does not exist as a directory.".format(project))

        self.config = RMSConfig()
        self.project = os.path.abspath(project)
        self.workflow = workflow
        self.run_path = run_path
        self.version = version
        self.readonly = readonly
        self.import_path = import_path
        self.export_path = export_path
        if target_file is None:
            self.target_file = None
        else:
            if os.path.isabs(target_file):
                self.target_file = target_file
            else:
                self.target_file = os.path.join(os.getcwd(), target_file)

            if os.path.isfile(self.target_file):
                st = os.stat(self.target_file)
                self.target_file_mtime = st.mtime
            else:
                self.target_file_mtime = None


        self.init_seed(iens)

    def init_seed(self, iens):
        if "RMS_SEED" in os.environ:
            seed = int(os.getenv("RMS_SEED"))
            for x in range(iens):
                seed *= RMSRun._seed_factor
        else:
            single_seed_file = os.path.join( self.run_path , RMSRun._single_seed_file)
            multi_seed_file = os.path.join( self.run_path , RMSRun._multi_seed_file)

            if os.path.exists( single_seed_file ):
                # Using existing single seed file
                with open(single_seed_file) as fileH:
                    seed  = int( float(fileH.readline( )) )
            elif os.path.exists( multi_seed_file ):
                with open( multi_seed_file) as fileH:
                    seed_list = [ int(x) for x in fileH.readlines() ]
                seed = seed_list[iens + 1]
            else:
                random.seed( )
                seed = random.randint( 0 , RMSRun._max_seed )

        self.seed = seed % RMSRun._max_seed



    def run(self):
        child_pid = os.fork( )
        if child_pid == 0:
            self._run_child()
        else:
            self.wait( child_pid )


    def _run_child(self):
        if not os.path.exists( self.run_path ):
            os.makedirs( self.run_path )

        self_exe, _ = os.path.splitext(os.path.basename(sys.argv[0]))
        exec_env_file = "%s_exec_env.json" % self_exe
        exec_env = { 'PATH': os.environ['PATH'] }
        if os.path.isfile(exec_env_file):
            exec_env.update(json.load(open(exec_env_file)))

        # TODO: This is to fix a bug in RMS 2013, should be removed when this
        # version is no longer to be supported in Equinor. The
        # _remove_python_from_path function, and adding PATH to the exec_env by
        # default, should then also be removed..
        if self.version is not None and self.version.startswith('2013'):
            exec_env['PATH'] = _remove_python_from_path(exec_env['PATH'])

        with pushd(self.run_path):
            fileH = open("RMS_SEED_USED", "a+")
            fileH.write("%s ... %d\n" % (time.strftime("%d-%m-%Y %H:%M:%S" , time.localtime(time.time())) , self.seed))
            fileH.close()

            if not os.path.exists( self.export_path ):
                os.makedirs( self.export_path )

            if not os.path.exists( self.import_path ):
                os.makedirs( self.import_path )

            self.exec_rms(exec_env)



    def wait(self, child_pid):
        _, wait_status = os.waitpid(child_pid, 0)
        exit_status = os.WEXITSTATUS(wait_status)

        if exit_status != 0:
            raise Exception("The RMS run failed with exit status: {}".format(exit_status))

        if self.target_file is None:
            return

        if not os.path.isfile(self.target_file):
            raise Exception("The RMS run did not produce the expected  file: {}".format(self.target_file))

        if self.target_file_mtime is None:
            return

        st = os.stat(self.target_file)
        if st.mtime == self.target_file_mtime:
            raise Exception("The target file:{} is unmodified - interpreted as failure".format(self.target_file))


    def exec_rms(self, exec_env):
        args = [self.config.executable,
                "-project", self.project,
                "-seed", str(self.seed),
                "-nomesa",
                "-export_path", self.export_path,
                "-import_path", self.import_path,
                "-batch", self.workflow]

        if self.version:
            args += ["-v", self.version]

        if self.readonly:
            args += ["-readonly"]

        if self.config.threads:
            args += ["-threads", str(self.config.threads)]

        if exec_env:
            env = os.environ.copy()
            for key, value in exec_env.items():
                if v is None:
                    env.pop(key)
                    continue

                v = v.strip()
                if len(v) == 0:
                    env.pop(v)
                    continue

                env[key] = str(value)

            os.execve( self.config.executable , args, env )
        else:
            os.execv(self.config.executable, args )



