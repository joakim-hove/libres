from ecl.util import ArgPack, CThreadPool

from res.job_queue import JobQueueManager

from res.enkf import ENKF_LIB
from res.enkf.ert_run_context import ErtRunContext
from res.enkf.run_arg import RunArg
from res.enkf.enums import EnkfRunType


class SimulationContext(object):
    def __init__(self, ert, init_fs, target_fs, size, itr , verbose=False):
        self._ert = ert
        """ :type: res.enkf.EnKFMain """
        self._size = size
        max_runtime = ert.analysisConfig().get_max_runtime()
               
        job_queue = ert.get_queue_config().create_job_queue()
        self._queue_manager = JobQueueManager(job_queue)
        self._queue_manager.startQueue(size, verbose=verbose)
        self._run_args = {}
        """ :type: dict[int, RunArg] """

        self._thread_pool = CThreadPool(8)
        self._thread_pool.addTaskFunction("submitJob", ENKF_LIB, "enkf_main_isubmit_job__")

        subst_list = self._ert.getDataKW( )
        path_fmt = self._ert.getModelConfig().getRunpathFormat()
        self._run_context = ErtRunContext( EnkfRunType.ENSEMBLE_EXPERIMENT, init_fs, target_fs, mask, path_fmt, subst_list, ier)


    def addSimulation(self, iens):
        if iens >= self._size:
            raise UserWarning("Realization number out of range: %d >= %d" % (iens, self._size))

        if iens in self._run_args:
            raise UserWarning("Realization number: '%d' already queued" % iens)

        run_arg = self._run_context.iensGet( iens )
        self._ert.createRunPath(run_arg)
        queue = self._queue_manager.get_job_queue()
        self._run_args[iens] = run_arg
        self._thread_pool.submitJob(ArgPack(self._ert, run_arg, queue))


    def isRunning(self):
        return self._queue_manager.isRunning()


    def getNumRunning(self):
        return self._queue_manager.getNumRunning()


    def getNumSuccess(self):
        return self._queue_manager.getNumSuccess()


    def getNumFailed(self):
        return self._queue_manager.getNumFailed()

    def getNumWaiting(self):
        return self._queue_manager.getNumWaiting()


    def didRealizationSucceed(self, iens):
        queue_index = self._run_args[iens].getQueueIndex()
        return self._queue_manager.didJobSucceed(queue_index)


    def didRealizationFail(self, iens):
        # For the purposes of this class, a failure should be anything (killed job, etc) that is not an explicit success.
        return not self.didRealizationSucceed(iens)


    def isRealizationQueued(self, iens):
        return iens in self._run_args


    def isRealizationFinished(self, iens):
        run_arg = self._run_args[iens]

        if run_arg.isSubmitted():
            queue_index = run_arg.getQueueIndex()
            return self._queue_manager.isJobComplete(queue_index)
        else:
            return False

    def __repr__(self):
        running = 'running' if self.isRunning() else 'not running'
        numRunn = self.getNumRunning()
        numSucc = self.getNumSuccess()
        numFail = self.getNumFailed()
        numWait = self.getNumWaiting()
        fmt = '%s, #running = %d, #success = %d, #failed = %d, #waiting = %d'
        fmt =  fmt % (running, numRunn, numSucc, numFail, numWait)
        return 'SimulationContext(%s)' % fmt
