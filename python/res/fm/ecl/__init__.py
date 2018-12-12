from .ecl_config import Ecl100Config
from .ecl_run import EclRun


def simulate(simulator, version, data_file, num_cpu = 1, check = True):
    run_cmd = "run_{}".format(simulator)
    if not check:
        run_cmd += "_nocheck"

    # The EclRun class should take a simulator instance as argument
    ecl_run = EclRun( [run_cmd, version, data_file, num_cpu])
    ecl_run.runEclipse( )



def run_ecl100(data_file, version = None, num_cpu = 1, check=True):
    config = Ecl100Config()
    sim = config.sim(version)
    ecl_run = EclRun(data_file, sim, num_cpu = num_cpu, check=check)
    ecl_run.runEclipse()
