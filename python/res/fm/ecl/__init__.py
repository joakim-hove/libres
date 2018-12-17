from .ecl_config import Ecl100Config, FlowConfig
from .ecl_run import EclRun
from .script import run

def simulate(simulator, version, data_file, num_cpu = 1, check = True):
    if simulator == "ecl100":
        return run_ecl100(data_file, version=version, num_cpu=num_cpu, check=check)
    elif simulator == "flow":
        return run_flow(data_file, version=version, num_cpu=num_cpu, check=check)
    elif simulator == "ecl300":
        return run_ecl300(data_file, version=version, num_cpu=num_cpu, check=check)
    else:
        raise Exception("No such simulator: {}".format(simulator))


