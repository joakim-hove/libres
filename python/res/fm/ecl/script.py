from argparse import ArgumentParser
from .ecl_run import EclRun




def run(config, argv):
    parser = ArgumentParser()
    parser.add_argument("ecl_case")
    parser.add_argument("--version", dest="version", type=str)
    parser.add_argument("--num-cpu", dest="num_cpu", type=int, default=1)
    parser.add_argument("--ignore-errors", dest="ignore_errors", action="store_true")

    options = parser.parse_args(argv)
    sim = config.sim(version = options.version)
    run = EclRun(options.ecl_case, sim, num_cpu = options.num_cpu, check_status = not options.ignore_errors)
    run.runEclipse()
