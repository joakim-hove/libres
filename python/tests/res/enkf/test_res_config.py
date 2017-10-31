#  Copyright (C) 2017  Statoil ASA, Norway.
#
#  The file 'test_res_config.py' is part of ERT - Ensemble based Reservoir Tool.
#
#  ERT is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ERT is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
#  for more details.
import os, os.path
from copy import deepcopy
from datetime import date

from ecl.test import ExtendedTestCase, TestAreaContext
from ecl.util import CTime
from ecl.util.enums import RngAlgTypeEnum, MessageLevelEnum

from res.sched import HistorySourceEnum

from res.enkf import ResConfig, SiteConfig, AnalysisConfig
from res.test import ErtTestContext

# The res_config object should set the environment variable
# 'DATA_ROOT' to the root directory with the config
# file. Unfortunately the python methods to get environment variable,
# os.getenv() and os.environ[] do not reflect the:
#
#    setenv( "DATA_ROOT" , ...)
#
# call in the res_config C code. We therefor create a wrapper to the
# underlying libc getenv() function to be used for testing.

from cwrap import Prototype
from cwrap import load
clib = load( None )
clib_getenv = Prototype(clib, "char* getenv( char* )" , bind = False)


config_defines = {
        "<USER>"            : "TEST_USER",
        "<SCRATCH>"         : "scratch/ert",
        "<CASE_DIR>"        : "the_extensive_case",
        "<ECLIPSE_NAME>"    : "XYZ"
        }

config_data = {
        "RUNPATH"           : "<SCRATCH>/<USER>/<CASE_DIR>/realization-%d/iter-%d",
        "NUM_REALIZATIONS"  : 10,
        "MAX_RUNTIME"       : 23400,
        "MIN_REALIZATIONS"  : "50%",
        "MAX_SUBMIT"        : 13,
        "QUEUE_SYSTEM"      : "LSF",
        "UMASK"             : int("007", 8),
        "MAX_RUNNING"       : "100",
        "DATA_FILE"         : "../../eclipse/model/SNAKE_OIL.DATA",
        "START"             : date(2017, 1, 1),
        "SUMMARY"           : ["WOPR:PROD", "WOPT:PROD", "WWPR:PROD", "WWCT:PROD",
                               "WWPT:PROD", "WBHP:PROD", "WWIR:INJ", "WWIT:INJ",
                               "WBHP:INJ", "ROE:1"],
        "GEN_KW"            : ["SIGMA"],
        "ENSPATH"           : "../output/storage/<CASE_DIR>",
        "PLOT_PATH"         : "../output/results/plot/<CASE_DIR>",
        "UPDATE_LOG_PATH"   : "../output/update_log/<CASE_DIR>",
        "LOG_FILE"          : "../output/log/ert_<CASE_DIR>.log",
        "RUNPATH_FILE"      : "../output/run_path_file/.ert-runpath-list_<CASE_DIR>",
        "REFCASE"           : "../input/refcase/SNAKE_OIL_FIELD",
        "SIGMA"             : {
                                  "TEMPLATE"  : "../input/templates/sigma.tmpl",
                                  "RESULT"    : "coarse.sigma",
                                  "PARAMETER" : "../input/distributions/sigma.dist"
                              },
        "JOBNAME"           : "SNAKE_OIL_STRUCTURE_%d",
        "INSTALL_JOB"       : {
                                  "SNAKE_OIL_SIMULATOR" : {
                                      "CONFIG"     : "../../snake_oil/jobs/SNAKE_OIL_SIMULATOR",
                                      "STDOUT"     : "snake_oil.stdout",
                                      "STDERR"     : "snake_oil.stderr",
                                      "EXECUTABLE" : "snake_oil_simulator.py"
                                      },
                                  "SNAKE_OIL_NPV" : {
                                      "CONFIG"     : "../../snake_oil/jobs/SNAKE_OIL_NPV",
                                      "STDOUT"     : "snake_oil_npv.stdout",
                                      "STDERR"     : "snake_oil_npv.stderr",
                                      "EXECUTABLE" : "snake_oil_npv.py"
                                      },
                                  "SNAKE_OIL_DIFF" : {
                                      "CONFIG"     : "../../snake_oil/jobs/SNAKE_OIL_DIFF",
                                      "STDOUT"     : "snake_oil_diff.stdout",
                                      "STDERR"     : "snake_oil_diff.stderr",
                                      "EXECUTABLE" : "snake_oil_diff.py"
                                      }
                              },
        "FORWARD_MODEL"     : ["SNAKE_OIL_SIMULATOR", "SNAKE_OIL_NPV", "SNAKE_OIL_DIFF"],
        "HISTORY_SOURCE"    : HistorySourceEnum.REFCASE_HISTORY,
        "OBS_CONFIG"        : "../input/observations/obsfiles/observations.txt",
        "LOAD_WORKFLOW"     : {
                                  "MAGIC_PRINT" : "../bin/workflows/MAGIC_PRINT"
                              },
        "LOAD_WORKFLOW_JOB" : {
                                  "UBER_PRINT"  : "../bin/workflows/workflowjobs/bin/uber_print.py"
                              },
        "LOG_LEVEL"         : MessageLevelEnum.LOG_INFO,
        "RNG_ALG_TYPE"      : RngAlgTypeEnum.MZRAN,
        "STORE_SEED"        : "../input/rng/SEED",
        "LOAD_SEED"         : "../input/rng/SEED",
        "GRID"              : "../../eclipse/include/grid/CASE.EGRID",
        "RUN_TEMPLATE"      : {
                                  "seed_template" : {
                                      "TEMPLATE_FILE" : "../input/templates/seed_template.txt",
                                      "TARGET_FILE"   : "seed.txt"
                                      }
                              }
        }

def expand_config_data():
    '''
    Expands all strings in config_data according to config_defines. This is to
    enable one to copy the expected data directly from the configuration file
    without having to do manual expantion (that is what we have computers for
    anyway).
    '''
    for define_key in config_defines:
        for data_key in config_data:
            if isinstance(config_data[data_key], str):
                config_data[data_key] = config_data[data_key].replace(
                                                        define_key,
                                                        config_defines[define_key]
                                                        )

class ResConfigTest(ExtendedTestCase):

    def set_up_simple(self):
        self.case_directory = self.createTestPath("local/simple_config/")

    def set_up_snake_oil_structure(self):
        self.case_directory = self.createTestPath("local/snake_oil_structure")
        self.config_file = "snake_oil_structure/ert/model/user_config.ert"
        expand_config_data()

    def test_invalid_user_config(self):
        self.set_up_simple()

        with TestAreaContext("void land"):
            with self.assertRaises(IOError):
                ResConfig("this/is/not/a/file")

    def test_init(self):
        self.set_up_simple()

        with TestAreaContext("res_config_init_test") as work_area:
            cwd = os.getcwd()
            work_area.copy_directory(self.case_directory)

            config_file = "simple_config/minimum_config"
            res_config = ResConfig(user_config_file=config_file)

            self.assertEqual( res_config.model_config.data_root( ) , os.path.join( cwd , "simple_config"))
            self.assertEqual( clib_getenv("DATA_ROOT") , os.path.join( cwd , "simple_config"))

            # This fails with an not-understandable Python error:
            #-----------------------------------------------------------------
            # res_config.model_config.set_data_root( "NEW" )
            # self.assertEqual( res_config.model_config.data_root( ) , "NEW")
            # self.assertEqual( clib_getenv("DATA_ROOT") , "NEW")

            self.assertIsNotNone(res_config)
            self.assert_same_config_file(config_file, res_config.user_config_file, os.getcwd())

            self.assertIsNotNone(res_config.site_config)
            self.assertTrue(isinstance(res_config.site_config, SiteConfig))

            self.assertIsNotNone(res_config.analysis_config)
            self.assertTrue(isinstance(res_config.analysis_config, AnalysisConfig))

            self.assertEqual( res_config.config_path , os.path.join( cwd , "simple_config"))

            config_file = os.path.join( cwd, "simple_config/minimum_config")
            res_config = ResConfig(user_config_file=config_file)
            self.assertEqual( res_config.config_path , os.path.join( cwd , "simple_config"))

            os.chdir("simple_config")
            config_file = "minimum_config"
            res_config = ResConfig(user_config_file=config_file)
            self.assertEqual( res_config.config_path , os.path.join( cwd , "simple_config"))

            subst_config = res_config.subst_config
            for t in subst_config:
                print t
            self.assertEqual( subst_config["<CONFIG_PATH>"], os.path.join( cwd , "simple_config"))


    def assert_same_config_file(self, expected_filename, filename, prefix):
        prefix_path = lambda fn: fn if os.path.isabs(fn) else os.path.join(prefix, fn)
        canonical_path = lambda fn: os.path.realpath(os.path.abspath(prefix_path(fn)))

        self.assertEqual(
                    canonical_path(expected_filename),
                    canonical_path(filename)
                    )


    def assert_model_config(self, model_config, config_data, working_dir):
        self.assertEqual(
                config_data["RUNPATH"],
                model_config.getRunpathAsString()
                )

        self.assertEqual(
                config_data["ENSPATH"],
                model_config.getEnspath()
                )

        self.assertEqual(
                config_data["JOBNAME"],
                model_config.getJobnameFormat()
                )

        self.assertEqual(
                config_data["FORWARD_MODEL"],
                model_config.getForwardModel().joblist()
                )

        self.assertEqual(
                config_data["HISTORY_SOURCE"],
                model_config.get_history_source()
                )

        self.assertEqual(
                config_data["NUM_REALIZATIONS"],
                model_config.num_realizations
                )

        self.assert_same_config_file(
                config_data["OBS_CONFIG"],
                model_config.obs_config_file,
                working_dir
                )


    def assert_analysis_config(self, analysis_config, config_data):
        self.assertEqual(
                config_data["MAX_RUNTIME"],
                analysis_config.get_max_runtime()
                )

        self.assertEqual(
                config_data["UPDATE_LOG_PATH"],
                analysis_config.get_log_path()
                )

    def assert_site_config(self, site_config, config_data, working_dir):
        self.assertEqual(
                config_data["MAX_SUBMIT"],
                site_config.queue_config.max_submit
                )

        self.assertEqual(
                config_data["QUEUE_SYSTEM"],
                site_config.queue_config.queue_name
                )

        self.assertEqual(
                config_data["QUEUE_SYSTEM"],
                site_config.queue_config.driver.name
                )

        self.assertEqual(
                config_data["MAX_RUNNING"],
                site_config.queue_config.driver.get_option("MAX_RUNNING")
                )

        self.assertEqual(
                config_data["UMASK"],
                site_config.umask
                )

        job_list = site_config.get_installed_jobs()
        for job_name in config_data["INSTALL_JOB"]:
            self.assertTrue(job_name in job_list)

            exp_job_data = config_data["INSTALL_JOB"][job_name]

            self.assert_same_config_file(
                    exp_job_data["CONFIG"],
                    job_list[job_name].get_config_file(),
                    working_dir
                    )

            self.assertEqual(
                    exp_job_data["STDERR"],
                    job_list[job_name].get_stderr_file()
                    )

            self.assertEqual(
                    exp_job_data["STDOUT"],
                    job_list[job_name].get_stdout_file()
                    )

    def assert_ecl_config(self, ecl_config, config_data, working_dir):
        self.assert_same_config_file(
                config_data["DATA_FILE"],
                ecl_config.getDataFile(),
                working_dir
                )

        self.assertEqual(
                CTime(config_data["START"]),
                ecl_config.getStartDate()
                )


        for extension in ["SMSPEC", "UNSMRY"]:
            self.assert_same_config_file(
                    config_data["REFCASE"] + "." + extension,
                    ecl_config.getRefcaseName() + "." + extension,
                    working_dir
                    )

        self.assert_same_config_file(
                config_data["GRID"],
                ecl_config.get_gridfile(),
                working_dir
                )

    def assert_ensemble_config(self, ensemble_config, config_data, working_dir):
        self.assertEqual(
            set(config_data["SUMMARY"] + config_data["GEN_KW"]),
            set(ensemble_config.alloc_keylist())
            )

        loaded_template_file = ensemble_config["SIGMA"].getKeywordModelConfig().getTemplateFile()
        self.assert_same_config_file(
                config_data["SIGMA"]["TEMPLATE"],
                loaded_template_file,
                working_dir
                )

        loaded_parameter_file = ensemble_config["SIGMA"].getKeywordModelConfig().getParameterFile()
        self.assert_same_config_file(
                config_data["SIGMA"]["PARAMETER"],
                loaded_parameter_file,
                working_dir
                )

        self.assert_same_config_file(
                config_data["SIGMA"]["RESULT"],
                ensemble_config["SIGMA"]._get_enkf_outfile(),
                working_dir
                )


    def assert_plot_config(self, plot_config, config_data):
        self.assertEqual(
                config_data["PLOT_PATH"],
                plot_config.getPath()
                )

    def assert_hook_manager(self, hook_manager, config_data, working_dir):
        self.assert_same_config_file(
                    config_data["RUNPATH_FILE"],
                    hook_manager.getRunpathList().getExportFile(),
                    working_dir
                    )

    def assert_ert_workflow_list(self, ert_workflow_list, config_data, working_dir):
        self.assertEqual(
                len(config_data["LOAD_WORKFLOW"]),
                len(ert_workflow_list.getWorkflowNames())
                )

        for w_name in config_data["LOAD_WORKFLOW"]:
            self.assertTrue(w_name in ert_workflow_list)

            self.assert_same_config_file(
                    config_data["LOAD_WORKFLOW"][w_name],
                    ert_workflow_list[w_name].src_file,
                    working_dir
                    )

        for wj_name in config_data["LOAD_WORKFLOW_JOB"]:
            self.assertTrue(ert_workflow_list.hasJob(wj_name))
            job = ert_workflow_list.getJob(wj_name)

            self.assertEqual(wj_name, job.name())
            self.assert_same_config_file(
                    config_data["LOAD_WORKFLOW_JOB"][wj_name],
                    job.executable(),
                    working_dir
                    )

    def assert_rng_config(self, rng_config, config_data, working_dir):
        self.assertEqual(
                config_data["RNG_ALG_TYPE"],
                rng_config.alg_type
                )

        self.assert_same_config_file(
                config_data["STORE_SEED"],
                rng_config.store_filename,
                working_dir
                )

        self.assert_same_config_file(
                config_data["LOAD_SEED"],
                rng_config.load_filename,
                working_dir
                )


    def assert_ert_templates(self, ert_templates, config_data, working_dir):
        self.assertEqual(
                    config_data["RUN_TEMPLATE"].keys(),
                    ert_templates.getTemplateNames()
                    )

        for template_name in ert_templates.getTemplateNames():
            ert_template = ert_templates.get_template(template_name)
            config_template = config_data["RUN_TEMPLATE"][template_name]

            self.assert_same_config_file(
                    config_template["TEMPLATE_FILE"],
                    ert_template.get_template_file(),
                    working_dir
                    )

            self.assertEqual(
                    config_template["TARGET_FILE"],
                    ert_template.get_target_file(),
                    working_dir
                    )

    def assert_log_config(self, log_config, config_data, working_dir):
        self.assert_same_config_file(
                config_data["LOG_FILE"],
                log_config.log_file,
                working_dir
                )

        self.assertEqual(
                config_data["LOG_LEVEL"],
                log_config.log_level
                )


    def test_extensive_config(self):
        self.set_up_snake_oil_structure()

        with TestAreaContext("enkf_test_other_area") as work_area:
            work_area.copy_directory(self.case_directory)

            # Move to another directory
            run_dir = "i/ll/camp/here"
            os.makedirs(run_dir)
            os.chdir(run_dir)

            rel_config_file = "/".join(
                                   [".."] * len(run_dir.split("/")) +
                                   [self.config_file]
                                   )

            res_config = ResConfig(rel_config_file)

            work_dir = os.path.split(rel_config_file)[0]
            rel2workdir = lambda path : os.path.join(
                                             work_dir,
                                             path
                                             )

            self.assert_model_config(res_config.model_config, config_data, work_dir)
            self.assert_analysis_config(res_config.analysis_config, config_data)
            self.assert_site_config(res_config.site_config, config_data, work_dir)
            self.assert_ecl_config(res_config.ecl_config, config_data, work_dir)
            self.assert_ensemble_config(res_config.ensemble_config, config_data, work_dir)
            self.assert_plot_config(res_config.plot_config, config_data)
            self.assert_hook_manager(res_config.hook_manager, config_data, work_dir)
            self.assert_ert_workflow_list(res_config.ert_workflow_list, config_data, work_dir)
            self.assert_rng_config(res_config.rng_config, config_data, work_dir)
            self.assert_ert_templates(res_config.ert_templates, config_data, work_dir)
            self.assert_log_config(res_config.log_config, config_data, work_dir)


            # TODO: Not tested
            # - MIN_REALIZATIONS

    def test_missing_directory(self):
        config = {
            "INTERNALS" :
            {
                "CONFIG_DIRECTORY" : "does_not_exist",
            },
            "SIMULATION" :
            {
                "QUEUE_SYSTEM" :
                {
                    "JOBNAME" : "Job%d",
                },
                "RUNPATH"            : "/tmp/simulations/run%d",
                "NUM_REALIZATIONS"   : 1,
                "JOB_SCRIPT"         : "script.sh",
                "ENSPATH"            : "Ensemble"
            }
        }

        with self.assertRaises(IOError):
            ResConfig( config = config )
