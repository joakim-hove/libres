import os
import sys
from ecl.test import ExtendedTestCase, TestAreaContext
from ecl.util import BoolVector

from res.test import ErtTestContext
from res.enkf import EnkfFs, EnkfConfigNode, NodeId, EnkfNode
from res.enkf import EnKFMain, ErtRunContext
from res.enkf.enums import EnKFFSType
from res.enkf.enums import RealizationStateEnum

class EverestTest(ExtendedTestCase):
    def setUp(self):
        pass


    def test_run(self):
        ens_size = 2
        config = {"SIMULATION" : {"NUM_REALIZATIONS" : ens_size,
                                  "QUEUE_SYSTEM" :
                                  {"JOBNAME" : "job%d",
                                   "QUEUE_SYSTEM" : "LOCAL"}}}

        with ErtTestContext( "simple_everest", config_dict = config, store_area = True) as ctx:
            ert = ctx.getErt( )
            ens_config = ert.ensembleConfig( )
            sys.stderr.write("CWD: %s \n" % os.getcwd())


            # Add control nodes
            order_control = EnkfConfigNode.create_ext_param("WELL_ORDER", ["W1","W2","W3"])
            injection_control = EnkfConfigNode.create_ext_param("WELL_INJECTION", ["W1","W4"])
            ens_config.addNode( order_control )
            ens_config.addNode( injection_control )

            # Add result nodes
            order_result = EnkfConfigNode.create_gen_data( "ORDER", "order_%d" )
            injection_result = EnkfConfigNode.create_gen_data( "INJECTION", "injection_%d" )
            ens_config.addNode( order_result )
            ens_config.addNode( injection_result )


            order_node = EnkfNode( order_control ) ; order_node_ext = order_node.as_ext_param( )
            injection_node = EnkfNode( injection_control ) ; injection_node_ext = injection_node.as_ext_param( )

            fs_manager = ert.getEnkfFsManager( )
            sim_fs = fs_manager.getFileSystem("sim_fs")
            state_map = sim_fs.getStateMap( )
            batch_size = ens_size + 2
            for iens in range(batch_size):
                node_id = NodeId( 0, iens)

                order_node_ext["W1"] = iens
                order_node_ext["W2"] = iens * 10
                order_node_ext["W3"] = iens * 100
                order_node.save( sim_fs, node_id)

                injection_node_ext["W1"] = iens + 1
                injection_node_ext["W4"] = 3*(iens + 1)
                injection_node.save( sim_fs, node_id)
                state_map[iens] = RealizationStateEnum.STATE_INITIALIZED

            mask = BoolVector( default_value = True, initial_size = batch_size )
            model_config = ert.getModelConfig( )
            runpath_fmt = model_config.getRunpathFormat( )
            jobname_fmt = model_config.getJobnameFormat( )
            subst_list = ert.getDataKW( )
            itr = 0

            run_context = ErtRunContext.ensemble_experiment( sim_fs, mask, runpath_fmt, jobname_fmt, subst_list, itr)
            ert.getEnkfSimulationRunner().createRunPath( run_context )

            job_queue = ert.get_queue_config().create_job_queue()
            #num = ert.getEnkfSimulationRunner().runEnsembleExperiment(job_queue, run_context)
