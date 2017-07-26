from itertools import product
from random import random

import sys
import time

from res.enkf.enums import ErtImplType
from res.enkf.export.custom_kw_collector import CustomKWCollector
from res.server import ErtRPCClient
from ecl.test import ExtendedTestCase
from tests.res.server import RPCServiceContext
from res.enkf import EnkfFsManager, EnkfFs, EnkfVarType


def initializeCase(ert, name, size):
    """
    @type ert: res.enkf.enkf_main.EnKFMain
    @type name: str
    @type size: int
    @rtype:
    """
    current_fs = ert.getEnkfFsManager().getCurrentFileSystem()
    fs = ert.getEnkfFsManager().getFileSystem(name)
    ert.getEnkfFsManager().switchFileSystem(fs)
    parameters = ert.ensembleConfig().getKeylistFromVarType(EnkfVarType.PARAMETER)
    ert.getEnkfFsManager().initializeFromScratch(parameters, 0, size - 1)

    ert.getEnkfFsManager().switchFileSystem(current_fs)
    return fs



class RPCStorageTest(ExtendedTestCase):
    def test_rpc_storage(self):
        config = self.createTestPath("local/snake_oil_no_data/snake_oil.ert")

        with RPCServiceContext("ert/server/rpc_storage/storage", config) as server:
            client = ErtRPCClient("localhost", server.port)

            group_name = "Test"
            storage_definition = {
                "PI": float,
                "DakotaVersion": str,
                "Gradient": float,
                "GradientDirection": str
            }
            client.prototypeStorage(group_name, storage_definition)

            with self.assertRaises(UserWarning):
                client.prototypeStorage(group_name, {"value": float})

            with self.assertRaises(TypeError):
                client.prototypeStorage(group_name, {"value": bool})


            ensemble_config = server.ert.ensembleConfig()

            self.assertIn(group_name, ensemble_config.getKeylistFromImplType(ErtImplType.CUSTOM_KW))

            custom_kw_config = ensemble_config.getNode(group_name).getCustomKeywordModelConfig()

            self.assertIn("PI", custom_kw_config)
            self.assertTrue(custom_kw_config.keyIsDouble("PI"))

            self.assertIn("DakotaVersion", custom_kw_config)
            self.assertFalse(custom_kw_config.keyIsDouble("DakotaVersion"))

            self.assertIn("Gradient", custom_kw_config)
            self.assertTrue(custom_kw_config.keyIsDouble("Gradient"))

            self.assertIn("GradientDirection", custom_kw_config)
            self.assertFalse(custom_kw_config.keyIsDouble("GradientDirection"))


            simulation_count = 10
            initializeCase(server.ert, "default", simulation_count)

            client.storeGlobalData("default", group_name, "PI", 3.1415)
            client.storeGlobalData("default", group_name, "DakotaVersion", "DAKOTA 6.2.0")

            gradients = [random() * 20.0 - 10.0 for _ in range(simulation_count)]
            gradient_directions = [("POSITIVE" if gradient >= 0.0 else "NEGATIVE") for gradient in gradients]
            for sim_id in range(simulation_count):
                gradient = gradients[sim_id]
                gradient_direction = gradient_directions[sim_id]
                client.storeSimulationData("default", group_name, "Gradient", gradient, sim_id)
                client.storeSimulationData("default", group_name, "GradientDirection", gradient_direction, sim_id)


            data = CustomKWCollector.loadAllCustomKWData(server.ert, "default")
            for sim_id in range(simulation_count):
                self.assertEqual(data["Test:PI"][sim_id], 3.1415)
                self.assertEqual(data["Test:DakotaVersion"][sim_id], "DAKOTA 6.2.0")
                self.assertEqual(data["Test:Gradient"][sim_id], gradients[sim_id])
                self.assertEqual(data["Test:GradientDirection"][sim_id], gradient_directions[sim_id])

                
    def test_rpc_storage_with_simulation(self):
        config = self.createTestPath("local/snake_oil_no_data/snake_oil.ert")

        with RPCServiceContext("ert/server/rpc_storage/simulation_and_storage", config, store_area=True) as server:
            client = ErtRPCClient("localhost", server.port)
            realization_count = 2

            client.prototypeStorage("SNAKEX", {"SNAKE_ID": float, "SNAKE_RUN": str})

            case_names = ["test_run_0", "test_run_1"]

            for batch_id, case_name in enumerate(case_names):
                self.runSimulation(client, realization_count, batch_id, case_name)

            for (batch_id, case_name), iens in product(enumerate(case_names), range(realization_count)):
                result = client.getCustomKWResult(case_name, iens, "SNAKEX")
                self.assertEqual(result["SNAKE_RUN"], "batch_%d" % batch_id)
                snake_id = realization_count * batch_id + iens
                self.assertEqual(result["SNAKE_ID"], snake_id)



    def runSimulation(self, client, realization_count, batch_id, case_name):
        client.startSimulationBatch("default", case_name, 2)
        kw = {"SNAKE_OIL_PARAM": [0.50, 6, 1.750, 0.250, 0.990, 2, 1.770, 0.330, 0.550, 0.770]} # identical runs

        for iens in range(realization_count):
            client.addSimulation(0,0,iens,  kw)

        while client.isRunning():
            time.sleep(0.2)

        for iens in range(realization_count):
            self.assertTrue(client.didRealizationSucceed(iens))


        client.storeGlobalData(case_name, "SNAKEX", "SNAKE_RUN", "batch_%d" % batch_id)

        for iens in range(realization_count):
            client.storeSimulationData(case_name, "SNAKEX", "SNAKE_ID", realization_count * batch_id + iens, iens)







