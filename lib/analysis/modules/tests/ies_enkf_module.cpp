#include <ert/util/test_util.hpp>
#include <ert/util/rng.h>

#include <ert/analysis/analysis_module.hpp>
#include <ert/analysis/std_enkf.hpp>

#include "es_testdata.hpp"

void test_steplength1(const char * module_lib, const char * path_testdata) {
  res::es_testdata testdata(path_testdata);
  rng_type * rng = rng_alloc( MZRAN, INIT_DEFAULT );

  analysis_module_type * std_module = analysis_module_alloc_internal("STD_ENKF");
  analysis_module_type * ies_module = analysis_module_alloc_external(module_lib);

  test_assert_true( analysis_module_set_var(std_module, ENKF_TRUNCATION_KEY_, "0.95") );
  test_assert_true( analysis_module_set_var(ies_module, ENKF_TRUNCATION_KEY_, "0.95") );

  analysis_module_init_update(std_module,
                              testdata.ens_mask,
                              testdata.obs_mask,
                              testdata.S,
                              testdata.R,
                              testdata.dObs,
                              testdata.E,
                              testdata.D,
                              rng);

  rng_free( rng );
  analysis_module_free(std_module);
  analysis_module_free(ies_module);
}


void test_load(const char * module_lib) {
  analysis_module_type * module = analysis_module_alloc_external( module_lib );
  test_assert_not_NULL( module );
  analysis_module_free( module );
}


int main(int argc, char ** argv) {
  const char * module_lib = argv[1];
  const char * path_testdata = argv[2];

  test_load(module_lib);
  test_steplength1(module_lib, path_testdata);
}
