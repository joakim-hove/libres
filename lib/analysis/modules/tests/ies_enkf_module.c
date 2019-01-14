#include <ert/util/test_util.hpp>

#include <ert/analysis/analysis_module.hpp>
#include <ert/analysis/std_enkf.hpp>

void test_steplength1(const char * module_lib) {
  analysis_module_type * std_module = analysis_module_alloc_internal("STD_ENKF");
  analysis_module_type * ies_module = analysis_module_alloc_external(module_lib);

  test_assert_true( analysis_module_set_var(std_module, ENKF_TRUNCATION_KEY_, "0.95") );
  test_assert_true( analysis_module_set_var(ies_module, ENKF_TRUNCATION_KEY_, "0.95") );


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

  test_load(module_lib);
  test_steplength1(module_lib);
}
