#include <ert/util/test_util.hpp>

#include <ert/analysis/analysis_module.hpp>



int main(int argc, char ** argv) {
  const char * module_lib = argv[1];
  analysis_module_type * module = analysis_module_alloc_external( module_lib );
  test_assert_not_NULL( module );
  analysis_module_free( module );
}
