#ifndef ES_TESTDATA_HPP
#define ES_TESTDATA_HPP

#include <string>

#include <ert/util/bool_vector.hpp>

#include <ert/res_util/matrix.hpp>


namespace res {
class es_testdata {
public:
  std::string path;

  matrix_type * S;
  matrix_type * E;
  matrix_type * R;
  matrix_type * D;
  matrix_type * dObs;
  bool_vector_type * ens_mask;
  bool_vector_type * obs_mask;
  int obs_size;
  int ens_size;


  es_testdata(const char * path);
  ~es_testdata();
};

}

#endif
