#ifndef ES_TESTDATA_HPP
#define ES_TESTDATA_HPP

#include <string>

#include <ert/util/bool_vector.hpp>


namespace res {
class es_testdata {
public:
  std::string path;

  matrix_type * S;
  matrix_type * E;
  matrix_type * R;
  matrix_type * D;
  matrix_type * dObs;
  int active_obs_size;
  int active_ens_size;


  es_testdata(const char * path);
  ~es_testdata();
};

}

#endif
