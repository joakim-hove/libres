#include <stdexcept>
#include <string>

#include <ert/res_util/matrix.hpp>

#include "es_testdata.hpp"

namespace res {

namespace {

class pushd {
public:
  pushd(const std::string& path) {
    if (!util_is_directory(path.c_str()))
      throw std::runtime_error("The path: " + path + " does not exist - can not proceed");

    util_chdir(path.c_str());
    this->org_cwd = util_alloc_cwd();
  }


  ~pushd() {
    util_chdir(this->org_cwd);
    free(this->org_cwd);
  }

private:
  char * org_cwd;
};

matrix_type * alloc_load(const char * name, int rows, int columns) {
  if (!util_file_exists(name))
    return NULL;

  bool row_major = true;
  FILE * stream = util_fopen(name, "r");
  matrix_type * m = matrix_alloc(rows, columns);
  matrix_fscanf_data(m, row_major, stream);
  fclose(stream);

  return m;
}


std::pair<int, int> load_size() {
  FILE * stream = fopen("size", "r");
  int active_ens_size, active_obs_size;

  fscanf(stream, "%d", &active_ens_size);
  fscanf(stream, "%d", &active_obs_size);

  fclose(stream);

  return {active_ens_size, active_obs_size};
}


}

es_testdata::es_testdata(const char * path) :
  path(path),
  S(nullptr),
  E(nullptr),
  R(nullptr),
  D(nullptr),
  dObs(nullptr)
{
  pushd tmp_path(path);

  auto size_pair = load_size();
  this->active_ens_size = size_pair.first;
  this->active_obs_size = size_pair.second;

  this->S = alloc_load("S", this->active_obs_size, this->active_ens_size);
  this->E = alloc_load("E", this->active_obs_size, this->active_ens_size);
  this->R = alloc_load("R", this->active_obs_size, this->active_obs_size);
  this->D = alloc_load("D", this->active_obs_size, this->active_ens_size);
  this->dObs = alloc_load("dObs", this->active_obs_size, 2);
}

es_testdata::~es_testdata() {
  if (this->S)
    matrix_free(this->S);

  if (this->E)
    matrix_free(this->E);

  if (this->R)
    matrix_free(this->R);

  if (this->D)
    matrix_free(this->D);

  if (this->dObs)
    matrix_free(this->dObs);
}
}
