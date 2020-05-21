/*
   Copyright (C) 2020  Equinor ASA, Norway.

   The file 'job_slurm_runtest.cpp' is part of ERT - Ensemble based Reservoir Tool.

   ERT is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ERT is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.

   See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
   for more details.
*/

#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>

#include <vector>

#include <ert/util/test_util.hpp>
#include <ert/util/util.hpp>
#include <ert/util/test_work_area.hpp>

#include <ert/job_queue/queue_driver.hpp>

void make_sleep_job(const char * fname, int sleep_time) {
  FILE * stream = util_fopen(fname, "w");
  fprintf(stream, "sleep %d \n", sleep_time);
  fclose(stream);

  mode_t fmode = S_IRWXU;
  chmod( fname, fmode);
}


void make_failed_job(const char * fname, int sleep_time) {
  FILE * stream = util_fopen(fname, "w");
  fprintf(stream, "sleep %d \n", sleep_time);
  fprintf(stream, "exit 1\n", sleep_time);
  fclose(stream);

  mode_t fmode = S_IRWXU;
  chmod( fname, fmode);
}



void run() {
  ecl::util::TestArea ta("slurm_submit", true);
  queue_driver_type * driver = queue_driver_alloc_slurm();
  std::vector<void *> jobs;
  const char * long_cmd = util_alloc_abs_path("long_run.sh");
  const char * ok_cmd = util_alloc_abs_path("ok_run.sh");
  const char * fail_cmd = util_alloc_abs_path("failed_run.sh");
  int num_jobs = 6;

  make_sleep_job(long_cmd, 10);
  make_sleep_job(ok_cmd, 1);
  make_failed_job(fail_cmd, 1);

  for (int i = 0; i < num_jobs; i++) {
    std::string run_path = ta.test_cwd() + "/" + std::to_string(i);
    std::string job_name = "job" + std::to_string(i);
    util_make_path(run_path.c_str());
    if (i == 0)
      jobs.push_back( queue_driver_submit_job(driver, long_cmd, 1, run_path.c_str(), job_name.c_str(), 0, nullptr));
    else if (i == num_jobs - 1)
      jobs.push_back( queue_driver_submit_job(driver, fail_cmd, 1, run_path.c_str(), job_name.c_str(), 0, nullptr));
    else
      jobs.push_back( queue_driver_submit_job(driver, ok_cmd, 1, run_path.c_str(), job_name.c_str(), 0, nullptr));
  }


  while (true) {
    int active_count = 0;
    for (auto * job_ptr : jobs) {
      auto status = queue_driver_get_status(driver, job_ptr);
      if (status == JOB_QUEUE_RUNNING || status == JOB_QUEUE_PENDING || status == JOB_QUEUE_WAITING)
        active_count += 1;
    }

    if (active_count == 0)
      break;

    auto * long_job = jobs[0];
    auto long_status = queue_driver_get_status(driver, long_job);
    if (long_status != JOB_QUEUE_IS_KILLED)
      queue_driver_kill_job( driver, long_job );

    usleep(100000);
  }


  for (int i = 0; i < num_jobs; i++) {
    auto * job_ptr = jobs[i];
    if (i == 0)
      test_assert_int_equal( queue_driver_get_status(driver, job_ptr), JOB_QUEUE_IS_KILLED );
    else if (i == num_jobs - 1)
      test_assert_int_equal( queue_driver_get_status(driver, job_ptr), JOB_QUEUE_EXIT );
    else
      test_assert_int_equal( queue_driver_get_status(driver, job_ptr), JOB_QUEUE_DONE );
  }

  for (auto * job_ptr : jobs)
    queue_driver_free_job(driver, job_ptr);

  queue_driver_free(driver);
}


int main( int argc , char ** argv) {
  run();
  exit(0);
}
