/*
   Copyright (C) 2020  Equinor ASA, Norway.

   The file 'slurm_driver.cpp' is part of ERT - Ensemble based Reservoir Tool.

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
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <dlfcn.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/types.h>

#include <cstddef>
#include <ctime>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <ert/util/util.hpp>
#include <ert/util/stringlist.hpp>
#include <ert/res_util/log.hpp>
#include <ert/res_util/res_log.hpp>
#include <ert/res_util/res_env.hpp>

#include <ert/job_queue/queue_driver.hpp>
#include <ert/job_queue/slurm_driver.hpp>


struct SlurmJob {
  SlurmJob(int job_id) :
    job_id(job_id),
    string_id(std::to_string(job_id))
  {
  }

  int job_id;
  std::string string_id;
};



class SlurmStatus {
public:

  void update(int job_id, job_status_type status) {
    pthread_rwlock_wrlock(&this->lock);
    this->jobs[job_id] = status;
    pthread_rwlock_unlock(&this->lock);
  }


  void new_job(int job_id) {
    this->update(job_id, JOB_QUEUE_PENDING);
  }

  static bool active_status(job_status_type status) {
    if (status == JOB_QUEUE_RUNNING)
      return true;

    if (status == JOB_QUEUE_PENDING)
      return true;

    return false;
  }


  std::vector<int> squeue_update(const std::unordered_map<int, job_status_type>& squeue_jobs) {
    std::vector<int> active_jobs;

    pthread_rwlock_wrlock(&this->lock);
    for (auto& job_pair : this->jobs) {
      auto job_id = job_pair.first;
      auto job_status = job_pair.second;

      auto squeue_pair = squeue_jobs.find( job_id );
      if (squeue_pair == squeue_jobs.end()) {
        if (this->active_status(job_status))
          active_jobs.push_back(job_id);
      } else
        this->jobs[job_id] = squeue_pair->second;
    }
    pthread_rwlock_unlock(&this->lock);

    return active_jobs;
  }


  job_status_type get(int job_id) const {

    pthread_rwlock_rdlock(&this->lock);
    auto status = this->jobs.at(job_id);
    pthread_rwlock_unlock(&this->lock);

    return status;
  }

private:
  std::unordered_map<int, job_status_type> jobs;
  mutable pthread_rwlock_t lock = PTHREAD_RWLOCK_INITIALIZER;
};


#define SLURM_DRIVER_TYPE_ID  70555081
#define DEFAULT_SBATCH_CMD    "sbatch"
#define DEFAULT_SCANCEL_CMD   "scancel"
#define DEFAULT_SQUEUE_CMD    "sqeueue"
#define DEFAULT_SCONTROL_CMD  "scontrol"

#define SLURM_PENDING_STATUS    "PENDING"
#define SLURM_COMPLETED_STATUS  "COMPLETED"
#define SLURM_RUNNING_STATUS    "RUNNING"
#define SLURM_FAILED_STATUS     "FAILED"
#define SLURM_CANCELED_STATUS   "CANCELLED"
#define SLURM_COMPLETING_STATUS "COMPLETING"


struct slurm_driver_struct {
  UTIL_TYPE_ID_DECLARATION;

  std::string sbatch_cmd;
  std::string scancel_cmd;
  std::string squeue_cmd;
  std::string scontrol_cmd;
  std::string partition;
  std::string username;
  mutable SlurmStatus status;
  mutable std::time_t status_timestamp;
  double status_timeout = 10;
};


static std::string load_file(const char * fname) {
    char * buffer = util_fread_alloc_file_content(fname, nullptr);
    std::string s = buffer;
    free(buffer);
    return s;
}


static std::string load_stdout(const char * cmd, int argc, const char ** argv) {
    std::string fname = std::string(cmd) + "-stdout";
    char * stdout = (char*) util_alloc_tmp_file("/tmp", fname.c_str(), true);

    util_spawn_blocking(cmd, argc, argv, stdout, nullptr);
    auto file_content = load_file(stdout);

    util_unlink_existing(stdout);
    free(stdout);
    return file_content;
}


static std::string load_stdout(const char * cmd, const std::vector<std::string>& args) {
    const char ** argv = static_cast<const char **>(util_calloc( args.size(), sizeof * argv ));
    for (std::size_t i=0; i < args.size(); i++)
        argv[i] = args[i].c_str();

    auto file_content = load_stdout(cmd, args.size(), argv);
    free(argv);
    return file_content;
}



UTIL_SAFE_CAST_FUNCTION( slurm_driver , SLURM_DRIVER_TYPE_ID)
static UTIL_SAFE_CAST_FUNCTION_CONST( slurm_driver , SLURM_DRIVER_TYPE_ID)

void * slurm_driver_alloc() {
  slurm_driver_type * driver = new slurm_driver_type();
  UTIL_TYPE_ID_INIT(driver, SLURM_DRIVER_TYPE_ID);
  driver->sbatch_cmd = DEFAULT_SBATCH_CMD;
  driver->scancel_cmd = DEFAULT_SCANCEL_CMD;
  driver->squeue_cmd = DEFAULT_SQUEUE_CMD;
  driver->scontrol_cmd = DEFAULT_SCONTROL_CMD;


  auto pwname = getpwuid( geteuid() );
  driver->username = pwname->pw_name;
  return driver;
}

void slurm_driver_free(slurm_driver_type * driver) {
  delete driver;
}

void slurm_driver_free__(void * __driver ) {
  slurm_driver_type * driver = slurm_driver_safe_cast( __driver );
  slurm_driver_free( driver );
}


const void * slurm_driver_get_option( const void * __driver, const char * option_key) {
  const slurm_driver_type * driver = slurm_driver_safe_cast_const( __driver );
  if (strcmp(option_key, SLURM_SBATCH_OPTION) == 0)
    return driver->sbatch_cmd.c_str();

  if (strcmp(option_key, SLURM_SCANCEL_OPTION) == 0)
    return driver->scancel_cmd.c_str();

  if (strcmp(option_key, SLURM_SCONTROL_OPTION) == 0)
    return driver->scontrol_cmd.c_str();

  if (strcmp(option_key, SLURM_SQUEUE_OPTION) == 0)
    return driver->squeue_cmd.c_str();

  if (strcmp(option_key, SLURM_PARTITION_OPTION) == 0)
    return driver->partition.c_str();

  return nullptr;
}


bool slurm_driver_set_option( void * __driver, const char * option_key, const void * value) {
  slurm_driver_type * driver = slurm_driver_safe_cast( __driver );
  if (strcmp(option_key, SLURM_SBATCH_OPTION) == 0) {
    driver->sbatch_cmd = static_cast<const char*>(value);
    return true;
  }

  if (strcmp(option_key, SLURM_SCANCEL_OPTION) == 0) {
    driver->scancel_cmd = static_cast<const char*>(value);
    return true;
  }

  if (strcmp(option_key, SLURM_SQUEUE_OPTION) == 0) {
    driver->squeue_cmd = static_cast<const char*>(value);
    return true;
  }

  if (strcmp(option_key, SLURM_SCONTROL_OPTION) == 0) {
    driver->scontrol_cmd = static_cast<const char*>(value);
    return true;
  }

  if (strcmp(option_key, SLURM_PARTITION_OPTION) == 0) {
    driver->partition = static_cast<const char*>(value);
    return true;
  }

  return false;
}


void slurm_driver_init_option_list(stringlist_type * option_list) {
  stringlist_append_copy(option_list, SLURM_PARTITION_OPTION);
  stringlist_append_copy(option_list, SLURM_SBATCH_OPTION);
  stringlist_append_copy(option_list, SLURM_SCONTROL_OPTION);
  stringlist_append_copy(option_list, SLURM_SQUEUE_OPTION);
  stringlist_append_copy(option_list, SLURM_SCANCEL_OPTION);
}


static std::string make_submit_script(const char * cmd, int num_cpu, int argc, const char ** argv) {
  char * submit        = (char*) util_alloc_tmp_file("/tmp" , "slurm-submit" , true);

  FILE * submit_stream = util_fopen(submit, "w");
  fprintf(submit_stream, "#!/bin/sh\n");

  fprintf(submit_stream, "%s", cmd);  // Without srun?
  for (int iarg=0; iarg < argc; iarg++)
    fprintf(submit_stream, " %s", argv[iarg]);
  fprintf(submit_stream, "\n");

  fclose(submit_stream);
  chmod(submit, S_IRWXU + S_IRGRP + S_IROTH);

  std::string submit_script = submit;
  free( submit );
  return submit_script;
}


/*
 The slurm jobs are submitted by first creating a submit script, which is a
 small shell which contains the command to run along with possible slurm
 options, and then this script is submitted with the 'sbatch' command.
*/

void * slurm_driver_submit_job( void * __driver, const char * cmd, int num_cpu, const char * run_path, const char * job_name, int argc, const char ** argv) {
  slurm_driver_type * driver = slurm_driver_safe_cast( __driver );

  auto submit_script = make_submit_script( cmd, num_cpu, argc, argv);
  std::vector<std::string> sbatch_argv = {"--workdir=" + std::string(run_path), "--job-name=" + std::string(job_name), "--parsable"};
  if (!driver->partition.empty())
    sbatch_argv.push_back( "--partition=" + driver->partition );
  sbatch_argv.push_back( submit_script );


  auto file_content = load_stdout(driver->sbatch_cmd.c_str(), sbatch_argv);
  util_unlink_existing( submit_script.c_str() );

  int job_id = 0;
  util_sscanf_int(file_content.c_str(), &job_id);
  if (job_id == 0)
    return nullptr;

  driver->status.new_job(job_id);
  return new SlurmJob(job_id);
}


static job_status_type slurm_driver_translate_status(const std::string& status_string, const std::string& string_id) {
  if (status_string == SLURM_PENDING_STATUS)
    return JOB_QUEUE_PENDING;

  if (status_string == SLURM_COMPLETED_STATUS)
    return JOB_QUEUE_DONE;

  if (status_string == SLURM_COMPLETING_STATUS)
    return JOB_QUEUE_RUNNING;

  if (status_string == SLURM_RUNNING_STATUS)
    return JOB_QUEUE_RUNNING;

  if (status_string == SLURM_FAILED_STATUS)
    return JOB_QUEUE_EXIT;

  if (status_string == SLURM_CANCELED_STATUS)
    return JOB_QUEUE_IS_KILLED;

  res_log_fwarning("The job status: \'%s\' for job:%s is not recognized", status_string.c_str(), string_id.c_str());
  return JOB_QUEUE_UNKNOWN;
}



static std::unordered_map<std::string, std::string> load_scontrol(const slurm_driver_type * driver, const std::string& string_id) {
  auto file_content = load_stdout(driver->scontrol_cmd.c_str(), {"show", "jobid", string_id});

  std::unordered_map<std::string, std::string> options;
  std::size_t offset = 0;
  while(true) {
    auto new_offset = file_content.find_first_of("\n ", offset);
    if (new_offset == std::string::npos)
      break;

    std::string key_value = file_content.substr(offset, new_offset - offset);
    auto split_pos = key_value.find('=');
    if (split_pos != std::string::npos) {
      std::string key = key_value.substr(0, split_pos);
      std::string value = key_value.substr(split_pos + 1);

      options.insert({key,value});
    }
    offset = file_content.find_first_not_of("\n ", new_offset);
  }
  return options;
}


static job_status_type slurm_driver_get_job_status_scontrol(const slurm_driver_type * driver, const std::string& string_id) {
  auto values = load_scontrol(driver, string_id);
  const auto status_iter = values.find("JobState");

  /*
    When a job has finished running it quite quickly - the order of minutes -
    falls out of the slurm database, and the scontrol command will not give any
    output. In this situation we guess that the job has completed succesfully
    and return status JOB_QUEUE_DONE. If the job has actually not succeded this
    should be picked up the libres post run checking.
  */
  if (status_iter == values.end()) {
    res_log_fwarning("The command \'scontrol show jobid %s\' gave no output for job:%s - assuming it is COMPLETED", string_id.c_str(), string_id.c_str());
    return JOB_QUEUE_DONE;
  }

  const auto& status_string = status_iter->second;
  auto status = slurm_driver_translate_status(status_string, string_id);

  if (status == JOB_QUEUE_UNKNOWN) {
    res_log_fwarning("The job status: \'%s\' for job:%s is not recognized - assuming it is RUNNING", status_string.c_str(), string_id.c_str());
    status = JOB_QUEUE_RUNNING;
  }

  return status;
}

static job_status_type slurm_driver_get_job_status_scontrol(const slurm_driver_type * driver, int job_id) {
  return slurm_driver_get_job_status_scontrol(driver, std::to_string(job_id));
}

static void slurm_driver_update_status_cache(const slurm_driver_type * driver) {
  driver->status_timestamp = time( nullptr );
  const std::string space = " \n";
  auto squeue_output = load_stdout(driver->squeue_cmd.c_str(), {"-h", "--user=" + driver->username, "--format=%i %T"});
  auto offset = squeue_output.find_first_not_of(space);

  std::unordered_map<int, job_status_type> squeue_jobs;
  while (offset != std::string::npos) {
    auto id_end = squeue_output.find_first_of(space, offset);
    auto job_string = squeue_output.substr( offset, id_end - offset);
    int job_id = std::stoi(squeue_output.substr( offset, id_end - offset));

    auto status_start = squeue_output.find_first_not_of(space, id_end + 1);
    auto status_end = squeue_output.find_first_of(space, status_start);
    auto status = slurm_driver_translate_status(squeue_output.substr( status_start, status_end - status_start ), std::to_string(job_id));

    squeue_jobs.insert( {job_id, status} );
    offset = squeue_output.find_first_not_of(space, status_end);
  }

  const auto& active_jobs = driver->status.squeue_update(squeue_jobs);
  for (const auto& job_id : active_jobs) {
    auto status = slurm_driver_get_job_status_scontrol(driver, job_id);
    driver->status.update(job_id, status);
  }
}



job_status_type slurm_driver_get_job_status(void * __driver , void * __job) {
  slurm_driver_type * driver = slurm_driver_safe_cast( __driver );
  const auto * job = static_cast<const SlurmJob*>(__job);
  auto update_cache = difftime(time(nullptr), driver->status_timestamp) > driver->status_timeout;
  if (update_cache)
    slurm_driver_update_status_cache(driver);

  return driver->status.get( job->job_id );
}


void slurm_driver_kill_job(void * __driver , void * __job ) {
  slurm_driver_type * driver = slurm_driver_safe_cast( __driver );
  const auto * job = static_cast<const SlurmJob*>(__job);
  const char ** argv = static_cast<const char **>(util_calloc( 1, sizeof * argv ));

  argv[0] = job->string_id.c_str();
  util_spawn_blocking(driver->scancel_cmd.c_str(), 1, argv, nullptr, nullptr);
  free(argv);
}

void slurm_driver_free_job(void * __job) {
  SlurmJob * job = static_cast<SlurmJob *>(__job);
  delete job;
}
