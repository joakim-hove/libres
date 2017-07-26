/*
   Copyright (C) 2017  Statoil ASA, Norway.

   The file 'queue_config.c' is part of ERT - Ensemble based Reservoir Tool.

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

#include <ert/util/util.h>

#include <ert/job_queue/job_queue.h>
#include <ert/job_queue/job_node.h>
#include <ert/job_queue/ext_job.h>
#include <ert/job_queue/ext_joblist.h>
#include <ert/job_queue/lsf_driver.h>
#include <ert/job_queue/rsh_driver.h>
#include <ert/job_queue/local_driver.h>
#include <ert/job_queue/queue_driver.h>

#include <ert/config/config_parser.h>
#include <ert/config/config_content_item.h>
#include <ert/config/config_content_node.h>
#include <ert/config/config_schema_item.h>

#include <ert/enkf/queue_config.h>
#include <ert/enkf/config_keys.h>
#include <ert/enkf/enkf_defaults.h>
#include <ert/enkf/site_config.h>
#include <ert/enkf/model_config.h>

struct queue_config_struct {
    job_driver_type driver_type;
    char * job_script;
    hash_type * queue_drivers;
    bool user_mode;
    int max_submit;
};

static void queue_config_add_queue_driver(queue_config_type * queue_config, const char * driver_name, queue_driver_type * driver);
static bool queue_config_init(queue_config_type * queue_config, const config_content_type * config_content);

static queue_config_type * queue_config_alloc_empty() {
    queue_config_type * queue_config = util_malloc(sizeof * queue_config);
    queue_config->queue_drivers = hash_alloc();
    queue_config->job_script = NULL;
    queue_config->driver_type = NULL_DRIVER;
    queue_config->user_mode = false;
    queue_config->max_submit = 2; // Default value

    return queue_config;
}

static queue_config_type * queue_config_alloc() {
  queue_config_type * queue_config = queue_config_alloc_empty();

  config_parser_type * config = config_alloc();
  config_content_type * content = site_config_alloc_content(config);

  queue_config_init(queue_config, content);

  config_free(config);
  config_content_free(content);

  return queue_config;
}

queue_config_type * queue_config_alloc_load(const char * user_config_file) {
  queue_config_type * queue_config = queue_config_alloc();

  if(user_config_file) {
    queue_config->user_mode = true;

    config_parser_type * config = config_alloc();
    config_content_type * content = model_config_alloc_content(
                                                    user_config_file,
                                                    config
                                                    );

    queue_config_init(queue_config, content);

    config_free(config);
    config_content_free(content);
  }

  return queue_config;
}

queue_config_type * queue_config_alloc_local_copy( queue_config_type * queue_config) {
    queue_config_type * queue_config_copy = util_malloc(sizeof * queue_config_copy);
    queue_config_copy->queue_drivers = hash_alloc();
    queue_config_add_queue_driver(queue_config_copy, LOCAL_DRIVER_NAME, queue_driver_alloc_local());
    queue_config_copy->user_mode = queue_config->user_mode;

    if (queue_config_has_job_script(queue_config)) {
        queue_config_copy->job_script = util_alloc_string_copy( queue_config->job_script );
    }
    else
        queue_config_copy->job_script = NULL;

    if (queue_config->driver_type == NULL_DRIVER)
        queue_config_copy->driver_type = NULL_DRIVER;
    else
        queue_config_copy->driver_type = LOCAL_DRIVER;
    
    queue_config_copy->max_submit = queue_config->max_submit;

    return queue_config_copy;
}

job_queue_type * queue_config_alloc_job_queue(const queue_config_type * queue_config) {
  job_queue_type * job_queue = job_queue_alloc(DEFAULT_MAX_SUBMIT, "OK", "STATUS", "ERROR");
  const char * driver_name = queue_config_get_queue_name(queue_config);
  if (driver_name != NULL)
    {
        queue_driver_type * driver = queue_config_get_queue_driver(queue_config, driver_name);
        job_queue_set_driver(job_queue, driver);
    }

  job_queue_set_max_submit(job_queue, queue_config->max_submit);

  return job_queue;
}

void queue_config_free(queue_config_type * queue_config) {
  hash_free(queue_config->queue_drivers);
  util_safe_free(queue_config->job_script);
  free(queue_config);
}

const char * queue_config_get_queue_name(const queue_config_type * queue_config) {
    switch(queue_config->driver_type) {
        case LSF_DRIVER:
            return LSF_DRIVER_NAME;
        case RSH_DRIVER:
            return RSH_DRIVER_NAME;
        case LOCAL_DRIVER:
            return LOCAL_DRIVER_NAME;
        case TORQUE_DRIVER:
            return TORQUE_DRIVER_NAME;
        default:
            return NULL;
    }
    return NULL;
}


int queue_config_get_max_submit(queue_config_type * queue_config) {
    return queue_config->max_submit;
}

const char * queue_config_get_job_script(const queue_config_type * queue_config) {
  return queue_config->job_script;
}


bool queue_config_has_job_script( const queue_config_type * queue_config ) {
  if (queue_config->job_script)
    return true;
  else
    return false;
}


bool queue_config_set_job_script(queue_config_type * queue_config, const char * job_script) {
  if (util_is_executable(job_script)) {
    char * job_script_full_path = util_alloc_realpath(job_script);
    queue_config->job_script = util_realloc_string_copy(queue_config->job_script, job_script_full_path);
    free(job_script_full_path);
    return true;
  }
  else
    return false;
}


static void queue_config_set_queue_option(queue_config_type * queue_config, const char * driver_name, const char * option_key, const char * option_value) {
  if (queue_config_has_queue_driver(queue_config, driver_name)) {
    queue_driver_type * driver = queue_config_get_queue_driver(queue_config, driver_name);
    if (!queue_driver_set_option(driver, option_key, option_value))
      fprintf(stderr, "** Warning: Option:%s or its value is not recognized by driver:%s- ignored \n", option_key, driver_name);
  } else
    fprintf(stderr, "** Warning: Driver:%s not recognized - ignored \n", driver_name);
}


static void queue_config_add_queue_driver(queue_config_type * queue_config, const char * driver_name, queue_driver_type * driver) {
  hash_insert_hash_owned_ref(queue_config->queue_drivers, driver_name, driver, queue_driver_free__);
}


queue_driver_type * queue_config_get_queue_driver(const queue_config_type * queue_config, const char * driver_name) {
  return hash_get(queue_config->queue_drivers, driver_name);
}


void queue_config_create_queue_drivers(queue_config_type * queue_config) {
  queue_config_add_queue_driver(queue_config, LSF_DRIVER_NAME, queue_driver_alloc_LSF(NULL, NULL, NULL));
  queue_config_add_queue_driver(queue_config, TORQUE_DRIVER_NAME, queue_driver_alloc_TORQUE());
  queue_config_add_queue_driver(queue_config, RSH_DRIVER_NAME, queue_driver_alloc_RSH(NULL, NULL));
  queue_config_add_queue_driver(queue_config, LOCAL_DRIVER_NAME, queue_driver_alloc_local());
}


bool queue_config_has_queue_driver(const queue_config_type * queue_config, const char * driver_name) {
  return hash_has_key(queue_config->queue_drivers, driver_name);
}


static bool queue_config_init(queue_config_type * queue_config, const config_content_type * config_content)
{


  if (!queue_config->user_mode)
    queue_config_create_queue_drivers(queue_config);

  if (config_content_has_item(config_content, QUEUE_SYSTEM_KEY)) {
    const char * queue_system = config_content_get_value(config_content, QUEUE_SYSTEM_KEY);

    if (strcmp(queue_system, LSF_DRIVER_NAME) == 0) {
      queue_config->driver_type = LSF_DRIVER;
    } else if (strcmp(queue_system, RSH_DRIVER_NAME) == 0)
      queue_config->driver_type = RSH_DRIVER;
    else if (strcmp(queue_system, LOCAL_DRIVER_NAME) == 0)
      queue_config->driver_type = LOCAL_DRIVER;
    else if (strcmp(queue_system, TORQUE_DRIVER_NAME) == 0)
      queue_config->driver_type = TORQUE_DRIVER;
    else {
      util_abort("%s: queue system :%s not recognized \n", __func__, queue_system);
      queue_config->driver_type = NULL_DRIVER;
    }
  }


  if (config_content_has_item(config_content, JOB_SCRIPT_KEY))
    queue_config_set_job_script(queue_config, config_content_get_value_as_abspath(config_content, JOB_SCRIPT_KEY));

  if (config_content_has_item(config_content, MAX_SUBMIT_KEY))
    queue_config->max_submit = config_content_get_value_as_int(config_content, MAX_SUBMIT_KEY);

  /* Setting QUEUE_OPTIONS */
  for (int i = 0; i < config_content_get_occurences(config_content, QUEUE_OPTION_KEY); i++) {
    const stringlist_type * tokens = config_content_iget_stringlist_ref(config_content, QUEUE_OPTION_KEY, i);
    const char * driver_name = stringlist_iget(tokens, 0);
    const char * option_key = stringlist_iget(tokens, 1);
    char * option_value = stringlist_alloc_joined_substring(tokens, 2, stringlist_get_size(tokens), " ");

    // If it is essential to keep the exact number of spaces in the
    // option_value, it should be quoted with "" in the configuration
    // file.
    queue_config_set_queue_option(queue_config, driver_name, option_key, option_value);
    free( option_value );
  }

  return true;
}



job_driver_type queue_config_get_driver_type(const queue_config_type * queue_config) {
    return queue_config->driver_type;
}

void queue_config_add_queue_config_items(config_parser_type * parser, bool site_mode) {
  
}



void queue_config_add_config_items(config_parser_type * parser, bool site_mode) {

  {
    config_schema_item_type * item = config_add_schema_item(parser, MAX_SUBMIT_KEY, false);
    config_schema_item_set_argc_minmax(item, 1, 1);
    config_schema_item_iset_type(item, 0, CONFIG_INT);
  }

  {
    config_schema_item_type * item = config_add_schema_item(parser, QUEUE_SYSTEM_KEY, site_mode);
    config_schema_item_set_argc_minmax(item, 1, 1);
  }

  {
    config_schema_item_type * item = config_add_schema_item(parser, QUEUE_OPTION_KEY, false);
    config_schema_item_set_argc_minmax(item, 3, CONFIG_DEFAULT_ARG_MAX);
  }

  {
    config_schema_item_type * item = config_add_schema_item(parser, JOB_SCRIPT_KEY, false);
    config_schema_item_set_argc_minmax(item, 1, 1);
    config_schema_item_iset_type(item, 0, CONFIG_EXISTING_PATH);
  }
}
