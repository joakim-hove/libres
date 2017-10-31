/*
   Copyright (C) 2011  Statoil ASA, Norway.

   The file 'enkf_state.c' is part of ERT - Ensemble based Reservoir Tool.

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

#include <sys/types.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <pthread.h>

#include <ert/util/path_fmt.h>
#include <ert/util/thread_pool.h>
#include <ert/util/hash.h>
#include <ert/util/util.h>
#include <ert/util/arg_pack.h>
#include <ert/util/stringlist.h>
#include <ert/util/node_ctype.h>
#include <ert/util/timer.h>
#include <ert/util/time_t_vector.h>
#include <ert/util/rng.h>
#include <ert/res_util/subst_list.h>

#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_io_config.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_endian_flip.h>

#include <ert/sched/sched_file.h>

#include <ert/job_queue/forward_model.h>
#include <ert/job_queue/job_queue.h>
#include <ert/job_queue/queue_driver.h>
#include <ert/job_queue/ext_joblist.h>

#include <ert/enkf/enkf_node.h>
#include <ert/enkf/enkf_state.h>
#include <ert/enkf/enkf_types.h>
#include <ert/enkf/field.h>
#include <ert/enkf/field_config.h>
#include <ert/enkf/gen_kw.h>
#include <ert/enkf/summary.h>
#include <ert/enkf/gen_data.h>
#include <ert/enkf/enkf_fs.h>
#include <ert/enkf/ensemble_config.h>
#include <ert/enkf/model_config.h>
#include <ert/enkf/site_config.h>
#include <ert/enkf/ecl_config.h>
#include <ert/enkf/ert_template.h>
#include <ert/enkf/member_config.h>
#include <ert/enkf/enkf_defaults.h>
#include <ert/enkf/state_map.h>
#include <ert/res_util/res_log.h>
#include <ert/enkf/run_arg.h>
#include <ert/enkf/summary_key_matcher.h>
#include <ert/enkf/forward_load_context.h>
#include <ert/enkf/enkf_config_node.h>

#define  ENKF_STATE_TYPE_ID 78132





/**
   This struct contains various objects which the enkf_state needs
   during operation, which the enkf_state_object *DOES NOT* own. The
   struct only contains pointers to objects owned by (typically) the
   enkf_main object.

   If the enkf_state object writes to any of the objects in this
   struct that can be considered a serious *BUG*.

   The elements in this struct should not change during the
   application lifetime?
*/

typedef struct shared_info_struct {
  model_config_type           * model_config;      /* .... */
  ext_joblist_type            * joblist;           /* The list of external jobs which are installed - and *how* they should be run (with Python code) */
  const site_config_type      * site_config;
  ert_templates_type          * templates;
  const ecl_config_type       * ecl_config;
} shared_info_type;






/*****************************************************************/

struct enkf_state_struct {
  UTIL_TYPE_ID_DECLARATION;
  hash_type             * node_hash;
  subst_list_type       * subst_list;              /* This a list of key - value pairs which are used in a search-replace
                                                      operation on the ECLIPSE data file. Will at least contain the key INIT"
                                                      - which will describe initialization of ECLIPSE (EQUIL or RESTART).*/
  ensemble_config_type  * ensemble_config;         /* The config nodes for the enkf_node objects contained in node_hash. */

  shared_info_type      * shared_info;             /* Pointers to shared objects which is needed by the enkf_state object (read only). */
  member_config_type    * my_config;               /* Private config information for this member; not updated during a simulation. */
  rng_type              * rng;                     /* This is owned and managed by the external rng_manager. */
};

/*****************************************************************/


static void enkf_state_internalize_eclipse_state(enkf_state_type * enkf_state ,
						 forward_load_context_type * load_context,
						 const model_config_type * model_config ,
                                                 int report_step ,
                                                 bool store_vectors);

static void enkf_state_fread(enkf_state_type * enkf_state , enkf_fs_type * fs , int mask , int report_step );

static enkf_node_type * enkf_state_get_node(const enkf_state_type * enkf_state , const char * node_key);

/*****************************************************************/

static UTIL_SAFE_CAST_FUNCTION( enkf_state , ENKF_STATE_TYPE_ID )


static shared_info_type * shared_info_alloc(const site_config_type * site_config , model_config_type * model_config, const ecl_config_type * ecl_config , ert_templates_type * templates) {
  shared_info_type * shared_info = util_malloc(sizeof * shared_info );
  shared_info->joblist      = site_config_get_installed_jobs( site_config );
  shared_info->site_config  = site_config;
  shared_info->model_config = model_config;
  shared_info->templates    = templates;
  shared_info->ecl_config   = ecl_config;
  return shared_info;
}


static void shared_info_free(shared_info_type * shared_info) {
  /**
      Adding something here is a BUG - this object does
      not own anything.
  */
  free( shared_info );
}





/*****************************************************************/
/** Helper classes complete - starting on the enkf_state proper object. */
/*****************************************************************/


/*
  This function does not acces the nodes of the enkf_state object.
*/
void enkf_state_initialize(enkf_state_type * enkf_state , enkf_fs_type * fs , const stringlist_type * param_list, init_mode_type init_mode) {
  if (init_mode != INIT_NONE) {
    int iens = enkf_state_get_iens( enkf_state );
    state_map_type * state_map = enkf_fs_get_state_map( fs );
    realisation_state_enum current_state = state_map_iget(state_map, iens);
    if ((current_state == STATE_PARENT_FAILURE) && (init_mode != INIT_FORCE))
      return;
    else {
      const ensemble_config_type * ensemble_config = enkf_state->ensemble_config;
      for (int ip = 0; ip < stringlist_get_size(param_list); ip++) {
        const enkf_config_node_type * config_node = ensemble_config_get_node( ensemble_config , stringlist_iget(param_list, ip));
        enkf_node_type * param_node = enkf_node_alloc( config_node );
        node_id_type node_id = { .report_step = 0, .iens = iens };
        bool has_data = enkf_node_has_data(param_node, fs, node_id);

        if ((init_mode == INIT_FORCE) || (has_data == false) || (current_state == STATE_LOAD_FAILURE)) {
          if (enkf_node_initialize(param_node, iens, enkf_state->rng))
            enkf_node_store(param_node, fs, true, node_id);
        }

        enkf_node_free( param_node );
      }
      state_map_update_matching(state_map , iens , STATE_UNDEFINED | STATE_LOAD_FAILURE , STATE_INITIALIZED);
      enkf_fs_fsync(fs);
    }
  }
}







/*
  void enkf_state_set_iens(enkf_state_type * enkf_state , int iens) {
  enkf_state->my_iens = iens;
  }
*/

int enkf_state_get_iens(const enkf_state_type * enkf_state) {
  return member_config_get_iens( enkf_state->my_config );
}

member_config_type * enkf_state_get_member_config(const enkf_state_type * enkf_state) {
  return enkf_state->my_config;
}



subst_list_type * enkf_state_get_subst_list( enkf_state_type * enkf_state ) {
  return enkf_state->subst_list;
}

void enkf_state_add_subst_kw(enkf_state_type * enkf_state , const char * kw , const char * value , const char * doc_string) {
  char * tagged_key = util_alloc_sprintf( INTERNAL_DATA_KW_TAG_FORMAT , kw );
  subst_list_append_owned_ref(enkf_state->subst_list , tagged_key , util_alloc_string_copy(value) , doc_string);
  free(tagged_key);
}





/**
   Sets all the static subst keywords which will not change during the simulation.
*/
static void enkf_state_set_static_subst_kw(enkf_state_type * enkf_state) {

  {
    int    iens        = member_config_get_iens( enkf_state->my_config );
    char * iens_s      = util_alloc_sprintf("%d"   , iens);
    char * iens4_s     = util_alloc_sprintf("%04d" , iens);
    char * iensp1_s    = util_alloc_sprintf("%d"   , iens + 1);

    enkf_state_add_subst_kw(enkf_state , "IENS"        , iens_s      , NULL);
    enkf_state_add_subst_kw(enkf_state , "IENSP1"      , iensp1_s    , NULL);
    enkf_state_add_subst_kw(enkf_state , "IENS4"       , iens4_s     , NULL);

    free(iensp1_s);
    free(iens_s);
    free(iens4_s);
  }
}


static void enkf_state_add_nodes( enkf_state_type * enkf_state, const ensemble_config_type * ensemble_config) {
  stringlist_type * container_keys = stringlist_alloc_new();
  stringlist_type * keylist  = ensemble_config_alloc_keylist(ensemble_config);
  int keys        = stringlist_get_size(keylist);

  // 1: Add all regular nodes
  for (int ik = 0; ik < keys; ik++) {
    const char * key = stringlist_iget(keylist, ik);
    const enkf_config_node_type * config_node = ensemble_config_get_node(ensemble_config , key);
    if (enkf_config_node_get_impl_type( config_node ) == CONTAINER) {
      stringlist_append_ref( container_keys , key );
    } else
      enkf_state_add_node(enkf_state , key , config_node);
  }

  // 2: Add container nodes - must ensure that all other nodes have
  //    been added already (this implies that containers of containers
  //    will be victim of hash retrieval order problems ....

  for (int ik = 0; ik < stringlist_get_size( container_keys ); ik++) {
    const char * key = stringlist_iget(container_keys, ik);
    const enkf_config_node_type * config_node = ensemble_config_get_node(ensemble_config , key);
    enkf_state_add_node( enkf_state , key , config_node );
  }

  stringlist_free(keylist);
  stringlist_free( container_keys );
}


enkf_state_type * enkf_state_alloc(int iens,
                                   rng_type                  * rng ,
                                   const char                * casename ,
                                   model_config_type         * model_config,
                                   ensemble_config_type      * ensemble_config,
                                   const site_config_type    * site_config,
                                   const ecl_config_type     * ecl_config,

                                   ert_templates_type        * templates,
                                   subst_list_type           * subst_parent) {

  enkf_state_type * enkf_state  = util_malloc(sizeof *enkf_state );
  UTIL_TYPE_ID_INIT( enkf_state , ENKF_STATE_TYPE_ID );

  enkf_state->ensemble_config   = ensemble_config;
  enkf_state->shared_info       = shared_info_alloc(site_config , model_config , ecl_config , templates);

  enkf_state->node_hash         = hash_alloc();
  enkf_state->subst_list        = subst_list_alloc( subst_parent );
  enkf_state->rng               = rng;
  /*
    The user MUST specify an INIT_FILE, and for the first timestep the
    <INIT> tag in the data file will be replaced by an

    INCLDUE
    EQUIL_INIT_FILE

    statement. When considering the possibility of estimating EQUIL this
    require a real understanding of the treatment of paths:

    * If not estimating the EQUIL contacts, all members should use the
    same init_file. To ensure this the user must specify the ABSOLUTE
    PATH to a file containing initialization information.

    * If the user is estimating initial contacts, the INIT_FILE must
    point to the ecl_file of the EQUIL keyword, this must be a pure
    filename without any path component (as it will be generated by
    the EnKF program, and placed in the run_path directory). We could
    let the EnKF program use the ecl_file of the EQUIL keyword if it
    is present.

    The <INIT> key is actually initialized in the
    enkf_state_set_dynamic_subst_kw() function.
  */

  /**
     Adding all the subst_kw keywords here, with description. Listing
     all of them here in one go guarantees that we have control over
     the ordering (which is interesting because the substititions are
     done in a cascade like fashion). The user defined keywords are
     added first, so that these can refer to the built in keywords.
  */

  enkf_state_add_subst_kw(enkf_state , "RUNPATH"       , "---" , "The absolute path of the current forward model instance. ");
  enkf_state_add_subst_kw(enkf_state , "IENS"          , "---" , "The realisation number for this realization.");
  enkf_state_add_subst_kw(enkf_state , "IENS4"         , "---" , "The realization number for this realization - formated with %04d.");
  enkf_state_add_subst_kw(enkf_state , "ECLBASE"       , "---" , "The ECLIPSE basename for this realization.");
  enkf_state_add_subst_kw(enkf_state , "ECL_BASE"      , "---" , "Depreceated - use ECLBASE instead.");
  enkf_state_add_subst_kw(enkf_state , "SMSPEC"        , "---" , "The ECLIPSE SMSPEC file for this realization.");
  enkf_state_add_subst_kw(enkf_state , "TSTEP1"        , "---" , "The initial report step for this simulation.");
  enkf_state_add_subst_kw(enkf_state , "TSTEP2"        , "---" , "The final report step for this simulation.");
  enkf_state_add_subst_kw(enkf_state , "TSTEP1_04"     , "---" , "The initial report step for this simulation - formated with %04d.");
  enkf_state_add_subst_kw(enkf_state , "TSTEP2_04"     , "---" , "The final report step for this simulation - formated withh %04d.");
  enkf_state_add_subst_kw(enkf_state , "RESTART_FILE1" , "---" , "The ECLIPSE restart file this simulation starts with.");
  enkf_state_add_subst_kw(enkf_state , "RESTART_FILE2" , "---" , "The ECLIPSE restart file this simulation should end with.");
  enkf_state_add_subst_kw(enkf_state , "RANDINT"       , "---" , "Random integer value (depreceated: use __RANDINT__() instead).");
  enkf_state_add_subst_kw(enkf_state , "RANDFLOAT"     , "---" , "Random float value (depreceated: use __RANDFLOAT__() instead).");
  enkf_state_add_subst_kw(enkf_state , "INIT"          , "---" , "The code which will be inserted at the <INIT> tag");
  if (casename != NULL)
    enkf_state_add_subst_kw(enkf_state , "CASE" , casename , "The casename for this realization - as loaded from the CASE_TABLE file.");
  else
    enkf_state_add_subst_kw(enkf_state , "CASE" , "---" , "The casename for this realization - similar to ECLBASE.");

  enkf_state->my_config = member_config_alloc( iens, casename );
  enkf_state_set_static_subst_kw( enkf_state );
  enkf_state_add_nodes( enkf_state , ensemble_config );

  return enkf_state;
}









static void enkf_state_log_GEN_DATA_load( const enkf_node_type * enkf_node , int report_step , forward_load_context_type * load_context) {
  if (forward_load_context_accept_messages(load_context)) {
    char * load_file = enkf_config_node_alloc_infile(enkf_node_get_config( enkf_node ) , report_step);
    int data_size = gen_data_get_size( enkf_node_value_ptr( enkf_node ));
    char * msg = util_alloc_sprintf("Loaded GEN_DATA:%s instance for step:%d from file:%s size:%d" ,
                                    enkf_node_get_key( enkf_node ) ,
                                    report_step ,
                                    load_file ,
                                    data_size);

    forward_load_context_add_message(load_context, msg);

    free( msg );
    free( load_file );
  }
}


static void enkf_state_log_custom_kw_load(const enkf_node_type * enkf_node, int report_step, forward_load_context_type * load_context) {
  if (forward_load_context_accept_messages(load_context)) {
    char * load_file = enkf_config_node_alloc_infile(enkf_node_get_config(enkf_node), report_step);
    char * msg = util_alloc_sprintf("Loaded CUSTOM_KW: %s instance for step: %d from file: %s",
                                    enkf_node_get_key(enkf_node),
                                    report_step,
                                    load_file);

    forward_load_context_add_message(load_context, msg);

    free(msg);
    free(load_file);
  }
}



static int_vector_type * __enkf_state_get_time_index(enkf_fs_type * sim_fs, const ecl_sum_type * summary) {
  time_map_type * time_map = enkf_fs_get_time_map( sim_fs );
  time_map_summary_update( time_map , summary );
  return time_map_alloc_index_map( time_map , summary );
}


/*
 * Check if there are summary keys in the ensemble config that is not found in Eclipse. If this is the case, AND we
 * have observations for this key, we have a problem. Otherwise, just print a message to the log.
 */
static void enkf_state_check_for_missing_eclipse_summary_data(const summary_key_matcher_type * matcher, const ecl_smspec_type * smspec,
							      const enkf_state_type * enkf_state, forward_load_context_type * load_context, const int iens ) {

  stringlist_type * keys = summary_key_matcher_get_keys(matcher);

  for (int i = 0; i < stringlist_get_size(keys); i++) {

    const char *key = stringlist_iget(keys, i);

    if (ecl_smspec_has_general_var(smspec, key) || !summary_key_matcher_summary_key_is_required(matcher, key))
      continue;

    if (!ensemble_config_has_key(enkf_state->ensemble_config, key))
      continue;

    const enkf_config_node_type *config_node = ensemble_config_get_node(enkf_state->ensemble_config, key);
    if (enkf_config_node_get_num_obs(config_node) == 0) {
      res_log_add_fmt_message(LOG_INFO, NULL, "[%03d:----] Unable to find Eclipse data for summary key: %s, but have no observations either, so will continue.",
                              iens, key);
    } else {
      res_log_add_fmt_message(LOG_ERROR, NULL, "[%03d:----] Unable to find Eclipse data for summary key: %s, but have observation for this, job will fail.",
                              iens, key);
      forward_load_context_update_result(load_context, LOAD_FAILURE);
      if (forward_load_context_accept_messages(load_context)) {
        char *msg = util_alloc_sprintf("Failed to load vector: %s", key);
        forward_load_context_add_message(load_context, msg);
        free(msg);
      }
    }
  }

  stringlist_free(keys);
}

static bool enkf_state_internalize_dynamic_eclipse_results(enkf_state_type * enkf_state ,
                                                           forward_load_context_type * load_context ,
                                                           const model_config_type * model_config) {

  bool load_summary = ensemble_config_has_impl_type(enkf_state->ensemble_config, SUMMARY);
  const run_arg_type * run_arg = forward_load_context_get_run_arg( load_context );
  ensemble_config_type * ens_config = enkf_state->ensemble_config;
  const summary_key_matcher_type * matcher = ensemble_config_get_summary_key_matcher(ens_config);
  const ecl_sum_type * summary = forward_load_context_get_ecl_sum( load_context );
  int matcher_size = summary_key_matcher_get_size(matcher);

  if (load_summary || matcher_size > 0 || summary) {
    int load_start = run_arg_get_load_start( run_arg );

    if (load_start == 0) { /* Do not attempt to load the "S0000" summary results. */
      load_start++;
    }

    {
      enkf_fs_type * sim_fs = run_arg_get_sim_fs( run_arg );
      /** OK - now we have actually loaded the ecl_sum instance, or ecl_sum == NULL. */
      if (summary) {
        int_vector_type * time_index = __enkf_state_get_time_index(sim_fs, summary);

        /*
	  Now there are two related / conflicting(?) systems for
	  checking summary time consistency, both internally in the
	  time_map and also through the
	  model_config_report_step_compatible() function.
        */

        /*Check the loaded summary against the reference ecl_sum_type */
        if (!model_config_report_step_compatible(model_config, summary))
          forward_load_context_update_result(load_context, REPORT_STEP_INCOMPATIBLE);


        /* The actual loading internalizing - from ecl_sum -> enkf_node. */
        const int iens   = run_arg_get_iens( run_arg );
        const int step2  = ecl_sum_get_last_report_step( summary );  /* Step2 is just taken from the number of steps found in the summary file. */

        int_vector_iset_block( time_index , 0 , load_start , -1 );
        int_vector_resize( time_index , step2 + 1);

        const ecl_smspec_type * smspec = ecl_sum_get_smspec(summary);

        for(int i = 0; i < ecl_smspec_num_nodes(smspec); i++) {
	  const smspec_node_type * smspec_node = ecl_smspec_iget_node(smspec, i);
	  if (smspec_node_is_valid( smspec_node )) {
	    const char * key = smspec_node_get_gen_key1(smspec_node);

              if(summary_key_matcher_match_summary_key(matcher, key)) {
                summary_key_set_type * key_set = enkf_fs_get_summary_key_set(sim_fs);
                summary_key_set_add_summary_key(key_set, key);

                enkf_config_node_type * config_node = ensemble_config_get_or_create_summary_node(ens_config, key);
                enkf_node_type * node = enkf_node_alloc( config_node );

                enkf_node_try_load_vector( node , sim_fs , iens );  // Ensure that what is currently on file is loaded before we update.

                enkf_node_forward_load_vector( node , load_context , time_index);
                enkf_node_store_vector( node , sim_fs , iens );
		enkf_node_free( node );
              }
            }
        }

        int_vector_free( time_index );

        /*
        Check if some of the specified keys are missing from the Eclipse data, and if there are observations for them. That is a problem.
        */
        enkf_state_check_for_missing_eclipse_summary_data(matcher, smspec, enkf_state, load_context, iens);

        return true;
      } else {
        res_log_add_fmt_message(LOG_WARNING, NULL, "Could not load ECLIPSE summary data from %s - this will probably fail later ... ",run_arg_get_runpath( run_arg ));
        return false;
      }
    }
  } else {
    return true;
  }
}






static void enkf_state_internalize_custom_kw(enkf_state_type * enkf_state,
                                             forward_load_context_type * load_context ,
                                             const model_config_type * model_config) {

  const run_arg_type * run_arg     = forward_load_context_get_run_arg( load_context );
  enkf_fs_type *sim_fs             = run_arg_get_sim_fs(run_arg);
  const int iens                   = run_arg_get_iens( run_arg );
  stringlist_type * custom_kw_keys = ensemble_config_alloc_keylist_from_impl_type(enkf_state->ensemble_config, CUSTOM_KW);
  const int report_step            = 0;

  custom_kw_config_set_type * config_set = enkf_fs_get_custom_kw_config_set(sim_fs);
  custom_kw_config_set_reset(config_set);

  for (int ikey=0; ikey < stringlist_get_size(custom_kw_keys); ikey++) {
    const char* custom_kw_key = stringlist_iget(custom_kw_keys, ikey);
    enkf_node_type * node = enkf_state_get_node(enkf_state, custom_kw_key);

    if (enkf_node_vector_storage(node))
      util_abort("%s: Vector storage not correctly implemented for CUSTOM_KW\n", __func__);

    if (enkf_node_internalize(node, report_step) && enkf_node_has_func(node, forward_load_func)) {
      if (enkf_node_forward_load(node, load_context)) {
        node_id_type node_id = {.report_step = report_step, .iens = iens };

        enkf_node_store(node, sim_fs, false, node_id);

        const enkf_config_node_type * config_node = enkf_node_get_config(node);
        const custom_kw_config_type * custom_kw_config = (custom_kw_config_type*) enkf_config_node_get_ref(config_node);
        custom_kw_config_set_add_config(config_set, custom_kw_config);
        enkf_state_log_custom_kw_load(node, report_step, load_context);
      } else {
        forward_load_context_update_result(load_context, LOAD_FAILURE);
        res_log_add_fmt_message(LOG_ERROR, stderr,
                                "[%03d:%04d] Failed load data for CUSTOM_KW node: %s.",
                                iens, report_step, enkf_node_get_key(node));

        if (forward_load_context_accept_messages(load_context)) {
          char * msg = util_alloc_sprintf("Failed to load: %s at step: %d", enkf_node_get_key(node), report_step);
          forward_load_context_add_message(load_context , msg);
          free( msg );
        }
      }
    }
  }

  stringlist_free(custom_kw_keys);
}



static void enkf_state_internalize_GEN_DATA(const ensemble_config_type * ens_config, 
                                            forward_load_context_type * load_context ,
                                            const model_config_type * model_config ,
                                            int last_report) {

  stringlist_type * keylist_GEN_DATA = ensemble_config_alloc_keylist_from_impl_type(ens_config, GEN_DATA);


  if (stringlist_get_size( keylist_GEN_DATA) > 0) 
    if (last_report <= 0)
      res_log_add_message( LOG_WARNING, NULL , "Trying to load GEN_DATA without properly set last_report - will only look for step 0 data.", false);

  const run_arg_type * run_arg            = forward_load_context_get_run_arg( load_context );
  enkf_fs_type * sim_fs                   = run_arg_get_sim_fs( run_arg );
  const int  iens                         = run_arg_get_iens( run_arg );

  for (int ikey=0; ikey < stringlist_get_size( keylist_GEN_DATA ); ikey++) {
    const enkf_config_node_type * config_node = ensemble_config_get_node( ens_config , stringlist_iget( keylist_GEN_DATA , ikey));

    /*
      This for loop should probably be changed to use the report
      steps configured in the gen_data_config object, instead of
      spinning through them all.
    */
    for (int report_step = run_arg_get_load_start( run_arg ); report_step <= util_int_max(0, last_report); report_step++) {
      if (!enkf_config_node_internalize(config_node , report_step))
        continue;

      forward_load_context_select_step(load_context, report_step);
      enkf_node_type * node = enkf_node_alloc( config_node );

      if (enkf_node_forward_load(node , load_context )) {
        node_id_type node_id = {.report_step = report_step ,
                                .iens = iens };

        enkf_node_store( node , sim_fs, false , node_id );
        enkf_state_log_GEN_DATA_load( node , report_step , load_context);
      } else {
        forward_load_context_update_result(load_context, LOAD_FAILURE);
        res_log_add_fmt_message(LOG_ERROR, stderr,
                                "[%03d:%04d] Failed load data for GEN_DATA node:%s.",
                                iens, report_step, enkf_node_get_key(node));

        if (forward_load_context_accept_messages(load_context)) {
          char * msg = util_alloc_sprintf("Failed to load: %s at step:%d",
                                          enkf_node_get_key(node), report_step);
          forward_load_context_add_message(load_context, msg);
          free( msg );
        }
      }
      enkf_node_free( node );
    }
  }
  stringlist_free( keylist_GEN_DATA );
}





static forward_load_context_type * enkf_state_alloc_load_context(const enkf_state_type * state,
                                                                 run_arg_type * run_arg,
                                                                 stringlist_type * messages) {
  bool load_summary = false;
  const ensemble_config_type * ens_config = state->ensemble_config;
  const summary_key_matcher_type * matcher = ensemble_config_get_summary_key_matcher(ens_config);
  if (summary_key_matcher_get_size(matcher) > 0)
    load_summary = true;

  if (ensemble_config_has_GEN_DATA(ens_config))
    load_summary = true;

  if (ensemble_config_has_impl_type(ens_config, SUMMARY))
    load_summary = true;

  forward_load_context_type * load_context;
  const ecl_config_type * ecl_config = state->shared_info->ecl_config;

  load_context = forward_load_context_alloc(run_arg,
                                            load_summary,
                                            ecl_config,
                                            messages);
  return load_context;

}


/**
   This function loads the results from a forward simulations from report_step1
   to report_step2. The details of what to load are in model_config and the
   spesific nodes for special cases.

   Will mainly be called at the end of the forward model, but can also
   be called manually from external scope.
*/
static int enkf_state_internalize_results(enkf_state_type * enkf_state , run_arg_type * run_arg , stringlist_type * msg_list) {
  const ensemble_config_type * ens_config = enkf_state->ensemble_config;
  const model_config_type * model_config = enkf_state->shared_info->model_config;
  forward_load_context_type * load_context = enkf_state_alloc_load_context( enkf_state , run_arg , msg_list);
  int report_step;

  /*
    The timing information - i.e. mainly what is the last report step
    in these results are inferred from the loading of summary results,
    hence we must load the summary results first.
  */

  enkf_state_internalize_dynamic_eclipse_results(enkf_state ,
                                                 load_context ,
                                                 model_config);

  enkf_fs_type * sim_fs = run_arg_get_sim_fs( run_arg );
  int last_report = time_map_get_last_step( enkf_fs_get_time_map( sim_fs ));
  if (last_report < 0)
    last_report = model_config_get_last_history_restart( model_config );

  /* Ensure that the last step is internalized? */
  if (last_report > 0)
    model_config_set_internalize_state( model_config , last_report);

  for (report_step = run_arg_get_load_start(run_arg); report_step <= last_report; report_step++) {
    bool store_vectors = (report_step == last_report);
    if (model_config_load_state( model_config , report_step))
      enkf_state_internalize_eclipse_state(enkf_state ,
                                           load_context ,
                                           model_config ,
                                           report_step ,
                                           store_vectors);
  }

  enkf_state_internalize_GEN_DATA(ens_config , load_context , model_config , last_report);
  enkf_state_internalize_custom_kw(enkf_state, load_context , model_config);

  int result = forward_load_context_get_result(load_context);
  forward_load_context_free( load_context );
  return result;
}


int enkf_state_forward_init(const ensemble_config_type * ens_config, 
                            run_arg_type * run_arg) {

  int result = 0;
  if (run_arg_get_step1(run_arg) == 0) {
    int iens = run_arg_get_iens( run_arg );
    hash_iter_type * iter = ensemble_config_alloc_hash_iter( ens_config );
    while ( !hash_iter_is_complete(iter) ) {
      enkf_config_node_type * config_node = hash_iter_get_next_value(iter);
      if (enkf_config_node_use_forward_init(config_node)) {
	enkf_node_type * node = enkf_node_alloc( config_node );
        enkf_fs_type * sim_fs = run_arg_get_sim_fs( run_arg );
        node_id_type node_id = {.report_step = 0 ,
                                .iens = iens };


        /*
           Will not reinitialize; i.e. it is essential that the
           forward model uses the state given from the stored
           instance, and not from the current run of e.g. RMS.
        */

        if (!enkf_node_has_data( node , sim_fs , node_id)) {
          if (enkf_node_forward_init(node , run_arg_get_runpath( run_arg ) , iens ))
            enkf_node_store( node , sim_fs , false , node_id );
          else {
            char * init_file = enkf_config_node_alloc_initfile( enkf_node_get_config( node ) , run_arg_get_runpath(run_arg) , iens );

            if (init_file && !util_file_exists( init_file ))
              fprintf(stderr,"File not found: %s - failed to initialize node: %s\n", init_file , enkf_node_get_key( node ));
            else
              fprintf(stderr,"Failed to initialize node: %s\n", enkf_node_get_key( node ));

            util_safe_free( init_file );
            result |= LOAD_FAILURE;
          }
        }
	enkf_node_free( node );
      }
    }
    hash_iter_free( iter );
  }
  return result;
}



int enkf_state_load_from_forward_model(enkf_state_type * enkf_state ,
                                       run_arg_type * run_arg ,
                                       stringlist_type * msg_list) {

  int result = 0;

  if (ensemble_config_have_forward_init( enkf_state->ensemble_config ))
    result |= enkf_state_forward_init( enkf_state->ensemble_config , run_arg );

  result |= enkf_state_internalize_results( enkf_state , run_arg , msg_list );
  state_map_type * state_map = enkf_fs_get_state_map( run_arg_get_sim_fs( run_arg ) );
  int iens = run_arg_get_iens( run_arg );
  if (result & LOAD_FAILURE)
    state_map_iset( state_map , iens , STATE_LOAD_FAILURE);
  else
    state_map_iset( state_map , iens , STATE_HAS_DATA);

  return result;
}


/**
   Observe that this does not return the loadOK flag; it will load as
   good as it can all the data it should, and be done with it.
*/

void * enkf_state_load_from_forward_model_mt( void * arg ) {
  arg_pack_type * arg_pack     = arg_pack_safe_cast( arg );
  enkf_state_type * enkf_state = enkf_state_safe_cast(arg_pack_iget_ptr( arg_pack  , 0 ));
  run_arg_type * run_arg       = arg_pack_iget_ptr( arg_pack  , 1 );
  stringlist_type * msg_list   = arg_pack_iget_ptr( arg_pack  , 2 );
  bool manual_load             = arg_pack_iget_bool( arg_pack , 3 );
  int * result                 = arg_pack_iget_ptr( arg_pack  , 4 );
  int iens                     = run_arg_get_iens( run_arg );

  if (manual_load)
    state_map_update_undefined(enkf_fs_get_state_map( run_arg_get_sim_fs(run_arg) ) , iens , STATE_INITIALIZED);

  *result = enkf_state_load_from_forward_model( enkf_state , run_arg , msg_list );
  if (*result & REPORT_STEP_INCOMPATIBLE) {
    // If refcase has been used for observations: crash and burn.
    fprintf(stderr,"** Warning the timesteps in refcase and current simulation are not in accordance - something wrong with schedule file?\n");
    *result -= REPORT_STEP_INCOMPATIBLE;
  }

  if (manual_load) {
    printf(".");
    fflush(stdout);
  }
  return NULL;
}





void enkf_state_free(enkf_state_type *enkf_state) {
  hash_free(enkf_state->node_hash);
  subst_list_free(enkf_state->subst_list);
  member_config_free(enkf_state->my_config);
  shared_info_free(enkf_state->shared_info);
  free(enkf_state);
}



/**
   This function will set all the subst_kw key=value pairs which
   change with report step.
*/

static void enkf_state_set_dynamic_subst_kw(enkf_state_type * enkf_state , const run_arg_type * run_arg) {
  const char * run_path = run_arg_get_runpath(run_arg);
  const char * job_name = run_arg_get_job_name( run_arg );
  int step1 = run_arg_get_step1(run_arg);
  int step2 = run_arg_get_step2(run_arg);
  
  const ecl_config_type * ecl_config = enkf_state->shared_info->ecl_config;
  const bool fmt_file  = ecl_config_get_formatted( ecl_config );


  if (run_path != NULL) {
    /** Make absolutely sure the path available as <RUNPATH> is absolute. */
    char * abs_runpath = util_alloc_realpath( run_path );
    enkf_state_add_subst_kw(enkf_state , "RUNPATH"       , abs_runpath      , NULL);
    free( abs_runpath );
  }

  enkf_state_add_subst_kw(enkf_state , "ECL_BASE"    , job_name , NULL);
  enkf_state_add_subst_kw(enkf_state , "ECLBASE"     , job_name , NULL);


  /* Time step */
  char * step1_s           = util_alloc_sprintf("%d" , step1);
  char * step2_s           = util_alloc_sprintf("%d" , step2);
  char * step1_s04         = util_alloc_sprintf("%04d" , step1);
  char * step2_s04         = util_alloc_sprintf("%04d" , step2);

  enkf_state_add_subst_kw(enkf_state , "TSTEP1"        , step1_s       , NULL);
  enkf_state_add_subst_kw(enkf_state , "TSTEP2"        , step2_s       , NULL);
  enkf_state_add_subst_kw(enkf_state , "TSTEP1_04"     , step1_s04     , NULL);
  enkf_state_add_subst_kw(enkf_state , "TSTEP2_04"     , step2_s04     , NULL);

  free(step1_s);
  free(step2_s);
  free(step1_s04);
  free(step2_s04);


  /* Restart file names and RESTART keyword in datafile. */
  if (ecl_config_have_eclbase( ecl_config )) {
    printf("run_mde:%d \n",run_arg_get_run_mode(run_arg));
    if (run_arg_get_run_mode(run_arg) != INIT_ONLY) {
      const bool fmt_file  = ecl_config_get_formatted( ecl_config );
      char * restart_file1 = ecl_util_alloc_filename(NULL , job_name , ECL_RESTART_FILE , fmt_file , step1);
      char * restart_file2 = ecl_util_alloc_filename(NULL , job_name , ECL_RESTART_FILE , fmt_file , step2);

      enkf_state_add_subst_kw(enkf_state , "RESTART_FILE1" , restart_file1 , NULL);
      enkf_state_add_subst_kw(enkf_state , "RESTART_FILE2" , restart_file2 , NULL);

      free(restart_file1);
      free(restart_file2);

      if (step1 > 0) {
        char * data_initialize = util_alloc_sprintf("RESTART\n   \'%s\'  %d  /\n" , job_name , step1);
        enkf_state_add_subst_kw(enkf_state , "INIT" , data_initialize , NULL);
        free(data_initialize);
      }
    }
  }


  /**
     The <INIT> magic string:
  */
  if (step1 == 0) {
    const char * init_file = ecl_config_get_equil_init_file(ecl_config);
    if (init_file != NULL) {
      char * tmp_include = util_alloc_sprintf("INCLUDE\n   \'%s\' /\n",init_file);
      enkf_state_add_subst_kw(enkf_state , "INIT" , tmp_include , NULL);
      free(tmp_include);
    } /*
         if init_file == NULL that means the user has not supplied the INIT_SECTION keyword,
         and the EQUIL (or whatever) info to initialize the model is inlined in the datafile.
      */
  }


  /**
     Adding keys for <RANDINT> and <RANDFLOAT> - these are only
     added for backwards compatibility, should be replaced with
     prober function callbacks.
  */
  char * randint_value    = util_alloc_sprintf( "%u"      , rng_forward( enkf_state->rng ));
  char * randfloat_value  = util_alloc_sprintf( "%12.10f" , rng_get_double( enkf_state->rng ));

  enkf_state_add_subst_kw( enkf_state , "RANDINT"   , randint_value   , NULL);
  enkf_state_add_subst_kw( enkf_state , "RANDFLOAT" , randfloat_value , NULL);

  free( randint_value );
  free( randfloat_value );
}








/**
   init_step    : The parameters are loaded from this EnKF/report step.
   report_step1 : The simulation should start from this report step;
                  dynamic data are loaded from this step.
   report_step2 : The simulation should stop at this report step. (unless run_mode == ENSEMBLE_PREDICTION - where it just runs til end.)

   For a normal EnKF run we well have init_step == report_step1, but
   in the case where we want rerun from the beginning with updated
   parameters, they will be different. If init_step != report_step1,
   it is required that report_step1 == 0; otherwise the dynamic data
   will become completely inconsistent. We just don't allow that!
*/
void enkf_state_init_eclipse(enkf_state_type *enkf_state, const run_arg_type * run_arg ) {
  const ecl_config_type * ecl_config = enkf_state->shared_info->ecl_config;

  util_make_path(run_arg_get_runpath( run_arg ));
  if (ecl_config_get_schedule_target( ecl_config ) != NULL) {

    char * schedule_file_target = util_alloc_filename(run_arg_get_runpath( run_arg ),
                                                      ecl_config_get_schedule_target( ecl_config ),
                                                      NULL);
    char * schedule_file_target_path = util_split_alloc_dirname(schedule_file_target);
    util_make_path(schedule_file_target_path);
    free(schedule_file_target_path);

    sched_file_fprintf( ecl_config_get_sched_file( ecl_config ) , schedule_file_target);

    free(schedule_file_target);
  }


  /**
     For reruns of various kinds the parameters and the state are
     generally loaded from different timesteps:
  */
  enkf_fs_type * sim_fs = run_arg_get_sim_fs( run_arg );
  /* Loading parameter information: loaded from timestep: run_arg->init_step_parameters. */
  enkf_state_fread(enkf_state , sim_fs , PARAMETER , 0);


  enkf_state_set_dynamic_subst_kw(  enkf_state , run_arg );
  ert_templates_instansiate( enkf_state->shared_info->templates , run_arg_get_runpath( run_arg ) , enkf_state->subst_list );
  enkf_state_ecl_write( enkf_state , run_arg , sim_fs);

  if (ecl_config_have_eclbase( ecl_config )) {

    /* Writing the ECLIPSE data file. */
    if (ecl_config_get_data_file( ecl_config ) != NULL) {
      char * data_file = ecl_util_alloc_filename(run_arg_get_runpath( run_arg ),
                                                 run_arg_get_job_name( run_arg ),
                                                 ECL_DATA_FILE,
                                                 true,
                                                 -1);
      subst_list_filter_file(enkf_state->subst_list,
                             ecl_config_get_data_file(ecl_config),
                             data_file);
      free( data_file );
    }
  }

  mode_t umask = site_config_get_umask(enkf_state->shared_info->site_config);

  /* This is where the job script is created */
  forward_model_formatted_fprintf( model_config_get_forward_model( enkf_state->shared_info->model_config ) ,
                                   run_arg_get_run_id( run_arg ),
                                   run_arg_get_runpath( run_arg ) ,
                                   model_config_get_data_root( enkf_state->shared_info->model_config ) ,
                                   enkf_state->subst_list,
                                   umask);

}



/**
    Observe that if run_arg == false, this routine will return with
    job_completeOK == true, that might be a bit misleading.

    Observe that if an internal retry is performed, this function will
    be called several times - MUST BE REENTRANT.
*/

static bool enkf_state_complete_forward_modelOK(enkf_state_type * enkf_state , run_arg_type * run_arg) {
  const int iens = run_arg_get_iens( run_arg );
  int result;


  /**
     The queue system has reported that the run is OK, i.e. it has
     completed and produced the targetfile it should. We then check
     in this scope whether the results can be loaded back; if that
     is OK the final status is updated, otherwise: restart.
  */
  res_log_add_fmt_message(LOG_INFO, NULL,
                          "[%03d:%04d-%04d] Forward model complete - starting to load results.",
                          iens , run_arg_get_step1(run_arg), run_arg_get_step2(run_arg));
  result = enkf_state_load_from_forward_model(enkf_state , run_arg , NULL);

  if (result & REPORT_STEP_INCOMPATIBLE) {
    // If refcase has been used for observations: crash and burn.
     fprintf(stderr,"** Warning the timesteps in refcase and current simulation are not in accordance - something wrong with schedule file?\n");
     result -= REPORT_STEP_INCOMPATIBLE;
  }


  if (result == 0) {
    /*
      The loading succeded - so this is a howling success! We set
      the main status to JOB_QUEUE_ALL_OK and inform the queue layer
      about the success. In addition we set the simple status
      (should be avoided) to JOB_RUN_OK.
    */
    run_arg_set_run_status( run_arg , JOB_RUN_OK);
    res_log_add_fmt_message(LOG_INFO, NULL,
                            "[%03d:%04d-%04d] Results loaded successfully.",
                            iens , run_arg_get_step1(run_arg), run_arg_get_step2(run_arg));

    run_arg_complete_run(run_arg);              /* free() on runpath */
  }

  return (result == 0) ? true : false;
}


bool enkf_state_complete_forward_modelOK__(void * arg ) {
  arg_pack_type * arg_pack = arg_pack_safe_cast( arg );
  enkf_state_type * enkf_state = enkf_state_safe_cast( arg_pack_iget_ptr( arg_pack , 0 ));
  run_arg_type * run_arg = run_arg_safe_cast( arg_pack_iget_ptr( arg_pack , 1 ));

  return enkf_state_complete_forward_modelOK( enkf_state , run_arg);
}



static bool enkf_state_complete_forward_model_EXIT_handler__(enkf_state_type * enkf_state,
                                                             run_arg_type * run_arg) {
  const int iens = run_arg_get_iens( run_arg );
  res_log_add_fmt_message(LOG_ERROR, NULL,
                          "[%03d:%04d-%04d] FAILED COMPLETELY.",
                          iens, run_arg_get_step1(run_arg), run_arg_get_step2(run_arg));

  if (run_arg_get_run_status(run_arg) != JOB_LOAD_FAILURE)
    run_arg_set_run_status( run_arg , JOB_RUN_FAILURE);

  state_map_type * state_map = enkf_fs_get_state_map(run_arg_get_sim_fs( run_arg ));
  state_map_iset(state_map, iens, STATE_LOAD_FAILURE);
  return false;
}


static bool enkf_state_complete_forward_model_EXIT_handler(void * arg) {
  arg_pack_type * arg_pack = arg_pack_safe_cast( arg );

  enkf_state_type * enkf_state = enkf_state_safe_cast( arg_pack_iget_ptr( arg_pack , 0 ) );
  run_arg_type * run_arg = run_arg_safe_cast( arg_pack_iget_ptr( arg_pack , 1 ) );

  return enkf_state_complete_forward_model_EXIT_handler__( enkf_state , run_arg);
}


bool enkf_state_complete_forward_modelEXIT__(void * arg ) {
  return enkf_state_complete_forward_model_EXIT_handler(arg);
}


/**
    This function is called when:

     1. The external queue system has said that everything is OK; BUT
        the ert layer failed to load all the data.

     2. The external queue system has seen the job fail.

    The parameter and state variables will be resampled before
    retrying. And all random elements in templates+++ will be
    resampled.
*/



static void enkf_state_internal_retry(enkf_state_type * enkf_state , run_arg_type * run_arg) {
  const int iens                          = run_arg_get_iens( run_arg );

  res_log_add_fmt_message(LOG_ERROR, NULL,
                          "[%03d:%04d - %04d] Forward model failed.",
                          iens, run_arg_get_step1(run_arg) , run_arg_get_step2(run_arg));
  if (run_arg_can_retry( run_arg ) ) {
    res_log_add_fmt_message(LOG_ERROR, NULL , "[%03d] Resampling and resubmitting realization." ,iens);
    /* Reinitialization of the nodes */
    stringlist_type * init_keys = ensemble_config_alloc_keylist_from_var_type( enkf_state->ensemble_config , PARAMETER );
    for (int ikey=0; ikey < stringlist_get_size( init_keys ); ikey++) {
      enkf_node_type * node = enkf_state_get_node( enkf_state , stringlist_iget( init_keys , ikey) );
      enkf_node_initialize( node , iens , enkf_state->rng );
    }
    stringlist_free( init_keys );

    /* Possibly clear the directory and do a FULL rewrite of ALL the necessary files. */
    enkf_state_init_eclipse( enkf_state , run_arg  );
    run_arg_increase_submit_count( run_arg );
  }
}


bool enkf_state_complete_forward_modelRETRY__(void * arg ) {
  arg_pack_type * arg_pack = arg_pack_safe_cast( arg );
  enkf_state_type * enkf_state = enkf_state_safe_cast( arg_pack_iget_ptr( arg_pack , 0 ) );
  run_arg_type * run_arg = run_arg_safe_cast( arg_pack_iget_ptr( arg_pack , 1 ) );

  if (run_arg_can_retry(run_arg)) {
    enkf_state_internal_retry(enkf_state, run_arg);
    return true;
  }

  return false;
}



/*****************************************************************/


rng_type * enkf_state_get_rng( const enkf_state_type * enkf_state ) {
  return enkf_state->rng;
}

unsigned int enkf_state_get_random( enkf_state_type * enkf_state ) {
  return rng_forward( enkf_state->rng );
}



const ensemble_config_type * enkf_state_get_ensemble_config( const enkf_state_type * enkf_state ) {
  return enkf_state->ensemble_config;
}


/**
  This function writes out all the files needed by an ECLIPSE simulation, this
  includes the restart file, and the various INCLUDE files corresponding to
  parameters estimated by EnKF.

  The writing of restart file is delegated to enkf_state_write_restart_file().
*/

void enkf_state_ecl_write(enkf_state_type * enkf_state, const run_arg_type * run_arg , enkf_fs_type * fs) {
  /**
     This iteration manipulates the hash (thorugh the enkf_state_del_node() call)

     -----------------------------------------------------------------------------------------
     T H I S  W I L L  D E A D L O C K  I F  T H E   H A S H _ I T E R  A P I   I S   U S E D.
     -----------------------------------------------------------------------------------------
  */

  const shared_info_type * shared_info   = enkf_state->shared_info;
  const model_config_type * model_config = shared_info->model_config;
  int iens                               = enkf_state_get_iens( enkf_state );
  const char * base_name                 = model_config_get_gen_kw_export_name(model_config);
  value_export_type * export             = value_export_alloc( run_arg_get_runpath( run_arg ), base_name );

  stringlist_type * key_list = ensemble_config_alloc_keylist_from_var_type( enkf_state->ensemble_config , PARAMETER );
  for (int ikey = 0; ikey < stringlist_get_size( key_list ); ikey++) {
    enkf_config_node_type * config_node = ensemble_config_get_node( enkf_state->ensemble_config, stringlist_iget( key_list , ikey));
    enkf_node_type * enkf_node = enkf_node_alloc( config_node );
    bool forward_init = enkf_node_use_forward_init( enkf_node );
    node_id_type node_id = {.report_step = run_arg_get_step1(run_arg),
                            .iens = iens };

    if ((run_arg_get_step1(run_arg) == 0) && (forward_init)) {

      if (enkf_node_has_data( enkf_node , fs , node_id))
        enkf_node_load(enkf_node, fs, node_id);
      else
        continue;
    } else
      enkf_node_load(enkf_node, fs, node_id);

    enkf_node_ecl_write(enkf_node , run_arg_get_runpath( run_arg ) , export , run_arg_get_step1(run_arg));
    enkf_node_free(enkf_node);
  }
  value_export( export );

  value_export_free( export );
  stringlist_free( key_list );
}


#include "enkf_state_nodes.c"
