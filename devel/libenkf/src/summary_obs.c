/**
   See the overview documentation of the observation system in enkf_obs.c
*/
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <util.h>
#include <stdio.h>
#include <summary_obs.h>
#include <obs_data.h>
#include <meas_matrix.h>
#include <summary.h>
#include <active_list.h>


#define SUMMARY_OBS_TYPE_ID 66103
#define OBS_SIZE 1

struct summary_obs_struct {
  UTIL_TYPE_ID_DECLARATION;
  char    * summary_key;    /** The observation, in summary.x syntax, e.g. GOPR:FIELD.    */
  char    * obs_key;
  
  double    value;          /** Observation value. */
  double    std;            /** Standard deviation of observation. */
};





/**
  This function allocates a summary_obs instance. The summary_key
  string should be of the format used by the summary.x program.
  E.g., WOPR:P4 would condition on WOPR in well P4.

  Observe that this format is currently *not* checked before the actual
  observation time.

  TODO
  Should check summary_key on alloc.
*/
summary_obs_type * summary_obs_alloc(const char   * summary_key,
                                     const char   * obs_key , 
                                     double value ,
                                     double std)
{
  summary_obs_type * obs = util_malloc(sizeof * obs , __func__);
  UTIL_TYPE_ID_INIT( obs , SUMMARY_OBS_TYPE_ID )

    obs->summary_key   = util_alloc_string_copy( summary_key );
  obs->obs_key       = util_alloc_string_copy( obs_key );
  obs->value         = value;
  obs->std           = std;
  
  return obs;
}


static UTIL_SAFE_CAST_FUNCTION_CONST(summary_obs   , SUMMARY_OBS_TYPE_ID);
static UTIL_SAFE_CAST_FUNCTION(summary_obs   , SUMMARY_OBS_TYPE_ID);
UTIL_IS_INSTANCE_FUNCTION(summary_obs , SUMMARY_OBS_TYPE_ID);


void summary_obs_free(summary_obs_type * summary_obs) {
  free(summary_obs->summary_key);
  free(summary_obs);
}







const char * summary_obs_get_summary_key(const summary_obs_type * summary_obs)
{
  return summary_obs->summary_key;
}


/**
   Hardcodes an assumption that the size of summary data|observations
   is always one; i.e. PARTLY_ACTIVE and ALL_ACTIVE are treated in the
   same manner.
*/
void summary_obs_get_observations(const summary_obs_type * summary_obs,
				  int                      restart_nr,
				  obs_data_type          * obs_data,
				  const active_list_type * __active_list) {

  int active_size              = active_list_get_active_size( __active_list , OBS_SIZE );
  if (active_size == 1) {
    obs_block_type * obs_block   = obs_data_add_block( obs_data , summary_obs->obs_key , OBS_SIZE );
    obs_block_iset( obs_block , 0 , summary_obs->value , summary_obs->std );
  }
}



void summary_obs_measure(const summary_obs_type * obs, const summary_type * summary, int report_step , int iens , meas_matrix_type * meas_matrix , const active_list_type * __active_list) {
  int active_size = active_list_get_active_size( __active_list , OBS_SIZE );
  if (active_size == 1) {
    meas_block_type * meas_block = meas_matrix_add_block( meas_matrix , obs->obs_key , active_size );
    meas_block_iset( meas_block , iens , 0 , summary_get(summary));
  }
}

 

double summary_obs_chi2(const summary_obs_type * obs,
			const summary_type     * summary) {
  double x = (summary_get(summary) - obs->value) / obs->std;
  return x*x;
}



void summary_obs_user_get(const summary_obs_type * summary_obs , const char * index_key , double * value , double * std, bool * valid) {
  *valid = true;
  *value = summary_obs->value;
  *std   = summary_obs->std;
}



/*****************************************************************/

VOID_FREE(summary_obs)
VOID_GET_OBS(summary_obs)
VOID_USER_GET_OBS(summary_obs)
VOID_MEASURE(summary_obs , summary)
VOID_CHI2(summary_obs , summary)
