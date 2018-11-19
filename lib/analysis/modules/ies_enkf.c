/*
   Copyright (C) 2011  Statoil ASA, Norway.

   The file 'ies_enkf.c' is part of ERT - Ensemble based Reservoir Tool.

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


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <time.h>

#include <ert/util/util.hpp>
#include <ert/util/type_macros.hpp>
#include <ert/util/rng.hpp>
#include <ert/util/bool_vector.hpp>

#include <ert/res_util/matrix.hpp>
#include <ert/res_util/matrix_blas.hpp>


#include <ert/analysis/analysis_module.hpp>
#include <ert/analysis/analysis_table.hpp>
#include <ert/analysis/enkf_linalg.hpp>
#include <ert/analysis/std_enkf.hpp>

#include <ies_enkf_config.h>

typedef struct ies_enkf_data_struct ies_enkf_data_type;

#define IES_ENKF_TYPE_ID 19640202

#define ENKF_SUBSPACE_DIMENSION_KEY      "ENKF_SUBSPACE_DIMENSION"
#define ENKF_TRUNCATION_KEY              "ENKF_TRUNCATION"
#define IES_STEPLENGTH_KEY               "IES_STEPLENGTH"
#define GAUSS_NEWTON_CONV_KEY            "GAUSS_NEWTON_CONV"
#define ITER_KEY                         "ITER"

#define IES_SUBSPACE_KEY                 "IES_SUBSPACE"
#define IES_INVERSION_KEY                "IES_INVERSION"
#define IES_LOGFILE_KEY                  "IES_LOGFILE"
#define IES_DEBUG_KEY                    "IES_DEBUG"

//#define DEFAULT_ANALYSIS_SCALE_DATA true


//**********************************************
// IES "object" data definition
//**********************************************
/*
  The configuration data used by the ies_enkf module is contained in a
  ies_enkf_data_struct instance. The data type used for the ies_enkf
  module is quite simple; with only a few scalar variables, but there
  are essentially no limits to what you can pack into such a datatype.

  All the functions in the module have a void pointer as the first
  argument, this will immediately be casted to an ies_enkf_data_type
  instance, to get some type safety the UTIL_TYPE_ID system should be
  used.

  The data structure holding the data for your analysis module should
  be created and initialized by a constructor, which should be
  registered with the '.alloc' element of the analysis table; in the
  same manner the desctruction of this data should be handled by a
  destructor or free() function registered with the .freef field of
  the analysis table.
*/

struct ies_enkf_data_struct {
   UTIL_TYPE_ID_DECLARATION;
   int       iteration_nr;            // Keep track of the outer iteration loop
   double    truncation;              // Controlled by config key: ENKF_TRUNCATION_KEY
   int       subspace_dimension;      // Controlled by config key: ENKF_SUBSPACE_DIMENSION_KEY (-1: use Truncation instead)
   double    ies_steplength;          // Step length in Gauss Newton iteration
   double    gauss_newton_conv;       // NOT USED
   bool      ies_subspace;            // NOT USED
   int       ies_inversion;                // 1 for exact R, 2 for R=EE', 3 for E
   char    * ies_logfile;             // Filename for ies logfile
   bool      ies_debug;               // True or false for storing debug information
   int       max_gauss_newton_it;     // NOT USED
   int state_size;                    // Initial state_size used for checks in subsequent calls
   bool_vector_type * ens_mask0;      // Initial ensemble mask of active realizations
   bool_vector_type * ens_mask;       // Ensemble mask of active realizations
   bool_vector_type * obs_mask0;      // Initial observation mask for active measurements
   bool_vector_type * obs_mask;       // Current observation mask
   matrix_type * W;                   // Coefficient matrix used to compute Omega = I + W (I -11'/N)/sqrt(N-1)
   matrix_type * A0;                  // Prior ensemble used in Ei=A0 Omega_i
   matrix_type * E;                   // Prior ensemble of measurement perturations (should be the same for all iterations)
   bool      converged;               // GN has converged
   ies_enkf_config_type * config;     // This I don't understand but I assume I include data from the ies_enkf_config_type defined in ies_enkf_config.c
   bool      AAprojection;            // For including the AAprojection of Y
};

static UTIL_SAFE_CAST_FUNCTION( ies_enkf_data , IES_ENKF_TYPE_ID )
static UTIL_SAFE_CAST_FUNCTION_CONST( ies_enkf_data , IES_ENKF_TYPE_ID )

void * ies_enkf_data_alloc( rng_type * rng) {
  ies_enkf_data_type * data = util_malloc( sizeof * data);
  UTIL_TYPE_ID_INIT( data , IES_ENKF_TYPE_ID );
  data->iteration_nr         = 0;
  data->truncation           = 0.0;
  data->subspace_dimension   = 0;
  data->ies_steplength       = 0.0;
  data->gauss_newton_conv    = 0.0;
  data->ies_subspace         = false;
  data->ies_inversion        = false;
  data->ies_logfile          = NULL ;
  data->ies_debug            = false;
  data->max_gauss_newton_it  = 0;
  data->state_size           = 0;
  data->ens_mask0            = NULL;
  data->ens_mask             = NULL;
  data->obs_mask0            = NULL;
  data->obs_mask             = NULL;
  data->W                    = NULL;
  data->A0                   = NULL;
  data->E                    = NULL;
  data->converged            = false;
  data->config               = ies_enkf_config_alloc();
  data->AAprojection         = true;
  return data;
}

void ies_enkf_data_free( void * arg ) {
  ies_enkf_data_type * data = ies_enkf_data_safe_cast( arg );
  ies_enkf_config_free( data->config );
  free( data );
}

void teccost(matrix_type * W,
             matrix_type * D,
             const char * fname,
             int ens_size,
             int izone)
{
   FILE *fp2;
   float costJ[100];

   if (izone == 1){
      fp2 = fopen(fname, "w");
      fprintf(fp2, "TITLE = \"%s\"\n",fname);
      fprintf(fp2, "VARIABLES = \"i\" \"Average J\" ");
      for (int i = 0; i < ens_size; i++)
         fprintf(fp2, "\"%d\" ",i);
      fprintf(fp2,"\n");
   } else {
      fp2 = fopen(fname, "a");
   }

   float costf=0.0;
   for (int i = 0; i < ens_size; i++){
      costJ[i]=matrix_column_column_dot_product(W,i,W,i)+matrix_column_column_dot_product(D,i,D,i);
      costf += costJ[i];
   }
   costf= costf/ens_size;

   fprintf(fp2,"%2d %10.3f ", izone-1, costf) ;
   for (int i = 0; i < ens_size; i++)
      fprintf(fp2,"%10.3f ", costJ[i] ) ;
   fprintf(fp2,"\n") ;
   fclose(fp2);
}

void teclog(matrix_type * W,
            matrix_type * D,
            matrix_type * DW,
            const char * fname,
            int ens_size,
            int itnr,
            double rcond,
            int nrsing,
            int nrobs)
{

   double diff1W = 0.0;
   double diff2W = 0.0;
   double diffW;
   double costfp,costfd,costf;
   FILE *fp;

   if (itnr == 1){
      fp = fopen(fname, "w");
      fprintf(fp, "TITLE = \"%s\"\n",fname);
      fprintf(fp, "VARIABLES = \"it\" \"Diff1W\" \"Diff2W\" \"J-prior\" \"J-data\" \"J\" \"rcond\" \"Sing\" \"nrobs\"\n");
   } else {
      fp = fopen(fname, "a");
   }

   for (int i = 0; i < ens_size; i++){
       diffW=matrix_column_column_dot_product(DW,i,DW,i)/ens_size;
       diff2W+=diffW;
       if (diffW > diff1W) diff1W=diffW ;
   }
   diff2W=diff2W/ens_size;

   costfp=0.0;
   costfd=0.0;
   for (int i = 0; i < ens_size; i++){
       costfp += matrix_column_column_dot_product(W,i,W,i);
       costfd += matrix_column_column_dot_product(D,i,D,i);
   }
   costfp=costfp/ens_size;
   costfd=costfd/ens_size;
   costf= costfp + costfd;
   fprintf(fp," %d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %d %d\n", itnr-1, diff1W, diff2W, costfp, costfd, costf, rcond, nrsing , nrobs ) ;
   fclose(fp);
}

void tecfld(matrix_type * W,
           const char * fname,
           const char * vname,
           int ii,
           int jj,
           int izone)
{
   FILE *fp2;
   int ic;

   if (izone == 1){
      fp2 = fopen(fname, "w");
      fprintf(fp2, "TITLE = \"%s\"\n",fname);
      fprintf(fp2, "VARIABLES = \"i-index\" \"j-index\" \"%s\"\n",vname);
      fprintf(fp2, "ZONE F=BLOCK, T=\"Iteration %d\", I=%d, J=%d, K=1\n",izone, ii, jj);
      ic=0;
      for (int j=0;j<ii;j++){
         for (int i=0;i<jj;i++){
            ic=ic+1;
            fprintf(fp2,"%d ",i);
            if (ic==30){
               ic=0;
               fprintf(fp2,"\n");
            }
         }
      }

      ic=0;
      for (int j=0;j<ii;j++){
         for (int i=0;i<jj;i++){
            ic=ic+1;
            fprintf(fp2,"%d ",j);
            if (ic==30){
               ic=0;
               fprintf(fp2,"\n");
            }
         }
      }
   } else {
      fp2 = fopen(fname, "a");
      fprintf(fp2, "ZONE F=BLOCK, T=\"Iteration %d\", I=%d, J=%d, K=1, VARSHARELIST = ([1,2]=1)",izone, ii, jj);
   }

   ic=0;
   for (int j=0;j<ii;j++){
      for (int i=0;i<jj;i++){
         ic=ic+1;
         fprintf(fp2,"%12.5e ",matrix_iget(W , i,j));
         if (ic==30){
            ic=0;
            fprintf(fp2,"\n");
         }
      }
   }
   fclose(fp2);
}



/***************************************************************************************************************
*  Set / Get iteration number
****************************************************************************************************************/
void ies_enkf_set_iteration_nr( ies_enkf_data_type * data , int iteration_nr) {
  data->iteration_nr = iteration_nr;
}

int ies_enkf_get_iteration_nr( const ies_enkf_data_type * data ) {
  return data->iteration_nr;
}

/*
 * -------------------------------------------------------------------------------------------------------------
 * I E n K S
 * IEnKS initialization (getting the obs_mask and ens_mask for active observations and realizations)
 * -------------------------------------------------------------------------------------------------------------
*/
void ies_enkf_init_update(void * arg ,
                          const bool_vector_type * ens_mask ,
                          const bool_vector_type * obs_mask ,
                          const matrix_type * S ,
                          const matrix_type * R ,
                          const matrix_type * dObs ,
                          const matrix_type * E ,
                          const matrix_type * D,
                          rng_type * rng) {
  ies_enkf_data_type * module_data = ies_enkf_data_safe_cast( arg );

/* Store current ens_mask in module_data->ens_mask for each iteration */
  if (module_data->ens_mask)
    bool_vector_free( module_data->ens_mask );

  module_data->ens_mask = bool_vector_alloc_copy( ens_mask );

/* Store obs_mask for initial iteration in module_data->obs_mask0,
*  for each subsequent iteration we store the current mask in module_data->obs_mask */
  if (!module_data->obs_mask0)
     module_data->obs_mask0 = bool_vector_alloc_copy( obs_mask );

  if (module_data->obs_mask)
    bool_vector_free( module_data->obs_mask );

  module_data->obs_mask = bool_vector_alloc_copy( obs_mask );
}



/*
 * -------------------------------------------------------------------------------------------------------------
 * I E n K S
 * IEnKS update (IES searching the solution in ensemble subspace)
 * -------------------------------------------------------------------------------------------------------------
*/

void ies_enkf_updateA( void * module_data,
                       matrix_type * A ,      // Updated ensemble A retured to ERT.
                       matrix_type * Yin ,    // Ensemble of predicted measurements
                       matrix_type * Rin ,    // Measurement error covariance matrix (not used)
                       matrix_type * dObs ,   // Actual observations (not used)
                       matrix_type * Ein ,    // Ensemble of observation perturbations
                       matrix_type * Din ,    // (d+E-Y) Ensemble of perturbed observations - Y
                       const module_info_type * module_info,
                       rng_type * rng) {

   ies_enkf_data_type * data = ies_enkf_data_safe_cast( module_data );

   int nrobs_msk     =bool_vector_size( data->obs_mask ); // Total number of observations
   int nrobs_inp     =matrix_get_rows( Yin );             // Number of active observations input in current iteration
   int nrobs         =nrobs_inp;                          // Number of selected active observations


   int ens_size_msk  =bool_vector_size(data->ens_mask);   // Total number of realizations
   int ens_size      =matrix_get_columns( Yin );          // Number of active realizations in current iteration

   int state_size    = matrix_get_rows( A );

   double rcond;
   int nrsing=0;

   if (data->state_size == 0) data->state_size=state_size;


   data->iteration_nr         = data->iteration_nr+1 ;
   data->ies_steplength       = ies_enkf_config_get_ies_steplength( data->config );
   data->ies_subspace         = ies_enkf_config_get_ies_subspace( data->config );
   data->ies_inversion        = ies_enkf_config_get_ies_inversion( data->config );
   data->ies_logfile          = ies_enkf_config_get_ies_logfile( data->config );
   data->ies_debug            = ies_enkf_config_get_ies_debug( data->config );

   data->truncation           = ies_enkf_config_get_truncation( data->config );
   data->subspace_dimension   = ies_enkf_config_get_enkf_subspace_dimension( data->config );

   data->gauss_newton_conv    = ies_enkf_config_get_gauss_newton_conv( data->config );

   if (data->iteration_nr > 3) data->ies_steplength=data->ies_steplength/2;
   if (data->iteration_nr > 6) data->ies_steplength=data->ies_steplength/4;


   FILE *fp;
   if (data->iteration_nr == 1){
      fp=freopen(data->ies_logfile, "w", stdout);
   } else {
      fp=freopen(data->ies_logfile, "a", stdout);
   }

   fprintf(stdout,"\n\n\n***********************************************************************\n");
   fprintf(stdout,"IES Iteration   = %d\n", data->iteration_nr);
   fprintf(stdout,"----ies_steplength  = %f\n", data->ies_steplength);
   fprintf(stdout,"----ies_inversion   = %d\n", data->ies_inversion);
   fprintf(stdout,"----ies_debug       = %d\n", data->ies_debug);
   fprintf(stdout,"----truncation      = %f %d\n", data->truncation, data->subspace_dimension);
   fprintf(stdout,"----ies_logfile     = %s\n", data->ies_logfile);
   bool dbg = ies_enkf_config_get_ies_debug( data->config ) ;

/***************************************************************************************************************/
/* Counting number of active observations for current iteration. The number requires that
   the observations were included in the initial call for storage in data->E as well as in
   in the current call. Thus, it is possible to remove observations but not include new ones. */
   nrobs=0;
   for (int i = 0; i < nrobs_msk; i++){
      if ( bool_vector_iget(data->obs_mask0,i) && bool_vector_iget(data->obs_mask,i) ){
         nrobs=nrobs+1;
      }
   }

   int nrmin  = util_int_min( ens_size , nrobs);
   double nsc = 1.0/sqrt(ens_size - 1.0);

/* dimensions for printing */
   int m_nrobs     =util_int_min(nrobs     -1,7);
   int m_ens_size  =util_int_min(ens_size  -1,16);
   int m_state_size=util_int_min(state_size-1,3);

/***************************************************************************************************************/
   if (!data->E){
      // We store the initial observation perturbations corresponding to data->obs_mask0.
      fprintf(stdout,"Allocating and assigning data->E \n");
      data->E=matrix_alloc( nrobs , ens_size  );
      matrix_assign(data->E,Ein);
      if (dbg) matrix_pretty_fprint_submat(data->E,"data->E","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;
   }

   if (!data->W){
      // We initialize data-W which will store W for use in next iteration                    (Line 9)
      fprintf(stdout,"Allocating data->W\n");
      data->W=matrix_alloc( ens_size , ens_size  );
      matrix_set(data->W , 0.0) ;
      if (dbg) matrix_pretty_fprint_submat(data->W,"Ini data->W","%11.5f",stdout,0,m_ens_size,0,m_ens_size) ;
   }

   if (!data->A0){
      // We store the initial ensemble to use it in final update equation                     (Line 11)
      fprintf(stdout,"Allocating and assigning data->A0 \n");
      data->A0=matrix_alloc( state_size , ens_size  );
      matrix_assign(data->A0,A);
      if (dbg) matrix_pretty_fprint_submat(data->A0,"Ini data->A0","%11.5f",stdout,0,m_state_size,0,m_ens_size) ;
   }



/***************************************************************************************************************/
   fprintf(stdout,"----active ens_size  = %d, total ens_size_msk = %d\n", ens_size,ens_size_msk);
   fprintf(stdout,"----active nrobs     = %d, nrobs_inp= %d, total nrobs_msk= %d\n", nrobs, nrobs_inp, nrobs_msk );
   fprintf(stdout,"----active state_size= %d\n", state_size );

/***************************************************************************************************************/
/* Print initial observation mask */
   fprintf(stdout,"obsmask_0:") ;
   for (int i = 0; i < nrobs_msk; i++){
      fprintf(stdout,"%d",bool_vector_iget(data->obs_mask0,i));
      if ((i+1)%10 == 0) fprintf(stdout," ") ;
      if ((i+1)%100 == 0) fprintf(stdout,"\nobsmask_0:") ;
   }
   fprintf(stdout,"\n");

/***************************************************************************************************************/
/* Print Current observation mask */
   fprintf(stdout,"obsmask_i:") ;
   for (int i = 0; i < nrobs_msk; i++){
      fprintf(stdout,"%d",bool_vector_iget(data->obs_mask,i));
      if ((i+1)%10 == 0) fprintf(stdout," ") ;
      if ((i+1)%100 == 0) fprintf(stdout,"\nobsmask_i:") ;
   }
   fprintf(stdout,"\n");

/***************************************************************************************************************/
/* Print Current ensemble mask */
   fprintf(stdout,"ensmask_i:");
   for (int i = 0; i < ens_size_msk; i++){
      fprintf(stdout,"%d",bool_vector_iget(data->ens_mask,i));
      if ((i+1)%10 == 0) fprintf(stdout," ") ;
      if ((i+1)%100 == 0) fprintf(stdout,"\n") ;
   }
   fprintf(stdout,"\n\n");

/***************************************************************************************************************
* Re structure input matrices according to new active obs_mask and ens_size.
*     Allocates the local matrices to be used.
*     Copies the initial measurement perturbations for the active observations into the current E matrix.
*     Copies the inputs in D, Y and R into their local representations
*/
   matrix_type * Y   = matrix_alloc( nrobs    , ens_size );
   matrix_type * E   = matrix_alloc( nrobs    , ens_size );
   matrix_type * D   = matrix_alloc( nrobs    , ens_size );
   matrix_type * Rtmp= matrix_alloc( nrobs    , nrobs_inp );
   matrix_type * R   = matrix_alloc( nrobs    , nrobs );

/* Subtract new measurement perturbations              D=D-E    */
   matrix_inplace_sub(Din,Ein);

/* E=data->E but only using the active obs also stored in data->E */
   int j=-1;  // counter for initial mask0
   int k=-1;  // counter for current mask
   int m=-1;  // counter for currently active measurements
   for (int iobs = 0; iobs < nrobs_msk; iobs++){
      if ( bool_vector_iget(data->obs_mask0,iobs) )
         j=j+1 ;

      if ( bool_vector_iget(data->obs_mask,iobs) )
         k=k+1 ;

      if ( bool_vector_iget(data->obs_mask0,iobs) && bool_vector_iget(data->obs_mask,iobs) ){
         m=m+1;

         {
            int i=-1;
            for (int iens = 0; iens < ens_size_msk; iens++){
               if ( bool_vector_iget(data->ens_mask,iens) ){
                  i=i+1 ;
                  matrix_iset_safe(E,m,i,matrix_iget(data->E,j,iens)) ;
               }
            }
         }

         matrix_copy_row(D,Din,m,k);
         matrix_copy_row(Y,Yin,m,k);
         matrix_copy_row(Rtmp,Rin,m,k);
         matrix_copy_column(R,Rtmp,m,k);
      }
   }
   if (ens_size_msk == ens_size && nrobs == nrobs_inp){
      fprintf(stdout,"data->E copied exactly to E: %d\n",matrix_equal(data->E,E)) ;
   }

   fprintf(stdout,"Input matrices\n");
   if (dbg) matrix_pretty_fprint_submat(E,"E","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;

   if (dbg) matrix_pretty_fprint_submat(Din,"Din","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;
   if (dbg) matrix_pretty_fprint_submat(D,"D","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;

   if (dbg) matrix_pretty_fprint_submat(Yin,"Yin","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;
   if (dbg) matrix_pretty_fprint_submat(Y,"Y","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;

   if (dbg) matrix_pretty_fprint_submat(Rin,"Rin","%11.5f",stdout,0,m_nrobs,0,m_nrobs) ;
   if (dbg) matrix_pretty_fprint_submat(R,"R","%11.5f",stdout,0,m_nrobs,0,m_nrobs) ;

   matrix_inplace_add(D,E);          // Add old measurement perturbations



   matrix_type * A0  = matrix_alloc( state_size, ens_size );  // Temporary ensemble matrix
   matrix_type * W0  = matrix_alloc( ens_size , ens_size  );  // Coefficient matrix
   matrix_type * W   = matrix_alloc( ens_size , ens_size  );  // Coefficient matrix
   matrix_type * DW  = matrix_alloc( ens_size , ens_size  );  // Coefficient matrix W - data->W
   matrix_type * H   = matrix_alloc( nrobs    , ens_size  );  // Innovation vector "H= S*W+D-Y"
   matrix_type * S   = matrix_alloc( nrobs    , ens_size );   // Predicted ensemble anomalies scaled with inv(Omeaga)

   matrix_type * YT  = matrix_alloc( ens_size, nrobs     );   // Y^T used in linear solver
   matrix_type * ST  = matrix_alloc( ens_size, nrobs     );   // current S^T used in linear solver
   matrix_type * STO = matrix_alloc( ens_size, nrobs     );   // previous S^T used in linear solver
   matrix_type * SD  = matrix_alloc( ens_size, nrobs     );   // difference between ST and STO in linear solver
   matrix_type * X   = matrix_alloc( ens_size, ens_size  );   // Used for Omega and transform matrix

   double      * eig = (double*)util_calloc( ens_size , sizeof * eig);


   if (dbg) matrix_pretty_fprint_submat(A,"Ain","%11.5f",stdout,0,m_state_size,0,m_ens_size) ;
   if (dbg) fprintf(stdout,"Computed matrices\n");

/***************************************************************************************************************
*  Subtract mean of predictions to generate predicted ensemble anomaly matrix                 (Line 5)
*/
   matrix_subtract_row_mean( Y );   // Y=Y*(I-(1/ens_size)*11)
   matrix_scale(Y,nsc);             // Y=Y / sqrt(ens_size-1)
   if (dbg) matrix_pretty_fprint_submat(Y,"Y","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;


/***************************************************************************************************************
*  COMPUTING THE PROJECTION Y= Y * (Ai^+ * Ai) (only used when state_size < ens_size-1)    */
   if (data->AAprojection && state_size <= ens_size -1){
      fprintf(stdout,"Activating AAi projection for Y\n");
      matrix_type * Ai    = matrix_alloc_copy( A );
      matrix_type * AAi   = matrix_alloc( ens_size, ens_size  );
      matrix_subtract_row_mean(Ai);
      matrix_type * VT    = matrix_alloc( state_size, ens_size  );
      matrix_dgesvd(DGESVD_NONE , DGESVD_MIN_RETURN , Ai , eig , NULL , VT);
      if (dbg) matrix_pretty_fprint_submat(VT,"VT","%11.5f",stdout,0,state_size-1,0,m_ens_size) ;
      matrix_dgemm(AAi,VT,VT,true,false,1.0,0.0);
      if (dbg) matrix_pretty_fprint_submat(AAi,"AAi","%11.5f",stdout,0,m_ens_size,0,m_ens_size) ;
      matrix_inplace_matmul(Y,AAi);
      matrix_free(Ai);
      matrix_free(AAi);
      matrix_free(VT);
      if (dbg) matrix_pretty_fprint_submat(Y,"Yprojected","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;
   }

/***************************************************************************************************************
*  COPY ACTIVE REALIZATIONS FROM data->W to W0 */
   {
      int i=-1;
      int j;
      for (int iens=0; iens < ens_size_msk; iens++){
         if ( bool_vector_iget(data->ens_mask,iens) ){
            i=i+1;
            j=-1;
            for (int jens=0; jens < ens_size_msk; jens++){
               if ( bool_vector_iget(data->ens_mask,jens) ){
                  j=j+1;
                  matrix_iset_safe(W0,i,j,matrix_iget(data->W,iens,jens)) ;
               }
            }
         }
      }
   }

   if (ens_size_msk == ens_size){
      fprintf(stdout,"data->W copied exactly to W0: %d\n",matrix_equal(data->W,W0)) ;
   }
   if (dbg) matrix_pretty_fprint_submat(data->W,"data->W","%11.5f",stdout,0,m_ens_size,0,m_ens_size) ;
   if (dbg) matrix_pretty_fprint_submat(W0,"W0","%11.5f",stdout,0,m_ens_size,0,m_ens_size) ;


/***************************************************************************************************************
* COMPUTE  X= I + W (I-11'/sqrt(ens_size))    from Eq. (36).                                   (Line 6)
*  When solving the system S = Y inv(Omega) we write
*     X^T S^T = Y^T
*  Here we compute the W (I-11'/N) / sqrt(N-1)  and transpose it).
*/

   matrix_assign(X,W0) ;            // X=data->W (from previous iteration used to solve for S)
   matrix_subtract_row_mean(X);     // X=X*(I-(1/N)*11')
   matrix_scale(X,nsc);             // X/sqrt(N-1)
   matrix_inplace_transpose(X);     // X=transpose(X)
   for (int i = 0; i < ens_size; i++){  // X=X+I
      matrix_iadd(X,i,i,1.0);
   }
   if (dbg) matrix_pretty_fprint_submat(X,"OmegaT","%11.5f",stdout,0,m_ens_size,0,m_ens_size) ;
   if (dbg) tecfld( X, "tecOmega.dat" , "Omega", ens_size, ens_size , data->iteration_nr);
   matrix_transpose(Y,YT);         // RHS stored in YT

/* Solve system and return S in YT                                                             (Line 7)   */
   fprintf(stdout,"Solving X' S' = Y' using LU factorization:\n");
   matrix_dgesvx(X,YT,&rcond);
   fprintf(stdout,"dgesvx condition number= %12.5e\n",rcond);

   matrix_transpose(YT,S);          // Copy solution to S


   if (data->iteration_nr == 1){
         fprintf(stdout,"dgesvx: Y exactly equal to S: %d\n",matrix_equal(Y,S)) ;
   }


   if (dbg) matrix_pretty_fprint_submat(S,"S","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;

/***************************************************************************************************************
*  INNOVATION H = S*W + D - Y   from Eq. (47)                                                  (Line 8)    */
   matrix_assign(H,D) ;                            // H=D=dobs + E - Y
   matrix_dgemm(H,S,W0,false,false,1.0,1.0);       // H=S*W + H

   if (dbg) matrix_pretty_fprint_submat(H,"H","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;

/* Store previous W for convergence test */
   matrix_assign(W,W0);

/***************************************************************************************************************
* COMPUTE NEW UPDATED W                                                                        (Line 9)
*  We first compute the expression
*          S'*(S*S'+R)^{-1} H           (a)
*  which in the case when R=I can be rewritten as
*          (S'*S + I)^{-1} * S' * H     (b)
*
*  With R=I the subspace inversion (ies_inversion=1) solving Eq. (a) with singular value
*  trucation=1.000 gives exactly the same solution as the exact inversion (ies_inversion=0).
*
*  Using ies_inversion=1, and a step length of 1.0, one update gives identical result to ENKF_STD
*  as long as the same SVD truncation is used.
*
*  With very large data sets it is likely that the inversion becomes poorly conditioned
*  and a trucation=1.000 is not a good choice. In this case the ies_inversion > 0 and EnKF_truncation set
*  to 0.999 or so, should stabelize the algorithm.
*
*  Using ies_inversion=2 and ies_inversion=3 gives identical results but ies_inversion=3 is much faster (N^2m)
*  than ies_inversion=2 (Nm^2).

   ies_inversion=0  -> exact inversion from (b) with exact R=I
   ies_inversion=1  -> subspace inversion from (a) with exact R
   ies_inversion=2  -> subspace inversion from (a) with R=EE
   ies_inversion=3  -> subspace inversion from (a) with R represented by E
*/

   if (data->ies_inversion > 0){
      fprintf(stdout,"Subspace inversion. (ies_inversion=%d)\n",data->ies_inversion);
      matrix_type * X1  = matrix_alloc( nrobs   , nrmin     );   // Used in subspace inversion
      matrix_type * X3  = matrix_alloc( nrobs   , ens_size  );   // Used in subspace inversion
      if (data->ies_inversion == 3){
         fprintf(stdout,"Subspace inversion using E to represent errors. (ies_inversion=%d)\n",data->ies_inversion);
         matrix_scale(E,nsc);
         enkf_linalg_lowrankE( S , E , X1 , eig , data->truncation , data->subspace_dimension);
      } else if (data->ies_inversion == 2){
         fprintf(stdout,"Subspace inversion using ensemble generated full R=EE. (ies_inversion=%d)'\n",data->ies_inversion);
         matrix_scale(E,nsc);
         matrix_type * Et = matrix_alloc_transpose( E );
         matrix_type * Cee = matrix_alloc_matmul( E , Et );
         matrix_scale(Cee,nsc*nsc); // since enkf_linalg_lowrankCinv solves (SS' + (N-1) R)^{-1}
         if (dbg) matrix_pretty_fprint_submat(Cee,"Cee","%11.5f",stdout,0,m_nrobs,0,m_nrobs) ;
         enkf_linalg_lowrankCinv( S , Cee , X1 , eig , data->truncation , data->subspace_dimension);
         matrix_free( Et );
         matrix_free( Cee );
      } else if (data->ies_inversion == 1){
         fprintf(stdout,"Subspace inversion using 'exact' full R. (ies_inversion=%d)\n",data->ies_inversion);
         matrix_scale(R,nsc*nsc); // since enkf_linalg_lowrankCinv solves (SS' + (N-1) R)^{-1}
         if (dbg) matrix_pretty_fprint_submat(R,"R","%11.5f",stdout,0,m_nrobs,0,m_nrobs) ;
         enkf_linalg_lowrankCinv( S , R , X1 , eig , data->truncation , data->subspace_dimension);
      }

      nrsing=0;
      fprintf(stdout,"\nEig:\n");
      for (int i=0;i<nrmin;i++){
         fprintf(stdout," %f ", eig[i]);
         if ((i+1)%20 == 0) fprintf(stdout,"\n") ;
         if (eig[i] < 1.0) nrsing+=1;
      }
      fprintf(stdout,"\n");

/*    X3 = X1 * diag(eig) * X1' * H (Similar to Eq. 14.31, Evensen (2007))                                  */
      enkf_linalg_genX3(X3 , X1 , H , eig);

      if (dbg) matrix_pretty_fprint_submat(X1,"X1","%11.5f",stdout,0,m_nrobs,0,util_int_min(m_nrobs,nrmin-1)) ;
      if (dbg) matrix_pretty_fprint_submat(X3,"X3","%11.5f",stdout,0,m_nrobs,0,m_ens_size) ;

/*    Update data->W = (1-ies_steplength) * data->W +  ies_steplength * S' * X3                          (Line 9)    */
      matrix_dgemm(W0 , S , X3 , true , false , data->ies_steplength , 1.0-data->ies_steplength);

      matrix_free( X1 );
      matrix_free( X3 );

   } else if (data->ies_inversion == 0) {
      fprintf(stdout,"Exact inversion using diagonal R=I. (ies_inversion=%d)\n",data->ies_inversion);
      matrix_type * Z      = matrix_alloc( ens_size , ens_size  );  // Eigen vectors of S'S+I
      matrix_type * StH    = matrix_alloc( ens_size , ens_size );
      matrix_type * StS    = matrix_alloc( ens_size , ens_size );
      matrix_type * ZtStH  = matrix_alloc( ens_size , ens_size );

      matrix_diag_set_scalar(StS,1.0);
      matrix_dgemm(StS,S,S,true,false,1.0,1.0);
      matrix_dgesvd(DGESVD_ALL , DGESVD_NONE , StS , eig , Z , NULL);

      matrix_dgemm(StH,S,H,true,false,1.0,0.0);
      matrix_dgemm(ZtStH,Z,StH,true,false,1.0,0.0);

      for (int i=0;i<ens_size;i++){
         eig[i]=1.0/eig[i] ;
         matrix_scale_row(ZtStH,i,eig[i]);
      }

      fprintf(stdout,"\nEig:\n");
      for (int i=0;i<ens_size;i++){
         fprintf(stdout," %f ", eig[i]);
         if ((i+1)%20 == 0) fprintf(stdout,"\n") ;
      }
      fprintf(stdout,"\n");

/*    Update data->W = (1-ies_steplength) * data->W +  ies_steplength * Z * (Lamda^{-1}) Z' S' H         (Line 9)    */
      matrix_dgemm(W0 , Z , ZtStH , false , false , data->ies_steplength , 1.0-data->ies_steplength);

      matrix_free(Z);
      matrix_free(StH);
      matrix_free(StS);
      matrix_free(ZtStH);
   }

   if (dbg) matrix_pretty_fprint_submat(W0,"Updated W","%11.5f",stdout,0,m_ens_size,0,m_ens_size) ;



/* Store active realizations from W0 to data->W */
   {
      int i=-1;
      int j;
      matrix_set(data->W , 0.0) ;
      for (int iens=0; iens < ens_size_msk; iens++){
         if ( bool_vector_iget(data->ens_mask,iens) ){
            i=i+1;
            j=-1;
            for (int jens=0; jens < ens_size_msk; jens++){
               if ( bool_vector_iget(data->ens_mask,jens) ){
                  j=j+1;
                  matrix_iset_safe(data->W,iens,jens,matrix_iget(W0,i,j)) ;
               }
            }
         }
      }
   }
   if (ens_size_msk == ens_size){
      fprintf(stdout,"W0 copied exactly to data->W: %d\n",matrix_equal(data->W,W0)) ;
   }

   if (dbg) tecfld( W, "tecW.dat" , "W", ens_size, ens_size , data->iteration_nr);


/***************************************************************************************************************
*  CONSTRUCT TRANFORM MATRIX X FOR CURRENT ITERATION                                         (Line 10)
*     X= I + W/sqrt(N-1)          */
   matrix_assign(X,W0);
   matrix_scale(X,nsc);
   for (int i = 0; i < ens_size; i++){
      matrix_iadd(X,i,i,1.0);
   }
   if (dbg) matrix_pretty_fprint_submat(X,"X","%11.5f",stdout,0,m_ens_size,0,m_ens_size) ;

/***************************************************************************************************************
*  COMPUTE NEW ENSEMBLE SOLUTION FOR CURRENT ITERATION  Ei=A0*X                              (Line 11)   */
   matrix_pretty_fprint_submat(data->A0,"data->A0","%11.5f",stdout,0,m_state_size,0,m_ens_size) ;
   matrix_pretty_fprint_submat(A,"A^f","%11.5f",stdout,0,m_state_size,0,m_ens_size) ;

   {
      int i=-1;
      for (int iens=0; iens < ens_size_msk; iens++){
         if ( bool_vector_iget(data->ens_mask,iens) ){
            i=i+1;
            matrix_copy_column(A0,data->A0,i,iens);
         }
      }
   }
   if (ens_size_msk == ens_size){
      fprintf(stdout,"data->A0 copied exactly to A0: %d\n",matrix_equal(data->A0,A0)) ;
   }

   if (dbg) tecfld( X, "tecX.dat" , "X", ens_size, ens_size , data->iteration_nr);
   matrix_matmul(A,A0,X);
   matrix_pretty_fprint_submat(A,"A^a","%11.5f",stdout,0,m_state_size,0,m_ens_size) ;


/***************************************************************************************************************
*  COMPUTE ||W0 - W|| AND EVALUATE COST FUNCTION FOR PREVIOUS ITERATE                        (Line 12)   */
   matrix_sub(DW,W0,W);
   if (dbg) tecfld( DW, "tecDW.dat" , "Delta W", ens_size, ens_size , data->iteration_nr);
   teclog(W,D,DW,"iesteclog.dat",ens_size,data->iteration_nr, rcond, nrsing, nrobs);
   teccost(W,D,"costf.dat",ens_size,data->iteration_nr);

/* DONE *********************************************************************************************************/


   fflush(stdout);
   fclose(fp);

   matrix_free( Y  );
   matrix_free( D  );
   matrix_free( E  );
   matrix_free( Rtmp);
   matrix_free( R  );
   matrix_free( A0 );
   matrix_free( W0 );
   matrix_free( W );
   matrix_free( DW );
   matrix_free( H  );
   matrix_free( S  );
   matrix_free( YT );
   matrix_free( ST );
   matrix_free( STO);
   matrix_free( SD );
   matrix_free( X  );
}





//**********************************************
// Set / Get basic types
//**********************************************
bool ies_enkf_set_int( void * arg , const char * var_name , int value) {
  ies_enkf_data_type * module_data = ies_enkf_data_safe_cast( arg );
  {
    bool name_recognized = true;

    if (strcmp( var_name , ENKF_SUBSPACE_DIMENSION_KEY) == 0)
      ies_enkf_config_set_enkf_subspace_dimension(module_data->config , value);
    else if (strcmp( var_name , ITER_KEY) == 0)
      ies_enkf_set_iteration_nr( module_data , value );
    else if (strcmp( var_name , IES_INVERSION_KEY) == 0)
      ies_enkf_config_set_ies_inversion( module_data->config , value );
    else
      name_recognized = false;

    return name_recognized;
  }
}

int ies_enkf_get_int( const void * arg, const char * var_name) {
  const ies_enkf_data_type * module_data = ies_enkf_data_safe_cast_const( arg );
  {
    if (strcmp(var_name , ITER_KEY) == 0)
      return ies_enkf_get_iteration_nr(module_data);
    else if (strcmp(var_name , ENKF_SUBSPACE_DIMENSION_KEY) == 0)
      return ies_enkf_config_get_enkf_subspace_dimension(module_data->config);
    else if (strcmp(var_name , IES_INVERSION_KEY) == 0)
      return ies_enkf_config_get_ies_inversion(module_data->config);
    else
      return -1;
  }
}

bool ies_enkf_set_string( void * arg , const char * var_name , const char * value) {
  ies_enkf_data_type * module_data = ies_enkf_data_safe_cast( arg );
  {
    bool name_recognized = true;

    if (strcmp( var_name , IES_LOGFILE_KEY) == 0)
      ies_enkf_config_set_ies_logfile( module_data->config , value );
    else
      name_recognized = false;

    return name_recognized;
  }
}

bool ies_enkf_get_string( void * arg , const char * var_name ) {
  ies_enkf_data_type * module_data = ies_enkf_data_safe_cast( arg );
  {
    if (strcmp( var_name , IES_LOGFILE_KEY) == 0)
      return ies_enkf_config_get_ies_logfile( module_data->config );

    else
       return false;
  }
}

bool ies_enkf_set_bool( void * arg , const char * var_name , bool value) {
  ies_enkf_data_type * module_data = ies_enkf_data_safe_cast( arg );
  {
    bool name_recognized = true;

    if (strcmp( var_name , IES_SUBSPACE_KEY) == 0)
      ies_enkf_config_set_ies_subspace( module_data->config , value);
    else if (strcmp( var_name , IES_DEBUG_KEY) == 0)
      ies_enkf_config_set_ies_debug( module_data->config , value );
    else
      name_recognized = false;

    return name_recognized;
  }
}

bool ies_enkf_get_bool( const void * arg, const char * var_name) {
  const ies_enkf_data_type * module_data = ies_enkf_data_safe_cast_const( arg );
  {
    if (strcmp(var_name , IES_SUBSPACE_KEY) == 0)
      return ies_enkf_config_get_ies_subspace( module_data->config );
    else if (strcmp(var_name , IES_DEBUG_KEY) == 0)
      return ies_enkf_config_get_ies_debug( module_data->config );
    else
       return false;
  }
}




bool ies_enkf_set_double( void * arg , const char * var_name , double value) {
  ies_enkf_data_type * module_data = ies_enkf_data_safe_cast( arg );
  {
    bool name_recognized = true;

    if (strcmp( var_name , ENKF_TRUNCATION_KEY) == 0)
      ies_enkf_config_set_truncation( module_data->config , value );
    else if (strcmp( var_name , IES_STEPLENGTH_KEY) == 0)
      ies_enkf_config_set_ies_steplength( module_data->config , value );
    else if (strcmp( var_name , GAUSS_NEWTON_CONV_KEY) == 0)
      ies_enkf_config_set_gauss_newton_conv( module_data->config , value );
    else
      name_recognized = false;

    return name_recognized;
  }
}

double ies_enkf_get_double( const void * arg, const char * var_name) {
  const ies_enkf_data_type * module_data = ies_enkf_data_safe_cast_const( arg );
  {
    if (strcmp(var_name , ENKF_TRUNCATION_KEY) == 0)
      return ies_enkf_config_get_truncation( module_data->config );

    if (strcmp(var_name , IES_STEPLENGTH_KEY) == 0)
      return ies_enkf_config_get_ies_steplength(module_data->config);

    if (strcmp(var_name , GAUSS_NEWTON_CONV_KEY) == 0)
      return ies_enkf_config_get_gauss_newton_conv(module_data->config);

    return -1;
  }
}


long ies_enkf_get_options( void * arg , long flag ) {
  ies_enkf_data_type * module_data = ies_enkf_data_safe_cast( arg );
  {
    return ies_enkf_config_get_option_flags( module_data->config );
  }
}

bool ies_enkf_has_var( const void * arg, const char * var_name) {
  {
    if (strcmp(var_name , ITER_KEY) == 0)
      return true;
    else if (strcmp(var_name , IES_STEPLENGTH_KEY) == 0)
      return true;
    else if (strcmp(var_name , GAUSS_NEWTON_CONV_KEY) == 0)
      return true;
    else if (strcmp(var_name , IES_SUBSPACE_KEY) == 0)
      return true;
    else if (strcmp(var_name , IES_INVERSION_KEY) == 0)
      return true;
    else if (strcmp(var_name , IES_LOGFILE_KEY) == 0)
      return true;
    else if (strcmp(var_name , IES_DEBUG_KEY) == 0)
      return true;
    else if (strcmp(var_name , ENKF_TRUNCATION_KEY) == 0)
      return true;
    else if (strcmp(var_name , ENKF_SUBSPACE_DIMENSION_KEY) == 0)
      return true;
    else
      return false;
  }
}

void * ies_enkf_get_ptr( const void * arg , const char * var_name ) {
  const ies_enkf_data_type * module_data = ies_enkf_data_safe_cast_const( arg );
  {
    if (strcmp(var_name , IES_LOGFILE_KEY) == 0)
      return (void *) ies_enkf_config_get_ies_logfile( module_data->config );
    else
      return NULL;
  }
}


//**********************************************
// Symbol table
//**********************************************
#ifdef INTERNAL_LINK
#define LINK_NAME IES_ENKF
#else
#define LINK_NAME EXTERNAL_MODULE_SYMBOL
#endif


analysis_table_type LINK_NAME = {
  .name            = "IES_ENKF",
  .initX           = NULL,
  .updateA         = ies_enkf_updateA,
  .init_update     = ies_enkf_init_update,
  .complete_update = NULL,
  .alloc           = ies_enkf_data_alloc,
  .freef           = ies_enkf_data_free,
  .has_var         = ies_enkf_has_var,
  .set_int         = ies_enkf_set_int ,
  .set_double      = ies_enkf_set_double ,
  .set_bool        = ies_enkf_set_bool ,
  .set_string      = ies_enkf_set_string ,
  .get_options     = ies_enkf_get_options ,
  .get_int         = ies_enkf_get_int,
  .get_double      = ies_enkf_get_double,
  .get_bool        = ies_enkf_get_bool ,
  .get_ptr         = ies_enkf_get_ptr ,
};

