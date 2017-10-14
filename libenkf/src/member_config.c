/*
   Copyright (C) 2011  Statoil ASA, Norway.

   The file 'member_config.c' is part of ERT - Ensemble based Reservoir Tool.

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
#include <time.h>

#include <ert/util/util.h>
#include <ert/util/path_fmt.h>
#include <ert/res_util/subst_list.h>

#include <ert/enkf/ecl_config.h>
#include <ert/enkf/member_config.h>
#include <ert/enkf/enkf_types.h>
#include <ert/enkf/ensemble_config.h>


/**
   This struct contains information which is private to this
   member. It is initialized at object boot time, and (typically) not
   changed during the simulation. [In principle it could change during
   the simulation, but the current API does not support that.]
*/


struct member_config_struct {
  char                  * casename;            /* The name of this case - will mostly be NULL. */
  int iens;
};


/*****************************************************************/


/******************************************************************/
/** Implementation of the member_config struct. All of this implementation
    is private - however some of it is exported through the enkf_state object,
    and it should be perfectly safe to export more of it.
*/


int member_config_get_iens( const member_config_type * member_config ) {
  return member_config->iens;
}





void member_config_free(member_config_type * member_config) {
  util_safe_free(member_config->casename );
  free(member_config);
}



const char * member_config_get_casename( const member_config_type * member_config ) {
  return member_config->casename;
}


member_config_type * member_config_alloc(int iens, const char * casename) {
  member_config_type * member_config = util_malloc( sizeof * member_config );
  member_config->casename            = util_alloc_string_copy( casename );
  member_config->iens = iens;
  return member_config;
}
