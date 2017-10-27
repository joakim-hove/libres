/*
   Copyright (C) 2017  Statoil ASA, Norway.

   The file 'callback_arg.c' is part of ERT - Ensemble based Reservoir Tool.

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

#include <ert/enkf/callback_arg.h>

#define CALLBACK_ARG_TYPE_ID 7814509


void callback_arg_free(callback_arg_type * cb_arg) {
  free( cb_arg );
}



callback_arg_type * callback_arg_alloc(run_arg_type * run_arg, enkf_state_type * enkf_state)
{
  callback_arg_type * cb = util_malloc( sizeof * cb );
  cb->run_arg = run_arg;
  cb->enkf_state = enkf_state;
  return cb;
}



UTIL_IS_INSTANCE_FUNCTION( callback_arg, CALLBACK_ARG_TYPE_ID )
UTIL_SAFE_CAST_FUNCTION( callback_arg, CALLBACK_ARG_TYPE_ID )


