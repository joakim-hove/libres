/*
   Copyright (C) 2013  Statoil ASA, Norway.

   The file 'enkf_site_config.c' is part of ERT - Ensemble based Reservoir Tool.

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
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <ert/util/test_util.h>
#include <ert/util/util.h>
#include <ert/util/arg_pack.h>
#include <ert/util/test_work_area.h>

#include <ert/config/config_parser.hpp>
#include <ert/config/config_content.hpp>

#include <ert/ecl/ecl_sum.h>

#include <ert/enkf/site_config.hpp>


#include <ert/enkf/site_config.hpp>




#define INCLUDE_KEY "INCLUDE"
#define DEFINE_KEY  "DEFINE"


void test_init(const char * config_file) {
  site_config_type * site_config = site_config_alloc_load_user_config(NULL);
  site_config_free( site_config );
}


void test_job_script() {
  test_work_area_type * test_area = test_work_area_alloc("site-config");
  {
    site_config_type * site_config = site_config_alloc_load_user_config(NULL);


    test_assert_false( site_config_set_job_script( site_config , "/does/not/exist" ));


    {
      FILE * job_script = util_fopen("Script.sh" , "w");
      fclose( job_script );
    }
    test_assert_false( site_config_set_job_script( site_config , "Script.sh" ));

    chmod("Script.sh" , S_IRWXU );
    test_assert_true( site_config_set_job_script( site_config , "Script.sh" ));


    test_assert_false( site_config_set_job_script( site_config , "DoesNotExits"));

    {
      char * full_path = util_alloc_realpath( "Script.sh" );
      queue_config_type * queue_config = site_config_get_queue_config( site_config );
      test_assert_string_equal( full_path , queue_config_get_job_script( queue_config));
      free( full_path );
    }
    site_config_free( site_config );
  }
  test_work_area_free( test_area );
}



int main(int argc , char ** argv) {
  const char * site_config_file = argv[1];

  util_install_signals();

  test_init( site_config_file );
  test_job_script();

  exit(0);
}

