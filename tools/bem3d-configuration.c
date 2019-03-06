/* bem3d-configuration.c
 * 
 * Copyright (C) 2017, 2018 Michael Carley
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

gint main(gint argc, gchar **argv)

{
  BEM3DConfiguration *config ;
  gchar ch, *progname ;
  GLogLevelFlags loglevel ;
  gboolean write_config, write_description ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  bem3d_configuration_init() ;
  config = bem3d_configuration_new() ;
  write_config = FALSE ;
  write_description = FALSE ;
  
  loglevel = G_LOG_LEVEL_MESSAGE ;

  while ( (ch = getopt(argc, argv, "hC:dg")) != EOF ) {
    switch (ch) {
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: generate and manipulate BEM3D configuration files\n\n",
	      progname) ;
      fprintf(stderr, "Usage: %s [options] > file.cfg\n", 
	      progname) ;
      fprintf(stderr, 
	      "Options:\n\n"
	      "  -h print this message and exit\n"
	      "  -C <name of configuration file to read>\n"
	      "  -d if writing configuration, add descriptions of settings "
	      "as comments\n"
	      "  -g write generic configuration file to stdout\n") ;
      return 0 ;
      break ;
    case 'C': g_assert_not_reached() ;
      bem3d_configuration_read(config, optarg) ;
      break ;
    case 'd': write_description = TRUE ; break ;
    case 'g': write_config = TRUE ; break ;
    }
  }

  bem3d_logging_init(stderr, "", loglevel, NULL) ;
  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  if ( write_config )
    bem3d_configuration_write_generic(write_description, 70, "#", "## ",
				      stdout) ;

  return 0 ;
}
