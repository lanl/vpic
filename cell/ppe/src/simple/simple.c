/* --------------------------------------------------------------- */
/* (C)Copyright 2006,2007,                                         */
/* International Business Machines Corporation,                    */
/* All Rights Reserved.                                            */
/*                                                                 */
/* This program is made available under the terms of the           */
/* Common Public License v1.0 which accompanies this distribution. */
/* --------------------------------------------------------------- */
/* PROLOG END TAG zYx                                              */
#include <stdlib.h>
#include <stdio.h>
#include <libspe2.h>
#include <sys/wait.h>
#include <string.h>

extern spe_program_handle_t simple;

int main() 
{
  spe_context_ptr_t speid;
  int rc;
  spe_stop_info_t stopinfo;
  unsigned int entry = SPE_DEFAULT_ENTRY;

  /* Create context */
  if ((speid = spe_context_create (0, NULL)) == NULL)
    {
      fprintf (stderr, "Failed spe_context_create(errno=%d strerror=%s)\n", errno, strerror(errno));
      exit (1);
    }
  /* Load program */
  if ((rc = spe_program_load (speid, &simple)) != 0)
    {
      fprintf (stderr, "Failed spe_program_load(errno=%d strerror=%s)\n", errno, strerror(errno));
      exit (1);
    }
  /* Run context */
  if ((rc = spe_context_run(speid, &entry, 0, NULL, NULL, &stopinfo)) != 0)
    {
      fprintf (stderr, "Failed spe_context_run(errno=%d strerror=%s)\n", errno, strerror(errno));
      exit (1);
    }
  /* Destroy context */
  if ((rc = spe_context_destroy (speid)) != 0)
    {
      fprintf (stderr, "Failed spe_context_destroy(rc=%d, errno=%d, strerror=%s)\n", rc, errno, strerror(errno));
      exit (1);
    }

  if (stopinfo.stop_reason == SPE_EXIT) {
        return (stopinfo.result.spe_exit_code);
  }
  else {
	fprintf(stderr, "stopinfo.stop_reason=%x, stopinfo.spe_exit_code=%x \n", 
		stopinfo.stop_reason,
		stopinfo.result.spe_exit_code);
        return -1;
  }

}
