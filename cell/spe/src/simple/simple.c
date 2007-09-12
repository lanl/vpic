/* --------------------------------------------------------------- */
/* (C)Copyright 2006                                               */
/* International Business Machines Corporation,                    */
/* All Rights Reserved.                                            */
/*                                                                 */
/* This program is made available under the terms of the           */
/* Common Public License v1.0 which accompanies this distribution. */
/* --------------------------------------------------------------- */
/* PROLOG END TAG zYx                                              */
#include <stdlib.h>
#include <stdio.h>

/* overlay functions */
extern int o1_test1(int);
extern int o1_test2(int);
extern int o2_test1(int);
extern int o2_test2(int);

#define TEST(x, y)							\
  if ((x) != (y)) {							\
    printf("olay func call failed! expecting %d, got %d\n", x, y);	\
    return (1);								\
  }

int
main (unsigned long long spuid __attribute__ ((__unused__)),
      unsigned long long argp  __attribute__ ((__unused__)),
      unsigned long long envp  __attribute__ ((__unused__)))
{
  int rc;

  /* bring in overlay 1 */
  rc = o1_test1(101);
  TEST(1, rc);
		
  rc = o1_test2(102);
  TEST(2, rc);

  /* bring in overlay 2 */
  rc = o2_test1(201);
  TEST(11, rc);

  rc = o2_test2(202);
  TEST(12, rc);

  /* back to overlay 1 */
  rc = o1_test1(301);
  TEST(1, rc);
		
  rc = o1_test2(302);
  TEST(2, rc);

  return 0;
}
