#ifndef DSLA_HPP_
#define DSLA_HPP_

#include "dsla.hpp"

#include <cstdlib>
#include <cstdio>
#include <string>
#include <unistd.h>


#define die(x) dieLoc(x,__FILE__,__LINE__)
void dieLoc(const char* msg, const char* fff, int ll){
  printf("\nFatal error:\n\n");
  printf("  %s\n\n",msg);
  printf("  in file '%s', line %d.\n\n",fff,ll);
  #ifdef QCPARALLEL
  printf("  [rank: %i]\n\n",parControl->get_rank());
  #endif
  fflush(stdout);
  char cdebug[1000];
  sprintf(cdebug, "gdb -batch -ex 'bt' -p %d", (int) getpid());
  int sys = system(cdebug);
  printf("sys %i \n", sys);
  exit(1);
}

#endif
