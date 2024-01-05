#ifndef DSLA_HPP_
#define DSLA_HPP_

#include <cstdlib>
#include <cstdio>
#include <string>
#include <unistd.h>


#define die(x) dieLoc(x,__FILE__,__LINE__)
void dieLoc(const char* msg, const char* fff, int ll);

#endif
