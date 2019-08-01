#ifndef _GLOBAL_DEFS_H_
#define _GLOBAL_DEFS_H_

#include <cmath>

#ifdef DOUBLE_PRECISION
typedef double real_type;
#else
typedef float real_type;
#endif

const real_type sqrt2 = sqrt(2.0);
const real_type oneoversqrt2 = 1.0 / sqrt2;

#endif
