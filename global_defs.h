#ifndef _GLOBAL_DEFS_H_
#define _GLOBAL_DEFS_H_

#include <vector>
#include <cmath>

#ifdef DOUBLE_PRECISION
typedef double real_type;
#else
typedef float real_type;
#endif

typedef std::vector<real_type> real_vector;

#define PI 3.14159265358979
#define degrad 0.01745329251994330   // Convert degrees to radians; e.g. 180 * degrad = 3.14159..
const real_type sqrt2 = sqrt(2.0);
const real_type oneoversqrt2 = 1.0 / sqrt2;

#endif
