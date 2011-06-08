#ifndef ylm_h
#define ylm_h ylm_h
#include <stdio.h>
#include <math.h>
#include <complex>
#include <factorial.h>

#define complex std::complex<double>

complex ylm(long l, long m, double theta, double phi);

#endif // ylm_h





