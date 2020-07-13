#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <float.h>

#include "sexp_macros.h"

#define ALLOC_REAL_VECTOR(S, D, N)                                             \
SEXP S;                                                                        \
PROTECT(S = allocVector(REALSXP, N));                                          \
double *D = REAL(S);

#ifndef MAX
#define MAX(A, B) ((A > B) ? (A) : (B))
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>

#ifdef DEBUG
#define SAMPLER_DEBUG(N, A, B) Rprintf("%8s(%f, %f)\n", N, A, B)
#else
#define SAMPLER_DEBUG(N, A, B)
#endif

extern "C"{
  double do_rtruncnorm_study(double a, double b, double mean, double sd);
}
