#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ggm_rjmcmc_map(
		void *, void *, void *, void *, void *, void *, void *, void *, 
		void *, void *, void *, void *, void *, void *, void *, void *, 
		void *, void *);
extern void gcgm_rjmcmc_map( 
		void *, void *, void *, void *, void *, void *, void *, void *, 
		void *, void *, void *, void *, void *, void *, void *, void *, 
		void *, void *, void *, void *, void *, void *, void *);
extern void copula( void *, void *, void *, void *, void *);
extern void check_nthread( void *);
extern void check_os( void * );
extern void omp_set_num_cores( void *);
extern void rgwish_c( void *, void *, void *, void *, void *);
extern void rwish_c( void *, void *, void *, void *);


static const R_CMethodDef CEntries[] = {
    {"ggm_rjmcmc_map", (DL_FUNC) &ggm_rjmcmc_map, 18},
	{"gcgm_rjmcmc_map", (DL_FUNC) &gcgm_rjmcmc_map,23},
	{"copula", (DL_FUNC) &copula,5},
	{"check_nthread", (DL_FUNC) &check_nthread, 1},
	{"check_os", (DL_FUNC) &check_os,1},
	{"omp_set_num_cores", (DL_FUNC) &omp_set_num_cores, 1},
	{"rgwish_c", (DL_FUNC) &rgwish_c, 5},
	{"rwish_c", (DL_FUNC) &rwish_c, 4},
    {NULL, NULL, 0}
};

void R_init_bdgraphStudy(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
