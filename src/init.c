#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mcmcloop_leroux_gauss(void *, void *,void *,void *,void *,void *,void *,void *,void *,
                                  void *, void *,void *,void *,void *,void *,void *,void *,void *,
                                  void *, void *,void *,void *,void *,void *);

extern void mcmcloop_leroux_GLM(void *,void *,void *,void *,void *,void *,void *,void *,void *,
                                void *,void *,void *,void *,void *,void *,void *,void *,void *,
                                void *,void *,void *,void *,void *,void *,void *,void *,void *,
                                void *);


static const R_CMethodDef CEntries[] = {
    {"mcmcloop_leroux_gauss", (DL_FUNC) &mcmcloop_leroux_gauss, 24},
    {"mcmcloop_leroux_GLM", (DL_FUNC) &mcmcloop_leroux_GLM, 28},
    {NULL, NULL, 0}
};


void R_init_eCAR(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
