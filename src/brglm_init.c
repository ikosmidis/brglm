#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void hatsc(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"hatsc", (DL_FUNC) &hatsc, 6},
    {NULL, NULL, 0}
};

void R_init_brglm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
