// Dependencies and variables
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "umfpk.h"
#include "src/time.h"
#include "src/prob.h"
#include "src/funcTools.h"
#include "src/multiGridMethod.h"

#define DL 1.0

// Direct solvers
int primme(int primme_n, int primme_m, int *primme_ia, int *primme_ja, double *primme_a,
           int nev, double *evals, double *evecs, int typePrecond);
