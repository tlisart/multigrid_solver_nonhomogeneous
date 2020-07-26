#include <stdlib.h>
#include "funcTools.h"
void smoothingJacobi(int n, double *r, double *d, double **buff);
void smoothingGSG(int n, int *ia, int *ja, double *a, double *r, double *d, double **buff);
void solveTrigSup(int n, int *ia, int *ja, double *a, double *r, double **buff);
void solveTrigInf(int n, int *ia, int *ja, double *a, double *r, double **buff);
int stationnaryItMethod(int n, int *ia, int *ja, double *a, double *b, double **u, int ptype, int max_it_smoothing);
