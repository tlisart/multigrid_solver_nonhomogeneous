// Dependencies and variables
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void residualTransfert(int n, int n_c, int m, int m_c, int *bounds, int *bounds_c, double *r, double **r_c);
int matching_fine(int m, int n, int *bounds, int x, int y, int version);
