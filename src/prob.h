#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define LS 1.0  // Length of the initial square

int prob(int m, int *n, int **ia, int **ja, double **a, double **b);
int probStat(int m, int *n, int **ia, int **ja, double **a, double *dom);
void defineLine(double x1,double x2,double y1,double y2, int x0, int y0, double *nx0, double *ny0);
