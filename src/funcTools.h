// Header file to concatenate all functions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "time.h"
#include "../umfpk.h"
#include "time.h"
#include "prob.h"
#include "../problem.h"

// Matrix-vector and Matrix-vector operations
double norm(int dim, double *v);
void matVec(int n, int *ia, int *ja, double *a, double *v, double **res);
void addVector(int size, double *v1, double *v2, double **res);
void subVector(int size, double *v1, double *v2, double **res);
void computeRes(int n, int *ia, int *ja, double *a, double *u, double *b, double **r);

// Tools, file managemment and display management
void matrixOutpout(int n, int *ia, int *ja, double *a);
void denseDisplay(int n, int *ia, int *ja, double *a);
void residualDisplay(int amount_iterations, double *resNorm);
void convergenceTimeAnalysis(int size, double *timeVector);
void printVector(int size, double *d);
void displayVector(int m, int size, double *vect);

// Plotting
void showVectorHomogeneious(int m, int size, double *vect, int typeDisplay);   // DESUET
void plot(int m, int size, double *vect, int typeDisplay);
