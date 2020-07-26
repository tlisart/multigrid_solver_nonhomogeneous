#include "prob.h"
#include "smoothing.h"
#include "geometry_vec.h"
#include "reduction_mat.h"
#include "prolongation_mat.h"

int multiGridMethod(int m, int n, int *ia, int *ja, double *a, double *b, double *u, int max_it_smoothing, double *resNorm, int ptype, int m_depth, int m_first);
