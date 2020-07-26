#include "smoothing.h"

// =============================================================================
//       ## Implementation of a stationnary iterative linear solver ##
// =============================================================================


int stationnaryItMethod(int n, int *ia, int *ja, double *a, double *b, double **u, int ptype, int max_it_smoothing){
  /*
    INPUT : see multiGridMethod.cs
    OUTPUT : approached u per reference

    Desc :  Implementing a stationnary iterative resolution with Jacobi or Sym Gauss-Seidel (SGS)
            Preconditionning of B. Solving iteratively the problem Au = b.

            Algorithm :
                A = La + Ua - Da

                _init_ :
                  u0 = arbitrary (zeroes)
                  r0 = b - A*u0
                loop :
                  d_m = B^(-1)*r_m

                  --> Solving for B*d_m = r_m where d_m unknown (B = B1/B2)

                  u_(m+1) = u_m + \tau*d_m
                  r_(m+1) = b - A*u_(m+1)

                  --> Acceleration (\tau corr) implemented in main.c
  */

  double *d, *r, *buff;

  r = malloc(n*sizeof(double));    // Residual
  d = malloc(n*sizeof(double));   // Diagonal

  // Will be used as a moving vector
  buff = malloc(n*sizeof(double));

  if(r == NULL || d == NULL || buff == NULL){
    printf("Error : impossible to allocate memory smoothing.c -- mem allocation r/d/buff\n");
    return 1;
  }

  // Saving the diagonal for later use
  for(int i = 0; i < n; i++){
    for(int k =  ia[i]; k < ia[i + 1]; k++){   //  Runs through the whole CSR mat
      if(ja[k] == i){  // Selecting only diagonal elements
        d[i] = a[k];
      }
    }
  }

  // Smoothing -- Stationnary iterative method (Jacoby and SGS precond)
  int count = 0;
  while (count < max_it_smoothing){

      // __init__
      matVec(n, ia, ja, a, *u, &r);      // A*u0
      subVector(n, b, r, &r);            // b - A*u0, as u0 = 0 ==> b

      // Smoothing method -- Choice of B matrix
      if(ptype == 1){
        smoothingGSG(n, ia, ja, a, r, d, &buff);
        // Stationnary method
        addVector(n, *u, r, &(*u));         // u_(m+1) = u_m + \tau*d_m
      }

      if(ptype == 2){
        smoothingJacobi(n, r, d, &buff);    // Approaching B by its diagonal
        // Stationnary method
        addVector(n, *u, buff, &(*u));     // u_(m+1) = u_m + \tau*d_m
      }
      count ++;  // Counting iterations
  }
  free(r); free(d); free(buff);
}


void smoothingGSG(int n, int *ia, int *ja, double *a, double *r, double *d, double **buff){
  /*
    INPUT : n, A, diagonnal, buff pointer
    OUTPUT : buff (approached solution after one GS swipe)

    Desc : triangular matrix swiping, summetric Gauss-Seidel expect one solve
           triangular superior one solve triangular inferior after the other to
           preserve the symmetric properties of the matrix.

           Solves the triangular inf, then triangular sup.
  */

    solveTrigInf(n, ia, ja, a, r, &(*buff));
    for(int i = 0; i < n; i++){(*buff)[i] *= d[i];}

    solveTrigSup(n, ia, ja, a, *buff, &r);
  }


void solveTrigInf(int n, int *ia, int *ja, double *a, double *r, double **buff){
  /*
    INPUT : n, A, diagonnal, buff pointer
    OUTPUT : buff (approached solution after one GS swipe)

    Desc : solve triangular inferior problem after the other to preserve the
           symmetric properties of the matrix.
  */

  for(int i = 0; i <n; i++){
    (*buff)[i] = r[i];
    for(int k = ia[i]; k < ia[i + 1]; k++){   //  Runs through the whole CSR mat

      if(ja[k] < i){  // Selecting only trig inf elements
        (*buff)[i] -= (*buff)[ja[k]]*a[k];
      }
      else if(ja[k] == i){
        (*buff)[i] = (*buff)[i]/a[k];
      }
    }
  }
}

void solveTrigSup(int n, int *ia, int *ja, double *a, double *r, double **buff){
  /*
    INPUT : n, A, diagonnal, buff pointer
    OUTPUT : buff (approached solution after one GS swipe)

    Desc : solve triangular superior problem after the other to preserve the
           symmetric properties of the matrix.
  */

  for(int i = n - 1; i >= 0; i--){
    (*buff)[i] = r[i];
    for(int k = ia[i + 1] - 1; k >= ia[i]; k--){   //  Runs through the whole CSR mat

      if(ja[k] > i){  // Selecting only trig inf elements
        (*buff)[i] -= (*buff)[ja[k]]*a[k];
      }
      else if(ja[k] == i){   // Diagonal
        (*buff)[i] = (*buff)[i]/a[k];
      }
    }
  }
}


void smoothingJacobi(int n, double *r, double *d, double **buff){
  /*
    INPUT : d : diagonal vector
            r : initial residual
    OUTPUT : buff : solution vector approached

    Desc : The simplest preconditionner, B is approached by the diagonnal of A
   */

  for(int i = 0; i < n; i++){
    (*buff)[i] = (1/d[i]) * r[i];
  }
}
