#include "multiGridMethod.h"

int multiGridMethod(int m, int n, int *ia, int *ja, double *a, double *b, double *u, int max_it_smoothing, double *resNorm, int ptype, int m_depth, int m_first){
  /*
    INPUT : m : finest discretization (int)
            n : dimension of the problem (int)
            A : discretization Laplacian (ia, ja, a) --> CSR format
            b : Independant boundary conditions (double)*n
            u : Pointer to the solution (double) *n

            max_it_smoothing : amount of Gauss-Seidel/Jacobi smoothing it (int)
            resNorm : pointer to the residual norm (double)
            ptype : smoothing method -> 1 GS / 2 Jacobi
            m_depth : maximum recursion depth (int)
            m_first : Saves the first m (thus the relative position in the recursion)

    OUTPUT : resNorm :  per reference (needs to be saved only once)
             u : Approached solution to the problem (double)*n

    Desc : This is the implementation of the multi-grid method solver for linear
           systems. The methods works as follow :

    multiGridMethod  --> | smoothing u0               | reduction matrix
                         | coarse_grid correction --> | solving coarse problem -->   |multiGridMethod --> UMFPACK
                         | smoothing                  | prolongating coarse to fine

  */

  // Coarce matrix  m_c --------------------------------------------------------
  int m_c = ((m - 1) / 2) + 1;
  int n_c, *ia_c, *ja_c;
  double *a_c, *b_c;
  prob(m_c, &n_c, &ia_c, &ja_c, &a_c, &b_c);

  // Moving vectors and residuals (coarse and fine) ---------------------------------------------
  double *x, *x_c, *r, *r_c;

  // Mem allocation
    // Unknowns
  x = calloc(n, sizeof(double));
  x_c = calloc(n_c, sizeof(double));   // calloc

    //Residuals
  r = malloc(n*sizeof(double));
  r_c = malloc(n_c*sizeof(double));

    // Geometric parameters, storing the coordinates of the bounday points
  int *bounds, *bounds_c;
  bounds = malloc(4*sizeof(int));
  bounds_c = malloc(4*sizeof(int));

  boundaries_matching(n, m, bounds);
  boundaries_matching(n_c, m_c, bounds_c);


  if(x == NULL || x_c == NULL || r == NULL || r_c == NULL){
    printf("Error : impossible to allocate memory multiGridMethod.c -- mem allocation d/d_c/r/r_c\n");
    return 1;
  }

  // =============================================================================
  //                   ## Using the multi-grid method ##
  // =============================================================================

    /*
      -- Pre-smoothing, iterative stationnary method
        ptype == 1 -> Gauss-seidel
        ptype == 2 -> Jacobi
    */


  // ---------------------------------------------------------------------------
  // Smoothing (pre)
  // ---------------------------------------------------------------------------
    stationnaryItMethod(n, ia, ja, a, b, &u, ptype, max_it_smoothing);

    // Computing the initial residual
    computeRes(n, ia, ja, a, u, b, &r);

    // Saving the residual
    if(m == m_first)
      *resNorm = norm(n, r);

    //printf("Multigrille : m = %d | m_c = %d --  résidu = %1e\n",m, m_c, norm(n, r));    //DEBUG


  // ---------------------------------------------------------------------------
  // Reduction matrix
  // ---------------------------------------------------------------------------

    residualTransfert(n, n_c, m, m_c, bounds, bounds_c, r, &r_c);

    // Debug visualization
    //plot(m, n, r, 2);
    //plot(m_c, n_c, r_c, 2);


  // ---------------------------------------------------------------------------
  // Solving the coarse Problem
  // ---------------------------------------------------------------------------

        /*We need to solve the coarced problem associated with

        ---- A_c * x_c = r_c

        As it is a multigrid solver, we either increase the depth or
        we are at the bottom of the recursion.
        */

      // Resetting the x_c vector (debug)
      if(m_c > m_depth){
        if(multiGridMethod(m_c, n_c, ia_c, ja_c, a_c, r_c, x_c, max_it_smoothing, &(*resNorm), ptype, m_depth, m_first))
          return 1;
      }else{
        if(solve_umfpack(n_c, ia_c, ja_c, a_c, r_c, x_c))
          return 1;
      }

  // ---------------------------------------------------------------------------
  // Prolongation matrix
  // ---------------------------------------------------------------------------

  // Setting vector to zero
   prolongationMatrix(n, n_c, m, m_c, bounds, bounds_c, x_c, &x);

   // Debug visualization
   //plot(m, n, x, 2);
   //plot(m_c, n_c, x_c, 2);


  // ---------------------------------------------------------------------------
  // Correcting u
  // ---------------------------------------------------------------------------

   addVector(n, u, x, &u);

/*
   showVectorHomogeneious(m_c, n_c, x_c);
   showVectorHomogeneious(m, n, u);
*/

  // ---------------------------------------------------------------------------
  // Smoothing (post)
  // ---------------------------------------------------------------------------

   stationnaryItMethod(n, ia, ja, a, b, &u, ptype, max_it_smoothing);

    //  coarse problem
    free(ia_c); free(a_c); free(ja_c); free(b_c);
    // unknown
    free(x); free(x_c); free(r); free(r_c);
    // others
    free(bounds); free(bounds_c);
    return 0;
}
