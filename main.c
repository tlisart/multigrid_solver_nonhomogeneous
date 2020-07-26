#include "main.h"
#include "problem.h"

int main(int argc, char *argv[])   {

  /* Declaring problem's variables*/
  int m = 513, nev = 1;               /*m the amount of discretized points in one direction (fine)*/
  int m_depth = ((m - 1) / 2) + 1;    // Two grid
  //m_depth = 9;                      // Multi-grid

  double *evals, *evecs;               /*eigenvalues and vectors*/

  double t1, t2;                       /*Time values */
  double *u, *b;
  int n;

// =============================================================================
//                     ## Problem's parameters ##
// =============================================================================

  double tau = 1.232827;      // When stable, show optimal relaxation has to be found exerimentaly, \tau = 2.0 / (2.0 - rho)

  int ptype = 1;              // 1 --> Gauss-Seidel, 2 --> Jacobi method
  int relaxationOn = 1;       // 1 --> Searching optimal \tau for the problem
  int stabCheck = 1;          // 0 -> Stability check off || 1 --> Stability check on
  int primmePrecond = 1;      // 0 -> Skips primme
  int primmePrecondType = 1;  // 1 --> multigrid; 2--> twogrid

  // Graphical solutions
  int showDiffusionMap = 1;    // Displays the results (non homogeneous BC)
  int showDiffusionSurface = 1;    // Displays the results (non homogeneous BC)
  int showEigenvalues = 1;  // Dispays PRIMME results eigenvalues
  int showResGraph = 1;


  //rates and relaxation parameters
  double resNorm = 10.0;
  double rho;            // Final/current convergence rate
  double rho_as = 0.0;   // Asymptotical rate (usual use of rho)
  double rho_th = 0.0;   // Theoritical rho if new \tau chosen
  double *rate;

  int stagn = 0;         // 0 --> rho doesn't sagnates || 1 --> rho stagnates
  double *res;


// =============================================================================
//                   ## Q1 - Non homogeneious Dirichlet BC   ##
// =============================================================================

  int *ia, *ja;
  double *a;

  if (prob(m, &n, &ia, &ja, &a, &b))
     return 1;

  // Initializing solution vector  u0 = vect(0)
  u = calloc(n, sizeof(double));
  if(u == NULL)
    printf("Error : impossible to allocate memory u in main.c\n");

/*
// For comparing to the solution --  Direct solver
    if(solve_umfpack(n, ia, ja, a, b, u))
      return 1;

    for(int i = 0; i < n; i++){u[i] = 0.0;}
*/

// =============================================================================
//                     ## Initializing ##
// =============================================================================

      printf("\n// ================================================================================================== \n");
      printf("\nGenerating multi-grid problem : maximum coarse depth m_c = %d\n", m_depth);
      printf("fine matrix m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n]);
      printf("\n");

      if(ptype == 1){printf("Smoothing method : Gauss-Seidel || ");}
      if(ptype == 2){printf("Smoothing method : Jacobi || ");}

      if(stabCheck == 0){
        if(relaxationOn == 0){
          printf("Iteration / smoothing : %d || Direct stability check OFF || Relaxation OFF || \n", MAXIT_PRECOND);
        }else if(relaxationOn == 1){
          printf("Iteration / smoothing : %d || Direct stability check OFF || Relaxation ON - tau = %f || \n", MAXIT_PRECOND, tau);
        }
      }
      if(stabCheck == 1){
        if(relaxationOn == 0){
          printf("Iteration / smoothing : %d || Direct stability check ON || Relaxation OFF ||\n", MAXIT_PRECOND);
        }else if(relaxationOn == 1){
          printf("Iteration / smoothing : %d || Direct stability check ON || Relaxation  ON - tau = %f ||\n", MAXIT_PRECOND, tau);
        }
      }

      printf("\n// ================================================================================================== \n");



// =============================================================================
//                       ## Q2 - Two-grid method   ##
// =============================================================================

  /* The two-grid method solves the linar system b =/= 0 in two steps :
      1 -- Stationnary iteration resolution method with GS preconditionning
      2 -- Coarsed grid correction
          -- Transfering the residual to coarse dimension
          -- Solving the coarse problem
          -- Prolongating the coarse to fine  x_c --> x
      */

// =============================================================================
//                       ## Q3 - Multi-grid method   ##
// =============================================================================


int count;
int foundRho = 0;

res = malloc(MAXIT_MET*sizeof(double));
rate = malloc(MAXIT_MET*sizeof(double));

// Solution buffer
double *u_buff;
u_buff = malloc(n*sizeof(double));

// Clock
t1 = mytimer();

while (count < MAXIT_MET && resNorm > TOL && stagn == 0){

  // Creating buffer to save the previous solution (relaxation)
  for(int i = 0; i < n; i++){
    u_buff[i] = u[i];
  }

  multiGridMethod(m, n, ia, ja, a, b, u, MAXIT_PRECOND, &resNorm, ptype, m_depth, m);

  // Recording the residual
  res[count] = resNorm;

  if(count > 0){
    rate[count] = res[count]/res[count - 1];
    rho = res[count]/res[count - 1];
    if(count > 1){
      if((rate[count] < 1.05 * rate[count -1]) && foundRho == 0){
        rho_as = rate[count];
        foundRho = 1;
      }
    }
  }


  // Relaxation (could be implemented in the stationnary function)
  for(int i = 0; i < n; i++){
    u[i] = u_buff[i] + tau*(u[i] - u_buff[i]);
  }
  // Checking convergenc rate
  if(1.2*res[count] > res[count - 1] && count > 5){stagn = 1;}

  count++;
  printf("|||  Residual norm : %1e  ||   Convergence rate : %f  ||   at iteration : %d ||\n",resNorm, rho, count);
  fflush(stdout);
}

t2 = mytimer();

// Graphical analysis
if(showDiffusionMap == 1)
  plot(m, n, u, 1);

if(showDiffusionSurface == 1)
  showVectorHomogeneious(m, n, u, 2);
if(showResGraph == 1){
  if(showEigenvalues == 1){
    printf("\n// ================================================================================================== \n");
    printf("Please close the residual graph to carry on with the eigenvalue problem\n");
  }
  residualDisplay(count, res);
}

// Relaxation acceleration
if(relaxationOn == 1){
  tau = 2.0/(2.0 - rho_as);
}

// =============================================================================
//                       ## Q4 - Results analysis ##
// =============================================================================

  // Direct Stability --> setting the directStabilityCheck (Q4.1)

    printf("\n// ================================================================================================== \n");
    if(resNorm <= TOL){printf("\nMinimum tol reached a it : %d\n", count);}
    if(count == MAXIT_MET){printf("Max iteration reached : %d\n", count);}
    if(stagn == 1){printf("The residual is stagnating at iteration : %d\n", count);}
    printf("Resolution Time multi grid method (CPU): %5.1f sec\n",t2-t1);
    printf("Final residual norm : %1e\n",resNorm);

    if(relaxationOn == 0){printf("Final convergence rate at stop : %f\n",rate[count - 1]);}
    if(relaxationOn == 1){

      printf("Final convergence rate at stop : %f\n",rate[count -1]);

      // Asymptotic convergence rate (Q4.2)
      if(rho_as == 0.0){
        printf("Asymptotic rate not computed (not enough iterations or instability)\n");
      }else{
        printf("Asymptotic rate rho = %1f , optimal results for tau = 2/(lambda_max - lambda_min) ~= (2/(2 - rho)) = %f\n", rho_as, tau);
        rho_th = rho_as/(2 - rho_as);
        printf("If chosen, the new relaxation coefficient should give a theoritical rho = %f\n", rho_th);
      }
    }

    if(stabCheck == 1){
      double eps_rel = res[count - 1]/norm(n, b);        // Relative error ||r||/||b||
      double h2 = pow(LS, 2)/((m-1)*(m-1));          // Discretization step

      // From the definitions ||b|| simplifies -->
      double eps = h2*res[count - 1]/(LS*norm(n, u)*8);  // Convergence number


      printf("\nConvergence quality analysis :\n");
      printf("     --->   Direct stability criterion : eps = ||r||/||b|| < O(u)*||A||_2*||u||/||b|| ?\n");
      printf("     --->   Relative error : eps_err = %1e || convergence number : eps = %1e\n", eps_rel, eps);
      printf("     --->   eps = %1e < u = 1.1e-6 ? \n",eps);

      /*The relative error is given by ||b||/||u||*/
      if(eps < 1.1e-6){
        printf("\n ----- The stability criterion is FULFILLED -----\n");
      }else{printf("The stability criterion is NOT fulfilled\n");}
    }


// =============================================================================
//                           ## PRIMME assessment ##
// =============================================================================

// Homogeneous problem with acceleration
if(primmePrecond == 1){
  printf("\n// ================================================================================================== \n");
  printf("Eigenvalue problem with primme \n");

  if(primmePrecondType == 1)
    printf("Preconditionner to PRIMME --> MULTIGRID \n");

  if(primmePrecondType == 2)
    printf("Preconditionner to PRIMME --> TWO-GRID\n");

  printf("\n// ================================================================================================== \n");

  evals = malloc(nev * sizeof(double));
  evecs = malloc(nev * n * sizeof(double));
  if (evals == NULL || evecs == NULL) {
    printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
    return 1;
  }

  t1 = mytimer();
  if(primme(n, m, ia, ja, a, nev, evals, evecs, primmePrecondType))
     return 1;
  t2 = mytimer();

  // temps de solution
  printf("\n  Temps de solution (CPU) (PRIMME): %5.1f sec\n",t2-t1);
  printf("  Valeur propre minimale calculée (PRIMME): %5.2f \n",evals[0]);

  if(showEigenvalues == 1){
    plot(m, n, evecs, 2);
  }
}

printf("\n// ================================================================================================== \n");

free(ia); free(ja); free(a); free(b); free(u); free(u_buff);
free(evals); free(evecs); free(res); free(rate);
}
