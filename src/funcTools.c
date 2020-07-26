#include "funcTools.h"

// =============================================================================
//                           # Residual analysis #
// =============================================================================

void computeRes(int n, int *ia, int *ja, double *a, double *u, double *b, double **r){
  /*
    INPUT : A matrix (ia, ja, a), b independant member, u approched solution, size n
    OUTPUT : residual r target

    Desc : Computes the residual vector b - A*u by reference r.
           r must be allocated outside of the function (size n)
  */

  double *vec;
  vec = malloc(n*sizeof(double));
  if(vec == NULL)
    printf("Error : impossible to allocate memory vec in computeRes() --> functools.c\n");

  matVec(n, ia, ja, a, u, &vec);   // A*u

  for(int i = 0; i < n; i++){
    (*r)[i] = b[i] - vec[i];   // b - Au
  }
  free(vec);
}

// =============================================================================
//                     # Matrices and vector operations #
// =============================================================================


void matVec(int n, int *ia, int *ja, double *a, double *v, double **res){
  /*
    INPUT : A : matrix (ia, ja, a) CSR format (n x n dimension)
            v : vector of size n
    OUTPUT : res per referencing

    Desc : Computes the matrix / vector multiplication (matrix linear transform)
           to a vector v of size n. memory allocation outside of the function
  */

  for(int i = 0; i < n; i++){
    (*res)[i] = 0.0;
    for(int k = ia[i]; k < ia[i + 1]; k++){
      (*res)[i] = (*res)[i] + a[k]* v[ja[k]];
    }
  }
}

double norm(int n, double *v){
  /*
    INPUT : n : dimension of the vector
            v : pointer to the array
    OUTPUT : (double) prod : cartesian norm of v

    Desc : Self explanatory
  */

  double prod = 0.0;
  for(int i = 0; i < n; i++){
    prod += v[i]*v[i];     // Multiplication euclidienne
  }
  return sqrt(prod);
}

void addVector(int n, double *v1, double *v2, double **res){
  /*
    INPUT : v1, v2 vector of size n
    OUTPUT : res solution target

    Desc : Adds two vectors of same size together
           res memory must be allocated outside of the function
  */

  for(int i = 0; i < n; i++){
    (*res)[i] = v1[i] + v2[i];
  }
}

void subVector(int size, double *v1, double *v2, double **res){
  /*
    INPUT : v1, v2 vector of size n
    OUTPUT : res solution target

    Desc : Substract two vectors of same size together
           res memory must be allocated outside of the function
  */

  for(int i = 0; i < size; i++){
    (*res)[i] = v1[i] - v2[i];
  }
}

// =============================================================================
//             # File management, data piping and graphical tools #
// =============================================================================


void denseDisplay(int n, int *ia, int *ja, double *a){
  /*
    INPUT : n : array size
            A : problem matrix (ia, ja, a)
    OUTPUT : pipes data to python

    Desc : Calls a python routine to display the dense matrix associated
           to the problem
  */

  matrixOutpout(n, ia, ja, a);

  // Piping into Python
  char cmd[MAX_CMD_LINE] = "", **p;
  strcat(cmd, "python2.7 data/matrices/display.py");
  if(system(cmd) == 1){printf("An error has occured");};
}

void matrixOutpout(int n, int *ia, int *ja, double *a){
  /*
    INPUT : n : size of the matrix
            A : problem matrix (ia, ja, a)
    OUTPUT : data text files a_txt, ja_txt, ia_txt

    Desc : Saves to files the matrix informations as plain text
  */

  FILE *a_txt;
  FILE *ja_txt;
  FILE *ia_txt;

  a_txt = fopen("data/matrices/a.txt", "w+");
  ja_txt = fopen("data/matrices/ja.txt", "w+");
  ia_txt = fopen("data/matrices/ia.txt", "w+");

  for (int i = 0; i <= n; i++){fprintf(ia_txt, "%d\n",ia[i]);}
  fclose(ia_txt);

  for (int i = 0; i < ia[n]; i++){
    fprintf(ja_txt, "%d\n", ja[i]);
    fprintf(a_txt, "%f\n",a[i]);
    }

  fclose(a_txt);
  fclose(ja_txt);
}

void residualDisplay(int amount_iterations, double *resNorm){
  /*
    INPUT : amount_iterations : vector size
            resNorm : Array with all residual norms
    OUTPUT : textfile : resData

    Desc : Takes all residual norms per iterations and pipes the results
           into python to be displayed as a graph
  */

  FILE *data;
  data = fopen("data/plotRes.txt", "w+");

  for (int i = 0; i < amount_iterations; i++){fprintf(data, "%1e\n",resNorm[i]);}
  fclose(data);

  // Piping into Python
  char cmd[MAX_CMD_LINE] = "", **p;
  strcat(cmd, "python2.7 data/plotRes.py");
  if(system(cmd) == 1){printf("An error has occured");};
}

void convergenceTimeAnalysis(int amount_iterations, double *timeVector){
  /*
    INPUT : amount_iterations : vector size
            resNorm : Array with all residual norms
    OUTPUT : textfile : timeDate

    Desc : Takes all time spent per iterations and pipes the results
           into python to be displayed as a graph
  */

  FILE *data;
  data = fopen("data/timeData.txt", "w+");

  for (int i = 0; i < amount_iterations; i++){fprintf(data, "%f\n",timeVector[i]);}
  fclose(data);

  // Piping into Python
  char cmd[MAX_CMD_LINE] = "", **p;
  strcat(cmd, "python2.7 data/plotTime.py");
  if(system(cmd) == 1){printf("An error has occured");};
}


void printVector(int n, double *d){
  /*
    INPUT : n : vector dimention
            d : pointer to vector

    Desc : prints a vector to the consol
  */

  for(int i = 0; i < n; i++){
    printf("%1e\n",d[i]);
  }
}

void plot(int m, int size, double *vect, int typeDisplay){
  /*
    INPUT : m : amount of discretization in one direction
            size : dimension of the solution
            vect : pointer to the vector
            typeDisplay : Diffusion heatmap -> 1 / pm3d surface --> 2
    OUTPUT : data.txt and cmd.txt files

    Desc : Pipes the results data to GNUplot and GNUplot configuration
           Formats 1:2:3 (x, y, z) for GNUplot, scaling necessary by h.
  */

  FILE *data;
  FILE *cmd;

  data = fopen("data.txt", "w+");
  cmd = fopen("cmd.txt", "w+");

  if(data == NULL){printf("File could not be created");}

  int  nnz, nx, x0, y0, imax, imin, jmax, jmin;    /*Problem's parameters -- coarse*/

  nx = m - 2; /* noeuds de Dirichlet ne sont pas pris en compte == bords */
  nnz = 0;

  double h = LS/(m-1);

  double L = fabs(X2 - X1);
  double l = fabs(Y2 - Y1);

  x0 = round(L*(m - 1));
  y0 = round(l*(m - 1));

  double nx0, nx0F;
  double ny0, ny0F;

  defineLine(X1 ,X2 ,Y1 ,Y2 ,x0 ,y0 ,&nx0 ,&ny0);

  imin = ny0 - 3;
  imax = 2*ny0 - 2;
  jmin = 0;
  jmax = nx0;

  // Let us create a x - y - z coordinates system to be displayed
  // Value "hack" to obtain the right scale --DEBUG should be corrected

  int ind = 0;
  double value = 0;

  fprintf(data," ");
  for(int iy = 0; iy < nx; iy++){

    for(int ix = 0; ix < nx; ix++){
      if( !((iy > imin && iy < imax) &&  ix < jmax)){  // In the problem
        if(iy <= imin){
          ind = iy*nx + ix;
          value = vect[ind];
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
        else if(iy > imin && iy < imax){
          ind = (iy - imin)*(nx - nx0) + imin*nx + ix;
          value = vect[ind];
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }else if(iy >= imax){
          ind = imin*nx +(imax - imin)*(nx - nx0) + (iy - imax)*nx + jmax + ix;
          value = vect[ind];
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
      }
      else{
        // Superior border
        if((iy == imax) && ((ix < jmax) && (ix > jmin))){
          value = 1.0;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
        // Inferior border
        if((iy == imin) && ((ix < jmax) && (ix > jmin))){
          value = 1.0;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
        // Right hand side border
        if(ix == jmax && ( (iy < imax) && (iy > imin))){
          value = 1.0;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        } // Zeroes
        else{
          value = -INFINITY;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
      }
    }
    fprintf(data," \n ");
    ind++;
  }

  // Piping into GNUplot
  if(typeDisplay == 1){       // Heatmap
    fprintf(cmd, "set view 50,55\n");
    fprintf(cmd, "set ticslevel 0\n");
    fprintf(cmd, "set terminal wxt size 800,600\n");

    fprintf(cmd, "set title 'Graph of the solution for m = %d (diffusion problem)' \n", m);
    fprintf(cmd, "set view map \n");
    fprintf(cmd, "set size ratio -1 \n");
    fprintf(cmd, "splot 'data.txt' using 1:2:3 with pm3d\n");
  }

  else if(typeDisplay == 2){   // Surface
    fprintf(cmd, "set view 50,55\n");
    fprintf(cmd, "set ticslevel 0\n");
    fprintf(cmd, "set terminal wxt size 800,600\n");

    fprintf(cmd, "set title 'Graph of the first eigenvector m = %d' \n", m);
    fprintf(cmd, "set size ratio -1 \n");
    fprintf(cmd, "splot 'data.txt' using 1:2:3 with pm3d\n");
  }

  fclose(cmd);
  fclose(data);


  if(system("gnuplot -persistent 'cmd.txt' "))
    printf("ERREUR GNUplot");
}


void showVectorHomogeneious(int m, int size, double *vect, int typeDisplay){
  /*
    INPUT : m : amount of discretization in one direction
            size : dimension of the solution
            vect : pointer to the vector
            typeDisplay : Diffusion heatmap -> 1 / pm3d surface --> 2
    OUTPUT : data.txt and cmd.txt files

    Desc : Pipes the results data to GNUplot and GNUplot configuration
           Formats 1:2:3 (x, y, z) for GNUplot, scaling necessary by h.
  */

  FILE *data;
  FILE *cmd;

  data = fopen("data.txt", "w+");
  cmd = fopen("cmd.txt", "w+");

  if(data == NULL){printf("File could not be created");}

  int  nnz, nx, x0, y0, imax, imin, jmax, jmin;    /*Problem's parameters -- coarse*/

  nx = m - 2; /* noeuds de Dirichlet ne sont pas pris en compte == bords */
  nnz = 0;

  double h = LS/(m-1);

  double L = fabs(X2 - X1);
  double l = fabs(Y2 - Y1);

  x0 = round(L*(m - 1));
  y0 = round(l*(m - 1));

  double nx0, nx0F;
  double ny0, ny0F;

  defineLine(X1 ,X2 ,Y1 ,Y2 ,x0 ,y0 ,&nx0 ,&ny0);

  imin = ny0 - 3;
  imax = 2*ny0 - 2;
  jmin = 0;
  jmax = nx0;

  // Let us create a x - y - z coordinates system to be displayed
  // Value "hack" to obtain the right scale --DEBUG should be corrected

  int ind = 0;
  double value = 0;

  fprintf(data," ");
  for(int iy = 0; iy < nx; iy++){

    for(int ix = 0; ix < nx; ix++){
      if( !((iy > imin && iy < imax) &&  ix < jmax)){  // In the problem
        if(iy <= imin){
          ind = iy*nx + ix;
          value = vect[ind];
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
        else if(iy > imin && iy < imax){
          ind = (iy - imin)*(nx - nx0) + imin*nx + ix;
          value = vect[ind];
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }else if(iy >= imax){
          ind = imin*nx +(imax - imin)*(nx - nx0) + (iy - imax)*nx + jmax + ix;
          value = vect[ind];
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
      }
      else{
        // Superior border
        if((iy == imax) && ((ix < jmax) && (ix > jmin))){
          value = 1.0;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
        // Inferior border
        if((iy == imin) && ((ix < jmax) && (ix > jmin))){
          value = 1.0;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
        // Right hand side border
        if(ix == jmax && ( (iy < imax) && (iy > imin))){
          value = 1.0;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        } // Zeroes
        else{
          value = -INFINITY;
          fprintf(data," %f %f %f \n",(ix+ 1)*h, (iy + 1)*h, value * 1e6 / 4);
        }
      }
    }
    fprintf(data," \n ");
    ind++;
  }

  // Piping into GNUplot
  if(typeDisplay == 1){       // Heatmap
    fprintf(cmd, "set view 50,55\n");
    fprintf(cmd, "set ticslevel 0\n");
    fprintf(cmd, "set terminal wxt size 800,600\n");

    fprintf(cmd, "set title 'Graph of the solution for m = %d (diffusion problem)' \n", m);
    fprintf(cmd, "set view map \n");
    fprintf(cmd, "set size ratio -1 \n");
    fprintf(cmd, "splot 'data.txt' using 1:2:3 with pm3d\n");
  }

  else if(typeDisplay == 2){   // Surface
    fprintf(cmd, "set view 50,55\n");
    fprintf(cmd, "set ticslevel 0\n");
    fprintf(cmd, "set terminal wxt size 800,600\n");

    fprintf(cmd, "set title 'Graph of the solution for m = %d (diffusion problem)' \n", m);
    fprintf(cmd, "set size ratio -1 \n");
    fprintf(cmd, "splot 'data.txt' using 1:2:3 with pm3d\n");
  }

  fclose(cmd);
  fclose(data);


  if(system("gnuplot -persistent 'cmd.txt' "))
    printf("ERREUR GNUplot");

}
