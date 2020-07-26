#include "geometry_vec.h"

void boundaries_matching(int n, int m, int *boundaries){
  /*
    INPUT : n : dimension
            m : discretisation
            boundaries : pointer to 4 size (int)
    OUTPUT : Fills boundaries vector.

    Desc : saves jmin,jmax, imin, imax to compute the geometry of the problem
           only once.
  */

  int  nx, x0, y0;
  double *dom = malloc(4*sizeof(double));
  if(dom == NULL)
    printf("Error : impossible to allocate memory dom in boundaries_matching() --> geometry_vec.c\n");

  dom[0] = 0.0; dom[1] = 1.0 / 4.0; dom[2] = 5.0 / 8.0; dom[3] = 1.0 / 2.0;

//############################# Problème carré ##############################

  nx = m - 2; /* noeuds de Dirichlet ne sont pas pris en compte == bords */


  double L = fabs(dom[2] - dom[0]);
  double l = fabs(dom[3]- dom[1]);

  x0 = round(L*(m-1));
  y0 = round(l*(m-1));

  // On converti en points discrétisés, avec prise en compte des bords

  double nx0;
  double ny0;

  defineLine(dom[0] ,dom[2] ,dom[1] ,dom[3] ,x0 ,y0 ,&nx0 ,&ny0);

  // ############### Détermination des indices du "trou" ########################
  // On cherche les indices juste " à coté" du trou, pour déterminer les cas limites

  /*
    imin = 0
    imax = 1
    jmin = 2
    jmax = 3
  */

  boundaries[0] = ny0 - 3;
  boundaries[1] = 2*ny0 - 2;
  boundaries[2] = 0;
  boundaries[3] = nx0;

  free(dom);
}
