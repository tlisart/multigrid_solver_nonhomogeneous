#include "reduction_mat.h"

void residualTransfert(int n, int n_c, int m, int m_c, int *bounds, int *bounds_c, double *r, double **r_c){
  /* INPUT : Dimensions, boundaries (geometry), starting vector u
     OUTPUT : reference u_c (coarse reduced vector)

     Desc : restricting the residual from r_f to r_c, using the R (restriction) matrix
            (nc x n). Using full weighting  R_fw = average(matching points + neighbour
            not present on the coarse grid).

            Scaling factor changes in function of the relative position, the scaling is
            0.25 the matching point
            0.125 the neighbouring points
            0.0625 the external neighbouring points
   */

   // geometry of the coarse matrix
   int nx_c = m_c - 2;

   int imin = bounds_c[0];
   int imax = bounds_c[1];
   int jmin = bounds_c[2];
   int jmax = bounds_c[3];

   int x, y;

   // Running through the coarse grid and finding matching points on the fine matrix
   int ind_c = 0;

   for (int y_c = 0; y_c < nx_c; y_c++) {
       for (int x_c = 0; x_c < nx_c; x_c++) {

         // Condition to appear in the coarse matrix
         if( !((y_c > imin && y_c < imax) &&  x_c < jmax) ){

           // Fine matrix coordinates
           x = 2*x_c + 1;
           y = 2*y_c + 1;

           // Exact matching on the fine grid

           (*r_c)[ind_c] = 0.25 * r[matching_fine(m, n, bounds, x, y, 0)];

           // Close neighbours
           (*r_c)[ind_c] += 0.125 * r[matching_fine(m, n, bounds, x + 1, y, 0)];
           (*r_c)[ind_c] += 0.125 * r[matching_fine(m, n, bounds, x, y + 1, 0)];
           (*r_c)[ind_c] += 0.125 * r[matching_fine(m, n, bounds, x - 1, y, 0)];
           (*r_c)[ind_c] += 0.125 * r[matching_fine(m, n, bounds, x, y - 1, 0)];

           // Extended neighbours
           (*r_c)[ind_c] += 0.0625 * r[matching_fine(m, n, bounds, x - 1, y + 1, 0)];
           (*r_c)[ind_c] += 0.0625 * r[matching_fine(m, n, bounds, x - 1, y - 1, 0)];
           (*r_c)[ind_c] += 0.0625 * r[matching_fine(m, n, bounds, x + 1, y + 1, 0)];
           (*r_c)[ind_c] += 0.0625 * r[matching_fine(m, n, bounds, x + 1, y - 1, 0)];


           //printf("A l'indice COARSE : %d -- A l'indice FINE : %d\n",ind_c, matching_fine(m, n, bounds, x, y));
           ind_c++;
       }
     }
   }
 }



  int matching_fine(int m, int n, int *bounds, int x, int y, int version){
    /* INPUT : boundaries (discretized domain info), sizes
       OUTPUT : CSR format / indice of the point (x, y) on the fine matrix

       DESC : crutial function associating the ind on the 1D vector on a grid
              of any discretization giving the coordinates of the point.

              Re-use of the geometrical data in bounds and the same indice methods
              from prob.c
    */

    // geometry of the coarse matrix
    int ind = 0;

    if(version == 0){
      int nx = m - 2;

      int imin = bounds[0];
      int imax = bounds[1];
      int jmin = bounds[2];
      int jmax = bounds[3];

      int indSafe = 0;

      if(y <= imin){
        ind = y*nx + x;
      }else if(y > imin && y < imax){
        ind = (y - imin)*(nx - jmax) + imin*nx + x;
      }else if(y >= imax){
        ind = imin*nx +(imax - imin)*(nx - jmax) + (y - imax)*nx + jmax + x;
      }
    }

    else if(version == 1){
      int nx = m - 2;
      ind = y*nx + x;
    }
    return ind;
  }
