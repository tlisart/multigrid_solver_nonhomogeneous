#include "prolongation_mat.h"


void prolongationMatrix(int n, int n_c, int m, int m_c, int *bounds, int *bounds_c, double *u_c, double **u){
  /* INPUT : Dimensions, boundaries (geometry), starting vector u_c
     OUTPUT : reference u (fine prolongated vector)

     Desc : Interpolating the coarse problem solution to fit the extended dimension
            (n X n_c) matrix interpolating the coarse problem results. As the problem
            is in 2D --> bilinear interpolation for nodes not on coarse grid.

            Going throught the coarse grid as we do not need to remove non-existant
            points.
   */

   // Variables
   int x, y;
   double curr_value;
   int indC = 0;
   int nx_c = m_c - 2;

   int imin = bounds_c[0];
   int imax = bounds_c[1];
   int jmin = bounds_c[2];
   int jmax = bounds_c[3];

   for(int y_c = 0; y_c < nx_c; y_c++){
     for(int x_c = 0; x_c < nx_c; x_c++){

       // Condition to appear in the coarse matrix
       if( !((y_c > imin && y_c < imax) &&  x_c < jmax)){

         // Fine matrix coordinates
         x = 2*x_c + 1;
         y = 2*y_c + 1;

         curr_value = u_c[indC];

         // Interpolating on the fine grid

         (*u)[matching_fine(m, n, bounds, x , y, 0)] = curr_value;

         // Direct neighbours
         (*u)[matching_fine(m, n, bounds,x , y + 1, 0)] = curr_value*0.5;
         (*u)[matching_fine(m, n, bounds,x , y - 1, 0)] = curr_value*0.5;
         (*u)[matching_fine(m, n, bounds,x + 1, y, 0)] = curr_value*0.5;
         (*u)[matching_fine(m, n, bounds,x - 1, y, 0)] = curr_value*0.5;

         // Far neighbours
         (*u)[matching_fine(m, n, bounds,x + 1 , y + 1, 0)] = curr_value*0.25;
         (*u)[matching_fine(m, n, bounds,x + 1 , y - 1, 0)] = curr_value*0.25;
         (*u)[matching_fine(m, n, bounds,x - 1 , y + 1, 0)] = curr_value*0.25;
         (*u)[matching_fine(m, n, bounds,x - 1 , y - 1, 0)] = curr_value*0.25;

         indC ++;
       }
     }
   }
}
