#include "prob.h"

int prob(int m, int *n, int **ia, int **ja, double **a, double **b){

/*
   But
   ===
   Génère la matrice n x n qui correspond à la disrétisation sur une grille
   cartesienne regulière m x m de l'operateur de Laplace à deux dimensions

            d    d        d    d
         - == ( == u ) - == ( == u )        sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec la fonction u qui satisfait les conditions aux limites de Dirichlet

         u = 0  sur (0,y), (1,y), (x,0) et (x,1), avec 0 <= x,y <= 1 .

  Arguments
  =========
  m (input)   - nombre de points par direction dans la grille
  n  (output) - pointeur vers le nombre d'inconnues dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A

*/

  int  nnz, nx, x0, y0, imax, imin, jmax, jmin;  /*nnz == n non-zero*/

  double *dom = malloc(4*sizeof(double));
  dom[0] = 0.0; dom[1] = 1.0 / 4.0; dom[2] = 5.0 / 8.0; dom[3] = 1.0 / 2.0;

  /*
  nx0 est la longueur suivant x du trou
  ny0 est la longueur suivant y du trou
  les imin, imax, jmin, jmax sont les indices maximum aux bords
  */
  long double invh2;


  //############################# Problème carré ##############################

    nx = m - 2; /* noeuds de Dirichlet ne sont pas pris en compte == bords */
    invh2 = ((m-1)*(m-1)); /* h^-2 pour L=1, le pas au carré */

    *n  = nx * nx; /* nombre d'inconnues pour le carré */
    nnz = 5 * nx * nx - 4 * nx;  /*prédiction du nbr d'élements non nuls de la matrice*/

  //############################## Problème adapté ##############################

  /*Avant d'allouer la mémoire, il s'agit de savoir quelle sera les inconnues à
  enlever avec le rectangle.

  Si on écrit le domaine de la forme [x1, y1]x[x2, y2], on adapte le domaine
  aux conventions lexicographiques*/

  double L = fabs(dom[2] - dom[0]);
  double l = fabs(dom[3]- dom[1]);

    // Si les sections ne sont pas parfaitement divisibles, il faut prendre
    // la meilleure approximation

  x0 = round(L*(m-1));
  y0 = round(l*(m-1));

  // On converti en points discrétisés, avec prise en compte des bords

  double nx0;
  double ny0;

  defineLine(dom[0] ,dom[2] ,dom[1] ,dom[3] ,x0 ,y0 ,&nx0 ,&ny0);


  /*
  Le nouveau n est le nouveau nombre d'inconnues,
  nnz = 5*nombre inconnues - nbr de zéros de Dirichlet - nbr de nouveaux zéros*/


  //################ Determination du nombre d'inconnues et de zéros ############

  *n -= nx0*ny0;  // On enlève les ponts sur la surface du bord
  nnz = 5*(*n) - 3*nx - ((nx - ny0) + ny0 + 2*nx0);


  // ############### Détermination des indices du "trou" ########################
  // On cherche les indices juste " à coté" du trou, pour déterminer les cas limites

  imin = ny0 - 3;
  imax = 2*ny0 - 2;
  jmin = 0;
  jmax = nx0;

  //printf("indices du trou imin  : %d imax : %d jmin : %d jmax : %d\n",imin, imax, jmin, jmax);

  /*Pour que l'algorithme fonctionne, il faut que les indices au dessus et en
  dessous du trou soient connus et référencés. Il faut 4 indices, imax, imin, jmax et jmin,
  respectivement pour la position relative du trou suivant y et suivant x*/


  //#################################################################################

  /* allocation des tableaux (calloc pour les zéros) */

  *ia  = malloc((*n + 1) * sizeof(int));
  *ja  = malloc(nnz * sizeof(int));
  *a   = malloc(nnz * sizeof(double));
  *b   = malloc(nnz * sizeof(double));   // Terme indépendant CB

  if(ia == NULL || ja == NULL || a == NULL || b == NULL){
    printf("Error : impossible to allocate memory prob.c -- mem allocation ia/ja/a/b\n");
    return 1;
  }


  // #############################################################################
  // ####################### Algorithme de remplissage ###########################
  // #############################################################################


  /* Le principe est de donner une équation numérotée par point de discrétisation, on
  parcours la matrice suivant iy (matrice ia), et à chaque ia on calcule chaque élément
  de la ligne par l'indice colonne (élement ja).

  Pour former la matrice du problème, le principe est alors d'ignorer les zones "à trous"
  c'est à dire vérifier si on est dans la zone en iy, et décaller le début du remplissage
  jusque la dimension en ix.
  */

  //################# Algorithme de remplissage : indices #######################

  int ind = 0;
  nnz = 0;

  for (int iy = 0; iy < nx; iy++) {
  int count = 0;
    for (int ix = 0; ix < nx; ix++) {
      // Conditions pour apparaitre dans la matrice
      if( !((iy > imin && iy < imax) &&  ix < jmax) ){
        /*On est en dehors du trou, on marque l'équation et on rempli*/
        (*ia)[ind] = nnz;

  //###################### Remplissage de la matrice #############################


        // Voisin sud ##########################################################
        if(iy > 0 && !((iy == imax) && (ix < jmax))){

          // -- > editer les indices d'équation, première zone
          if(iy <= imin){
            (*ja)[nnz] = ind - nx;
            (*a)[nnz] = -invh2;
          }

          // Deuxième zone
          else if(iy < imax && iy > imin){
            (*a)[nnz] = -invh2;
            (*ja)[nnz] =  ind - (nx - nx0);
          }

          // Passé le trou
          else {
            (*a)[nnz] = -invh2;
            (*ja)[nnz] = ind - nx;
          }
          nnz++;
        }

        // Voisin Ouest ########################################################
        if(ix > 0 && !(ix == jmax && (iy < imax && iy > imin))){
          (*a)[nnz] = -invh2;
          (*ja)[nnz] = ind - 1;
          nnz++;
        }

        // Element diagonal ####################################################
        (*a)[nnz] = 4.0*(invh2);
        (*ja)[nnz] = ind;

        // --  Non homogenious Dirichlet conditions-----------------------------

            // Superior border
            if((iy == imax) && ((ix < jmax) && (ix > jmin))){
              (*b)[ind] = 1.0;
            }

            // Inferior border
            if((iy == imin) && ((ix < jmax) && (ix > jmin))){
              (*b)[ind] = 1.0;
            }

            // Right hand side border
            if(ix == jmax && ( (iy < imax) && (iy > imin))){
              (*b)[ind] = 1.0;
            }
        nnz++;

        // Element Est #########################################################
        if(ix < nx - 1){
          (*a)[nnz] = -invh2;
          (*ja)[nnz] = ind + 1;
          nnz++;
        }

        // Element Nord ########################################################
        if((iy < nx - 1) && !(iy == imin && ix < jmax)) {

          // Première zone
          if(iy < imin){
            (*a)[nnz] = -invh2;
            (*ja)[nnz] = ind + nx;
          }

          // Au niveau du trou
          else if(iy >= imin && iy < imax - 1){
            (*a)[nnz] = -invh2;
            (*ja)[nnz] =  ind + (nx - nx0);
          }

          // Passé le trou
          else{
            (*a)[nnz] = -invh2;
            (*ja)[nnz] = ind + nx;
          }
          nnz++;
        }
        ind++;
      }
    }
  }

  /*Il reste à déterminer les derniers éléments des matrices*/
  (*ia)[ind] = nnz;

  free(dom);


return 0;
}

void defineLine(double x1,double x2,double y1,double y2, int x0, int y0, double *nx0, double *ny0){
  /*
    INPUT : relative coordinates of the hole in the grid
    OUTPUT : nx0 and ny0

    Desc : computes the amount of point on grid per lines
  */


  if(x1 == 0.0 || x2 == LS){
    *nx0 = x0; // On ne prend pas le bord confondu avec le coté gauche
    if(y1 == 0.0 || y2 == LS){
      *ny0 = y0;
    }
    else{
      *ny0 = y0 + 1;
    }
  }

  else if(y1 ==0.0 || y2 == LS){
    *ny0 = y0;
    if(x1 == 0.0 || x2 == LS){
      *nx0 = x0;
    }
    else{
      *nx0 = x0 + 1;
    }
  }

  else{
    *nx0 = x0 + 1;
    *ny0 = y0 + 1;
  }
}
