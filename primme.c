#include <stdlib.h>
#include "primme.h"
#include "umfpk.h"
#include "main.h"


// Those needs to be static for it to work
static double *a;
static int n, m, type, *ia, *ja;
double resNorm = 10.0;


//void precond_primme(void *vx, void *vy, int *blockSize, primme_params *primme);
void matvec_primme(void *vx, void *vy, int *blockSize, primme_params *primme)

/*
   But
   ===
   Calcule le produit matrice-vecteur
                              vy = A*vx
   pour le solveur aux valeurs propres PRIMME. La matrice A doit être
   "stoquée" au préalable dans les variables statiques 'n', 'ia', 'ja' et 'a'
   en utilisant le format CSR (Compressed Sparse Rows). Par "stoquer"
   on veut dire ici stoquer la valeur de 'n' et les pointeurs vers les
   tableaux 'ia', 'ja' et 'a'.

   Arguments
   =========
   vx        (input) - vecteur(s) d'entrée
   vy       (output) - vecteur(s) de produit A*vx
   blockSize (input) - nombre de vecteurs d'entrée
   primme    (input) - paramètres fournis par primme pour optimiser le calcul
                       (pas utilisé)
*/

{
  int i, j, b;
  double *x = vx, *y=vy;

  for(b = 0; b < (*blockSize)*n; b+= n)
      for(i = 0; i < n; i++){
          y[b+i] = 0;
          for (j = ia[i]; j < ia[i + 1]; j++){
            y[b+i] += a[j] * x[b+ja[j]];
          }
      }
}

void precond_primme(void *vx, void *vy, int *blockSize, primme_params *primme)
/*
  Custom preconditioning method following matvec_primme conventions.
  Will use the two-grid method as preconditioner adapted to the square problem.
  The function is called by multiblock calling.

  As stated in the documentation, we search to approximate y = M^{-1}*x where the
  matrix is an approximate of A - \sigma*B

  */
  {
  double *x = vx, *y = vy;    //+
  int i, j;

  // Initial approximation
  for(i = 0; i < (*blockSize)*n; i += n){

    // Setting first approximate of u
    for(j = 0; j < n; j++){y[i + j] = 0.0;}

    if(type == 1)
      if(multiGridMethod(m, n, ia, ja, a, &x[i], &y[i], MAXIT_PRECOND, &resNorm, 1, 9, m)){printf("ERROR multigrid précond\n");}

    if(type == 2)
      if(multiGridMethod(m, n, ia, ja, a, &x[i], &y[i], MAXIT_PRECOND, &resNorm, 1, ((m - 1) / 2) + 1, m) ){printf("ERROR twogrid précond\n");}
  }
}

int primme(int primme_n, int primme_m, int *primme_ia, int *primme_ja, double *primme_a,
           int nev, double *evals, double *evecs, int typePrecond)

/*
   But
   ===
   Calcule les nve valeurs propres les plus basses de la matrice A stoquée
   dans le format CSR à l'aide du scalaire primme_n et des vecteurs
   primme_ia, primme_ja et primme_ia.

  Arguments
  =========
  primme_n  (input) - le nombre d'inconus dans le système
  primme_ia (input) - le tableau 'ia' de la matrice A
  primme_ja (input) - le tableau 'ja' de la matrice A
  primme_a  (input) - le tableau 'a' de la matrice A
  nev       (input) - le nombre de valeurs propres recherchées
  evals    (output) - le tableau des valeurs propres
  evecs    (output) - le tableau des vecteurs propres

  Retourne 0 si le caclul s'est bien déroulé, 1 si non.
*/
{
    int err;

    /* norme des résidus */
    double *resn = malloc(nev * sizeof(double));
    if (resn == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour un vecteur auxilière dans la fonction primme\n\n");
        return 1;
    }

    /* sauvgarder les pointeurs dans des variables statiques */
    n = primme_n;
    m = primme_m;
    a = primme_a;
    ja = primme_ja;
    ia = primme_ia;
    type = typePrecond;

    /* encoder les paramètres de PRIMME */
    primme_params primme;
    primme_initialize (&primme);
    primme.matrixMatvec = matvec_primme; /* MV product */
    // init the preconditionner

    primme.n = primme_n;         /* Matrix dimensions */
    primme.numEvals = nev; /* Number of wanted eigenpairs */
    primme.printLevel = 2; /* 1-4 */

    primme.applyPreconditioner = precond_primme;   // Using the same conventions than matrixMatvec...
    primme.correctionParams.precondition = 1;


    // Changing from the largest to the smallest
    primme.target = primme_smallest;

    //primme.target = primme_largest;
    //primme.target = primme_closest_geq;




    if(err = primme_set_method (DEFAULT_MIN_TIME, &primme)){
        printf("\nPRIMME: erreur N % dans le choix de la methode \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }

    /* afficher les papramètres de PRIMME */
    primme_display_params (primme);

    /* Caclul des valeurs et vecteurs propres */
    if(err = dprimme (evals, evecs, resn, &primme)){
        printf("\nPRIMME: erreur N %d dans le calcul des valeurs propres \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }

    /* libérer la mémoire */
    primme_Free (&primme); free(resn);

    return 0;
}
