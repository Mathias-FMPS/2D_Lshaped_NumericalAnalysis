#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int calcul_residu(int n, int *ia, int *ja, double *a, double *b,double *x, double *residu)
{
/*
	 But
     ===
	                              ||Ax-b||
	 Calculer la norme du résidu  --------
	                                ||b||
	                               
    en minimisant le temps d'accès aux données en mémoire
    
    Arguments
    =========
   n  (input) - nombre d'inconnues
   ia (input) - tableau 'ia' de la matrice A
   ja (input) - tableau 'ja' de la matrice A
   a  (input) - tableau 'a' de la matrice A
   b  (input) - tableau 'b'
   x  (input) - tableau 'x', tableau de solution du système
   residu (output) - pointeur vers la valeur du résidu
	 
*/
	 
	double S1=0,S2=0, Ax=0;
	

    /*
     Les deux bouclent for suivantes permettent de calculer:
     
     
     1) Le produit Ax*b avec A en format CSR:
     
     Pour chaque élément i du vecteur résultant Ax (ici seulement une variable double et pas un vecteur), on effectue la multiplication a *x et la somme des a*x uniquement pour les éléments non nuls, 
     c-a-d les éléments définies par ia[i]-->ia[i+1]
     
     
     
     2) La norme ||Ax-b|| et ||b|| en réutilisant Ax pour stocker Ax-b
     
    */
    
    
    for(int i=0;i<n;i++)
    {
		
      for(int nnz=ia[i];nnz<ia[i+1];nnz++)
      {

         Ax+= a[nnz]*x[ja[nnz]];

      } 
      /* Ax[i]-b[i]*/
      Ax-=b[i];

      /*Norme de (Ax-b)*/
      S1+=Ax*Ax;

      /*Norme de b*/
      S2+=b[i]*b[i];
      
      /*Réinitialisation de Ax pour i suivan */
      Ax=0.0;
    }

    /*
     L'opération racine carré étant couteuse en ressource, on préfère l'utiliser 
     sur le quotient directement pour obtenir une petite optimisation (opération effectuée 1 fois au lieu de 2)
    */
    *residu= sqrt(S1/S2);
     
	return 0;

}
