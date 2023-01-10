#include <math.h>
#include <stdio.h>

double ecart_type(double *x, int n)
{
	/*
	  But
	  ===
	  
	  Calcul de l'écart type. De facto, la température moyenne est calculée mais comme elle n'est pas utilisé on ne la passe pas en argument (peut être changer si besoin est)
	  Fonction de type double qui retourne l'écart-type
	 
	 Arguments
	 =======
	 
	 x (input)   -pointeur vers le vecteur solution
	 n (input)   -nombre d'inconnues
	 
	 */
	
	/*Déclaration des variables */	
	double S=0, ecart_t, moy;
	
	/*calcul de la moyenne */
	for(int i=0;i<n;i++)
    {
	   S=S+x[i];
	}
	
	moy=S/n;
	S=0;
	
	/* calcul de l'écart-type */
	for(int i=0;i<n;i++)
	{
		S=S+(x[i]-moy)*(x[i]-moy);
	}
	
	ecart_t=sqrt(S/n);
	
	/* on retourne l'écart-type */
	return ecart_t;
}
