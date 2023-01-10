#include <stdio.h>
#include <stdlib.h>
#include "Rho.h"
#include "ecart_type.h"
#include "umfpk.h"
#include "math.h"
#include "plotting.h"
#include "struct_variables.h"
#include "modif_b.h"
#include "Calcul_Flux.h"


int Mesure_Rho(int m, int n, double *Rho,double *et_opt,double *flux, double I1, double I2,double tol, double (*pf)(int, int,double, double,struct variables v),int para, void *Numeric, int *ia, int *ja, double *a,struct variables v)
{
	/*
	 But
     ===
	 
	 Trouver le Rho optimal pour une certaine configuration du radiateur (horizontale, verticale et position)
	 Pour ce faire on utilise la méthode du nombre d'or qui permet de n'évaluer la fonction qu'une seule fois par itération
	 Les triplets ia, ja, a ansi que la factorisation LU (Numeric) sont passés en argument et ont déjà été évalués ce qui permet une réduction de 
	 coût de la résolution puisqu'il faut seulement modifier b en fonction de la valeur de Rho et résoudre un système Ax=b où la factorisation LU est déjà 
	 calculée
	 
	
	 
    Arguments
    =========
	
     m       (input)      	- nombre de points par direction dans la grille
     n       (input)      	- nombre d'inconnus dans le système
     Rho     (output)    	- pointeur vers le Rho_optimal qui sera retenu
     et_opt  (output) 	    - pointeur vers l'écart type optimal qui sera retenu
     flux    (ouput)        - pointeur vers le flux pour Rho_optimal 
     I1      (input)        - borne inférieure de l'intervalle de recherche
     I2      (input)        - borne supérieure de l'intervalle de recherche
     tol     (input)        - tolérance sur le résultat final 
     (*pf)   (input)        - pointeur vers une fonction pour modifier le b pour les configurations du radiateur
     para    (input)        - 0,1 ou 2 option d'affichage du résulat
     Numeric (input)        - Factorisation LU de la matrice A
     ia      (input)        - tableau 'ia' de la matrice A
     ja      (input)        - tableau 'ja' de la matrice A
     a       (input)        - tableau 'a' de la matrice A
	
	*/
	
	
	/* Déclaration des variables */
	
	int max_iter=50;
	
	/*tableau qui vont acceuilir b et les solutions des deux bornes de l'intervalle */
    double  *b,*x_I1, *x_I2;
    
    /* ET1, ET2 pour l'écart type et permettre la comparaison des bornes de l'intervalle */
    double ET1,ET2, pas_h=v.L/(m-1);
    
    
    /*Nombre d'or */
	double c=(-1+sqrt(5))/2;
	
	
	/* L'espacement des bornes d'évaluation est basé sur le nombre d'or */
	double x1=c*I1+(1-c)*I2, x2=(1-c)*I1+c*I2;

	/*INITIALISATION:
	 
	 C'est la seule fois où on devra évaluer 2 fois la fonction
	  
	 */
	

	
  /* Allocation mémoire pour b, x_I1 et x_I2 */
	b=malloc(n*sizeof(double));

	x_I1 = malloc(n*sizeof(double));
	x_I2 = malloc(n*sizeof(double));
	
	if(x_I1==NULL || x_I2==NULL)
	{
		printf("\n ERREUR : pas de mémoire pour les vecteurs b et x \n\n");
		return 1;
		 
    }
    
   /* Résolution pour la borne inférieure de l'intervalle */
	 if(modif_b(m,pas_h,&b,x1, pf,v))
		return 1;

	 if(solve_umfpack(n, ia, ja, a, b, x_I1,Numeric) )
			return 1;
	
	 ET1=ecart_type(x_I1, n);
	 
	
	
   /* Résolution pour la borne supérieure de l'intervalle */
	
	 if(modif_b(m,pas_h,&b,x2, pf,v))
		return 1;
		
	 if(solve_umfpack(n, ia, ja, a, b, x_I2,Numeric) )
		return 1;
			
	ET2=ecart_type(x_I2, n);
	
	//FIN INITIALISATION
	
  /*début recherche */
	int i=0;
	
	while(fabs(I1-I2)>tol && i<max_iter)
	{
	   
	
	   if(ET1 <ET2)
	   {
			I2=x2;
			x2=x1;
			ET2=ET1;
			x1=c*I1+(1-c)*I2;
			
			if (modif_b(m,pas_h,&b,x1, pf,v)) 
				return 1;
	        if(solve_umfpack(n, ia, ja, a, b, x_I1,Numeric) )
				return 1;
			ET1=ecart_type(x_I1, n);
		  
	   }
	   
	   else
	   {
		   I1=x1;
		   x1=x2;
		   ET1=ET2;
		   x2=(1-c)*I1+c*I2;
		   
		   if (modif_b(m,pas_h,&b,x2, pf,v)) 
			   return 1;
	       if(solve_umfpack(n, ia, ja, a, b, x_I2,Numeric) )
			   return 1;
		   ET2=ecart_type(x_I2, n);
	 
	   }
	 
	   /*nombre d'itérations*/
	   i++;
			
	}
	
	printf("\t%i itérations : le minimum se trouve entre [%lf, %lf]\n",i, I1, I2);
	
	/*recherche du min*/
	if(ET1<ET2)
	{
		*et_opt=ET1;
		*Rho=I1;
		printf("\tValeur optimale sous hypothèses: %lf\n",*Rho);
		printf("\tEcart Type de la solution : %lf\n", *et_opt);
		*flux=Calcul_Flux(m,pas_h,x_I1,v.f1,v.f2,v.Tf,v.k);
		printf("\n\n\n");
		plotting(x_I1,m,para,v);
		printf("\n\n\n");
	
	}
	else
	{
		*et_opt=ET2;
		*Rho=I2;
		printf("\tValeur optimale sous hypothèses: %lf\n",*Rho);
		printf("\tEcart Type de la solution : %lf\n", *et_opt);
		*flux=Calcul_Flux(m,pas_h,x_I2,v.f1,v.f2,v.Tf,v.k);
		printf("\n\n\n");
		plotting(x_I2,m,para,v);
		printf("\n\n\n");
	}

  /* libération mémoire */
   free(b);free(x_I1);free(x_I2);
    
    
  return 0;

}


