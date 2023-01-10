#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prob.h"
#include "calcul_residu.h"
#include "plotting.h"
#include "umfpk.h"
#include "time.h"
#include "struct_variables.h"
#include "Rho.h"
#include "ecart_type.h"
#include "Comparaison_Solveur.h"
#include "plot_comp.h"
#include "Mesure_Rho.h"
#include "Calcul_Flux.h"
#include "factoLU_umfpack.h"
#include "generation_m.h"
#include "umfpack.h"
#include "verification_geometrique.h"




int main(int argc, char *argv[])
{
/*	But
    ===
    Appel des différentes fonction pour répondre aux Q1--Q5
    
    
    Q1,Q2,Q3: Résolution, mesure de temps de solution, calcul résidu, écart-type et affichage pour une pièce sans radiateur
    
    Q4: Résolution avec radiateur horizontal et vertical, détermination de la puissance optimale afin de donner un écart-type minimal (fct Mesure_Rho())
        Comparaison de l'ecart-type horizontal et vertical pour déterminer la configuration la plus homogène afin de calculer la puissance et la comparer au flux sortant par la fenêtre
        
        
    Q5: Comparaison solveur direct UMFPACK avec le solveur itératif AGMG
    
   	!!!!variables d'utilisations importantes : 
   
    rang --> définit m, plus le rang est élevé, plus m le sera.
    precision--> définit la précision géométrique des éléments dans v. Ex si p1=3.1 alors précision =0.1
    para --> option d'affichage
    
    NB: plus le degré de précision est fort plus des rangs faibles donneront des m élevés
  
    Ex pour L=4.5: precision=0.5, rang 1-->m=10
                 precision=0.01, rang 1-->m=451 
   
    par défaut la précision a été mise sur 10cm
    
*/
  /* déclaration des variables */
  
  /*Initialisation des variables géométriques dans une structure */
                        
    struct variables v={ 
	.L=4.5,
	.l=3,
	.h=1.5,
	.f1=1,
	.f2=3,
	.p1=3,
	.p2=4,
	.r_v1=1.5,
	.r_v2=3.5,
	.r_h1=1,
	.r_h2=3,
	.dist_rad=0.5,
	.tickness=0.25,
	.Tf=0.0,
	.Tp=20.0,
	.k=0.026};
   
    /*m: 
      n: nombre d'inconnus
   rang: variable intervenant dans la génération de m (multiplication par rang de la solution) doit être >0
   		 plus rang est élevé, plus m est élevé.
     */  
	int m, n,rang=11;
	int  *ia, *ja, para=0; //para (int) permet de choisir l'affichage (para=0: HeatMap, para=1: affichage 3D avec couverture en heatmap, para=2, les deux (x2 nbr de fenêtres))
	
    /*precision géométrique sur les données 
	 Ex: définiton au cm prés --> generation_m(rang, 0.01,v.L)*/
    double precision=0.1;

	double *a, *b, *x;

	/*Rho_optH et Rho_optV = résultat du calcul de la puissance optimale pour le radiateur */
	double t1, t2, Rho_optH,Rho_optV;
	/*et_h et et_v sont les écart-types pour les puissances optimales du radiateur
	flux_h et flux_v pour les valeurs de flux
	*/
	double et_h,et_v,flux_h,flux_v, res;
	double puissance;
  
	/*Permet de récupérer la factorisation LU faite par UMFPACK*/
	void *Numeric;
	
	
	/*Verification geometrique de la structure variable */
	if(verification_geometrique(v))
		return 1;
		
   /* Génération des m convenables en fonction de la tolérance sur les dispositions géométriques*/
	m=generation_m(rang,precision,v.L);
	
  /*Résolution pour Rho=0 (sans radiateur)*/
  
	if (prob(m,&n, &ia, &ja, &a, &b,0,Rho_h,v)) 
		return 1;
     
	printf("\nPROBLEM: ");
	printf("m = %5d   n = %8d  nnz = %9d\n\n", m, n, ia[n]);
	
  /* allouer la mémoire pour le vecteur de solution */

	x = malloc(n* sizeof(double));
  
	if(x == NULL)
	{
		printf("\n ERREUR : pas de mémoire pour vecteur des solutions 1\n\n");
		return 1;
	}
   
  /*résoudre et mesurer le temps de solution*/

	t1 = mytimer_cpu(); 
	if(factoLU_umfpack(n,ia,ja,a,&Numeric))
		return 1;
	if( solve_umfpack(n, ia, ja, a, b, x,Numeric) )
		return 1;   
	t2 = mytimer_cpu();
  
    /*Calcul du résidu*/
	if(calcul_residu(n, ia, ja, a, b, x,&res))
		return 1;
	
   /*Calcul écart type*/	
	et_h=ecart_type(x, n);
	
	printf("\tRESIDU: %e\n",res);
	printf("\tEcart-type de la solution %lf\n",et_h); 
	printf("\tTemps de solution (CPU): %lf sec\n",t2-t1);
	printf("\n\n");
	
	/*Affichage de la solution*/
	plotting(x,m,para,v);
	printf("\n\n\n\n\n");
	
  /*Libérer la mémoire alloué pour b et x*/
  
	 free(b); free(x);

  //==================Q4=======================//

  /*Calcul du Rho_opt pour radiateur horizontal */
	printf("Calcul du Rho optimal pour un radiateur horizontal\n\n");
	if(Mesure_Rho(m,n, &Rho_optH,&et_h,&flux_h,0,300,1,Rho_h,para,Numeric,ia, ja, a,v))
		return 1;


  /*Calcul du Rho_opt pour radiateur vertical */
    
    printf("Calcul du Rho optimal pour un radiateur vertical\n\n");
	if(Mesure_Rho(m,n,&Rho_optV,&et_v,&flux_v,0,10,0.5,Rho_v,para,Numeric, ia, ja, a,v))
		return 1;
	
	/*Libérer la mémoire alloué pour a, ia, ja et Numeric*/
	umfpack_di_free_numeric(&Numeric);
	free(a);free(ia);free(ja);	
		
  /*Determination du min*/
  
	printf("\n Position du radiateur optimale et puissance du radiateur:\n\n");
  
  /* Comparaison des résultats (écart type) pour les deux positions de radiateurs */
   
	if(et_h<et_v)
	{
		puissance=Rho_optH*v.k*(v.r_h2-v.r_h1)*v.tickness;
		printf("\tConfiguration Horizontale, puissance= %lf W/m\n", puissance);
		printf("\tDifférence de %.3lf %% (on néglige le flux de la porte)\n", (1-(puissance/flux_h))*100);
	
	}
	else
	{
		puissance=Rho_optV*v.k*(v.r_v2-v.r_v1)*v.tickness;
		printf("\tConfiguration Verticale, puissance= %lf W/m\n", puissance);
		printf("\tDifférence de %.3lf %% (on néglige le flux de la porte)\n",(1-(puissance/flux_v))*100);
	}
	
  //==================Q5=======================//
  
  /*Comparaison_Solveur et plot des résultats*/
    printf("\n\nComparaison du solveur direct UMFPACK et du solveur itératif AGMG\n\n");
	Comparaison_Solveur();
	plot_comp();

  /*retour normal de fonction */
	return 0;
  
}

