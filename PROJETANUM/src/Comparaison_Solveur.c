#include <stdlib.h>
#include <stdio.h>
#include "calcul_residu.h"
#include "prob.h"
#include "Rho.h"
#include "time.h"
#include "umfpk.h"
#include "affichage.h"
#include "factoLU_umfpack.h"
#include "generation_m.h"
#include "umfpack.h"
#include "struct_variables.h"
#include "modif_b.h"

void dagmg_(int*,double*,int*,int*,double*,double*,int*,int*,int*,int*,double*);

int Comparaison_Solveur(void)
{
  /*
	 But
     ===
	 
	 Comparer le solveur direct de systèmes linéaires pour matrices creuses (UMFPACK) 
	 au solveur itératif de systèmes linéaires pour matrices creuses (AGMG).
	 
	 La comparaison s'effectue de la manière suivante:
	 
	 1) Résolution pour m croissant
	 
	 2) AU sein d'un même m, résolution sous AGMG pour différentes tolérance (1.e-5, 1.e-10,1.e-15)
	 
	 Pour les deux solveurs: évaluation des résidus et du temps de solution (cpu)
	 
	 Finalement les données UMFPACK(résidus, temps) et AGMG(tol) (résidus, temps, nbr d'itération nécessaire pour converger) sont reprises dans tableau récapitulatif
	 + affichage graphique du temps de résolution des différents solveurs en fonction du pas de discrétisation m.
	 
    Si l'on souhaite obtenir un échantillon de comparaison plus petit ou plus grand, il suffit de modifier la variable taille
	 
   */ 
	 
	/*déclaration des variables */
	
    /* Initialisation des paramètres de AGMG (nbr max d'iter=100, pas d'affichage (iprint=-1), faire ijob, nrest)*/
	 int iter=100, iprint=-1, ijob=0,nrest=1;
	 int n, *ia, *ja, nnz;
	 /*Nombre de m différents à tester pour les deux solveurs*/
	 int taille=9;
	 /*On teste AGMG pour différents m mais aussi pour différentes tolérances, ici 3 */
	 int nbr_mesure_AGMG=3;
	 /*Déclaration du tableau contenant les m à tester*/
	 int m_adm [taille];
	 
	 double pas_h;
	 
	 /*Déclaration du tableau contenant les différents tolérances à tester pour AGMG*/
     double tol_[3]={1.e-5,1.e-10,1.e-15}; 
     /*a,b,x pour la résolution et t1, t2 pour la mesure du temps de solution*/
	 double *a, *b, *x, t1, t2;
	 /* déclaration des différents tableaux qui receuillent les résultats de la comparaison*/
	 double times_umfpack[taille],residus_umfpack[taille],times_AGMG[3*taille],residus_AGMG[3*taille];
	 int iter_AGMG[3*taille];
	 
	 /*Définition de données géométriques avec une précision 0.5 (sauf tickness, 0.25 mais sans incidence)*/
	 struct variables v2={4.5, 3, 1.5, 1, 3, 3, 4, 1.5, 3.5, 1, 3, 0.5, 0.25, 0.0, 20.0, 0.026};

	 /*Génération des m :
	  
	  On génère des m avec une précision de 0.5 pour avoir un maximum de m disponibles.
	  Pour avoir des m de différents ordres de grandeur en taille résolution,
	  on utilise la variable j pour passer les m similaires en terme d'ordre de grandeur
	  
	  L'échantillon ainsi obtenu s'étend de 10 à 586
	   
	  */
	 for(int i=0;i<taille;i++)
	 {
		 int j=i*i+1;
		  m_adm[i]=generation_m(j,0.5,v2.L);
	 }
	 
	 
     /*Résolution et comparaison pour les différents m*/
	 for(int ix=0;ix<taille;ix++)
	 { 
		 /*déclaration du pointeur type void pour la résolution UMFPACK*/
		   void *Numeric;
		   
		   /*Variable locale m */
		   int m=m_adm[ix];
		  
		  pas_h=v2.L/(m-1);

		   /*génération du système*/
		   
		   /*Génération de ia, ja, a et b  pour une configuration sans radiateur*/
	       if(prob(m, &n, &ia, &ja, &a, &b,0, Rho_h,v2))
				return 1;
	       
	       /*Allocation du vecteur solution x */
	       
	       x = malloc(n* sizeof(double)); 
             
           if(x == NULL)
           {
				printf("\n ERREUR : pas de mémoire pour vecteur des solutions avec m= %i\n Réduire la taille de l'échantillon\n",m_adm[ix]);
				return 1;
		   }
           
           nnz=ia[n];
             
             
	       /*Résolution par umfpack (mesure du temps de résolution + calcul résidu) */
	      
	       t1 = mytimer_cpu();
	       if(factoLU_umfpack(n,ia,ja,a,&Numeric))
				return 1;
	       if(solve_umfpack(n, ia, ja, a, b, x,Numeric))
				return 1;
	       t2 = mytimer_cpu(); 
	       times_umfpack[ix]=t2-t1;        
	       calcul_residu(n, ia, ja, a, b, x,&residus_umfpack[ix]); 
	       
	       /*Libérer la mémoire de la factorisation LU */
	       umfpack_di_free_numeric(&Numeric);
	      
	     /*résolution par AGMG */
	 
		   int indice;
	      
		   for(int i=0;i<nbr_mesure_AGMG;i++)
		   {      
			
			/*Redéfintion du nombre d'iter à chaque itération car iter = IN/OUT dans AGMG*/
				iter=100;
	     
			/*décalage d'indice Fortran */
	     
				for(int iy=0;iy<n+1;iy++)
				{
					ia[iy]++; ja[iy]++;
				}
				for(int iy=n+1;iy<nnz;iy++) 
				{
					ja[iy]++;
				}
			   
		  /*résolution par AGMG: mesure du temps de résolution et résidu */
		  
				t1 = mytimer_cpu();
				dagmg_(&n,a,ja,ia,b,x,&ijob,&iprint, &nrest,&iter,&tol_[i]);
				t2 =mytimer_cpu();   
				indice=3*ix+i;
				times_AGMG[indice]= t2-t1;
				iter_AGMG[indice]= iter;
	       
	      /*décalage d'indice Fortran après passage par AGMG */
	     
				for(int iy=0;iy<n+1;iy++)
				{
					ia[iy]--; ja[iy]--;
				}
				for(int iy=n+1;iy<nnz;iy++) 
				{
					ja[iy]--;
				}
		   
	      /*regénération du vecteur b car  b est un paramètre IN/OUT dans AGMG, utilisation de la fonction modif_b*/
	       
				if(modif_b(m,pas_h,&b,0,Rho_h,v2))
					return 1;
			
				calcul_residu(n, ia, ja, a, b, x,&residus_AGMG[indice]); 
     	       
        } 
        
        
         /*libération mémoire pour la résolution à taille m_i */
	   free(ia); free(ja); free(a); free(b);free(x);
    
      
	 }
	 
	 
	 /*Transcription des données dans un fichier texte 'data_comp.txt' en vu du plot*/
	 /*Cela aurait pu être fait dans la boucle précédente juste après obtention des données,
	  mais pour plus de clarté, on y consacre une autre boucle */
	 FILE *data=fopen("datafiles/data_comp.txt","w");
	 
	 for(int ix=0;ix<taille;ix++)
	 {
		
		  fprintf(data,"%i\t\t", m_adm[ix]);
		  fprintf(data,"%f\t\t",times_umfpack[ix]);

		  for(int i=0;i<nbr_mesure_AGMG;i++)
		  {
			  fprintf(data,"%f\t\t",times_AGMG[nbr_mesure_AGMG*ix+i]);
			  
		  }
		 
		fprintf(data,"\n");
	 }
	 fclose(data);
	  
	/*Affichage du tableau synthétisant les résultats*/
	  
	 affichage(taille,m_adm, times_umfpack,residus_umfpack,times_AGMG, residus_AGMG, iter_AGMG);

	 return 0;
	
}
