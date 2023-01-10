#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void affichage(int taille,int m_adm[], double times_umfpack[], double residus_umfpack[], double times_AGMG[],  double residu_AGMG[],  int iter_AGMG[])
{ 
	/*
	 
	 BUT
	 ===
	 Simple affichage des résultats de la comparaison de UMFPACK et AGMG
	 
     ARGUMENT
     =======
     
     taille (input)      -Nombre de résolution et donc de comparaison
     m_adm[] (input)     -Les m utilisés pour les différentes résolution
     Tableaux[] (input)  -Différents tableaux contenants les résultats de comparaison en terme de temps de résolution, résidu, 
                          nombre d'itération (AGMG seulement)
	 
	 
	 */
	 
	
	 //En-tête
	printf("\n\n");
	printf("%*s",15,"");
	printf("%s","UMFPACK");
	printf("%*s",15,"AGMG(tol)");   
	printf("%*s",13,"E10-5");
	printf("%*s",15,"E10-10");   
	printf("%*s",16,"E10-15");  
	printf("\n\n");
	
	
	printf("----------------------------------------------------------------------------------\n");
	printf("%-*s",15,"");
	printf("%-*s\n",15,"residu");
	printf("%-*s",15,"");
	printf("%-*s\n",15,"temps(CPU)");
	printf("%-*s",15,"");
	printf("%-*s\n",15,"Iter(AGMG)");
	printf("----------------------------------------------------------------------------------\n");
	
	
	for(int ix=0;ix<taille;ix++)
	{ 
		
	   
	   printf("  m=%-*i",10,m_adm[ix]);
	   printf("%.3e",residus_umfpack[ix]);
	   printf("%*s",15,"");
	     
	   for(int iy=0;iy<3;iy++)
	   {
		   printf("%*.3e",15,residu_AGMG[3*ix+iy]);
	   }
	   printf("\n");
	   
	   printf("%*.3e",23,times_umfpack[ix]);
	   printf("%*s",15,"");
	   
	   for(int iy=0;iy<3;iy++)
	   {
	   
			printf("%*.3e",15,times_AGMG[3*ix+iy]);
	   
	   }
	   
	   printf("\n");
	   printf("%*s",37,"");
	   for(int iy=0;iy<3;iy++)
	   {
	   
			printf("%*i",16-iy,iter_AGMG[3*ix+iy]);
	   
	   }
	   printf("\n");
	   printf("----------------------------------------------------------------------------------\n");
	   printf("\n");
	
	
	}
	printf("\n");printf("\n");printf("\n");printf("\n");
	
	
}
