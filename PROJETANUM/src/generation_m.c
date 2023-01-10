#include <stdio.h>

/*
	 But
     ===
	 
	 1)Génération_m --> générer les m qui menent à une bonne discrétisation du problème. Ces m sont aussi fonctions de la précision géometriques souhaitées 

    Arguments
    =========
	 
	 rang (input)        -correspond aux nombre de divisions du pas fondamental 
	 
     precision (input)   -précision géométrique désirée dans la modélisation du problème
   
   
  EX: L=4.5, précision 0.05
	 
  rang 1-->m=91 (pas_h=0.05)
  rang 2-->m=181 (pas_h=0.025)
  rang 5-->m=451 (pas_h=0.01)
*/

int generation_m(int rang,double precision, double L)
{
	
	return (L/precision)*(rang)+1;

}
