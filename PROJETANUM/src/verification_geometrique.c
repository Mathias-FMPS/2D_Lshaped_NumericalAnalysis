#include <stdio.h>
#include <stdlib.h>
#include "struct_variables.h"

int verification_geometrique(struct variables v)
{
	
  /*
   
   BUT
   ===
   
   Vérifier que les données géométriques du problème sont cohérentes (liste de conditions non exhaustives)
   
   
   ARGUMENT
   ========
   
   v (input)    -structure contenant l'ensemble des données géométriques et physiques du problèmes.

*/
    double *p;  
    p=(double *)&v;
    int i;
    
    //On parcours tout les éléments de la structure (sauf les paramètres physiques Tp, Tf, k)
    for(i=0;i<13;i++){
		if(*(p+i)<0){ printf("\nERREUR: donnée(s) géométrique(s) négative(s)\n"); return 1;}
	    if(i>0){
			if(*(p) <= *(p+i)){ printf("\nERREUR: L doit être strictement supérieur à tout autre paramètre géométrique\n"); return 1;}} 
	}
	if(v.f1>v.f2) { printf("\nERREUR: f1 > f2\n"); return 1;}
	if(v.p1>v.p2) { printf("\nERREUR: p1 > p2\n"); return 1;}
	if(v.r_h1>v.r_h2){ printf("\nERREUR: r_h1 > r_h2\n"); return 1;}
	if(v.r_v1>v.r_v2) { printf("\nERREUR: r_v1 > r_v2\n"); return 1;}
	if(v.r_h2 > v.l && v.dist_rad+v.tickness >v.h) { printf("\nERREUR: radiateur horizontal en dehors du domaine\n"); return 1;}
		
		return 0;
}


