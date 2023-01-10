#include <stdio.h>
#include "struct_variables.h"

/*
   But
   ===
   Modifier le vecteur b grâce à une fonction passée en argument de prob
   
   On évalue via la position dans le domaine, si le point se trouve dans le domaine du radiateur
  
  

  Arguments
  =========
  ix(input)    	 		Position horizontale dans le domaine
  iy(input)     		Position verticale dans le domaine
  pas_h(input)  		Pas de discrétisation
  rho (input)   		Valeur à imposer aux points dans le domaine du radiateur
  variables v (input)	structure contenant toutes les informations du problème



 */

double Rho_h(int ix,int iy,double invpas, double rho,struct variables v)
{
    /*On initialise la valeur à retourner à 0, en cas où le point n'est pas dans le domaine du radiateur*/
	double Rho=0;
	
	/*On considère un radiateur à une distance dist_rad du mur et !d'épaisseur unidirectionnelle! dans la direction orthogonale et opposée au mur*/
	if(iy>=v.dist_rad*invpas && iy<=(v.dist_rad+v.tickness)*invpas)	
	{
		if(ix>=v.r_h1*invpas && ix<=v.r_h2*invpas)
		{	
			/* Si dans le domaine du radiateur Rho différent de 0 et égal au rho imposé*/
			Rho=rho;	
		}	
	}
	/*On retourne la valeur qui va venir s'ajouter dans b*/
	return Rho;
}

//même fonctionnement que Rho_h
double Rho_v(int ix,int iy,double invpas, double rho,struct variables v)
{
    double Rho=0;
  
	if(ix<=(v.l-v.dist_rad)*invpas && ix>=(v.l-v.dist_rad-v.tickness)*invpas)
	{
		if(iy>=v.r_v1*invpas && iy<=v.r_v2*invpas)
		{	
			Rho=rho;	
		}		
	}
	return Rho;
}

