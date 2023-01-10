#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "struct_variables.h"

int modif_b(int m,double pas_h, double **b,double iRho, double (*pf)(int, int,double, double,struct variables v),struct variables v)
{
/*
   But
   ===
   Permettre la modification de b sans regénération du triplet ia, ja, a qui sont constants pour certaines
   exécutions notamment la mesure du rho optimal

   La fonction est une version simplifiée de prob.c


  Arguments
  =========
  m     (input)    - nombre de points par direction dans la grille
  pas_h (input)    - pas de discrétisation
  n     (output)   - pointeur vers le nombre d'inconnus dans le système
  b     (output)   - pointeur vers le tableau 'b'
  iRho  (input)    - "imposed Rho", Rho à prendre en compte dans le pointeur de fonction pf
  				      et qui est lié à la puissance du radiateur
  *pf   (input)    - pointeur de fonction vers Rho_h ou Rho_v

*/

/*Déclaration variables */

	int passe=0;
    int ind;
    int m2=m;
    double invh2,alpha,invpas=(1/pas_h);

    invh2 = ((m-1)*(m-1))/(v.L* v.L); 
	alpha = (v.h/pas_h)+ 1;
   
    ind=0;
    
    
    for (int  iy = 0; iy < m; iy++)
	{
		if (iy == alpha)
		{
			//Le domaine peut être décomposé en 2 rectangles, le changement de m se fait quand on passe du rectangle primaire au secondaire

			m2 = (v.l / pas_h) + 1;
			
		}
	
		for (int ix = 0; ix < m2; ix++)
		{
			

			/*LISTE DES CONDITIONS
			On va utiliser un booléan 0, 1 pour voir si l'on est sur la porte ou la fenêtre.
			
			*/
			if (iy == 0 && (ix>= v.f1*invpas && ix <= v.f2*invpas)) //teste si l'on est sur la fenêtre
			{
				passe = 1;
			}
			if (ix == 0 && (iy>= v.p1*invpas && iy<=v.p2*invpas)) //teste si l'on est sur la porte
			{
				passe = 1;
			}

			//Début bloc principal
			if (passe == 0)
			{
				
				//appel des fonctions Rho_h ou Rho_v pour le vecteur b
				(*b)[ind] = (*pf)(ix, iy,invpas,iRho,v);
				

				//TRAITEMENT VOISIN SUD 
			 
				if (  (iy == m - 1  || (iy == alpha - 1 && ix> v.l*invpas)))
				{
					
					//Condition plus générale dans le cas où la porte se situe à l'extremité supérieure du mur gauche
					if ( (iy-1)==v.p2*invpas && ix==0)
					{
						(*b)[ind] += 2*v.Tp * invh2;
							
					}

				}
				else if (iy > 0)
				{

					if (iy == 1 && (ix>= v.f1*invpas && ix <= v.f2*invpas)) //au desssus de la fenêtre
					{
						(*b)[ind] += v.Tf * invh2;
						
					}

					else if ( ((iy-1)==v.p2*invpas) && ix==0) //Au dessus de la porte
					{
						(*b)[ind] += v.Tp * invh2;
						
					}
		
				}
	
				if (ix > 0)
				{
					/*On Est à Droite de la Fenêtre*/
					if (iy == 0 && (ix-1)==v.f2*invpas)
					{
						(*b)[ind] += v.Tf * invh2;
						
					}
					/*On est à droite de la porte */
					else if (ix==1 && (iy >=v.p1*invpas && iy<=v.p2*invpas))
					{
						(*b)[ind] += v.Tp * invh2;
					}
				

				}
				
				if (ix == 0 && (ix+1)==v.f1*invpas)
				{
					
						(*b)[ind] += 2*v.Tf * invh2;
					

				}

			    /*On Est à Gauche de la Fenêtre */
				if (iy == 0 && (ix + 1) == v.f1*invpas)  
				{
						(*b)[ind] += v.Tf * invh2;
						
				}

				/*On est En-Dessous de la porte*/
				if ( (iy + 1)== v.p1*invpas && ix==0)
				{
						(*b)[ind] += v.Tp * invh2;
						
				}
				
                /*Incrémentation de l'indice*/
				ind++;
			}
			
		    passe = 0;
			
		}
	
	}
	
    /* retour habituel de fonction */
    return 0;

}
