#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "struct_variables.h"


int prob(int m, int *n, int **ia, int **ja, double **a, double **b,double iRho, double (*pf)(int, int,double, double,struct variables v),struct variables v)
/*
   But
   ===
   Générer le système linéaire n x n 
                          
                             Au = b                                   

   qui correspond à la discrétisation sur une grille cartésienne 
   régulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )  = 0     sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         
         u = 20     sur (0,y)                 , avec 0 <=  y  <= 1
         du/dn = 1  sur (1,y), (x,0) et (x,1) , avec 0 <= x,y <= 1 .

  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice A est retournée dans le format 
  CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.
  
  

  Arguments
  =========
  m  (input)  - nombre de points par direction dans la grille
                (les valeurs m inférieures à 2 ne sont pas valides) 
  n  (output) - pointeur vers le nombre d'inconnus dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  b  (output) - pointeur vers le tableau 'b'

*/
{
	/* Déclaration des variables:
	 
	 
		1) nnz (int) -->  (nombre non nul)
		
		2) ind (int) --> Détermine l'indice de l'équation associée à un point
		
		3) var (int) --> Variable qui va modifier l'indice en fonction des spécificités des points, 
		 
		                 va permettre de traiter des conditions simultanées indépendamment les unes des autres.
		4) m2 (int)  --> Va permettre de modifier la valeur maximale de la boucle for en fonction de la largueur de la pièce
		                 (passage du domaine rectangulaire primaire au secondaire)
		
		5) invh2 (double) -->
		
		6)Tp, Tf (float) --> Température Porte et Fenêtre
		
		7) alpha, beta(double) --> Variable fonction de la géométrie de la pièce intervenant dans le calcul de n
		                   alpha=nbr de points nécessaires dans la direction iy pour changer la largueur du rectangle
		                   beta= nbr de points  dans la direction ix du rectangle supérieur
	     
	    8) Theta, Lambda (double) --> Nbre de points sur la fenêtre et sur la porte
	    
	    utilisation de double pour éviter les erreurs d'arrondies dans le calcul mais les valeurs seront entières (m calibré pour qu'elles le soient)

        9)pas_h et invpas (double) --> pas de discrétisation et son inverse utilisés dans les calculs
	    

	    On désignera "rectangle supếrieur" et "rectangle inférieur" de la manière suivante:
	     ______
	    |     |--> rectangle supérieur
	    |     |
        |_____|__ 	    
	    |       | --> rectangle inférieur
	    |_______|
	 
*/
	 
	double alpha, beta,theta,lambda;
	int  passe=0;
    int  nnz=0,var=0;
    int m2=m;
    double invh2=0, pas_h=v.L/(m-1), invpas=(1/pas_h);
    // Nbre de pas de discrétisation à faire avant d'être à hauteur du changement de largueur de la pièce
    alpha = (v.h/pas_h)+ 1;
	beta = (v.l/pas_h)+ 1;
	theta=(v.f2-v.f1)/pas_h +1;
	lambda=(v.p2-v.p1)/pas_h+1;
	
    
    invh2 = ((m-1)*(m-1))/(v.L* v.L); /* h^-2 pour L=1 */
    *n  =alpha * m + (m - alpha) * beta - (theta+lambda); /* nombre d'inconnues */
    nnz =  5*(*n)-4*m-4; /* nombre d'éléments non nuls */


    /* allocation des tableaux */

    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));

    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }



    /* partie principale : remplissage de la matrice */

   /*Initialisation de l'indice et du nnz qui sont incrémentés */
   int ind=0;
   nnz=0;
   
   for (int  iy = 0; iy < m; iy++)
	{
	
		if(iy==alpha)
		{
			/*Le domaine peut être décomposé en 2 rectangles, le changement de m se fait quand on passe du rectangle inférieur au supérieur*/

			m2 = beta;
			
		}
	
		for (int ix = 0; ix < m2; ix++)
		{
              
			/*LISTE DES CONDITIONS
			On va utiliser un booléan 0, 1 pour voir si l'on est sur la porte ou la fenêtre.*/
			
			if (iy == 0 && (ix>=v.f1*invpas && ix<= v.f2*invpas)) //teste si l'on est sur la fenêtre //OK
			{
				passe = 1;
			}
			if (ix == 0 && (iy>=v.p1*invpas && iy<= v.p2*invpas)) //teste si l'on est sur la porte
			{
				passe = 1;
			}

			/*Début bloc principal*/
			
			if (passe == 0)
			{
				(*ia)[ind] = nnz;
				
				/*appel des fonctions Rho_h ou Rho_v pour le vecteur b*/
				
				(*b)[ind] = (*pf)(ix, iy,invpas,iRho,v);


				/*TRAITEMENT VOISIN SUD*/
			
                
				if (  (iy == m - 1  || (iy == alpha - 1 && ix > v.l*invpas)))
				{
					
					/*porte se situe à l'extremité supérieure du mur gauche*/
					
					if ( (iy-1)==v.p2*invpas && ix==0) //à peut être modifié
					{
						(*b)[ind] += 2*v.Tp * invh2;		
					}
					
					/*Si l'on est sur l'extremité supérieure de la pièce y compris dans le rectangle inférieure*/
					 
					else
					{
						(*a)[nnz] = -2 * invh2;
						(*ja)[nnz] = ind - m2;
						nnz++;
					}

				}
				
				else if (iy > 0)
				{
                    /*Au-dessus de la fenêtre*/
                    
					if (iy == 1 && (ix>=v.f1*invpas && ix<= v.f2*invpas)) 
					{
						(*b)[ind] += v.Tf * invh2;
					
						
					}
					
					/*Au-dessus de la porte*/
					
					else if ((iy-1)==v.p2*invpas && ix==0) 
					{
						(*b)[ind] += v.Tp * invh2;
	
					}

					else
					{
						
						(*a)[nnz] = -invh2;
						
						/*l'utilisation de var permet de traiter plusieurs conditions indépendamment les unes des autres ex:
						 ex: diminution de la largueur de la pièce + porte */
						var = ind - m2;

						if (iy == 1 && (ix < v.f1*invpas)) // au-dessus mais à gauche de la fenêtre
						{
							var = var + theta;
						
						}

						if (iy == alpha)  //Changement dans la manière de compter car diminution de la largueur de la pièce
						{
							var = ind-m;
							
						}
						if (iy >= v.p1*invpas && iy<= v.p2*invpas) //Changement dans la manière de compter car points en moins dû à la porte
						{
							var++;
						}

						(*ja)[nnz] = var;		
						nnz++;

					}
				}

				/*TRAITEMENT VOISIN OUEST*/

				/*Condition de bord Est, m2 varie en fonction de la largeur du rectangle (inférieur ou supérieur)*/
				
				if (ix == m2 - 1) 
				{
					(*a)[nnz] = -2 * invh2;
					(*ja)[nnz] = ind - 1;
					nnz++;
				}
	
				else if (ix > 0)
				{
					/*On Est à Droite de la fenêtre (lambda occurence)*/
					
					if (iy == 0 && (ix - 1)==v.f2*invpas)
					{
						(*b)[ind] += v.Tf * invh2;
						
					}
					//On est à droite de la porte 
					else if (ix==1 && (iy>=v.p1*invpas && iy<=v.p2*invpas))
					{
						
						(*b)[ind] += v.Tp * invh2;
					}
					
					else
					{
						
						(*a)[nnz] = -invh2;
						(*ja)[nnz] = ind - 1;
						nnz++;
					}

				}

				/*TRAITEMENT ELEMENT DIAG*/

				(*a)[nnz] = 4.0 * invh2;
				(*ja)[nnz] = ind;
				nnz++;



				/*TRAITEMENT VOISIN EST*/
				//Condition Bord Ouest
				if (ix == 0)
				{
					(*a)[nnz] = -2 * invh2;
					(*ja)[nnz] = ind + 1;
					nnz++;
				}
				
				/* On traite tous les cas sauf celui du Bord Est */
				
				else if (ix < m2 - 1)
				{
					/*On Est à Gauche de la Fenêtre*/
					
					if (iy == 0 && (ix + 1)==v.f1*invpas)  
					{
						(*b)[ind] += v.Tf * invh2;
						
					}
					else
					{
						(*a)[nnz] = -invh2;
						(*ja)[nnz] = ind + 1;
						nnz++;
					}

				}

				/*TRAITEMENT VOISIN NORD*/
				
				if (iy == 0)
				{
					
					(*a)[nnz] = -2 * invh2;
					var = ind + m2;
					
					/*On Est à Gauche de la Fenêtre-->entraine un décalage d'indice de theta pour le voisin supếrieur*/
					
					if (ix < v.f1*invpas)
					{
						var = var - theta; 
						
					}

					(*ja)[nnz] = var;
					nnz++;

				}
				
				/*On traite les cas différents de Bord Nord, Intersection du rectangle primaire et secondaire ou ...*/
				
				else if (iy < m - 1 && (iy!=alpha-1  || ix<=v.l*invpas))
				{
					(*a)[nnz] = -1 * invh2;
					var = ind + m2;

				    //On est En-Dessous de la porte
						if ((iy+1)==v.p1*invpas && ix==0)
					{
						(*b)[ind] += v.Tp * invh2;
						
					}
					
					//On est dans la zone de hauteur de la porte--> provoque une modification de l'indice de l'équation du voisin nord
					else if (  (iy + 1) >= v.p1*invpas && iy<v.p2*invpas)
					{
						//On retranche 1 à var à cause des points déjà connus de la porte
						var = var - 1;
						(*ja)[nnz] = var;
						nnz++;
					}
					/*Cas général sans spécifité d'indice ou du vecteur b*/
					
					else
					{
						
						(*ja)[nnz] = var;
						nnz++;
					
					}
					
				}
						
                /*Incrémentation de l'indice*/
				ind++;

			}


			passe = 0;		

		}
		
	}

	(*ia)[ind] = nnz;  //comme on a utilisé une incrémentation de l'indice, on modifie le ind-1 en ind

    /*retour habituel de fonction */
    return 0;
}




