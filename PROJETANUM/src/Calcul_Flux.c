#include <stdio.h>


double Calcul_Flux(int m,double pas_h, double *x,double f1, double f2, double Tf, const double k)
{ 

    
    /*	But
    ===
    Mesurer le flux sortant par la fenêtre. On utilise la loi de Fourier (J=-k*∇(T)) avec k, la conductivité thermique de l'air  [W*m⁻1*K⁻1] 
            (sa valeur est définie dans la structure variables)
    
    ∇(T) peut être développé comme ∇(T)=(T2-T1)/pas_h. Ici T2 est donné par un élément du vecteur x et T1=Tf
    
    un élément "infinitésimal" du flux dq(i) est donc donné par dq(i)=-k*(x[i]-Tf)/pas_h
    
    
    Pour obtenir le flux q, il suffit de sommer toutes les contributions infinitésimales dq(i) sur la longueur de la fênetre
    
    A ceci près, que dq(i) sera considéré constant sur toute la zone centré sur T(i) de pas_h si bien que le pas_h, dans la somme se simplifie.
    
    Illustration:
    
    
    T(i)----T(i+1)----T(i+2)---...----T(i+n)
       
     -->|<------->  |<------->
        |zone T(i+1)|zone T(i+2)|
        | pas_h     |pas_h
    <-------------------------------> longueur fenêtre
    
   NB : pour le premier et dernier point, la zone considéré est de longueur pas_h/2.
   
    Arguments
    =========
    
    m (input)      -nombre de points par direction dans la grille   
    pas_h (input)  -pas de discrétisation
    x (input)    -pointeur vers le vecteur de solution donné par la fonction Mesure_Rho() et qui correspond à la solution conduisant à la 
                   répartition de température la plus homogène
*/
   
   /* Nombre de points sur la fenêtre */
   double theta=(f2-f1)/pas_h+1;
    
 
   /* On se repère dans le vecteur solution pour trouver les solutions correspondant aux points juste au dessus de la porte */
	 int f = m-theta+f1/pas_h;
	 int g = m-1+f1/pas_h;


   /* q = variable qui donnera le flux */
    double q=0;
  
  /* Evaluation premier point et première zone voir NB*/
    q=q+(k*(x[f]-Tf)/2);
  
  /* Somme sur tous les éléments */
	for(int i=f+1;i<g;i++)
	{
	    q=q+k*(x[i]-Tf); 
	}

   /* Evaluation dernier point et dernière zone voir NB*/
	q=q+(k*(x[g]-Tf)/2);

	printf("\tLe flux linéique sortant par la fenêtre est: %lf W/m \n",q);

	return q;
}
