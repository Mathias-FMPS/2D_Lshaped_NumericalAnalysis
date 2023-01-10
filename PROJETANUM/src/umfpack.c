#include <stdio.h>
#include "umfpack.h"

int solve_umfpack(int n, int *ia, int *ja, double *a, double *b, double *x, void *Numeric)
{  
	
   /*
    
     BUT
     ===
     
     Résoudre le système Ax=b constitué de la matrice A en format CSR(triplet ia, ja, a), et du vecteur b. 
     
   Arguments
   =========
   n  (output)     - le nombre d'inconus dans le système
   ia (output)     - le tableau 'ia' de la matrice A
   ja (output)     - le tableau 'ja' de la matrice A
   a  (output)     - le tableau 'a' de la matrice A
   b  (output)     - le tableau des seconds membres B
   x  (output)     - le tableau des solutions X
   Numeric (input) - type void *, contient la factorisation LU de A
     
    */
	
	int status;
    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
   

    /* solution */
    /* status = umfpack_di_solve (UMFPACK_A, ia, ja, a, x, b, Numeric, Control, Info) ;*/
    status = umfpack_di_solve (UMFPACK_At, ia, ja, a, x, b, Numeric, Control, Info) ;
    /* UMFPACK utilise CSC et pas CSR =>  */
    umfpack_di_report_info (Control, Info) ;
    umfpack_di_report_status (Control, status) ;
    if (status < 0){
	printf("\n ERREUR : umfpack_di_solve a échoué\n\n");
        return 1;
    }
return 0;
}
