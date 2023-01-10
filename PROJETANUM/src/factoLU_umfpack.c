#include <stdio.h>
#include "umfpack.h"

int factoLU_umfpack(int n, int *ia, int *ja, double *a, void ** Numeric)
 {
   /*
    
     BUT
     ===
     
     Récupérer la factorisation LU pour optimiser le temps d'exécution lorsque seul b change, on peut donc réutiliser cette 
     facorisation à m constant pour b variant
     
     
     ARGUMENT
     ========
     
    n  (input)         -nombre d'inconnues
    ia (input)         - tableau 'ia' de la matrice A
    ja (input)         - tableau 'ja' de la matrice A
    a  (input)         - tableau 'a' de la matrice A
    Numeric (output)   - pointeur vers le type void *Numeric qui sera réutilisé pour résoudre le système (accès à la factorisation LU)
     
    */
    int status;
    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
    void *Symbolic;
    
    /* initialisation des paramètres par défaut */

    umfpack_di_defaults (Control) ;


    /* factorization symbolique */

    status = umfpack_di_symbolic (n, n, ia, ja, a, &Symbolic, Control, Info) ;
    if (status < 0)
    {
	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
        printf("\n ERREUR : umfpack_di_symbolic a échoué\n\n");
        return 1;
    }

    /* factorization symbolique - affichage */

    (void) umfpack_di_report_symbolic (Symbolic, Control) ;

    /* factorization numérique */

    status = umfpack_di_numeric (ia, ja, a, Symbolic, &(*Numeric), Control, Info);
    if (status < 0)
    {
	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
        printf("\n ERREUR : umfpack_di_numeric a échoué\n\n");
        return 1;
    }

    /* factorization numérique - affichage */

    (void) umfpack_di_report_numeric ((*Numeric), Control) ;
    
    
    return 0;
    
 }
