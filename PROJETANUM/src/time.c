#include <stdio.h>
#include <time.h>
#include "time.h"

/* retourne le temps CPU depuis le début de l'execution du programme */
double mytimer_cpu(void){
	
	
    return (double) clock() /CLOCKS_PER_SEC;
}

#include <sys/time.h>
double mytimer_wall(void){
	
	struct timeval dummy;
	gettimeofday(&dummy,0);
	return (double) ((double) dummy.tv_sec + ((double) dummy.tv_usec)/1000000.0);
	
}
