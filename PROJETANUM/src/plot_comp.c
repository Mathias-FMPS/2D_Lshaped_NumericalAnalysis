#include <stdio.h>
#include <stdlib.h>


int plot_comp()
{  
  /*
    But
    ===
    Afficher le graphique de comparaison entre le solveur UMFPACK et AGMG (comparaison du temps de résolution en fonction du pas de discrétisation et pour plusieurs tolérances soumises à AGMG))

*/ 
	/*Tableau de char contenant les instructions de plot */

	 char *GnuCommandesPlot[]= {"set title 'Comparaison Temps de Solution AGMG et UMFPACK'","set size 1,1","set logscale x","set logscale y","set xrange [10:1200]","set style data lines","set key left top","set border 3","set xtics nomirror","set ytics nomirror"
		 ,"set xlabel'Pas de Discrétisation-m'","set ylabel 'Temps de Resolution [sec]'"
		 ,"plot 'datafiles/data_comp.txt' u 1:2 lw 1 t'Umfpack',\'datafiles/data_comp.txt' u 1:3 lw 1 t'AGMG E-5',\'datafiles/data_comp.txt' u 1:4 lw 1 t 'AGMG E-10',\'datafiles/data_comp.txt' u 1:5 lw 1 t 'AGMG E-15'"};
	 
	

	 FILE*plot=fopen("datafiles/PlotInstruc.txt","w");
	 
	 for(int i=0;i<13;i++)
	 {
	 
	   fprintf(plot,"%s\n",GnuCommandesPlot[i]);
	 
	 }
	 fclose(plot);
	 
	 system("gnuplot -persistent datafiles/PlotInstruc.txt");
	  
	 return 0;

}
