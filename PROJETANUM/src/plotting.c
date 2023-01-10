#include <stdio.h>
#include <stdlib.h>
#include "struct_variables.h"

void plotting(double *b, int m, int para, struct variables v)
{
	
	/*
    But
    ===
    Proposer un affichage graphique de la solution pour la distribution de la température 
    dans la pièce


    Arguments
    ========
    b(input)  -tableau de double qui contient la solution à afficher

   */ 
   char *GnuCommandesHeatMap []={"unset key","set title 'HeatMap: Température [°C] en fonction de la position dans la pièce [m,m]'","set xlabel 'Distance horizontale [m]'","set ylabel 'Distance verticale [m]'",
	"set cblabel 'Température [°C]'","set size ratio 1","set tics out"
   ,"set tics nomirror","set palette rgb 33,13,10","plot 'datafiles/data.txt' using 2:1:3 with image"};
   
   
   char *GnuCommandesDGraph [] ={"unset key","set title 'Température [°C] en fonction de la position dans la pièce [m,m]'","set xlabel 'Distance horizontale [m]'","set ylabel 'Distance verticale [m]'","set zlabel 'Température [°C]' rotate by 90",
   "set pm3d","set palette rgb 33,13,10","set cblabel 'Température [°C]'","splot 'datafiles/data.txt' using 2:1:3 with linespoints lt nodraw pt 0"};
   
   int ind=0;
   double pas_h= v.L/(m-1),invpas=(1/pas_h);
   
  /* creation d'un fichier data et écriture des données en vu du plot*/
   FILE *data = fopen("datafiles/data.txt", "w");

   
   for(int iy=0;iy<m;iy++) 
   {
  

      
      for(int ix=0;ix<m;ix++)
      {
         
         
        if(ix>v.l*invpas && iy>v.h*invpas)
		{
			 
			 fprintf(data,"%f\t%f\t%s\n",iy*pas_h,ix*pas_h,"NaN");
		}
		else
		{
         if(iy >=v.p1*invpas && iy<=v.p2*invpas &&ix==0)
         { 	
			fprintf(data,"%f\t%f\t%f\n",iy*pas_h,ix*pas_h,v.Tp);
         }
         else if(iy==0 &&(ix >=v.f1*invpas && ix <=v.f2*invpas))
         {
			  
			 fprintf(data,"%f\t%f\t%f\n",iy*pas_h,ix*pas_h,v.Tf);
		 }
         else
         {
			; 		
			fprintf(data,"%f\t%f\t%f\n",iy*pas_h,ix*pas_h,b[ind]);
          
           
			ind++;
          
          }
         }
         
      } 
      
      fprintf(data,"\n");
      
   } 
         
   fclose(data);
   

   /*Switch en fonction de para(int) pour différentes options d'affichage 
    
    para=0 --> HeatMap
    
    para=1 --> 3D
    
    
    para=2 --> Les deux
    
    
    */
    
   FILE *HeatMap;
   FILE *DGraph ;
   switch(para)
   {
	   
	 case 0:    HeatMap = fopen("datafiles/HeatMap.txt", "w");
				
				//Permet d'avoir un afichage correct pour tout L
				fprintf(HeatMap,"set xrange[0:%f]\n",v.L);
				fprintf(HeatMap,"set yrange[0:%f]\n",v.L);
				//
			    for(int i=0;i<10;i++)
			    {
					fprintf(HeatMap,"%s\n",GnuCommandesHeatMap[i]);
	
			    } 
     		
  
				fclose(HeatMap);
				system("gnuplot -persistent datafiles/HeatMap.txt");
				break;
					
					
	case 1: 	DGraph = fopen("datafiles/DGraph.txt","w");
				
				fprintf(DGraph,"set xrange[0:%f]\n",v.L);
				fprintf(DGraph,"set yrange[0:%f]\n",v.L);
				for(int i=0;i<9;i++)
				{
					fprintf(DGraph,"%s\n",GnuCommandesDGraph[i]);
	
				} 
   
				fclose(DGraph);
				system("gnuplot -persistent datafiles/DGraph.txt");
				break;
	   
    case 2:  	HeatMap = fopen("datafiles/HeatMap.txt", "w");
    
                fprintf(HeatMap,"set xrange[0:%f]\n",v.L);
				fprintf(HeatMap,"set yrange[0:%f]\n",v.L);
				for(int i=0;i<10;i++)
				{
						fprintf(HeatMap,"%s\n",GnuCommandesHeatMap[i]);
		
				} 
     		
  
				fclose(HeatMap);
   
  
				DGraph = fopen("datafiles/DGraph.txt","w");
				fprintf(DGraph,"set xrange[0:%f]\n",v.L);
				fprintf(DGraph,"set yrange[0:%f]\n",v.L);
				for(int i=0;i<9;i++)
				{
						fprintf(DGraph,"%s\n",GnuCommandesDGraph[i]);
	
				} 
   
				fclose(DGraph);
   
				// executer les commandes dans le fichier 
   
				system("gnuplot -persistent datafiles/HeatMap.txt");
				system("gnuplot -persistent datafiles/DGraph.txt");
				break;
  
   } 


  
}
