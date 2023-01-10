#include "struct_variables.h"

int Mesure_Rho(int m,int n,double *Rho,double *et_opt,double *flux,double I1, double I2,double tol, double (*pf)(int, int,double, double,struct variables v),int para, void *Numeric, int *ia, int *ja, double *a,struct variables v);
