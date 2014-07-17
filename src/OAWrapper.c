#include <stdlib.h>
#include "hungarian2.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
        

void OAWrapper(double* Rpre, int* m, int* n, int *mod, double *res);
  

// void R_init_GOSim(DllInfo *info){
//  /* Register routines, allocate resources. */
//  R_CMethodDef cMethods[] = {{"OAWrapper", &OAWrapper, -1,NULL},
// // 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP}}, 
//          {NULL, NULL, 0}};
//  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
// }
//           
// void R_unload_GOSim(DllInfo *info){        
// }

void OAWrapper(double* Rpre, int* m, int* n, int* mod, double* res){  
  hungarian_problem_t prob;   
  int i, j;
  double** R = (double**)R_alloc(*m,sizeof(double*));
  double maxweight = -999999;
  for(i = 0; i < *m; i++){
    R[i] = (double*)R_alloc(*n,sizeof(double));
    for(j = 0; j < *n; j++){
      R[i][j] = (Rpre[j* *m + i]);
      if(R[i][j] > maxweight)
        maxweight = R[i][j];
    } 
  }     
  int mode = (*mod == 1)? HUNGARIAN_MODE_MAXIMIZE_UTIL: HUNGARIAN_MODE_MINIMIZE_COST; 
  int dim = hungarian_init(&prob,R,*m,*n, mode);      
    double cost=hungarian_solve(&prob);   
    if(mode == HUNGARIAN_MODE_MAXIMIZE_UTIL)  
    cost = dim * maxweight - cost;  
  hungarian_free(&prob);  
  res[0] = cost;
}
