/*
  author: Shiva Rudraraju
 */
#if !defined(SOLVER_H_)
#define SOLVER_H_


//snes convegence test
PetscErrorCode SNESConverged_Interactive(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx){
  BVPStruct *bvp  = (BVPStruct*) ctx;
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);  
  //custom test
  if ((it>199) || (fnorm<1.0e-10)){
    *reason = SNES_CONVERGED_ITS;
    return(0);
  }

  //default test
  PetscFunctionReturn(SNESConvergedDefault(snes,it,xnorm,snorm,fnorm,reason,ctx));
}

#endif
