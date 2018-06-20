/*
  Evaluator functions for Residual and Jacobian (through Automatic Differentiation) computations
  author: Shiva Rudraraju
 */
#if !defined(EVALUATORS_H_)
#define EVALUATORS_H_

#undef  __FUNCT__
#define __FUNCT__ "Residual"
PetscErrorCode Residual(IGAPoint p,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx)
{
  PetscFunctionBegin;
  
  ResidualFunction<PetscReal>(p, shift, V, t, U, t0, U0, R, ctx);
  /*
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  std::vector<doubleAD> U_AD(nen*dof);
  for(int i=0; i<nen*dof; i++){
    U_AD[i]=U[i];
    U_AD[i].diff(i, dof*nen);
  }
  std::vector<doubleAD> tempR(nen*dof);
  Function<doubleAD> (p, shift, V, t, &U_AD[0], t0, U0, &tempR[0], ctx);
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      R[n1*dof+d1]= tempR[n1*dof+d1].val();
    }
  }
  */

  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Jacobian"
PetscErrorCode Jacobian(IGAPoint p,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscScalar *U,
			PetscReal t0,const PetscScalar *U0,
			PetscScalar *K,void *ctx)
{
  PetscFunctionBegin;

  AppCtx *user = (AppCtx *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  //
  std::vector<doubleAD> U_AD(nen*dof);
  for(int i=0; i<nen*dof; i++){
    U_AD[i]=U[i];
    U_AD[i].diff(i, dof*nen);
  } 
  std::vector<doubleAD> R(nen*dof);
  Function<doubleAD> (p, shift, V, t, &U_AD[0], t0, U0, &R[0], ctx);
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
      	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] = R[n1*dof+d1].dx(n2*dof+d2);
	}
      }
    }				
  }

  PetscFunctionReturn(0);
}

#endif
