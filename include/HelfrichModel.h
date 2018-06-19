/*
  Helfrich energy material model
  author: Shiva Rudraraju
 */
#if !defined(HELFRICH_H_)
#define HELFRICH_H_

typedef struct {
  PetscReal l;
  PetscReal kMean, kGaussian, mu, epsilon, delta;
  //stress terms
  T sigma_contra[2][2], M_contra[2][2];
} HelfrichModel;

#undef  __FUNCT__
#define __FUNCT__ "computeStress"
template <class T>
PetscErrorCode computeStress(IGAPoint p,
			     PetscReal shift,const PetscScalar *V,
			     PetscReal t,const T * tempU,
			     PetscReal t0,const PetscScalar * tempU0,
			     T *R,void *ctx)
{
  PetscFunctionBegin;
   
  //contra-variant stress and bending moment tensors
  T sigma_contra[2][2], M_contra[2][2];
  //For Helfrich energy formulation
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
#ifdef LagrangeMultiplierMethod
      sigma_contra[i][j]=(q+ K*dH*dH - KGaussian*Kappa)*a_contra[i][j]-2*K*dH*b_contra[i][j];
#else
      sigma_contra[i][j]=(Lambda*(J-1.0) + K*dH*dH - KGaussian*Kappa)*a_contra[i][j]-2*K*dH*b_contra[i][j];
#endif
      sigma_contra[i][j]+=(Mu/(J*J))*(A_contra[i][j]-0.5*I1*a_contra[i][j]); //stabilization term
      M_contra[i][j]=(K*dH + 2*KGaussian*H)*a_contra[i][j]-KGaussian*b_contra[i][j];
    }
  }
  
  PetscFunctionReturn(0);
}

#endif
