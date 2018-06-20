/*
  Helfrich energy material model
  author: Shiva Rudraraju
 */
#if !defined(HELFRICH_H_)
#define HELFRICH_H_

template <class T>
typedef struct {
  PetscReal kMean, kGaussian, mu; //material modulus
  PetscReal lambda;   //penalty parameter for enforcing incompressibility
  //model or history variables
  double H0;        //instantaneous curvature
  T q;              //Lagrange multiplier 
  //stress terms
  T sigma_contra[2][2], M_contra[2][2];
} HelfrichModel;

#undef  __FUNCT__
#define __FUNCT__ "computeStress"
template <class T>
PetscErrorCode computeStress(void* _kinematics, void* _materialModel)
{
  PetscFunctionBegin;

  KinematicsStruct<T> *k = (KinematicsStruct<T> *) _kinematics;
  HelfrichModel<T>    *m = (HelfrichModel<T> *) _materialModel;
  
  //contra-variant stress and bending moment tensors
  T (*sigma_contra)[2] = (T (*)[2]) m->sigma_contra;
  T (*moment_contra)[2] = (T (*)[2]) m->moment_contra;

  //model variables
  double K=m->kMean;
  double KGaussian=m->kGaussian;
  double mu=m->mu;
  double lambda=m->lambda;
  T dH=(k->H) - (m->H0);
  T Kappa=k->Kappa;
  T J=k->J;
  T q=m->q;
  
  //For Helfrich energy formulation
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
#ifdef LagrangeMultiplierMethod
      sigma_contra[i][j]=(q+ K*dH*dH - KGaussian*Kappa)*k->a_contra[i][j]-2*K*dH*k->b_contra[i][j];
#else
      sigma_contra[i][j]=(lambda*(J-1.0) + K*dH*dH - KGaussian*Kappa)*k->a_contra[i][j]-2*K*dH*k->b_contra[i][j];
#endif
      sigma_contra[i][j]+=(mu/(J*J))*(k->A_contra[i][j]-0.5*I1*k->a_contra[i][j]); //stabilization term
      moment_contra[i][j]=(K*dH + 2*KGaussian*H)*k->a_contra[i][j]-KGaussian*k->b_contra[i][j];
    }
  }
  
  PetscFunctionReturn(0);
}

#endif
