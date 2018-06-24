/*
  Helfrich energy material model
  author: Shiva Rudraraju
 */
#if !defined(HELFRICH_H_)
#define HELFRICH_H_

struct BVPStruct;

template <class T>
struct HelfrichModel{
  BVPStruct *bvp;
  //stress and moment terms
  T sigma_contra[2][2], moment_contra[2][2];
};

#undef  __FUNCT__
#define __FUNCT__ "computeStress"
template <class T>
PetscErrorCode computeStress(KinematicsStruct<T>& k, HelfrichModel<T>& m)
{
  PetscFunctionBegin;
  
  //contra-variant stress and bending moment tensors
  T (*sigma_contra)[2] = (T (*)[2]) m.sigma_contra;
  T (*moment_contra)[2] = (T (*)[2]) m.moment_contra;

  //model variables
  double K=m.bvp->kMean;
  double KGaussian=m.bvp->kGaussian;
  double mu=m.bvp->mu;
  double lambda=m.bvp->lambda;
  T H=k.H;
  T dH=(k.H) - (k.H0);
  T Kappa=k.Kappa;
  T I1=k.I1;
  T J=k.J;
  T q=k.q;
  
  //For Helfrich energy formulation
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
#ifdef LagrangeMultiplierMethod
      sigma_contra[i][j]=(q+ K*dH*dH - KGaussian*Kappa)*k.a_contra[i][j]-2*K*dH*k.b_contra[i][j];
#else
      sigma_contra[i][j]=(lambda*(J-1.0) + K*dH*dH - KGaussian*Kappa)*k.a_contra[i][j]-2*K*dH*k.b_contra[i][j];
#endif
      moment_contra[i][j]=(K*dH + 2*KGaussian*H)*k.a_contra[i][j]-KGaussian*k.b_contra[i][j];
    }
  }
  
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "computeEnergy"
PetscErrorCode computeEnergy(KinematicsStruct<PetscScalar>& k, HelfrichModel<PetscScalar>& m, PetscScalar *energy){
  PetscFunctionBegin;
  
  //model variables
  double K=m.bvp->kMean;
  double KGaussian=m.bvp->kGaussian;
  double mu=m.bvp->mu;
  double lambda=m.bvp->lambda;
  double H=k.H;
  double dH=(k.H) - (k.H0);
  double Kappa=k.Kappa;
  double I1=k.I1;
  double J=k.J;
  double q=k.q;
  
  //For Helfrich energy formulation
  energy[0]=(K*dH*dH+KGaussian*Kappa)*J; //bending energy
#ifdef LagrangeMultiplierMethod
  energy[1]=q*(J-1)*J; //tension energy for Lagrange multiplier method
#else
  energy[1]=0.5*lambda*(J-1)*(J-1); //tension energy for penalty method
#endif
  PetscFunctionReturn(0);
}

#endif
