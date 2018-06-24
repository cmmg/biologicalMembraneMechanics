/*
  Residual implementation for Kirchhoff-Love shell formulation
  author: Shiva Rudraraju
 */
#if !defined(RESIDUAL_H_)
#define RESIDUAL_H_
#include "kinematics.h"
#include "HelfrichModel.h"

#define PI 3.14159

struct BVPStruct{
  IGA iga;
  //
  PetscReal l;       //characteristic length scale of the domain
  PetscReal epsilon; //penalty parameter for enforcing rotation BC's
  //
  PetscReal kMean, kGaussian, mu; //material modulus
  PetscReal lambda;   //penalty parameter for enforcing incompressibility
  PetscReal surfaceTensionAtBase;
  //
  PetscInt stabilization;
  //B.C's
  PetscInt type;
  PetscReal uDirichlet;
  Vec xDirichlet;
  bool angleConstraints[2];
  PetscReal angleConstraintValues[2];
  //load stepping
  PetscReal c_time;
  PetscInt load_increment;
  //Output
  bool isProc0;
  FILE * fileForUROutout;
};

#undef  __FUNCT__
#define __FUNCT__ "ResidualFunction"
template <class T>
PetscErrorCode ResidualFunction(IGAPoint p,
				PetscReal shift,const PetscScalar *V,
				PetscReal t,const T * U,
				PetscReal t0,const PetscScalar * U0,
				T *R,void *ctx)
{
  PetscFunctionBegin;
  BVPStruct *bvp = (BVPStruct *)ctx;
  //get number of shape functions (nen) and dof's
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //shape functions: value, grad, hess
  const PetscReal (*N) = (const PetscReal (*)) p->shape[0];
  const PetscReal (*N1)[2] = (const PetscReal (*)[2]) p->basis[1];
  const PetscReal (*N2)[2][2] = (const PetscReal (*)[2][2]) p->basis[2];
  
  //get kinematic quantities
  KinematicsStruct<T> k;
  getKinematics<T>(p, U, U0, k);
  
  //get stress
  HelfrichModel<T> m;
  m.bvp=bvp;
  computeStress(k, m);

  //Add stabilization
  PetscReal mu=bvp->mu;
  T J=k.J;
  T sigma_contra_inPlane[2][2], sigma_contra_outPlane[2][2];
  T sigma_contra_StabilizationTerm;
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      sigma_contra_inPlane[i][j]=m.sigma_contra[i][j];
      sigma_contra_outPlane[i][j]=m.sigma_contra[i][j];
      //stabization terms
      switch (bvp->stabilization) {
      case 0: //A method
	sigma_contra_StabilizationTerm=(mu/(J))*(k.A_contra[i][j]-k.a_contra[i][j]); 
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      case 1: //A-t method
	sigma_contra_StabilizationTerm=(mu/(J))*(k.A_contra[i][j]-k.a_contra[i][j]); 
	sigma_contra_inPlane[i][j]+=sigma_contra_StabilizationTerm;
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      case 2: //A-s method
	sigma_contra_StabilizationTerm=(mu/(J*J))*(k.A_contra[i][j]-0.5*k.I1*k.a_contra[i][j]); 
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      case 3: //A-st method
	sigma_contra_StabilizationTerm=(mu/(J*J))*(k.A_contra[i][j]-0.5*k.I1*k.a_contra[i][j]); 
	sigma_contra_inPlane[i][j]+=sigma_contra_StabilizationTerm;
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      case 4: //a method
	sigma_contra_StabilizationTerm=(mu/(J))*(k.aPre_contra[i][j]-k.a_contra[i][j]); 
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      case 5: //a-t method
	sigma_contra_StabilizationTerm=(mu/(J))*(k.aPre_contra[i][j]-k.a_contra[i][j]); 
	sigma_contra_inPlane[i][j]+=sigma_contra_StabilizationTerm;
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      case 6: //a-s method
	sigma_contra_StabilizationTerm=(mu/(k.JPre*k.JPre))*(k.aPre_contra[i][j]-0.5*k.I1Pre*k.a_contra[i][j]); 
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      case 7: //a-st method
	sigma_contra_StabilizationTerm=(mu/(k.JPre*k.JPre))*(k.aPre_contra[i][j]-0.5*k.I1Pre*k.a_contra[i][j]); 
	sigma_contra_inPlane[i][j]+=sigma_contra_StabilizationTerm;
	sigma_contra_outPlane[i][j]+=sigma_contra_StabilizationTerm; break;
      }
    }
  }
  
  //
  bool surfaceFlag=p->atboundary;
  //PetscReal *boundaryNormal = p->normal;
  PetscReal pCoords[3]; IGAPointFormGeomMap(p,pCoords);
  PetscReal tempR=std::sqrt(pCoords[0]*pCoords[0]+pCoords[2]*pCoords[2]);
  if (tempR==0){tempR=1.0;}
  PetscReal rVec[3]={pCoords[0]/tempR, 0, pCoords[2]/tempR}; //Radial vector used for applying radial forces

  //rotational constraints
  bool hasRotationalConstraint = false; 
  PetscReal theta=0.0;
  if (std::abs(pCoords[1])<1.0e-2*bvp->l){ //bottom surface
    if (bvp->angleConstraints[0]){theta=bvp->angleConstraintValues[0]; hasRotationalConstraint=true;}
  }
  else{ //top surface
    if (bvp->angleConstraints[1]){theta=bvp->angleConstraintValues[1]; hasRotationalConstraint=true;}
  }
  PetscReal nVec[3]={0,0,0};
  if (hasRotationalConstraint){
    if (theta==0){ nVec[1]=1.0; }
    else if (theta==90){ nVec[0]=-rVec[0];  nVec[2]=-rVec[2]; }
  }
  
  //Residual
  PetscReal L=bvp->l;
  PetscReal K=bvp->kMean;
  PetscReal Epsilon=bvp->epsilon;
  if (!surfaceFlag) {
    for (unsigned int n=0; n<(unsigned int)nen; n++) {
      for (unsigned int i=0; i<3; i++){
	T sigma_in=0.0, sigma_out=0.0, moment=0.0;
	for (unsigned int a=0; a<2; a++){
	  for (unsigned int b=0; b<2; b++){
	    //sigma*grad(Na)*dxdR*J
	    sigma_in+=sigma_contra_inPlane[a][b]*N1[n][a]*k.dxdR[i][b];
	    T temp_ab=0.0;
	    for (unsigned int j=0; j<2; j++){
	      temp_ab+=k.normal[j]*k.dxdR2[j][a][b];
	    }
	    sigma_in+=sigma_contra_inPlane[a][b]*N[n]*k.normal[i]*temp_ab; //in plane components
	    sigma_out+=-sigma_contra_outPlane[a][b]*N[n]*k.normal[i]*temp_ab; //out of plane components
	    //M*(hess(Na)-Gamma*grad(Na))*n*J
	    moment += m.moment_contra[a][b]*(N2[n][a][b])*k.normal[i];
	    for (unsigned int l=0; l<2; l++){
	      moment += -m.moment_contra[a][b]*(k.Gamma[l][a][b]*N1[n][l])*k.normal[i];
	    }
	  }
	}
	T Ru_i = ((L*L)/K)*(sigma_in+sigma_out)*J+(L/K)*moment*J;
	//
	/*
	bool isCollar=false;
	//collar implementation
	if (std::abs(pCoords[1]-CollarZ)<=CollarDepth) {isCollar=true;}
	if (isCollar){ //axis aligned along Y-axis
	  Ru_i+=((L*L*L)/K)*N[n]*(t*CollarForce)*(k.normal[i])*J;
	}
	*/
	R[n*dof+i] = Ru_i; 
      }
#ifdef LagrangeMultiplierMethod
      //Lagrange multiplier residual, J-1
      R[n*dof+3] = N[n]*(L*L/K)*(J-1.0)*J;
#endif
    }
  }
  else{
    if (hasRotationalConstraint){
      if (bvp->angleConstraints[0]==true){
	//std::cout << pCoords[1] << "v: (" << nVec[0] << ", " << nVec[1] << ", " << nVec[2] << ")\n";
	//std::cout << pCoords[1] << "n: (" << k.normal[0].val() << ", " << k.normal[1].val() << ", " << k.normal[2].val() << ")\n";
      }
    }
    for (unsigned int n=0; n<(unsigned int)nen; n++) {
      for (unsigned int i=0; i<3; i++){
	T Ru_i=0.0;
	//Rotational constraints
	//change of curve length not currently accounted as the corresponding Jacobian not yet implemented.
	//Ru_i+=epsilonBar*N1[n][a]*normal[i]*(nValue[j]-normal[j])*dxdR_contra[j][a]
	if (hasRotationalConstraint){
	  for (unsigned int j=0; j<3; j++){  
	    for (unsigned int a=0; a<2; a++){  
	      Ru_i+=-(L/K)*Epsilon*k.normal[i]*(k.normal[j]-nVec[j])*(N1[n][a]*k.dxdR_contra[j][a]);
	    }
	  }
	}
	Ru_i+=-((L*L)/K)*N[n]*rVec[i]*bvp->surfaceTensionAtBase;
	R[n*dof+i] = Ru_i; 
      }
#ifdef LagrangeMultiplierMethod
      R[n*dof+3] = 0.0;
#endif
    }
  }
  //
  PetscFunctionReturn(0);
}


#include "evaluators.h"
#endif
