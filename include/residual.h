/*
  Residual implementation for Kirchhoff-Love shell formulation
  author: Shiva Rudraraju
 */
#if !defined(RESIDUAL_H_)
#define RESIDUAL_H_
#include "include/kinematics.h"
#include "include/HelfrichModel.h"

template <class T>
typedef struct {
  IGA iga;
  HelfrichModel<T> m;
  //
  PetscReal l;       //characteristic length scale of the domain
  PetscReal epsilon; //penalty parameter for enforcing rotation BC's 
} BVPStruct;

#undef  __FUNCT__
#define __FUNCT__ "Residual"
template <class T>
PetscErrorCode Residual(IGAPoint p,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * tempU,
			PetscReal t0,const PetscScalar * tempU0,
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
  
  //get X
  const PetscReal (*X)[3] = (const PetscReal (*)[3]) p->geometry;

  //get x, q
  T x[nen][3];
#ifdef LagrangeMultiplierMethod
  T (*u)[3+1] = (T (*)[3+1])tempU;
  T q=0.0; //q is the Lagrange multiplier value at this point
#else
  T (*u)[3] = (T (*)[3])tempU;
#endif
  for(unsigned int n=0; n<(unsigned int) nen; n++){
    for(unsigned int d=0; d<3; d++){
      x[n][d]= X[n][d]+ u[n][d];
    }
#ifdef LagrangeMultiplierMethod
    q+=N[n]*u[n][3];
#endif
  }
  
  //kinematics
  KinematicsStruct<T> k;

  //basis vectors, dXdR and dxdR, gradient of basis vectors, dXdR2 and dxdR2
  double (*dXdR)[2] = (double (*)[2]) &k.dXdR[0];
  double (*dXdR2)[2][2] = (double (*)[2][2]) &k.dXdR2[0];
  T (*dxdR)[2] = (double (*)[2]) &k.dxdR[0];
  T (*dxdR2)[2][2] = (double (*)[2]) &k.dxdR2[0];
  for(unsigned int d=0; d<3; d++){
    dXdR[d][0]=dXdR[d][1]=0.0;
    dXdR2[d][0][0]=dXdR2[d][0][1]=0.0;
    dXdR2[d][1][0]=dXdR2[d][1][1]=0.0;
    dxdR[d][0]=dxdR[d][1]=0.0;
    dxdR2[d][0][0]=dxdR2[d][0][1]=0.0;
    dxdR2[d][1][0]=dxdR2[d][1][1]=0.0;
    for(unsigned int n=0; n<(unsigned int) nen; n++){
      dXdR[d][0]+=N1[n][0]*X[n][d];
      dXdR[d][1]+=N1[n][1]*X[n][d];
      dXdR2[d][0][0]+=N2[n][0][0]*X[n][d];
      dXdR2[d][0][1]+=N2[n][0][1]*X[n][d];
      dXdR2[d][1][0]+=N2[n][1][0]*X[n][d];	    
      dXdR2[d][1][1]+=N2[n][1][1]*X[n][d];
      //
      dxdR[d][0]+=N1[n][0]*x[n][d];
      dxdR[d][1]+=N1[n][1]*x[n][d];
      dxdR2[d][0][0]+=N2[n][0][0]*x[n][d];
      dxdR2[d][0][1]+=N2[n][0][1]*x[n][d];
      dxdR2[d][1][0]+=N2[n][1][0]*x[n][d];	    
      dxdR2[d][1][1]+=N2[n][1][1]*x[n][d];
    }
  }

  //get all kinematic quantities
  getKinematics(&k);

  //get stress
  computeStress(&k, &bvp->m);
    
  //
  bool surfaceFlag=p->atboundary;
  PetscReal *boundaryNormal = p->normal;
  PetscReal pCoords[3]; IGAPointFormGeomMap(p,pCoords);
  
  //Residual
  PetscReal L=bvp->l;
  PetscReal K=bvp->m.KMean;
  T J=k.J;
  if (!surfaceFlag) {
    for (unsigned int n=0; n<(unsigned int)nen; n++) {
      for (unsigned int i=0; i<3; i++){
	T Ru_i=0.0;
	for (unsigned int j=0; j<2; j++){
	  for (unsigned int k=0; k<2; k++){
	    //sigma*grad(Na)*dxdR*J
	    Ru_i += (L*L/K)*sigma_contra[j][k]*N1[n][j]*dxdR[i][k]*J;
	    //M*(hess(Na)-Gamma*grad(Na))*n*J
	    Ru_i += (L/K)*M_contra[j][k]*(N2[n][j][k])*normal[i]*J;
	    for (unsigned int l=0; l<2; l++){
	      Ru_i += -(L/K)*M_contra[j][k]*(Gamma[l][j][k]*N1[n][l])*normal[i]*J;
	    }
	  }
	}
	//
	bool isCollar=false;
	//collar implementation
	if (std::abs(pCoords[1]-CollarZ)<=CollarDepth) {isCollar=true;}
	if (isCollar){ //axis aligned along Y-axis
	  Ru_i+=((L*L*L)/K)*N[n]*(t*CollarForce)*(normal[i])*J;
	}
	R[n*dof+i] = Ru_i; 
      }
#ifdef LagrangeMultiplierMethod
      //Lagrange multiplier residual, J-1
      R[n*dof+3] = N[n]*(L*L/K)*(J-1.0)*J;
#endif
    }
  }
  else{
    for (unsigned int n=0; n<(unsigned int)nen; n++) {
      for (unsigned int i=0; i<3; i++){
	T Ru_i=0.0;
	//Rotational constraints 
	for (unsigned int j=0; j<2; j++){
	  //change of curve length not currently accounted as the corresponding Jacobian not yet implemented.
	  //Ru_i+=epsilonBar*std::abs(nValue[0])*(N1[n][j]*normal[i]*nValue[d]*dxdR_contra[d][j]);
	  if (std::abs(pCoords[0])<1.0e-8){
	    PetscReal nValue[3]={0.0, 0.0, -1.0};  //normal along -Z for \eta_2=1
	    Ru_i+=-epsilonBar*normal[0]*normal[0]*(N1[n][j]*dxdR_contra[i][j]); 
	  }
	  else{
	    PetscReal nValue[3]={1.0, 0.0, 0.0}; //normal along +X for \eta_2=0
	    Ru_i+=-epsilonBar*normal[2]*normal[2]*(N1[n][j]*dxdR_contra[i][j]);
	  }  
	}
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


#include "include/evaluators.h"
#endif
