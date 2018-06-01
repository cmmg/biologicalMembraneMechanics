/*
  Kirchhoff-Love shell implementation in PetIGA
  author: rudraa
 */
#include <math.h> 
//extern "C" {
#include "petiga.h"
//}

//include automatic differentiation library
#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> doubleAD;

typedef struct {
  IGA iga;
  PetscReal l;
  PetscReal kMean, kGaussian, delta, mu, epsilon;
} AppCtx;

#undef  __FUNCT__
#define __FUNCT__ "Function"
template <class T>
PetscErrorCode Function(IGAPoint p,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * tempU,
			PetscReal t0,const PetscScalar * tempU0,
			T *R,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;
  PetscReal L=user->l; //Length scale for normalization
  PetscReal K=user->kMean; //bending modulus
  PetscReal KGaussian=user->kGaussian; //bending modulus
  PetscReal Lambda=user->delta; //penalty parameter for incompressibility
  PetscReal Mu=user->mu; //stabilzation parameter
  PetscReal Epsilon=user->epsilon;
  
  //normalization
  PetscReal kBar=K/K, kGaussianBar=KGaussian/K, lambdaBar=Lambda*(L*L)/K;
  PetscReal muBar=Mu*(L*L)/K, epsilonBar=Epsilon*(L/K); 
  
  //instantaneous curvature
  double H0=0.0; 

  //get number of shape functions (nen) and dof's
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //shape functions: value, grad, hess
  const PetscReal (*N) = (const PetscReal (*)) p->shape[0];
  const PetscReal (*N1)[2] = (const PetscReal (*)[2]) p->basis[1];
  const PetscReal (*N2)[2][2] = (const PetscReal (*)[2][2]) p->basis[2];
  
  //get X
  const PetscReal (*X)[3] = (const PetscReal (*)[3]) p->geometry;

  //get x
  T x[nen][3];
  T (*u)[3] = (T (*)[3])tempU;
  for(unsigned int n=0; n<(unsigned int) nen; n++){
    for(unsigned int d=0; d<3; d++){
      x[n][d]= X[n][d]+ u[n][d];
      //normalization
      //X[n][d]*=1.0/L;
      //x[n][d]*=1.0/L;
    }
  }

  //compute basis vectors, dxdR and dXdR, gradient of basis vectors, dxdr2, and co-variant metric tensors, a and A. 
  T dxdR[3][2], dxdR2[3][2][2], a[2][2];
  double A[2][2], dXdR[3][2];
  for(unsigned int d=0; d<3; d++){
    dXdR[d][0]=dXdR[d][1]=0.0;
    dxdR[d][0]=dxdR[d][1]=0.0;
    dxdR2[d][0][0]=dxdR2[d][0][1]=0.0;
    dxdR2[d][1][0]=dxdR2[d][1][1]=0.0;
    for(unsigned int n=0; n<(unsigned int) nen; n++){
      dXdR[d][0]+=N1[n][0]*X[n][d];
      dXdR[d][1]+=N1[n][1]*X[n][d];	 
      dxdR[d][0]+=N1[n][0]*x[n][d];
      dxdR[d][1]+=N1[n][1]*x[n][d];
      dxdR2[d][0][0]+=N2[n][0][0]*x[n][d];
      dxdR2[d][0][1]+=N2[n][0][1]*x[n][d];
      dxdR2[d][1][0]+=N2[n][1][0]*x[n][d];	    
      dxdR2[d][1][1]+=N2[n][1][1]*x[n][d];
    }
  }
  //
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      A[i][j]=0.0;
      a[i][j]=0.0;
      for(unsigned int d=0; d<3; d++){
	A[i][j]+=dXdR[d][i]*dXdR[d][j];
	a[i][j]+=dxdR[d][i]*dxdR[d][j];
      }
    }
  }

  //compute Jacobians
  T J, J_a; double J_A;
  J_A=std::sqrt(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
  J_a=std::sqrt(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
  J=J_a/J_A;
  if (J_A<0.0) {std::cout << "negative jacobian\n";  exit(-1);}
  //std::cout << J.val() << ", ";

  //compute normal
  T normal[3];
  normal[0]=(dxdR[1][0]*dxdR[2][1]-dxdR[2][0]*dxdR[1][1])/J_a;
  normal[1]=(dxdR[2][0]*dxdR[0][1]-dxdR[0][0]*dxdR[2][1])/J_a;
  normal[2]=(dxdR[0][0]*dxdR[1][1]-dxdR[1][0]*dxdR[0][1])/J_a;
  //std::cout << normal[0].val() << ", " << normal[1].val() << ", " << normal[2].val() << "\n";
  
  //compute curvature tensor, b
  T b[2][2];
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      b[i][j]=0.0;
      for(unsigned int d=0; d<3; d++){
	b[i][j]+=normal[d]*dxdR2[d][i][j];
      }
    }
  }

  //compute determinants of metric tensors and curvature tensor
  T det_a=a[0][0]*a[1][1]-a[0][1]*a[1][0];
  T det_b=b[0][0]*b[1][1]-b[0][1]*b[1][0];
  double det_A=A[0][0]*A[1][1]-A[0][1]*A[1][0];
  
  //compute contra-variant metric tensors, a_contra, A_contra.
  //needed for computing the contra-variant tanget vectors dxdR_contra, which are needed for computing the Christoffel symbols
  T a_contra[2][2], dxdR_contra[3][2];
  double A_contra[2][2];
  a_contra[0][0]=a[1][1]/det_a; a_contra[1][1]=a[0][0]/det_a;
  a_contra[0][1]=-a[0][1]/det_a; a_contra[1][0]=-a[1][0]/det_a;
  A_contra[0][0]=A[1][1]/det_A; A_contra[1][1]=A[0][0]/det_A;
  A_contra[0][1]=-A[0][1]/det_A; A_contra[1][0]=-A[1][0]/det_A;
  for(unsigned int d=0; d<3; d++){
    dxdR_contra[d][0]=a_contra[0][0]*dxdR[d][0]+a_contra[0][1]*dxdR[d][1];
    dxdR_contra[d][1]=a_contra[1][0]*dxdR[d][0]+a_contra[1][1]*dxdR[d][1];
  }

  //compute invariants of Green-Lagrange strain, C
  T I1=0.0;
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      I1+=A_contra[i][j]*a[i][j];
    }
  }
  
  //compute contra-variant curvature tensor, b_contra.
  T b_contra[2][2];
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      b_contra[i][j]=0.0;
      for(unsigned int k=0; k<2; k++){
	for(unsigned int l=0; l<2; l++){
	  b_contra[i][j]+=a_contra[i][k]*b[k][l]*a_contra[l][j]; //*
	}
      }
    }
  }
  
  //compute Christoffel symbols
  T Gamma[2][2][2];
  for(unsigned int i=0; i<2; i++){
    for(unsigned int j=0; j<2; j++){
      for(unsigned int k=0; k<2; k++){
	Gamma[i][j][k]=0.0;
	for(unsigned int d=0; d<3; d++){
	  Gamma[i][j][k]+=dxdR_contra[d][i]*dxdR2[d][j][k];
	}
      }
    }
  }

  //compute mean curvature, H
  T H=0;
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      H+=0.5*a[i][j]*b_contra[i][j];
    }
  }
  T dH=H-H0;
  //std::cout << dH.val() << ", ";
  
  //compute Gaussian curvature, Kappa
  T Kappa=det_b/det_a;
    
  //compute contra-variant stress and bending moment tensors
  T sigma_contra[2][2], M_contra[2][2];
  //For Helfrich energy formulation
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      sigma_contra[i][j]=(L*L/K)*((Lambda*(J-1.0)+K*dH*dH - KGaussian*Kappa)*a_contra[i][j]-2*K*dH*b_contra[i][j]);
      sigma_contra[i][j]+=(L*L/K)*(Mu/(J*J))*(A_contra[i][j]-0.5*I1*a_contra[i][j]); //stabilization term
      M_contra[i][j]=(L/K)*(K*dH + 2*KGaussian*H)*a_contra[i][j]-KGaussian*b_contra[i][j];
    }
  }
  //std::cout << sigma[0][0].val() << ", " << sigma[0][1].val() << ", " << sigma[1][0].val() << ", " << sigma[1][1].val() << "  :  ";
  //std::cout << M[0][0].val() << ", " << M[0][1].val() << ", " << M[1][0].val() << ", " << M[1][1].val() << "\n";
  
  bool surfaceFlag=p->atboundary;
  PetscReal *boundaryNormal = p->normal;
  //Residual
  if (!surfaceFlag) {
    for (unsigned int n=0; n<(unsigned int)nen; n++) {
      for (unsigned int i=0; i<3; i++){
	T Ru_i=0.0;
	for (unsigned int j=0; j<2; j++){
	  for (unsigned int k=0; k<2; k++){
	    //sigma*grad(Na)*dxdR*J
	    Ru_i += sigma_contra[j][k]*N1[n][j]*dxdR[i][k]*J;
	    //M*(hess(Na)-Gamma*grad(Na))*n*J
	    Ru_i += M_contra[j][k]*(N2[n][j][k])*normal[i]*J;
	    for (unsigned int l=0; l<2; l++){
	      Ru_i += -M_contra[j][k]*(Gamma[l][j][k]*N1[n][l])*normal[i]*J;
	    }
	  }
	}
	R[n*3+i] = Ru_i;
      }
    }
  }
  else{
    PetscReal theta=t*3.14159/6.0;
    PetscReal pCoords[3];
    IGAPointFormGeomMap(p,pCoords);
    if (pCoords[0]>0.5*(3*user->l)){theta*=-1;}
    PetscReal nValue[3]={std::sin(theta), 0.0, std::cos(theta)};
    //
    for (unsigned int n=0; n<(unsigned int)nen; n++) {
      for (unsigned int i=0; i<3; i++){
	T Ru_i=0.0;
	for (unsigned int j=0; j<2; j++){
	  for (unsigned int d=0; d<3; d++){
	    Ru_i+=epsilonBar*std::abs(nValue[0])*(N1[n][j]*normal[i]*nValue[d]*dxdR_contra[d][j]);
	  }
	}
	R[n*3+i] = Ru_i;
      }
    }
  }
  return 0; 
}


#undef  __FUNCT__
#define __FUNCT__ "Residual"
PetscErrorCode Residual(IGAPoint p,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx)
{
  Function<PetscReal>(p, shift, V, t, U, t0, U0, R, ctx);
  /*
  const PetscInt nen=p->nen, dof=3;
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
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "Jacobian"
PetscErrorCode Jacobian(IGAPoint p,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscScalar *U,
			PetscReal t0,const PetscScalar *U0,
			PetscScalar *K,void *ctx)
{
  //std::cout << "J";
  AppCtx *user = (AppCtx *)ctx;
  const PetscInt nen=p->nen, dof=3;
  //const PetscReal (*U2)[DIM] = (PetscReal (*)[DIM])U;
  /*
  if (dof*nen!=numVars) {
    PetscPrintf(PETSC_COMM_WORLD,"\ndof*nen!=numVars.... Set numVars = %u\n",dof*nen); exit(-1);
  }
  */
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
      	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] = R[n1*dof+d1].fastAccessDx(n2*dof+d2);
	}
      }
    }				
  }
  return 0;    
}


//snes convegence test
PetscErrorCode SNESConverged_Interactive(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx){
  AppCtx *user  = (AppCtx*) ctx;
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);  
  //custom test
  if ((it>199) || (fnorm<1.0e-10)){
    *reason = SNES_CONVERGED_ITS;
    return(0);
  }

  //default test
  PetscFunctionReturn(SNESConvergedDefault(snes,it,xnorm,snorm,fnorm,reason,ctx));
}

#undef __FUNCT__
#define __FUNCT__ "OutputMonitor"
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  char           filename[256];
  sprintf(filename,"./outU%d.vts",it_number);
  ierr = IGADrawVecVTK(user->iga,U,filename);CHKERRQ(ierr);
  //std::cout << c_time << "\n";
  //
  ierr = IGASetBoundaryValue(user->iga,1,0,1,-user->l*0.5*(c_time));CHKERRQ(ierr); //Y=lambda on \eta_2=0
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {

  char           filename[PETSC_MAX_PATH_LEN] = "mesh.dat";
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  AppCtx user;
  user.l=1.0;
  user.kMean=1.0;
  user.kGaussian=-0.5*user.kMean;
  user.delta=2.5;
  user.mu=0.01;
  user.epsilon=100*user.kMean/user.l;
  
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDof(iga,3);CHKERRQ(ierr); // dofs = {ux,uy,uz}
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetGeometryDim(iga,3);CHKERRQ(ierr);
  ierr = IGARead(iga,filename);CHKERRQ(ierr);  

  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  user.iga = iga;
  
  //Print knots to output
  IGAAxis axisX;  
  ierr = IGAGetAxis(iga,0,&axisX);CHKERRQ(ierr);
  PetscInt mX; PetscReal* UX;
  IGAAxisGetKnots(axisX, &mX, &UX);
  std::cout << mX << " knotsX: ";
  for (unsigned int i=0; i<(mX+1); i++){
    std::cout << UX[i] << ", ";
  }
  std::cout << std::endl;
  IGAAxis axisY;  
  ierr = IGAGetAxis(iga,1,&axisY);CHKERRQ(ierr);
  PetscInt mY; PetscReal* UY;
  IGAAxisGetKnots(axisY, &mY, &UY);
  std::cout << mY << " knotsY: ";
  for (unsigned int i=0; i<(mY+1); i++){
    std::cout << UY[i] << ", ";
  }
  std::cout << std::endl;
  
  //Dirichlet BC's u = 0, v = [0:1]
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr); //X=0 on \eta_1=0
  ierr = IGASetBoundaryValue(iga,0,0,2,0.0);CHKERRQ(ierr); //Z=0 on \eta_1=0
  ierr = IGASetBoundaryValue(iga,0,1,2,0.0);CHKERRQ(ierr); //Z=0 on \eta_1=1
  ierr = IGASetBoundaryValue(iga,1,1,1,0.0);CHKERRQ(ierr); //Y=0 on \eta_2=1
  
  //Boundary form for Neumann BC's
  
  IGAForm form;
  ierr = IGAGetForm(iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_FALSE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_FALSE);CHKERRQ(ierr);
  
  // // //
  Vec U,U0;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&U0);CHKERRQ(ierr);
  //ierr = FormInitialCondition(iga, U0, &user); //set initial conditions
  ierr = VecSet(U0,0.0);CHKERRQ(ierr);
  ierr = VecCopy(U0, U);CHKERRQ(ierr);
  //
  ierr = IGASetFormIEFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(iga,Jacobian,&user);CHKERRQ(ierr);
  //
  ierr = IGADrawVecVTK(iga,U,"mesh.vts");CHKERRQ(ierr);
  //
  TS ts;
  PetscInt timeSteps=100;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  //ierr = TSSetMaxSteps(ts,timeSteps+1);CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts, 1.0);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,1.0/timeSteps);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  
  //
  SNES snes;
  TSGetSNES(ts,&snes);
  SNESSetConvergenceTest(snes,SNESConverged_Interactive,(void*)&user,NULL);
  //SNESLineSearch ls;
  //SNESGetLineSearch(snes,&ls);
  //SNESLineSearchSetType(ls,SNESLINESEARCHBT);
  //ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  //SNESLineSearchView(ls,NULL);
 
#if PETSC_VERSION_LE(3,3,0)
  ierr = TSSolve(ts,U,NULL);CHKERRQ(ierr);
#else
  ierr = TSSolve(ts,U);CHKERRQ(ierr);
#endif

  //
  PetscMPIInt rank,size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if(rank == size-1){
    //PetscScalar value;
    //PetscInt index;
    //ierr = VecGetValues(U0,1,&index,&value);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF,"x[%d]=%g\n",index,(double)value);CHKERRQ(ierr);
  }

  ierr = IGAWriteVec(iga,U,"mesh.out");CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  ierr = IGADrawVecVTK(iga,U,"mesh.vts");CHKERRQ(ierr);
#endif

  //PetscBool draw = IGAGetOptBool(NULL,"-draw",PETSC_FALSE);
  //if (draw) {ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);}
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&U0);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
