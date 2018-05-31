/*
  Kirchhoff-Love shell implementation in PetIGA
  author: rudraa
 */
#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
#include "fields.h"

//include automatic differentiation library
#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> doubleAD;

typedef struct {
  IGA iga;
  PetscReal k, delta;
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
  PetscReal K  = user->k; //bending modulus
  PetscReal delta  = user->delta; //penalty parameter for incompressibility
  
  //instantaneous curvature
  double H0=5.0;
  
  //get number of shape functions (nen) and dof's
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //shape functions: value, grad, hess
  const PetscReal (*N) = (const PetscReal (*)) p->shape[0];
  const PetscReal (*N1)[2] = (const PetscReal (*)[2]) p->basis[1];
  const PetscReal (*N2)[2][2] = (const PetscReal (*)[2][2]) p->basis[2];
  
  //get X
  const PetscReal (*X)[3] = (const PetscReal (*)[3]) p->geometry;;

  //get x
  T x[nen][3];
  T (*u)[3] = (T (*)[3])tempU;
  for(unsigned int n=0; n<(unsigned int) nen; n++){
    for(unsigned int d=0; d<3; d++){
      x[n][d]= X[n][d]+ u[n][d];
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
  //if (J_A<0.0) exit(-1);
  J=J_a/J_A;
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

  //compute contra-variant metric tensor, a_contra.
  //needed for computing the contra-variant tanget vectors dxdR_contra, which are needed for computing the Christoffel symbols
  T a_contra[2][2], dxdR_contra[3][2];
  T det_a=a[0][0]*a[1][1]-a[0][1]*a[1][0];
  a_contra[0][0]=a[1][1]/det_a; a_contra[1][1]=a[0][0]/det_a;
  a_contra[0][1]=-a[0][1]/det_a; a_contra[1][0]=-a[1][0]/det_a;
  for(unsigned int d=0; d<3; d++){
    dxdR_contra[d][0]=a_contra[0][0]*dxdR[d][0]+a_contra[0][1]*dxdR[d][1];
    dxdR_contra[d][1]=a_contra[1][0]*dxdR[d][0]+a_contra[1][1]*dxdR[d][1];
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

  //compute Gaussian curvature, H
  T H=0;
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      H+=0.5*a_contra[i][j]*b[i][j];
    }
  }
  T dH=H-H0;
  //std::cout << H.val() << ", ";
  
  //compute contra-variant stress and bending moment tensors
  T sigma[2][2], M[2][2];
  //For Helfrich energy formulation
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      sigma[i][j]=((delta*(J-1.0)+K*dH*dH)*a_contra[i][j]-2*K*dH*b_contra[i][j]);
      M[i][j]=(K*dH)*a[i][j];
    }
  }
  
  //Residual
  for (unsigned int n=0; n<(unsigned int)nen; n++) {
    for (unsigned int i=0; i<3; i++){
      T Ru_i=0.0;
      for (unsigned int j=0; j<2; j++){
	for (unsigned int k=0; k<2; k++){
	  //sigma*grad(Na)*dxdR*J
	  Ru_i += sigma[j][k]*N1[n][j]*dxdR[i][k]*J;
	  //M*(hess(Na)-Gamma*grad(Na))*n*J
	  Ru_i += M[j][k]*(N2[n][j][k])*normal[i]*J;
	  for (unsigned int l=0; l<2; l++){
	    Ru_i += -M[j][k]*(Gamma[l][j][k]*N1[n][l])*normal[i]*J;
	  }
	}
      }
      R[n*3+i] = Ru_i;
    }
  }
  //R[0]+=1.0;
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
  //std::cout << "R";
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
  
  //Function<PetscReal>(p, shift, V, t, U, t0, U0, R, ctx);
  //PetscInt nen, dof;
  //IGAPointGetSizes(p,0,&nen,&dof);
  //for(int n1=0; n1<nen*dof; n1++) std::cout << U[n1] << ", ";
  //std::cout << "\n";
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
  //std::cout << "cccccccccccccccccccccccccccccccccccccccccccccc";
  return 0;    
}


/*
//Initial conditions for U
typedef struct {
  PetscReal ux, uy, uz;
} Field;

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::srand(5);
  DM da;
  ierr = IGACreateNodeDM(iga,DIM,&da);CHKERRQ(ierr);
  Field ***u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscInt i,j,k;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(k=info.zs;k<info.zs+info.zm;k++){
	u[k][j][i].ux=0.0;
	u[k][j][i].uy=0.0;
	u[k][j][i].uz=0.0;
      }
    }

  }
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}
*/

//snes convegence test
PetscErrorCode SNESConverged_Interactive(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx){
  AppCtx *user  = (AppCtx*) ctx;
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);  
  //custom test
  if ((it>500) || (fnorm<1.0e-9)){
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
  //
  //ierr = IGASetBoundaryValue(user->iga,1,1,2,0.01*(it_number+1));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {

  char           filename[PETSC_MAX_PATH_LEN] = "mesh.dat";
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  AppCtx user;
  user.k  = 1.0;
  user.delta  = 1.0;
  
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDof(iga,3);CHKERRQ(ierr); // dofs = {ux,uy,uz}
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetGeometryDim(iga,3);CHKERRQ(ierr);
  ierr = IGARead(iga,filename);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);

  // Boundary conditions on u = 0, v = [0:1]
  //Symmetric BCs
  ierr = IGASetBoundaryValue(iga,1,0,0,0.0);CHKERRQ(ierr); //X=0 on \eta_1=0
  ierr = IGASetBoundaryValue(iga,1,0,1,0.0);CHKERRQ(ierr); //Z=0 on \eta_2=0
  ierr = IGASetBoundaryValue(iga,1,0,2,0.0);CHKERRQ(ierr); //Z=0 on \eta_2=0  
  ierr = IGASetBoundaryValue(iga,1,1,2,0.01);CHKERRQ(ierr); //Y=0 on \eta_2=1
  //ierr = IGASetBoundaryValue(iga,0,1,0,0.1);CHKERRQ(ierr);

  //
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  user.iga = iga;
  
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
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetMaxSteps(ts,10);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,0.1);CHKERRQ(ierr);
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
