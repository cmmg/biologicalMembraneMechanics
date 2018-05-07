/*
  Kirchhoff-Love shell implementation in PetIGA
  author: rudraa
 */
#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
#include "fields.h"

#define DIM 3
#define DOF 3

#define numVars 27
//include automatic differentiation library
#include <Sacado.hpp>
//typedef Sacado::Fad::SFad<double,numVars> doubleAD;
typedef Sacado::Fad::DFad<double> doubleAD;

typedef struct {
  PetscReal nu,E,t,k;
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
  PetscReal nu = user->nu;
  PetscReal E  = user->E;
  PetscReal mu=user->E;
  PetscReal thick  = user->t;
  PetscReal k  = user->k;

  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //shape functions: value, grad, hess
  //PetscReal (*N) = (PetscReal (*)) p->shape[0];
  PetscReal (*N1)[2] = (PetscReal (*)[2]) p->basis[1];
  //PetscReal (*N2)[2][2] = (PetscReal (*)[2][2]) p->shape[2];
  
  //get X
  const PetscReal *tempX = p->geometry;
  PetscReal (*X)[3] = (PetscReal (*)[3])tempX;
  //get x
  T x[nen][3];
  T (*u)[3] = (T (*)[3])tempU;
  for(unsigned int n=0; n<(unsigned int) nen; n++){
    for(unsigned int d=0; d<3; d++) x[n][d]= X[n][d]+ u[n][d];
  }

  //compute metric tensors, a and A
  T a[2][2], basis_x[3][2];
  double A[2][2], basis_X[3][2];;
  for(unsigned int d=0; d<3; d++){
    basis_X[d][0]=basis_X[d][1]=0.0;
    basis_x[d][0]=basis_x[d][1]=0.0;
    for(unsigned int n=0; n<(unsigned int) nen; n++){
      basis_X[d][0]+=N1[n][0]*X[n][d];
      basis_X[d][1]+=N1[n][1]*X[n][d];	 
      basis_x[d][0]+=N1[n][0]*x[n][d];
      basis_x[d][1]+=N1[n][1]*x[n][d];	 
    }
  }

  for(unsigned int d1=0; d1<2; d1++){
    for(unsigned int d2=0; d2<2; d2++){
      A[d1][d2]=0.0;
      a[d1][d2]=0.0;
      for(unsigned int n=0; n<3; n++){
	A[d1][d2]+=basis_X[n][d1]*basis_X[n][d2];
	a[d1][d2]+=basis_x[n][d1]*basis_x[n][d2];
      }
    }
  }

  //compute Jacobian
  T J, J_a; double J_A;
  J_A=(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
  J_a=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
  //std::cout << J_A << ", " << A[0][0] <<"\n";
  //if (J_A<0.0) exit(-1);
  J=J_a/J_A;
  
  //Stress
  T Tau[2][2];
  //For incompressible Neo-Hookean material
  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++){
      Tau[i][j]=mu*(A[i][j]-a[i][j]/(J*J));
    }
  }

  //Residual
  for (unsigned int n=0; n<(unsigned int)nen; n++) {
    for (unsigned int i=0; i<3; i++){
      T Ru_i=0.0;
      for (unsigned int j=0; j<2; j++){
	for (unsigned int k=0; k<2; k++){
	  //grad(Na)*Tau
	  Ru_i += N1[n][j]*Tau[j][k]*basis_x[i][k];
	}
      }
      R[n*3+i] = Ru_i;
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
  //std::cout << "R";
  Function<PetscReal>(p, shift, V, t, U, t0, U0, R, ctx);
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
  const PetscInt nen=p->nen, dof=DIM;
  const PetscReal (*U2)[DIM] = (PetscReal (*)[DIM])U;
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

int main(int argc, char *argv[]) {

  char           filename[PETSC_MAX_PATH_LEN] = "mesh.dat";
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  AppCtx user;
  user.nu = 0.3;
  user.E  = 3.e7;
  user.t  = 1.;
  user.k  = 5./6.;

  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDof(iga,3);CHKERRQ(ierr); // dofs = {ux,uy,uz}
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetGeometryDim(iga,3);CHKERRQ(ierr);
  ierr = IGARead(iga,filename);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);

  // Boundary conditions on u = 0, v = [0:1]
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,0,2,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,1,2,0.0001);CHKERRQ(ierr);
  //
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);

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
  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetMaxSteps(ts,2);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,0.1);CHKERRQ(ierr);
  //ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  
  //
  //SNES snes;
  //TSGetSNES(ts,&snes);
  //SNESSetConvergenceTest(snes,SNESConverged_Interactive,(void*)&user,NULL);
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
