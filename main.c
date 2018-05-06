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

#define numVars 27
//include automatic differentiation library
#include <Sacado.hpp>
typedef Sacado::Fad::SFad<double,numVars> doubleAD;
//typedef Sacado::Fad::DFad<double> doubleAD;

typedef struct {
  PetscReal nu,E,t,k;
} AppCtx;

#undef  __FUNCT__
#define __FUNCT__ "Function"
template <class T>
PetscErrorCode Function(IGAPoint p,
				PetscReal shift,const PetscScalar *V,
				PetscReal t,const T * U,
				PetscReal t0,const PetscScalar * U0,
				T *R,void *ctx)
{
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
  //std::cout << "s";
  Function(p, shift, V, t, U, t0, U0, R, ctx);
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
  return 0;
  AppCtx *user = (AppCtx *)ctx;
  const PetscInt nen=p->nen, dof=DIM+1;
  const PetscReal (*U2)[DIM+1] = (PetscReal (*)[DIM+1])U;
  if (dof*nen!=numVars) {
    PetscPrintf(PETSC_COMM_WORLD,"\ndof*nen!=numVars.... Set numVars = %u\n",dof*nen); exit(-1);
  }
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
  ierr = IGASetBoundaryValue(iga,0,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,0,2,0.0);CHKERRQ(ierr);
  // Boundary conditions on u = 1, v = [0:1]
  ierr = IGASetBoundaryValue(iga,0,1,0,0.0);CHKERRQ(ierr);
  // Boundary conditions on u = [0:1], v = 0
  ierr = IGASetBoundaryValue(iga,1,0,1,0.0);CHKERRQ(ierr);
  // Boundary conditions on u = [0:1], v = 1
  ierr = IGASetBoundaryValue(iga,1,1,1,0.0);CHKERRQ(ierr);

  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);

  // // //
  Vec U,U0;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&U0);CHKERRQ(ierr);
  //ierr = FormInitialCondition(iga, U0, &user); //set initial conditions
  ierr = VecSet(U0,1.0);CHKERRQ(ierr);
  ierr = VecCopy(U0, U);CHKERRQ(ierr);
  //
  ierr = IGASetFormIEFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(iga,Jacobian,&user);CHKERRQ(ierr);
  //
  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetMaxSteps(ts,10);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,0.1);CHKERRQ(ierr);
  //ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);
  //
  SNES snes;
  TSGetSNES(ts,&snes);
  //SNESSetConvergenceTest(snes,SNESConverged_Interactive,(void*)&user,NULL);
  //SNESLineSearch ls;
  //SNESGetLineSearch(snes,&ls);
  //SNESLineSearchSetType(ls,SNESLINESEARCHBT);
  //ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  //SNESLineSearchView(ls,NULL);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
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

  ierr = IGAWriteVec(iga,U0,"mesh.out");CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  ierr = IGADrawVecVTK(iga,U0,"mesh.vts");CHKERRQ(ierr);
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
