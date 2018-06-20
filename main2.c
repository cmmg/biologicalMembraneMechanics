/*
  Kirchhoff-Love shell implementation in PetIGA
  author: Shiva Rudraraju
 */

#include "petiga.h"
#include <math.h> 

#define LagrangeMultiplierMethod

//include automatic differentiation library
#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> doubleAD;

#include "include/residual.h"
#include "include/output.h"

//main
int main(int argc, char *argv[]) {

  char           filename[PETSC_MAX_PATH_LEN] = "mesh.dat";
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  AppCtx user;
  user.l=1.0;
  user.kMean=1.0;
  user.kGaussian=0*-0.5*user.kMean;
  user.mu=0.01;
  user.epsilon=0*user.kMean/user.l;
#ifndef LagrangeMultiplierMethod
  user.delta=1000.0;
#endif
  user.c_time=0.0;
  
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
#ifdef LagrangeMultiplierMethod
  ierr = IGASetDof(iga,4);CHKERRQ(ierr); // dofs = {ux,uy,uz, lambda}
#else
  ierr = IGASetDof(iga,3);CHKERRQ(ierr); // dofs = {ux,uy,uz}
#endif
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetGeometryDim(iga,3);CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(iga->axis[1],PETSC_TRUE);CHKERRQ(ierr);
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
  
  //Boundary form for Neumann BC's
  IGAForm form;
  ierr = IGAGetForm(iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_FALSE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_FALSE);CHKERRQ(ierr);
  
  //
  Vec U,U0;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&U0);CHKERRQ(ierr);
  ierr = VecSet(U0,0.0);CHKERRQ(ierr);
  ierr = VecCopy(U0, U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&bvp.x);CHKERRQ(ierr);
  
  //
  ierr = IGASetFormIEFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(iga,Jacobian,&user);CHKERRQ(ierr);

  //
  TS ts;
  PetscInt timeSteps=1000;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  //ierr = TSSetMaxSteps(ts,timeSteps+1);CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts, 1.01);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,1.0/timeSteps);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  SNES snes;
  TSGetSNES(ts,&snes);
  SNESSetConvergenceTest(snes,SNESConverged_Interactive,(void*)&user,NULL);
  ierr = TSSolve(ts,U);CHKERRQ(ierr);

  //
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&U0);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
