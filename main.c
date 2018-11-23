/*
  Kirchhoff-Love shell implementation in PetIGA
  author: Shiva Rudraraju
 */

#include "petiga.h"
#include <math.h> 
#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> doubleAD;

#define LagrangeMultiplierMethod
#define enableFastResidualComputation
#define enableForceControl //default is displacement control for some BVPs

#include "include/residual.h"
#include "include/project.h"
#include "include/output.h"
#include "include/solvers.h"

//parameters
#define bvpType 3
#define stabilizationMethod 8 //Note: Method 8 will make the solution a bit time step dependent as previous time step solution (dx0dR, aPre terms, etc) are used.
#define numLoadSteps 1000

#undef  __FUNCT__
#define __FUNCT__ "setBCs"
PetscErrorCode setBCs(BVPStruct& bvp, PetscInt it_number, PetscReal c_time)
{
  PetscFunctionBegin; 
  PetscErrorCode ierr;

  ierr = IGASetFixTable(bvp.iga,NULL);CHKERRQ(ierr); /* Clear vector to read BCs from */

  //clear old BC's
  IGAForm form;
  ierr = IGAGetForm(bvp.iga,&form);CHKERRQ(ierr);
  for (PetscInt dir=0; dir<2; dir++){
    for (PetscInt side=0; side<2; side++){
      ierr =   IGAFormClearBoundary(form,dir,side);
    }
  }
  
  //Boundary form for Neumann BC's
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_FALSE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_FALSE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_FALSE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_FALSE);CHKERRQ(ierr);
  bvp.angleConstraints[0]= bvp.angleConstraints[1]=false;
  
  //Dirichlet and Neumann BC's
  switch (bvp.type) {
  case 0: //cap BVP
    bvp.uDirichlet=0.8*c_time*bvp.l*1.0; //X=Z=uDirichlet at the bottom of the cap (displacement control)
    ProjectL2(&bvp);
  
    //Dirichlet
    ierr = IGASetBoundaryValue(bvp.iga,0,1,0,0.0);CHKERRQ(ierr); //X=0 at the top of the cap
    ierr = IGASetBoundaryValue(bvp.iga,0,1,2,0.0);CHKERRQ(ierr); //Z=0 at the top of the cap
    ierr = IGASetBoundaryValue(bvp.iga,0,0,1,0.0);CHKERRQ(ierr); //Y=0 at the bottom of the cap
    ierr = IGASetBoundaryValue(bvp.iga,0,0,0,/*dummy*/0.0);CHKERRQ(ierr); //init for X=uDirichlet at the bottom of the cap
    ierr = IGASetBoundaryValue(bvp.iga,0,0,2,/*dummy*/0.0);CHKERRQ(ierr); //init for Z=uDirichlet at the bottom of the cap
   
    //Neumann
    ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr); //phi=90 at the bottom of the cap
    ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_FALSE);CHKERRQ(ierr); //phi=0  at the top of the cap
    bvp.angleConstraints[0]=true; bvp.angleConstraintValues[0]=90;
    bvp.angleConstraints[1]=false; bvp.angleConstraintValues[1]=0;
    bvp.epsilon=bvp.kMean;
    
    //Non-homogeneous Dirichlet BC values
    ierr = IGASetFixTable(bvp.iga,bvp.xDirichlet);CHKERRQ(ierr);    /* Set vector to read BCs from */
    break;
    
  case 1: //tube BVP
#ifdef enableForceControl
    bvp.isCollar=true;
    bvp.CollarLocation=bvp.l*0.0; //At the bottom
    bvp.CollarHeight=bvp.l*0.025; //0.5nm (20*0.025)
    bvp.CollarPressure=c_time*55;
#else
    bvp.uDirichlet=0.9*c_time*bvp.l*1.0; //X=Z=uDirichlet at the bottom of the tube (displacement control)
    ProjectL2(&bvp);
#endif
    
    //Dirichlet
    ierr = IGASetBoundaryValue(bvp.iga,0,1,0,0.0);CHKERRQ(ierr); //X=0 at the top of the tube
    ierr = IGASetBoundaryValue(bvp.iga,0,1,2,0.0);CHKERRQ(ierr); //Z=0 at the top of the tube
    ierr = IGASetBoundaryValue(bvp.iga,0,1,1,0.0);CHKERRQ(ierr); //Y=0 at the top of the tube
    //ierr = IGASetBoundaryValue(bvp.iga,0,0,1,0.0);CHKERRQ(ierr); //Y=0 at the bottom of the tube
#ifndef  enableForceControl
    ierr = IGASetBoundaryValue(bvp.iga,0,0,0,/*dummy*/0.0);CHKERRQ(ierr); //init for X=uDirichlet at the bottom of the tube
    ierr = IGASetBoundaryValue(bvp.iga,0,0,2,/*dummy*/0.0);CHKERRQ(ierr); //init for Z=uDirichlet at the bottom of the tube
#endif
    //Neumann. Comment out for Asymmetric mode.
    ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr); //phi=90 at the bottom of the tube
    ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr); //phi=90 at the top of the tube
    bvp.angleConstraints[0]=true; bvp.angleConstraintValues[0]=90;
    bvp.angleConstraints[1]=true; bvp.angleConstraintValues[1]=90;
    bvp.epsilon=bvp.kMean;
    
    
#ifndef  enableForceControl
    //Non-homogeneous Dirichlet BC values
    ierr = IGASetFixTable(bvp.iga,bvp.xDirichlet);CHKERRQ(ierr);    /* Set vector to read BCs from */
#endif
    break;
    
  case 2: //base BVP
#ifdef enableForceControl
    bvp.isCollar=true;
    bvp.CollarLocation=bvp.l*0.45;
    bvp.CollarHeight=bvp.l*0.1; //2.0nm 
    bvp.CollarPressure=c_time*3.7;
#else
    bvp.uDirichlet=c_time*bvp.l*1.0; //X=Z=uDirichlet at the bottom of the base (displacement control)
    ProjectL2(&bvp);
#endif
    
    //Dirichlet
    ierr = IGASetBoundaryValue(bvp.iga,0,1,0,0.0);CHKERRQ(ierr); //X=0 at the top of the base
    ierr = IGASetBoundaryValue(bvp.iga,0,1,2,0.0);CHKERRQ(ierr); //Z=0 at the top of the base
    //ierr = IGASetBoundaryValue(bvp.iga,0,1,1,0.0);CHKERRQ(ierr); //Y=0 at the top of the base
    ierr = IGASetBoundaryValue(bvp.iga,0,0,1,0.0);CHKERRQ(ierr); //Y=0 at the bottom of the base
#ifndef  enableForceControl
    ierr = IGASetBoundaryValue(bvp.iga,0,0,0,/*dummy*/0.0);CHKERRQ(ierr); //init for X=uDirichlet at the bottom of the base
    ierr = IGASetBoundaryValue(bvp.iga,0,0,2,/*dummy*/0.0);CHKERRQ(ierr); //init for Z=uDirichlet at the bottom of the base
#endif

    //Neumann
    bvp.angleConstraints[0]=false; 
    bvp.angleConstraints[1]=false; 
    bvp.epsilon=bvp.kMean*0.0;

#ifndef  enableForceControl
    //Non-homogeneous Dirichlet BC values
    ierr = IGASetFixTable(bvp.iga,bvp.xDirichlet);CHKERRQ(ierr);    /* Set vector to read BCs from */
#endif
    break;
  case 3: //pullout BVP
    //properties
    bvp.kGaussian=-0.7*bvp.kMean; //Gaussian curvature modulus
    //bottom surface
    ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
    bvp.surfaceTensionAtBase=1.0;
    //topsurface
    ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
    bvp.tractionOnTop=c_time*300; //600;
#ifndef  enableForceControl
    bvp.uDirichlet= (c_time)*bvp.l*0.05; //pull out height
    ierr = IGASetBoundaryValue(bvp.iga,0,1,1,bvp.uDirichlet);CHKERRQ(ierr); //Y at the top of the base
#endif

    //Dirichlet
    ierr = IGASetBoundaryValue(bvp.iga,0,1,0,0.0);CHKERRQ(ierr); //X=0 at the top of the base
    ierr = IGASetBoundaryValue(bvp.iga,0,1,2,0.0);CHKERRQ(ierr); //Z=0 at the top of the base
    ierr = IGASetBoundaryValue(bvp.iga,0,0,1,0.0);CHKERRQ(ierr); //Y=0 at the bottom of the base
#ifndef  enableForceControl
    //ierr = IGASetBoundaryValue(bvp.iga,0,0,0,0.0);CHKERRQ(ierr); 
    //ierr = IGASetBoundaryValue(bvp.iga,0,0,2,0.0);CHKERRQ(ierr); 
#endif

    //Neumann. Comment out for Asymmetric mode.
    ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr); //phi=0 at the bottom of the tube
    bvp.angleConstraints[0]=true; bvp.angleConstraintValues[0]=0;
    bvp.epsilon=10*bvp.kMean;

    break;
  }
  //
  PetscFunctionReturn(0);
}

//main
int main(int argc, char *argv[]) {
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  //
  PetscMPIInt rank,size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  //IGA and BVP setup
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  
  //setup BVP and material model parameters
  BVPStruct bvp;
  //non-dimentionalzation factors used for the pinching problems
  bvp.lengthFactor=1.0;
  bvp.kFactor=1.0; 
  bvp.forceFactor=1.0; 
  bvp.energyFactor=1.0;
  //material constants (in actual units)
  bvp.l=20.0;              //20nm
  bvp.kMean=320.0*16.0;         //320pN-nm, mean curvature modulus
  bvp.kGaussian=0.0;       //Gaussian curvature modulus
  bvp.mu=1.0*bvp.kMean;       //shear modulus for stabilization terms
  bvp.lambda=10*bvp.kMean;        //penalty parameter
  bvp.surfaceTensionAtBase=0.0;
  bvp.tractionOnTop=0.0;
  bvp.epsilon=0.0;       //penalty parameter for rotational constraints
  bvp.xMin=bvp.l;
  //
  bvp.type=bvpType;
  bvp.stabilization=stabilizationMethod;
  bvp.c_time=0.0;
  bvp.isCollar=false;
  bvp.isCollarHelix=false;
  bvp.CollarLocation=0.0;
  bvp.CollarHeight=0.0;
  bvp.CollarHelixPitch=0.0;
  bvp.CollarRadius=0.0;
  bvp.CollarPressure=0.0;
  
  //processor zero for printing output in MPI jobs
  if(rank == 0){bvp.isProc0=true;}
  else {bvp.isProc0=false;}
  
#ifdef LagrangeMultiplierMethod
  ierr = IGASetDof(iga,4);CHKERRQ(ierr); // dofs = {ux,uy,uz,q}
#else
  ierr = IGASetDof(iga,3);CHKERRQ(ierr); // dofs = {ux,uy,uz}
#endif
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetGeometryDim(iga,3);CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(iga->axis[1],PETSC_TRUE);CHKERRQ(ierr);
  switch (bvp.type) {
  case 0: //cap BVP
    ierr = IGARead(iga,"meshes/capMeshr80h40C1.dat"); CHKERRQ(ierr);
    break;
  case 1: //tube BVP
    //ierr = IGARead(iga,"meshes/tubeMeshr160h80C1.dat"); CHKERRQ(ierr);
    ierr = IGARead(iga,"meshes/tubeMeshr80h80C1H2R.dat"); CHKERRQ(ierr);
    break;
  case 2: //base BVP
    ierr = IGARead(iga,"meshes/base90DegMeshr80h40C1H2R.dat"); CHKERRQ(ierr);
    break;
  case 3: //pullout BVP
    ierr = IGARead(iga,"meshes/baseCircleMeshr40h80C1.dat"); CHKERRQ(ierr);
    //ierr = IGARead(iga,"meshes/baseCircleMeshr60h40C1.dat"); CHKERRQ(ierr);
    break;
  }
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  bvp.iga = iga;
  
  //Print knots to output
  IGAAxis axisX;  
  ierr = IGAGetAxis(bvp.iga,0,&axisX);CHKERRQ(ierr);
  PetscInt mX; PetscReal* UX;
  IGAAxisGetKnots(axisX, &mX, &UX);
  IGAAxis axisY;  
  ierr = IGAGetAxis(bvp.iga,1,&axisY);CHKERRQ(ierr);
  PetscInt mY; PetscReal* UY;
  IGAAxisGetKnots(axisY, &mY, &UY);
  if(bvp.isProc0){
    std::cout << mX << " knotsX: ";
    for (unsigned int i=0; i<(mX+1); i++){
      std::cout << UX[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << mY << " knotsY: ";
    for (unsigned int i=0; i<(mY+1); i++){
      std::cout << UY[i] << ", ";
    }
    std::cout << std::endl;
  }
  
  //solutions vectors
  Vec U;
  ierr = IGACreateVec(bvp.iga,&U);CHKERRQ(ierr);
  ierr = VecSet(U,0.0);CHKERRQ(ierr);
  ierr = IGADrawVecVTK(iga,U,"mesh.vts");CHKERRQ(ierr);
  ierr = IGACreateVec(bvp.iga,&bvp.xDirichlet);CHKERRQ(ierr);
  //
  ierr = IGASetFormIEFunction(bvp.iga,Residual,&bvp);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(bvp.iga,Jacobian,&bvp);CHKERRQ(ierr);
  //
  setBCs(bvp, 0, 0.0);
  
  //open file for U,R output
  if (bvp.isProc0){
#ifdef enableForceControl
    bvp.fileForUROutout=fopen ("URbyForceControl.txt","w");
    if (bvp.type!=3){ 
      fprintf (bvp.fileForUROutout, "%12s, %12s, %12s, %12s\n", "Radius[nm]", "Pressure[pN/nm]", "E1[pN-nm]", "E2[pN-nm]");
    }
    else{
      fprintf (bvp.fileForUROutout, "%12s, %12s, %12s, %12s\n", "Height[nm]", "Force[pN]", "E1[pN-nm]", "E2[pN-nm]");
    }
#else
    bvp.fileForUROutout=fopen ("URbyDisplacementControl.txt","w");
    fprintf (bvp.fileForUROutout, "%12s, %12s, %12s, %12s\n", "Radius[nm]", "Reaction[pN]", "E1[pN-nm]", "E2[pN-nm]");
#endif
  }
  
  //load stepping
  TS ts;
  PetscInt timeSteps=numLoadSteps;
  ierr = IGACreateTS(bvp.iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  //ierr = TSSetMaxSteps(ts,timeSteps+1);CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts, 1.01);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,1.0/timeSteps);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor,&bvp,NULL);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  SNES snes;
  TSGetSNES(ts,&snes);
  SNESSetConvergenceTest(snes,SNESConverged_Interactive,(void*)&bvp,NULL);
#if PETSC_VERSION_LE(3,3,0)
  ierr = TSSolve(ts,U,NULL);CHKERRQ(ierr);
#else
  ierr = TSSolve(ts,U);CHKERRQ(ierr);
#endif

  //
  if (bvp.isProc0){fclose(bvp.fileForUROutout);}
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&bvp.iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
