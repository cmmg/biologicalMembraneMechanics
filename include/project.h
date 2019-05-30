/*
  L2 Projection functions
  author: Shiva Rudraraju
 */

#if !defined(PROJECT_H_)
#define PROJECT_H_

#undef __FUNCT__
#define __FUNCT__ "FunctionEnergy"
PetscErrorCode FunctionEnergy(IGAPoint p,const PetscScalar *U,PetscInt n,PetscScalar *energy,void *ctx){
  BVPStruct *bvp = (BVPStruct *)ctx;
  PetscFunctionBegin;

  //get kinematic quantities
  KinematicsStruct<PetscScalar> k;
  getKinematics<PetscScalar>(p, U, U, k);

  //energy
  HelfrichModel<PetscScalar> m;
  m.bvp=bvp;
  computeEnergy(k,m,energy);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FunctionFields"
PetscErrorCode FunctionFields(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  BVPStruct *bvp = (BVPStruct *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //get kinematic quantities
  KinematicsStruct<PetscScalar> k;
  getKinematics<PetscScalar>(p, U, U, k);

  //get u at point
#ifdef LagrangeMultiplierMethod
  double (*tempU)[3+1] = (double (*)[3+1])U;
#else
  double (*tempU)[3] = (double (*)[3])U;
#endif
  
  double u[3];
  const PetscReal (*N) = (const PetscReal (*)) p->shape[0];;  
  for(unsigned int d=0; d<3; d++){
    u[d]=0.0;
    for(unsigned int n=0; n<(unsigned int) nen; n++){
      u[d]+=N[n]*tempU[n][d];
    }
  }
  
  //L2 projection residual
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      PetscReal val=0.0;
      switch (d1) {
      case 0:
	val=u[0]; break;
      case 1:
	val=u[1]; break;
      case 2:
	val=u[2]; break;
      case 3:
	val= k.H; break;
      }
      R[n1*dof+d1] = N[n1]*val*k.J_a;
    }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FunctionDirichletL2"
PetscErrorCode FunctionDirichletL2(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  BVPStruct *bvp = (BVPStruct *)ctx;
  
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal x[3];
  IGAPointFormGeomMap(p,x);

  //get kinematic quantities
  KinematicsStruct<PetscScalar> k;
  getKinematics<PetscScalar>(p, U, U, k);
  
  //normal direction in XZ plane
  PetscReal n[3];
  if ((x[0]*x[0]+x[2]*x[2])>0.0){
    n[0]=x[0]/sqrt(x[0]*x[0]+x[2]*x[2]);
    n[1]=0.0;
    n[2]=x[2]/sqrt(x[0]*x[0]+x[2]*x[2]);
  }
  else{
    n[0]=n[1]=n[2]=0.0;
  }
	
  PetscReal uDirichletVal=-bvp->uDirichlet;  
  if (x[1]>(0.5*bvp->l)){ uDirichletVal=0.0;}//for top surface

  //L2 projection residual
  const PetscReal (*N) = (const PetscReal (*)) p->shape[0];;  
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      PetscReal val=0.0;
      switch (d1) {
      case 0:
	val=uDirichletVal*n[0]; break;
      case 1:
	val=0.0; break;
      case 2:
	val=uDirichletVal*n[2]; break;
      case 3:
	val=0.0; break;
      }
      R[n1*dof+d1] = N[n1]*val*k.J_a;
    }
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FunctionReactions"
PetscErrorCode FunctionReactions(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx)
{
  PetscFunctionBegin;
  BVPStruct *bvp = (BVPStruct *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal pCoords[3]; IGAPointFormGeomMap(p,pCoords);

   //get kinematic quantities
  KinematicsStruct<PetscScalar> k;
  getKinematics<PetscScalar>(p, U, U, k);
  const PetscReal (*N) = (const PetscReal (*)) p->shape[0];
  
  if (bvp->type==4){ //helix BVP
    //body forces (force collar)
    bool isCollar=false;
    double CollarPressure=1.0; 
    if (bvp->isCollarHelix){
      unsigned int numTracePoints=100;
      for (unsigned int t=0; t<numTracePoints; t++){
	PetscReal theta=(((double)t)/numTracePoints)*2*PI;
	PetscReal z=bvp->CollarLocation+theta*bvp->CollarHelixPitch/(2*PI);
	PetscReal ri=bvp->CollarRadius;
	PetscReal ro=bvp->CollarHelixBaseRadius+bvp->CollarRadius;
	PetscReal r=bvp->CollarHelixBaseRadius;
	if (z<r){
	  ri=ro-(r*std::sqrt(1.0-std::pow(1.0-z/r,2.0)));
	}
	PetscReal x=ri*std::cos(theta);
	PetscReal y=ri*std::sin(theta);
	if (std::sqrt(std::pow(pCoords[0]-x,2)+std::pow(pCoords[1]-z,2)+std::pow(pCoords[2]-y,2))<=0.5*bvp->CollarHeight) {
	  isCollar=true; break;
	}
      }
    }
    //
    if (isCollar) {
      for (unsigned int n=0; n<(unsigned int)(nen*dof); n++){
	for (unsigned int i=0; i<3; i++){
	  double Ru_i=0.0;
	  if (i!=1){ //remove Y component
	    Ru_i+=-N[n]*CollarPressure*k.normal[i]; 
	  }
	  R[n*dof+i] = Ru_i*k.J_a;
	}
      }
    }
  }
  else{
    //
    if (bvp->type!=3){ //non-pullout BVP
      if (std::abs(pCoords[1])<0.5*bvp->l){ 
	Residual(p, 0, 0, 0, U, 0, U, R, ctx);
	for (unsigned int n=0; n<(unsigned int)(nen); n++) R[n*dof+1]=0.0; //null the y component 
      }
      else{
	for (unsigned int n=0; n<(unsigned int)(nen*dof); n++) R[n]=0.0; 
      }
    }
    else{
      Residual(p, 0, 0, 0, U, 0, U, R, ctx); //This is for pullout BVP
    }
  }
  //
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FunctionUAtBase"
PetscErrorCode FunctionUAtBase(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx)
{
  PetscFunctionBegin;
  BVPStruct *bvp = (BVPStruct *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal pCoords[3]; IGAPointFormGeomMap(p,pCoords);

   //get kinematic quantities
  KinematicsStruct<PetscScalar> k;
  getKinematics<PetscScalar>(p, U, U, k);

  //get Xmin (the coordinate at the maximum value of pinching)
  double x=std::sqrt(k.x0[0]*k.x0[0]+k.x0[2]*k.x0[2]); //radial distance from the central cylindrical axis
  bvp->xMin=std::min(bvp->xMin, x);
  
  //
#ifdef LagrangeMultiplierMethod
  double (*u)[3+1] = (double (*)[3+1])U;
#else
  double (*u)[3] = (double (*)[3])U;
#endif
  
  //
  const PetscReal (*N) = (const PetscReal (*)) p->shape[0];
  if (std::abs(pCoords[1])<0.5*bvp->l){ 
    for (unsigned int n=0; n<(unsigned int)(nen); n++){
      	R[n*dof+0] = N[n]*u[n][0]*k.J_a;
	R[n*dof+1] = N[n]*u[n][1]*0; //null the y component as this by itself can be higher than the x,z components
	R[n*dof+2] = N[n]*u[n][2]*k.J_a;
#ifdef LagrangeMultiplierMethod   
      R[n*dof+3] = 0.0;
#endif
    }
  }
  else{
    for (unsigned int n=0; n<(unsigned int)(nen*dof); n++) R[n]=0.0; 
  }
  //
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "JacobianL2"
PetscErrorCode JacobianL2(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  BVPStruct *bvp = (BVPStruct *)ctx;
  
  PetscInt dim = p->dim;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //get kinematic quantities
  KinematicsStruct<PetscScalar> k;
  getKinematics<PetscScalar>(p, U, U, k);
  
  const PetscReal *N = (const PetscReal (*)) p->shape[0];

  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
	  PetscReal val2=0.0;
	  if (d1==d2) {val2 = N[n1] * N[n2];}
	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] =val2*k.J_a;
	}
      }
    }
  }

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "ProjectL2"
PetscErrorCode ProjectL2(void *ctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  BVPStruct *bvp = (BVPStruct *)ctx;
  
  //Solve L2 projection problem
  Mat A;
  Vec b;
  ierr = VecSet(bvp->xDirichlet,0.0);CHKERRQ(ierr);
  ierr = IGACreateMat(bvp->iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(bvp->iga,&b);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGASetFormFunction(bvp->iga,FunctionDirichletL2,bvp);CHKERRQ(ierr);
  ierr = IGASetFormJacobian(bvp->iga,JacobianL2,bvp);CHKERRQ(ierr);
  ierr = IGAComputeFunction(bvp->iga,bvp->xDirichlet,b);CHKERRQ(ierr);
  ierr = IGAComputeJacobian(bvp->iga,bvp->xDirichlet,A);CHKERRQ(ierr);
  //Solver
  {
    KSP ksp;
    ierr = IGACreateKSP(bvp->iga,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,bvp->xDirichlet);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  }
  /*
  char           filename[256];
  sprintf(filename,"./UU%d.vts",bvp->load_increment);
  ierr = IGADrawVecVTK(bvp->iga,bvp->xDirichlet,filename);CHKERRQ(ierr);
  */
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ProjectFields"
PetscErrorCode ProjectFields(Vec& U, void *ctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  BVPStruct *bvp = (BVPStruct *)ctx;

  //
  ierr = IGASetFixTable(bvp->iga,NULL);CHKERRQ(ierr); /* Clear vector to read BCs from */
  //clear old BC's
  IGAForm form;
  ierr = IGAGetForm(bvp->iga,&form);CHKERRQ(ierr);
  for (PetscInt dir=0; dir<2; dir++){
    for (PetscInt side=0; side<2; side++){
      ierr =   IGAFormClearBoundary(form,dir,side);
    }
  }
  
  //Solve L2 projection problem
  Mat A;
  Vec X, R, UAtBase;
  Vec bx, br, bUAtBase;
  ierr = IGACreateMat(bvp->iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(bvp->iga,&X);CHKERRQ(ierr);
  ierr = IGACreateVec(bvp->iga,&R);CHKERRQ(ierr);
  ierr = IGACreateVec(bvp->iga,&UAtBase);CHKERRQ(ierr);
  ierr = IGACreateVec(bvp->iga,&bx);CHKERRQ(ierr);
  ierr = IGACreateVec(bvp->iga,&br);CHKERRQ(ierr);
  ierr = IGACreateVec(bvp->iga,&bUAtBase);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  //Jacobian
  ierr = IGASetFormJacobian(bvp->iga,JacobianL2,bvp);CHKERRQ(ierr);
  ierr = IGAComputeJacobian(bvp->iga,U,A);CHKERRQ(ierr);
  //Fields
  ierr = IGASetFormFunction(bvp->iga,FunctionFields,bvp);CHKERRQ(ierr);
  ierr = IGAComputeFunction(bvp->iga,U,bx);CHKERRQ(ierr);
  //Reactions
  ierr = IGASetFormFunction(bvp->iga,FunctionReactions,bvp);CHKERRQ(ierr);
  ierr = IGAComputeFunction(bvp->iga,U,br);CHKERRQ(ierr);
  //UAtBase
  ierr = IGASetFormFunction(bvp->iga,FunctionUAtBase,bvp);CHKERRQ(ierr);
  ierr = IGAComputeFunction(bvp->iga,U,bUAtBase);CHKERRQ(ierr);
  //Solver
  {
    KSP ksp;
    ierr = IGACreateKSP(bvp->iga,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,bx,X);CHKERRQ(ierr); //Project fields
    ierr = KSPSolve(ksp,br,R);CHKERRQ(ierr); //Project reactions
    ierr = KSPSolve(ksp,bUAtBase,UAtBase);CHKERRQ(ierr); //Project UAtBase
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  }
  //write to vtk files
  char           filename[256];
  //write fields
  sprintf(filename,"./fields%d.vts",bvp->load_increment);
  DMDASetFieldName(bvp->iga->draw_dm,0,"Ux"); VecStrideScale(X, 0, 1.0/bvp->lengthFactor);
  DMDASetFieldName(bvp->iga->draw_dm,1,"Uy"); VecStrideScale(X, 1, 1.0/(bvp->lengthFactor*bvp->lengthFactor));
  DMDASetFieldName(bvp->iga->draw_dm,2,"Uz"); VecStrideScale(X, 2, 1.0);
#ifdef LagrangeMultiplierMethod
  DMDASetFieldName(bvp->iga->draw_dm,3,"H"); VecStrideScale(X, 3, 1.0);
#endif
  ierr = IGADrawVecVTK(bvp->iga,X,filename);CHKERRQ(ierr);
  //write reactions
  sprintf(filename,"./reactions%d.vts",bvp->load_increment);
  DMDASetFieldName(bvp->iga->draw_dm,0,"Rx");
  DMDASetFieldName(bvp->iga->draw_dm,1,"Ry");
  DMDASetFieldName(bvp->iga->draw_dm,2,"Rz");
#ifdef LagrangeMultiplierMethod
  DMDASetFieldName(bvp->iga->draw_dm,3,"none");
#endif
  VecScale(R, bvp->forceFactor);
  ierr = IGADrawVecVTK(bvp->iga,R,filename);CHKERRQ(ierr);
  
  //
  double uVal=0.0, rVal=0.0;
  VecNorm(UAtBase,NORM_INFINITY,&uVal);
  VecNorm(R,NORM_INFINITY,&rVal);

  //
#ifdef enableForceControl
  if (bvp->type!=3){  //non-pullout BVP
    PetscReal xMin=bvp->xMin;
    MPIU_Allreduce(&xMin,&xMin,1,MPIU_REAL,MPIU_MIN,PetscObjectComm((PetscObject)bvp->iga->draw_dm));
    uVal=xMin;
    rVal=bvp->CollarPressure*bvp->CollarHeight; //Force per unit length (pressure*thickness)
  }
  else{ //pullOut BVP
    VecStrideMax(U,1,NULL,&uVal); 
    rVal=bvp->tractionOnTop*2*3.142*(bvp->l/100.0); //Force (traction*circumference) 
  }
#else
  if (bvp->type!=3){
    uVal=bvp->l-bvp->uDirichlet;
    VecStrideMax(R,0,NULL,&rVal); //non-pullout BVP
  }
  else{
    uVal=bvp->uDirichlet;
    VecStrideMax(R,1,NULL,&rVal); //pull out BVP
  }
#endif
  
  //Get energy values
  PetscScalar    energy[2];// = {0,0};
  ierr = IGAComputeScalar(bvp->iga,U,2,energy,FunctionEnergy,ctx);CHKERRQ(ierr);

  //
  if (bvp->isProc0){
    printf ("Stats: Radius: %14.8e nm, Pressure: %14.8e pN/nm, E1: %14.8e pN-nm, E2: %14.8e pN-nm\n", bvp->lengthFactor*std::abs(uVal), bvp->forceFactor*std::abs(rVal), (double)bvp->energyFactor*energy[0], (double)bvp->energyFactor*energy[1]);
    fprintf (bvp->fileForUROutout, "%12.5e, %12.5e, %12.5e, %12.5e\n", bvp->lengthFactor*std::abs(uVal), bvp->forceFactor*std::abs(rVal), (double)bvp->energyFactor*energy[0], (double)bvp->energyFactor*energy[1]);
    fflush(bvp->fileForUROutout);
  }
   
  ierr = VecDestroy(&bx);CHKERRQ(ierr);
  ierr = VecDestroy(&br);CHKERRQ(ierr);
  ierr = VecDestroy(&bUAtBase);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&R);CHKERRQ(ierr);
  ierr = VecDestroy(&UAtBase);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#endif
