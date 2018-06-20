/*
  Kirchhoff-Love shell implementation in PetIGA
  author: rudraa
 */

#undef __FUNCT__
#define __FUNCT__ "FunctionL2"
PetscErrorCode FunctionL2(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal x[3];
  IGAPointFormGeomMap(p,x);

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
	
  PetscReal uDirichletVal=-1.0*user->l*user->c_time;  
  if (x[1]>(10.5*user->l)){ uDirichletVal=0.0;}//for top surface
  //std::cout << "(," << uDirichletVal << ", " << x[1] << "), ";

  if (!user->projectBC){
    Residual(p,0.0,0,0.0,U,0.0,0,R,mctx);
  }
  else{
    //store L2 projection residual
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
	  val=1.0; break;
	}
	R[n1*dof+d1] = N[n1]*val;
      }
    }
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "JacobianL2"
PetscErrorCode JacobianL2(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  
  PetscInt dim = p->dim;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  
  const PetscReal *N = (const PetscReal (*)) p->shape[0];

  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
	  PetscReal val2=0.0;
	  if (d1==d2) {val2 = N[n1] * N[n2];}
	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] =val2;
	}
      }
    }
  }

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "ProjectL2"
PetscErrorCode ProjectL2(IGA iga, PetscInt step, Vec U, void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;

  //clear old BC's
  IGAForm form;
  ierr = IGAGetForm(user->iga,&form);CHKERRQ(ierr);
  for (PetscInt dir=0; dir<2; dir++){
    for (PetscInt side=0; side<2; side++){
      ierr =   IGAFormClearBoundary(form,dir,side);
    }
  }
  /*
  char           filename[256];
  sprintf(filename,"./project%d.vts",step);
  user->projectBC=false;
  Vec x2;
  ierr = IGACreateVec(user->iga,&x2);CHKERRQ(ierr);
  ierr = VecSet(x2,0.0);CHKERRQ(ierr);
  */
  /* Solve L2 projection problem */
  Mat A;
  Vec b;
  ierr = VecSet(user->x,0.0);CHKERRQ(ierr);
  ierr = IGACreateMat(user->iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(user->iga,&b);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGASetFormFunction(user->iga,FunctionL2,user);CHKERRQ(ierr);
  ierr = IGASetFormJacobian(user->iga,JacobianL2,user);CHKERRQ(ierr);
  ierr = IGAComputeJacobian(user->iga,user->x,A);CHKERRQ(ierr);
  /*
  //Solver
  {
    KSP ksp;
    ierr = IGACreateKSP(user->iga,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = IGAComputeFunction(user->iga,U,b);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,x2);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  }
  //
  ierr = IGADrawVecVTK(user->iga,x2,filename);CHKERRQ(ierr);
  */
  //
  user->projectBC=true;
  //
  //Solver
  {
    KSP ksp;
    ierr = IGACreateKSP(user->iga,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = IGAComputeFunction(user->iga,user->x,b);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,user->x);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  }
  
  //
  //PetscReal xVal;
  //VecNorm(user->x,NORM_INFINITY,&xVal);
  //std::cout << "xVal: " << xVal << "\n";
  //
  ierr = IGASetBoundaryValue(user->iga,0,0,0,0.0);CHKERRQ(ierr); 
  ierr = IGASetBoundaryValue(user->iga,0,0,2,0.0);CHKERRQ(ierr);
  //ierr = IGASetBoundaryValue(user->iga,0,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,0,1,0,/*dummy*/0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,0,1,1,0.0);CHKERRQ(ierr); 
  ierr = IGASetBoundaryValue(user->iga,0,1,2,/*dummy*/0.0);CHKERRQ(ierr);
  ierr = IGASetFixTable(user->iga,user->x);CHKERRQ(ierr);    /* Set vector to read BCs from */
  //
  //ierr = VecDestroy(&x2);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  PetscFunctionReturn(0);
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
  ierr = IGASetFixTable(user->iga,NULL);CHKERRQ(ierr); /* Clear vector to read BCs from */
  user->c_time=c_time;
  ProjectL2(user->iga, it_number, U, mctx);
  //
  //ierr = IGASetBoundaryValue(user->iga,0,0,2,user->l*(c_time));CHKERRQ(ierr); //Y=t on \eta_2=0
  PetscFunctionReturn(0);
}

#endif
