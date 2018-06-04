#include <math.h> 
#include "petiga.h"
#include <Sacado.hpp>

PetscErrorCode System(IGAPoint p,PetscScalar *K,PetscScalar *F,void *ctx)
{
  PetscInt nen = p->nen;
  PetscInt dim = p->dim;
  PetscInt dof=1;
  //
  PetscReal l=1.0;
  PetscReal CollarRadius=l;
  PetscReal CollarHeight=10*l;
  //PetscReal CollarZ=1.005*CollarHeight; //1.015*CollarHeight; //Cap
  //PetscReal CollarZ=0.50*CollarHeight;  //Tube
  PetscReal CollarZ=0.035*CollarHeight;  //Base
  PetscReal CollarDepth=0.0025*CollarHeight;
  PetscReal CollarHelixHeight=2*CollarDepth;
  //
  PetscReal pCoords[3];
  IGAPointFormGeomMap(p,pCoords);
  //
  const PetscReal *N = (PetscReal *) p->shape[0];
  PetscReal *normal = p->normal;
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
      	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] = N[n1] * N[n2];
	}
      }
      PetscReal FVal=0.0;
      //
      bool isCollar=false;
      //collar implementation
      if (std::abs(pCoords[1]-CollarZ)<=CollarDepth) {
	isCollar=true;
      }
      //helix implementation
      /*
      PetscReal cc=CollarHelixHeight/(2*3.1415);
      PetscReal tt=(pCoords[1]-CollarZ)/cc;
      if ((tt>=0) && (tt<=2*3.1415)){ 
	PetscReal xx=CollarRadius*std::cos(tt);
	PetscReal yy=CollarRadius*std::sin(tt);
	if (std::sqrt(std::pow(pCoords[0]-xx,2)+std::pow(pCoords[2]-yy,2))<=CollarDepth) {isCollar=true;}
      }
      */
      if (isCollar){ 
	FVal=N[n1]*1.0; //normal[d1];
	//std::cout << FVal << ", ";
      }
      F[n1*dof+d1]= FVal;
    }
  }
  return 0;
}


int main(int argc, char *argv[]) {

  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDof(iga,1);CHKERRQ(ierr); 
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetGeometryDim(iga,3);CHKERRQ(ierr);
  //ierr = IGARead(iga,filename);CHKERRQ(ierr);  
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  
  Mat A;
  Vec x,b;
  ierr = IGACreateMat(iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&x);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&b);CHKERRQ(ierr);
  ierr = IGASetFormSystem(iga,System,0);CHKERRQ(ierr);
  ierr = IGAComputeSystem(iga,A,b);CHKERRQ(ierr);
  ierr = VecSet(x,0.0);CHKERRQ(ierr);
  
  KSP ksp;
  ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  PetscReal bNorm; VecNorm(b,NORM_2,&bNorm);
  PetscReal xNorm; VecNorm(x,NORM_2,&xNorm);
  std::cout << bNorm << ", " << xNorm << std::endl;
  
  ierr = IGADrawVecVTK(iga, x,"reactionsX.vts");CHKERRQ(ierr);
  ierr = IGADrawVecVTK(iga, b,"reactionsB.vts");CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
