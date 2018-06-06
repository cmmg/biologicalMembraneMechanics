/*
  Write mesh in dat format to vtk format
  author: rudraa
 */
#include <math.h> 
//extern "C" {
#include "petiga.h"
//}

#define DIM 3
#define DOF 3

int main(int argc, char *argv[]) {

  char           filename[PETSC_MAX_PATH_LEN] = "mesh.dat";
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDof(iga,3);CHKERRQ(ierr); // dofs = {ux,uy,uz}
  ierr = IGASetDim(iga,2);CHKERRQ(ierr);
  ierr = IGASetGeometryDim(iga,3);CHKERRQ(ierr);
  ierr = IGARead(iga,filename);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  Vec U;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = VecSet(U,0.0);CHKERRQ(ierr);
  ierr = IGAWriteVec(iga,U,"mesh.out");CHKERRQ(ierr);
  ierr = IGADrawVecVTK(iga,U,"mesh.vts");CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
