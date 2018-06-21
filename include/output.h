/*
  Output functions
  author: Shiva Rudraraju
 */

#if !defined(OUTPUT_H_)
#define OUTPUT_H_

PetscErrorCode setBCs(BVPStruct& bvp, PetscInt it_number, PetscReal c_time);

#undef __FUNCT__
#define __FUNCT__ "OutputMonitor"
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *ctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  BVPStruct *bvp = (BVPStruct *)ctx;
  char           filename[256];
  sprintf(filename,"./solution%d.vts",it_number);
  ierr = IGADrawVecVTK(bvp->iga,U,filename);CHKERRQ(ierr);
  //
  bvp->c_time=c_time;
  bvp->load_increment=it_number;
  setBCs(*bvp, it_number, c_time);
  //
  PetscFunctionReturn(0);
}

#endif
