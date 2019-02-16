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
  DMDASetFieldName(bvp->iga->draw_dm,0,"Ux");
  DMDASetFieldName(bvp->iga->draw_dm,1,"Uy");
  DMDASetFieldName(bvp->iga->draw_dm,2,"Uz");
#ifdef LagrangeMultiplierMethod
  DMDASetFieldName(bvp->iga->draw_dm,3,"q");
#endif
  ierr = IGADrawVecVTK(bvp->iga,U,filename);CHKERRQ(ierr);
  //
  bvp->c_time=c_time;
  bvp->load_increment=it_number;

  //project and output fields and reaction forces
  ProjectFields(U, ctx);

  //
  //PetscViewer    viewer;
  std::string filehistory("history");
  std::ostringstream temp;
  temp << it_number;
  filehistory.append(temp.str());
  filehistory.append(".dat");
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,filehistory.c_str(),FILE_MODE_WRITE,&viewer);
  //VecView(U,viewer);
  //PetscViewerDestroy(&viewer);
  
  //setup BCs for next load increment 
  setBCs(*bvp, it_number, c_time);
  //
  PetscFunctionReturn(0);
}

#endif
