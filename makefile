fileName=main
rm *.o $fileName

TRILINOS_DIR=/Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5
#TRILINOS_DIR=/home/krudraraju/software/trilinos/trilinos-12.12.1-Source/install
#TRILINOS_DIR=/Users/rudraa/workspace/software/trilinos/trilinos-12.12.1-Source/installDir

#Compile
mpicxx -c -O2 -fPIC  -std=c++11 -Wall -Wwrite-strings -Wno-unused-variable -Wno-unused-value  -Wno-uninitialized -Wno-strict-aliasing -Wno-unknown-pragmas -I/$PETSC_DIR/$PETSC_ARCH/include -I/$PETSC_DIR/include -I$PETIGA_DIR/$PETSC_ARCH/include -I$PETIGA_DIR/include -I$TRILINOS_DIR/include/ $fileName.c 
mpicxx -o $fileName $fileName.o -Wl,-rpath,$PETIGA_DIR/$PETSC_ARCH/lib -L$PETIGA_DIR/$PETSC_ARCH/lib -lpetiga -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib  -lpetsc -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib $TRILINOS_DIR/lib//libteuchoscomm.dylib	 $TRILINOS_DIR/lib//libteuchoscore.dylib $TRILINOS_DIR/lib//libteuchosnumerics.dylib $TRILINOS_DIR/lib//libteuchosparameterlist.dylib $TRILINOS_DIR/lib//libteuchosremainder.dylib

#mpicxx -c -O2 -fPIC -Wall -Wwrite-strings -Wno-unused-variable -Wno-unused-value  -Wno-uninitialized -Wno-strict-aliasing -Wno-unknown-pragmas -I/$PETSC_DIR/$PETSC_ARCH/include -I/$PETSC_DIR/include -I$PETIGA_DIR/$PETSC_ARCH/include -I$PETIGA_DIR/include -I$TRILINOS_DIR/include/ $fileName.c
#mpicxx -o $fileName $fileName.o -Wl,-rpath,$PETIGA_DIR/$PETSC_ARCH/lib -L$PETIGA_DIR/$PETSC_ARCH/lib -lpetiga -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib -lpetsc -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib $TRILINOS_DIR/lib//libteuchoscomm.so $TRILINOS_DIR/lib//libteuchoscore.so $TRILINOS_DIR/lib//libteuchosnumerics.so $TRILINOS_DIR/lib//libteuchosparameterlist.so $TRILINOS_DIR/lib//libteuchosremainder.so

mpiexec -np 4 ./$fileName -iga_view -ts_monitor -snes_monitor -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -ts_adapt_type none -ts_max_snes_failures 500 -snes_max_it 200 -snes_max_funcs 50000 -snes_type newtontr
#  -snes_linesearch_type nleqerr  
#-snes_atol 1.0e-12 -snes_rtol 1.0e-10 -snes_stol 1.0e-12 

