fileName=main
rm *.o $fileName

TRILINOS_DIR=/Applications/deal.II.app/Contents/Resources/opt/trilinos-0b08cd5

#Compile according to AD type
#Sacado
mpicxx -c -O2 -fPIC  -std=c++11 -Wall -Wwrite-strings -Wno-unused-variable -Wno-unused-value  -Wno-uninitialized -Wno-strict-aliasing -Wno-unknown-pragmas -I/$PETSC_DIR/$PETSC_ARCH/include -I/$PETSC_DIR/include -I$PETIGA_DIR/$PETSC_ARCH/include -I$PETIGA_DIR/include -I$TRILINOS_DIR/include/ $fileName.c 
mpicxx -o $fileName $fileName.o -Wl,-rpath,$PETIGA_DIR/$PETSC_ARCH/lib -L$PETIGA_DIR/$PETSC_ARCH/lib -lpetiga -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib  -lpetsc -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib $TRILINOS_DIR/lib//libteuchoscomm.dylib	 $TRILINOS_DIR/lib//libteuchoscore.dylib $TRILINOS_DIR/lib//libteuchosnumerics.dylib $TRILINOS_DIR/lib//libteuchosparameterlist.dylib $TRILINOS_DIR/lib//libteuchosremainder.dylib

#mpiexec -np 4 ./$fileName -iga_view -ts_monitor -snes_monitor -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -ts_adapt_type none -iga_rule_size 4 -iga_geometry mesh.dat 
#-iga_periodic 0,1
mpiexec -np 4 ./$fileName -iga_view -ts_monitor -snes_monitor -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -ts_adapt_type none -ts_max_snes_failures 500 -snes_max_it 100 -snes_max_funcs 50000 -snes_type newtontr 
#-iga_geometry mesh.dat 
#-iga_periodic 0,1
#-ts_max_snes_failures 500 -snes_max_it 100 -snes_max_funcs 50000 -snes_type newtontr
#-ksp_monitor
#-file_prefix "test" -N 100 -dt 1.0e-5 -ch_monitor -ts_monitor -snes_monitor -snes_converged_reason -log_summary -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps 
#-snes_type newtontr
