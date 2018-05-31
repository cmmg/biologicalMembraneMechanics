import sys, numpy as np
from igakit.nurbs import NURBS

N = 15
p = 2
if len(sys.argv) >= 2:
    N = int(sys.argv[1])
if len(sys.argv) >= 3:
    p = int(sys.argv[2])

U = [0,0,     1,1]
V = [0,0,0, 1,1,1]
C = np.zeros((2,3,4))
val = np.sqrt(2)*0.5
C[0,0,:] = [  0,-100,  0,1]
C[1,0,:] = [100,-100,  0,1]
C[0,1,:] = [  0,-100,100,1]
C[1,1,:] = [100,-100,100,1]
C[0,2,:] = [  0,   0,100,1]
C[1,2,:] = [100,   0,100,1]
C[:,1,:] *= val

geom = NURBS([U,V],C)
geom.elevate(0,max(p-1,0)).elevate(1,max(p-2,0))

h = 1./N
insert = np.linspace(h,1.-h,N-1)
geom.refine(0,insert).refine(1,insert)

if True:
    from igakit.io import PetIGA
    PetIGA().write("mesh.dat", geom, nsd=3)
