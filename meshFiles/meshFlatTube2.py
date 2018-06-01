import sys, numpy as np
from igakit.nurbs import NURBS
from igakit.cad import circle, line, join, revolve, ruled
from igakit.plot import plt
from math import pi as Pi
from igakit.io import PetIGA
from numpy import linspace

#Geometry
R=0.01;

c1=line((0, 0, 0), (100*R, 0, 0));
S = revolve(c1, (0,0), 2, angle=[0,2*Pi])
S.elevate(0,1)

#refine along X
to_insertX = np.setdiff1d(linspace(0,1.0,11)[1:-1],S.knots[0]);
S.refine(0,to_insertX)

#refine along Y
to_insertY = np.setdiff1d(linspace(0,1.0,21)[1:-1],S.knots[1]);
S.refine(1,to_insertY)

#S.elevate (1,1);

#periodicity
S.unclamp(1, continuity=1)

#plt.plot(S,color='g')
#plt.show()
PetIGA().write("mesh.dat",S,nsd=3)
