import sys, numpy as np
from igakit.nurbs import NURBS
from igakit.cad import circle, line, join, revolve, ruled
from igakit.plot import plt
from math import pi as Pi
from igakit.io import PetIGA
from numpy import linspace

#Geometry
R=0.01;

c1=circle(R, (0,0,1))
c2=circle(100*R, (0,0,1))
#c1 = circle(radius=R)
#c2 = circle(radius=20*R)
S = ruled(c1,c2).transpose().elevate(0,1)

#refine along X
to_insertX = np.setdiff1d(linspace(0,0.25,11)[1:-1],S.knots[0]);
S.refine(0,to_insertX)
to_insertX = np.setdiff1d(linspace(0.25,1.0,5)[1:-1],S.knots[0]);
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
