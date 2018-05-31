import sys, numpy as np
from igakit.nurbs import NURBS
from igakit.cad import circle, line, join, revolve
from igakit.plot import plt
from math import pi as Pi
from igakit.io import PetIGA
from numpy import linspace

#Geometry
R=0.1;
H=R;

C1=line((R, 0, 0), (R, 0, H));
S = revolve(C1, (0,0), 2, angle=[0,0.5*Pi])

#refine along X
to_insertX1 = linspace(0,1.0,5)[1:-1]
S.refine(0,to_insertX1)
S.refine(1,to_insertX1)

#plt.plot(S,color='g')
#plt.show()
PetIGA().write("mesh.dat",S,nsd=3)
