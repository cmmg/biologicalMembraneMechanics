import sys, numpy as np
from igakit.nurbs import NURBS
from igakit.cad import circle, line, join, revolve, ruled
from igakit.plot import plt
from math import pi as Pi
from igakit.io import PetIGA
from numpy import linspace

#Geometry
L=1.0

C1=line((0, 0), (0, L));
C2=line((L, 0), (L, L));
S = ruled(C1,C2).transpose()
to_insert = linspace(0,1.0,11)[1:-1];
S.refine(0,to_insert)
S.refine(1,to_insert)

#plt.plot(S,color='g')
#plt.show()
PetIGA().write("mesh.dat",S,nsd=3)
