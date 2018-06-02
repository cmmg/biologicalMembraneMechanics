import sys, numpy as np
from igakit.nurbs import NURBS
from igakit.cad import circle, line, join, revolve
from igakit.plot import plt
from math import pi as Pi
from igakit.io import PetIGA
from numpy import linspace

#Geometry
R=0.1;
H=10*R;
r=0.5*R;

C1=circle(R, (0,H), (Pi/2,0))
C2=line((R, H), (R, r));
C3=circle(r, (r+R,r), (Pi,3*Pi/2));
C4=line((R+r,0), (R+2*r,0));
C5 = join(C1, C2, 0);
C6 = join(C5, C3, 0);
C7 = join(C6, C4, 0);
S = revolve(C7, (0,0), 1, angle=[0,2*Pi])

#refine along X
to_insertX1 = np.setdiff1d(linspace(0,1.0,9)[1:-1],S.knots[0]);
S.refine(0,to_insertX1)
to_insertX2 = np.setdiff1d(linspace(1.0,2.0,21)[1:-1],S.knots[0]);
S.refine(0,to_insertX2)
to_insertX3 = np.setdiff1d(linspace(2.0,3.0,5)[1:-1],S.knots[0]);
S.refine(0,to_insertX3)
to_insertX4 = np.setdiff1d(linspace(3.0,4.0,5)[1:-1],S.knots[0]);
S.refine(0,to_insertX4)

#refine along Y
to_insertY = np.setdiff1d(linspace(0,1.0,21)[1:-1],S.knots[1]);
S.refine(1,to_insertY)

#periodicity
#S.unclamp(0)
S.unclamp(1, continuity=1)

#plt.plot(S,color='g')
#plt.show()
PetIGA().write("mesh.dat",S,nsd=3)
