import sys, numpy as np
from igakit.nurbs import NURBS
from igakit.cad import circle, line, join, revolve, ruled, extrude
from igakit.plot import plt
from math import pi as Pi
from igakit.io import PetIGA
from numpy import linspace

#Geometry
Lx=3.0;
Ly=1.0; #This is the length scale of the problem

c1=line((0, 0, 0), (Lx, 0, 0));
S=extrude(c1, displ=Ly,axis=1);
S.elevate(0,1)
S.elevate(1,1)

#refine along X
to_insertX = np.setdiff1d(linspace(0,1.0,9)[1:-1],S.knots[0]);
S.refine(0,to_insertX)

#output
PetIGA().write("mesh.dat",S,nsd=3)
