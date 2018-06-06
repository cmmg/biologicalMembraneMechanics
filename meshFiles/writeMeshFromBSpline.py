import sys, numpy as np
from igakit.nurbs import NURBS
#from igakit.cad import circle, line, join, revolve
from cad2 import circle, line, join, revolve
import scipy.io
from math import pi as Pi
from numpy import linspace

#Read from mat file
mat = scipy.io.loadmat('bsplineBase40.mat')

order = np.array(mat['order'])
order = order.tolist()[0][0];

knots = np.array(mat['knots'])
knots = knots.tolist()[0];

C = np.transpose(np.array(mat['controlPoints']))
U=knots;
Curve = NURBS([U],C)
S = revolve(Curve, (0,0), 1, angle=[0,2*Pi])

#refine along X
#to_insertX = np.setdiff1d(linspace(0.25,0.5,11)[1:-1],S.knots[0]);
#S.refine(0,to_insertX)

#refine along Y
to_insertY = np.setdiff1d(linspace(0,1.0,21)[1:-1],S.knots[1]);
S.refine(1,to_insertY)

#N = 10;
#p = order-1;
#geom.elevate(0,max(p-1,0)).elevate(1,max(p-2,0))
#h = 1./N
#insert = np.linspace(h,1.-h,N-1)
#geom.refine(0,insert).refine(1,insert)

from igakit.io import PetIGA
PetIGA().write("mesh.dat", S, nsd=3)

