import sys, numpy as np
from igakit.nurbs import NURBS
#from igakit.cad import circle, line, join, revolve
from cad3 import circle, line, join, revolve
import scipy.io
import math
from math import pi as Pi
from numpy import linspace

#BSpline files
circleFileName='circle160.mat';
splineFileName='tube640.mat'
outputFileName='tubeMeshr160h640C1.dat'
C2Continuity=False;

#Read from mat file
mat = scipy.io.loadmat(splineFileName)

order = np.array(mat['order'])
order = order.tolist()[0][0];
knots = np.array(mat['knots'])
knots = knots.tolist()[0];

C = np.transpose(np.array(mat['controlPoints']))
U=knots;
Curve = NURBS([U],C)
S = revolve(Curve, circleFileName, (0,0), 1)
print S.knots[0];
print S.knots[1];

#refine along X
#to_insertX = np.setdiff1d(linspace(0.25,0.5,11)[1:-1],S.knots[0]);
#S.refine(0,to_insertX)

#refine along Y
#to_insertY = np.setdiff1d(linspace(0,1.0,21)[1:-1],S.knots[1]);
#S.elevate(0,1);

#S.elevate(1,1);
#S.refine(1,kX);
#S.refine(1,kX);
#S.insert(1,0.5,1)
#S.rotate(1,0.5*Pi)
if C2Continuity:
    S.elevate(0,1); S.elevate(1,1);
    S.unclamp(1, continuity=2);
else:
    S.unclamp(1, continuity=1);

#
from igakit.io import PetIGA
PetIGA().write(outputFileName, S, nsd=3)


