import sys, numpy as np
from igakit.nurbs import NURBS
#from igakit.cad import circle, line, join, revolve
from cad3 import circle, line, join, revolve
import scipy.io
import math
from math import pi as Pi
from numpy import linspace

numPointsX=10;
tol=0.001;
kX=[];
C = circle();
L=linspace(0,Pi,numPointsX)[1:-1];
L=np.append(L.tolist(),(L+Pi).tolist());
for t in L:
    print t*(180/Pi);
    xR=math.cos(t); yR=math.sin(t);
    for i in linspace(0,1.0,10001):
        x=C([i])[0][0]; y=C([i])[0][1];
        dist = math.hypot(x-xR, y-yR);
        if (dist<=tol):
            kX.append(i);
            #print t*(180/Pi), i, x, y;
            break;
print "Num knots:", len(kX);
print kX;

#Read from mat file
mat = scipy.io.loadmat('bsplineTube.mat')

order = np.array(mat['order'])
order = order.tolist()[0][0];

knots = np.array(mat['knots'])
knots = knots.tolist()[0];

C = np.transpose(np.array(mat['controlPoints']))
U=knots;
Curve = NURBS([U],C)
S = revolve(Curve, (0,0), 1)

#refine along X
#to_insertX = np.setdiff1d(linspace(0.25,0.5,11)[1:-1],S.knots[0]);
#S.refine(0,to_insertX)

#refine along Y
#to_insertY = np.setdiff1d(linspace(0,1.0,21)[1:-1],S.knots[1]);
#S.elevate(0,1);
S.remove(0, 0.052631658, 1, 1.0e-3);
print S.knots[0];
#S.elevate(1,1);
#S.refine(1,kX);
#S.refine(1,kX);
S.unclamp(1, continuity=1);
#S.insert(1,0.5,1)
#S.rotate(1,0.5*Pi)


from igakit.io import PetIGA
PetIGA().write("mesh.dat", S, nsd=3)


