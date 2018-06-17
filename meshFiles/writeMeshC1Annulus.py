import sys, numpy as np
from igakit.nurbs import NURBS
#from igakit.cad import circle, line, join, revolve
from cad2 import circle, line, join, revolve, ruled
import scipy.io
import math
from math import pi as Pi
from numpy import linspace
from igakit.io import PetIGA

c1 = circle()
c2 = circle(radius=2)
srf = ruled(c1,c2).transpose().elevate(0,1)
to_insert = linspace(0,0.25,5)[1:-1]
for i in range(4):
    srf.refine(0,i*0.25+to_insert)
    srf.refine(1,i*0.25+to_insert)
srf.unclamp(1) # <----------------------- added this line
PetIGA().write("annulus.dat",srf)

