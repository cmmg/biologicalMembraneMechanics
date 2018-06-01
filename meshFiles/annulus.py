from igakit import cad
from igakit.io import PetIGA
from numpy import linspace
c1 = cad.circle()
c2 = cad.circle(radius=2)
srf = cad.ruled(c1,c2).transpose().elevate(0,1)
to_insert = linspace(0,0.25,11)[1:-1]
for i in range(4):
    srf.refine(0,i*0.25+to_insert)
    srf.refine(1,i*0.25+to_insert)
PetIGA().write("annulus.dat",srf)
