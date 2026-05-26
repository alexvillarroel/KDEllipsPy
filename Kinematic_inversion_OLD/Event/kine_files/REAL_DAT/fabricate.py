#!/usr/bin/env python

from pylab import *
from matplotlib.patches import Ellipse


sismog=open("fake_disp_x",'w')
dt=0.06
x=arange(512)*dt
desplaz=zeros(512)
cero=0.
for i in range(12):
	for j in range(512):
		print >>sismog,cero

