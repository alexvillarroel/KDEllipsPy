#!/usr/bin/env python3

# Main Python imports:
import numpy as np
import matplotlib.pyplot as plt



#############################################################
#   SCRIPT TO PLOT THE SLIP DISTRIBUTION AND RUPTURE TIME   #
#            INSIDE THE ELLIPTICAL RUPTURE ZONE             #
#############################################################



# ======================================================================
#  FAULT DATA
# ======================================================================
# Load fault discretization data:
input = open("../Evento/input.dat",'r')
(nx,ny,nz) = input.readline().split()[0:3]
nx = int(nx)
ny = int(ny)
(nt,) = input.readline().split()[0:1]
nt = int(nt)
(dx,) = input.readline().split()[0:1]
dx = float(dx)/1000.

xmax = int(nx*dx)
ymax = int(ny*dx)

print("Grid: %d x %d points" % (nx, ny))
print("dx = %i (m)" % (dx*1000))

# Read hypocenter location on the fault:
input2 = open("../Faille/faille.in",'r')
input2.readline()
(xhypo,yhypo) = input2.readline().split()[0:2]




# ======================================================================
#  SLIP DISTRIBUTION FIGURE
# ======================================================================
# Load slip data:
A = np.loadtxt("../Evento/slip.res")
max_slip = A.max()
print("Maximum slip = %0.2f (m)" % (max_slip))

# Calculate the average slip:
sum_slip = []
for i in range(0,len(A)):
   sample = A[i]
   # The threshold of 0.1 was obtained from visual inspection (points 
   # outside ellipse are close to zero, but not zero):
   if sample > 0.1:
      sum_slip = np.append(sum_slip, sample)
average_slip = np.mean(sum_slip)
print("Average slip = %0.2f (m)" % (average_slip))

# Fracture energy (according to Poli & Prieto (2016)):
Gc_1=7.63-(0.03*average_slip**(1.72-0.03))
Gc_2=7.63+(0.03*average_slip**(1.72+0.03))

# Plot figure:
plt.figure(1,figsize=(8,7))
A.shape = nx,ny
plt.imshow(A.T,vmin=0., vmax=max_slip+0.1,interpolation='bilinear',cmap='binary')
plt.xticks(np.arange(2)*nx,('0',xmax), fontsize=27)
plt.yticks(np.arange(2)*ny,('',ymax), fontsize=27)
plt.xlim(0.,nx)
plt.title("Coseismic slip", fontsize=27)
plt.ylabel("Along strike (km)", fontsize=27)
plt.xlabel("Along dip (km)", fontsize=27)
cbar = plt.colorbar(ticks=[0, 0.3, 0.6, 0.9, 1.2, 1.5])
cbar.ax.set_title('(m)', fontsize=27)
cbar.ax.tick_params(labelsize=27)
plt.gca().set_aspect('equal')
plt.gca().invert_yaxis()
plt.savefig("../Figures/Source_Slip.pdf")



# ======================================================================
#  RUPTURE TIME FIGURE
# ======================================================================
# Load rupture time data:
A = np.loadtxt("../Evento/ruptime.res")
max_rt=A.max()

# Plot figure:
plt.figure(2,figsize=(8,7))
print("Rupture time = %0.2f (s)" % (max_rt))
A.shape = nx,ny
plt.imshow(A.T,vmin=0.,vmax=max_rt+0.1,interpolation='bilinear',cmap='binary')
plt.xticks(np.arange(2)*nx,('0',xmax), fontsize=12)
plt.yticks(np.arange(2)*ny,('',ymax), fontsize=12)
plt.xlim(0.,nx)
plt.title("Rupture time")
plt.ylabel("Along strike (km)")
plt.xlabel("Along dip (km)")
cbar = plt.colorbar()
cbar.ax.set_title('(s)')
plt.gca().set_aspect('equal')
plt.gca().invert_yaxis()
plt.savefig("../Figures/Source_Rupture_Time.png", dpi=300)



plt.show()
