#!/usr/bin/env python3

# Main Python imports:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse



#############################################################
#   SCRIPT TO PLOT THE SLIP DISTRIBUTION AND RUPTURE TIME   #
#            INSIDE THE ELLIPTICAL RUPTURE ZONE             #
#############################################################



# ===============================================
#  LOAD MODEL DATA:
# ===============================================
# Read kine.in file:
file_in = open("../Event/kine_files/kine.in",'r')
for i in range(5):file_in.readline()
(nsource,) = file_in.readline().split()[0:1]
(xhypo,yhypo,nstr,ndip) = file_in.readline().split()[0:4]
nstr = int(nstr)
ndip = int(ndip)
nsource = int(nsource)

# Read faille.in file to retrieve fault size:
faille_file = open("../Faille/faille.in",'r')
(xsize,ysize) = faille_file.readline().split()[0:2]
nx1 = int(xsize)/1000
ny1 = int(ysize)/1000
print("Fault size:", nx1, "x", ny1, "(km)")
print('Grid:', nstr, "x", ndip, "sub-faults")

# Read axi.hist file to adjust hypocenter coordinates:
hist = np.loadtxt("../Event/axi.hist")
hist.shape = nsource,8
dx = float(hist[1,5])/4.
dy = float(hist[1,6])/4.
xhypo = float(xhypo)/dx
yhypo = float(yhypo)/dy

# Loading slip data:
A = np.loadtxt("../Event/kine_files/slip_fin.dat")
max_slip = A.max()
print("Maximum slip = %0.2f (m)" % (max_slip))

# Arrange data to get the spatial distribution:
nx = nstr*4
ny = ndip*4
A.shape = ny,nx

# Loading rupture time data:
B = np.loadtxt("../Event/kine_files/tr_fin.dat")
max_rt = B.max()
print("Rupture time = %0.2f (s)" % (max_rt))

# Arrange data to get the spatial distribution:
B.shape = ny,nx



# ======================================================
#  FIGURES
# ======================================================
# Plot slip distribution:
fig = plt.figure(1,figsize=(8,7))
ax = fig.add_subplot(111)
img = plt.imshow(A,vmin=0,vmax=A.max(),interpolation='bilinear',cmap='binary')
el = Ellipse((xhypo, yhypo), 2, 2)
el.set_facecolor('white')
ax.add_patch(el)
plt.xticks(np.arange(2)*nx,('0',nx1),fontsize=12)
plt.yticks(np.arange(2)*ny,(ny1,''),fontsize=12)
plt.xlim(0,nx)
cbar = plt.colorbar()
cbar.ax.set_title('(m)')
plt.title("Slip distribution")
plt.xlabel("Along strike (km)")
plt.ylabel("Along dip (km)")
plt.subplots_adjust(left=0.08,right=1.00,bottom=0.08,top=0.92)
plt.savefig("../Figures/Source_slip_distribution.png", dpi=300)


# Plot rupture time:
fig = plt.figure(2,figsize=(8,7))
ax = fig.add_subplot(111)
img = plt.imshow(B,vmin=0.,vmax=B.max(),interpolation='bilinear',cmap='binary')
el = Ellipse((xhypo, yhypo), 2, 2)
el.set_facecolor('white')
ax.add_patch(el)
plt.xticks(np.arange(2)*nx,('0',nx1),fontsize=12)
plt.yticks(np.arange(2)*ny,(ny1,''),fontsize=12)
plt.xlim(0,nx)
cbar = plt.colorbar()
cbar.ax.set_title('(s)')
plt.title("Rupture time")
plt.xlabel("Along strike (km)")
plt.ylabel("Along dip (km)")
plt.subplots_adjust(left=0.08,right=1.00,bottom=0.08,top=0.92)
plt.savefig("../Figures/Source_rupture_time.png", dpi=300)



plt.show()
