#!/usr/bin/env python3

# Main Python imports:
import numpy as np
import matplotlib.pyplot as plt



#########################################################
#   SCRIPT TO PLOT THE SLIP RATE AND STRESS DROP RATE   #
#########################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
dt = 0.25   # Time increment for time snaapshots in plots
x1 = 2.5    # x-point on fault to obtain dynamic parameters as function of time
y1 = 2.2    # y-point on fault to obtain dynamic parameters as function of time
cols = 5    # Number of columns for time-snapshots plots



# ======================================================================
#  DYNAMIC MODEL DATA
# ======================================================================
# Load fault discretization data:
input = open("../Evento/input.dat",'r')
(nx,ny,nz) = input.readline().split()[0:3]
nx = int(nx)
ny = int(ny)

# Load fault information:
input_file = open("../Faille/faille.in",'r')
(nx1,ny1) = input_file.readline().split()[0:2]
nx1 = int(nx1)/1000
ny1 = int(ny1)/1000

# Load dynamic model data:
A = np.loadtxt("../Evento/result/sliprate.res")
B = np.loadtxt("../Evento/result/disp.res")
S = np.loadtxt("../Evento/result/ssx3d.res")

print("Fault discretization =", nx, "x", ny, "points")
print("Slip samples =", A.shape)
print("Stress samples =", S.shape)

Amin = A.min()
Amax = A.max()
print("Slip = [%0.2f, %0.2f] (m)" % (Amin, Amax))
Smin = S.min()
Smax = S.max()
print("Stress = [%0.2f, %0.2f] (MPa)" % (Smin/1.E6, Smax/1.E6))

x = np.arange(nx+1)
y = np.arange(ny+1)
X,Y = np.meshgrid(x,y)




# ======================================================================
#  FIGURES
# ======================================================================
# Find the iteration when the slip becomes insignificant (rupture ends):
for i in range(1,40,1):
   data = A[nx*ny*(i-1):nx*ny*i]
   max_data = data.max()
   if max_data < 0.1:
      final_iter = i
      break

# Set the subplots according to the number of iterations:
if ((final_iter-1) % cols == 0):
   nrow_plot = final_iter/cols
else:
   fn = (final_iter-1)/cols
   nrow_plot = int(fn) + 1

# Plot slip distribution in time:
plt.figure(1,figsize=(10,8))
for i in range(1,final_iter,1):
   plt.subplot(nrow_plot,cols,i)
   time = (i-1)*dt
   Picture = A[nx*ny*(i-1):nx*ny*i]
   Picture.shape = ny, nx
   img = plt.imshow(Picture,vmin=0., vmax=0.5*Amax, cmap='binary')
   plt.plot([nx/2],[ny/2],'rx',ms=7)
   # Plot point of measure for Figure 3:
   plt.plot([nx/x1],[ny/y1],'bd',ms=5)
   # Subplot properties:
   plt.xticks(np.arange(2)*nx,('0',nx1), fontsize=10)
   plt.yticks(np.arange(2)*ny,('0',ny1), fontsize=10)
   plt.xlim(0,nx)
   plt.ylim(ny,0)
   plt.title("t = %0.2f (s)" % (time), fontsize=12)
   # Axis labels:
   if (i >= (final_iter-cols) and i < final_iter):
      plt.xlabel("(km)",fontsize=10)
   for j in range(1,final_iter,int(cols)):
      if (i == j):
         plt.ylabel("(km)",fontsize=10)
# Add title, set subplots, and save figure:
t = plt.gcf().text(0.5, 0.95, "Slip distribution in time (m)",horizontalalignment='center',fontsize=14)
plt.subplots_adjust(wspace=0.45,hspace=0.45,left=0.05,right=0.97,bottom=0.06,top=0.89)
plt.savefig("../Figures/Rate_slip.png", dpi=300)


# Plot stress drop in time:
plt.figure(2,figsize=(10,8))
for i in range(1,final_iter,1):
   plt.subplot(nrow_plot,cols,i)
   time = (i-1)*dt
   Picture = S[nx*ny*(i-1):nx*ny*i]/1.E6
   Picture.shape = ny, nx
   img = plt.imshow(Picture, vmin=-0.5*(Smax/1.E6), vmax=0.5*(Smax/1.E6), cmap='binary')
   # Subplot properties:
   plt.xticks(np.arange(2)*nx, ('0',nx1), fontsize=10)
   plt.yticks(np.arange(2)*ny, ('0',ny1), fontsize=10)
   plt.xlim(0,nx)
   plt.ylim(ny,0)
   plt.title("t = %0.2f [s]" % (time),fontsize=12)
   # Axis labels:
   if (i >= (final_iter-cols) and i < final_iter):
      plt.xlabel("[km]",fontsize=10)
   for j in range(1,final_iter,int(cols)):
      if (i == j):
         plt.ylabel("[km]",fontsize=10)
# Add title, set subplots, and save figure:
plt.colorbar()
t = plt.gcf().text(0.5, 0.95, "Stress variations in time [MPa]",horizontalalignment='center',fontsize=14)
plt.subplots_adjust(wspace=0.45,hspace=0.45,left=0.05,right=0.97,bottom=0.06,top=0.89)
plt.savefig("../Figures/Rate_stress_drop.png", dpi=300)


# Reshape data matrices and create time vector:
T = S/1.E7
T.shape = 40, ny, nx
A.shape = 40, ny, nx
B.shape = 40, ny, nx
time = [i*dt for i in range(0,40)]
# Plot dynamic parameters as function of time:
plt.figure(3)
plt.plot(time,T[:,int(nx/x1),int(ny/y1)],'r',label="Stress variations (x10 [MPa])")
plt.plot(time,A[:,int(nx/x1),int(ny/y1)],'b',label="Static field derivative")
plt.plot(time,B[:,int(nx/x1),int(ny/y1)],'k',label="Static field")
plt.legend(prop={'size':10})
plt.title("Parameters")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.savefig("../Figures/Rate_3.png", dpi=300)



plt.show()
