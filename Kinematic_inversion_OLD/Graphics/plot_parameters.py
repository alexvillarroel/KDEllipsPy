#!/usr/bin/env python3

# Main Python imports:
import math
import numpy as np
import matplotlib.pyplot as plt



#########################################################
#   SCRIPT TO PLOT THE CONVERGENCE OF EACH ONE OF THE   #
#                  INVERTED PARAMETERS                  #
#########################################################



# ======================================================================
#  INVERSION AND MODEL DATA
# ======================================================================
# Retrieve the number of model parameters:
param_file = open("../Event/kine_files/kine_param", 'r')
i = 0
for row in param_file:
    cols = row.split()
    if i == 0:
        num_params = int(cols[0])
        i = i+1

# Load models and parameters tested in the inversion:
file = "../Event/kine_files/models.dat"
print("File =", file)
models = np.loadtxt(file)
param_size = np.size(models)
num_models = int(param_size/num_params)
iter = np.arange(num_models)
print("There are %i models tested in total" % (num_models))
models.shape = num_models,num_params

# Load misfit values:
misfits = np.loadtxt("../Event/kine_files/misfits.dat")
res_min = misfits.min()
index_min = misfits.argmin()

print("First model =", models[0,:])
print("Best model =", models[index_min,:])
print("Last model =", models[num_models-1,:])
print("Minimum misfit = %0.3f" % (res_min))
print("Index of minimum misfit =", index_min)



# ======================================================================
#  FIGURES
# ======================================================================
# Plot length of ellipse axes and ellipse center coordinates:
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(models[:,0],'ro',label="Axis 1")
plt.plot(models[:,1],'bo',label="Axis 2")
plt.title("Ellipse axes")
plt.xlabel("Models")
plt.ylabel("Length (km)")
plt.legend(prop={'size':10})
plt.grid()

plt.subplot(2,1,2)
x0 = models[:,0]*models[:,3]*np.cos(2*math.pi*models[:,4])
y0 = models[:,1]*models[:,3]*np.sin(2*math.pi*models[:,4])
plt.plot(x0[:],'bo',label='x')
plt.plot(y0[:],'ro',label='y')
plt.legend(prop={'size':10})
plt.title("Ellipse center with respect to hypocenter")
plt.xlabel("Models")
plt.ylabel("Location (km)")
plt.grid()
plt.subplots_adjust(hspace=0.50)
plt.savefig("../Figures/Convergence_axes_center.png", dpi=300)


# Plot ellipse rotation angle:
plt.figure(2)
plt.plot(models[:,2]*math.pi*(180/math.pi),'ro')
plt.title("Ellipse rotation angle")
plt.xlabel("Models")
plt.ylabel("Angle (degrees)")
plt.grid()
plt.savefig("../Figures/Convergence_angle.png", dpi=300)


# Plot maximum slip:
plt.figure(3)
plt.plot(models[:,5],'ro')
plt.title("Maximum slip")
plt.xlabel("Models")
plt.ylabel("Maximum slip (m)")
plt.grid()
plt.savefig("../Figures/Convergence_Dmax.png", dpi=300)


# Plot rupture velocity:
plt.figure(4)
plt.plot(models[:,6],'ro')
plt.title("Rupture velocity")
plt.xlabel("Models")
plt.ylabel("Rupture velocity (km/s)")
plt.grid()
plt.savefig("../Figures/Convergence_Vr.png", dpi=300)


# Plot global convergence:
plt.figure(5)
plt.plot(misfits,'ro')
plt.ylim(0,5)
plt.title("Global convergence")
plt.xlabel("Models")
plt.ylabel("Misfit")
plt.grid()
text_str = '\n'.join(("Min. misfit = %0.3f" % (res_min), "Last misfit = %0.3f" % (misfits[-1])))
plt.text(0.98, 0.97, text_str, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.8, edgecolor='k'))
plt.savefig("../Figures/Convergence_misfit.png", dpi=300)



plt.show()
