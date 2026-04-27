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
# USER PARAMETERS
# ======================================================================
# Start of model range (reference for histograms):
startmod = 10000



# ======================================================================
#  INVERSION AND MODEL DATA
# ======================================================================
# Load fault discretization data:
input = open("../Evento/input.dat",'r')
(nx,ny,nz) = input.readline().split()[0:3]
(nt,) = input.readline().split()[0:1]
(dx,) = input.readline().split()[0:1]
dx = float(dx)/1000

# Retrieve the number of model parameters:
param_file = open("../Evento/fd3d_param", 'r')
i = 0
for row in param_file:
    cols = row.split()
    if i == 0:
        num_params = int(cols[0])
        i = i+1

# Load models tested in the inversion:
file = "../Evento/models.dat"
print("File =", file)
models = np.loadtxt(file)
param_size = np.size(models)
num_models = int(param_size/num_params)
itera = np.arange(num_models)
print("There are %i models tested in total" % (num_models))
models.shape = num_models,num_params

# Loading misfit values:
misfits = np.loadtxt("../Evento/misfits.dat")
min_misfit = misfits.min()
index_min = misfits.argmin()
model_best = models[index_min,:]

print("First model =", models[0,:])
print("Best model =", models[index_min,:])
print("Last model =", models[num_models-1,:])
print("Minimum misfit = %0.3f" % (min_misfit))
print("Index of minimum misfit =", index_min)

# Calculate Kappa for each of the models:
mu = 7.768456e10   # Shear modulus or rigidity of rocks (N/m2)
kappa = ["" for j in range(0,len(models[:,0]))]
itera = ["" for j in range(0,len(models[:,0]))]
for j in range(0,len(models[:,0])):
    kappa[j] = (((models[j,5]*1e6)**2)*(np.mean([models[j,0],models[j,1]])*dx*1000))/(mu*(models[j,5]*models[j,6]*1e6)*models[j,9]*2)
    itera[j] = j

# Calculate Kappa of the best model:   
kappa_best = (((model_best[5]*1e6)**2)*(np.mean([model_best[0],model_best[1]])*dx*1000))/(mu*(model_best[5]*model_best[6]*1e6)*model_best[9]*2)
    
    
    
# ======================================================================
#  FIGURES
# ======================================================================
# Plot length of ellipse axes and ellipse center coordinates:
plt.figure(1, figsize=(9,5))
plt.subplot(2,3,1)
plt.scatter(models[:,0]*dx, models[:,1]*dx, c=misfits, vmin=0, vmax=1, edgecolor='none', marker='o', s=8)
plt.scatter([model_best[0]*dx], [model_best[1]*dx], s=200, edgecolor='k', facecolor='w', lw=1, marker='*')
plt.xlabel("a (km)")
plt.ylabel("b (km)")

plt.subplot(2,3,2)
plt.scatter(models[:,2]*dx, models[:,3]*dx, c=misfits, vmin=0, vmax=1, edgecolor='none', marker='o', s=8)
plt.scatter([model_best[2]*dx], [model_best[3]*dx], s=200, edgecolor='k', facecolor='w', lw=1, marker='*')
plt.xlabel("x0 (km)")
plt.ylabel("y0 (km)")

plt.subplot(2,3,3)
plt.scatter(models[:,5], models[:,5]*models[:,6]*models[:,7], c=misfits, vmin=0, vmax=1, edgecolor='none', marker='o', s=8)
plt.scatter([model_best[5]], [model_best[5]*model_best[6]*model_best[7]], s=200, edgecolor='k', facecolor='w', lw=1, marker='*')
plt.xlabel("Te (MPa)")
plt.ylabel("Tu' (MPa)")

plt.subplot(2,3,4)
plt.scatter(models[:,5]*models[:,6], models[:,9]*2, c=misfits, vmin=0, vmax=1, edgecolor='none', marker='o', s=8)
plt.scatter([model_best[5]*model_best[6]], [model_best[9]*2], s=200, edgecolor='k', facecolor='w', lw=1, marker='*')
plt.xlabel("Tu (MPa)")
plt.ylabel("Dc (m)")

plt.subplot(2,3,5)
plt.scatter(models[:,8]*dx, models[:,5]*models[:,6]*models[:,7], c=misfits, vmin=0, vmax=1, edgecolor='none', marker='o', s=8)
plt.scatter([model_best[8]*dx], [model_best[5]*model_best[6]*model_best[7]], s=200, edgecolor='k', facecolor='w', lw=1, marker='*')
plt.xlabel("Nucleation\nradius (km)")
plt.ylabel("Tu' (MPa)")
#cbar = plt.colorbar()
#cbar.ax.set_title("Misfit", fontsize=10)

plt.subplots_adjust(left=0.08,right=0.97,bottom=0.15,top=0.97,wspace=0.43,hspace=0.41)
plt.savefig("../Figures/2D_marginals.jpg", dpi=300)
plt.show()
