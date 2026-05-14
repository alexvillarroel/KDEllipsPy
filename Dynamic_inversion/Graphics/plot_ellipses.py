#!/usr/bin/env python3

# Main Python imports:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse



########## PARAMETERS ##########
mis0 = 0.3   # Maximum misfit for model plotting



########## DATA LOAD AND PROCESSING ##########
# Retrieve the number of model parameters:
param_file = open("../Evento/fd3d_param", 'r')
i = 0
for row in param_file:
    cols = row.split()
    if i == 0:
        num_params = int(cols[0])
        i = i+1

# Read output data from the inversion:
models = np.loadtxt("../Evento/models.dat")
misfits = np.loadtxt("../Evento/misfits.dat")

# Reshape the data from the "models" file:
all_params = len(models)
num_models = int(all_params/num_params)
models.shape = (num_models,num_params)

# Retrieve parameters of selected models:
a = []
b = []
x0 = []
y0 = []
phi = []
for i in range(len(misfits)):
   # Select only models with misfits equal or less than a certain threshold:
   if misfits[i] <= mis0:
      a.append(models[i,0]*0.2)
      b.append(models[i,1]*0.2)
      x0.append(models[i,2]*0.2)
      y0.append(models[i,3]*0.2)
      phi.append(np.degrees(models[i,4]))

# Find the model with minimum misfit (best model) and convert units to "km"
# and degrees when appropriate:
min_misfit = misfits.min()
index_min = misfits.argmin()
a_best = models[index_min,0]*0.2
b_best = models[index_min,1]*0.2
x0_best = models[index_min,2]*0.2
y0_best = models[index_min,3]*0.2
phi_best = np.degrees(models[index_min,4])

# Find the average model within the selection and convert units to "km"
# and degrees when appropriate:
a_ave = np.mean(a)
b_ave = np.mean(b)
x0_ave = np.mean(x0)
y0_ave = np.mean(y0)
phi_ave = np.mean(phi)



########## CREATE FIGURE ##########
# Create ellipses of every selected model:
ells = [Ellipse((x0[i],y0[i]),2*a[i],2*b[i],phi[i], fill=None, edgecolor=(0.5,0.5,0.5)) for i in range(len(x0))]

# Ellipse of the best model:
ebest = Ellipse((x0_best,y0_best),2*a_best,2*b_best,phi_best, fill=None, edgecolor=(0,0,1), linewidth=2)

# Ellipse of average model:
eave = Ellipse((x0_ave,y0_ave),2*a_ave,2*b_ave,phi_ave, fill=None, edgecolor=(1,0,0), linewidth=2)

# Plot ellipses based on misfit threshold:
fig, ax = plt.subplots()
for e in ells:
    ax.add_artist(e)
ax.add_artist(ebest)
ax.add_artist(eave)

# Add rupture plane center and misfit legend:
plt.plot(80*0.2,80*0.2,'*r', markersize=12)
plt.text(110*0.2,10*0.2,"Misfit $\leq$ %0.2f" % (mis0), fontsize=11, bbox=dict(facecolor='w', alpha=0.8, edgecolor='k'))
plt.xlim(0,160*0.2)
plt.ylim(0,160*0.2)
plt.title("Earthquake source models", fontsize=12)
plt.xlabel("Along strike (km)", fontsize=11)
plt.ylabel("Alond dip (km)", fontsize=11)
plt.xticks([0, 80*0.2, 160*0.2], ('0', '16', '32'))
plt.yticks([0, 80*0.2, 160*0.2], ('0', '16', '32'))
ax.xaxis.set_tick_params(labelsize=11)
ax.yaxis.set_tick_params(labelsize=11)
plt.gca().invert_yaxis()
plt.gca().set_aspect('equal')
plt.savefig("../Figures/Ellipses.png", dpi=300)



plt.show()
