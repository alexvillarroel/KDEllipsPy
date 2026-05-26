#!/usr/bin/env python3

# Main Python imports:
import math
import numpy as np
from scipy.stats import norm, sem
import matplotlib.pyplot as plt



# ======================================================================
# USER PARAMETERS
# ======================================================================
# Histograms are plotted starting at this model number:
startmod = 10000

# Number of standard deviations to be plotted along horizontal axis and number 
# of bins for histograms:
nstd = 5
nbin = 30

# Shear modulus or rigidity of rocks (N/m2)
mu = 7.768456e10



# ======================================================================
# DATA
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
print("There are %i model parameters" % (num_params))

# Load models tested in the inversion:
models = np.loadtxt("../Evento/models.dat")
param_size = np.size(models)
num_models = int(param_size/num_params)
print("There are %i models tested in total" % (num_models))
models.shape = num_models,num_params

# Load misfit values and find the best model:
misfits = np.loadtxt("../Evento/misfits.dat")
min_misfit = misfits.min()
index_min = misfits.argmin()
model_best = models[index_min,:]

# Trim range of models:
models0 = models[startmod:,:]

# Calculate Kappa for each of the models:
kappa0 = ["" for j in range(0,len(models0[:,0]))]
itera0 = ["" for j in range(0,len(models0[:,0]))]
for j in range(0,len(models0[:,0])):
    kappa0[j] = (((models0[j,5]*1e6)**2)*(np.mean([models0[j,0],models0[j,1]])*dx*1000))/(mu*(models0[j,5]*models0[j,6]*1e6)*models0[j,9]*2)
    itera0[j] = j

# Calculate Kappa of the best model:   
kappa_best = (((model_best[5]*1e6)**2)*(np.mean([model_best[0],model_best[1]])*dx*1000))/(mu*(model_best[5]*model_best[6]*1e6)*model_best[9]*2)



# ======================================================================
#  FIGURES
# ======================================================================
# Create plot:
plt.figure(1, figsize=(8,7))

plt.subplot(3,4,1)
# Gaussian distribution:
ave0 = np.mean(models0[:,0]*dx)
std0 = np.std(models0[:,0]*dx, ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)

#h = sem(models0[:,0]*dx)*norm.ppf((1+0.95)/2, len(models0[:,0])-1)
#print(ave0-h, ave0, ave0+h)
#print(norm.interval(0.2, loc=ave0, scale=sem(models[:,0]*dx)))
#print(1.645*(std0/np.sqrt(len(models[:,0]))))

text_str = '\n'.join(("$S_{0}$ = %0.3f" % (model_best[0]*dx), "$\mu$ = %0.3f" % (ave0), "$\sigma$ = %0.3f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,0]*dx, bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[0]*dx, color='r', linewidth=2)
plt.xlabel("a (km)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,2)
# Gaussian distribution:
ave0 = np.mean(models0[:,1]*dx)
std0 = np.std(models0[:,1]*dx, ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.3f" % (model_best[1]*dx), "$\mu$ = %0.3f" % (ave0), "$\sigma$ = %0.3f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,1]*dx, bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[1]*dx, color='r', linewidth=2)
plt.xlabel("b (km)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,3)
# Gaussian distribution:
ave0 = np.mean(models0[:,2]*dx)
std0 = np.std(models0[:,2]*dx, ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.3f" % (model_best[2]*dx), "$\mu$ = %0.3f" % (ave0), "$\sigma$ = %0.3f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,2]*dx, bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[2]*dx, color='r', linewidth=2)
plt.xlabel("x0 (km)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,4)
# Gaussian distribution:
ave0 = np.mean(models0[:,3]*dx)
std0 = np.std(models0[:,3]*dx, ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.3f" % (model_best[3]*dx), "$\mu$ = %0.3f" % (ave0), "$\sigma$ = %0.3f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,3]*dx, bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[3]*dx, color='r', linewidth=2)
plt.xlabel("y0 (km)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,5)
# Gaussian distribution:
ave0 = np.mean(models0[:,4]*(180/math.pi))
std0 = np.std(models0[:,4]*(180/math.pi), ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.1f" % (model_best[4]*(180/math.pi)), "$\mu$ = %0.1f" % (ave0), "$\sigma$ = %0.1f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,4]*(180/math.pi), bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[4]*(180/math.pi), color='r', linewidth=2)
plt.xlabel("Rotation\nangle (deg)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,6)
# Gaussian distribution:
ave0 = np.mean(models0[:,5])
std0 = np.std(models0[:,5], ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.2f" % (model_best[5]), "$\mu$ = %0.2f" % (ave0), "$\sigma$ = %0.2f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,5], bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[5], color='r', linewidth=2)
plt.xlabel("Te (MPa)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,7)
# Gaussian distribution:
ave0 = np.mean(models0[:,5]*models0[:,6])
std0 = np.std(models0[:,5]*models0[:,6], ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.2f" % (model_best[5]*model_best[6]), "$\mu$ = %0.2f" % (ave0), "$\sigma$ = %0.2f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,5]*models0[:,6], bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[5]*model_best[6], color='r', linewidth=2)
plt.xlabel("Tu (MPa)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,8)
# Gaussian distribution:
ave0 = np.mean(models0[:,5]*models0[:,6]*models0[:,7])
std0 = np.std(models0[:,5]*models0[:,6]*models0[:,7], ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.2f" % (model_best[5]*model_best[6]*model_best[7]), "$\mu$ = %0.2f" % (ave0), "$\sigma$ = %0.2f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,5]*models0[:,6]*models0[:,7], bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[5]*model_best[6]*model_best[7], color='r', linewidth=2)
plt.xlabel("Tu' (MPa)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,9)
# Gaussian distribution:
ave0 = np.mean(models0[:,8]*dx)
std0 = np.std(models0[:,8]*dx, ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.3f" % (model_best[8]*dx), "$\mu$ = %0.3f" % (ave0), "$\sigma$ = %0.3f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,8]*dx, bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[8]*dx, color='r', linewidth=2)
plt.xlabel("Nucleation\nradius (km)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,10)
# Gaussian distribution:
ave0 = np.mean(models0[:,9]*2)
std0 = np.std(models0[:,9]*2, ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.3f" % (model_best[9]*2), "$\mu$ = %0.3f" % (ave0), "$\sigma$ = %0.3f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(models0[:,9]*2, bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=model_best[9]*2, color='r', linewidth=2)
plt.xlabel("Dc (m)")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplot(3,4,11)
# Gaussian distribution:
ave0 = np.mean(kappa0)
std0 = np.std(kappa0, ddof=1)
x = np.linspace(ave0-std0*nstd,ave0+std0*nstd,50)
pdf0 = norm.pdf(x,ave0,std0)
text_str = '\n'.join(("$S_{0}$ = %0.3f" % (kappa_best), "$\mu$ = %0.3f" % (ave0), "$\sigma$ = %0.3f" % (std0)))
plt.text(0.99, 0.93, text_str, fontsize=8, horizontalalignment='right', verticalalignment='top', transform = plt.gca().transAxes, bbox=dict(facecolor='w', alpha=0.0))
bins = np.linspace(ave0-std0*nstd,ave0+std0*nstd,nbin)
plt.hist(kappa0, bins, density=True, facecolor='#4242ff')
plt.plot(x,pdf0,'c',linewidth=2)
plt.axvline(x=kappa_best, color='r', linewidth=2)
plt.xlabel("$\kappa$")
plt.xlim(ave0-std0*(nstd+2),ave0+std0*(nstd+2))
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.subplots_adjust(left=0.03,right=0.97,bottom=0.1,top=0.97,hspace=0.5)
plt.savefig("../Figures/param_hist.jpg", dpi=300)
plt.show()
