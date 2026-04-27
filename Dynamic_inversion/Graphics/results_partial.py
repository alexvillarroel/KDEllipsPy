#!/usr/bin/env python3

# Main Python imports:
import math, glob
import numpy as np
import matplotlib.pyplot as plt



###################################################################
#   SCRIPT TO PLOT THE PARTIAL RESULTS OF THE INVERSION PROCESS   #
###################################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
mu = 7.768456e10   # Shear modulus or rigidity of rocks (N/m2)



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

# Load misfit values:
misfits = np.loadtxt("../Evento/misfits.dat")
num_misfit = np.size(misfits)
min_misfit = misfits.min()
index_min = misfits.argmin()

# Load models tested in the inversion:
file = "../Evento/models.dat"
print("File =", file)
models = np.loadtxt(file)
print("There are", num_misfit, "models tested in total")
print("The minimum misfit is", "%0.3f" % (min_misfit), "at the model number", index_min)
models.shape = num_misfit,num_params

# Print the best model and its parameters:
best_par = models[index_min,:]
print("Best model is:")
print(best_par)
print("Axis 1 = %0.2f km" % (best_par[0]*dx))
print("Axis 2 = %0.2f km" % (best_par[1]*dx))
print("Te = %0.2f MPa" % (best_par[5]))
print("Tu = %0.2f MPa" % (best_par[5]*best_par[6]))
print("Dc = %0.2f m" % (best_par[9]*2))

# Export the best model into a file:
out_model = open("../Evento/model_best.dat",'w')
out_model.write("%i\n" % (num_params))
for j in range(num_params):
	out_model.write("%f\n" % float(models[index_min,j]))
out_model.close()



# ======================================================================
#  SIMILARITY PARAMETER (KAPPA)
# ======================================================================
# Calculate Kappa of the best model:
kappa1 = (((best_par[5]*1e6)**2)*(best_par[0]*dx*1000))/(mu*(best_par[5]*best_par[6]*1e6)*best_par[9]*2)
kappa2 = (((best_par[5]*1e6)**2)*(best_par[1]*dx*1000))/(mu*(best_par[5]*best_par[6]*1e6)*best_par[9]*2)
print("Kappa Axis 1 =", "%0.2f" % (kappa1))
print("Kappa Axis 2 =", "%0.2f" % (kappa2))

# Calculate Kappa for each one of the models:
kappa = ["" for j in range(0,len(models))]
itera = ["" for j in range(0,len(models))]
for j in range(0,len(models)):
   param = models[j,:]
   # Force to use the Kappa of the shorter axis:
   if param[0] < param[1]:
      kappa[j] = (((param[5]*1e6)**2)*(param[0]*dx*1000))/(mu*(param[5]*param[6]*1e6)*param[9]*2)
   else:
      kappa[j] = (((param[5]*1e6)**2)*(param[1]*dx*1000))/(mu*(param[5]*param[6]*1e6)*param[9]*2)
   itera[j] = j



# ======================================================================
#  MOMENT MAGNITUDE
# ======================================================================
# Identify the number of existent momentos.dat files and load their data:
m_files = glob.glob("../Evento/momentos*.dat")
all_m = {}
for i in range(len(m_files)):
   all_m["momentos_file{0}".format(i)] = np.loadtxt("../Evento/momentos%i.dat" % (i))

# Retrieve the slip values from the momentos.dat files:
slip = np.zeros(num_misfit)
k = 0
for i in range(int(num_misfit/len(m_files))):
   for j in range(len(m_files)):
      momentos_file = all_m["momentos_file%d" % (j)]
      slip[k] = momentos_file[i,2]
      k = k+1

# Calculate seismic moment (Mo) and moment magnitude (Mw) for each of 
# the models:
Mo = mu*slip*((dx*1000)*(dx*1000))
Mw = np.zeros(len(Mo))
for i in range(len(Mo)):
    if Mo[i] != 0:
        Mw[i] = (math.log10(Mo[i])-9.1)/1.5



# ======================================================================
#  FIGURES
# ======================================================================
plt.figure(1)
plt.scatter(itera,misfits,c=kappa,vmin=0,vmax=3,marker='o',edgecolor='none',cmap='jet')
plt.title("Misfit convergence")
plt.xlabel("Number of models")
plt.ylim(0,5)
plt.xlim(0,len(itera)+1000)
plt.ylabel("Misfit")
cbar = plt.colorbar()
cbar.ax.set_title('$\kappa$')
plt.savefig("../Figures/Partial_results_Kappa.png", dpi=300)

plt.figure(2)
Mw_best = round(Mw[index_min],1)
Mw_ticks = np.linspace(Mw_best-0.5,Mw_best+0.5,11)
plt.scatter(itera,misfits,c=Mw,vmin=Mw_best-0.5,vmax=Mw_best+0.5,marker='o',edgecolor='none',cmap='jet')
plt.title("Misfit convergence")
plt.xlabel("Number of models")
plt.ylim(0,5)
plt.xlim(0,len(itera)+1000)
plt.ylabel("Misfit")
cbar = plt.colorbar(ticks=Mw_ticks)
cbar.ax.set_title('$M_{w}$')
plt.savefig("../Figures/Partial_results_Mw.png", dpi=300)

plt.figure(3)
plt.scatter(itera[0:34000],misfits[0:34000],c=models[0:34000,5],vmin=best_par[5]-2,vmax=best_par[5]+2,marker='o',edgecolor='none',cmap='jet')
plt.title("Misfit convergence")
plt.xlabel("Number of models")
plt.ylim(0,5)
plt.xlim(0,len(itera)+1000)
plt.ylabel("Misfit")
cbar = plt.colorbar()
cbar.ax.set_title('(MPa)')
plt.savefig("../Figures/Partial_results_Te.png", dpi=300)



plt.show()
