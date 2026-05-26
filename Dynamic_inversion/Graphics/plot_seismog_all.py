#!/usr/bin/env python3

# Main Python imports:
import math
import numpy as np
import matplotlib.pyplot as plt



##################################################
#   SCRIPT TO PLOT OBSERVED AND MODELED TRACES   #
#                 IN ONE FIGURE                  #
##################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
units_param=1



# ======================================================================
#  SEISMOGRAM DATA
# ======================================================================
# Retrieve time increment and number of points of the seismograms:
f = open("../Source/fd3d_subs/dimension.inc", 'r')
for row in f:
    cols = row.split(',')
    if "parameter" in cols[0]:
        subcols1 = cols[0].split('=')
        dt = float(subcols1[1])
        subcols2 = cols[4].split('=')
        npts = int(subcols2[1])
        
# Retrieve the number of stations:
f = open("../Stations/stationn", 'r')
nsta = len(f.readlines())

# Load observed and modeled data (original data is in (m/s) or (m) for
# velocity or displacement, respectively):
real_x = 100*(np.loadtxt("../Evento/real_disp_x"))
sismog_x = 100*(np.loadtxt("../Evento/dyn_disp_x"))
real_y = 100*(np.loadtxt("../Evento/real_disp_y"))
sismog_y = 100*(np.loadtxt("../Evento/dyn_disp_y"))
real_z = 100*(np.loadtxt("../Evento/real_disp_z"))
sismog_z = 100*(np.loadtxt("../Evento/dyn_disp_z"))



# ======================================================================
#  PLOT PARAMETERS
# ======================================================================
# Set time array and station names:
time = np.arange(npts)*dt
fd_titulos = open("../Stations/stationn")
titulos = ["" for i in range(0,nsta+1)]
j = 0
for line in fd_titulos.readlines():
   line.strip()
   titulos[j] = line[:-1]
   j+=1

# Physical units of seismograms:
if (units_param == 1):
   label_units = "Displacement (cm)"
   units = "(cm)"
elif (units_param == 2):
   label_units = "Velocity (cm/s)"
   units = "(cm/s)"
else:
   label_units = ""
   units = ""

# Function to control vertical limits:
def round_to_1(x):
   return round(x, -int(math.floor(math.log10(x))))



# ======================================================================
#  FIGURES
# ======================================================================
# Plot all components and stations:
plt.figure(figsize=(7,9))
for j in range(0,nsta):
   # North component:
   plt.subplot(nsta,3,3*j+1)	
   plt.plot(time,real_x[npts*j:npts*(j+1)],'b')
   plt.plot(time,sismog_x[npts*j:npts*(j+1)],'r')
   plt.xlim(0,npts*dt)
   # Set vertical limits:
   max_obs = max(abs(real_x[npts*j:npts*(j+1)]))
   max_mod = max(abs(sismog_x[npts*j:npts*(j+1)]))
   if max_obs > max_mod:
      maxv = max_obs
   elif max_obs < max_mod:
      maxv = max_mod
   limit = round_to_1(maxv*1.5)
   plt.ylim(-limit,limit)
   plt.yticks([-limit, 0, limit])
   if 3*j+1 != 3*nsta-2:
      plt.tick_params(axis='x', labelbottom=False)
   plt.text(0.,limit/4.0," %s" % (titulos[j]))


   # East component:
   plt.subplot(nsta,3,3*j+2)
   plt.plot(time,real_y[npts*j:npts*(j+1)],'b')
   plt.plot(time,sismog_y[npts*j:npts*(j+1)],'r')
   plt.xlim(0,npts*dt)
   # Set vertical limits:
   max_obs = max(abs(real_y[npts*j:npts*(j+1)]))
   max_mod = max(abs(sismog_y[npts*j:npts*(j+1)]))
   if max_obs > max_mod:
      maxv = max_obs
   elif max_obs < max_mod:
      maxv = max_mod
   limit = round_to_1(maxv*1.5)
   plt.ylim(-limit,limit)
   plt.yticks([-limit, 0, limit])
   if 3*j+2 != 3*nsta-1:
      plt.tick_params(axis='x', labelbottom=False)
   plt.text(0.,limit/4.0," %s" % (titulos[j]))


   # Vertical component:
   plt.subplot(nsta,3,3*j+3)
   plt.plot(time,real_z[npts*j:npts*(j+1)],'b')
   plt.plot(time,sismog_z[npts*j:npts*(j+1)],'r')
   plt.xlim(0,npts*dt)
   # Set vertical limits:
   max_obs = max(abs(real_z[npts*j:npts*(j+1)]))
   max_mod = max(abs(sismog_z[npts*j:npts*(j+1)]))
   if max_obs > max_mod:
      maxv = max_obs
   elif max_obs < max_mod:
      maxv = max_mod
   limit = round_to_1(maxv*1.5)
   plt.ylim(-limit,limit)
   plt.yticks([-limit, 0, limit])
   if 3*j+3 != 3*nsta:
      plt.tick_params(axis='x', labelbottom=False)
   plt.text(0.,limit/4.0," %s" % (titulos[j]))
# Add title, labels, set subplots and save figure:
t = plt.gcf().text(0.22, 0.98, "N-S",horizontalalignment='center')
t = plt.gcf().text(0.55, 0.98, "E-W",horizontalalignment='center')
t = plt.gcf().text(0.88, 0.98, "Z",horizontalalignment='center')
t = plt.gcf().text(0.5, 0.02, "Time (s)",horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units,verticalalignment='center', rotation='vertical')
plt.subplots_adjust(wspace=0.55,hspace=0.55,left=0.12,right=0.98,bottom=0.06,top=0.97)
plt.savefig("../Figures/Seismog_all.pdf")



plt.show()
