#!/usr/bin/env python3

# Main Python imports:
import math
import numpy as np
import matplotlib.pyplot as plt



##################################################
#   SCRIPT TO PLOT OBSERVED AND MODELED TRACES   #
##################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
units_param=2



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
time = [i*dt for i in range(0,npts)]
f = open("../Stations/stationn", 'r')
stations = []
for row in f:
    row1 = row.strip()
    stations.append(row1)

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

# Set the subplots according to the number of stations:
if (nsta % 2 == 0):
   # Even number of stations:
   nrow_plot = nsta/2
else:
   # Odd number of stations:
   nrow_plot = (nsta/2)+0.5

# Function to control vertical limits:
def round_to_1(x):
   return round(x, -int(math.floor(math.log10(x))))



# ======================================================================
#  FIGURES
# ======================================================================
# Plot North components and absolute minimums and maximums:
print("Observed North = [%0.3f, %0.3f] %s" % (min(real_x), max(real_x), units))
plt.figure(1, figsize=(7,9))
for j in range(0,nsta):
   plt.subplot(nrow_plot,2,j+1)	
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
   plt.yticks([-limit, -limit/2, 0, limit/2, limit])
   plt.text(0.,limit/4.0," %s" % (stations[j]))
# Add title, labels, set subplots and save figure:
t = plt.gcf().text(0.5, 0.98, "North component",horizontalalignment='center')
t = plt.gcf().text(0.5, 0.01, "Time [s]",horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units,verticalalignment='center',rotation='vertical')
plt.subplots_adjust(wspace=0.35,hspace=0.50,bottom=0.06,top=0.95,left=0.15,right=0.96)
plt.savefig("../Figures/Seismog_north.png", dpi=300)


# Plot East components and absolute minimums and maximums:
print("Observed East = [%0.3f, %0.3f] %s" % (min(real_y), max(real_y), units))
plt.figure(2, figsize=(7,9))
for j in range(0,nsta):
   plt.subplot(nrow_plot,2,j+1)
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
   plt.yticks([-limit, -limit/2, 0, limit/2, limit])
   plt.text(0.,limit/4.0," %s" % (stations[j]))
# Add title, labels, set subplots and save figure:
t = plt.gcf().text(0.5, 0.98, "East component",horizontalalignment='center')
t = plt.gcf().text(0.5, 0.01, "Time [s]",horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units,verticalalignment='center',rotation='vertical')
plt.subplots_adjust(wspace=0.35,hspace=0.50,bottom=0.06,top=0.95,left=0.15,right=0.96)
plt.savefig("../Figures/Seismog_east.png", dpi=300)


# Plot Vertical components and absolute minimums and maximums:
print("Observed Vertical = [%0.3f, %0.3f] %s" % (min(real_z), max(real_z), units))
plt.figure(3, figsize=(7,9))
for j in range(0,nsta):
   plt.subplot(nrow_plot,2,j+1)
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
   plt.yticks([-limit, -limit/2, 0, limit/2, limit])
   plt.text(0.,limit/4.0," %s" % (stations[j]))
# Add title, labels, set subplots and save figure:
t = plt.gcf().text(0.5, 0.98, "Vertical component",horizontalalignment='center')
t = plt.gcf().text(0.5, 0.01, "Time [s]",horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units,verticalalignment='center',rotation='vertical')
plt.subplots_adjust(wspace=0.35,hspace=0.50,bottom=0.06,top=0.95,left=0.15,right=0.96)
plt.savefig("../Figures/Seismog_vertical.png", dpi=300)



plt.show()
