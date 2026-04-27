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
units_param = 1



# ======================================================================
#  SEISMOGRAM DATA
# ======================================================================
# Retrieve time increment and number of points of the seismograms:
f = open("../src-covm_inkm/kine_param.inc", 'r')
for row in f:
    cols = row.split(',')
    if "parameter" in cols[0] and "dt=" in cols[0]:
        subcols1 = cols[0].split('=')
        dt = float(subcols1[1])
        subcols2 = cols[4].rstrip().split('=')
        subcols2 = subcols2[1].split(')')
        npts = int(subcols2[0])

# Retrieve the number of stations:
f = open("../Stations/stationn", 'r')
nsta = len(f.readlines())

# Loading observed and modeled data:
real_x = 100*(np.loadtxt("../Event/kine_files/real_disp_x"))
sismog_x = 100*(np.loadtxt("../Event/kine_files/best_disp_x"))
real_y = 100*(np.loadtxt("../Event/kine_files/real_disp_y"))
sismog_y = 100*(np.loadtxt("../Event/kine_files/best_disp_y"))
real_z = 100*(np.loadtxt("../Event/kine_files/real_disp_z"))
sismog_z = 100*(np.loadtxt("../Event/kine_files/best_disp_z"))



# ======================================================================
#  PLOT PARAMETERS
# ======================================================================
# Set time array, titles and names of the stations:
time = np.arange(npts)*dt
f = open("../Stations/stationn")
stations = ["" for i in range(0,nsta+1)]
j = 0
for line in f.readlines():
    line.strip()
    stations[j] = line[:-1]
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
plt.figure(figsize=(9,10))
for j in range(0,nsta):
    # Set common vertical limits for the current station:
    mobs_x = max(abs(real_x[npts*j:npts*(j+1)]))
    mmod_x = max(abs(sismog_x[npts*j:npts*(j+1)]))
    mobs_y = max(abs(real_y[npts*j:npts*(j+1)]))
    mmod_y = max(abs(sismog_y[npts*j:npts*(j+1)]))
    mobs_z = max(abs(real_z[npts*j:npts*(j+1)]))
    mmod_z = max(abs(sismog_z[npts*j:npts*(j+1)]))
    maxamps = [mobs_x, mmod_x, mobs_y, mmod_y, mobs_z, mmod_z]
    maxv = max(maxamps)
    limit = round_to_1(maxv*1.5)
   
   
    # North component:
    plt.subplot(nsta,3,3*j+1)	
    plt.plot(time,real_x[npts*j:npts*(j+1)],'b')
    plt.plot(time,sismog_x[npts*j:npts*(j+1)],'r')
    plt.xlim(0,npts*dt)
    plt.ylim(-limit,limit)
    plt.yticks([-limit, -limit/2, 0, limit/2, limit])
    plt.text(0.,limit/4.0," %s" % (stations[j]))


    # East component:
    plt.subplot(nsta,3,3*j+2)	
    plt.plot(time,real_y[npts*j:npts*(j+1)],'b')
    plt.plot(time,sismog_y[npts*j:npts*(j+1)],'r')
    plt.xlim(0,npts*dt)
    plt.ylim(-limit,limit)
    plt.yticks([-limit, -limit/2, 0, limit/2, limit])
    plt.text(0.,limit/4.0," %s" % (stations[j]))


    # Vertical component:
    plt.subplot(nsta,3,3*j+3)	
    plt.plot(time,real_z[npts*j:npts*(j+1)],'b')
    plt.plot(time,sismog_z[npts*j:npts*(j+1)],'r')
    plt.xlim(0,npts*dt)
    plt.ylim(-limit,limit)
    plt.yticks([-limit, -limit/2, 0, limit/2, limit])
    plt.text(0.,limit/4.0," %s" % (stations[j]))

# Add titles, labels, set subplots, and save figure:
t = plt.gcf().text(0.5, 0.98, "All Seismograms",horizontalalignment='center')
t = plt.gcf().text(0.22, 0.95, "N-S",horizontalalignment='center')
t = plt.gcf().text(0.55, 0.95, "E-W",horizontalalignment='center')
t = plt.gcf().text(0.88, 0.95, "Z",horizontalalignment='center')
t = plt.gcf().text(0.5, 0.02, "Time [s]",horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units,verticalalignment='center', rotation='vertical')
plt.subplots_adjust(wspace=0.53,hspace=0.55,left=0.12,right=0.98,bottom=0.06,top=0.94)
plt.savefig("../Figures/Seismog_all.png", dpi=300)



plt.show()
