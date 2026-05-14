#!/usr/bin/env python3

# Main Python imports:
import math
import numpy as np
import matplotlib.pyplot as plt



##########################################################
#   SCRIPT TO PLOT OBSERVED AND MODELED ROTATED TRACES   #
#                     IN ONE FIGURE                      #
##########################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
units_param = 2
window = 20



# ======================================================================
#  SEISMOGRAM DATA
# ======================================================================
# Azimuth and window parameters:
rw_data = np.loadtxt("../Evento/azi_times.txt")
azi = rw_data[:,0]
tP = rw_data[:,1]
tS = rw_data[:,2]

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

# Loading observed and modeled data:
real_x = 100*(np.loadtxt("../Evento/real_disp_x"))
sismog_x = 100*(np.loadtxt("../Evento/dyn_disp_x"))
real_y = 100*(np.loadtxt("../Evento/real_disp_y"))
sismog_y = 100*(np.loadtxt("../Evento/dyn_disp_y"))
real_z = 100*(np.loadtxt("../Evento/real_disp_z"))
sismog_z = 100*(np.loadtxt("../Evento/dyn_disp_z"))



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
    # Calculate radial and transverse components:
    radial_obs = (real_x[npts*j:npts*(j+1)]*math.cos(azi[j]))+(real_y[npts*j:npts*(j+1)]*math.sin(azi[j]))
    radial_mod = (sismog_x[npts*j:npts*(j+1)]*math.cos(azi[j]))+(sismog_y[npts*j:npts*(j+1)]*math.sin(azi[j]))
    trans_obs = (real_y[npts*j:npts*(j+1)]*math.cos(azi[j]))-(real_x[npts*j:npts*(j+1)]*math.sin(azi[j]))
    trans_mod = (sismog_y[npts*j:npts*(j+1)]*math.cos(azi[j]))-(sismog_x[npts*j:npts*(j+1)]*math.sin(azi[j]))
    
    # Set common vertical limits for the current station:
    mobs_r = max(abs(radial_obs))
    mmod_r = max(abs(radial_mod))
    mobs_t = max(abs(trans_obs))
    mmod_t = max(abs(trans_mod))
    mobs_z = max(abs(real_z[npts*j:npts*(j+1)]))
    mmod_z = max(abs(sismog_z[npts*j:npts*(j+1)]))
    maxamps = [mobs_r, mmod_r, mobs_t, mmod_t, mobs_z, mmod_z]
    maxv = max(maxamps)
    limit = round_to_1(maxv*1.5)
    
    # Radial component:
    plt.subplot(nsta,3,3*j+1)
    plt.axvspan(tP[j], tP[j]+window, color='yellow', alpha=0.5)
    plt.plot(time,radial_obs,'b')
    plt.plot(time,radial_mod,'r')
    if 3*j+1 != 3*(nsta-1)+1:
        plt.gca().xaxis.set_ticklabels([])
    else:
        plt.xticks(fontsize=15)
    plt.gca().yaxis.set_ticks([])
    plt.xlim(0,npts*dt)
    plt.ylim(-limit,limit)
    plt.text(0.,limit/2.0," %s" % (stations[j]), fontsize=14, weight='bold')
    plt.text(0.,-limit/1.2," %0.3f" % (max(maxamps[0],maxamps[1])), fontsize=14)


    # Transverse component:
    plt.subplot(nsta,3,3*j+2)
    plt.axvspan(tS[j], tS[j]+window, color='yellow', alpha=0.5)
    plt.plot(time,trans_obs,'b')
    plt.plot(time,trans_mod,'r')
    if 3*j+2 != 3*(nsta-1)+2:
        plt.gca().xaxis.set_ticklabels([])
    else:
        plt.xticks(fontsize=15)
    plt.gca().yaxis.set_ticks([])
    plt.xlim(0,npts*dt)
    plt.ylim(-limit,limit)
    plt.text(0.,limit/2.0," %s" % (stations[j]), fontsize=14, weight='bold')
    plt.text(0.,-limit/1.2," %0.3f" % (max(maxamps[2],maxamps[3])), fontsize=14)


    # Vertical component:
    plt.subplot(nsta,3,3*j+3)
    plt.axvspan(tP[j], tP[j]+window, color='yellow', alpha=0.5)
    plt.plot(time,real_z[npts*j:npts*(j+1)],'b')
    plt.plot(time,sismog_z[npts*j:npts*(j+1)],'r')
    if 3*j+3 != 3*(nsta-1)+3:
        plt.gca().xaxis.set_ticklabels([])
    else:
        plt.xticks(fontsize=15)
    plt.gca().yaxis.set_ticks([])
    plt.xlim(0,npts*dt)
    plt.ylim(-limit,limit)
    plt.text(0.,limit/2.0," %s" % (stations[j]), fontsize=14, weight='bold')
    plt.text(0.,-limit/1.2," %0.3f" % (max(maxamps[4],maxamps[5])), fontsize=14)

# Add titles, labels, set subplots, and save figure:
t = plt.gcf().text(0.18, 0.97, "Radial", horizontalalignment='center', fontsize=15)
t = plt.gcf().text(0.50, 0.97, "Transverse", horizontalalignment='center', fontsize=15)
t = plt.gcf().text(0.85, 0.97, "Vertical", horizontalalignment='center', fontsize=15)
t = plt.gcf().text(0.5, 0.01, "Time (s)", horizontalalignment='center', fontsize=15)
plt.subplots_adjust(wspace=0.2,hspace=0.17,left=0.04,right=0.98,bottom=0.07,top=0.96)
plt.savefig("../Figures/Seismog_all_rot.pdf")



plt.show()
