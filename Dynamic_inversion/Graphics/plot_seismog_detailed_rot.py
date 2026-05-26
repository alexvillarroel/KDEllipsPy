#!/usr/bin/env python3

# Main Python imports:
import math
import numpy as np
import matplotlib.pyplot as plt



##########################################################
#   SCRIPT TO PLOT OBSERVED AND MODELED ROTATED TRACES   #
##########################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
units_param=2
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
    nrow_plot = int(nsta/2)
else: 
    # Odd number of stations:
    nrow_plot = int((nsta/2)+0.5)

# Function to control vertical limits:
def round_to_1(x):
    return round(x, -int(math.floor(math.log10(x))))



# ======================================================================
#  FIGURES
# ======================================================================
# Plot Radial components:
plt.figure(1, figsize=(7,9))
for j in range(0,nsta):
    plt.subplot(nrow_plot,2,j+1)
    # Calculate radial components:
    radial_obs = (real_x[npts*j:npts*(j+1)]*math.cos(azi[j]))+(real_y[npts*j:npts*(j+1)]*math.sin(azi[j]))
    radial_mod = (sismog_x[npts*j:npts*(j+1)]*math.cos(azi[j]))+(sismog_y[npts*j:npts*(j+1)]*math.sin(azi[j]))
    # Plot waveforms:
    plt.axvspan(tP[j], tP[j]+window, color='yellow', alpha=0.5)
    plt.plot(time,radial_obs,'b')
    plt.plot(time,radial_mod,'r')
    plt.xlim(0,npts*dt)
    # Set vertical limits:
    max_obs = max(abs(radial_obs))
    max_mod = max(abs(radial_mod))
    if max_obs > max_mod:
        maxv = max_obs
    elif max_obs < max_mod:
        maxv = max_mod
    limit = round_to_1(maxv*1.5)
    plt.ylim(-limit,limit)
    plt.yticks([-limit, -limit/2, 0, limit/2, limit])
    plt.text(0.,limit/4.0," %s" % (stations[j]))
# Add title, labels, set subplots, and save figure:
t = plt.gcf().text(0.5, 0.98, "Radial component", horizontalalignment='center')
t = plt.gcf().text(0.5, 0.01, "Time (s)", horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units, verticalalignment='center', rotation='vertical')
plt.subplots_adjust(wspace=0.35,hspace=0.50,bottom=0.06,top=0.95,left=0.15,right=0.96)
plt.savefig("../Figures/Seismog_radial.png", dpi=300)


# Plot Transverse components:
plt.figure(2, figsize=(7,9))
for j in range(0,nsta):
    plt.subplot(nrow_plot,2,j+1)
    # Calculate transverse components:
    trans_obs = (real_y[npts*j:npts*(j+1)]*math.cos(azi[j]))-(real_x[npts*j:npts*(j+1)]*math.sin(azi[j]))
    trans_mod = (sismog_y[npts*j:npts*(j+1)]*math.cos(azi[j]))-(sismog_x[npts*j:npts*(j+1)]*math.sin(azi[j]))
    # Plot waveforms:
    plt.axvspan(tS[j], tS[j]+window, color='yellow', alpha=0.5)
    plt.plot(time,trans_obs,'b')
    plt.plot(time,trans_mod,'r')
    plt.xlim(0,npts*dt)
    # Set vertical limits:
    max_obs = max(abs(trans_obs))
    max_mod = max(abs(trans_mod))
    if max_obs > max_mod:
        maxv = max_obs
    elif max_obs < max_mod:
        maxv = max_mod
    limit = round_to_1(maxv*1.5)
    plt.ylim(-limit,limit)
    plt.yticks([-limit, -limit/2, 0, limit/2, limit])
    plt.text(0.,limit/4.0," %s" % (stations[j]))
# Add title, labels, set subplots, and save figure:
t = plt.gcf().text(0.5, 0.98, "Transverse component", horizontalalignment='center')
t = plt.gcf().text(0.5, 0.01, "Time (s)", horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units, verticalalignment='center', rotation='vertical')
plt.subplots_adjust(wspace=0.35,hspace=0.50,bottom=0.06,top=0.95,left=0.15,right=0.96)
plt.savefig("../Figures/Seismog_transverse.png", dpi=300)
	       

# Plot Vertical components:
plt.figure(3, figsize=(7,9))
for j in range(0,nsta):
    plt.subplot(nrow_plot,2,j+1)
    plt.axvspan(tP[j], tP[j]+window, color='yellow', alpha=0.5)
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
# Add title, labels, set subplots, and save figure:
t = plt.gcf().text(0.5, 0.98, "Vertical component", horizontalalignment='center')
t = plt.gcf().text(0.5, 0.01, "Time (s)", horizontalalignment='center')
t = plt.gcf().text(0.02, 0.5, label_units, verticalalignment='center', rotation='vertical')
plt.subplots_adjust(wspace=0.35,hspace=0.50,bottom=0.06,top=0.95,left=0.15,right=0.96)
plt.savefig("../Figures/Seismog_vertical.png", dpi=300)



plt.show()
