#!/usr/bin/env python3

# Main Python imports:
import math
import numpy as np
import matplotlib.pyplot as plt
import os

##################################################
#   SCRIPT TO PLOT OBSERVED AND MODELED TRACES   #
#               IN ONE FIGURE                    #
##################################################

# ======================================================================
#   USER PARAMETERS
# ======================================================================
# Parameters:
# 1 = Displacement (cm), 2 = Velocity (cm/s)
units_param = 2

# ======================================================================
#   SEISMOGRAM DATA
# ======================================================================
# Retrieve time increment (dt) and number of points (npts) from the file
kine_param_path = os.path.join("..", "src-covm_inkm", "kine_param.inc")
with open(kine_param_path, 'r') as f:
    dt = None
    npts = None
    for row in f:
        cols = row.split(',')
        if "parameter" in cols[0] and "dt=" in cols[0]:
            dt = float(cols[0].split('=')[1])
            npts = int(cols[4].strip().split('=')[1].replace(')', ''))
    
    if dt is None or npts is None:
        raise ValueError("Could not find dt or npts in kine_param.inc")

# Retrieve the number of stations and their names
station_file_path = os.path.join("..", "Stations", "stationn")
with open(station_file_path, 'r') as f:
    stations = [line.strip() for line in f.readlines()]
nsta = len(stations)

# Loading observed and modeled data
data_path = os.path.join("..", "Event", "kine_files")
real_x = 100 * np.loadtxt(os.path.join(data_path, "real_disp_x"))
sismog_x = 100 * np.loadtxt(os.path.join(data_path, "best_disp_x"))
real_y = 100 * np.loadtxt(os.path.join(data_path, "real_disp_y"))
sismog_y = 100 * np.loadtxt(os.path.join(data_path, "best_disp_y"))
real_z = 100 * np.loadtxt(os.path.join(data_path, "real_disp_z"))
sismog_z = 100 * np.loadtxt(os.path.join(data_path, "best_disp_z"))

# ======================================================================
#   PLOT PARAMETERS
# ======================================================================
# Set time array, titles and names of the stations:
time = np.arange(npts) * dt

# Physical units of seismograms:
if units_param == 1:
    label_units = "Displacement (cm)"
    units = "(cm)"
elif units_param == 2:
    label_units = "Velocity (cm/s)"
    units = "(cm/s)"
else:
    label_units = ""
    units = ""

# Function to control vertical limits:
def round_to_1(x):
    if x == 0:
        return 0
    return round(x, -int(math.floor(math.log10(x))))

# ======================================================================
#   FIGURES
# ======================================================================
# Plot all components and stations:
plt.figure(figsize=(9, 10))
for j in range(nsta):
    # Set common vertical limits for the current station:
    start_idx = npts * j
    end_idx = npts * (j + 1)
    
    mobs_x = np.max(np.abs(real_x[start_idx:end_idx]))
    mmod_x = np.max(np.abs(sismog_x[start_idx:end_idx]))
    mobs_y = np.max(np.abs(real_y[start_idx:end_idx]))
    mmod_y = np.max(np.abs(sismog_y[start_idx:end_idx]))
    mobs_z = np.max(np.abs(real_z[start_idx:end_idx]))
    mmod_z = np.max(np.abs(sismog_z[start_idx:end_idx]))
    
    maxamps = [mobs_x, mmod_x, mobs_y, mmod_y, mobs_z, mmod_z]
    maxv = np.max(maxamps)
    limit = round_to_1(maxv * 1.5)
    
    # North component:
    ax1 = plt.subplot(nsta, 3, 3 * j + 1) 
    ax1.plot(time, real_x[start_idx:end_idx], 'b')
    ax1.plot(time, sismog_x[start_idx:end_idx], 'r')
    ax1.set_xlim(0, npts * dt)
    ax1.set_ylim(-limit, limit)
    ax1.set_yticks([-limit, -limit / 2, 0, limit / 2, limit])
    ax1.text(0.05, 0.95, stations[j], transform=ax1.transAxes, va='top', ha='left')

    # East component:
    ax2 = plt.subplot(nsta, 3, 3 * j + 2) 
    ax2.plot(time, real_y[start_idx:end_idx], 'b')
    ax2.plot(time, sismog_y[start_idx:end_idx], 'r')
    ax2.set_xlim(0, npts * dt)
    ax2.set_ylim(-limit, limit)
    ax2.set_yticks([-limit, -limit / 2, 0, limit / 2, limit])
    ax2.text(0.05, 0.95, stations[j], transform=ax2.transAxes, va='top', ha='left')

    # Vertical component:
    ax3 = plt.subplot(nsta, 3, 3 * j + 3) 
    ax3.plot(time, real_z[start_idx:end_idx], 'b')
    ax3.plot(time, sismog_z[start_idx:end_idx], 'r')
    ax3.set_xlim(0, npts * dt)
    ax3.set_ylim(-limit, limit)
    ax3.set_yticks([-limit, -limit / 2, 0, limit / 2, limit])
    ax3.text(0.05, 0.95, stations[j], transform=ax3.transAxes, va='top', ha='left')

# Add titles, labels, set subplots, and save figure:
plt.suptitle("All Seismograms", fontsize=16)
plt.figtext(0.22, 0.95, "N-S", ha='center', fontsize=12)
plt.figtext(0.55, 0.95, "E-W", ha='center', fontsize=12)
plt.figtext(0.88, 0.95, "Z", ha='center', fontsize=12)
plt.figtext(0.5, 0.02, "Time [s]", ha='center', fontsize=12)
plt.figtext(0.02, 0.5, label_units, va='center', rotation='vertical', fontsize=12)
plt.subplots_adjust(wspace=0.53, hspace=0.55, left=0.12, right=0.98, bottom=0.06, top=0.94)

# Create the output directory if it doesn't exist
output_dir = os.path.join("..", "Figures")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# plt.savefig(os.path.join(output_dir, "Seismog_all.png"), dpi=300)

plt.show()