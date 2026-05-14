#!/usr/bin/env python3

# Main Python imports:
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Cartopy for 2D map plotting
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Rockhound for fetching Slab 2.0 data
import rockhound

# Import for 3D plotting
from mpl_toolkits.mplot3d import Axes3D


#####################################################################
#   SCRIPT TO PLOT THE FAULT GEOMETRY ON A MAP AND CROSS-SECTIONS   #
#####################################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
mu = 4.68075e10   # Shear modulus or rigidity of rocks (N/m2)


# ======================================================================
#  SLAB 2.0 DATA FETCHING
# ======================================================================
print("Fetching Slab 2.0 data for South America...")
slab2 = rockhound.fetch_slab2(zone="south_america")
slab2['longitude'] = slab2['longitude'] - 360 # Adjust longitude to -180 to 180 range
slab2['depth'] = slab2['depth'] / 1000       # Convert depth from meters to kilometers
print("Slab 2.0 data fetched.")


# ======================================================================
#  STATIONS AND FAULT GEOMETRY DATA
# ======================================================================
# Read station_in (station distribution):
eq_dat = open("../Stations/station_in",'r')
linea = eq_dat.readline().strip()
(xlat,ylon,zdepth,) = linea.split()[0:3]
print(xlat,ylon,zdepth)
titulo = linea[29:]
print(titulo)
xlat = float(xlat)
ylon = float(ylon) # Cartopy handles negative lon well.

# Read faille.in (fault geometry):
input_file = open("../Faille/faille.in",'r')
input_file.readline()
(xhypo,yhypo) = input_file.readline().split()[0:2]
(zhypo,) = input_file.readline().split()[0:1]
(nx,ny) = input_file.readline().split()[0:2]
xhypo = float(xhypo)/1000
yhypo = float(yhypo)/1000
zhypo = float(zhypo)/1000

faille_file = open("../Faille/faille.in",'r')
(xsize,ysize) = faille_file.readline().split()[0:2]
nx1 = int(xsize)/1000
ny1 = int(ysize)/1000

print("Hypocenter on the fault plane (x, y, depth) in (km) =", xhypo, yhypo, zhypo)
print('Sub-faults on the fault plane =:', nx, "x", ny)
print("Fault size:", nx1, "x", ny1, "(km)")

# Read source.dat (source parameters):
fault = np.loadtxt("../AXITRA_F90/source.dat")
size_fault = int(nx)*int(ny)
fault_min = min(fault[0:size_fault,3])
fault_max = max(fault[0:size_fault,3])

# Read axi.hist file, calculate Mo and Mw, and retrieve strike, dip and rake:
hist = np.loadtxt("../Event/axi.hist")
hist.shape = size_fault,8
slip = 0
for i in range(size_fault):
   slip += hist[i,1]
Mo = mu*slip*(hist[1,5]*hist[1,6])
Mw = (np.log10(Mo)-9.1)/1.5
print("Mo = %0.4e (Nm)" % (Mo))
print("Mw = %0.2f" % (Mw))
strike = hist[1,2]
dip = hist[1,3]
rake = hist[1,4]

# Retrieve the maximum slip on the fault:
hist_max = max(hist[0:size_fault,1])
print("Maximum slip = %0.2f (m)" % (hist_max))

# Transform slip into shade:
shade = hist[0:size_fault,1]/hist_max

# Open stationn and station (station names and positions):
estaciones = np.loadtxt("../Stations/station")
nombres = open("../Stations/stationn",'r')
nstac = int(np.size(estaciones)/3)
print("There are %i stations:" % (nstac))

estaciones.shape = nstac,3
fault.shape = size_fault,4
fault[:,1] /= 1000. # Convert fault 'y' (N-S) to km
fault[:,2] /= 1000. # Convert fault 'x' (E-W) to km
fault[:,3] /= 1000. # Convert fault 'z' (depth) to km
estaciones[:,0] /= 1000. # Convert station 'y' (N-S) to km
estaciones[:,1] /= 1000. # Convert station 'x' (E-W) to km


# --- Helper function to convert km offsets to latitude and longitude ---
def km_to_lonlat(ref_lon, ref_lat, dx_km, dy_km):
    """
    Converts km offsets (East, North) from a reference point to lon/lat.
    - dx_km: Distance East in km.
    - dy_km: Distance North in km.
    """
    lat_rad = np.radians(ref_lat)
    delta_lat = dy_km / 111.132
    delta_lon = dx_km / (111.320 * np.cos(lat_rad))
    return ref_lon + delta_lon, ref_lat + delta_lat

# Convert fault and station coordinates from local km to lon/lat
# Note: Original plot used fault[:,2] for E-W (x-axis) and fault[:,1] for N-S (y-axis)
fault_lons, fault_lats = km_to_lonlat(ylon, xlat, fault[:,2], fault[:,1])

# Note: Original plot used estaciones[:,1] for E-W (x-axis) and estaciones[:,0] for N-S (y-axis)
station_lons, station_lats = km_to_lonlat(ylon, xlat, estaciones[:,1], estaciones[:,0])

# Prepare data for scatter plot (more efficient for many points)
fault_colors = np.zeros((size_fault, 3))
fault_colors[:,0] = shade # Red channel corresponds to slip
fault_sizes_overview = 0.5 * shade + 1
fault_sizes_zoom = 10 * shade + 1

# Define the map projection centered on the event
projection = ccrs.AzimuthalEquidistant(central_longitude=ylon, central_latitude=xlat)
geodetic_proj = ccrs.Geodetic()


# ======================================================================
#  FIGURE 1: MAP VIEWS (with Cartopy)
# ======================================================================
fig1 = plt.figure(1, figsize=(14, 6))

# Subplot 1: Wider "Study Area" view
ax1 = fig1.add_subplot(1, 2, 1, projection=projection)
ax1.set_title("Study Area")

extent_m = 500 * 1000
ax1.set_extent([-extent_m, extent_m, -extent_m, extent_m], crs=projection)

ax1.add_feature(cfeature.LAND, edgecolor='black', zorder=0)
ax1.add_feature(cfeature.OCEAN, zorder=0)
ax1.add_feature(cfeature.COASTLINE, zorder=1)
gl = ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, zorder=2)
gl.top_labels = False
gl.right_labels = False

ax1.plot(station_lons, station_lats, 'ro', ms=7, transform=geodetic_proj, label='Stations')
ax1.scatter(fault_lons, fault_lats, c=fault_colors, s=fault_sizes_overview, marker='.', transform=geodetic_proj, zorder=3)

station_names_list = nombres.readlines()
for i in range(nstac):
    name = station_names_list[i].strip()
    if name:
        ax1.text(station_lons[i] + 0.01, station_lats[i] - 0.01, name, transform=geodetic_proj, clip_on=True)


# Subplot 2: Zoomed-in "Fault Plane" view
ax2 = fig1.add_subplot(1, 2, 2, projection=projection)
ax2.set_title("Fault Plane")

xfault = nx1
yfault = ny1
ax2.set_extent([-xfault*1000, xfault*1000, -yfault*1000, yfault*1000], crs=projection)

ax2.add_feature(cfeature.LAND, edgecolor='black', zorder=0)
ax2.add_feature(cfeature.OCEAN, zorder=0)
ax2.add_feature(cfeature.COASTLINE, zorder=1)
gl2 = ax2.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, zorder=2)
gl2.top_labels = False
gl2.right_labels = False

ax2.scatter(fault_lons, fault_lats, c=fault_colors, s=fault_sizes_zoom, marker='.', transform=geodetic_proj, zorder=3)
ax2.plot(ylon, xlat, 'bd', ms=12, transform=geodetic_proj, zorder=4)

ax2.text(0.02, 0.95, f"Max slip = {hist_max:.2f} (m)", transform=ax2.transAxes, va='top')
ax2.text(0.02, 0.90, f"Mo  = {Mo:0.4e} (Nm)", transform=ax2.transAxes, va='top')
ax2.text(0.02, 0.85, f"Mw  = {Mw:0.2f}", transform=ax2.transAxes, va='top')
ax2.text(0.98, 0.95, f"Fault size = {int(xfault)} x {int(yfault)} (km)", transform=ax2.transAxes, ha='right', va='top')
ax2.text(0.98, 0.90, f"Sub-faults = {int(float(nx))} x {int(float(ny))}", transform=ax2.transAxes, ha='right', va='top')
ax2.text(0.5, 0.02, titulo, color=[0,0,0.4], transform=ax2.transAxes, ha='center', va='bottom')


plt.subplots_adjust(top=0.94,bottom=0.11,left=0.06,right=0.99, wspace=0.3)
plt.savefig("../Figures/Geometry_map_cartopy.png", dpi=300)


# ======================================================================
#  FIGURE 2: CROSS SECTIONS (unchanged)
# ======================================================================
plt.figure(2, figsize=(12,5))
plt.subplot(1,2,1)
for i in range(size_fault):
   plt.plot([fault[i,1]],[-fault[i,3]],marker='.',color=[shade[i],0,0],ms=10*shade[i]+1)
depth=-float(zhypo)
plt.plot([0.],[depth],'bd',ms=12)
plt.xlabel("S <-- Distance (km) --> N")
plt.ylabel("Depth (km)")
plt.xlim(-xfault,xfault)
plt.ylim(depth-yfault,depth+yfault)
plt.title("Cross-section along longitude")

plt.subplot(1,2,2)
for i in range(size_fault):
   plt.plot([fault[i,2]],[-fault[i,3]],marker='.',color=[shade[i],0,0],ms=10*shade[i]+1)
depth=-float(zhypo)
plt.plot([0.],[depth],'bd',ms=12)
plt.xlabel("W <-- Distance (km) --> E")
plt.ylabel("Depth (km)")
plt.xlim(-xfault,xfault)
plt.ylim(depth-yfault,depth+yfault)
plt.title("Cross-section along latitude")
plt.subplots_adjust(top=0.94,bottom=0.11,left=0.06,right=0.99)
plt.savefig("../Figures/Geometry_cross_section.png", dpi=300)


# ======================================================================
#  FIGURE 3: 3D Visualization of Fault within Slab
# ======================================================================
print("Generating 3D visualization...")
fig3 = plt.figure(3, figsize=(10, 8))
ax3 = fig3.add_subplot(111, projection='3d')
ax3.set_title(f"3D View: Fault within {titulo} Slab")

# Plot the Slab 2.0 surface
# Reshape the slab data for surf plotting
slab_lons_grid = slab2.longitude.values.reshape(slab2.depth.shape)
slab_lats_grid = slab2.latitude.values.reshape(slab2.depth.shape)
slab_depth_grid = slab2.depth.values.reshape(slab2.depth.shape)

# Create a colormap for slab depth
cmap_slab = plt.cm.viridis_r
norm_slab = mcolors.Normalize(vmin=slab_depth_grid.min(), vmax=slab_depth_grid.max())
slab_colors = cmap_slab(norm_slab(slab_depth_grid))

# Plot the slab surface
ax3.plot_surface(slab_lons_grid, slab_lats_grid, -slab_depth_grid,
                 facecolors=slab_colors, rstride=5, cstride=5,
                 linewidth=0, antialiased=False, alpha=0.6)


# Plot the fault plane (points) with slip values
# Depth for fault points: fault[:,3] is already in km
fault_depths = -fault[:,3] # Negative for plotting depth downwards

ax3.scatter(fault_lons, fault_lats, fault_depths,
            c=shade, cmap='Reds', s=fault_sizes_zoom * 2, marker='o', depthshade=False,
            label='Fault Slip (m)', edgecolors='k', linewidths=0.5)

# Plot the epicenter
hypocenter_depth = -zhypo
ax3.scatter(ylon, xlat, hypocenter_depth,
            color='blue', marker='*', s=300, label='Hypocenter', edgecolors='k', linewidths=1)

# Set labels
ax3.set_xlabel("Longitude")
ax3.set_ylabel("Latitude")
ax3.set_zlabel("Depth (km)")

# Invert z-axis to make depth increase downwards
ax3.invert_zaxis()

# Set reasonable limits based on fault and local slab extent
lon_min = min(fault_lons.min(), slab_lons_grid.min()) - 0.5
lon_max = max(fault_lons.max(), slab_lons_grid.max()) + 0.5
lat_min = min(fault_lats.min(), slab_lats_grid.min()) - 0.5
lat_max = max(fault_lats.max(), slab_lats_grid.max()) + 0.5
depth_min = min(-slab_depth_grid.max(), hypocenter_depth - yfault) - 10
depth_max = max(-slab_depth_grid.min(), hypocenter_depth + yfault) + 10

ax3.set_xlim(lon_min, lon_max)
ax3.set_ylim(lat_min, lat_max)
ax3.set_zlim(depth_max, depth_min) # Note inversion for depth

# Adjust camera angle for better view
ax3.view_init(elev=30, azim=-120)

plt.colorbar(plt.cm.ScalarMappable(norm=norm_slab, cmap=cmap_slab), ax=ax3, shrink=0.5, aspect=10, label='Slab Depth (km)')
plt.legend()
plt.tight_layout()
plt.savefig("../Figures/Geometry_3D_Slab_Fault.png", dpi=300)
print("3D visualization saved to ../Figures/Geometry_3D_Slab_Fault.png")

plt.show()