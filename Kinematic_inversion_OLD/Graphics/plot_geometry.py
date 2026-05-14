#!/usr/bin/env python3

# Main Python imports:
import numpy as np
import matplotlib.pyplot as plt



#####################################################################
#   SCRIPT TO PLOT THE FAULT GEOMETRY ON A MAP AND CROSS-SECTIONS   #
#####################################################################



# ======================================================================
#  USER PARAMETERS
# ======================================================================
# Parameters:
mu = 4.68075e10   # Shear modulus or rigidity of rocks (N/m2)



# ======================================================================
#  STATIONS AND FAULT GEOMETRY DATA
# ======================================================================
# Read station_in (station distribution):
eq_dat = open("../Stations/station_in",'r')
linea = eq_dat.readline().strip()
(xlat,ylon,zdepth,) = linea.split()[0:3]
titulo = linea[29:]
print(titulo)
xlat = float(xlat)
ylon = 360.+float(ylon)

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
shade = (fault[0:size_fault,3]-fault_min)/(fault_max-fault_min)/2.+0.5
shade = hist[0:size_fault,1]/hist_max

# Open stationn and station (station names and positions):
estaciones = np.loadtxt("../Stations/station")
nombres = open("../Stations/stationn",'r')
nstac = int(np.size(estaciones)/3)
print("There are %i stations:" % (nstac))

estaciones.shape = nstac,3
fault.shape = size_fault,4
fault[:,1] /= 1000.
fault[:,2] /= 1000.
fault[:,3] /= 1000.
estaciones[:,0] /= 1000.
estaciones[:,1] /= 1000.



# ======================================================================
#  FIGURES
# ======================================================================
# Map views of the rupture area:
plt.figure(1, figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(estaciones[:,1],estaciones[:,0],'ro',ms=7)
for i in range(size_fault):
   plt.plot([fault[i,2]],[fault[i,1]],marker='.',color=[shade[i],0,0],ms=0.5*shade[i]+1)
plt.text(0,30,'Fault plane', color='green')
# Add station names:
i=0
for name in nombres.readlines():
   print(i+1, name[0:-1])
   if(name==''): continue
   plt.text(estaciones[i,1]+0.500,estaciones[i,0]-35,name, clip_on=True)
   i+=1
   if(i>=nstac): break
# Plot limits:
plt.xlabel("W <-- Distance (km) --> E")
plt.ylabel("S <-- Distance (km) --> N")
plt.xlim(-500,500)
plt.ylim(-500,500)


plt.title("Study area")

plt.subplot(1,2,2)
for i in range(size_fault):
   plt.plot([fault[i,2]],[fault[i,1]],marker='.',color=[shade[i],0,0],ms=10*shade[i]+1)
# Add rupture details to the plot:
faultg = open("../Faille/faille.in",'r')
(xfault,yfault) = faultg.readline().split()[0:2]
xfault = int(xfault)/1000
yfault = int(yfault)/1000
plt.plot([0.],[0.],'bd',ms=12)
plt.text(-xfault+3,yfault-3,"Max slip = %1.2f (m)" % (hist_max))
plt.text(-xfault+3,yfault-6,"Mo  = %0.4e (Nm)" % (Mo))
plt.text(-xfault+3,yfault-9,"Mw  = %0.2f" % (Mw))
plt.text(xfault-3,yfault-3,"Fault size = %i x %i (km)" % (xfault, yfault), horizontalalignment='right')
plt.text(xfault-3,yfault-6,"Sub-faults = %i x %i" % (float(nx), float(ny)), horizontalalignment='right')
# Add event name to the plot:
eq_dat=open("../Stations/station_in",'r')
linea=eq_dat.readline().strip()
ev_name=linea[29:]
plt.text(0,-xfault+3,ev_name,color=[0,0,0.4],horizontalalignment='center')
# Plot limits:
plt.xlabel("W <-- Distance (km) --> E")
plt.ylabel("S <-- Distance (km) --> N")
plt.xlim(-xfault,xfault)
plt.ylim(-yfault,yfault)
plt.title("Fault plane")
plt.subplots_adjust(top=0.94,bottom=0.11,left=0.06,right=0.99)
plt.savefig("../Figures/Geometry_map.png", dpi=300)



# Cross sections:
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



plt.show()
