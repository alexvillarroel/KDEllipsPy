# Main Python imports:
import math
# Obspy imports:
from obspy.taup import TauPyModel
from obspy.clients.iris import Client



# Import station coordinates:
infile = open("../Stations/station_in", 'r')
count = 1
stinfo = []
for row in infile:
    cols = row.split()
    
    # Retrieve event location:
    if count == 1:
        evlat = float(cols[0])
        evlon = float(cols[1])
        evdep = float(cols[2])
        
    # Retrieve station locations:
    if count >= 3:
        stlat = float(cols[0])
        stlon = float(cols[1])
        stname = str(cols[3])
        stinfo.append([stlat, stlon, stname])
        
    count = count+1
    
    
# Main calculations:
client = Client()
model = TauPyModel(model="iasp91")
out = open("../Event/azi_times.txt", 'w')
for i in range(len(stinfo)):
    # Calculate azimuth and distance:
    result = client.distaz(stinfo[i][0], stinfo[i][1], evlat, evlon)
    distdeg = float(result['distance'])
    azi = float(result['azimuth'])

    # Calculate theoretical arrival times of P and S waves:
    try:
        P_time = model.get_travel_times(evdep, distdeg, phase_list="P")
        Pt = P_time[0].time
    except IndexError as ex:
        P_time = model.get_travel_times(evdep, distdeg, phase_list="p")
        Pt = P_time[0].time
    try:
        S_time = model.get_travel_times(evdep, distdeg, phase_list="S")
        St = S_time[0].time
    except IndexError as ex:
        S_time = model.get_travel_times(evdep, distdeg, phase_list="s")
        St = S_time[0].time
        
    # Output file:
    out.write("%0.4f %0.1f %0.1f\n" % (math.radians(azi), Pt-1, St-1))
    
out.close()
