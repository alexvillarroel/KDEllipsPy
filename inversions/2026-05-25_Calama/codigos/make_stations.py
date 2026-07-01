#!/usr/bin/env python
"""Lee SAC/{ACC,VEL,DISP} y escribe stations.txt: lat lon elev(0) nombre."""
import sys, glob, os
from obspy import read

base = sys.argv[1] if len(sys.argv) > 1 else "SAC"
sub = next((d for d in ("ACC", "VEL", "DISP") if glob.glob(f"{base}/{d}/*.sac")), None)
assert sub, f"sin .sac en {base}/{{ACC,VEL,DISP}}"

st = {}
for f in glob.glob(f"{base}/{sub}/*.sac"):
    h = read(f)[0].stats.sac
    st[read(f)[0].stats.station] = (h.stla, h.stlo)  # ponytail: re-read, files son pocos

with open("stations.txt", "w") as out:
    for name in sorted(st):
        lat, lon = st[name]
        out.write(f"{lat:.5f} {lon:.5f} 0.0000 {name}\n")

print(f"stations.txt: {len(st)} estaciones desde {base}/{sub}")
