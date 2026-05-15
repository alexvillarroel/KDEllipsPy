#!/usr/bin/env python
# %% esto es un python interactivo

import obspy

data = obspy.read("./RAW/*.sac")
for tr in data:
    print(tr.stats.starttime)
