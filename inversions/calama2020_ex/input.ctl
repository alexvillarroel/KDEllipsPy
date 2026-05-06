#1. Observed Data Parameters
Time window start (t1)          :  0.0
Time window end (t2)            :  128.0
Number of points (Npts)         :  512
Delta / Time step               :  0.25
Units (1:disp, 2:vel)           :  2

#2. Source Position & Focal Mechanism
Event Name                      :  Calama 2020 Intraplate
Event Origin Date (UTC, YYYY-MM-DD)           :     2020-06-03
Event Origin Time (UTC, HH:MM:SS[.sss])       :     07:35:34.000
Latitude                        :  -23.247
Longitude                       :  -68.53
Depth                           :  123.4
Strike                          :  333.0
Dip                             :  60.0
Rake                            :  -91.0

#3. Fault Plane Parameters
Length along strike (Lx)        :  32000.0
Length along dip (Ly)           :  32000.0
Hypocenter position strike (Hx) :  16000.0
Hypocenter position dip (Hy)    :  16000.0
Number of subfaults along strike (Nx) :  30
Number of subfaults along dip (Ny)    :  30

#4. Ellipse Parameters & Frequency Band
Number of ellipses              :  1
Initial slip                    :  0
Slip shape                      :  1
Frequency 1 (Freq1)             :  0.02
Frequency 2 (Freq2)             :  0.1
Time shift (T0)                 :  3.0

#5. Parameters to Invert
Param 1 : Length of axis 1                 (km)    :  2.0   10.0  1
Param 2 : Length of axis 2                 (km)    :  2.0   10.0  1
Param 3 : Rotation angle                 (x pi)    :  0.0   2.0  1
Param 4 : Position of the center np                :  0.0   1.0  1
Param 5 : Position of the center tp     (x 2pi)    :  0.0   1.0  1
Param 6 : Maximum slip (Dmax)               (m)    :  0.5   4.0  1
Param 7 : Rupture velocity (Vr)          (km/s)    :  0.5   3.5  1

#6. Inversion Process Parameters
Algorithm type                  :  0
Number of iterations            :  7
Sample size for first iteration (SS1) :  100
Sample size for other iterations      :  30
Cells to resample               :  4
Misfit time window              :  0.0
MCMC total steps                :  500
MCMC burn-in                    :  0
MCMC proposal scale             :  0.08
MCMC thinning                   :  1
MCMC chains                     :  1

#7. Moment Tensor (Full MT)
Moment Tensor Flag              :  0
MT Scaling Mode                 :  mt_factored
Mrr                             :  -2.486
Mtt                             :  7.932
Mpp                             :  0.083
Mrt                             :  0.038
Mrp                             :  -0.739
Mtp                             :  0.473
Exponent (iexp)                 :  18.0

#8. Station Parameters
# lat          lon         height   name   use_N  use_E  use_Z
    -27.3612      -70.3390      0.000  A05C      1  1  1
    -26.1479      -70.5987      0.000  AC01      1  1  1
    -26.8355      -69.1291      0.000  AC02      1  1  1
    -22.8528      -70.2023      0.000  PB05      1  1  1
    -22.7058      -69.5717      0.000  PB06      1  1  1
    -20.1411      -69.1535      0.000  PB08      1  1  1
    -23.9048      -69.2900      0.000  PB19      1  1  1
    -20.2393      -70.0540      0.000  T15A      1  1  1

#9. Velocity Model 1D
# thickness    vp        vs       rho      qp      qs
       0.0     5210.0     2990.0     2500.0     600.0     400.0
    2500.0     5370.0     3090.0     2500.0     600.0     400.0
    4500.0     5550.0     3190.0     2500.0     600.0     400.0
    6500.0     5720.0     3290.0     2600.0     600.0     400.0
    8500.0     5890.0     3390.0     2600.0     600.0     400.0
   10500.0     5980.0     3440.0     2600.0     600.0     400.0
   15000.0     6800.0     3750.0     2800.0     600.0     400.0
   20000.0     6810.0     3880.0     2800.0     600.0     400.0
   25000.0     6950.0     3940.0     3000.0     600.0     400.0
   30000.0     6980.0     4050.0     3000.0     600.0     400.0
   35000.0     7110.0     4110.0     3100.0     600.0     400.0
   40000.0     7410.0     4180.0     3300.0     600.0     400.0
   45000.0     7690.0     4300.0     3300.0     600.0     400.0
   50000.0     8050.0     4390.0     3300.0     600.0     400.0
   60000.0     8480.0     4730.0     3400.0     600.0     400.0
   70000.0     8480.0     4780.0     3400.0     600.0     400.0
