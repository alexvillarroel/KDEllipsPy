c dimension file to run convm
c contains basic info like delta t, nb signals etc ...
c here dt = 0.47 (60sec/128 sples)
	parameter   (dt=0.25, nsp=5000,ncp=16,nrp=8,ntp=512,
     *  nrtp=nrp*ntp, fr2=0.1,fr1=0.02, bp_order=2)
