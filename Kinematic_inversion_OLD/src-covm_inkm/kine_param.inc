c
c	Parameter file for SMA subroutines
c
        parameter       (maxlayer = 30,
     &                   maxsublayer = 100,
     &                   max_nd = maxlayer * 2,
     &                   maxmoddim = max_nd,
     &                   maxdata = 500000,
     &                   maxwave = 200)
        
	parameter  (dt=0.25,nsp=1000,ncp=40,nrp=40,ntp=512)
        parameter (nrtp=nrp*ntp,bp_order=1)
