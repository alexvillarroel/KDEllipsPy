c @(#) fsource.F       FSOURCE 1.5      03/19/96 1
c**********************************************************
c      FSOURCE
c      
c      Different kind of source function defined in the 
c      frequency domain
c      Source function for dislocation (step, ramp, haskell)
c      are normalized so that in far-field the low-frequency 
c      level is proportionnal to the seismic moment with a factor
c      equals to: Rad/(4.PI.rho.beta^3) * 1./r 
c      where Rad= Radiation coefficient with (possibly) 
c      free surface effect
c
c      input:      
c            type       ->       see below
c            omega      ->      angular frequency
c            t0,t1      ->      time constant when needed
c            dt      ->      sampling rate
c**********************************************************

      function      fsource (type, omega, t0, t1, dt)

      implicit       none
      integer            type
      real            pi,pi2,dt,t0,t1
      real*8            uur,uui,trise,trupt
      complex*16      fsource,uu,uex,uxx,omega,shx,ai

      common      /par/      ai,pi,pi2

c TYPE=0               Source = Dirac en deplacement
       if (type.eq.0) then
          fsource=1
       endif
 
c TYPE=1        Source = Ricker en deplacement
      if (type.eq.1) then
          uu=omega*t0
          uu=uu*uu/pi2/pi2
          uu=exp(-uu)
          uu=omega*omega*uu*dt
        fsource= uu
        endif
 
c TYPE=2        Source = step en deplacement
c            2 steps possibles (1) real=1/(ai*omega)
c                          (2) bouchon's
        if (type.eq.2) then
          shx=exp(omega*pi*t0/2.)      !Bouchon's
          shx=1./(shx-1./shx)
          uu=-ai*t0*pi*shx
          fsource= uu
        endif
c TYPE=7        Source = step en deplacement
      if (type.eq.7) then
        uu=1./ai/omega
          fsource= uu
      endif

c TYPE=3        Source = fichier 'axi.sou'
c               sismogramme dans l'unite choisie dans le fichier
        if (type.eq.3) then
          read(13,*) uur,uui
          fsource=cmplx(uur,uui)
        endif
 
c TYPE=4        Source = triangle en deplacement
        if (type.eq.4) then
        trise=t0
        trupt=t0
        uu=ai*omega*trise
          uu=(1.-exp(-uu))/uu            ! ramp
          uxx=ai*omega*trupt/2.            ! finite fault
          uex=exp(uxx)
          uxx=(uex-1./uex)/uxx/2.
          fsource=uu*uxx/(ai*omega)
        endif
 
c TYPE=5        Source = rampe causale
c               rise time T=t0
        if (type.eq.5) then
        trise=t0
          uu=ai*omega*trise
          uu=(1.-exp(-uu))/uu
          fsource=uu/(ai*omega)
        endif
      
c TYPE=6,8        Source = modele d'haskell, trapezoide
c      1 ere cste de temps rise time: riset
c      2 eme cste de temps, duree de la rupture
c        trupt = Length/2/rupt_velocity (Haskell)
        if ((type.eq.6).or.(type.eq.8)) then
        trise=t0
        trupt=t1
        uu=ai*omega*trise
          uu=(1.-exp(-uu))/uu            ! ramp
          uxx=ai*omega*trupt/2.            ! finite fault
          uex=exp(uxx)
          uxx=(uex-1./uex)/uxx/2.
          fsource=uu*uxx/(ai*omega)
        endif
    
        return
      end
        

