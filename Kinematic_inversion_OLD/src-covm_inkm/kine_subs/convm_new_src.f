      subroutine readaxires
      include 'kine_param.inc'
      character headfile*50,sourcefile*10,statfile*10,
     &       histfile*50, filename*50, axires*50

      complex*16 uxf(6),uyf(6),uzf(6),ai
      
      complex uuxf(6,nrp,nsp,ntp/4),uuyf(6,nrp,nsp,ntp/4),
     $     uuzf(6,nrp,nsp,ntp/4) 
      common/hist/uuxf,uuyf,uuzf
        namelist /input/ nc,nfreq,tl,aw,nr,ns,xl,t0,fr1,fr2,
     &     ikmax,uconv,sourcefile,statfile
c
c  Lecture de axi.res 
c      
      headfile='axi.head'
      axires='axi.res'
      open (10,form='formatted',file=headfile)

      read(10,input)
      if ((ns.gt.nsp).or.(nr.gt.nrp).or.(nc.gt.ncp)) then
         write(*,*) "parametres sous-dimensionnes"
         stop
      endif
      close(10)
c
c ========================================
c       read axi.res
c ========================================
c
      open (12,form='unformatted',file='axi.res')
      print *, 'XXXXX read axi.res with nfreq,ns,nr = ',nfreq,ns,nr
      do jf=1,nfreq
         do is=1,ns
            do ir=1,nr
c       
c       FAUT lire par block de six complexes a cause du format 
c       binaire fortran
c       
           read(12,end=200)(uxf(it),it=1,6)
           read(12,end=200)(uyf(it),it=1,6)
           read(12,end=200)(uzf(it),it=1,6)
           do it=1,6
              uuxf(it,ir,is,jf)=uxf(it)
              uuyf(it,ir,is,jf)=uyf(it)
              uuzf(it,ir,is,jf)=uzf(it)
           enddo
            enddo
         enddo
      enddo   
      goto 201
 200      print *,'Merde, axi.res is not big eneugh'
      nfreq=jf-1
 
 201      continue
      print *,' read nfreq = ',nfreq,' axitra lines' 
      close(12)
        return
        end
c
c =============================================
c       seismograms
c =============================================
c
      subroutine seismograms(ndata,nwave,predicted_data,ncomp)       
      include 'kine_param.inc'
      real         lat0,lon0
      parameter   (lon0=-118,lat0=34,stdlat1=33,stdlat2=45)

      real tax(ntp),tay(ntp),taz(ntp)

      character   filein*50,sourcefile*10,statfile*10,
     1     histfile*50, filename*50, axires*50, axistationn*50
      integer     iwk(ntp),jf,ir,is,it,nc,ns,nr,nfreq,ikmax,mm,nt,
     1     io,isc(nsp)
      real  hc(ncp),vp(ncp),vs(ncp),rho(ncp),delay(nsp),
     1     xr(nrp),yr(nrp),zr(nrp),a(6,nsp),qp(ncp),qs(ncp),
     2     tl,xl,uconv,mu(nsp),strike(nsp),dip(nsp),
     3     rake(nsp),disp(nsp),hh,zsc,dfreq,freq,rvel(nsp),
     4     pi,pi2,xs(nsp),ys(nsp),zs(nsp),aw,ck,xmm,cc,
     5     uout(ntp),pas,width(nsp),length(nsp),xref,yref,
     7     lat,long,t1,momo(nsp),moto
      complex*16  ux(ntp,nrp),uy(ntp,nrp),uz(ntp,nrp),omega,
     1     ai,deriv,fsource,us,uux,uuy,uuz
      
      real*4      predicted_data(maxdata,maxwave)
      common      /par/ai,pi,pi2
      namelist    /input/ nc,nfreq,tl,aw,nr,ns,xl,t0,fr1,fr2,
     &     ikmax,uconv,sourcefile,statfile
      data  ux,uy,uz/nrtp*(0.,0.),nrtp*(0.,0.),
     &                     nrtp*(0.,0.)/
      complex uuxf(6,nrp,nsp,ntp/4),uuyf(6,nrp,nsp,ntp/4),
     $     uuzf(6,nrp,nsp,ntp/4) 
      common/hist/uuxf,uuyf,uuzf
      integer filt_order, nf_ord, bb_channels
!      integer     is(1600)
!      real        phi_pi(1600), slip(1600), delta_pi(1600)
!      real        rake(1600), ddip(1600), dstk(1600),ts(1600) 
!      common /axihist/is, slip, phi_pi, delta_pi, rake,
!     &        ddip, dstk, ts
 
      pi=3.1415926535
      pi2=2.*pi
      ai=(0.,1.)
      dt0=0.
      histfile="axi.hist"
      filein='axi.head'

c      print *,'in seismograms'
c      print *, 'lee axi.res with nfreq,ns,nr = ',nfreq,ns,nr
 
      open (10,form='formatted',file=filein)

      read(10,input)
      if ((ns.gt.nsp).or.(nr.gt.nrp).or.(nc.gt.ncp)) then
         write(0,*) "parametres sous-dimensionnes"
         stop
      endif
c      print *,'leyo input'
      open (13,form='formatted',file=sourcefile)
      open (14,form='formatted',file=statfile)
      open (15,form='formatted',file=histfile)
      ics=5
!      print *,'T0,fr1,fr2,dt = ',t0,fr1,fr2,dt
c Why is it 1 and then 0. 
c icc = 1 DISPLACEMENT. 2 VELOCITY
      icc=2
      icc=icc-1

c++++++++++++
c     Milieu, stations et sources
c++++++++++++
      read(10,*)
      do ic=1,nc
         read(10,*)hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)
!         print *,'leyo 10 ic, hc = ',ic, hc(ic)
      enddo
!      read(10,*)
      do iso=1,ns
         read(13,*) index,xs(iso),ys(iso),zs(iso)
         indexin=-1
         rewind(15)
         do while (indexin.ne.index)
            read(15,*) indexin,disp(iso),strike(iso),dip(iso),
     1           rake(iso),width(iso),length(iso),delay(iso)
         enddo
         delay(iso)=delay(iso)+dt0
      enddo
      close(13)
      close(15)
c     write(*,*) 'nr is   ', nr
      do ir=1,nr
         read(14,*) xr(ir),yr(ir),zr(ir)
c     reference: x = nord; y = est
         
      enddo
      close(14)
c++++++++++++
c     Parametres sources
c++++++++++++
c     convertion depth -> thickness
      if (hc(1).eq.0.) then
         do ic=1,nc-1
            hc(ic)=hc(ic+1)-hc(ic)
         enddo
      endif
      
      moto=0
      do is=1,ns
         hh=0.
         isc(is)=1
         zsc=zs(is)
         do ic=1,nc-1
            hh=hc(ic)
            if (zsc.gt.hh) then
               zsc=zsc-hh
               isc(is)=ic+1
            else
               goto 91
            endif
         enddo
 91      continue
c
c         print *,'la source est dans la couche ',isc(is)
         mu(is)=vs(isc(is))*vs(isc(is))*rho(isc(is))
         call cmoment (mu(is),strike(is),dip(is),
     &        rake(is),disp(is),width(is)*length(is),a(1,is))
      enddo
      
c++++++++++++
c     Lecture fonctions de transfert
c++++++++++++
      xmm=log(real(nfreq))/log(2.)
!      mm=int(xmm+0.05)+2
      mm=int(xmm+0.05)+2
      nt=2**mm
!      print *,'nt,ndata,nr,nwave, ncomp =',nt,ndata,nr,nwave,ncomp
!      print *,'maxdata,maxwave =',maxdata,maxwave
      if (nt.gt.ntp) then
         write(0,*) 'tableau sous dimensionnes !!!'
         stop
      endif
      dfreq=1./tl
      pas=tl/nt
      aw=-pi*aw/tl
c      print *,'tl,nt, pas = ',tl,nt,pas
      do jf=1,nt
          do ir=1,nr
            ux(jf,ir)=cmplx(0.,0.)       
            uy(jf,ir)=cmplx(0.,0.)       
            uz(jf,ir)=cmplx(0.,0.)       
         enddo
      enddo   
  
      do  jf=1,nfreq
         freq=float(jf-1)/tl
         omega=cmplx(pi2*freq,aw)
         if(icc==0)then 
            deriv=1.
         else
            deriv=(ai*omega)**icc
         endif
         
         do is=1,ns
            us=fsource(ics,omega,t0,pas) *
     1           deriv * exp(-ai*omega*delay(is))
            do ir=1,nr
               uux=0.
               uuy=0.
               uuz=0.
               do it=1,6
                  uux = uux + uuxf(it,ir,is,jf)*a(it,is)
                  uuy = uuy + uuyf(it,ir,is,jf)*a(it,is)
                  uuz = uuz + uuzf(it,ir,is,jf)*a(it,is)
               enddo
               ux(jf,ir)=ux(jf,ir) + uux * us 
               uy(jf,ir)=uy(jf,ir) + uuy * us
               uz(jf,ir)=uz(jf,ir) + uuz * us
            enddo      !fin boucle recepteur
         enddo                  !fin boucle source
      enddo 
      close(10)

c++++++++++++
c                Calcul des sismogrammes
c++++++++++++

      iorder=1  
      do  ir=1,nr
      
c     on complete le spectre pour les hautes frequences
c     avec inversion du signe de la partie imaginaire pour
c     la FFT inverse qui n-existe pas avec fft2cd

         do ifr=nt+2-nfreq,nt
            ux(ifr,ir)=dconjg(ux(nt+2-ifr,ir))
            uy(ifr,ir)=dconjg(uy(nt+2-ifr,ir))
            uz(ifr,ir)=dconjg(uz(nt+2-ifr,ir))
         enddo
         
         call fft2cd(ux(1,ir),mm,iwk)
         call fft2cd(uy(1,ir),mm,iwk)
         call fft2cd(uz(1,ir),mm,iwk)
         
         do it=1,nt
            ck=float(it-1)/nt
            cc=exp(-aw*tl*ck)/tl
            ux(it,ir)=ux(it,ir)*dble(cc)
            uy(it,ir)=uy(it,ir)*dble(cc)
            uz(it,ir)=uz(it,ir)*dble(cc)
         enddo
         
         do it=1,nt
            tax(it)=real(real(ux(it,ir)))
            tay(it)=real(real(uy(it,ir)))
            taz(it)=real(real(uz(it,ir)))
         enddo
         
c     20/11/04               
c     now let s filter 
c     on filtre les signaux de fr1 a fr2 Hertz
c     filtre obtenu par Sophie Peyrat via Anthony Sladen
c     les param dt et fr1 et fr2 sont ocntenus dans
c     dimension.inc
         call XAPIIR(taz,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,pas,1)
         call XAPIIR(tay,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,pas,1)
         call XAPIIR(tax,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,pas,1)


         do iter=1,nt
            predicted_data(iorder,1)=tax(iter)
            predicted_data(iorder,2)=tay(iter)
            predicted_data(iorder,3)=taz(iter)
            iorder=iorder+1
         enddo
      
      enddo

      return
      end
c
c ==================================================
c  calcule du moment
c ==================================================
c
      subroutine cmoment (mu,strike,dip,
     1             rake,disp,surf,a)

      implicit none
    
      real xmoment,mu,strike,dip,rake,disp,surf,a(6),pi2,
     1          sd,cd,sp,cp,s2p,s2d,c2p,c2d,x1,x2,x3,x4,x5,x6,cl,sl,pi
      complex*16 ai
      common      /par/ai,pi,pi2

      if (surf.eq.0.) then
         xmoment=disp
      else
         xmoment=mu*disp*surf
      endif
c     write(6,*) "moment (Nm):",xmoment
c     write(6,*) "moment (Dyne.cm):",xmoment*1.e7
      strike=strike*pi/180.
      dip=dip*pi/180.
      rake=rake*pi/180.
      sd=sin(dip)
      cd=cos(dip)
      sp=sin(strike)
      cp=cos(strike)
      sl=sin(rake)
      cl=cos(rake)
      s2p=2.*sp*cp
      s2d=2.*sd*cd
      c2p=cp*cp-sp*sp
      c2d=cd*cd-sd*sd
      
c     Coefficient pour les sources Mxx,Mxy,Mxz,Myy,Myz,Mzz
      x1 =-(sd*cl*s2p + s2d*sl*sp*sp)*xmoment
      x2 = (sd*cl*c2p + s2d*sl*s2p/2.)*xmoment
      x3 =-(cd*cl*cp  + c2d*sl*sp)*xmoment
      x4 = (sd*cl*s2p - s2d*sl*cp*cp)*xmoment
      x5 =-(cd*cl*sp  - c2d*sl*cp)*xmoment
      x6 =             (s2d*sl)*xmoment
      
c     Coefficient pour les sources bis (5 dislocations elementaires
c     et une source isotrope)
      a(1) = x2
      a(2) = x3
      a(3) =-x5
      a(4) = (-2.*x1 + x4 + x6)/3.
      a(5) = (x1 -2*x4 + x6)/3.
      a(6) = 0.
      return
      end
cFSOURCE
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
c      type       ->       see below
c      omega      ->      angular frequency
c      t0,t1      ->      time constant when needed
c      dt      ->      sampling rate
c**********************************************************

      function fsource (type, omega, t0, dt)

      implicit none
      integer  type
      real            pi,pi2,dt,t0
      real*8            uur,uui,trise,trupt
      complex*16      fsource,uu,uex,uxx,omega,shx,ai

      common/par/ai,pi,pi2

 
c TYPE=5        Source = rampe causale
c               rise time T=t0
c        trise=t0
        uu=ai*omega*t0
        uu=(1.-exp(-uu))/uu
        fsource=uu/(ai*omega)
        return
      end
c
c  ===============================================================
c
c        
      subroutine  fft2cd (a,m,iwk)                                       
c                                  specifications for arguments         
      implicit real*8 (a-h,o-z)
      integer            m,iwk(*)                                       
      complex*16      a(*)                                           
c                                  specifications for local variables   
      integer     i,isp,j,jj,jsp,k,k0,k1,k2,k3,kb,
     &              kn,mk,mm,mp,n,n4,n8,n2,lm,nn,jk 
      real*8      rad,c1,c2,c3,s1,s2,s3,ck,sk,sq,a0,a1,a2,a3,    
     &              b0,b1,b2,b3,twopi,temp,                        
     &              zero,one,z0(2),z1(2),z2(2),z3(2)               
      complex*16  za0,za1,za2,za3,ak2                            
      equivalence  (za0,z0(1)),(za1,z1(1)),(za2,z2(1)),           
     &              (za3,z3(1)),(a0,z0(1)),(b0,z0(2)),(a1,z1(1)),  
     &              (b1,z1(2)),(a2,z2(1)),(b2,z2(2)),(a3,z3(1)),   
     &              (b3,z3(2))                                     
      data        sq,sk,ck,twopi/.7071068,.3826834,
     &              .9238795,6.283185/                                
      data        zero/0.0/,one/1.0/                             
c                   sq=sqrt2/2,sk=sin(pi/8),ck=cos(pi/8) 
c                   twopi=2*pi                           
c                                  first executable statement           
      mp = m+1                                                          
      n = 2**m                                                          
      iwk(1) = 1                                                        
      mm = (m/2)*2                                                      
      kn = n+1                                                          
c                                  initialize work vector               
      do 5  i=2,mp                                                      
         iwk(i) = iwk(i-1)+iwk(i-1)                                     
5      continue                                                          
      rad = twopi/n                                                     
      mk = m - 4                                                        
      kb = 1                                                            
      if (mm .eq. m) go to 15                                           
      k2 = kn                                                           
      k0 = iwk(mm+1) + kb                                               
10      k2 = k2 - 1                                                       
      k0 = k0 - 1                                                       
      ak2 = a(k2)                                                       
      a(k2) = a(k0) - ak2                                               
      a(k0) = a(k0) + ak2                                               
      if (k0 .gt. kb) go to 10                                          
15      c1 = one                                                          
      s1 = zero                                                         
      jj = 0                                                            
      k = mm - 1                                                        
      j = 4                                                             
      if (k .ge. 1) go to 30                                            
      go to 70                                                          
20      if (iwk(j) .gt. jj) go to 25                                      
      jj = jj - iwk(j)                                                  
      j = j-1                                                           
      if (iwk(j) .gt. jj) go to 25                                      
      jj = jj - iwk(j)                                                  
      j = j - 1                                                         
      k = k + 2                                                         
      go to 20                                                          
25      jj = iwk(j) + jj                                                  
      j = 4                                                             
30      isp = iwk(k)                                                      
      if (jj .eq. 0) go to 40                                           
c                                  reset trigonometric parameter(s       )
      c2 = jj * isp * rad                                               
      c1 = cos(c2)                                                      
      s1 = sin(c2)                                                      
35      c2 = c1 * c1 - s1 * s1                                            
      s2 = c1 * (s1 + s1)                                               
      c3 = c2 * c1 - s2 * s1                                            
      s3 = c2 * s1 + s2 * c1                                            
40      jsp = isp + kb                                                    
c                                  determine fourier coefficients       
c                                    in groups of 4                     
      do 50 i=1,isp                                                     
         k0 = jsp - i                                                   
         k1 = k0 + isp                                                  
         k2 = k1 + isp                                                  
         k3 = k2 + isp                                                  
         za0 = a(k0)                                                    
         za1 = a(k1)                                                    
         za2 = a(k2)                                                    
         za3 = a(k3)                                                    
         if (s1 .eq. zero) go to 45                                     
         temp = a1                                                      
         a1 = a1 * c1 - b1 * s1                                         
         b1 = temp * s1 + b1 * c1                                       
         temp = a2                                                      
         a2 = a2 * c2 - b2 * s2                                         
         b2 = temp * s2 + b2 * c2                                       
         temp = a3                                                      
         a3 = a3 * c3 - b3 * s3                                         
         b3 = temp * s3 + b3 * c3                                       
45         temp = a0 + a2                                                 
         a2 = a0 - a2                                                   
         a0 = temp                                                      
         temp = a1 + a3                                                 
         a3 = a1 - a3                                                   
         a1 = temp                                                      
         temp = b0 + b2                                                 
         b2 = b0 - b2                                                   
         b0 = temp                                                      
         temp = b1 + b3                                                 
         b3 = b1 - b3                                                   
         b1 = temp                                                      
         a(k0) = cmplx(a0+a1,b0+b1)                                     
         a(k1) = cmplx(a0-a1,b0-b1)                                     
         a(k2) = cmplx(a2-b3,b2+a3)                                     
         a(k3) = cmplx(a2+b3,b2-a3)                                     
50      continue                                                          
      if (k .le. 1) go to 55                                            
      k = k - 2                                                         
      go to 30                                                          
55      kb = k3 + isp                                                     
c                                  check for completion of final        
c                                    iteration                          
      if (kn .le. kb) go to 70                                          
      if (j .ne. 1) go to 60                                            
      k = 3                                                             
      j = mk                                                            
      go to 20                                                          
60      j = j - 1                                                         
      c2 = c1                                                           
      if (j .ne. 2) go to 65                                            
      c1 = c1 * ck + s1 * sk                                            
      s1 = s1 * ck - c2 * sk                                            
      go to 35                                                          
65      c1 = (c1 - s1) * sq                                               
      s1 = (c2 + s1) * sq                                               
      go to 35                                                          
70      continue                                                          
c                                  permute the complex vector in        
c                                    reverse binary order to normal     
c                                    order                              
      if(m .le. 1) go to 9005                                           
      mp = m+1                                                          
      jj = 1                                                            
c                                  initialize work vector               
      iwk(1) = 1                                                        
      do 75  i = 2,mp                                                   
         iwk(i) = iwk(i-1) * 2                                          
75      continue                                                          
      n4 = iwk(mp-2)                                                    
      if (m .gt. 2) n8 = iwk(mp-3)                                      
      n2 = iwk(mp-1)                                                    
      lm = n2                                                           
      nn = iwk(mp)+1                                                    
      mp = mp-4                                                         
c                                  determine indices and switch a       
      j = 2                                                             
80      jk = jj + n2                                                      
      ak2 = a(j)                                                        
      a(j) = a(jk)                                                      
      a(jk) = ak2                                                       
      j = j+1                                                           
      if (jj .gt. n4) go to 85                                          
      jj = jj + n4                                                      
      go to 105                                                         
85      jj = jj - n4                                                      
      if (jj .gt. n8) go to 90                                          
      jj = jj + n8                                                      
      go to 105                                                         
90      jj = jj - n8                                                      
      k = mp                                                            
95      if (iwk(k) .ge. jj) go to 100                                     
      jj = jj - iwk(k)                                                  
      k = k - 1                                                         
      go to 95                                                          
100      jj = iwk(k) + jj                                                  
105      if (jj .le. j) go to 110                                          
      k = nn - j                                                        
      jk = nn - jj                                                      
      ak2 = a(j)                                                        
      a(j) = a(jj)                                                      
      a(jj) = ak2                                                       
      ak2 = a(k)                                                        
      a(k) = a(jk)                                                      
      a(jk) = ak2                                                       
110      j = j + 1                                                         
c                                  cycle repeated until limiting number 
c                                    of changes is achieved             
      if (j .le. lm) go to 80                                           
c                                                                       
9005      return                                                            
      end                                                               
C
C XAPIIR -- SUBROUTINE:   IIR FILTER DESIGN AND IMPLEMENTATION 
C
C  AUTHOR:  Dave Harris
C  LAST MODIFIED:  September 12, 1990C
C  ARGUMENTS:                                                                  C  ----------                                                                   
C                                                                               
C    DATA           REAL ARRAY CONTAINING SEQUENCE TO BE FILTERED               
C                     ORIGINAL DATA DESTROYED, REPLACED BY FILTERED DATA        
C                                                                               
C    NSAMPS         NUMBER OF SAMPLES IN DATA                                   
C                                                                               
C                                                                               
C    APROTO         CHARACTER*8 VARIABLE, CONTAINS TYPE OF ANALOG               
C                     PROTOTYPE FILTER                                          
C                     '(BU)TTER  ' -- BUTTERWORTH FILTER                        
C                     '(BE)SSEL  ' -- BESSEL FILTER                             
C                     'C1      ' -- CHEBYSHEV TYPE I                            
C                     'C2      ' -- CHEBYSHEV TYPE II                           
C                                                                               
C    TRBNDW         TRANSITION BANDWIDTH AS FRACTION OF LOWPASS                 
C                   PROTOTYPE FILTER CUTOFF FREQUENCY.  USED                    
C                   ONLY BY CHEBYSHEV FILTERS.                                  
C                                                                               
C    A              ATTENUATION FACTOR.  EQUALS AMPLITUDE                       
C                   REACHED AT STOPBAND EDGE.  USED ONLY BY                     
C                   CHEBYSHEV FILTERS.                                          
C                                                                               
C    IORD           ORDER (#POLES) OF ANALOG PROTOTYPE                          
C                   NOT TO EXCEED 10 IN THIS CONFIGURATION.  4 - 5              
C                   SHOULD BE AMPLE.                                            
C                                                                               
C    TYPE           CHARACTER*8 VARIABLE CONTAINING FILTER TYPE                 
C                     'LP' -- LOW PASS                                          
C                     'HP' -- HIGH PASS                                         
C                     'BP' -- BAND PASS                                         
C                     'BR' -- BAND REJECT                                       
C                                                                               
C    FLO            LOW FREQUENCY CUTOFF OF FILTER (HERTZ)                      
C                   IGNORED IF TYPE = 'LP'                                      
C                                                                               
C    FHI            HIGH FREQUENCY CUTOFF OF FILTER (HERTZ)                     
C                   IGNORED IF TYPE = 'HP'                                      
C                                                                               
C    TS             SAMPLING INTERVAL (SECONDS)                                 
C                                                                              
C    PASSES           INTEGER VARIABLE CONTAINING THE NUMBER OF PASSES         
C                  1 -- FORWARD FILTERING ONLY                                 
C                   2 -- FORWARD AND REVERSE (I.E. ZERO PHASE) FILTERING       
C                                                                              
C                                                                              
C  SUBPROGRAMS REFERENCED:  BILIN2, BUROOTS, WARP, CUTOFFS, LPTHP, LPTBP,      
C    LP, LPTBR, BEROOTS, C1ROOTS, C2ROOTS, CHEBPARM, DESIGN, APPLY             
C                                                                              
      SUBROUTINE XAPIIR( DATA, NSAMPS, APROTO, TRBNDW, A, IORD,         
     +                   TYPE, FLO, FHI, TS, PASSES )
C    
        DIMENSION DATA(1)                                               
        CHARACTER*8 TYPE, APROTO                                        
        INTEGER NSAMPS, PASSES, IORD                                    
        REAL*4 TRBNDW, A, FLO, FHI, TS, SN(30), SD(30)                  
        LOGICAL ZP                                                      
C                                                                              
C  Filter designed
C
        CALL DESIGN( IORD, TYPE(1:2), APROTO(1:2), A, TRBNDW,           
     &               FLO, FHI, TS, SN, SD, NSECTS )                     
C
C  Filter data                                                                 
C                                                                              
        IF (   PASSES .EQ. 1 ) THEN                                     
          ZP = .FALSE.                                                  
        ELSE                                                            
          ZP = .TRUE.                                                   
        END IF                                                          
        CALL APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS )                  
C
      RETURN                                                            
      END                                                                       
C                                                               APPLY           
C  Subroutine to apply an iir filter to a data sequence.                        
C    The filter is assumed to be stored as second order sections.               
C    Filtering is in-place.                                                     
C    Zero-phase (forward and reverse) is an option.                             
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    DATA                           Array containing data                       
C                                                                               
C    NSAMPS                         Number of data samples                      
C                                                                               
C    ZP                             Logical variable, true for                  
C                                     zero phase filtering, false               
C                                     for single pass filtering                 
C                                                                               
C    SN                             Numerator polynomials for second            
C                                     order sections.                           
C                                                                               
C    SD                             Denominator polynomials for second          
C                                     order sections.                           
C                                                                               
C    NSECTS                         Number of second-order sections             
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    DATA                          Data array (same as input)                   
C                                                                               
C                                                                               
      SUBROUTINE APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS )              
C                                                                               
        REAL*4 SN(1), SD(1), DATA(1)                                    
        REAL*4 OUTPUT                                                   
        LOGICAL ZP                                                      
C                                                                               
        JPTR = 1                                                        
        DO    1 J = 1, NSECTS                                           
C                                                                               
          X1 = 0.0                                                      
          X2 = 0.0                                                      
          Y1 = 0.0                                                      
          Y2 = 0.0                                                      
          B0 = SN(JPTR)                                                 
          B1 = SN(JPTR+1)                                               
          B2 = SN(JPTR+2)                                               
          A1 = SD(JPTR+1)                                               
          A2 = SD(JPTR+2)                                               
C                                                                               
          DO    2 I = 1, NSAMPS                                         
C                                                                               
            OUTPUT = B0*DATA(I) + B1*X1 + B2*X2                         
            OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )                         
            Y2 = Y1                                                     
            Y1 = OUTPUT                                                 
            X2 = X1                                                     
            X1 = DATA(I)                                                
            DATA(I) = OUTPUT                                            
C                                                                               
    2     CONTINUE                                                      
C                                                                               
          JPTR = JPTR + 3                                               
C                                                                               
    1   CONTINUE                                                        
C                                                                               
        IF (   ZP ) THEN                                                
C                                                                               
          JPTR = 1                                                      
          DO    3 J = 1, NSECTS                                         
C                                                                               
            X1 = 0.0                                                    
            X2 = 0.0                                                    
            Y1 = 0.0                                                    
            Y2 = 0.0                                                    
            B0 = SN(JPTR)                                               
            B1 = SN(JPTR+1)                                             
            B2 = SN(JPTR+2)                                             
            A1 = SD(JPTR+1)                                             
            A2 = SD(JPTR+2)                                             
C                                                                               
            DO    4 I = NSAMPS, 1, -1                                   
C                                                                               
              OUTPUT = B0*DATA(I) + B1*X1 + B2*X2                       
              OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )                       
              Y2 = Y1                                                   
              Y1 = OUTPUT                                               
              X2 = X1                                                   
              X1 = DATA(I)                                              
              DATA(I) = OUTPUT                                          
C                                                                               
    4       CONTINUE                                                    
C                                                                               
            JPTR = JPTR + 3                                             
C                                                                               
    3     CONTINUE                                                      
C                                                                               
        END IF                                                          
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C  DESIGN -- Subroutine to design IIR digital filters from analog              
C    prototypes.
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    IORD                Filter order (10 MAXIMUM)                              
C                                                                               
C    TYPE                Character*2 variable containing filter type            
C                          LOWPASS (LP)                                         
C                          HIGHPASS (HP)                                        
C                          BANDPASS (BP)                                        
C                          BANDREJECT (BR)                                      
C                                                                               
C    APROTO              Character*2 variable designating analog prototype      
C                          Butterworth (BU)                                     
C                          Bessel (BE)                                          
C                          Chebyshev Type I (C1)                                
C                          Chebyshev Type II (C2)                               
C                                                                               
C    A                   Chebyshev stopband attenuation factor                  
C                                                                               
C    TRBNDW              Chebyshev transition bandwidth (fraction of            
C                          lowpass prototype passband width)                    
C                                                                               
C    FL                  Low-frequency cutoff                                   
C                                                                               
C    FH                  High-frequency cutoff                                  
C                                                                               
C    TS                  Sampling interval (in seconds)                         
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                  Array containing numerator coefficients of             
C                        second-order sections packed head-to-tail.             
C                                                                               
C    SD                  Array containing denominator coefficients              
C                        of second-order sections packed head-to-tail.          
C                                                                               
C    NSECTS              Number of second-order sections.                       
C                                                                               
C                                                                               
      SUBROUTINE DESIGN( IORD, TYPE, APROTO, A, TRBNDW,                 
     &                   FL, FH, TS, SN, SD, NSECTS )                   
C                                                                               
        COMPLEX P(10), Z(10)                                            
        CHARACTER*2 TYPE, APROTO                                        
        CHARACTER*3 STYPE(10)                                           
        REAL*4 SN(1), SD(1)                                             
C                                                                               
C  Analog prototype selection                                                   
C                                                                               
        IF (     APROTO .EQ. 'BU' ) THEN                                
C                                                                               
          CALL BUROOTS( P, STYPE, DCVALUE, NSECTS, IORD )               
C                                                                               
        ELSE IF (    APROTO .EQ. 'BE' ) THEN                            
C                                                                               
          CALL BEROOTS( P, STYPE, DCVALUE, NSECTS, IORD )               
C                                                                               
        ELSE IF (    APROTO .EQ. 'C1' ) THEN                            
C                                                                               
          CALL CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )                 
          CALL C1ROOTS( P, STYPE, DCVALUE, NSECTS, IORD, EPS )          
C                                                                               
        ELSE IF (    APROTO .EQ. 'C2' ) THEN                            
C                                                                               
          OMEGAR = 1. + TRBNDW                                          
          CALL C2ROOTS( P, Z, STYPE, DCVALUE, NSECTS, IORD, A, OMEGAR ) 
C                                                                               
        END IF                                                          
C                                                                               
C  Analog mapping selection                                                     
C                                                                               
        IF (     TYPE .EQ. 'BP' ) THEN                                  
C                                                                               
          FLW = WARP( FL*TS/2., 2. )                                    
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LPTBP( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )  
C                                                                               
        ELSE IF (   TYPE .EQ. 'BR' ) THEN                               
C                                                                               
          FLW = WARP( FL*TS/2., 2. )                                    
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LPTBR( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )  
C                                                                               
        ELSE IF (    TYPE .EQ. 'LP' ) THEN                              
C                                                                               
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )               
          CALL CUTOFFS( SN, SD, NSECTS, FHW )                           
C                                                                               
        ELSE IF (    TYPE .EQ. 'HP' ) THEN                              
C                                                                               
          FLW = WARP( FL*TS/2., 2. )                                    
          CALL LPTHP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )            
          CALL CUTOFFS( SN, SD, NSECTS, FLW )                           
C                                                                               
        END IF                                                          
C                                                                               
C  Bilinear analog to digital transformation                                    
C                                                                               
        CALL BILIN2( SN, SD, NSECTS )                                   
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C BUROOTS -- SUBROUTINE TO COMPUTE BUTTERWORTH POLES FOR                        
C   NORMALIZED LOWPASS FILTER                                                   
C                                                                               
C LAST MODIFIED:  SEPTEMBER 7, 1990                                             
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C                                                                               
      SUBROUTINE BUROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )             
C                                                                               
        COMPLEX P(1)                                                    
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            
C                                                                               
        PI=3.14159265                                                   
C                                                                               
        HALF = IORD/2                                                   
C                                                                               
C TEST FOR ODD ORDER, AND ADD POLE AT -1                                        
C                                                                               
        NSECTS = 0                                                      
        IF (    2*HALF .LT. IORD ) THEN                                 
          P(1) = CMPLX( -1., 0. )                                       
          RTYPE(1) = 'SP'                                               
          NSECTS = 1                                                    
        END IF                                                          
C       
        DO    1  K = 1, HALF
          ANGLE = PI * ( .5 + FLOAT(2*K-1) / FLOAT(2*IORD) )            
          NSECTS = NSECTS + 1   
          P(NSECTS) = CMPLX( COS(ANGLE), SIN(ANGLE) )                   
          RTYPE(NSECTS) = 'CP'                                          
    1   CONTINUE                                                        
C                                                                               
        DCVALUE = 1.0                                                   
C                                                                               
      RETURN                                                            
      END  
C                                                                               
C BEROOTS -- SUBROUTINE TO RETURN BESSEL POLES FOR                              
C   NORMALIZED LOWPASS FILTER                                                   
C                                                                               
C LAST MODIFIED:  April 15, 1992. Changed P and RTYPE to adjustable 
C                 array by using an "*" rather than a "1".     
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C                                                                               
      SUBROUTINE BEROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )             
C                                                                               
        COMPLEX P(*)                                                    
        INTEGER NSECTS, IORD                                            
        CHARACTER*3 RTYPE(*)                                            
C                                                                               
        IF (   IORD .EQ. 1 ) THEN                                       
C                                                                               
          P(1) = CMPLX( -1.0, 0.0 )                                     
          RTYPE(1) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 2 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -1.1016013,  0.6360098 )                        
          RTYPE(1) = 'CP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 3 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -1.0474091, 0.9992645 )                         
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3226758, 0.0 )                               
          RTYPE(2) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 4 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9952088,  1.2571058 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3700679, 0.4102497 )                         
          RTYPE(2) = 'CP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 5 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9576766,  1.4711244 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3808774,  0.7179096 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.5023160, 0.0 )                               
          RTYPE(3) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 6 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9306565,  1.6618633 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3818581,  0.9714719 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.5714904,  0.3208964 )                        
          RTYPE(3) = 'CP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 7 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9098678,  1.8364514 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3789032,  1.1915667 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.6120388,  0.5892445 )                        
          RTYPE(3) = 'CP'                                               
          P(4) = CMPLX( -1.6843682, 0.0 )                               
          RTYPE(4) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 8 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.8928710,  1.9983286 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3738431,  1.3883585 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.6369417,  0.8227968 )                        
          RTYPE(3) = 'CP'                                               
          P(4) = CMPLX( -1.7574108,  0.2728679 )                        
          RTYPE(4) = 'CP'                                               
C                                                                               
        END IF                                                          
C                                                                               
        NSECTS = IORD - IORD/2                                          
C                                                                               
        DCVALUE = 1.0                                                   
C                                                                               
C  DONE                                                                         
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C  CHEBPARM - Calculates Chebyshev type I and II design parameters              
C                                                                               
C                                                                               
C  INPUT ARGUMENTS                                                              
C  ---------------                                                              
C                                                                               
C       A                Desired stopband attenuation                           
C                          i.e. max stopband amplitude is 1/ATTEN               
C                                                                               
C       TRBNDW           Transition bandwidth between stop and passbands        
C                          as a fraction of the passband width                  
C                                                                               
C       IORD             Filter order (number of poles)                         
C                                                                               
C                                                                               
C  OUTPUT ARGUMENTS                                                             
C  ----------------                                                             
C                                                                               
C       EPS              Chebyshev passband parameter                           
C                                                                               
C       RIPPLE           Passband ripple                                        
C                                                                               
      SUBROUTINE CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )               
                                                                        
          OMEGAR  =  1. + TRBNDW                                        
          ALPHA = ( OMEGAR + SQRT( OMEGAR**2 - 1. ) ) ** IORD           
          G = ( ALPHA**2 + 1. ) / (2.*ALPHA)                            
          EPS = SQRT( A**2 - 1. ) / G                                   
          RIPPLE = 1. / SQRT( 1. + EPS**2 )                             
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C C1ROOTS -- SUBROUTINE TO COMPUTE CHEBYSHEV TYPE I POLES FOR                   
C   NORMALIZED LOWPASS FILTER                                                   
C                                                                               
C LAST MODIFIED:  SEPTEMBER 7, 1990                                             
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        RESPONSE OF FILTER AT ZERO FREQUENCY                      
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C      EPS            CHEBYSHEV PARAMETER RELATED TO PASSBAND RIPPLE            
C                                                                               
      SUBROUTINE C1ROOTS( P, RTYPE, DCVALUE, NSECTS, IORD, EPS )        
C                                                                               
        COMPLEX P(1)                                                    
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            
C                                                                               
        PI = 3.14159265                                                 
        HALF = IORD/2                                                   
C                                                                               
C  INTERMEDIATE DESIGN PARAMETERS                                               
C                                                                               
        GAMMA = ( 1. + SQRT( 1. + EPS*EPS ) ) / EPS                     
        GAMMA = ALOG(GAMMA) / FLOAT(IORD)                               
        GAMMA = EXP(GAMMA)                                              
        S = .5 * ( GAMMA - 1./GAMMA )                                   
        C = .5 * ( GAMMA + 1./GAMMA )                                   
C                                                                               
C  CALCULATE POLES                                                              
C                                                                               
        NSECTS = 0                                                      
        DO    1  I = 1 ,  HALF                                          
          RTYPE(I) = 'CP'                                               
          ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)                       
          SIGMA = -S * SIN(ANGLE)                                       
          OMEGA =  C * COS(ANGLE)                                       
          P(I) = CMPLX( SIGMA, OMEGA )                                  
          NSECTS = NSECTS + 1                                           
    1   CONTINUE                                                        
        IF (   2*HALF .LT. IORD ) THEN                                  
          RTYPE( HALF + 1 ) = 'SP'                                      
          P(HALF+1) = CMPLX( -S, 0.0 )                                  
          NSECTS = NSECTS + 1                                           
          DCVALUE = 1.0                                                 
        ELSE                                                            
          DCVALUE = 1./SQRT( 1 + EPS**2 )                               
        END IF                                                          
C                                                                               
C  DONE                                                                         
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C C2ROOTS -- SUBROUTINE TO COMPUTE ROOTS FOR NORMALIZED LOWPASS                 
C   CHEBYSHEV TYPE 2 FILTER                                                     
C                                                                               
C LAST MODIFIED:  SEPTEMBER 7, 1990                                             
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      Z              COMPLEX ARRAY CONTAINING ZEROS                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL ZEROS                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C      A              STOPBAND ATTENUATION FACTOR                               
C                                                                               
C      OMEGAR         CUTOFF FREQUENCY OF STOPBAND                              
C                     PASSBAND CUTOFF IS AT 1.0 HERTZ                           
C                                                                               
C                                                                               
      SUBROUTINE C2ROOTS( P, Z, RTYPE, DCVALUE, NSECTS, IORD, A, OMEGAR 
     1)                                                                 
C                                                                               
        COMPLEX P(1), Z(1)                                              
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            
C                                                                               
        PI = 3.14159265                                                 
        HALF = IORD/2                                                   
C                                                                               
C  INTERMEDIATE DESIGN PARAMETERS                                               
C                                                                               
        GAMMA = (A+SQRT(A*A-1.))                                        
        GAMMA = ALOG(GAMMA)/FLOAT(IORD)                                 
        GAMMA = EXP(GAMMA)                                              
        S = .5*(GAMMA-1./GAMMA)                                         
        C = .5*(GAMMA+1./GAMMA)                                         
C                                                                               
        NSECTS = 0                                                      
        DO    1 I = 1, HALF                                             
C                                                                               
C  CALCULATE POLES                                                              
C                                                                               
          RTYPE(I) = 'CPZ'                                              
C                                                                               
          ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)                       
          ALPHA = -S*SIN(ANGLE)                                         
          BETA = C*COS(ANGLE)                                           
          DENOM = ALPHA*ALPHA + BETA*BETA                               
          SIGMA = OMEGAR*ALPHA/DENOM                                    
          OMEGA = -OMEGAR*BETA/DENOM                                    
          P(I) = CMPLX( SIGMA, OMEGA )                                  
C                                                                               
C  CALCULATE ZEROS                                                              
C                                                                               
          OMEGA = OMEGAR/COS(ANGLE)                                     
          Z(I) = CMPLX( 0.0, OMEGA )                                    
C                                                                               
          NSECTS = NSECTS + 1                                           
C                                                                               
    1   CONTINUE                                                        
C                                                                               
C  ODD-ORDER FILTERS                                                            
C                                                                               
        IF (  2*HALF .LT. IORD ) THEN                                   
          RTYPE(HALF+1) = 'SP'                                          
          P(HALF+1) = CMPLX( -OMEGAR/S, 0.0 )                           
          NSECTS = NSECTS + 1                                           
        END IF                                                          
C                                                                               
C  DC VALUE                                                                     
C                                                                               
        DCVALUE = 1.0                                                   
C                                                                               
C  DONE                                                                         
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C WARP -- FUNCTION, APPLIES TANGENT FREQUENCY WARPING TO COMPENSATE             
C         FOR BILINEAR ANALOG -> DIGITAL TRANSFORMATION                         
C                                                                               
C ARGUMENTS:                                                                    
C ----------                                                                    
C                                                                               
C      F       ORIGINAL DESIGN FREQUENCY SPECIFICATION (HERTZ)                  
C      TS      SAMPLING INTERVAL (SECONDS)                                      
C                                                                               
C  LAST MODIFIED:  SEPTEMBER 20, 1990                                           
C                                                                               
      REAL FUNCTION WARP( F , TS )                                      
C                                                                               
        TWOPI = 6.2831853                                               
        ANGLE = TWOPI*F*TS/2.                                           
        WARP = 2.*TAN(ANGLE)/TS                                         
        WARP = WARP/TWOPI                                               
C                                                                               
      RETURN                                                            
      END      
C                                                                               
C  Subroutine to generate second order section parameterization                 
C    from an pole-zero description for lowpass filters.                         
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeros                             
C                                                                               
C    RTYPE                   Character array containing root type information   
C                              (SP)  Single real pole or                        
C                              (CP)  Complex conjugate pole pair                
C                              (CPZ) Complex conjugate pole and zero pairs      
C                                                                               
C    DCVALUE                 Zero-frequency value of prototype filter           
C                                                                               
C    NSECTS                  Number of second-order sections                    
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C                                                                               
      SUBROUTINE LP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )             
C                                                                               
        COMPLEX P(*), Z(*)                                              
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          IF (   RTYPE(I) .EQ. 'CPZ' ) THEN                             
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
     &            / REAL( Z(I) * CONJG( Z(I) ) )                        
            SN( IPTR )     = REAL( Z(I) * CONJG( Z(I) ) ) * SCALE       
            SN( IPTR + 1 ) = -2. * REAL( Z(I) ) * SCALE                 
            SN( IPTR + 2 ) = 1. * SCALE                                 
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )               
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
            SN( IPTR )     = SCALE                                      
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )               
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SCALE = -REAL( P(I) )                                       
            SN( IPTR )     = SCALE                                      
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = -REAL( P(I) )                              
            SD( IPTR + 1 ) = 1.                                         
            SD( IPTR + 2 ) = 0.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
        SN(1) = DCVALUE * SN(1)                                         
        SN(2) = DCVALUE * SN(2)                                         
        SN(3) = DCVALUE * SN(3)                                         
                                                                        
C                                                                               
      RETURN                                                            
      END
C                                                    LPTBP                      
C                                                                               
C  Subroutine to convert an prototype lowpass filter to a bandpass filter via   
C    the analog polynomial transformation.  The lowpass filter is               
C    described in terms of its poles and zeros (as input to this routine).      
C    The output consists of the parameters for second order sections.           
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeros                             
C                                                                               
C    RTYPE                   Character array containing type information        
C                              (SP) single real pole  or                        
C                              (CP) complex conjugate pole pair  or             
C                              (CPZ) complex conjugate pole/zero pairs          
C                                                                               
C    DCVALUE                 Zero frequency value of filter                     
C                                                                               
C    NSECTS                  Number of second-order sections upon input         
C                                                                               
C    FL                      Low-frequency cutoff                               
C                                                                               
C    FH                      High-frequency cutoff                              
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C    NSECTS                  Number of second order sections upon output        
C                              This subroutine doubles the number of            
C                              sections.                                        
C                                                                               
C                                                                               
      SUBROUTINE LPTBP( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )  
C                                                                               
        COMPLEX P(*), Z(*), CTEMP, P1, P2, Z1, Z2, S, H                 
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    
C                                                                               
        PI = 3.14159265                                                 
        TWOPI = 2.*PI                                                   
        A = TWOPI*TWOPI*FL*FH                                           
        B = TWOPI*( FH - FL )                                           
C                                                                               
        N = NSECTS                                                      
        NSECTS = 0                                                      
        IPTR = 1                                                        
        DO    1 I = 1, N                                                
C                                                                               
          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN                            
C                                                                               
            CTEMP = ( B*Z(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            Z1 = 0.5*( B*Z(I) + CTEMP )                                 
            Z2 = 0.5*( B*Z(I) - CTEMP )                                 
            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z1 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z2 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = A                                          
            SD( IPTR + 1 ) = -B*REAL( P(I) )                            
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 1                                         
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
C  Scaling - use the fact that the bandpass filter amplitude at sqrt( omega_l * 
C            equals the amplitude of the lowpass prototype at d.c.              
C                                                                               
        S = CMPLX( 0., SQRT(A) )                                        
        H = CMPLX( 1., 0. )                                             
C                                                                               
        IPTR = 1                                                        
        DO    2 I = 1, NSECTS                                           
          H = H * ( ( SN(IPTR+2)*S + SN(IPTR+1) )*S + SN(IPTR) )        
     &      / ( ( SD(IPTR+2)*S + SD(IPTR+1) )*S + SD(IPTR) )            
          IPTR = IPTR + 3                                               
    2   CONTINUE                                                        
        SCALE = DCVALUE / SQRT( REAL( H ) * CONJG( H ) )                
                                                                        
        SN(1) = SN(1) * SCALE                                           
        SN(2) = SN(2) * SCALE                                           
        SN(3) = SN(3) * SCALE                                           
C                                                                               
      RETURN                                                            
      END                                                               
C                                                    LPTBR                      
C                                                                               
C  Subroutine to convert a lowpass filter to a band reject filter               
C    via an analog polynomial transformation.  The lowpass filter is            
C    described in terms of its poles and zeros (as input to this routine).      
C    The output consists of the parameters for second order sections.           
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeros                             
C                                                                               
C    RTYPE                   Character array containing type information        
C                              (SP)  single real pole or                        
C                              (CP)  complex conjugate pole pair                
C                              (CPZ) complex conjugate pole/zero pairs          
C                                                                               
C    DCVALUE                 Zero-frequency value of prototype filter           
C                                                                               
C    NSECTS                  Number of second-order sections                    
C                              prior to transformation                          
C                                                                               
C    FL                      Low-frequency cutoff                               
C                                                                               
C    FH                      High-frequency cutoff                              
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C    NSECTS                  Number of second order sections following          
C                              transformation.  The number is doubled.          
C                                                                               
C                                                                               
      SUBROUTINE LPTBR( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )  
C                                                                               
        COMPLEX P(*), Z(*), CINV, CTEMP, P1, P2, Z1, Z2                 
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*)                                             
C                                                                               
        PI = 3.14159265                                                 
        TWOPI = 2.*PI                                                   
        A = TWOPI*TWOPI*FL*FH                                           
        B = TWOPI*( FH - FL )                                           
C                                                                               
        N = NSECTS                                                      
        NSECTS = 0                                                      
        IPTR = 1                                                        
        DO    1 I = 1, N                                                
C                                                                               
          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN                            
C                                                                               
            CINV = 1./Z(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            Z1 = 0.5*( B*CINV + CTEMP )                                 
            Z2 = 0.5*( B*CINV - CTEMP )                                 
            CINV = 1./P(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*CINV + CTEMP )                                 
            P2 = 0.5*( B*CINV - CTEMP )                                 
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z1 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z2 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            CINV = 1./P(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*CINV + CTEMP )                                 
            P2 = 0.5*( B*CINV - CTEMP )                                 
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = -A*REAL( P(I) )                            
            SD( IPTR + 1 ) = B                                          
            SD( IPTR + 2 ) = -REAL( P(I) )                              
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 1                                         
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
C  Scaling - use the fact that the bandreject filter amplitude  at d.c.         
C            equals the lowpass prototype amplitude at d.c.                     
C                                                                               
        H = 1.0                                                         
C                                                                               
        IPTR = 1                                                        
        DO    2 I = 1, NSECTS                                           
          H = H * SN(IPTR) / SD(IPTR)                                   
          IPTR = IPTR + 3                                               
    2   CONTINUE                                                        
        SCALE = DCVALUE / ABS(H)                                        
        SN(1) = SN(1) * SCALE                                           
        SN(2) = SN(2) * SCALE                                           
        SN(3) = SN(3) * SCALE                                           
C                                                                               
      RETURN                                                            
      END                                                               
C                                                    LPTHP                      
C                                                                               
C  Subroutine to convert a lowpass filter to a highpass filter via              
C    an analog polynomial transformation.  The lowpass filter is                
C    described in terms of its poles and zeroes (as input to this routine).     
C    The output consists of the parameters for second order sections.           
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeroes                            
C                                                                               
C    RTYPE                   Character array containing root type information   
C                              (SP) single real pole or                         
C                              (CP)  complex conjugate pair                     
C                              (CPZ) complex pole/zero pairs                    
C                                                                               
C    DCVALUE                 Zero-frequency value of prototype filter           
C                                                                               
C    NSECTS                  Number of second-order sections                    
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C                                                                               
      SUBROUTINE LPTHP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )          
C                                                                               
        COMPLEX P(*), Z(*)                                              
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          IF (     RTYPE(I) .EQ. 'CPZ' ) THEN                           
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
     &            / REAL( Z(I) * CONJG( Z(I) ) )                        
            SN( IPTR )     = 1.  *  SCALE                               
            SN( IPTR + 1 ) = -2. * REAL( Z(I) )  *  SCALE               
            SN( IPTR + 2 ) = REAL( Z(I) * CONJG( Z(I) ) )  *  SCALE     
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )               
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = SCALE                                      
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )               
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SCALE = -REAL( P(I) )                                       
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = SCALE                                      
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -REAL( P(I) )                              
            SD( IPTR + 2 ) = 0.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
        SN(1) = SN(1) * DCVALUE                                         
        SN(2) = SN(2) * DCVALUE                                         
        SN(3) = SN(3) * DCVALUE                                         
C                                                                               
      RETURN                                                            
      END             
C                                                    CUTOFFS                    
C                                                                               
C  Subroutine to alter the cutoff of a filter.  Assumes that the                
C    filter is structured as second order sections.  Changes                    
C    the cutoffs of normalized lowpass and highpass filters through             
C    a simple polynomial transformation.                                        
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    F                       New cutoff frequency                               
C                                                                               
C  Input/Output Arguments:                                                      
C  -----------------------                                                      
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C    NSECTS                  Number of second order sectionsects                
C                                                                               
C                                                                               
      SUBROUTINE CUTOFFS( SN, SD, NSECTS, F )                           
C                                                                               
        REAL*4 SN(1), SD(1)                                             
C                                                                               
        SCALE = 2.*3.14159265*F                                         
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          SN( IPTR + 1 ) = SN( IPTR + 1 ) / SCALE                       
          SN( IPTR + 2 ) = SN( IPTR + 2 ) / (SCALE*SCALE)               
          SD( IPTR + 1 ) = SD( IPTR + 1 ) / SCALE                       
          SD( IPTR + 2 ) = SD( IPTR + 2 ) / (SCALE*SCALE)               
          IPTR = IPTR + 3                                               
C                                                                               
    1   CONTINUE                                                        
C                                                                               
      RETURN                                                            
      END                                                               

C                                                                               
C  Transforms an analog filter to a digital filter via the bilinear transformati
C    Assumes both are stored as second order sections.  The transformation is   
C    done in-place.                                                             
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    SN                   Array containing numerator polynomial coefficients for
C                           second order sections.  Packed head-to-tail.        
C                                                                               
C    SD                   Array containing denominator polynomial coefficients f
C                           second order sections.  Packed head-to-tail.        
C                                                                               
C    NSECTS               Number of second order sections.                      
C                                                                               
C                                                                               
      SUBROUTINE BILIN2( SN, SD, NSECTS )                               
C                                                                               
        REAL*4 SN(1), SD(1)                                             
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          A0 = SD(IPTR)                                                 
          A1 = SD(IPTR+1)                                               
          A2 = SD(IPTR+2)                                               
C                                                                               
          SCALE = A2 + A1 + A0                                          
          SD(IPTR)   = 1.                                               
          SD(IPTR+1) = (2.*(A0 - A2)) / SCALE                           
          SD(IPTR+2) = (A2 - A1 + A0) / SCALE                           
C                                                                               
          A0 = SN(IPTR)                                                 
          A1 = SN(IPTR+1)                                               
          A2 = SN(IPTR+2)                                               
C                                                                               
          SN(IPTR)   = (A2 + A1 + A0) / SCALE                           
          SN(IPTR+1) = (2.*(A0 - A2)) / SCALE                           
          SN(IPTR+2) = (A2 - A1 + A0) / SCALE                           
C                                                                               
          IPTR = IPTR + 3                                               
C                                                                               
    1   CONTINUE                                                        
C                                                                               
      RETURN                                                            
      END                                                               


      
