c     @(#) convm.F       CONVM 1.19      12/8/93 
c******************************************************************************
c     *			PROGRAMME CONVM                                       *
c     *									      *
c     *                 Version du 15/11/09 
c     *	                                                                      *
c******************************************************************************
c     MODIF C.Francois-Holden + R. Madariaga  2009 + S. Di Carli 2007
c

      subroutine seismograms(ntfd3,dt,nskip,ioutput,iproc)
      
      parameter   (nsp=5000,ncp=33,nrp=19,ntp=512,
     *  nrtp=nrp*ntp, fr2ori=0.2,fr1ori=0.02, bp_order=4)
c
      parameter (nxx=200,nyy=200)
      include 'fd3d_param.inc'
c

      real*8         t_tot
      integer tt,ll,kk,ioutput
      real*8  slip(nxx,nyy),temp(nxx,nxx),slipr(nxx/4,nyy/4,ntp),
     *    temp1(512)
c
      integer   icc,ics,if,ic

      real tax(ntp),tay(ntp),taz(ntp),mask(nxx/4,nyy/4)
      character    filename*25
      character   sourcefile*10,statfile*10
      integer     iwk(ntp),jf,ir,is,it,nc,ns,nr,nfreq,ikmax,mm,nt,
     1     openxdr,wrxdr,io,isc(nsp)

      real  hc(ncp),vp(ncp),vs(ncp),rho(ncp),delay(nsp),
     1     xr(nrp),yr(nrp),zr(nrp),a(6,nsp),qp(ncp),qs(ncp),
     2     tl,xl,uconv,mu(nsp),strike(nsp),dip(nsp),
     3     rake(nsp),disp(nsp),hh,zsc,dfreq,freq,rvel(nsp),
     4     pi,pi2,xs(nsp),ys(nsp),zs(nsp),aw,ck,xmm,cc,
     5     uout(ntp),pas,width(nsp),length(nsp),xref,yref,
     7     t1,momo(nsp),ureal(ntp),fr1,fr2

      complex*16  ufft(ntp,nsp),ux(ntp,nrp),uy(ntp,nrp),uz(ntp,nrp),
     1     omega,uxf(6),uyf(6),uzf(6),ai,us,
     2     uux,uuy,uuz,ucmplx(2*ntp)
      common /pred_com/predicted_data(maxdata),iacc
     &                  ndata,
     &                  lu_obs, lu_dc, fname,
     &                  nummod,num_max_mod
c      common /taille/a(2,maxmoddim),b(2,maxmoddim)
      common/axitra/ nc,nfreq,tl,aw,nr,ns,xl,t0,
     &     ikmax,uconv,sourcefile,statfile,fr1,fr2

      integer            ndata

      complex uuxf(ntp/4,nsp,nrp,6),uuyf(ntp/4,nsp,nrp,6),
     $     uuzf(ntp/4,nsp,nrp,6) 
      common/hist/uuxf,uuyf,uuzf

      common /fd3d_com/nt,nxt,nyt,nz,dh,ndat,obs
      common      /par/ai,pi,pi2
      namelist    /input/ nc,nfreq,tl,aw,nr,ns,xl,t0,
     &     ikmax,uconv,sourcefile,statfile,fr1,fr2
      data  ux,uy,uz/nrtp*(0.,0.),nrtp*(0.,0.),
     &     nrtp*(0.,0.)/

      pi=3.1415926535
      pi2=2.*pi
      ai=(0.,1.)


cccccccccccccccc OJO cccccccccccccccccccc


c      fr1=0.02
c      fr2=0.12

c
c  read axitra header     
c
      open (10,form='formatted',file='axi.head')
      read(10,input)
c      dt0=0.
c    =================================================
c         parameters of numerical simulation
c    ================================================= 
c     nftd3  duration of numerical simulation
c     nskip  sortie de la vitesse tous les combien pas de temps
c     dt teps de la TdF
c     
c      
c++++++++++++
c   Determination de la longueur et du pas de temps 
c++++++++++++
c
      xmm=log(real(nfreq))/log(2.)
      mm=int(xmm+0.05)+2
      nt=2**mm
      dfreq=1./tl
      pas=tl/nt
      aw=-pi*aw/tl 
! PLOMBE
      nttot=ntfd3/nskip
      dt=dt*float(nskip)
c      if(dt.ne.pas) print *,'ATTENTION DT simulation dt, pas = ',dt,pas
      t_tot=nttot*dt
c      write(*,*) 'Time of dynamic simul ntfd3, nttot,ttot,nskip,dt = ',
c     *        ntfd3,nttot,t_tot,nskip,dt
c|      read(*,*) t_tot
c      write(*,*) ' aplique mask ? 1 (oui) or 0(no)'

      open (13,form='formatted',file=sourcefile)
      open (14,form='formatted',file=statfile)
      open (15,form='formatted',file='axi.hist')
      
c     
      
c++++++++++++
c     Milieu, stations et sources
c++++++++++++
      read(10,*)
      do  ic=1,nc
 3       read(10,*) hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)
      enddo
      read(10,*)
      do is=1,ns
         read(13,*) index,xs(is),ys(is),zs(is)
         indexin=-1
         do while (indexin.ne.index)
            read(15,*) indexin,disp(is),strike(is),dip(is),
     1           rake(is),width(is),length(is),delay(is)
         enddo
c         delay(is)=delay(is)+dt0
      enddo
c      write(*,*) 'Numero fuentes, nxt, nyt = ',ns, nxt, nyt

      do ir=1,nr
         read(14,*) xr(ir),yr(ir),zr(ir)
      enddo
      close(13)
      close(14)
      close(15)
      close(10)
c      
c++++++++++++
c    Calcul du tenseur de moments
c++++++++++++
c     convertion depth -> thickness
c
      if (hc(1).eq.0.) then
         do ic=1,nc-1
            hc(ic)=hc(ic+1)-hc(ic)
         enddo
      endif
c
c Boucle sur les sources
c      
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
         mu(is)=vs(isc(is))*vs(isc(is))*rho(isc(is))
         
         call cmoment (mu(is),strike(is),dip(is),
     &        rake(is),disp(is),width(is)*length(is),a(1,is))
      enddo
c
c++++++++++++
c     Lecture de la distribution de sliprate
c++++++++++++
c
c  Decimation
c
c     
      idec=4
      nx=nxt/idec
      ny=nyt/idec
      nxlee=nxt
      nylee=nyt                  
      if(ioutput.eq.1)then 
          open (52,form='formatted',file='result/sliprate.res')
c         print *, 'leyo sliprate.res'
      else
          ip=mod(iproc,10)
          write(filename,1001)ip
1001      format('/tmp/sliprateruiz3',i1,'.res')
c          write(*,*)'seismograms open file ',filename
  
          open (52,form='formatted',file=filename)
      endif
c
c
c
c      dt=0.2
      do tt=1,nt
         do kk=1,nx
            do jj=1,ny
               slipr(kk,jj,tt)=0.
            enddo
         enddo
      enddo      
      do jf=1,nt
         do ir=1,nr
                ux(jf,ir)=0
                uy(jf,ir)=0
                uz(jf,ir)=0
         enddo
      enddo  
c
c  Loop over time for reading sliprate
c    sliprate is defined over ntt points at steps of dt=0.2
c
      slipmx=0.

      do tt=1,nttot
         do ii=1,nxlee
            do jj=1,nylee
               read(52,*) slip(ii,jj)
               if(slipmx.lt.slip(ii,jj))slipmx=slip(ii,jj)
            enddo               ! enddo jj
         enddo                  ! enddo ll
c
c  decimation along x
c
         do kk=1,nylee
            do ii=1,nx
               ll=idec*(ii-1)
               temp(ii,kk)=0.
               do jj=1,idec
                  temp(ii,kk)=temp(ii,kk)+slip(ll+jj,kk)/float(idec)
               enddo 
            enddo       
         enddo
         
c
c  decimation along y
c
         do kk=1,nx
            do jj=1,ny
               ll=idec*(jj-1);
               slipr(kk,jj,tt)=0.
               do ii=1,idec
                  slipr(kk,jj,tt)=slipr(kk,jj,tt)
     $              +temp(kk,ll+ii)/float(idec)
               enddo 
            enddo
         enddo
      enddo                     ! enddo tt
      close(52)
c      print *, 'finish reading sliprate, slipmx =',slipmx
c
c  FFT en temps pour chaque element de failee ii, jj
c La ttf est numerote suivant s=1,ns=nx,ny
c
      do ii=1,nx
        do jj=1,ny
           is=(ii-1)*nx+jj
           do it=1,nt
              ck=float(it-1)/nt
              cc=exp(aw*tl*ck)*tl
              ucmplx(it)=cmplx(slipr(ii,jj,it),0)*cc
           enddo
         
           call fft2cd(ucmplx,mm,iwk)
         
           do it=1,nt
              ucmplx(it)=conjg(ucmplx(it))/nt
              ufft(it,is)=ucmplx(it)
           enddo                  !enddo nt
         enddo 
      enddo                     !enddo ns
c
c  Boucle sur les fréquences pour la convolution
c      
      do jf=1,nfreq
         freq=float(jf-1)/tl
         omega=cmplx(pi2*freq,aw)
c
c sources
c
         do is=1,ns
            us= ufft(jf,is)
c
c Fonction de Green pour chaque récepteur
c            
            do ir=1,nr
               uux=0.
               uuy=0.
               uuz=0.
               do it=1,6
                  uux = uux +  uuxf(jf,is,ir,it)*a(it,is)
                  uuy = uuy +  uuyf(jf,is,ir,it)*a(it,is)
                  uuz = uuz +  uuzf(jf,is,ir,it)*a(it,is)
               enddo
c              if(jf<2.and.is<10)print *,
c     *                   'primer bucle is, ir, uuxf(jf,is,ir,it) =',
c     *                     jf,is,ir,uuxf(jf,is,ir,1)
               ux(jf,ir)=ux(jf,ir) + uux * us 
               uy(jf,ir)=uy(jf,ir) + uuy * us
               uz(jf,ir)=uz(jf,ir) + uuz * us
            enddo               !fin boucle recepteur
         enddo                  !fin boucle source
 
c Boucle pour avoir les synth en deplacement   
          amax=0.         
c          do ir=1,nr
c               ux(jf,ir)=ux(jf,ir)/(ai * omega)
c               uy(jf,ir)=uy(jf,ir)/(ai * omega)
c               uz(jf,ir)=uz(jf,ir)/(ai * omega)
c               if(abs(ux(jf,ir)).gt.amax)amax=abs(ux(jf,ir))
c          enddo
c          if(jf<2)print *,'primer bucle jf, amax =',jf, amax
      enddo   
      
c
c++++++++++++
c     Calcul des sismogrammes par FFT inverse
c++++++++++++
c

      if(ioutput.eq.1)then       
        open (62,form='formatted',file='dyn_disp_x')
        open (63,form='formatted',file='dyn_disp_y')
        open (64,form='formatted',file='dyn_disp_z')
        open (83,form='formatted',file='synth.res')
c        print *, 'leyo synth.res'
      else
        ip=mod(iproc,10)
        write(filename,1000)ip
1000    format('/tmp/synth',i1,'.res')
c        write(*,*)'open synth file iproc = ',ip, filename
c        open (83,form='formatted',file=filename)
      endif
c
c  kd indice sur les observes
c
c      print *,'just before it fails nr, nfreq,nt =',nr, nfreq,nt 
c      print *,'fr1, fr2 = ',fr1,fr2
      kd=0
      do 30 ir=1,nr
         
c     on complete le spectre pour les hautes frequences
c     avec inversion du signe de la partie imaginaire pour
c     la FFT inverse qui n'existe pas avec fft2cd
         
         do ifr=nt+2-nfreq,nt
            ux(ifr,ir)=dconjg(ux(nt+2-ifr,ir))
            uy(ifr,ir)=dconjg(uy(nt+2-ifr,ir))
            uz(ifr,ir)=dconjg(uz(nt+2-ifr,ir))
         enddo
         
         call fft2cd(ux(1,ir),mm,iwk)
         call fft2cd(uy(1,ir),mm,iwk)
         call fft2cd(uz(1,ir),mm,iwk)
         
         
c         write(*,*) 'aw end = ', aw
         do it=1,nt
            ck=float(it-1)/nt
            cc=exp(-aw*tl*ck)/tl
            ux(it,ir)=ux(it,ir)*dble(cc)
            uy(it,ir)=uy(it,ir)*dble(cc)
            uz(it,ir)=uz(it,ir)*dble(cc)
         enddo
ccc
         do it=1,nt
            tax(it)=real(real(ux(it,ir)))
            tay(it)=real(real(uy(it,ir)))
            taz(it)=real(real(uz(it,ir)))
         enddo
         
c     20/11/04               
         
         call XAPIIR(taz,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,dt,1)
         call XAPIIR(tay,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,dt,1)
         call XAPIIR(tax,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,dt,1)
         if(ioutput.eq.1) then 
            do iter=1,nt
               write(64,*) taz(iter)
               write(62,*) tax(iter)
               write(63,*) tay(iter)
            enddo
         endif
         do iter=1,nt
            kd=kd+1
            predicted_data(kd)=tax(iter)
         enddo
         do iter=1,nt
            kd=kd+1
            predicted_data(kd)=tay(iter)
         enddo
         do iter=1,nt
            kd=kd+1
            predicted_data(kd)=taz(iter)
         enddo
         
 30   continue
      if(ioutput.eq.1)then       
        do kkd=1,kd
            write(83,*)predicted_data(kkd)
c            print *, predicted_data(kkd)
        enddo
      endif

      close(83)
c      print *,'wrote ',kd,'time points'      
      if(ioutput.eq.1)then
        close(62)
        close(63)
        close(64)
      endif
      return
      end
c
c   =================================================      
c     cmoment
c   =================================================      
c
      subroutine cmoment (mu,strike,dip,
     1     rake,disp,surf,a)
      
      implicit none
      
      real xmoment,mu,strike,dip,rake,disp,surf,a(6),pi2,
     1     sd,cd,sp,cp,s2p,s2d,c2p,c2d,x1,x2,x3,x4,x5,x6,cl,sl,pi
      complex*16 ai
      common      /par/ai,pi,pi2
      
      if (surf.eq.0.) then
         xmoment=disp
      else
c         xmoment=mu*disp*surf
         xmoment=mu*surf
         
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
c     If surf=-1. it is an explosion
      if (surf.eq.-1.) then
         xmoment=disp
         a(1)=0.
         a(2)=0.
         a(3)=0.
         a(4)=0.
         a(5)=0.
         a(6)=xmoment
c     If surf=-2. it is a tensile crack whose normal 
c     should be given
      else if (surf.eq.-2.) then
      endif   
      return
      end
c   
c ========================================
c   READ AXIRES
c ========================================
c

      subroutine readaxires

      parameter   (nsp=5000,ncp=33,nrp=19,ntp=512,
     *     nrtp=nrp*ntp, fr2ori=0.2,fr1ori=0.02, bp_order=4)

c      include 'kine_param.inc'

      character   headfile*50,sourcefile*10,statfile*10

      complex*16 uxf(6),uyf(6),uzf(6),ai
      common/axitra/ nc,nfreq,tl,aw,nr,ns,xl,t0,
     &     ikmax,uconv,sourcefile,statfile,fr1,fr2
      namelist    /input/ nc,nfreq,tl,aw,nr,ns,xl,t0,
     &     ikmax,uconv,sourcefile,statfile,fr1,fr2
      
      complex uuxf(ntp/4,nsp,nrp,6),uuyf(ntp/4,nsp,nrp,6),
     $     uuzf(ntp/4,nsp,nrp,6) 
      common/hist/uuxf,uuyf,uuzf
c     
c     Lecture de axi.head 
c     
      open (20,form='formatted',file='axi.head')

      read(20,input)
c      read(20,*)nc
c      read(20,*)nfreq
c      read(20,*)tl
c      read(20,*)aw
c      read(20,*)nr
c      read(20,*)ns
c      read(20,*)xl
c      read(20,*)t0
c      read(20,*)fr1
c      read(20,*)fr2
c      read(20,*)ikmax
c      read(20,*)uconv
c      read(20,*)sourcefile
c      read(20,*)statfile
c      if ((ns.gt.nsp).or.(nr.gt.nrp).or.(nc.gt.ncp)) then
c         write(0,*) "parametres sous-dimensionnes"
c         stop
c      endif
c      print *,'axi.head :'
c      print  *,'nc,nfreq,tl,aw,nr = ',nc,nfreq,tl,aw,nr
c      print  *,'xl,ikmax,uconv = ',xl,ikmax,uconv
c      print  *,'files = ',sourcefile,statfile
      close(20)
c     
c     ========================================
c     read axi.res
c     ========================================
c     
      open (12,form='unformatted',file='axi.res')
c     print *, 'read axi.res with nfreq,ns,nr = ',nfreq,ns,nr
      amax=0.
      amin =1.e+20
      do jf=1,nfreq
         do is=1,ns
            do ir=1,nr
c     
c     FAUT lire par block de six complexes a cause du format 
c     binaire fortran
c     
               read(12,end=200)(uxf(it),it=1,6)
               read(12,end=200)(uyf(it),it=1,6)
               read(12,end=200)(uzf(it),it=1,6)
               do it=1,6
                  if(amin.gt.abs(uxf(it)))amin=abs(uxf(it))
                  if(amax.lt.abs(uxf(it)))amax=abs(uxf(it))
                  uuxf(jf,is,ir,it)=uxf(it)
                  uuyf(jf,is,ir,it)=uyf(it)
                  uuzf(jf,is,ir,it)=uzf(it)
               enddo
            enddo
         enddo
      enddo   
      goto 201
 200  print *,'Merde, axi.res is not big enough'
      nfreq=jf-1
      
 201  continue
c      print *, 'leyo axi.res  amin, amax =',amin, amax
c      print *,' read nfreq = ',nfreq,' axitra lines' 
      close(12)

      return
      end
