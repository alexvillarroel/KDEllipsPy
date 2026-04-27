c
c   =====================================================
c    main routine to do fd     
c    modified spetember 2011 for dip slip
c    friction law is applied on sigma_{yz}
c  =====================================================
c
     
      subroutine indyna3d(nt,dt,ntiskp,ioutput,iproc,dmax,amom)
      parameter (ntt=10000,io=20,nl=2)

      include 'parstat'

      dimension vpe(nl),vse(nl),den(nl),src(ntt),src1(ntt),
     +     nbhx(100),nbhy(100),nbhz(100)
      integer broken(nx,ny),incrack(nx,ny)
      real ruptime(nx,ny),slip(nx,ny)
      real*4 amom
!
      common /fd3d_com/ntfd3d,nxt,nyt,nzt,dh,ndat,obs
c      print *, 'in indyna3d ioutput = ',ioutput
      nxt=nx
      nyt=ny
      nzt=nz
c      ioutput = 1
      call open3d(1,ioutput,iproc)
c
ccc   read parameters
c
      ifre=0
      iout=1
      call rddt3d(nt,nxt,nyt,nzt,dh,dt,domp,str,dip,rake,slip,isrc,
     *     ainc,az,
     +     nbgx,nedx,nskpx,nbgy,nedy,nskpy,
     +     ntiskp,ndamp,nsv,ifre,nxsc,nysc,nzsc,nholes,nbhx,nbhy,
     +     nbhz,idebug,ibc, dmax_old,vmax)

      call rdmd3d(vpe,vse,den,nxt,nyt,nzt,ifre)
    
c
c     print them out
c
c      print *, '  FD3D du 14 FEVRIER 2008'
c      print *
c      print *,' grid size :  nxt,nyt,nzt = ', nxt,nyt,nzt
c      print *,' center of fault: nxrsc,nysc,nzsc = ', nxsc,nysc,nzsc
c      print *,' debug = ',idebug
c      print *,' dmax = ',dmax
c      print *,' vmax = ',vmax
c
c   initialize the different arrays 
c
      do i=1,nxt
         do j=1,nyt
            broken(i,j)=0
            incrack(i,j)=0.
            ruptime(i,j)=0.
            slip(i,j)=0.
            do k=1,nzt
               u1(i,j,k)=0.
               v1(i,j,k)=0.
               w1(i,j,k)=0.
               xx(i,j,k)=0.
               yy(i,j,k)=0.
               zz(i,j,k)=0.
               xy(i,j,k)=0.
               yz(i,j,k)=0.
               xz(i,j,k)=0.
             enddo
         enddo
      enddo

      ntbegin=1 
      cfl=vpe(2)*dt/dh
c      print *,' dt, dh, cfl = ',dt,dh,cfl
c     
c     
c      call wrpm3d(nxt,nyt,nzt,nt,dh,dt,nsou,
c     +     domp,ndamp,ntiskp,nbgx,nedx,nskpx,nbgy,nedy,
c     +     nskpy,vpe,vse,den,str,dip,rake,slip,
c     +     ifre,axx,ayy,azz,axy,axz,ayz,nxsc,nysc,nzsc,
c     +     peak,ruptur,rupmax, strout, strasp,itype,
c     +     dmax,vmax, xmin,xmax,ymin,ymax )
c     
c     -------------------------------------------------------
c     Loop over time 
c     -------------------------------------------------------
c
      do 10 ntst = ntbegin,nt
         time =ntst*dt
c         if(mod(ntst,100).eq.1)
c     $     print *,' ntime step, time = ', ntst,time  
c
c  -------------------------------------------------------------
c    Apply 4th-order differencing to interior particle velocities
c             from  2:nxt;2:nyt;2:nzt
c  -------------------------------------------------------------
c

         call dvel(nxt,nyt,nzt,dh,dt)
c
ccc   compute absorbing boundary conditions in velocities
c
         call abc(nxt,nyt,nzt,dt,dh)
c
ccc   compute velocities of 2nd order accuracy at boundary
c
         call bnd2d(nxt,nyt,nzt,dt,dh)
c      
ccc   exponential damping along boundaries
c
c         call dahs3d(nxt,nyt,nzt,ndamp)
c
c  -------------------------------------------------------------
ccc   4th-order differencing of pressure
c  -------------------------------------------------------------
c
         call dstres(nxt,nyt,nzt,dh,dt)
c
cc   compute stress tensor of 2nd order accuracy at boundary
c
         call strbnd(nxt,nyt,nzt,dh,dt)
c
c   ======================================================================
c CRACK BOUNDARY CONDITIONS  (applied right after stress update)
c   ======================================================================
c
c   The boundary conditions are applied now in ABSOLUTE STRESS, not in
c     relative stress, iotherwise we ghet all sort of problems
c
 
        do i=1,nxt
           do k=1,nzt
c
              tabsx = xy(i,nysc,k)+striniy(i,k)
              tabsz = yz(i,nysc,k)+strinix(i,k)
c              phi = atan2(tabsy,tabsx)
              if(tabsz.gt.peak_xz(i,k)) then
                 broken(i,k)=1
                 incrack(i,k)=1
                 ruptime(i,k)=time
c                 if(idebug.eq.1)then
c                    print *,i,k, ruptime(i,k)
c                 endif
             endif
           enddo
        enddo
c
c  Apply Boundary Conditions
c
        j=nysc
        do i=1,nxt
           do k=1,nzt
c
c   Do boundary conditions only for points that are considered broken,
c   i.e. those that are slipping
c
              if (incrack(i,k).eq.1) then
                 if (broken(i,k).eq.1) then
                    if(w1(i,j+1,k).lt.0.) then
                       broken(i,k)=0.
                       u1(i,j,k)=0.
                       w1(i,j,k)=0.
                    else
c
c
                       sliprate = w1(i,nysc+1,k)
                       gliss=slip(i,k)
c                       gliss = u2(i,nysc+1,k)
c
c   Apply friction law to absolute slip, absolute sliprate and absolute
c       traction
c
c   Classical slip weakening law
c

                       if(dmax.eq.0.)then
                          dd=0.
                       else
                          if(gliss.le.dmax)then
                             dd=1-gliss/dmax
                          else
                             dd=0.
                          endif
                       endif
                       friction=peak_xz(i,k)*dd
c -strinix(i,j)
                       yz(i,j,k) = friction-strinix(i,k) 
                     endif
                 endif
              endif
              slip(i,k)=slip(i,k)+w1(i,nysc+1,k)*dt
           enddo
        enddo
c
ccc   write out seismograms
c
         call wrss3d(ntst,ntiskp,nxt,nyt,nzt,dt,dh,nholes,nbhx,nbhy,
     +        nbhz,nbgx,nedx,nskpx,nbgy,nedy,nskpy,ifre,ndamp,
     $        nxsc,nysc,nzsc,slip,ioutput)
c
 10   continue
      if(ioutput.eq.1) then
        open(88,file="ruptime.res")
        open(98,file="slip.res")
        open(99,file="stressdrop.res")
        do j=1,nzt
          do i=1,nxt
             write(88,*)ruptime(i,j)
             write(98,*)2.*slip(i,j)
             write(99,*)yz(i,nysc,j)
          enddo
c         write(88,'()')
        enddo
     
        close(88)
        close(98)
        close(99)
      endif
      amom=0
      do j = 1,nzt
        do i = 1,nxt
           amom=amom+2*slip(i,j)
        enddo
      enddo
c      print *,'in indyna  amom =',amom


c
c -------------------------------------------------------------
c   END main computational loop
c -------------------------------------------------------------
ccc   close files
c
      call open3d(0,ioutput,iproc)
c
      end
c
c   ====================================================================
c              READ MODEL 3d
c   ====================================================================
c
      subroutine rdmd3d(vpe,vse,den,nxt,nyt,nzt,ifre)

ccc   routine to read the model 
c     nxt     # of node points in x-direction (integer)(sent)
c     nyt     # of node points in y-direction (integer)(sent)
c     nzt     # of node points in z-direction (integer)(sent)
c     vpe(1)  lowest p-velocity               (real)   (returned)
c     vpe(2)  highest p-velocity              (real)   (returned)
c     vse(1)  lowest s-velocity               (real)   (returned)
c     vse(2)  highest s-velocity              (real)   (returned)
c     den(1)  lowest density                  (real)   (returned)
c     den(2)  highest density                 (real)   (returned)
c     ifre    flag; =0 flat fs, <>0  topo     (integer)(returned)

      include 'parstat'
      dimension vpe(*),vse(*),den(*)
      vpe(2)=0.
      vpe(1)=1.0E10 
      vse(2)=0.
      vse(1)=1.0E10 
      den(2)=0.
      den(1)=1.0E10

c
c LAYERED MEDIUM 
c

      open(10,file='model.dat')
c      write(21,*)' A ouvert model.dat en lecture'  
      do  k=nzt,1,-1
         read(10,*) vp,vs,dd
         if (vp.gt.vpe(2)) vpe(2)=vp
         if (vp.lt.vpe(1)) vpe(1)=vp
         if (vs.gt.vse(2)) vse(2)=vs
         if (vs.lt.vse(1)) vse(1)=vs
         if (dd.gt.den(2)) den(2)=dd
         if (dd.lt.den(1)) den(1)=dd
c
c distribue les valeurs sur chaque couche du modèle
c
         do j=1,nyt
            do i=1,nxt      
               lam1(i,j,k)=dd*(vp**2-2.*vs**2)
               mu1(i,j,k)=dd*vs**2
               d1(i,j,k)=dd
            enddo
         enddo
      enddo
      close(10)
      
      if (ifre.ne.0) then
         
         do  181 j=1,nyt
            do 181 i=1,nxt
               do 180 k=nzt,1,-1
                  if (lam1(i,j,k).gt.0.) then
                     jj1(i,j)=k
                     goto 181
                  endif
 180           continue
 181        continue
         endif
         
         goto 20
   15 stop 'end of model file' 
   20 continue
      return  
      end 
c
c   ====================================================================
c              READ DATA 3d
c   ====================================================================
c
      subroutine rddt3d(nt,nxt,nyt,nzt,dh,dt,domp,str,dip,rake,slip,
     +     isrc,ainc,az,
     +     nbgx,nedx,nskpx,nbgy,nedy,nskpy,
     +     ntiskp,ndamp,nsv,ifre,nxsc,nysc,nzsc,nholes,nbhx,nbhy,
     +     nbhz,idebug,ibc,dmax,vmax )
                                                  
ccc   routine to read model and source parameters
                                                                     
c     nt     # of time steps                       (integer)(returned)
c     dh     spatial discretization step           (real)   (returned)
c     dt     time discretization step              (real)   (returned)
c     domp   dominant period of source             (real)   (returned)
c     ndamp  # of points to be damped              (integer)(returned)  
c     ntiskp loop increment for time samples (t)   (integer)(returned)
c     nskpx  loop increment for time samples (x)   (integer)(returned)
c     nskpy  loop increment for time samples (y)   (integer)(returned)
c     nbgx   first line in x dir with receiver     (integer)(returned)
c     nedx   last  line in x dir with receiver     (integer)(returned)
c     nbgy   first line in y dir with receiver     (integer)(returned)
c     nedy   last  line in y dir with receiver     (integer)(returned)
c     nxsc   x node point for source               (integer)(returned)
c     nysc   y node point for source               (integer)(returned)
c     nzsc   z node point for source               (integer)(returned)
c     nsv    save model every nsv timestep         (integer)(returned)
c     ifre   flag; =0 flat fs, <>0 topo            (integer)(returned)
c     nholes # of boreholes to write seismograms   (integer)(returned)
c     nbhx   array containing x-nodes of boreholes (integer)(returned)
c     nbhy   array containing y-nodes of boreholes (integer)(returned)
                                                                    
      include 'parstat'
      dimension nbhx(1),nbhy(1),nbhz(1)
      character*80 aa
c      print *,'read from unit 22'
      read(22,*) nxt,nyt,nzt                 
      read(22,*) nt
c      print *,'nt =',nt
      read(22,*) dh
      read(22,*) dt
      read(22,*) nxsc
      read(22,*) nysc
      read(22,*) nzsc
      read(22,*) nbgx
      read(22,*) nedx
      read(22,*) nskpx
      read(22,*) nbgy
      read(22,*) nedy
      read(22,*) nskpy
      read(22,*) ntiskp
      read(22,*) ndamp
      read(22,*) ifre
c
c  read debug 1 yes 
c
      read(22,*) idebug
c
c  read type of boundary condition
c
      read(22,*) ibc
c 
c  read friction law parameters dmax and vmax
c
      read(22,*) dmax
      read(22,*) vmax
c
      close(22)

      str=str*3.1415927/180.
      dip=dip*3.1415927/180.
      rake=rake*3.1415927/180.

      return                                                       
      end                                                         
c
c======================================================================
c             OPEN FILES
c======================================================================
c 
      subroutine open3d(io,ioutput,iproc)
      character*25 filename

c     routine to open and close files

c     io   unit # (integer) (sent)
!      print *,'in open3d ioutput = ',ioutput
      if (io.eq.1) then

         open(21,file='/tmp/verif.res')
         open(22,file='input.dat')
         if(ioutput.eq.1) then
!              print *, 'in open3d with ioutput = ',ioutput
              open(23,file='result/ssx3d.res')
              open(35,file='result/disp.res')
              open(36,file='result/sliprate.res')
         else
!              print *, 'in open3d with ioutput = ',ioutput
              ip=mod(iproc,10)
              write(filename,1000)ip
1000          format('/tmp/sliprateruiz3',i1,'.res')
c              write(*,*)'open file iproc = ',ip, filename
              open(36,file=filename)
         endif
c
c section temps
c
         
      else
         
         close(21)
         close(22)
         close(23)
         close(26)
         close(27)

         close(35)
         close(36)
         close(37)
         close(38)
         close(39)
         close(40)
         
      endif
      
      return
      
      end
c
c   ====================================================================
c              WRITE MODEL 3d
c   ====================================================================
c
 
      subroutine wrpm3d(nxt,nyt,nzt,nt,dh,dt,nsou,               
     +     domp,ndamp,ntiskp,nbgx,nedx,nskpx,nbgy,nedy,
     +     nskpy,vpe,vse,den,str,dip,rake,slip,
     +     ifre,axx,ayy,azz,axy,axz,ayz,nxsc,nysc,nzsc,
     +     peak,ruptur,rupmax, strout, strasp,itype,
     +     dmax,vmax, xmin,xmax,ymin,ymax )
                                                                      
ccc   routine to write model and source parameters                   

c     nxt    # of node points in x-direction       (integer)(sent)
c     nyt    # of node points in y-direction       (integer)(sent)
c     nzt    # of node points in z-direction       (integer)(sent)
c     nt     # of time steps                       (integer)(sent)    
c     dh     spatial discretization step           (real)   (sent)    
c     dt     time discretization step              (real)   (sent)    
c     domp   dominant period of source             (real)   (sent)    
c     ndamp  # of points to be damped              (integer)(sent)      
c     ntiskp loop increment for time samples (t)   (integer)(sent)    
c     nskpx  loop increment for time samples (x)   (integer)(sent)    
c     nskpy  loop increment for time samples (y)   (integer)(sent)    
c     nbgx   first line in x dir with receiver     (integer)(sent)    
c     nedx   last  line in x dir with receiver     (integer)(sent)    
c     nbgy   first line in y dir with receiver     (integer)(sent)    
c     nedy   last  line in y dir with receiver     (integer)(sent)    
c     nxsc   x node point for source               (integer)(sent)    
c     nysc   y node point for source               (integer)(sent)    
c     nzsc   z node point for source               (integer)(sent)    
c     nsv    save model every nsv timestep         (integer)(sent)    
c     ifre   flag; =0 flat fs, <>0 topo            (integer)(sent)    
c     nsou   number of points in source            (integer)(sent)    
c     str   strike for rupture mechanism           (real)   (sent)
c     dip   dip    for rupture mechanism           (real)   (sent)
c     rake  rake   for rupture mechanism           (real)   (sent)
c     slip  slip                                   (real)   (sent)
c     axx   xx stress for moment tensor            (real)   (sent)
c     ayy   yy stress for moment tensor            (real)   (sent)
c     azz   zz stress for moment tensor            (real)   (sent)
c     axy   xy stress for moment tensor            (real)   (sent)
c     ayz   xz stress for moment tensor            (real)   (sent)
c     axz   xz stress for moment tensor            (real)   (sent)
c     vpe(1)  lowest p-velocity                    (real)   (returned)
c     vpe(2)  highest p-velocity                   (real)   (returned)
c     vse(1)  lowest s-velocity                    (real)   (returned)
c     vse(2)  highest s-velocity                   (real)   (returned)
c     den(1)  lowest density                       (real)   (returned)
c     den(2)  highest density                      (real)   (returned)
                                                                    
      include 'parstat'
                                                             
      dimension vpe(*),vse(*),den(*)
                                                                 
      write(21,70) vpe(2)*dt/dh
      write(21,55) nxt                                               
      write(21,56) nyt                                              
      write(21,57) nzt                                             
      write(21,58) nt                                             
      write(21,59) dh                                            
      write(21,60) dt                                           
      write(21,61) domp                                        
      write(21,62) ndamp                                      
      write(21,84) vpe(2)                                              
      write(21,85) vpe(1)                                             
      write(21,86) vse(2)                                            
      write(21,87) vse(1)                                           
      write(21,88) den(2)                                            
      write(21,89) den(1)                                           
      write(21,90) dip*180./3.1415927
      write(21,91) str*180./3.1415927
      write(21,92) rake*180./3.1415927
      write(21,93) slip
      write(21,94) axx
      write(21,95) ayy
      write(21,96) azz
      write(21,97) axy
      write(21,98) ayz
      write(21,99) axz
      write(21,100) nxsc
      write(21,101) nysc
      write(21,102) nzsc
      if (ifre.eq.0) then
      write(21,48)                                                     
      else
      write(21,49)                                               
      endif
                                                                  
      write(21,28) ntiskp
                                                               
      write(21,72) nsou                                       
                                                             
  106 format(50i6)  
   28 format('skip of seismograms in time (loop counter)',i5)          
   35 format(4i8)                                                     
   36 format(10i4)                                                   
   48 format('using flat free surface boundary condition')
   49 format('using irregular free surface boundary condition')
   51 format('receiver lines in x dir')                              
   52 format('receiver lines in x dir')                             
   53 format('receiver lines in x dir')                            
   55 format('# of x nodes', i5)                                       
   56 format('# of y nodes',i5)                                        
   57 format('# of z nodes',i5)                                       
   58 format('# of time steps',i5)                                   
   59 format('discretization in space',f10.4)                       
   60 format('discretization in time' ,f10.4)                       
   61 format('dominant source period',f8.3)                    
   62 format('exp damping # of points',i5)                    
   65 format('source direction',9i5)                            
   70 format('stability criteria .5 > vp(max)*dt/dx = ',e15.7)    
   72 format('# of samples in source',i10)                            
   84 format('highest p-velocity encountered',f10.2)              
   85 format('lowest  p-velocity encountered',f10.2)              
   86 format('highest s-velocity encountered',f10.2)            
   87 format('lowest  s-velocity encountered',f10.2)                    
   88 format('highest density    encountered',f10.2)                       
   89 format('lowest  density    encountered',f10.2)                     
   90 format('fault dip',f10.2)                     
   91 format('fault strike',f10.2)                     
   92 format('fault rake',f10.2)                     
   93 format('fault slip',f10.2)                     
   94 format('moment tensor element axx',f10.2)                     
   95 format('moment tensor element ayy',f10.2)                     
   96 format('moment tensor element azz',f10.2)                     
   97 format('moment tensor element axy',f10.2)                     
   98 format('moment tensor element ayz',f10.2)                     
   99 format('moment tensor element axz',f10.2)                     
  100 format('source node position (x)',i5)                     
  101 format('source node position (y)',i5)                     
  102 format('source node position (z)',i5)                     
  103 format('plane wave incidence angle',f10.2)                     
  104 format('plane wave azimuth',f10.2)                     
c
c     print them out
c
      write(21,*)
      write(21,*) '   SOURCE DATA '
      write(21,*) ' '
      write(21,*)' center of fault: nxrsc,nysc,nzsc = ', nxsc,nysc,nzsc
      write(21,*)' grid size :  nxt,nyt,nzt = ', nxt,nyt,nzt
      write(21,*)' slip weaken. zone dmax   =', dmax
      write(21,*)' velocity weaken. zone vmax,   =', vmax
c      call flush(21)
                                                                     
      return                                                        
      end

c
c   ====================================================================
c              WRITE SEISMOGRAMS
c   ====================================================================
c
 
      subroutine wrss3d(i,ntiskp,nxt,nyt,nzt,dt,dx,nholes,nbhx,nbhy,
     +     nbhz,nbgx,nedx,nskpx,nbgy,nedy,nskpy,ifre,ndamp,
     $     nxsc,nysc,nzsc,slip,ioutput)
                                                                     
c     routine to save seismograms
                                                                    
c     i      current timestep                          (int)(sent)
c     ntiskp Step length in time for seismograms       (int)(sent)     
c     nxt    nodal points in x direction               (int)(sent)    
c     nyt    nodal points in y direction               (int)(sent)   
c     nzt    nodal points in z direction               (int)(sent)  
c     nskpx  loop increment in x dir for receivers     (int)(sent) 
c     nskpy  loop increment in y dir for receivers     (int)(sent)
c     nbgx   first line in x dir to contain receivers  (int)(sent)    
c     nedx   last  line in x dir to contain receivers  (int)(sent)   
c     nbgy   first line in y dir to contain receivers  (int)(sent)  
c     nedy   last  line in y dir to contain receivers  (int)(sent) 
c     dx     discretization in space                  (real)(sent)    
c     dt     discretization in time                   (real)(sent)   
c     ifre   flag; =0 flat fs, <>0 topo            (integer)(sent)    
c     ndamp  # of points to be damped              (integer)(sent)
                                                                    
      include 'parstat'
      real slip(nx,ny)
      dimension nbhx(1),nbhy(1),nbhz(1)
      nbgz=1
      nedz=nzt
      nskpz=1
      itt=ntiskp*(i/ntiskp)                                        
      if(itt.ne.i) return          
c     
      nbgz=1
      nedz=nzt
      nskpz=1

c
c write cross sections
c    
      do k=nbgz,nedz,nskpz
         write(36,'(1g12.4)') (2.*w1(ii,nysc+1,k),
     +        ii=nbgx,nedx,nskpx)
      enddo
c    
c      print *,'i, ioutput =',i, ioutput
c
      if(ioutput.eq.1) then
         do k=nbgz,nedz,nskpz
            write(23,'(1g12.4)') (yz(ii,nysc,k),
     +        ii=nbgx,nedx,nskpx)
         enddo
         do k=nbgz,nedz,nskpz
            write(35,'(1g12.4)') (slip(ii,k),
     +        ii=nbgx,nedx,nskpx)
         enddo
      endif
c    
         
      return                                                   
      end                                                      














