c-------------------------------------------------------------------------
c
c    PROGRAM kine_na_direct
c
c-------------------------------------------------------------------------


c                                    initialize receiver
c                                    function forward modelling
      include 'kine_subs/kine_param.inc'
      real*4            misfitval
      real*4            model(100)
      
      real*8            rmodel(100)

      real*4            observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave),
     &                  weight(maxwave,3),
     &                  constant_a(maxwave),
     &                  constant_c(maxwave),
     &                  time_shift(maxwave),
     &                  time_begin(maxwave),
     &                  time_end(maxwave)
      real*4            slip0(1600)
      real*4            estk,edip
      
      integer            ndd,nstk,ndip,nell,islip0,ishape
      integer            ndata, nwave, nsource
      
      logical            verbose,debug,timing,summary
      logical            lroot
      
      character*40      chars,kname,fname(maxwave),out_name(maxwave)
c     Info and Logical unit common 
c     blocks used by NA routines
c     
      character   headfile*50,sourcefile*10,statfile*10,fileout*10,
     &     histfile*50, filename*50, axires*50
      
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &     lu_nad,verbose,debug,timing,summary
      common /NAMPI/iproc,nproc,lroot
      common /rfi_com/observed_data, predicted_data, weight,
     &     constant_a, constant_c, time_shift,
     &     time_begin, time_end, ndata,
     &     nwave, fs, lu_mod, fname, nd
      common /ellipse/ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape
      
c     
c     Set up logical units 
c     LUs for standard input and 
c     output
      lroot=.TRUE.
      lu_rfi = 15
      lu_out = 6
      lu_od = 28
c     
c     
c     Set up logical units for kine files
c     LU for input of velocity 
c     model
      lu_vel = 11
c                                    LU for output of model 
c                                    parameters
!
!   kine_files/kine.in   contiene todos los datos para la inversion.
!   na.in  contiene datos de control de inversion con NA
!
      open(lu_rfi,file='kine_files/kine.in',status='old')
      
      read(lu_rfi,*)
      read(lu_rfi,*)
      read(lu_rfi,*) 
      
      write(lu_out,*)
      write(lu_out,*)' Kine_direct output'  
      write(lu_out,*)
      
      read(lu_rfi,'(a)') chars
      lw=lofw(chars)
      write(kname,'("kine_files/",a,a1)') chars(1:lw),char(0)
      write(lu_out,*)
!      write(lu_out,*) '* Now open ... ',kname(1:lw+19)
      open(lu_vel,file=kname,status='old')
      read(lu_vel,*) nd

      read(lu_rfi,'(a)') chars
      lw=lofw(chars)
!      write(kname,'("kine_files/",a,a1)') chars(1:lw),char(0)
      write(lu_out,*)
      write(lu_out,*) '* Now open ... ',kname(1:lw+19)
      open(lu_mod, file=kname, status='unknown')
      
c     INFO SUR LES ELLIPSES
      
      read(lu_rfi,*) ndd
      read(lu_rfi,*) estk,edip,nstk,ndip
      read(lu_rfi,*) nell
      read(lu_rfi,*) islip0
      read(lu_rfi,*) ishape
      write(*,*) "ndd,estk,edip,nstk,ndip,nell,islip0,ishape,nd"
      write(*,*) ndd,estk,edip,nstk,ndip,nell,islip0,ishape,nd
      do i=1,ndd
            slip0(i)=0
      enddo
c     
      close(lu_vel)

      nsource=nd
c     
c     Read in waveforms for nwave stations
      read(lu_rfi,*) nwave
      
      
      do iw=1,3
         read(lu_rfi,'(a)') chars
         lw=lofw(chars)
         write(fname(iw),'("kine_files/",a,a1)') 
     &        chars(1:lw),char(0)
         write(lu_out,*) '* Read in data from ',
     &        fname(iw)(1:lw+20)
         
         read(lu_rfi,*) time_begin(iw), time_end(iw), ndata
         open(lu_od,file=fname(iw),status='old')
         do ii=1,ndata*nwave
            read(lu_od,*) observed_data(ii,iw)
         enddo
         close(lu_od)
         
      end do
      
      do iww=1,nwave
         read(lu_rfi,*)(weight(iww,iw),iw=1,3)
!         write(*,*)(weight(iww,iw),iw=1,3)
      enddo
      close(lu_rfi)
c     
c   read the Gren function data base axi.res  
c      
      call readaxires
      
      
c     
c    Read best model
c     output
      open(lu_vel,file='kine_files/model_best.dat',status='old')
c
      read(lu_vel,*)
      read(lu_vel,*)
      read(lu_vel,*)
c      nsource = nd
c      print *,'il y a nd = ',nd,' parameters'
      do i=1,nd
        read(lu_vel,*) indice, model(i)
!        print *,'i, model =',i,model(i)
      end do
       
      do j=1,nd
         rmodel(j) = dble(model(j))
      end do
c
c     Perform forward modelling 
c     
!      print *,'call forward_model nd,ndata,nwave = ',nd,ndata,nwave 
      call forward_modelling_new(rmodel,nd,ndata,nwave, 
     &     predicted_data)
c     
c    Calculate misfit function
c     
!      print *,'call misfit ndata,nwave = ',ndata,nwave 
      call calcmisfit(
     &     predicted_data, observed_data, ndata,
     &     nwave, misfitval,weight)
      
      print *,'misfit = ',misfitval 
c      print *,'call plot_model '
c
c          PLOT_MODEL
c
      call plot_model(rmodel,nd,ndata,nwave) 
c      print *,'nd = ',nd
c     
c     write out seismograms for best model
c      
      open(lu_rfi,file='kine_files/kine.in',status='old')
      do iw=1,11
          read(lu_rfi,*)
      enddo
       
      do iw=1,3
         read(lu_rfi,'(a)') chars
         lw=lofw(chars)
         write(out_name(iw),'("kine_files/best_disp_",a,a1)') 
     &               chars(lw:lw),char(0)
        write(*,*) 'file best seismograms ', out_name(iw)(1:lw+20)
        open(14,file=out_name(iw)(1:lw+20),status='unknown')
        do ii=1,ndata*nwave
            write(14,*) predicted_data(ii,iw)
        end do
        read(lu_rfi,*)
        close(14)
      end do

      stop
      end

c------------------------------------------------------------------------
c
c	Subroutine forward_modelling - calculates the array "predicted_
c				       data" from the model "rmodel"
c
c	Note_1: Calls are made to : theo
c
c       Note_2: rmodel(1:nlayer): thickness
c               rmodel(nlayer+1:2*nlayer): velocity at upper interface
c               rmodel(2*nlayer+1:3*nlayer): velocity at lower interface
c               rmodel(3*nlayer+1:4*nlayer): Vp/Vs ratio
c
c-------------------------------------------------------------------------

      subroutine forward_modelling_new(
     &              rmodel, nd,
     &              ndata, nwave,
     &              predicted_data)

      include 'kine_param.inc'

      character*40    chars, fname(maxwave)
      character       sourcefile*20,statfile*20, axihist*50, disp_name*40
      real            ts(1600), phi_pi(1600), slip(1600), delta_pi(1600)
      real            rake(1600), ddip(1600), dstk(1600)
      integer            nd,tt
      real*8            rmodel(max_nd)
      real        nlayer
      real*4      predicted_data(maxdata,maxwave)

      integer     ndata,is(1600),
     &     nwave,nsource

      real        pi,cte1,cte2
      real        a,b,alpha,np,tp,x01,y01,xe,ye,xpos,ypos,xx,yy
      real        gauss,d,smax,vr,rak,delay
      dimension    iis(1600),iid(1600)
      real*4      slip0(1600)
      real*4      estk,edip

      integer       ndd,nstk,ndip,nell,islip0,ishape

      namelist /input/ nc,nfreq,tl,aw,nr,ns,xl,t0,fr1,fr2,
     &                    ikmax,uconv,sourcefile,statfile
  
      common /ellipse/ ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape
      
      pi=3.14159265359
c      print *,' forward_modellin, ndata, nwave = ',ndata,nwave

c       logical

      in1=17
      in2=18
 
      nlayer=nc
      nsource=nd
c                         
c       les fichiers source et station sont les memes
c      dans ce cas particulier
c      seuls le fichier axi.hist change
c       D'ABORD'
c      ici une procedure pour effacer et recreer le fichier axi.hist
c      lire le input list

      axihist='axi.hist'
c      open(in1,file=axihist,status='unknown')
      open(in1,file=axihist,status='unknown')
 
 

      do i=1,ndd
        read(in1,*) is(i), slip(i), phi_pi(i), delta_pi(i), rake(i),
     &        ddip(i), dstk(i), ts(i)
c     slip(i)=rmodel(i)
      enddo
      
c     do i=1,nd
c     slip(i)=rmodel(i)
c     enddo
      
      rewind(in1)
      
      
c********************************************************************
c     
c     PREMIERE ELLIPSE LIEE A L'HYPOCENTRE
c
     
      a=rmodel(1)*1000
      b=rmodel(2)*1000
      alpha=rmodel(3)*pi
      np=rmodel(4)
      tp=rmodel(5)*2*pi
      smax=rmodel(6)
      vr=rmodel(7)
      rak=0

c      if (b > a) then
c         cte1=a
c         cte2=b
c         b=cte1
c         a=cte2
c      end if

      
c******
c     centre de l'ellipse
c     
      
      x01=a*np*cos(tp)
      y01=b*np*sin(tp)
      xe=x01*cos(alpha)+y01*sin(alpha)+estk
      ye=-x01*sin(alpha)+y01*cos(alpha)+edip
      
      i=1
      do idip=1,ndip
         do istk=1,nstk
            iis(i)=istk
            iid(i)=idip
            i=i+1
         enddo
      enddo
      
      do i=1,ndd
         xpos=dstk(i)*iis(i)
         ypos=ddip(i)*iid(i)
c     ts(i)=((xpos-dstk(i)/2-xe)**2+
c     &                    (ypos-ddip(i)-ye)**2)**0.5/vr/1000
         ts(i)=((xpos-dstk(i)/2-estk)**2+
     &        (ypos-ddip(i)-edip)**2)**0.5/vr/1000
         xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
         yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
         d=xx**2/a**2+yy**2/b**2
         if (ishape.eq.0) then
            gauss=smax
         elseif (ishape.eq.1) then
             gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
c            gauss=smax*exp(log(0.01)*(xx**2/a**2+yy**2/b**2))
         elseif (ishape.eq.2) then
            if (d.ge.1.) then
            gauss=0
            else
            gauss=smax*sqrt(1-d)
            endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
         endif
         
         if(d.le.1) then
            slip(i)=gauss+slip0(i)
         else
            slip(i)=0.+slip0(i)
c            ts(i)=0.
         endif
      enddo
c********************************************************************
c     
c      AUTRES ELLIPSES
      
      if (nell.gt.1) then
         
         do ie=2,2
            a=rmodel(6)*1000
            b=rmodel(7)*1000
            alpha=rmodel(8)*pi
c            alpha=rmodel(10)
            xe=rmodel(9)*1000
            ye=rmodel(10)*1000
            smax=rmodel(11)
            vr=rmodel(5)

      
            do i=1,ndd
               xpos=dstk(i)*iis(i)
               ypos=ddip(i)*iid(i)
               xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
               yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
               d=xx**2/a**2+yy**2/b**2
               if (ishape.eq.0) then
                  gauss=smax
               elseif (ishape.eq.1) then
                  gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
               elseif (ishape.eq.2) then
                if (d.ge.1.) then
                gauss=0
                else
                gauss=smax*sqrt(1-d)
               endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
              endif
               
               if(d.le.1) then
                  ts(i)=((xpos-dstk(i)/2-estk)**2+
     &                 (ypos-ddip(i)-edip)**2)**0.5/vr/1000
                  slip(i)=gauss+slip(i)
               else
                  ts(i)=ts(i)
                  slip(i)=0.+slip(i)
               endif
            enddo
         enddo
c         print *,' a fini'
         
      endif


c      3 ELIPSE AUTRES ELLIPSES
      
      if (nell.gt.2) then
         
         do ie=3,nell
            a=rmodel(12)*1000
            b=rmodel(12)*1000
c            alpha=rmodel(17)
            alpha=0
            xe=rmodel(13)*1000
            ye=rmodel(14)*1000
            smax=rmodel(15)
            vr=rmodel(5)

      
            do i=1,ndd
!               print *,'i, ndd = ',i,ndd
               xpos=dstk(i)*iis(i)
               ypos=ddip(i)*iid(i)
               xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
               yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
               d=xx**2/a**2+yy**2/b**2
               if (ishape.eq.0) then
                  gauss=smax
               elseif (ishape.eq.1) then
                  gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
               elseif (ishape.eq.2) then
            if (d.ge.1.) then
            gauss=0
            else
            gauss=smax*sqrt(1-d)
            endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
               endif
               
               if(d.le.1) then
                  ts(i)=((xpos-dstk(i)/2-estk)**2+
     &                 (ypos-ddip(i)-edip)**2)**0.5/vr/1000
                  slip(i)=gauss+slip(i)
               else
                  ts(i)=ts(i)
                  slip(i)=0.+slip(i)
               endif
            enddo
         enddo
         print *,' a fini'
         
      endif



c     
c     ECRITURE AXI.HIST SLIP
c         
      do i=1,ndd
         write(in1,831) is(i), slip(i), phi_pi(i), delta_pi(i),rake(i),
     &        ddip(i), dstk(i), ts(i) 
      enddo
      close(in1)
      
 831  format(i10,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,e8.1,1x,
     &     e8.1,1x,f8.3)
      
      
 1001 format(8f9.3)
c       
c      print *,' avant sismo ndata, nwave =',ndata,nwave
      call seismograms(ndata,nwave,predicted_data,3)       
      
      return
      end
      
c------------------------------------------------------------------------
c
c	Subroutine plot-model - calculates the array "predicted_
c				       data" from the model "rmodel"
c
c	Note_1: Calls are made to : theo
c
c       Note_2: rmodel(1:nlayer): thickness
c               rmodel(nlayer+1:2*nlayer): velocity at upper interface
c               rmodel(2*nlayer+1:3*nlayer): velocity at lower interface
c               rmodel(3*nlayer+1:4*nlayer): Vp/Vs ratio
c
c-------------------------------------------------------------------------
 
      subroutine plot_model(
     &              rmodel, nd,
     &              ndata, nwave)

      include 'kine_param.inc'

      character*40    chars, fname(maxwave)
      character       sourcefile*20,statfile*20, axihist*50, disp_name*40
      real            ts(1600), phi_pi(1600), slip(1600), delta_pi(1600)
      real            rake(1600), ddip(1600), dstk(1600)
      integer            nd,tt
      real*8            rmodel(max_nd)
      real        nlayer

      integer     ndata,is(1600),
     &     nwave,nsource

      real        pi,cte1,cte2
      real        a,b,alpha,np,tp,x01,y01,xe,ye,xpos,ypos,xx,yy
      real        gauss,d,smax,vr,rak,delay
      dimension    iis(10000),iid(10000)
      real*4      slip0(1600)
      real*4      displ(10000),retard(10000)
      real*4      estk,edip

      integer       ndd,nstk,ndip,nell,islip0,ishape

      namelist /input/ nc,nfreq,tl,aw,nr,ns,xl,t0,fr1,fr2,
     &                    ikmax,uconv,sourcefile,statfile

      common /ellipse/ ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape
      
      pi=3.14159265359
c        print *,' forward_modellin_ell'

c       logical

      in1=17
      in2=18
 
      nlayer=nc
      nsource=nd
c                         
c      lire le input list

      axihist='axi.hist'
c      open(in1,file=axihist,status='unknown')
      open(in1,form='formatted',file=axihist,status='unknown')


      do i=1,ndd
         read(in1,*) is(i), slip(i), phi_pi(i), delta_pi(i), rake(i),
     &        ddip(i), dstk(i), ts(i)
c     slip(i)=rmodel(i)
      enddo
      
      close(in1)
      
      open(99,file='kine_files/slip_fin.dat')
      open(88,file='kine_files/tr_fin.dat')
      
c********************************************************************
c     
c     PREMIERE ELLIPSE LIEE A L'HYPOCENTRE
c
      
      a=rmodel(1)*1000
      b=rmodel(2)*1000
      alpha=rmodel(3)*pi
      np=rmodel(4)
      tp=rmodel(5)*2*pi
      smax=rmodel(6)
      vr=rmodel(7)
      rak=0

c      if (b > a) then
c         cte1=a
c         cte2=b
c         b=cte1
c         a=cte2
c      end if

      
c******
c     centre de l'ellipse
c     
      
      x01=a*np*cos(tp)
      y01=b*np*sin(tp)
      xe=x01*cos(alpha)+y01*sin(alpha)+estk
      ye=-x01*sin(alpha)+y01*cos(alpha)+edip

      i=1
      ndp=ndip*4
      nst=nstk*4
      ddp=ddip(1)/4.
      dst=dstk(1)/4.
      npatches=ndp*nst
      do idip=1,ndp
         do istk=1,nst
            iis(i)=istk
            iid(i)=idip
            i=i+1
         enddo
      enddo
      print *,'first ellipse '
      print *, 'position = ', xe, ye
      print *,'axes, angle = ',a, b, rmodel(3)
      print *,'smax, vr =',smax, vr
      
      do i=1,npatches
         xpos=dst*iis(i)
         ypos=ddp*iid(i)
c     ts(i)=((xpos-dstk(i)/2-xe)**2+
c     &                    (ypos-ddip(i)-ye)**2)**0.5/vr/1000
         rupture_time=((xpos-dst/2-estk)**2+
     &        (ypos-ddp-edip)**2)**0.5/vr/1000
         xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
         yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
         d=xx**2/a**2+yy**2/b**2
         if (ishape.eq.0) then
            gauss=smax
         elseif (ishape.eq.1) then
             gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
c            gauss=smax*exp(log(0.01)*(xx**2/a**2+yy**2/b**2))
         elseif (ishape.eq.2) then
            if (d.ge.1.) then
            gauss=0
            else
            gauss=smax*sqrt(1-d)
            endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
         endif
         if(d.le.1) then
            displ(i)=gauss
            retard(i)=rupture_time
         else
            displ(i)=0.
            retard(i)=0.
         endif
      enddo
c********************************************************************
c     
c      AUTRES ELLIPSES
      
      if (nell.gt.1) then
         
         do ie=2,2
            a=rmodel(6)*1000
            b=rmodel(7)*1000
            alpha=rmodel(8)*pi
c            alpha=rmodel(10)
c            np=rmodel(4)
c            tp=rmodel(5)*2*pi
            xe=rmodel(9)*1000
            ye=rmodel(10)*1000
            smax=rmodel(11)
            vr=rmodel(5)
  
          
c            print *,'ellipse no ',ie
c            print *, 'position = ', xe, ye
c            print *,'axes, angle = ',a, b, alpha*180./pi
c            print *,'smax, vr =',smax, vr
      
            do i=1,npatches
               xpos=dst*iis(i)
               ypos=ddp*iid(i)
               xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
               yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
               d=xx**2/a**2+yy**2/b**2
               if (ishape.eq.0) then
                  gauss=smax
               elseif (ishape.eq.1) then
                  gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
               elseif (ishape.eq.2) then
                  if (d.ge.1.) then
                   gauss=0
                   else
                   gauss=smax*sqrt(1-d)
                  endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))

               endif
               if(d.le.1) then
                  retard(i)=((xpos-dst/2-estk)**2+
     &                 (ypos-ddp-edip)**2)**0.5/vr/1000
                  displ(i)=gauss+displ(i)
               else
                  retard(i)=retard(i)
                  displ(i)=displ(i)
               endif
              
            enddo
         enddo
         
      endif


c      3 ELIPSE AUTRES ELLIPSES
      
      if (nell.gt.2) then
         
         do ie=3,nell
            a=rmodel(12)*1000
            b=rmodel(12)*1000
c            alpha=rmodel(17)
            alpha=0
c            np=rmodel(4+(ie-1)*7)
c            tp=rmodel(5+(ie-1)+7)*2*pi
            xe=rmodel(13)*1000
            ye=rmodel(14)*1000
            smax=rmodel(15)
            vr=rmodel(5)


            do i=1,npatches
               xpos=dst*iis(i)
               ypos=ddp*iid(i)
               xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
               yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
               d=xx**2/a**2+yy**2/b**2
               if (ishape.eq.0) then
                  gauss=smax
               elseif (ishape.eq.1) then
                  gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
               elseif (ishape.eq.2) then
                  if (d.ge.1.) then
                   gauss=0
                   else
                   gauss=smax*sqrt(1-d)
                  endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))

               endif
               if(d.le.1) then
                  retard(i)=((xpos-dst/2-estk)**2+
     &                 (ypos-ddp-edip)**2)**0.5/vr/1000
                  displ(i)=gauss+displ(i)
               else
                  retard(i)=retard(i)
                  displ(i)=displ(i)
               endif
              
            enddo
         enddo
         
      endif
c     
c     ECRITURE AXI.HIST SLIP
         
      do i=1,npatches
         write(99,*)displ(i)
         write(88,*)retard(i)
      enddo
      close(99)
      close(88)
      
      return
      end
     
