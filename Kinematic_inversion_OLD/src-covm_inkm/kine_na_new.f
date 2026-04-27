c
c ----------------------------------------------------------------------------
c
c                              Call NA routine to do the work            

      Program kine_na

      call na
        
      stop
      end
c-------------------------------------------------------------------------
c
c      user_init - performs all user specific initialization tasks. 
c                In this case for kinematic source  inversion.
c
c      Input: - none
c
c      Output:
c            nd            :Number of dimensions in parameter space
c            ranges(2,nd)      :bounds on parameter space
c            scales(nd+1)      :scale factors in parameter space
c
c-------------------------------------------------------------------------

      subroutine user_init(nd,ranges,scales)
c                                    initialize receiver
c                                    function forward modelling
      include 'kine_subs/kine_param.inc'

      real*4              ranges(2,*),scales(*)
      real*4            observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave),
     &                  weight(maxwave,3)
      real*4            slip0(1600)
      real*4            estk,edip
      
      integer            ndd,nstk,ndip,nell,islip0,ishape
      integer            ndata, nwave, nsource
      
      logical            verbose,debug,timing,summary
      logical            lroot
      
      character*40      chars,kname,fname(maxwave)
c     Info and Logical unit common 
c     blocks used by NA routines
c     
      character   headfile*50,sourcefile*10,statfile*10,fileout*10,
     &     histfile*50, filename*50, axires*50
      
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &     lu_nad,verbose,debug,timing,summary
      common /NAMPI/iproc,nproc,lroot
      common /rfi_com/observed_data, predicted_data, weight
      common /kine_data/ndata,
     &     nwave, fs, lu_mod, fname
      common /ellipse/ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape
      
c     
c     Set up logical units 
c
      lroot=.TRUE. 
      lu_rfi = 65
      lu_out = 6
      lu_od = 28
      lu_params = 61
      lu_mod = 63
      
      open(66,file='kine_files/misfits.dat',status='unknown')
!      open(67,file='minmisfits.dat',status='unknown')
      open(68,file='kine_files/models.dat',status='unknown')
!
!   kine_files/kine.in   contiene todos los datos para la inversion.
!   na.in  contiene datos de control de inversion con NA
!
      open(lu_rfi,file='kine_files/kine.in',status='old')
      
      read(lu_rfi,*)
      read(lu_rfi,*)
      read(lu_rfi,*) 
      
      if(lroot)write(lu_out,*)
      if(lroot)write(lu_out,*)' User routines output'  
      if(lroot)write(lu_out,*)
      if(lroot)write(lu_out,*)' Opening kine files...'
      
      read(lu_rfi,'(a)') chars
      lw=lofw(chars)
      write(kname,'("kine_files/",a,a1)') chars(1:lw),char(0)
      if(lroot)write(lu_out,*)
      if(lroot)write(lu_out,*) '* Now open ... ',kname(1:lw+19)
      open(lu_params,file=kname,status='old')
!     
!     Read in parameter file and ranges of the parameters
!
      read(lu_params,*) nd,scales(1)
      do i=1,nd
        read(lu_params,*) ranges(1,i), ranges(2,i), scales(i+1)
      end do

!      call read_vmodelNA(lu_params, ranges, scales, nd) 
      close(lu_params)
      
      nsource=nd
!     
      read(lu_rfi,'(a)') chars
      lw=lofw(chars)
      write(kname,'("kine_files/",a,a1)') chars(1:lw),char(0)
      if(lroot)write(lu_out,*)
      if(lroot)write(lu_out,*) '* Now open ... ',kname(1:lw+19)
      open(lu_mod, file=kname, status='unknown')
      
c     INFO ABOUT THE SOURCE FILE TO BE USED WITH seismogram 
      
      read(lu_rfi,*) ndd
      read(lu_rfi,*) estk,edip,nstk,ndip
      read(lu_rfi,*) nell
      read(lu_rfi,*) islip0
      read(lu_rfi,*) ishape
      print *,'leyo nell, nd, islip0,ishape =',nell,nd, islip0,ishape
      write(*,*) "ndd,estk,edip,nstk,ndip,nell,islip0,ishape"
      write(*,*) ndd,estk,edip,nstk,ndip,nell,islip0,ishape
      if (islip0.eq.1) then
         open(69,file='slip0.dat',status='old')
         do i=1,ndd
            read(69,*)slip0(i)
         enddo
         close(69)
      else
         do i=1,ndd
            slip0(i)=0
         enddo
      endif
!
!   observed seismograms in a single line
!
      read(lu_rfi,*) nwave
      print *,'leyo nwave = ',nwave
      do iw=1,3
         read(lu_rfi,'(a)') chars
         lw=lofw(chars)
         write(fname(iw),'("kine_files/",a,a1)') 
     &        chars(1:lw),char(0)
         if(lroot)write(lu_out,*)
         if(lroot)write(lu_out,*) '* Read in data from ',
     &        fname(iw)(1:lw+20)
         
         read(lu_rfi,*) time_begin, time_end, ndata
         if(lroot)write(lu_out,*) '* Time window =',time_begin,
     &        ' - ',time_end
         if(lroot)write(lu_out,*)' '
         
         open(lu_od,file=fname(iw),status='old')
         do ii=1,ndata*nwave
            read(lu_od,*) observed_data(ii,iw)
         enddo
         close(lu_od)
         
      end do
      print *,'antes de weights nwave = ',nwave
      
      do iww=1,nwave
         read(lu_rfi,*)(weight(iww,iw),iw=1,3)
!         write(*,*)(weight(iww,iw),iw=1,3)
      enddo
      close(lu_rfi)
c     
c   read the Gren function data base axi.res  
c      
      call readaxires
      print *,'despues de readaixxres  nwave = ',nwave
     
      if(lroot)write(lu_out,*)' finished initialisation'
      
 100  format(1x,60("-")//)
 101  format(/1x,60("-")/)
      
      return
      end
c
c-------------------------------------------------------------------------
c
c      forward - performs forward modelling for user supplied problem.
c
c      Input: 
c            nd            :Number of dimensions in parameter space
c            model(nd)            :input velocity model
c
c      Output:
c            lppd            :negative log(ppd)
c
c-------------------------------------------------------------------------

      subroutine forward(nd,model,lppd)
      
      include 'kine_subs/kine_param.inc'
      
      real*4            lppd
      real*4            misfitval
      real*4            model(nd)
      
      real*8            rmodel(max_nd)
      
      real*4            observed_data(maxdata,maxwave),
     &     predicted_data(maxdata,maxwave),
     &     weight(maxwave,3)
      real*4            slip0(1600)
      real*4            estk,edip
      
      integer            ndd,nstk,ndip,nell,islip0,ishape
      integer            ndata, nd
      logical            verbose,debug,timing,summary
      character*40      fname(maxwave)
c     
c     Info and Logical unit common 
c     blocks used by NA routines
      
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &                lu_nad,verbose,debug,timing,summary
      
      common /rfi_com/observed_data, predicted_data, weight
      common /kine_data/ndata,
     &     nwave, fs, lu_mod, fname
      
      common /ellipse/ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape
c     
      do j=1,nd
         rmodel(j) = dble(model(j))
         write(68,*)model(j)
      end do
!
!     Perform forward modelling and compute misfit
!      
!     print *,'in forward nwave, ndata = ',nwave,ndata
      call forward_modelling(rmodel,nd,ndata,nwave, 
     &     predicted_data)
c     
      call calcmisfit(
     &     predicted_data, observed_data, ndata,
     &     nwave, misfitval,weight)
      
      lppd = misfitval
      write(66,*)misfitval
      write(*,*)'misfitval = ',misfitval
      
      call flush(66)
      return
      end
c     
c     *-------------------------------------------------------------------------
c     
c*      writemodels - user supplied routine to write out models produced
c*                  by Neighbourhood algorithm in user own format.
c
c*      Input: 
c            nd              :number of dimensions in parameter space
c            ntot              :number of models generated by NA 
c            models(nd,ntot)     :models generated by NA
c            misfit               :array of model misfits (-lppd)
c            ns1                :initial sample size used by NA
c            ns2                :normal sample size used by NA
c            itmax                :number of iterations
c            nh_max                :maximum length of nad file header 
c
c      Output:         
c            nh                :length of nad file header 
c            header(nh)          :character string containing nad file header
c
c-------------------------------------------------------------------------

      subroutine writemodels
     &     (nd, ntot, models, misfit, ns1, ns2, itmax,
     &     nh_max, nh, header)
      
      include 'kine_subs/kine_param.inc'
      
c     NA variables and arrays
      real*4            models(nd,*)
      real*4            misfit(ntot)
      real*4            mfitmin
      real*4            mfitminc
      real*4            mfitmean
      
      character*(*)   header
      
      logical            verbose,debug,timing,summary
      
      real*8            rmodel(max_nd)
      
      real*4            observed_data(maxdata,maxwave),
     &     predicted_data(maxdata,maxwave),
     &     weight(maxwave,3)
      
      real*4            slip0(1600)
      real*4            estk,edip
      
      integer            ndd,nstk,ndip,nell,islip0,ishape
      integer            ndata
      
      character*40      chars, fname(maxwave), out_name(maxwave),toto
      
c     
c     Info and Logical unit common 
c     block used by NA routines
      
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &     lu_nad,verbose,debug,timing,summary
      
      common /rfi_com/observed_data, predicted_data, weight
      common /kine_data/ndata,
     &     nwave, fs, lu_mod, fname
      
      common /ellipse/ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape
     c
 
      nsource= nd
      lu_rfi = 19
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     write out models for each iteration
      mfitmin = misfit(1)
      ns = ns1
      np = 0
      mopt = 1
      write(*,*) 'here we are in write kine_models'
      lu_out2 = lu_out
      lu_out2 = 0
      write(lu_mod,*)ns1,' Number of samples in starting pool'
      write(lu_mod,*)ns2,' Number of new samples per iteration'
      write(lu_mod,*)itmax,' Number of iterations'
c     
c   ===============================================
c     loop over iterations
c   ===============================================
c
      do it=1,itmax+1
         mfitminc = misfit(np+1)
         mfitmean = 0.0
c     find minimum and mean misfit
         do i=1,ns
            jj = np + i
            if(misfit(jj).lt.mfitmin)then
               mfitmin = misfit(jj)
               mopt = jj
            end if
            mfitminc = min(mfitminc,misfit(jj))
            mfitmean = mfitmean + misfit(jj)
         end do
         mfitmean = mfitmean/ns
         write(lu_mod,801) it-1, mfitmin, mfitmean, mfitminc
!
!  write residuals of each model in this iteration to kine_models
!
         do i=1,ns
            jj = np + i
            write(lu_mod,807) jj, misfit(jj)
  807       format( '  model:',i5,',   misfit value:',f10.5 )
         end do
         np = np + ns
         ns = ns2
c     
         
         call output_summary(
     &        lu_out2, lu_sum, it-1, models(1,mopt), nd,
     &        np, mfitmin, mfitmean, mfitminc, mopt, nsource)
      end do
c
c     Write out best model
c      
      call display_final(
     &     lu_out, models(1,mopt), nsource, nd)
      
      call display_final(
     &     lu_sum, models(1,mopt), nsource, nd)
      
      call display_final(
     &     lu_mod, models(1,mopt), nsource, nd)
      
      open(94,file="kine_files/model_best.dat")
      call display_final(
     &     94, models(1,mopt), nsource, nd)
      close(94)
      do j=1,nd
         rmodel(j) = dble(models(j,mopt))
      end do
c      
c     repeat forward modelling 
c     for optimum model
c
      call forward_modelling(
     &     rmodel, nd, ndata, nwave, predicted_data)
c
c  plot best model
c
      call plot_model(rmodel,nd,ndata,nwave) 
c
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
      close(lu_rfi)
      
 801  format( 'iteration:',i5,',  misfit: min=',f10.5,
     &     ', mean=',f10.5,', minc=',f10.5 )
c
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

      namelist /input/ nc,nfreq,tl,aw,nr,ns,xl,
     &                    ikmax,uconv,sourcefile,statfile

      common /ellipse/ ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape
      
      pi=3.14159265359

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
      ndp=80
      nst=80
      ddp=400.
      dst=400.
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
c            alpha=0.
            xe=rmodel(9)*1000
            ye=rmodel(10)*1000
            smax=rmodel(11)
            vr=rmodel(5)
            
            print *,'ellipse no ',ie
            print *, 'position = ', xe, ye
            print *,'axes, angle = ',a, b, alpha*180./pi
            print *,'smax, vr =',smax, vr
      
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
                  gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
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

c********************************************************************
c     
c      TERCERA ELLIPSES
      
      if (nell.gt.2) then
         
         do ie=3,nell
            a=rmodel(12)*1000
            b=rmodel(12)*1000
c            alpha=rmodel(17)
            alpha=0.
            xe=rmodel(13)*1000
            ye=rmodel(14)*1000
            smax=rmodel(15)
            vr=rmodel(5)
            
            print *,'ellipse no ',ie
            print *, 'position = ', xe, ye
            print *,'axes, angle = ',a, b, alpha*180./pi
            print *,'smax, vr =',smax, vr
      
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
                  gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
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
     
