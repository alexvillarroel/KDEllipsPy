c ----------------------------------------------------------------------------
c      INITIAL PARAMETERS INVERSION
c      S.PEYRAT 2002-2008
C      R. MADARIAGA 2009-2010
c
c      Inversion stressini or peakini or Dc
c
c      Forward : DYNAMIC RUPTURE MODELING by FD
c      from K.B. Olsen, R. Madariaga, C. Marcinkovich, S. Peyrat 1996-2003
c 
c      Inversion : Neighbourhood algorithm from M. Sambridge
c
c ----------------------------------------------------------------------------
c                                               M. Sambridge, RSES, ANU.
c                                    Last updated Sept. 1999.
c
c ----------------------------------------------------------------------------
      Program fd3d_na

c     MEMORY AND ARRAY SIZES
c     
c     Call NA routine to do the work            
      call na
      
      stop
      end
c     
c-------------------------------------------------------------------------
c
c     user_init - performs all user specific initialization tasks. 
c     
c
c     Modif S.Peyrat April 2002 - June 2003
c
c-------------------------------------------------------------------------

      subroutine user_init(nd,ranges,scales,nd1,nd2,iproc)
      include 'fd3d_subs/fd3d_param.inc'
c
      real*4            ranges(2,*),scales(*)
c
      logical            verbose,debug,timing,summary
c
      character*40      chars,
     &                  kname,
     &                  fname
      character*25 filename
c     Info and Logical unit common blocks used by NA routines
c
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &           lu_nad,verbose,debug,timing,summary,lu_misf,lu_modls

      common /pred_com/predicted_data(maxdata),iacc
      common /obs_com/observed_data(maxdata)
      common /fd3d_com/nt,nxt,nyt,nzt,dh,ndata,obs

      common /mis_com/nstaa
      common /taille/xll(2,maxmoddim),yll(2,maxmoddim)

      integer            ndata

      lu_fd3d  = 15
      lu_out = 6
c     LU for input of  model
      lu_vel = 11
c   
c     Open fd3d.in files NEEDS COMPLETE CHANGE
c
      open(lu_fd3d,file='fd3d.in',status='old')
      read(lu_fd3d,*)
      read(lu_fd3d,*)
      read(lu_fd3d,*)
c
c Read Parameter file fd3d_param
c
      read(lu_fd3d,'(a)') chars
      lw=lofw(chars)
      write(kname,'(a,a1)') chars(1:lw),char(0)
      open(lu_vel,file=kname,status='old')
c
      read(lu_vel,*) nd,nd1,nd2,scales(1)
c
      do i=1,nd
        read(lu_vel,*) ranges(1,i), ranges(2,i), scales(i+1)
      end do
c
      close(lu_vel)
c
c read input.dat first time to get grid size etc
c
      lu_input=67
      open(lu_input,file='input.dat',status='old')
      read(lu_input,*) nxt,nyt,nzt
      read(lu_input,*) nt
      read(lu_input,*) dh
      read(lu_input,*) dt
      close(lu_input)
c
c  GENERATE GREN FUNCTIONS
c 
      call readaxires()
c
c     Read in name of observed file TO_DELETE
c
c
c     Read nt number of time steps, nx, ny size of model
c     Read weight (?), iacc (old), nstaa
c
      read(lu_fd3d,*) ntobs
      read(lu_fd3d,*) nstaa
      close(lu_fd3d)
c
c  READ DATA 
c
      ndata=ntobs*nstaa*3
      open(51,file='real_disp_x')
      open(52,file='real_disp_y')
      open(53,file='real_disp_z')
      obs=0.
      k=1
      do j=1,nstaa
         do ii=1,ntobs
               read(51,*)observed_data(k)
               obs=obs+observed_data(k)**2
               k=k+1
         enddo
         do ii=1,ntobs
               read(52,*)observed_data(k)
               obs=obs+observed_data(k)**2
               k=k+1
         enddo
         do ii=1,ntobs
               read(53,*)observed_data(k)
               obs=obs+observed_data(k)**2
               k=k+1
         enddo
      enddo 
      close(51)
      close(52)
      close(53)
      print *,'rms observed = ',obs
c
      close(lu_fd3d)
      nummod=0
c
833   format("momentos",i2,".dat")
832   format("momentos",i1,".dat")
      if(iproc.lt.10)then
         write(filename,832)iproc
      else
         write(filename,833)iproc
      endif
      open(822,file=filename)
      return
      end
c
c-------------------------------------------------------------------------
c
c      forward - performs forward modelling for user supplied problem.
c              and calculates the misfit measure
c              between observation and prediction.
c
c      Input: 
c            nd            :Number of dimensions in parameter space
c            model(nd)            :input velocity model
c
c      Output:
c            lppd            :negative log(ppd)
c
c
c-------------------------------------------------------------------------
c
c       Modif S.Peyrat April 2002 - June 2003
c
c-------------------------------------------------------------------------

      subroutine forward(nd,model,lppd,nd1,nd2,iproc)

c                                    initialize receiver
c                                    function forward modelling
      include 'fd3d_subs/fd3d_param.inc'

      real*4            lppd
      real*4            misfitval
      real*4            model(nd)
      real*8            rmodel(max_nd)
c
      logical            verbose,debug,timing,summary

      character*40      fname
c
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &          lu_nad,verbose,debug,timing,summary,lu_misf,lu_modls


      common /pred_com/predicted_data(maxdata),iacc
      common /obs_com/observed_data(maxdata)
      common /fd3d_com/nt,nxt,nyt,nzt,dh,ndata,obs
      common /mis_com/nstaa
      common /taille/xll(2,maxmoddim),yll(2,maxmoddim)

      integer  ndata
c
      do j=1,nd
           rmodel(j) = dble(model(j))
      enddo
c     Perform forward modelling on model rmodel.
c
      nummod=nummod+1
c
c   compute forward model write into local files
c
      ioutput=0
      call forward_modelling(rmodel,nd,ndata,ioutput,iproc,amom)
c
c      print *,'in forward iproc, momento=', iproc, amom
      call calcmisfit(misfitval,ndata,obs,iproc,nstaa)
      lppd = misfitval
      write(822,'(i4,13g12.4)')iproc,misfitval,amom,model
      print *,iproc,misfitval
      return
      end
c
c-------------------------------------------------------------------------
c
c      writemodels - user supplied routine to write out models produced
c                  by Neighbourhood algorithm in user's own format.
c
c      Input: 
c            nd              :number of dimensions in parameter space
c            ntot              :number of models generated by NA 
c            models(nd,ntot)     :models generated by NA
c            misfit               :array of model misfits (-lppd's)
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
c
c       Modif S.Peyrat April 2002 - June 2003
c
c-------------------------------------------------------------------------
c
      subroutine writemodels
     &             (nd, ntot, models, misfit, ns1, ns2, itmax,
     &              nh_max, nh, header,nd1,nd2)

      include 'fd3d_subs/fd3d_param.inc'

c     NA variables and arrays
      real*4            models(nd,*)
      real*4            misfit(ntot)
      real*4            mfitmin
      real*4            mfitminc
      real*4            mfitmean

      character*(*)   header

      logical            verbose,debug,timing,summary

      real*8            rmodel(max_nd)
      real*4 misfitval
      integer            ndata

      character*40      fname

c
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &                lu_nad,verbose,debug,timing,summary,lu_mod,lu_mis
      
      common /pred_com/predicted_data(maxdata),iacc
      common /obs_com/observed_data(maxdata)
      common /fd3d_com/nt,nxt,nyt,nzt,dh,ndata,obs
      common /mis_com/nstaa
      common /taille/xll(2,maxmoddim),yll(2,maxmoddim)

c     write out models   at each iteration
      mfitmin = misfit(1)
      ns = ns1
      np = 0
      mopt = 1
c     turn off writing to standard out by setting lu to zero
      lu_out2 = lu_out
      lu_out2 = 0
c
c     ==========================================================
c     loop over iterations
c     ==========================================================
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
         np = np + ns
         ns = ns2
c     
         call output_summary(
     &        lu_out2, lu_sum, it-1, models(1,mopt), nd,
     &        np, mfitmin, mfitmean, mfitminc, mopt)
      end do
c
c     Write out final model
c     
      call display_final(lu_sum, models(1,mopt), nd, mfitmin)
c     
      do j=1,nd
         rmodel(j) = dble(models(j,mopt))
      end do
      
c     repeat forward modelling 
c     for optimum model
      print *,'call forward modelling for optimum model'
      print *,(rmodel(j),j=1,nd)
      open(121,file='model_best.dat')
      write(121,*)nd
      do j=1,nd
           write(121,*)rmodel(j)
      enddo
      close(121)
      ioutput=1
      call forward_modelling(rmodel, nd,ndata,ioutput,iproc,amom)
      call calcmisfit(misfitval,ndata,obs,iproc,nstaa)
      print *,'misfit = ',misfitval
c     
c     calculate size of header   TO DELETE ?
c
      do i=1,nd
         k1 = 2*(i-1)*6 + 25
         k2 = k1 + 11
      end do
      nh = k2
      print *, 'dans writemodels nh = ',nh
      if(nh.gt.nh_max)then
         write(*,*)
         write(*,*)' Error - header array too small'
         write(*,*)
         write(*,*)'         current size = ',nh_max
         write(*,*)'        required size = ',nh
         write(*,*)
         stop
      end if

c     write info into character string
      
      write(header(1:24),fmt='(4i6)')ns1,ns2,itmax,nd
      print *,' header = ',header(1:24)
        
 801  format( 'iteration:',i5,',  misfit: min=',e14.6,
     &       ', mean=',e14.6,', minc=',e14.6 )

      return
      end
