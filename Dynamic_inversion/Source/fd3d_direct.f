c
c  FD3D  02/2008 pronto per l'inversione
c  Madariaga Mars 2008 (DNS SRC)
c
      include 'fd3d_subs/fd3d_param.inc'
c
      logical            verbose,debug,timing,summary
c
      character*40      chars,
     &                  kname,
     &                  fname
      real*4            misfitval
      real*4            model(100)

      real*8            rmodel(100)
c
      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &           lu_nad,verbose,debug,timing,summary,lu_misf,lu_modls

      common /pred_com/predicted_data(maxdata),iacc
      common /obs_com/observed_data(maxdata)

      common /fd3d_com/nt,nxt,nyt,nzt,dh,ndat,obs
      common /mis_com/nstaa
      common /taille/xll(2,maxmoddim),yll(2,maxmoddim)

      integer            ndata

      lu_fd3d  = 15
      lu_vel = 11
      lu_obs=13
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
      open(lu_vel,file='model_best.dat',status='old')
c
      read(lu_vel,*)nd
      print *,'il y a nd = ',nd,' parameters'
      do i=1,nd
        read(lu_vel,*) model(i)
        print *,'i, model =',i,model(i)
      end do
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
c  Read Axitra Header and GREN FUNCTIONS
c 
      call readaxires()          
cc
c     Read nt number of time steps,
c     Read nstaa
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
      
      do j=1,nd
           rmodel(j) = dble(model(j))
      enddo
      ioutput=1
      call forward_modelling(rmodel,nd,ndata,ioutput,iproc,amom)
!      print *,'volvio de forward istop = ',istop 
      call calcmisfit(misfitval,ndata,obs,iproc,nstaa)
      print *,'misfitval = ',misfitval
      stop
      end










