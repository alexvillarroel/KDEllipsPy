c Modif S.Peyrat 2002-2003
c
c-------------------------------------------------------------------------
c
c      Subroutine calcmisfit - calculates a misfit value between 
c                        observed data and predicted data
c
c------------------------------------------------------------------------
c
      subroutine calcmisfit(misfitval,ndata,obs,iproc,nstaa)
      include 'fd3d_param.inc'

      common /pred_com/predicted_data(maxdata),iacc
      common /obs_com/observed_data(maxdata)



c     Defining variables:
      REAL, ALLOCATABLE :: azt_data(:,:)
      REAL, DIMENSION(nstaa) :: azi, tP, tS
      INTEGER :: i, j, k, m, ndata, nstaa, npts, 
     & start_data, end_data, time_window
      REAL :: sampling, 
     & r_obs, r_pred, t_obs, t_pred, z_obs, z_pred
      REAL*4 :: misfitval, num, den
      
c     Parameters:
      time_window = 20   ! Time window [s]
      sampling = 4       ! Sampling rate [Hz]

c     Loading azimuths and arrival time data data:
      allocate(azt_data(nstaa,3))
      open(2,file='azi_times.txt')
      do m=1,nstaa
         read(2,*) azt_data(m,:)
      end do
      close(2)
      
c     Creating vectors of azimuth and arrival times of each station:
      azi = azt_data(:,1)   ! Azimuths [radians]
      tP = azt_data(:,2)   ! P-wave arrival times [s]
      tS = azt_data(:,3)   ! S-wave arrival times [s]


c     --------------------------------------------
c     L2 norm misfit between observed and computed
c     --------------------------------------------

      npts=(ndata/nstaa)/3
      num=0.0
      den=0.0
      do i=0,2   ! Loop over components
         do j=0,nstaa-1   ! Loop over stations

            ! Considering P-wave in radial and vertical:
            if ((i == 0 ) .or. (i == 2)) then
               start_data = nint(tP(j+1)) * sampling
               end_data = start_data + (time_window * sampling)
            ! Considering SH-wave in tangential:
            else if (i == 1) then
               start_data = nint(tS(j+1)) * sampling
               end_data = start_data + (time_window * sampling)
            end if

            do k=start_data,end_data   ! Loop over samples

               ! Calculating radial component:
               if (i == 0) then
                  r_obs=(observed_data(k+(npts*3*j)+(npts*i))*
     &             cos(azi(j+1)))+
     &             (observed_data(k+(npts*3*j)+(npts*(i+1)))*
     &             sin(azi(j+1)))
                  r_pred=(predicted_data(k+(npts*3*j)+(npts*i))*
     &             cos(azi(j+1)))+
     &             (predicted_data(k+(npts*3*j)+(npts*(i+1)))*
     &             sin(azi(j+1)))
                  
                  num=num+(r_obs-r_pred)**2
                  den=den+(r_obs)**2
                  
               ! Calculating transverse component:
               else if (i == 1) then
                  t_obs=(observed_data(k+(npts*3*j)+(npts*i))*
     &             cos(azi(j+1)))-
     &             (observed_data(k+(npts*3*j)+(npts*(i-1)))*
     &             sin(azi(j+1)))
                  t_pred=(predicted_data(k+(npts*3*j)+(npts*i))*
     &             cos(azi(j+1)))-
     &             (predicted_data(k+(npts*3*j)+(npts*(i-1)))*
     &             sin(azi(j+1)))

                  num=num+(t_obs-t_pred)**2
                  den=den+(t_obs)**2   

               ! Calculating vertical component:
               else if (i == 2) then
                  z_obs=observed_data(k+(npts*3*j)+(npts*i))
                  z_pred=predicted_data(k+(npts*3*j)+(npts*i))

                  num=num+(z_obs-z_pred)**2
                  den=den+(z_obs)**2
               end if  

            end do
         end do
      end do

      ! Misfit equation. If doesn't nucleate (misfit=1.0), 
      ! set misfit to infinite:
      misfitval=(num/den)
      if (misfitval == 1.0) then
         misfitval=misfitval/0.0
      end if



      return
      end
