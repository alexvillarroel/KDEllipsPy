c-------------------------------------------------------------------------
c
c      Subroutine calcmisfit - calculates a misfit value between 
c                        obs data and pred data
c                        plus model roughness
c
c      Note: Calls no other routines
c
c------------------------------------------------------------------------
c
      subroutine calcmisfit(
     &            pred_data, obs_data, ndata, 
     &            nwave, misfitval )

      include 'kine_param.inc'



c     Defining variables:
      REAL, ALLOCATABLE :: azt_data(:,:)
      REAL, DIMENSION(nwave) :: azi, tP, tS
      INTEGER :: i, j, k, m, ndata, start_data, end_data,
     & time_window
      REAL :: dt, delta, sampling,
     & r_obs, r_pred, t_obs, t_pred, z_obs, z_pred
      REAL*4 :: pred_data(maxdata,maxwave), 
     & obs_data(maxdata,maxwave), misfitval, num, den
      
c     Parameters:
      time_window = 20   ! Time window [s]
      sampling = 1.0/dt  ! Sampling rate [Hz] 

c     Loading azimuths and arrival time data data:
      allocate(azt_data(nwave,3))
      open(2,file='azi_times.txt')
      do m=1,nwave
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

      num=0.0
      den=0.0 
      do i=1,3    ! Loop for channels
        do j=1,nwave   ! Loop for stations

           ! Considering P-wave in radial and vertical:
           if ((i == 1 ) .or. (i == 3)) then
              start_data = nint(tP(j)) * sampling
              end_data = start_data + (time_window * sampling)
           ! Considering SH-wave in transverse:
           else if (i == 2) then
              start_data = nint(tS(j)) * sampling
              end_data = start_data + (time_window * sampling)
           end if

           do k=start_data,end_data   ! Loop over samples

              ! Calculating radial component:
              if (i == 1) then
                 r_obs=(obs_data((j-1)*ndata+k,i)*cos(azi(j)))+
     &            (obs_data((j-1)*ndata+k,i+1)*sin(azi(j)))
                 r_pred=(pred_data((j-1)*ndata+k,i)*cos(azi(j)))+
     &            (pred_data((j-1)*ndata+k,i+1)*sin(azi(j)))

                 num=num+(r_obs-r_pred)**2
                 den=den+(r_obs)**2

              ! Calculating transverse component:
              else if (i == 2) then
                 t_obs=(obs_data((j-1)*ndata+k,i)*cos(azi(j)))-
     &            (obs_data((j-1)*ndata+k,i-1)*sin(azi(j)))
                 t_pred=(pred_data((j-1)*ndata+k,i)*cos(azi(j)))-
     &            (pred_data((j-1)*ndata+k,i-1)*sin(azi(j)))

                 num=num+(t_obs-t_pred)**2
                 den=den+(t_obs)**2

              ! Calculating vertical component:
              else if (i == 3) then
                 z_obs=obs_data((j-1)*ndata+k,i)
                 z_pred=pred_data((j-1)*ndata+k,i)

                 num=num+(z_obs-z_pred)**2
                 den=den+(z_obs)**2
              end if

           end do 
        end do
      end do

      ! Misfit equation:
      misfitval=num/den

      return
      end
