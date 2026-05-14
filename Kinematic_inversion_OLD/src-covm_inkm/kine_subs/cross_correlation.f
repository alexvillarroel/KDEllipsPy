c-------------------------------------------------------------------------
c
c      Subroutine cross_correlation - calculates the Pearson's correlation  
c                        coefficient (rho) between obs_data and pred_data
c                        and then delays pred_data according to the 
c                        location of the best value of rho.
c
c-------------------------------------------------------------------------
c

      subroutine cross_correlation(sobs,spred,ndata,
     &            nwave,new_spred)

      include 'kine_param.inc'
   
c     Define variables:
      INTEGER, PARAMETER :: tau = 250
      REAL*4 :: spred(maxdata,maxwave), sobs(maxdata,maxwave)
      INTEGER :: ndata, nwave, delta, i, j, k, m
      INTEGER, DIMENSION(1) :: pos_max
      REAL, DIMENSION(ndata) :: dobs, dpred, new_dpred
      REAL, DIMENSION((tau*2)+1) :: rho
      REAL*4 :: new_spred(maxdata,maxwave)



c     Numerical calculations:
      new_spred = 0.0
      do i=1,3   ! Loop for channels
         do j=1,nwave   ! Loop for stations

            ! Reordering data. Select seismograms of each station:
            dobs = sobs((j-1)*ndata+1:(j-1)*ndata+ndata,i)
            dpred = spred((j-1)*ndata+1:(j-1)*ndata+ndata,i)

            ! Calculating correlation coefficient:
            call cc(dobs,dpred,tau,size(dobs),rho)

            ! Search for the maximum correlation coefficient and define 
            ! the delay around the center of the cross-correlation array 
            ! (zero delay):
            pos_max = maxloc(rho)
            delta = abs((size(rho)/2.0)+0.5-pos_max(1))

            ! Delaying modeled signal according to calculated delay (delta):
            new_dpred = 0.0
            if (pos_max(1) < (size(rho)/2.0)+0.5) then
               ! If is a negative delay:
               do m=1,size(dpred)-delta
                  new_dpred(m) = dpred(delta+m)
               end do
            else if (pos_max(1) > (size(rho)/2.0)+0.5) then
               ! If is a positive delay:
               do m=1,size(dpred)-delta
                  new_dpred(delta+m) = dpred(m)
               end do
            end if

            ! Reordering delayed predicted data to its original dimension:
            new_spred((j-1)*ndata+1:(j-1)*ndata+ndata,i) = new_dpred

         end do
      end do
   

      return
      end


c     ----------------------------------------------------------
c     ------------ CROSS-CORRELATION SUBROUTINE ----------------
c     ----------------------------------------------------------

      subroutine cc(obs_data,pred_data,lag,lth,rho)

      ! Define variables:
      INTEGER :: lag, lth                          ! Input
      REAL, DIMENSION(lth) :: obs_data, pred_data  ! Input
      REAL, DIMENSION((lag*2)+1) :: rho            ! Output
      REAL, ALLOCATABLE :: s1(:), s2(:)
      REAL, ALLOCATABLE :: num1(:), num2(:), den1(:), den2(:)
      REAL :: k1, k2, k3, c, mean_s1, mean_s2
      INTEGER :: i



      ! -----------------------------------------------
      ! Cross-correlation between observed and computed   
      ! -----------------------------------------------

      ! Starting array of zeros:
      rho = 0

      ! For a negative delay:
      do i=1,lag
         allocate(s1(size(obs_data)-i))
         allocate(s2(size(pred_data)-i))
         allocate(num1(size(obs_data)-i))
         allocate(num2(size(pred_data)-i))
         allocate(den1(size(obs_data)-i))
         allocate(den2(size(pred_data)-i))

         ! Applying delay:
         s1 = obs_data(1:size(obs_data)-i)
         s2 = pred_data(i+1:size(pred_data))

         ! Calculating Pearson's correlation coefficient:
         mean_s1 = sum(s1)/size(s1)
         mean_s2 = sum(s2)/size(s2)

         num1 = s1-mean_s1
         num2 = s2-mean_s2
         k1 = sum(num1*num2)

         den1 = (s1-mean_s1)**2
         den2 = (s2-mean_s2)**2
         k2 = sum(den1)
         k3 = sum(den2)

         ! For cases when the function is undefined:
         if (k2 == 0 .or. k3 == 0) then
            c = 0
         else
            c = k1/(sqrt(k2)*sqrt(k3))
         end if

         ! Creating the array with correlation coefficient values:
         rho(lag+1-i)=c

         deallocate(s1)
         deallocate(s2)
         deallocate(num1)
         deallocate(num2)
         deallocate(den1)
         deallocate(den2)
      end do


      ! For a delay equal to zero or positive:
      do i=0,lag
         allocate(s1(size(obs_data)-i))
         allocate(s2(size(pred_data)-i))
         allocate(num1(size(obs_data)-i))
         allocate(num2(size(pred_data)-i))
         allocate(den1(size(obs_data)-i))
         allocate(den2(size(pred_data)-i))

         ! Applying delay:
         s1 = obs_data(i+1:size(obs_data))
         s2 = pred_data(1:size(pred_data)-i)

         ! Calculating Pearson's correlation coefficient:
         mean_s1 = sum(s1)/size(s1)
         mean_s2 = sum(s2)/size(s2)

         num1 = s1-mean_s1
         num2 = s2-mean_s2
         k1 = sum(num1*num2)
 
         den1 = (s1-mean_s1)**2
         den2 = (s2-mean_s2)**2
         k2 = sum(den1)
         k3 = sum(den2)

         ! For cases when the function is undefined:
         if (k2 == 0 .or. k3 == 0) then
            c = 0
         else
            c = k1/(sqrt(k2)*sqrt(k3))
         end if

         ! Creating the array with correlation coefficient values:
         rho(lag+1+i)=c

         deallocate(s1)
         deallocate(s2)
         deallocate(num1)
         deallocate(num2)
         deallocate(den1)
         deallocate(den2)
      end do
      end
