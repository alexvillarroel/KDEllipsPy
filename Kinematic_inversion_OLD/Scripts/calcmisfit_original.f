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



c     Define variables:
      REAL*4 :: pred_data(maxdata,maxwave),
     & obs_data(maxdata,maxwave), misfitval, num, den
      INTEGER :: ndata, nwave



c     --------------------------------------------
c     L2 norm misfit between observed and computed
c     --------------------------------------------

      num=0.0
      den=0.0 
      do i=1,3  ! Loop for channels
         do j=1,nwave   ! Loop for stations
            do k=1,ndata   ! Loop over samples

               dobs=obs_data((j-1)*ndata+k,i)
               dpred=pred_data((j-1)*ndata+k,i)

               num=num+(dobs-dpred)**2
               den=den+(dobs)**2

            end do 
         end do
      end do

      ! Misfit equation:
      misfitval=(num/den)



      return
      end
