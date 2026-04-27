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

      dimension predicted_data1(maxdata)
      integer ndata, nstaa, npts
      real*4 num,den,dobs,dpred,misfitval
      real*4 minmisfit(maxdata)
      dimension ki(maxdata)



c     --------------------------------------------
c     L2 norm misfit between observed and computed
c     --------------------------------------------
      npts=(ndata/nstaa)/3
      misfitval0=0.0
      num=0.0
      den=0.0
      do i=0,2   ! Loop over components
         do j=0,nstaa-1   ! Loop over stations
            do k=1,npts   ! Loop over points
               dobs=observed_data(k+(npts*3*j)+(npts*i))
               dpred=predicted_data(k+(npts*3*j)+(npts*i))
               
               num=num+(dobs-dpred)**2
               den=den+(dobs)**2
            end do
         end do
      end do
      
      misfitval=num/den
      
      return
      end
