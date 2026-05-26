c Modif S.Peyrat 2002-2003
c
c-------------------------------------------------------------------------
c
c      Subroutine display_final - writes out the final model contained in 
c                           array rmodel to the file assigned to 
c                                  logical unit lu
c
c      Note: Calls no other routines
c
c------------------------------------------------------------------------
c
      subroutine display_final(lu, rmodel, moddim, misfitval)
c
      include 'fd3d_param.inc'
c
c
!      common /fd3d_com/nt,nx,ny,weight,obs,
!     &                  ndata,
!     &                  lu_obs, lu_dc, fname,
!     &                  nummod

      integer            ndata

      real*4            rmodel(maxmoddim)
c
      real*4            misfitval
c
c                                    Write out final model
c

      write(lu,801)
      do j=1,moddim
          write(lu,811)j,rmodel(j)
      end do
c
  801      format(/'*** Final model ***' /
     &          3x,'#',4x,'a1',4x,'b1',4x,'a2',
     &          4x,'b2')

  811      format( i10,e14.6)
c
      return
      end
