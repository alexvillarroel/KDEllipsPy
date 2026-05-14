c Modif S.Peyrat 2002-2003
c
c-------------------------------------------------------------------------
c
c      Subroutine display_model - writes out the model contained in 
c                           array rmodel and misfit value to  lu
c
c      Note: Calls no other routines
c
c------------------------------------------------------------------------
c
      subroutine display_model(
     &            lu, imod, rmodel, moddim, misfitval)
c
      include      'fd3d_param.inc'
c
c
      real*4            rmodel(maxmoddim)
c
      real*4            misfitval
c
c                                    Write out model rmodel
c
      write(lu,801) imod, misfitval
c
      do j=1,moddim
        write(lu,811) j, rmodel(j)
      end do
c
  801      format( '  model:',i5,',   misfit value:',e14.6 )
  811      format( 5x,i3,1e14.6 )
c
      return
      end
