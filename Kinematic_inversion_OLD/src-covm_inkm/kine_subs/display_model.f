c-------------------------------------------------------------------------
c
c      Subroutine display_model - writes out the model contained in 
c                           array rmodel and misfit value to 
c                           the file assigned to logical 
c                           unit lu
c
c      Note: Calls no other routines
c
c------------------------------------------------------------------------
c
      subroutine display_model(
     &            lu,imod, nd, rmodel, nsource, misfitval )
c
      include      'kine_param.inc'
c
c
      real*4            rmodel(maxmoddim)
c
      real*4            misfitval
      integer         nsource, j2
c
c                                    Write out model `rmodel`
c
      write(lu,801) imod, misfitval
c
      do j=1,nd
        write(lu,821) j, rmodel(j)
      enddo 
c         write(lu,821)
c     &  1, rmodel(1), rmodel(nsource+1), rmodel(nd) 
c        do j=2,nsource
c          j2=j+nsource      
c          write(lu,811) j, rmodel(j), rmodel(j2)
c        end do
c
  801      format( '  model:',i5,',   misfit value:',f10.5 )
  811      format( 3x,i3,f10.3 )
  821   format( i10,f15.2,f15.2,f15.2)
c
        return
        end
