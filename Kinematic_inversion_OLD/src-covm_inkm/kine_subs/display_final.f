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
      subroutine display_final(
     &            lu, rmodel, nsource, nd)
c
      include 'kine_param.inc'
c
c
      real*4            rmodel(max_nd)
c
        integer         nsource, nd
c
c
c                                    Write out final model
c
      
      write(lu,801)

      do j=1,nd
        write(lu,821) j, rmodel(j)
      enddo

        
  801      format(/'*** Final model ***' /
     &          3x,'Parameter #',5x,'Value')
  821   format( i10,f15.2,f15.2,f15.2)
  811      format( i10,f15.2,f15.2)
c
      return
      end
