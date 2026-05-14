c-------------------------------------------------------------------------
c
c      Subroutine summary - writes out a summary of the current
c                       iteration of the Neighbourhood algorithm 
c                       together with the model `rmodel.
c
c--------------------------------------------------------------------------
c
      subroutine output_summary(
     &            lu_out, lu_sum, it, rmodel, moddim, ntot, 
     &            mfitmin, mfitmean, mfitminc, mopt)
c
        include 'kine_param.inc'
c
c
      real*4            rmodel(*)
c
        real*4          mfitmin,
     &                  mfitmean,
     &                  mfitminc
      logical            lw
        integer         nsource
c
      lw = .true.
      if(lu_out.eq.0)lw = .false.
      write(lu_sum,801)
      if(lw)write(lu_out,801)
c
      write(lu_sum,811) it,ntot,mfitmin,mfitmean,
     &                    mfitminc,mopt
      if(lw)write(lu_out,811) it,ntot,mfitmin,mfitmean,
     &                    mfitminc,mopt
c


      do j=1,nd
        write(lu_sum,812) j, rmodel(j) 
      enddo


c        write(lu_sum,821)
c     &  1, rmodel(1), rmodel(nsource+1), rmodel(moddim) 
c      do j=2,nsource
c        j2=nsource+j
c        write(lu_sum,812) j, rmodel(j), rmodel(j2) 
c          if(lw)write(lu_out,812) j, rmodel(j), rmodel(j2)
c      end do
c
      write(lu_sum,*)
      if(lw)write(lu_out,*)
c
  801      format( 3x,'It',1x,'Nsampled',3x,'Mfitmin',
     &          2x,'Mfitmean',1x,'Mfitmeanc',1x,
     &             3x,'Mopt')
  821 format( 5x,i3,3f10.3)
  811 format( i5,i9,3f10.5,1x,i7)
  812 format( 5x,i3,4f10.3,1x,i7)
c 
      return
      end
