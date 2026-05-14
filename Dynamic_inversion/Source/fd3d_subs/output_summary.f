c Modif S.Peyrat 2002-2003
c
c-------------------------------------------------------------------------
c
c	Subroutine summary - writes out a summary of the current
c			     iteration of the Neighbourhood algorithm 
c			     together with the model `rmodel'.
c
c--------------------------------------------------------------------------
c
	subroutine output_summary(
     &		lu_out, lu_sum, it, rmodel, moddim, ntot, 
     &		mfitmin, mfitmean, mfitminc, mopt)
c
        include 'fd3d_param.inc'
c
c
	real*4		rmodel(*)
c
        real*4          mfitmin,
     &                  mfitmean,
     &                  mfitminc
	logical		lw
c
c						Write out headers for files
c						with summary of optimization
c						performance.
	lw = .true.
	if(lu_out.eq.0)lw = .false.
c	write(lu_sum,801)
	if(lw)write(lu_out,801)
c
c 	write(lu_sum,811) it,ntot,mfitmin,mfitmean,
c     &                    mfitminc,mopt
 	if(lw)write(lu_out,811) it,ntot,mfitmin,mfitmean,
     &                    mfitminc,mopt
c
	do j=1,moddim
c	  write(lu_sum,812) j, rmodel(j)
	  if(lw)write(lu_out,812) j, rmodel(j)
	end do
c
c	write(lu_sum,*)
	if(lw)write(lu_out,*)
c
  801	format( 3x,'It',1x,'Nsampled',3x,'Mfitmin',
     &          2x,'Mfitmean',1x,'Mfitmeanc',1x,
     &             3x,'Mopt')
  811	format( i5,i9,3e14.6,1x,i7)
  812	format( 5x,i3,4e14.6,1x,i7)
c 
	return
	end
