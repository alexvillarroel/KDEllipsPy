c------------------------------------------------------------------------
c
c     Subroutine read_vmodelNA - reads in parameters associated with
c                     velocity model and 
c                     reads in all definitions to set up 
c                     the parameter space
c
      subroutine read_vmodelNA(lu_vel, range, scales, 
     &               moddim,moddim1,moddim2)
c
        include 'fd3d_param.inc'
c
      common /fd3d_com/nt,nx,ny,weight,obs,
     &                  ndata,
     &                  lu_obs, lu_dc, fname,
     &                  nummod,num_max_mod
c     common /taille/a(2,maxmoddim),b(2,maxmoddim)

      integer          ndata

      real*4          range(2,maxmoddim)
c
      real          scales(*)
c
      read(lu_vel,*) moddim,moddim1,moddim2,scales(1)
c
      do i=1,moddim
        read(lu_vel,*) range(1,i), range(2,i), scales(i+1)
      end do
c
      return
      end

