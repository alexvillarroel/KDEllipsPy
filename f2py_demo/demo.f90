! demo.f90 — mini-ejemplo de f2py (sandbox de aprendizaje)
! ----------------------------------------------------------------------
! Dos subrutinas para ver el flujo Fortran -> .so -> Python:
!   1) vecadd      : "hola mundo" — suma de dos vectores.
!   2) ellipse_patch : mini-mkstress — rellena una malla 2D con una
!                      elipse de alto valor (anticipo de la Fase 1).
!
! Las líneas "!f2py intent(...)" le dicen a f2py qué entra y qué sale.
! (En el Fortran fijo F77 del legacy, lo mismo se escribe "cf2py ...").
! ----------------------------------------------------------------------

      subroutine vecadd(a, b, c, n)
      implicit none
      integer n
      real*8 a(n), b(n), c(n)
!f2py intent(in)  :: a, b
!f2py intent(out) :: c
!f2py intent(hide):: n          ! n se deduce del tamaño de a
      integer i
      do i = 1, n
         c(i) = a(i) + b(i)
      end do
      end subroutine vecadd


      subroutine ellipse_patch(nx, ny, ax, bx, strin, strout, field)
!     Rellena field(nx,ny): 'strin' dentro de una elipse centrada en la
!     malla con semiejes (ax,bx) [en puntos], 'strout' fuera. Es la
!     versión bebé de mkstress.f del solver dinámico.
      implicit none
      integer nx, ny
      real*8 ax, bx, strin, strout
      real*8 field(nx, ny)
!f2py intent(in)  :: nx, ny, ax, bx, strin, strout
!f2py intent(out) :: field
      integer i, j
      real*8 x0, y0, rr
      x0 = 0.5d0 * (nx + 1)
      y0 = 0.5d0 * (ny + 1)
      do j = 1, ny
         do i = 1, nx
            rr = ((i - x0) / ax)**2 + ((j - y0) / bx)**2
            if (rr .lt. 1.0d0) then
               field(i, j) = strin
            else
               field(i, j) = strout
            end if
         end do
      end do
      end subroutine ellipse_patch
