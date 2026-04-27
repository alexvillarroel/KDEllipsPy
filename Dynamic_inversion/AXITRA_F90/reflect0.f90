
!******************************************************************************
!*                                                                            *
!*                     SUBROUTINE REFLECT0                                    *
!*                                                                            *
!*    Calcul de coefficients dependant du nombre d'onde kr et de fonctions    *
!*    de Bessel.                                                              *
!*                                                                            *
!******************************************************************************

subroutine reflect0 (ik,iklast)
  use axi1
  use axi2
  use axi3
  use axi4
  implicit none
  include 'precision_subroutine.f90'

  integer :: ir,ik,iklast,ic
  real(kind=pcal) :: fj0,arg

  real(kind=pcal) :: vy

  complex(kind=pcal) :: ccv,cvp,cvs,cc

  ! initialisations pour kr=zero

  if (ik .eq. 1) then
     u(:,:,:) = zero
  endif

  ! calcul des fonctions de Bessel J0,J1,k1,k2,k3,k4,k5,k0

  do ir = 1,nrs
     
     arg = rr(ir)*kr
     
     if (ik .gt. nkmax) then
        call ff01ad(fj0,vy,arg,0)
        call ff02ad(fj1(ir),vy,arg,0)
     else
        if (ik .gt. iklast) then
           call ff01ad(jj0(ik,ir),vy,arg,0)
           call ff02ad(jj1(ik,ir),vy,arg,0)
        endif
        fj0 = jj0(ik,ir)
        fj1(ir) = jj1(ik,ir)
     endif

     if (rr(ir) .ne. 0) then
        k0(ir) = kr*fj0
        k2(ir) = fj1(ir)/rr(ir)
        k1(ir) = k0(ir)-2.*k2(ir)
        k4(ir) = k1(ir)/rr(ir)
        k3(ir) = -(2*k4(ir)+kr2*fj1(ir))
     else
        ! lorsque rr=zero il faut utiliser un developpement limite
        k1(ir) = zero
        k2(ir) = kr/2
        k3(ir) = zero
        k4(ir) = zero
     endif
     k5(ir) = k0(ir)-k2(ir)

  end do

  if (ik .gt. iklast) iklast = ik

  ! calcul des nombres d'onde verticaux

  do ic = 1,nc
     ccv = one+ai/(qp(ic)+qs(ic))
     cvp = vp(ic)*ccv/(one+one/4/qp(ic)/qp(ic))/(one-xlnf/qp(ic))
     cvs = vs(ic)*ccv/(one+one/4/qs(ic)/qs(ic))/(one-xlnf/qs(ic))
     cka(ic) = omega/cvp
     ckb(ic) = omega/cvs
     ckb2(ic) = ckb(ic)*ckb(ic)
     cka2(ic) = cka(ic)*cka(ic)
     cc = cka2(ic)-kr2

     cnu(ic) = sqrt(cc)
     if (aimag(cnu(ic)) .gt. zero) cnu(ic) = -cnu(ic)

     cc = ckb2(ic)-kr2

     cgam(ic) = sqrt(cc)
     if (aimag(cgam(ic)) .gt. zero) cgam(ic) = -cgam(ic)

  end do

  do ic = 1,nc
     c2(ic) = kr*kr/ai/cnu(ic)
  end do

end subroutine reflect0
