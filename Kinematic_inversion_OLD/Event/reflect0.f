c @(#) reflect0.F	AXITRA 4.14     9/22/97 4
c******************************************************************************
c*                                                                            *
c*                     SUBROUTINE REFLECT0                                    *
c*                                                                            *
c*    Calcul de coefficients dependant du nombre d'onde kr et de fonctions    *
c*    de Bessel.                                                              *
c*                                                                            *
c******************************************************************************







      subroutine reflect0 (ik, iklast)

c Global
      include "parameter.inc"
      include "dimension1.inc"
      include "dimension2.inc"
c Local
      integer ir,ier
      real*8    fj0,arg

      real*8    jj0(nkmax,nrsp),jj1(nkmax,nrsp)



      real*8    vy


      save        jj0,jj1
      dimension   cu(11*nrsp)
      equivalence (cu,u)

c     initialisations pour kr=0.

      if (ik.eq.1) then
      do 5 i=1,11*nrsp
      cu(i)=0.
 5    continue
      endif


c     Calcul des fonctions de Bessel J0 et J1, k1,k2,k3,k4,k5,k0

      do 10 ir=1,nrs
      arg=rr(ir)*kr

      if (ik.gt.nkmax) then


        call ff01ad(fj0,vy,arg,0)
        call ff02ad(fj1(ir),vy,arg,0)


      else
       if (ik.gt.iklast) then


        call ff01ad(jj0(ik,ir),vy,arg,0)
        call ff02ad(jj1(ik,ir),vy,arg,0)


       endif
       fj0=jj0(ik,ir)
       fj1(ir)=jj1(ik,ir)
      endif

       
      if (rr(ir).ne.0) then
       k0(ir)=kr*fj0
       k2(ir)=fj1(ir)/rr(ir)
       k1(ir)=k0(ir)-2.*k2(ir)
       k4(ir)=k1(ir)/rr(ir)
       k3(ir)=-(2.*k4(ir)+kr2*fj1(ir))
      else
c		Lorsque rr=0. il faut utiliser
c		un developpement limite
       k1(ir)=0.
       k2(ir)=kr/2.
       k3(ir)=0.
       k4(ir)=0.
      endif
      k5(ir)=k0(ir)-k2(ir)


 10   continue
      if(ik.gt.iklast) iklast=ik

c               Calcul des nombres d'onde verticaux

      do 11 ic=1,nc
      ccv=1.+ai/(qp(ic)+qs(ic))
      cvp=vp(ic)*ccv/(1.+.25/qp(ic)/qp(ic))/(1.-xlnf/qp(ic))
      cvs=vs(ic)*ccv/(1.+.25/qs(ic)/qs(ic))/(1.-xlnf/qs(ic))
      cka(ic)=omega/cvp
      ckb(ic)=omega/cvs
      ckb2(ic)=ckb(ic)*ckb(ic)
      cka2(ic)=cka(ic)*cka(ic)
      cc=cka2(ic)-kr2

      cnu(ic)=cdsqrt(cc)
      if (dimag(cnu(ic)).gt.0.d0) cnu(ic)=-cnu(ic)

      cc=ckb2(ic)-kr2

      cgam(ic)=cdsqrt(cc)
      if (dimag(cgam(ic)).gt.0.d0) cgam(ic)=-cgam(ic)

 11   continue
      do 12 ic=1,nc
      c2(ic)=kr*kr/ai/cnu(ic)
 12   continue
      
      return
      end
