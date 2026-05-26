c @(#) reflect5.F	AXITRA 4.15     9/22/97 4
c***********************************************************
c*                                                         *
c*              SUBROUTINE REFLECT5                        *
c*                                                         *
c*        Calcul des deplacements avec diverses rotations  *
c*        et recombinaisons. Passage aux sources du tenseur*
c*        des moments sismiques (M1 a M6)                  *
c*        Multiplication par les termes frequentiel        *
c*        et angulaire :                                   *
c*                     u=u*C3(theta)*CFF(omega )           *
c***********************************************************




     
      subroutine reflect5

      include "parameter.inc"
      include "dimension1.inc"
      include "dimension2.inc"

      real*8   cor,cor2,co2r
      complex*16    urxx,utxx,uzxx,urxy,utxy,uzxy,urxz,utxz,uzxz,
     1           uryy,utyy,uzyy,uryz,utyz,uzyz,urzz,utzz,uzzz,
     2           ux(6),uy(6),uz(6)

      do is=1,ns
      do ir=1,nr
      do it=1,11
      u(ir,is,it)=u(ir,is,it)*a1*cff(is)
      enddo

c+++++++++++++	
c	Deplacement dus aux sources Mxx,Mxy,Mxz,Myy,Myz,Mzz avec
c	convention de signe inverse pour le deplacement vertical
c	(positif vers le haut)
c+++++++++++++

      cor=cosr(ir,is)
      sir=sinr(ir,is)
      cor2=cor*cor
      sir2=sir*sir
      co2r=cor2-sir2
      si2r=2.*cor*sir

c	Mxx
      
      urxx=-cor2 *  u(ir,is,1) - u(ir,is,2) - co2r*u(ir,is,5)
      utxx= si2r * (u(ir,is,2) + u(ir,is,6)/2.)
      uzxx= cor2 *  u(ir,is,3) + u(ir,is,4)

c	Mxy+Myx

      urxy=-si2r * (u(ir,is,1) + 2.*u(ir,is,5))
      utxy=-co2r * (2.*u(ir,is,2) + u(ir,is,6))
      uzxy= si2r *  u(ir,is,3)

c	Mxz+Mzx

      urxz=-cor * u(ir,is,7)
      utxz= sir * u(ir,is,8)
      uzxz= cor * u(ir,is,9) 

c	Myy

      uryy=-sir2 *  u(ir,is,1) - u(ir,is,2) + co2r*u(ir,is,5)
      utyy=-si2r * (u(ir,is,2) + u(ir,is,6)/2.) 
      uzyy= sir2 *  u(ir,is,3) + u(ir,is,4)

c	Myz+Mzy

      uryz=-sir * u(ir,is,7)
      utyz=-cor * u(ir,is,8)
      uzyz= sir * u(ir,is,9)

c	Mzz

      urzz=-u(ir,is,10)
      utzz= 0.
      uzzz=-u(ir,is,11)

c+++++++++++
c     Passage aux sources bis, 5 dislocations elementaires et
c     une source isotrope + rotation des composantes pour passer de 
c     radial/tangentiel a Ox/Oy
c+++++++++++

c	M1 = (Mxy + Myx)
      ux(1)=urxy*cor - utxy*sir
      uy(1)=urxy*sir + utxy*cor
      uz(1)=uzxy

c	M2 = (Mxz + Mzx)
      ux(2)=urxz*cor - utxz*sir
      uy(2)=urxz*sir + utxz*cor
      uz(2)=uzxz

c	M3 = -(Myz + Mzy)
      ux(3)=-uryz*cor + utyz*sir
      uy(3)=-uryz*sir - utyz*cor
      uz(3)=-uzyz

c	M4 = -Mxx + Mzz
      ux(4)=(-urxx + urzz)*cor - (-utxx + utzz)*sir
      uy(4)=(-urxx + urzz)*sir + (-utxx + utzz)*cor
      uz(4)=-uzxx + uzzz

c	M5 = -Myy + Mzz
      ux(5)=(-uryy + urzz)*cor - (-utyy + utzz)*sir
      uy(5)=(-uryy + urzz)*sir + (-utyy + utzz)*cor
      uz(5)=-uzyy + uzzz

c	M6 = Mxx + Myy + Mzz
      ux(6)=(urxx + uryy + urzz)*cor - (utxx + utyy + utzz)*sir
      uy(6)=(urxx + uryy + urzz)*sir + (utxx + utyy + utzz)*cor
      uz(6)=uzxx + uzyy + uzzz



      write(out2) (ux(it),it=1,6)
      write(out2) (uy(it),it=1,6)
      write(out2) (uz(it),it=1,6)

!          do it=1,6
!             if(amin.gt.abs(ux(it)))amin=abs(ux(it))
!             if(amax.lt.abs(ux(it)))amax=abs(ux(it))
!          enddo 


      enddo
      enddo
      
!      print *,'escribio axi.res  amin, amax =',amin, amax

       
      return
      end


