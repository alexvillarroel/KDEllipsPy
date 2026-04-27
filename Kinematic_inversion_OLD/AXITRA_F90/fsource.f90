
function fsource(type,omega,t0,t1,dt)
  use spec_numbers
  implicit none
  include 'precision_subroutine.f90'
  
  integer :: type
  real(kind=pcal) :: dt,t0,t1
  real(kind=pcal) :: uur,uui,trise,trupt
  complex(kind=pcal) :: fsource,uu,uex,uxx,omega,shx
  
  ! TYPE=0 Source = Dirac en deplacement
  if (type .eq. 0) then
     fsource = cmplx(one,zero,pcal)
  endif
  
  ! TYPE=1 Source = Ricker en deplacement
  if (type .eq. 1) then
     uu = omega*t0
     uu = uu*uu/twopi/twopi
     uu = exp(-uu)
     uu = omega*omega*uu*dt
     fsource = uu
  endif
  
  ! TYPE=2 Source = step en deplacement
  ! 2 steps possibles (1) real=1/(ai*omega)
  !                   (2) bouchon's
  if (type .eq. 2) then
     shx = exp(omega*pi*t0/two)	!Bouchon's
     shx = one/(shx-one/shx)
     uu = -ai*t0*pi*shx
     fsource = uu
  endif

  ! TYPE=7 Source = step en deplacement
  if (type .eq. 7) then
     uu = one/ai/omega
     fsource = uu
  endif

  ! TYPE=3 Source = fichier 'axi.sou' sismogramme dans l'unite choisie dans le fichier
  if (type .eq. 3) then
     read(13,*) uur,uui
     fsource = cmplx(uur,uui,pcal)
  endif
 
  ! TYPE=4 Source = triangle en vitesse
  if (type .eq. 4) then
     !	  trise=t0
     !	  trupt=t0
     !	  uu=ai*omega*trise
     !          uu=(1.-exp(-uu))/uu		! ramp
     !          uxx=ai*omega*trupt/2.		! finite fault
     !          uex=exp(uxx)
     !          uxx=(uex-1./uex)/uxx/2.
     !          fsource=uu*uxx/(ai*omega)
     uu = ai*omega*t0
     fsource = 4*t0*(exp(-uu)*(exp(uu/2)-1)**2)/(uu**3)
  endif
 
  ! TYPE=5 Source = rampe causale
  ! rise time T=t0
  if (type .eq. 5) then
     trise = t0
     uu = ai*omega*trise
     uu = (one-exp(-uu))/uu
     fsource = uu/(ai*omega)
  endif
  
! TYPE=6,8        Source = modele d'haskell, trapezoide
! 1 ere cste de temps rise time: riset
! 2 eme cste de temps, duree de la rupture
! trupt = Length/2/rupt_velocity (Haskell)

  if ((type .eq. 6) .or. (type .eq. 8)) then
     trise = t0
     trupt = t1
     uu = ai*omega*trise
     uu = (one-exp(-uu))/uu		! ramp
     uxx = ai*omega*trupt/two		! finite fault
     uex = exp(uxx)
     uxx = (uex-1./uex)/uxx/two
     fsource = uu*uxx/(ai*omega)
  endif

end function fsource
        

