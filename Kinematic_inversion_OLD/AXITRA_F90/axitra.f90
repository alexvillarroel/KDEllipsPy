
!******************************************************************************
!*									      *
!*			PROGRAMME AXITRA                                      *
!*									      *
!*   	Calcul de sismogrammes synthetiques en milieu stratifie a symetrie    *
!*      cylindrique.                                                          *
!*	Propagation par la methode de la reflectivite, avec coordonnees       *
!*      cylindriques (r, theta, z)				              *
!*      Attenuation sur les ondes P et S                                      *
!*									      *
!*      auteur : Coutant O. reecriture Favreau P.                             * 
!*	Bibliographie :                                                       *
!*                      Kennett GJRAS vol57, pp557R, 1979                     *
!*			Bouchon JGR vol71, n4, pp959, 1981                    *
!*									      *
!******************************************************************************

! les fonctions mathematiques usuelles doivent normalement calculer avec la precision de l'argument
! les fonctions ff01a et ff02a sont limitees a la double precision
! la fonction cmplx est appelee avec son 3eme argument pour eviter la conversion en complexe de precision inferieure
! pour chaque compilateur, verifier que les fonctions real et aimag ne convertissent pas en reel de precision inferieure
! allocation dynamique avec nr,nc,ns (nrp,ncp,nsp n'existent plus !)
! semble supporter les options de debuggage pointilleuses et les options d'optimisation agressives

program axitra
  use axi1
  use axi2
  use axi4
  implicit none
  include 'precision_subroutine.f90'
  
  character(len=10) :: sourcefile,statfile
  integer :: ic,ir,is
  real(kind=pcal) :: dfreq,freq,pil,zom,rw,phi
  logical,dimension(:,:),allocatable :: tconv
  
  integer :: ifreq,nfreq,ik,ikmax,iklast,lastik
  real(kind=pcal) :: xl,aw,tl,fr1,fr2

  namelist /input/ nc,nfreq,fr1,fr2,tl,aw,nr,ns,xl,ikmax,uconv,sourcefile,statfile
  
  zero = real(0.0q0,pcal)
  one = real(1.0q0,pcal)

  pi = real(4*atan(1.0q0),pcal)
  pi2 = 2*pi

  ai = cmplx(zero,one,pcal)

  open (in1,form='formatted',file='axi.data')
  
  open (out,form='formatted',file='axi.head')
  open (out2,form='unformatted',file='axi.res')
  
  rewind(out)
  rewind(out2)
  
  !c++++++++++
  !c           LECTURE DES PARAMETRES D'ENTREE
  !c              
  !c               sismogramme : nfreq,tl,xl
  !c               recepteurs  : nr,xr(),yr(),zr()
  !c               source      : xs,ys,zs
  !c               modele      : nc,hc(),vp(),vs(),rho()
  !c                            
  !c               si hc(1)=0 on donne les profondeurs des interfaces, sinon
  !c               on donne les epaisseurs des couches
  !c++++++++++
  
  read(in1,input)

  allocate(xr(1:nr),yr(1:nr),zr(1:nr))
  allocate(cosr(1:nr,1:ns),sinr(1:nr,1:ns))
  allocate(xs(1:ns),ys(1:ns),zs(1:ns))
  allocate(rr(1:nr*ns))
  allocate(irs(1:nr,1:ns))
  allocate(hc(1:nc),vp(1:nc),vs(1:nc),vp2(1:nc),vs2(1:nc),rho(1:nc),qp(1:nc),qs(1:nc))
  allocate(irc(nc),nzr(nc))
  allocate(nzrr(nr,nc))
  allocate(izrr(nr,nr,nc))
  allocate(isc(nc),nzs(nc))
  allocate(nzss(ns,nc))
  allocate(izss(ns,ns,nc))
  allocate(tconv(nr,ns))

  do ic = 1,nc
     read(in1,*) hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)
  enddo
  open (in2,form='formatted',file=sourcefile)
  open (in3,form='formatted',file=statfile)
  
  write(out,input)
  write(out,*) 'hc,vp,vs,rho,Qp,Qs'
  do ic = 1,nc
    write(out,'(8f9.3)') hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)
  end do

  call flush(out)

  !c++++++++++
  !c           INITIALISATIONS
  !c++++++++++

  uconv = uconv*uconv
  dfreq = one/tl
  aw = -pi*aw/tl
  freq = -dfreq
  pil = pi2/xl
  iklast = 0
  
  call initdata
  
  !c               ***************************
  !c               ***************************
  !c               **  BOUCLE EN FREQUENCE  **
  !c               ***************************
  !c               ***************************
  
  allocate(jj0(1:nkmax,1:nr*ns),jj1(1:nkmax,1:nr*ns))

  do ifreq = 1,nfreq
     
     freq = freq+dfreq
     rw = pi2*freq
     omega = cmplx(rw,aw,pcal)
     omega2 = omega*omega
     a1 = (one/omega2/xl)/2
     zom = sqrt(rw*rw+aw*aw)
     if (ifreq .eq. 1) then
        phi = -pi/2
     else
        phi = atan(aw/rw)
     endif
     do ir = 1,nr
        do is = 1,ns
  	   tconv(ir,is) = .false.
        enddo
     enddo
     ttconv = .false.
     
     xlnf = (ai*phi+log(zom))/pi
     
     
     !c            ******************************************
     !c            ******************************************
     !c            **  RESOLUTION PAR BOUCLE EXTERNE EN Kr **
     !c            ******************************************
     !c            ******************************************
     
     do ik = 0,ikmax
        
        kr = (ik+real(0.258q0,pcal))*pil
        kr2 = kr*kr
        
        
        !c+++++++++++++
        !c              Calcul de nombreux coefficients et des fonctions de Bessel
        !c+++++++++++++
        
        call reflect0 (ik+1,iklast)
        
        
        !c+++++++++++++
        !c              Calcul des coefficients de reflexion/transmission
        !c	       Matrice de Reflection/Transmission et Dephasage
        !c+++++++++++++
        
        call reflect1
        
        !c+++++++++++++
        !c              Calcul des matrices de reflectivite : mt(),mb(),nt(),nb()
        !c              (rapport des differents potentiels montant/descendant
        !c                        en haut et en bas de chaque couche)
        !c+++++++++++++
        
        call reflect2
        
        !c+++++++++++++
        !c	       Calcul des matrices de passage des vecteurs potentiel 
        !c		source, aux vecteurs potentiel PHI, PSI et KHI au sommet
        !c		de chaque couche
        !c+++++++++++++
        
        call reflect3
        
        !c+++++++++++++
        !c	       Calcul des potentiels et des deplacement dus aux sources du
        !c		tenseur, en chaque recepteur (termes en kr, r, z)
        !c+++++++++++++
        
        call reflect4 (ik .gt. ikmin,tconv)
        
        if (ttconv) goto 21
	
     end do
     
21   continue
     
     !c+++++++++++++
     !c               Calcul des deplacements aux recepteurs 
     !c		Sortie des resultats
     !c+++++++++++++
     
     lastik = ik-1
     write(out,*) 'freq =',freq,'iter =',lastik
     call flush(out)
     if (ifreq .eq. 1) lastik = 0
     
     call reflect5
     
     if (ik .ge. ikmax) then
        write(6,*) 'Depassement du nombre d iteration maximum'
        stop
     endif
     
  end do

end program axitra
