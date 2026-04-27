module axi1
  implicit none
  include 'precision_module.f90'

  integer :: in1,in2,in3,out,out2
  parameter(in1=10,in2=11,in3=12,out=13,out2=14)
  
  integer,parameter :: ikmin = 100
  integer,parameter :: nkmax = 10000

  real(kind=pcal) :: explim = -600.0d0
  
  complex(kind=pcal) :: ai,pi,pi2

  real(kind=pcal) :: zero,one

end module axi1

module axi2
  implicit none
  include 'precision_module.f90'
  
  integer :: ncr,nrs,ncs,nc,nr,ns
  
  complex(kind=pcal) :: omega,omega2,a1,xlnf

  real(kind=pcal) :: kr,kr2,uconv

  real(kind=pcal),dimension(:),allocatable :: xr,yr,zr
  real(kind=pcal),dimension(:,:),allocatable :: cosr,sinr
  real(kind=pcal),dimension(:),allocatable :: xs,ys,zs
  real(kind=pcal),dimension(:),allocatable :: hc,vp,vs,vp2,vs2,rho,qp,qs
  real(kind=pcal),dimension(:),allocatable :: rr
  integer,dimension(:,:),allocatable :: irs
  integer,dimension(:),allocatable :: irc,nzr
  integer,dimension(:,:),allocatable :: nzrr
  integer,dimension(:,:,:),allocatable :: izrr
  integer,dimension(:),allocatable :: isc,nzs
  integer,dimension(:,:),allocatable :: nzss
  integer,dimension(:,:,:),allocatable :: izss

  logical :: ttconv
  
end module axi2

module axi3
  implicit none
  include 'precision_module.f90'

  complex(kind=pcal),dimension(:),allocatable :: cka,ckb,cka2,ckb2,cnu,cgam,c2
  complex(kind=pcal),dimension(:),allocatable :: cff
  
  real(kind=pcal),dimension(:),allocatable :: fj1,k0,k1,k2,k3,k4,k5
  
  complex(kind=pcal),dimension(:,:,:),allocatable :: rd,ru,td,tu
  complex(kind=pcal),dimension(:),allocatable :: rdsh,rush,tdsh,tush,me1,me2
  
  complex(kind=pcal),dimension(:,:,:),allocatable :: nt,mt
  complex(kind=pcal),dimension(:),allocatable :: ntsh,mtsh
  
  complex(kind=pcal),dimension(:,:,:),allocatable :: fdo,fup
  complex(kind=pcal),dimension(:),allocatable :: fupsh,fdosh
  
  complex(kind=pcal),dimension(:,:),allocatable :: su1,sd1,su2,sd2,su3,sd3,su4,sd4
  complex(kind=pcal),dimension(:),allocatable :: su1sh,sd1sh,su2sh,sd2sh
  
  complex(kind=pcal),dimension(:,:,:),allocatable :: u

end module axi3

module axi4
  implicit none
  include 'precision_module.f90'

  real(kind=pcal),dimension(:,:),allocatable :: jj0,jj1

end module axi4
