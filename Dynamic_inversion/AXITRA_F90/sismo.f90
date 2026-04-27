module spec_numbers
  implicit none
  include 'precision_module.f90'
  real(kind=pcal) :: zero,one,two,pi
end module spec_numbers

program convm
  use spec_numbers
  implicit none
  include 'precision_subroutine.f90'

  integer :: icc,icu,icv,ica

  integer :: rec20

  integer :: id
  
  integer :: index,indexin
  
  character(len=50) :: filein,histfile,convmfile
  character(len=10) :: sourcefile,statfile

  integer :: ir,is,it,it1,nc,ns,nr,nfreq,ikmax,mm,nt
  real(kind=pcal),dimension(:),allocatable :: mu,rvel,delay,strike,dip,rake,disp,length,width,xs,ys,zs

  real(kind=pcal) :: tl,xl,aw,uconv,xmm,pas
  
  real(kind=pcal),dimension(:,:,:,:),allocatable :: fsu
  real(kind=pcal),dimension(:,:),allocatable :: fs
  real(kind=pcal),dimension(:,:,:),allocatable :: u

  namelist /input/ nc,nfreq,tl,aw,nr,ns,xl,ikmax,uconv,sourcefile,statfile
  
  two = real(2.0d0,pcal)
  one = real(1.0d0,pcal)
  zero = real(0.0d0,pcal)
  pi = real(4*atan(1.0d0),pcal)

  histfile = "axi.hist"
  filein = 'axi.head'
  convmfile = 'convm.data'
  open(10,form='formatted',file=filein)
  read(10,input)

  allocate(mu(1:ns),rvel(1:ns),delay(1:ns))
  allocate(strike(1:ns),dip(1:ns),rake(1:ns),disp(1:ns),length(1:ns),width(1:ns))
  allocate(xs(1:ns),ys(1:ns),zs(1:ns))

  open(13,form='formatted',file=sourcefile)
  open(15,form='formatted',file=histfile)
  open(16,form='formatted',file=convmfile)

  xmm = log(real(nfreq,pcal))/log(two)
  mm = int(xmm)+2
  nt = 2**mm

  pas = tl/nt

  allocate(fs(1:nt,1:ns))
  
  open(20,file='func_sources',form='unformatted',access='direct',recl=pcal)

  do is = 1,ns
     read(13,*) index,xs(is),ys(is),zs(is)
     indexin = -1
     rec20 = -nt
     rewind(15)
     do while (indexin .ne. index)
        rec20 = rec20 + nt
        read(15,*) indexin,disp(is),strike(is),dip(is),rake(is),width(is),length(is),delay(is)
     enddo
     do it = 1,nt
        read(20,rec=rec20+it) fs(it,is)
     end do
  enddo

  close(20)
  
  read(16,*)
  read(16,*) icu,icv,ica

  do icc = 0,2

     if (((icc .eq. 0) .and. (icu .eq. 1)) .or. ((icc .eq. 1) .and. (icv .eq. 1)) .or. ((icc .eq. 2) .and. (ica .eq. 1))) then

        select case(icc)
        case(0)
           open(60,file='convm.green_u',form='unformatted')
        case(1)
           open(60,file='convm.green_v',form='unformatted')
        case(2)
           open(60,file='convm.green_a',form='unformatted')
        end select
        
        allocate(fsu(1:nt,1:ns,1:3,1:nr))
        fsu(:,:,:,:) = zero
        
        do ir = 1,nr
           do id = 1,3
              do is = 1,ns
                 read(60) fsu(:,is,id,ir)
              end do
           end do
        end do
        
        close(60)
        
        allocate(u(1:nt,1:3,1:nr))
        u(:,:,:) = zero
        
        do ir = 1,nr
           do id = 1,3
              do it = 1,nt
                 do is = 1,ns
                    do it1 = 1,it
                       u(it,id,ir) = u(it,id,ir) + fsu(it-it1+1,is,id,ir)*fs(it1,is)
                    end do
                 end do
              end do
           end do
        end do
        
        select case(icc)
        case(0)
           open(51,form='formatted',file='disp_x')
           open(52,form='formatted',file='disp_y')
           open(53,form='formatted',file='disp_z')
        case(1)
           open(51,form='formatted',file='velo_x')
           open(52,form='formatted',file='velo_y')
           open(53,form='formatted',file='velo_z')
        case(2)
           open(51,form='formatted',file='acce_x')
           open(52,form='formatted',file='acce_y')
           open(53,form='formatted',file='acce_z')
        end select
        
        do ir = 1,nr
           do it = 1,nt
              write(51,'(e,a,$)') u(it,1,ir)*pas,' '
              write(52,'(e,a,$)') u(it,2,ir)*pas,' '
              write(53,'(e,a,$)') u(it,3,ir)*pas,' '
           end do
           write(51,*)
           write(52,*)
           write(53,*)
        end do
        
        close(51)
        close(52)
        close(53)

        deallocate(fsu)
        deallocate(u)

     end if

  end do
  
  open(54,form='formatted',file='sismo_t')
  do it = 1,nt
     write(54,'(e,a,$)') it*pas,' '
  end do
  close(54)

end program convm

