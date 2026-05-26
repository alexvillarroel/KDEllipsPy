module spec_numbers
  implicit none
  include 'precision_module.f90'
  real(kind=pcal) :: zero,one,two,pi,twopi
  complex(kind=pcal) :: ai,zerocmplx
end module spec_numbers

program convm
  use spec_numbers
  implicit none
  include 'precision_subroutine.f90'

  integer :: icc,ics,ic,id,icu,icv,ica
  
  real(kind=pcal) :: t0,dt0,rpvel
  
  integer :: index,indexin,ifr
  
  real(kind=pcal) :: xphi
  
  character(len=50) :: filein,histfile,axires,convmfile
  character(len=10) :: sourcefile,statfile

  integer,dimension(:),allocatable :: iwk

  integer,dimension(:),allocatable :: isc
  integer :: jf,ir,is,it,nc,ns,nr,nfreq,ikmax,mm,nt
  real(kind=pcal),dimension(:),allocatable :: hc,vp,vs,rho,qp,qs
  real(kind=pcal),dimension(:),allocatable :: mu,rvel,delay,strike,dip,rake,disp,length,width,xs,ys,zs
  real(kind=pcal),dimension(:,:),allocatable :: a
  real(kind=pcal),dimension(:),allocatable :: xr,yr,zr
  real(kind=pcal),dimension(:),allocatable :: uout

  real(kind=pcal) :: tl,xl,uconv,hh,zsc,dfreq,freq,aw,ck,xmm,cc,pas,t1
  
  complex(kind=pcal),dimension(:,:),allocatable :: ux,uy,uz
  complex(kind=pcal) :: omega
  complex(kind=pcal),dimension(1:6) :: uxf,uyf,uzf
  
  complex(kind=pcal),dimension(:,:,:,:),allocatable :: fsu

  complex(kind=pcal) :: deriv,fsource,us,uux,uuy,uuz
  
  namelist /input/ nc,nfreq,tl,aw,nr,ns,xl,ikmax,uconv,sourcefile,statfile
  
  two = real(2.0d0,pcal)
  one = real(1.0d0,pcal)
  zero = real(0.0d0,pcal)
  zerocmplx = cmplx(zero,zero,pcal)
  ai = cmplx(zero,one,pcal)

  pi = real(4*atan(1.0d0),pcal)
  twopi = 2*pi	

  dt0 = zero
  
  filein = 'axi.head'
  axires = 'axi.res'
  histfile = 'axi.hist'
  convmfile = 'convm.data'

  open(12,form='unformatted',file=axires)
  open(10,form='formatted',file=filein)
  read(10,input)

  allocate(hc(1:nc),vp(1:nc),vs(1:nc),rho(1:nc),qp(1:nc),qs(1:nc))
  allocate(mu(1:ns),rvel(1:ns),delay(1:ns))
  allocate(strike(1:ns),dip(1:ns),rake(1:ns),disp(1:ns),length(1:ns),width(1:ns))
  allocate(xs(1:ns),ys(1:ns),zs(1:ns))
  allocate(a(1:6,1:ns))
  allocate(xr(1:nr),yr(1:nr),zr(1:nr))
  allocate(isc(1:ns))

  open(13,form='formatted',file=sourcefile)
  open(14,form='formatted',file=statfile)
  open(15,form='formatted',file=histfile)
  open(16,form='formatted',file=convmfile)
  
  rpvel = one

  read(16,*) ics
  read(16,*) icu,icv,ica

  select case(ics)
  case(0) ! dirac
  case(1) ! ricker
     read(16,*) t0 ! Ricker pseudo period
  case(2) ! step
     read(16,*) t0 ! source duration
  case(3) ! function stored in file axi.sou
  case(4) ! triangle
     read(16,*) t0 ! triangle width
  case(5) ! ramp
     read(16,*) t0 ! rise time
  case(6) ! haskell not available
     !read(16,*) t0,repvel ! rise time, rupture velocity(fract. of S wave eg. 0.8)
  case(7) ! true step
  case(8) ! trapezoid
     read(16,*) t0,t1 ! rise time, source process time t1 (total time = t1+t0)
     t1 = t1+t0
  case(9) ! green functions for each source
  end select

  do icc = 0,2

     if (((icc .eq. 0) .and. (icu .eq. 1)) .or. ((icc .eq. 1) .and. (icv .eq. 1)) .or. ((icc .eq. 2) .and. (ica .eq. 1))) then
  
        rewind(10)
        read(10,input)
        rewind(12)
        rewind(13)
        rewind(14)
        rewind(15)

        read(10,*)
        do ic = 1,nc
           read(10,*) hc(ic),vp(ic),vs(ic),rho(ic),qp(ic),qs(ic)
        end do
        read(10,*)
        do is = 1,ns
           read(13,*) index,xs(is),ys(is),zs(is)
           indexin = -1
           rewind(15)
           do while (indexin .ne. index)
              read(15,*) indexin,disp(is),strike(is),dip(is),rake(is),width(is),length(is),delay(is)
           enddo
           delay(is) = delay(is)+dt0
        enddo
  
        do ir = 1,nr
           read(14,*) xr(ir),yr(ir),zr(ir)
        enddo
        
        if (hc(1) .eq. zero) then
           do ic = 1,nc-1
              hc(ic) = hc(ic+1)-hc(ic)
           enddo
        endif
        
        do is = 1,ns
           hh = zero
           isc(is) = 1
           zsc = zs(is)
           do ic = 1,nc-1
              hh = hc(ic)
              if (zsc .gt. hh) then
                 zsc = zsc-hh
                 isc(is) = ic+1
              else
                 goto 91
              endif
           enddo
           
91         continue
           
           mu(is) = vs(isc(is))*vs(isc(is))*rho(isc(is))
           
           call cmoment (mu(is),strike(is),dip(is),rake(is),disp(is),width(is)*length(is),a(:,is))
           
        end do
        
        xmm = log(real(nfreq,pcal))/log(two)
        mm = int(xmm)+2
        nt = 2**mm
        
        allocate(iwk(1:nt))
        allocate(uout(1:nt))
        allocate(ux(1:nt,1:nr),uy(1:nt,1:nr),uz(1:nt,1:nr))
        
        
        if (ics .ne. 9) then
           
           ux(:,:) = zerocmplx
           uy(:,:) = zerocmplx
           uz(:,:) = zerocmplx
           
           dfreq = one/tl
           pas = tl/nt
           aw = -pi*aw/tl
           
           do jf = 1,nfreq
              
              freq = real(jf-1,pcal)/tl
              omega = cmplx(twopi*freq,aw,pcal)
              deriv = (ai*omega)**icc
              
              do is = 1,ns
                 
                 xphi = atan2(yr(1)-ys(is),xr(1)-xs(is))
                 rvel(is) = rpvel * vs(isc(is))/(one-rpvel*cos(xphi-strike(is)))
                 
                 if (ics.ne.8) then
                    t1 = length(is)/rvel(is)
                 end if
                 
                 us = fsource(ics,omega,t0,t1,pas)*deriv*exp(-ai*omega*delay(is))
                 
                 do ir = 1,nr
                    read(12)(uxf(it),it=1,6)
                    read(12)(uyf(it),it=1,6)
                    read(12)(uzf(it),it=1,6)           
                    uux = zero
                    uuy = zero
                    uuz = zero
                    do it = 1,6
                       uux = uux + uxf(it)*a(it,is)
                       uuy = uuy + uyf(it)*a(it,is)
                       uuz = uuz + uzf(it)*a(it,is)
                    enddo
                    ux(jf,ir) = ux(jf,ir) + uux * us 
                    uy(jf,ir) = uy(jf,ir) + uuy * us
                    uz(jf,ir) = uz(jf,ir) + uuz * us
                 end do
              end do
              
           end do
           
           if (icc.eq.0) then
              open (52,form='formatted',file='disp_x')
              open (53,form='formatted',file='disp_y')
              open (54,form='formatted',file='disp_z')
           endif
           if (icc.eq.1) then
              open (52,form='formatted',file='velo_x')
              open (53,form='formatted',file='velo_y')
              open (54,form='formatted',file='velo_z')
           endif
           if (icc.eq.2) then
              open (52,form='formatted',file='acce_x')
              open (53,form='formatted',file='acce_y')
              open (54,form='formatted',file='acce_z')
           endif
           
           
           do ir = 1,nr
              
              do ifr = nt+2-nfreq,nt
                 ux(ifr,ir) = conjg(ux(nt+2-ifr,ir))
                 uy(ifr,ir) = conjg(uy(nt+2-ifr,ir))
                 uz(ifr,ir) = conjg(uz(nt+2-ifr,ir))
              enddo
              
              call fft2cd(nt,ux(:,ir),mm,iwk(:))
              call fft2cd(nt,uy(:,ir),mm,iwk(:))
              call fft2cd(nt,uz(:,ir),mm,iwk(:))
              
              do it = 1,nt
                 ck = real(it-1,pcal)/nt
                 cc = exp(-aw*tl*ck)/tl
                 ux(it,ir) = ux(it,ir)*real(cc)
                 uy(it,ir) = uy(it,ir)*real(cc)
                 uz(it,ir) = uz(it,ir)*real(cc)
              enddo
              
              do it = 1,nt
                 write(52,'(e,a,$)') real(ux(it,ir)),' '
                 write(53,'(e,a,$)') real(uy(it,ir)),' '
                 write(54,'(e,a,$)') real(uz(it,ir)),' '
              end do
              write(52,*)
              write(53,*)
              write(54,*)
              
           end do
           
           close(52)
           close(53)
           close(54)
           
        else
           
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
           
           dfreq = one/tl
           pas = tl/nt
           aw = -pi*aw/tl
           
           do jf = 1,nfreq
              
              freq = real(jf-1,pcal)/tl
              omega = cmplx(twopi*freq,aw,pcal)
              deriv = (ai*omega)**icc
              
              us = deriv
              
              do is = 1,ns
                 
                 do ir = 1,nr
                    read(12)(uxf(it),it=1,6)
                    read(12)(uyf(it),it=1,6)
                    read(12)(uzf(it),it=1,6)
                    uux = zero
                    uuy = zero
                    uuz = zero
                    do it = 1,6
                       uux = uux + uxf(it)*a(it,is)
                       uuy = uuy + uyf(it)*a(it,is)
                       uuz = uuz + uzf(it)*a(it,is)
                    enddo
                    fsu(jf,is,1,ir) = fsu(jf,is,1,ir) + uux * us 
                    fsu(jf,is,2,ir) = fsu(jf,is,2,ir) + uuy * us
                    fsu(jf,is,3,ir) = fsu(jf,is,3,ir) + uuz * us
                 end do
              end do
              
           end do
           
           do ir = 1,nr
              do id = 1,3
                 do is = 1,ns
                    
                    do ifr = nt+2-nfreq,nt
                       fsu(ifr,is,id,ir) = conjg(fsu(nt+2-ifr,is,id,ir))
                    enddo
                    
                    call fft2cd(nt,fsu(:,is,id,ir),mm,iwk(:))
                    
                    do it = 1,nt
                       ck = real(it-1,pcal)/nt
                       cc = exp(-aw*tl*ck)/tl
                       fsu(it,is,id,ir) = fsu(it,is,id,ir)*real(cc)
                    end do
                    
                 end do
              end do
           end do
           
           do ir = 1,nr
              do id = 1,3
                 do is = 1,ns
                    write(60) real(fsu(:,is,id,ir),pcal)
                 end do
              end do
           end do
           
           close(60)
           
           deallocate(fsu)

        end if

        deallocate(iwk)
        deallocate(uout)
        deallocate(ux,uy,uz)
        
     end if
        
  end do
     
end program convm


subroutine cmoment(mu,strike,dip,rake,disp,surf,a)
  use spec_numbers
  implicit none
  include 'precision_subroutine.f90'
  
  real(kind=pcal),dimension(1:6) :: a
  real(kind=pcal) :: xmoment,mu,strike,dip,rake,disp,surf,sd,cd,sp,cp,s2p,s2d,c2p,c2d,x1,x2,x3,x4,x5,x6,cl,sl

  if (surf .eq. zero) then
     xmoment = disp
  else
     xmoment = mu*disp*surf
  endif
  strike = strike*pi/180
  dip = dip*pi/180
  rake = rake*pi/180
  sd = sin(dip)
  cd = cos(dip)
  sp = sin(strike)
  cp = cos(strike)
  sl = sin(rake)
  cl = cos(rake)
  s2p = 2*sp*cp
  s2d = 2*sd*cd
  c2p = cp*cp-sp*sp
  c2d = cd*cd-sd*sd
  
  x1 = -(sd*cl*s2p + s2d*sl*sp*sp)*xmoment
  x2 = (sd*cl*c2p + s2d*sl*s2p/2.)*xmoment
  x3 = -(cd*cl*cp  + c2d*sl*sp)*xmoment
  x4 = (sd*cl*s2p - s2d*sl*cp*cp)*xmoment
  x5 = -(cd*cl*sp  - c2d*sl*cp)*xmoment
  x6 =  (s2d*sl)*xmoment
  
  a(1) = x2
  a(2) = x3
  a(3) = -x5
  a(4) = (-2*x1 + x4 + x6)/3
  a(5) = (x1 -2*x4 + x6)/3
  a(6) = zero
  
  if (surf .eq. -one) then
     xmoment = disp
     a(1) = zero
     a(2) = zero
     a(3) = zero
     a(4) = zero
     a(5) = zero
     a(6) = xmoment
  elseif (surf .eq. -2*one) then
  endif

end subroutine cmoment

