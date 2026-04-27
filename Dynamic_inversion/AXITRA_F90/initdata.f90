!******************************************************************************
!*                                                                            *
!*                     SUBROUTINE INITDATA                                    *
!*                                                                            *
!*    Initialisation de divers parametres                                     *
!*                                                                            *
!*    Input:
!*	hc,zr,zs,nc,nr
!*    Output:
!*     	ncr,irc,nzr,irzz,nzrr,rr
!*     	ncs,isc,nzs,iszz,nzss,rr,iss
!*    Modified:
!*   	hc
!******************************************************************************

subroutine initdata
  use axi1
  use axi2
  use axi3
  implicit none
  include 'precision_subroutine.f90'
  
  integer :: ir,ir1,ir2,ic,jr,jrr,js,jss,is,is1,is2
  integer,dimension(1:ns) :: index
  logical :: tc
  real(kind=pcal) :: hh,tmp
  real(kind=pcal),dimension(1:nr,1:ns) :: r

  allocate(cka(1:nc),ckb(1:nc),cka2(1:nc),ckb2(1:nc),cnu(1:nc),cgam(1:nc),c2(1:nc))
  allocate(cff(1:ns))

  allocate(fj1(1:nr*ns),k0(1:nr*ns),k1(1:nr*ns),k2(1:nr*ns),k3(1:nr*ns),k4(1:nr*ns),k5(1:nr*ns))

  allocate(rd(1:nc,1:2,1:2),ru(1:nc,1:2,1:2),td(1:nc,1:2,1:2),tu(1:nc,1:2,1:2))
  allocate(rdsh(1:nc),rush(1:nc),tdsh(1:nc),tush(1:nc),me1(1:nc),me2(1:nc))

  allocate(nt(1:nc,1:2,1:2),mt(1:nc,1:2,1:2))
  allocate(ntsh(1:nc),mtsh(1:nc))

  allocate(fdo(1:nc,1:2,1:2),fup(1:nc,1:2,1:2))
  allocate(fupsh(1:nc),fdosh(1:nc))

  allocate(su1(1:ns,1:2),sd1(1:ns,1:2),su2(1:ns,1:2),sd2(1:ns,1:2),su3(1:ns,1:2),sd3(1:ns,1:2),su4(1:ns,1:2),sd4(1:ns,1:2))
  allocate(su1sh(1:ns),sd1sh(1:ns),su2sh(1:ns),sd2sh(1:ns))

  allocate(u(1:nr,1:ns,1:11))


  !++++++++++++
  !        Lecture coordonnees stations et recepteurs
  !++++++++++++

  do is = 1,ns
     read(in2,*) index(is),xs(is),ys(is),zs(is)
  enddo
  do ir = 1,nr
     read(in3,*) xr(ir),yr(ir),zr(ir)
  enddo


  !++++++++++++
  !        conversion interface -> epaisseur des couches     
  !++++++++++++

  if (hc(1) .eq. zero) then
     do ic = 1,nc-1
        hc(ic) = hc(ic+1)-hc(ic)
     enddo
  endif

  !++++++++++++
  !        on reordonne les sources par profondeur croissante
  !++++++++++++

  do is1 = 1,ns-1
     do is2 = is1,ns
        if (zs(is1) .gt. zs(is2)) then
           tmp = xs(is1)
           xs(is1) = xs(is2)
           xs(is2) = tmp
           tmp = ys(is1)
           ys(is1) = ys(is2)
           ys(is2) = tmp
           tmp = zs(is1)
           zs(is1) = zs(is2)
           zs(is2) = tmp
           tmp = index(is1)
           index(is1) = index(is2)
           index(is2) = tmp
        endif
     enddo
  enddo
  rewind(in2)
  do is = 1,ns
     write(in2,'(i,a,e,a,e,a,e,$)') index(is),' ',xs(is),' ',ys(is),' ',zs(is)
     write(in2,*)
  enddo
  close(in2)

  nzss(:,:) = 0 ! ajout
  nzrr(:,:) = 0 ! ajout

  !++++++++++++
  !       on calcule :
  !       ncs: nombre de couches contenant un source
  !       isc(): liste des couches contenant un source
  !       nzs(i): nbre de sources de prof. differente dans couche i
  !       nzss(j,i): nbre de sources a la prof j, dans la couche i
  !       izss(,j,i): indice dans xr(),yr(),zr() des sources a la prof j
  !                   dans la couche i
  !++++++++++++

  do is = 1,ns
     ! compute ic,zc
     ic = 1    
     hh = hc(1)
     do while ((zs(is) .gt. hh) .and. (ic .lt. nc))
        zs(is) = zs(is)-hh
        ic = ic+1
        hh = hc(ic)
     enddo
     cff(is) = one/rho(ic)
     ! compute isc(),ncs,js
     if (is .eq. 1) then
        isc(1) = ic
        ncs = 1
        js = 1
     else
        is1 = 1
        tc = .true.
        do while (is1 .le. ncs)
           if (ic .eq. isc(is1)) then
              js = is1
              tc = .false.
           endif
           is1 = is1+1
        enddo
        if (tc) then
           ncs = ncs+1
           isc(ncs) = ic
           js = ncs
           nzs(js) = 0
        endif
     endif
     ! compute nzs(),jss
     if (is .eq. 1) then
        nzs(1) = 1
        jss = 1
        tc = .false.
     else
        is2 = 1
        tc = .true.
        do while (is2 .le. nzs(js))
           if (zs(is) .eq. zs(izss(1,is2,js))) then
              jss = is2
              tc = .false.
           endif
           is2 = is2+1
        enddo
     endif
     if (tc) then
        nzs(js) = nzs(js)+1
        jss = nzs(js)
     endif
     ! compute nzss(,),izss(,,)
     nzss(jss,js) = nzss(jss,js)+1
     izss(nzss(jss,js),jss,js) = is
  enddo


  !++++++++++++
  !	 on reordonne les stations par profondeur croissante
  !++++++++++++

  do ir1 = 1,nr-1
     do ir2 = ir1,nr
        if (zr(ir1) .gt. zr(ir2)) then
           tmp = xr(ir1)
           xr(ir1) = xr(ir2)
           xr(ir2) = tmp
           tmp = yr(ir1)
           yr(ir1) = yr(ir2)
           yr(ir2) = tmp
           tmp = zr(ir1)
           zr(ir1) = zr(ir2)
           zr(ir2) = tmp
        end if
     end do
  end do

  rewind(in3)
  do ir = 1,nr
     write(in3,'(e,a,e,a,e,$)') xr(ir),' ',yr(ir),' ',zr(ir)
     write(in3,*)
  end do
  close(in3)

  !++++++++++++
  !	on calcule :
  !	ncr: nombre de couches contenant un recepteur
  !	irc(): liste des couches contenant un recept
  !       nzr(i): nbre de recept. de prof. differente dans couche i
  !	nzrr(j,i): nbre de recept a la prof j, dans la couche i
  !	izrr(,j,i): indice dans xr(),yr(),zr() des recept a la prof j
  !		    dans la couche i
  !++++++++++++

  do ir = 1,nr
     ! compute ic,zc
     ic = 1	
     hh = hc(1)
     do while ((zr(ir) .gt. hh) .and. (ic .lt. nc))
	zr(ir) = zr(ir)-hh
	ic = ic+1
	hh = hc(ic)
     enddo
     ! compute irc(),ncr,jr
     if (ir .eq. 1) then 
	irc(1) = ic
	ncr = 1
	jr = 1
     else
	ir1 = 1
	tc = .true.
	do while (ir1 .le. ncr)
           if (ic.eq.irc(ir1)) then
              jr = ir1
              tc = .false.
           endif
           ir1 = ir1+1
        enddo
	if (tc) then
           ncr = ncr+1
           irc(ncr) = ic
           jr = ncr
           nzr(jr) = 0
        endif
     endif
     ! compute nzr(),jrr
     if (ir .eq. 1) then
	nzr(1) = 1
	jrr = 1
	tc = .false.
     else
	ir2 = 1
	tc = .true.
	do while (ir2 .le. nzr(jr))
           if (zr(ir) .eq. zr(izrr(1,ir2,jr))) then
              jrr = ir2
              tc = .false.
           endif
           ir2 = ir2+1
        enddo
     endif
     if (tc) then
	nzr(jr) = nzr(jr)+1
	jrr = nzr(jr)
     endif
     ! compute nzrr(,),izrr(,,)
     nzrr(jrr,jr) = nzrr(jrr,jr)+1
     izrr(nzrr(jrr,jr),jrr,jr) = ir
  enddo

  !++++++++++++
  !         distances radiales / source
  !         on ne garde que les distances differentes, stockees dans 
  !         rr(). tableau d'indirection irr().
  !++++++++++++

  nrs = 0		! calcule dist. rad.
  do is = 1,ns
     do ir = 1,nr
        nrs = nrs+1
        r(ir,is) = sqrt((xr(ir)-xs(is))*(xr(ir)-xs(is))+(yr(ir)-ys(is))*(yr(ir)-ys(is)))
        rr(nrs) = r(ir,is)
     enddo
  enddo

  ir1 = 1		! elimine dist. rad. egales
  do while (ir1 .lt. nrs)
     ir2 = ir1+1
     do while (ir2 .le. nrs)
	if (rr(ir1) .eq. rr(ir2)) then
           rr(ir2) = rr(nrs)
           nrs = nrs-1
        else
           ir2 = ir2+1
        endif
     enddo
     ir1 = ir1+1
  enddo

  ! tableau d'indirection
  do is = 1,ns
     do ir = 1,nr
        do ir2 = 1,nrs
           if (r(ir,is) .eq. rr(ir2)) irs(ir,is) = ir2
        enddo
     enddo
  enddo

  ! coef azimut.
  do is = 1,ns
     do ir = 1,nr
        if (r(ir,is) .ne. zero) then
           cosr(ir,is) = (xr(ir)-xs(is))/r(ir,is)
           sinr(ir,is) = (yr(ir)-ys(is))/r(ir,is)
        else
           cosr(ir,is) = one
           sinr(ir,is) = zero
        endif
     enddo
  enddo

  do ic = 1,nc
     vs2(ic) = vs(ic)*vs(ic)
     vp2(ic) = vp(ic)*vp(ic)
  end do

end subroutine initdata

