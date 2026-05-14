c******************************************************************************
c*                                                                            *
c*                     SUBROUTINE INITDATA                                    *
c*                                                                            *
c*    Initialisation de divers parametres                                     *
c*                                                                            *
c*    Input:
c*      hc,zr,zs,nc,nr
c*    Output:
c*           ncr,irc,nzr,irzz,nzrr,rr
c*           ncs,isc,nzs,iszz,nzss,rr,iss
c*    Modified:
c*         hc
c******************************************************************************

      subroutine initdata

c Global
      include "parameter.inc"
      include "dimension1.inc"
      include "dimension2.inc"
c Local
      integer      ir,ir1,ir2,ic,jr,jrr,js,jss,is,is1,is2,index(nsp)
      logical      tc
      real      hh,tmp,r(nrp,nsp)

c++++++++++++
c        Lecture coordonnees stations et recepteurs
c++++++++++++

      do is=1,ns
        read(in2,*) index(is),xs(is),ys(is),zs(is)
      enddo
      do ir=1,nr
        read(in3,*) xr(ir),yr(ir),zr(ir)
      enddo


c++++++++++++
c        conversion interface -> epaisseur des couches     
c++++++++++++

      if (hc(1).eq.0.) then
        do ic=1,nc-1
         hc(ic)=hc(ic+1)-hc(ic)
        enddo
      endif

c++++++++++++
c        on reordonne les sources par profondeur croissante
c++++++++++++
      do is1=1,ns-1
       do is2=is1,ns
        if (zs(is1).gt.zs(is2)) then
         tmp=xs(is1)
         xs(is1)=xs(is2)
         xs(is2)=tmp
         tmp=ys(is1)
         ys(is1)=ys(is2)
         ys(is2)=tmp
         tmp=zs(is1)
         zs(is1)=zs(is2)
         zs(is2)=tmp
         tmp=index(is1)
         index(is1)=index(is2)
         index(is2)=tmp
        endif
       enddo
      enddo
      rewind (in2)
      do is=1,ns
        write(in2,'(i5,3g15.5)') index(is),xs(is),ys(is),zs(is)
      enddo
      close(in2)

c++++++++++++
c       on calcule :
c       ncs: nombre de couches contenant un source
c       isc(): liste des couches contenant un source
c       nzs(i): nbre de sources de prof. differente dans couche i
c       nzss(j,i): nbre de sources a la prof j, dans la couche i
c       izss(,j,i): indice dans xr(),yr(),zr() des sources a la prof j
c                   dans la couche i
c++++++++++++
 
      do is=1,ns
c                       compute ic,zc
       ic=1    
       hh=hc(1)
       do while ((zs(is).gt.hh).and.(ic.lt.nc))
        zs(is)=zs(is)-hh
        ic=ic+1
        hh=hc(ic)
       enddo
       cff(is)=1./rho(ic)
c                       compute isc(),ncs,js
       if (is.eq.1) then
        isc(1)=ic
        ncs=1
        js=1
       else
        is1=1
        tc=.true.
        do while (is1.le.ncs)
         if (ic.eq.isc(is1)) then
          js=is1
          tc=.false.
         endif
         is1=is1+1
        enddo
        if (tc) then
         ncs=ncs+1
         isc(ncs)=ic
         js=ncs
         nzs(js)=0
        endif
       endif
c                       compute nzs(),jss
       if (is.eq.1) then
        nzs(1)=1
        jss=1
        tc=.false.
       else
        is2=1
        tc=.true.
        do while (is2.le.nzs(js))
         if (zs(is).eq.zs(izss(1,is2,js))) then
          jss=is2
          tc=.false.
         endif
         is2=is2+1
        enddo
       endif
       if (tc) then
        nzs(js)=nzs(js)+1
        jss=nzs(js)
       endif
c                       compute nzss(,),izss(,,)
       nzss(jss,js)=nzss(jss,js)+1
       izss(nzss(jss,js),jss,js)=is
      enddo
 

c++++++++++++
c       on reordonne les stations par profondeur croissante
c++++++++++++
      do ir1=1,nr-1
       do ir2=ir1,nr
        if (zr(ir1).gt.zr(ir2)) then
         tmp=xr(ir1)
         xr(ir1)=xr(ir2)
         xr(ir2)=tmp
         tmp=yr(ir1)
         yr(ir1)=yr(ir2)
         yr(ir2)=tmp
         tmp=zr(ir1)
         zr(ir1)=zr(ir2)
         zr(ir2)=tmp
        endif
       enddo
      enddo
 
      rewind(in3)
      do ir=1,nr
        write(in3,*) xr(ir),yr(ir),zr(ir)
      enddo
      close(in3)

c++++++++++++
c      on calcule :
c      ncr: nombre de couches contenant un recepteur
c      irc(): liste des couches contenant un recept
c       nzr(i): nbre de recept. de prof. differente dans couche i
c      nzrr(j,i): nbre de recept a la prof j, dans la couche i
c      izrr(,j,i): indice dans xr(),yr(),zr() des recept a la prof j
c                dans la couche i
c++++++++++++

      do ir=1,nr
c                   compute ic,zc
       ic=1      
       hh=hc(1)
       do while ((zr(ir).gt.hh).and.(ic.lt.nc))
      zr(ir)=zr(ir)-hh
      ic=ic+1
      hh=hc(ic)
       enddo
c                  compute irc(),ncr,jr
       if (ir.eq.1) then 
      irc(1)=ic
      ncr=1
      jr=1
       else
      ir1=1
      tc=.true.
      do while (ir1.le.ncr)
       if (ic.eq.irc(ir1)) then
        jr=ir1
        tc=.false.
       endif
       ir1=ir1+1
        enddo
      if (tc) then
       ncr=ncr+1
       irc(ncr)=ic
       jr=ncr
       nzr(jr)=0
        endif
       endif
c                  compute nzr(),jrr
       if (ir.eq.1) then
      nzr(1)=1
      jrr=1
      tc=.false.
       else
      ir2=1
      tc=.true.
      do while (ir2.le.nzr(jr))
       if (zr(ir).eq.zr(izrr(1,ir2,jr))) then
        jrr=ir2
        tc=.false.
       endif
       ir2=ir2+1
        enddo
       endif
       if (tc) then
      nzr(jr)=nzr(jr)+1
      jrr=nzr(jr)
       endif
c                  compute nzrr(,),izrr(,,)
       nzrr(jrr,jr)=nzrr(jrr,jr)+1
       izrr(nzrr(jrr,jr),jrr,jr)=ir
      enddo

c++++++++++++
c         distances radiales / source
c         on ne garde que les distances differentes, stockees dans 
c         rr(). tableau d-indirection irr().
c++++++++++++
      nrs=0            !calcule dist. rad.
      do is=1,ns
      do ir=1,nr
         nrs=nrs+1
         r(ir,is)=sqrt((xr(ir)-xs(is))*(xr(ir)-xs(is))+
     &                 (yr(ir)-ys(is))*(yr(ir)-ys(is)))
         rr(nrs)=r(ir,is)
      enddo
      enddo
 
      ir1=1            !elimine dist. rad. egales
      do while (ir1.lt.nrs)
        ir2=ir1+1
        do while (ir2.le.nrs)
      if (rr(ir1).eq.rr(ir2)) then
        rr(ir2)=rr(nrs)
          nrs=nrs-1
        else
          ir2=ir2+1
        endif
        enddo
      ir1=ir1+1
      enddo

c Tableau d-indirection
      do is=1,ns
       do ir=1,nr
        do ir2=1,nrs
          if (r(ir,is).eq.rr(ir2)) irs(ir,is)=ir2
        enddo
       enddo
      enddo

c coef azimut.
      do is=1,ns
       do ir=1,nr
        if (r(ir,is).ne.0.) then
         cosr(ir,is)=(xr(ir)-xs(is))/r(ir,is)
         sinr(ir,is)=(yr(ir)-ys(is))/r(ir,is)
        else
         cosr(ir,is)=1.
         sinr(ir,is)=0.
        endif
       enddo
      enddo

      do 2 ic=1,nc
      vs2(ic)=vs(ic)*vs(ic)
      vp2(ic)=vp(ic)*vp(ic)
 2    continue

      return
      end
