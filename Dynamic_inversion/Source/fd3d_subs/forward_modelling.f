c
c------------------------------------------------------------------------
c
c      Subroutine forward_modelling - calculates the array "predicted_
c                               data" from the model "rmodel"
!      This routine controls all the computations from
!          numerical modelling (indyna3d)
!          to axitra (seismograms) 
c
c-------------------------------------------------------------------------
c
      subroutine forward_modelling(rmodel,moddim,ndata,ioutput,
     & iproc,amom)
c
        include 'fd3d_param.inc'
c
c
      real*8            rmodel(maxmoddim)
      real*4            amom
c
      common /pred_com/predicted_data(maxdata),iacc
      common /fd3d_com/nt,nxt,nyt,nz,dh,ndat,obs


      integer            ndata
c 
      real xo(2),yo(2),a(2),b(2),phi(2),strin,strbin
      real r(2),cte(2),dmax(2)
      real x0,y0,x1,y1,xa,ya,strbout,straout,strain
      
      real*4            dc(maxnx,maxny),te(maxnx,maxny),tu(maxnx,maxny)
      real*4            te2(maxnx,maxny)
c
      character*40      kname
c
c      real*4   radius,rr
c
      factor=400./dh 
      factor=1
c      print *,"factor = ",factor
c------ construct stress and peak models from the parameters
c       all numbers are in units of grid points!!!
c
c      strbin=1E+7
c      strbout=-8.0E+7
c
c data for initial asperity
c
c      strain=6.5E+7
c      straout=-8E+7
      xa=nxt/2
      ya=nyt/2
c      r=4*factor

      nell=1
      a(1)=rmodel(1)*factor
      b(1)=rmodel(2)*factor
      xo(1)=rmodel(3)
      yo(1)=rmodel(4)
      phi(1)=rmodel(5)   
      strbin=rmodel(6)*1E+6
      strin=rmodel(7)
      strin=strin*strbin
      strain=rmodel(8)
      strain=strain*strin
      strbout=-10*1E+7
      straout=0
      r(1)=rmodel(9)*factor
c      cte(1)=rmodel(9)
      dmax(1)=rmodel(10)
c       dmax(1)=0.5

c
c  size of grid
c
      x0=1
      x1=nxt
      y0=1
      y1=nyt
c      print *, 'mkstress xa,yb, x1,y1, r =',xa,ya,x1,y1,r
c      print *, 'mkstre strbin, strbout,strain, strin =',strbin,strbout,
c     *      strain,strin
      call mkstress(nxt,nyt,nzt,x0,y0,xo,yo,x1,y1,xa,ya,r,
     $        a,b,phi,strbin,strbout,strain,straout,nell,ioutput)
c       call mkstress(x0,y0,x1,y1,xa,ya,r,
c     $        strbin,strbout,strain,straout,ioutput)

c
c     data for peak
c
c       strin=5E+7
c       strout=strbout
c
c  les deux ellipses utilisent les memes valeurs de strin et strout
c      
c      print *,'mkpeak ellipse ', xo(1),yo(1),a(1),b(1),phi(1)
c      print *,'strin, strout ', strin,strout
      call mkpeak(nxt,nyt,nzt,nell,xo,yo,a,b,phi,strin,strout,ioutput,r) 
c        call mkpeak(nell,xo,yo,a,b,phi,strin,strout,ioutput) 
      
      predicted_data(1:ndata)=0
C
c  NT, DT are the values used in the fd3d program not those of the 
!  observed or modelled data (ndata, dt=10*dt)
!  en axitra esta dado por el valor de tl
!  y de numero = nfreq*4=256
C
c      print *,'va a indyna  iproc =',iproc
      call indyna3d(nt,dt,nskip,ioutput,iproc,dmax,amom)
c      print *,'entre indyna y seismograms nt, dt, nskip =',nt,dt,nskip
      call seismograms(nt,dt,nskip,ioutput,iproc)     
      return
      end
