c------------------------------------------------------------------------
c
c	Subroutine forward_modelling - calculates the array "predicted_
c				       data" from the model "rmodel"
c
c-------------------------------------------------------------------------

      subroutine forward_modelling(
     &              rmodel, nd,
     &              ndata, nwave,
     &              predicted_data)

      include 'kine_param.inc'

      character*40    chars, fname(maxwave)
      character       sourcefile*20,statfile*20, axihist*50, disp_name*40
      integer            nd,tt
      real*8            rmodel(max_nd)
      real        nlayer
      real*4      predicted_data(maxdata,maxwave)
      real*4      predicted_bis(maxdata,maxwave)

      integer     ndata, nwave,nsource

      real        pi,cte1,cte2
      real        a,b,alpha,np,tp,x01,y01,xe,ye,xpos,ypos,xx,yy
      real        gauss,d,smax,vr,rak,delay
      dimension    iis(1600),iid(1600)
      real*4      slip0(1600)
      real*4      estk,edip

      integer       ndd,nstk,ndip,nell,islip0,ishape

      namelist /input/ nc,nfreq,tl,aw,nr,ns,xl,t0,fr1,fr2,
     &                    ikmax,uconv,sourcefile,statfile
      integer     is(1600)
      real        phi_pi(1600), slip(1600), delta_pi(1600)
      real        rake(1600), ddip(1600), dstk(1600),ts(1600) 
      common /axihist/is, slip, phi_pi, delta_pi, rake,
     &        ddip, dstk, ts
      common /ellipse/ ndd,nstk,ndip,nell,islip0,estk,edip,slip0,ishape

      
      pi=3.14159265359
!       print *,'in f_m_ell nwave, nd, ndata = ',nwave,nd,ndata
c       logical

      in1=17
      in2=18
 
      nlayer=nc
      nsource=nd

c                         
                
 
c       les fichiers source et station sont les memes
c      seuls le fichier axi.hist change

      axihist='axi.hist'
c      open(in1,file=axihist,status='unknown')
      open(in1,file=axihist,status='unknown')

      do i=1,ndd
         read(in1,*) is(i), slip(i), phi_pi(i), delta_pi(i), rake(i),
     &        ddip(i), dstk(i), ts(i)
      enddo
      
      rewind(in1)
      
c********************************************************************
c     
c     PREMIERE ELLIPSE LIEE A L'HYPOCENTRE
c
      
      a=rmodel(1)*1000
      b=rmodel(2)*1000
      alpha=rmodel(3)*pi
      np=rmodel(4)
      tp=rmodel(5)*2*pi
      smax=rmodel(6)
      vr=rmodel(7)
c     rak=rmodel(8)

c      if (b > a) then
c         cte1=a
c         cte2=b
c         b=cte1
c         a=cte2
c      end if
      
c******
c     centre de l'ellipse
c     
      
      x01=a*np*cos(tp)
      y01=b*np*sin(tp)
      xe=x01*cos(alpha)+y01*sin(alpha)+estk
      ye=-x01*sin(alpha)+y01*cos(alpha)+edip
      
      i=1
      do idip=1,ndip
         do istk=1,nstk
            iis(i)=istk
            iid(i)=idip
            i=i+1
         enddo
      enddo
      
      do i=1,ndd
         xpos=dstk(i)*iis(i)
         ypos=ddip(i)*iid(i)
c     ts(i)=((xpos-dstk(i)/2-xe)**2+
c     &                    (ypos-ddip(i)-ye)**2)**0.5/vr/1000
         ts(i)=((xpos-dstk(i)/2-estk)**2+
     &        (ypos-ddip(i)-edip)**2)**0.5/vr/1000
         xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
         yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
         d=xx**2/a**2+yy**2/b**2
         if (ishape.eq.0) then
            gauss=smax
         elseif (ishape.eq.1) then
            gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
c            gauss=smax*exp(log(0.01)*(xx**2/a**2+yy**2/b**2))
         elseif (ishape.eq.2) then
            if (d.ge.1.) then
            gauss=0
            else
            gauss=smax*sqrt(1-d)
            endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
         endif
         
         if(d.le.1) then
            slip(i)=gauss+slip0(i)
         else
            slip(i)=0.+slip0(i)
         endif
      enddo
c********************************************************************
c     
c      AUTRES ELLIPSES
      
      if (nell.gt.1) then
         
         do ie=2,2
            a=rmodel(6)*1000
            b=rmodel(7)*1000
            alpha=rmodel(8)*pi
c            alpha=0.
            xe=rmodel(9)*1000
            ye=rmodel(10)*1000
            smax=rmodel(11)
            vr=rmodel(5)
            
            do i=1,ndd
               xpos=dstk(i)*iis(i)
               ypos=ddip(i)*iid(i)
               
               xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
               yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
               d=xx**2/a**2+yy**2/b**2
               if (ishape.eq.0) then
                  gauss=smax
               elseif (ishape.eq.1) then
                  gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
               elseif (ishape.eq.2) then
                 if (d.ge.1.) then
                 gauss=0
                 else
                gauss=smax*sqrt(1-d)
               endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
             endif
               
               if(d.le.1) then
                  ts(i)=((xpos-dstk(i)/2-estk)**2+
     &                 (ypos-ddip(i)-edip)**2)**0.5/vr/1000
                  slip(i)=gauss+slip(i)
               else
                  ts(i)=ts(i)
                  slip(i)=0.+slip(i)
               endif
            enddo
         enddo
         
      endif


c********************************************************************
c     
c      TERCERA AUTRES ELLIPSES
      
      if (nell.gt.2) then
         
         do ie=3,nell
            a=rmodel(12)*1000
            b=rmodel(12)*1000
c            alpha=rmodel(17)
            alpha=0.
            xe=rmodel(13)*1000
            ye=rmodel(14)*1000
            smax=rmodel(15)
            vr=rmodel(5)
            
            do i=1,ndd
               xpos=dstk(i)*iis(i)
               ypos=ddip(i)*iid(i)
               
               xx=(xpos-xe)*cos(alpha)-(ypos-ye)*sin(alpha)
               yy=(xpos-xe)*sin(alpha)+(ypos-ye)*cos(alpha)
               d=xx**2/a**2+yy**2/b**2
               if (ishape.eq.0) then
                  gauss=smax
               elseif (ishape.eq.1) then
                  gauss=smax*exp(-((xx**2/a**2)+(yy**2/b**2)))
               elseif (ishape.eq.2) then
            if (d.ge.1.) then
            gauss=0
            else
            gauss=smax*sqrt(1-d)
            endif
!            gauss=smax*(1-(xx**2/a**2+yy**2/b**2))
               endif
               
               if(d.le.1) then
                  ts(i)=((xpos-dstk(i)/2-estk)**2+
     &                 (ypos-ddip(i)-edip)**2)**0.5/vr/1000
                  slip(i)=gauss+slip(i)
               else
                  ts(i)=ts(i)
                  slip(i)=0.+slip(i)
               endif
            enddo
         enddo
         
      endif



c     
c     ECRITURE AXI.HIST SLIP
         
!      print *, 'escribe axi.hist'
      do i=1,ndd
         write(in1,831) is(i), slip(i), phi_pi(i), delta_pi(i),rake(i),
     &        ddip(i), dstk(i), ts(i) 
         
      enddo
      close(in1)
      
 831  format(i10,1x,f8.2,1x,f8.3,1x,f8.3,1x,f8.3,1x,e8.1,1x,
     &     e8.1,1x,f8.3)
      
!      print *,' avant sismo ndata, nwave =',ndata,nwave

      call seismograms(ndata,nwave,predicted_data,3)       
      
      return
      end
