*  SUBROUTINE TO DO LAMBERT CONFORMAL CONIC PROJECTION
*  From USGS Bulletin no.1532, 1982, p. 105
      subroutine projlcc(rlon,rlat,xloc,yloc,isens)

      common /laccproj/ rlon0,radius,p0,fcon,cone,pi4
      common /trig/ pi,radpdeg

*  CONSTANTS SET IN EARLIER ROUTINE
*           cone=alog(cos(radpdeg*stdlat1)/cos(radpdeg*stdlat2))/
*    +         alog(tan(pi4+radpdeg*stdlat2/2.)/
*    +         tan(pi4+radpdeg*stdlat1/2.))
*           fcon=cos(radpdeg*stdlat1)*
*    +         (tan(pi4+radpdeg*stdlat1/2.)**cone)/cone
*           p0=radius*fcon/(tan(pi4+radpdeg*rlat0/2.)**cone)

      if (isens.eq.1) then
         p=radius*fcon/(tan(pi4+radpdeg*rlat/2.)**cone)
         theta=cone*radpdeg*(rlon-rlon0)
         xloc=p*sin(theta)
         yloc=p0-p*cos(theta)
      else
         theta=atan2(xloc, (p0-yloc))
         p=sqrt(xloc**2 + (p0-yloc)**2)
         rlon = theta/cone/radpdeg+rlon0
         rlat = (atan ((radius*fcon/p)**(1./cone))-pi4)*2./radpdeg
      endif
      
      return
      end
      
*     SUBROUTINE TO INITIALIZE LAMBERT CONF CONIC PROJECTION CONSTANTS
      subroutine init_lcc(rlat0,lon0,stdlat1,stdlat2)
      real lon0
      
      common /laccproj/ rlon0,radius,p0,fcon,cone,pi4
      common /trig/ pi,radpdeg
      rlon0=lon0
      pi4=atan(1.)
      pi=pi4*4.
      radpdeg=pi/180.
      radius=6378.
      
      cone=alog(cos(radpdeg*stdlat1)/cos(radpdeg*stdlat2))/
     +     alog(tan(pi4+radpdeg*stdlat2/2.)/
     +     tan(pi4+radpdeg*stdlat1/2.))
      fcon=cos(radpdeg*stdlat1)*
     +     (tan(pi4+radpdeg*stdlat1/2.)**cone)/cone
      p0=radius*fcon/(tan(pi4+radpdeg*rlat0/2.)**cone)
      
      return
      end
