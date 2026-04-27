************************************************************************
* convertion geographique (lon,lat) en Lambert(x,y) pour les coordonnees
* des stations
************************************************************************

      program stations

      implicit none

      real xlat0,xlon0,Y,X,Z,xlat,xlon
      integer k,nb
      character*4 name
      
      open(10,file='station')
      open(11,file='stationn')
      write(6,*) 'latitude longitude de l''epicentre'
      read(5,*) xlat0,xlon0
      write(6,*) xlat0,xlon0
      call InitGlobalLambert(xlat0,xlon0)
      write(6,*) 'nombre de stations'
      read(5,*) nb
      write(6,*) nb
      write(6,*) 'latitude, longitude, altitude(m), nom(4 lettres)'
      write(6,*) 'x,y,lon,lat'
      do k=1,nb
          read(5,*) xlat,xlon,Z,name
          call GeoToGlobalLambert(xlat,xlon,Y,X)
          write(6,*) name,X,Y,xlon,xlat
          write(10,*) 1000*X,1000*Y,1000*Z 
          write(11,'(a4)') name 
          enddo
          close(10)
          close(11)
      end

      subroutine InitGlobalLambert(lat0,lon0)

c       PROJECTION LAMBERT CALCULEE EN SPHERIQUE
c       AVEC LAT0 ET LON0 COMME CENTRE DE PROJECTION }

      implicit none

      real lat0,lon0

      common /GlobalLambert/L0,latitude0,longitude0,n,pi,r0,initialized
      logical initialized
      real L0,latitude0,longitude0,n,pi,r0

        if (lon0.lt.0) lon0=lon0+360
      if (abs(lon0).ge.360.or.abs(lat0).ge.90) then
        print *,'InitGlobalLambert: parameter values out of range'
        stop
      endif
      pi=atan(1.)*4
      latitude0=lat0*pi/180
      longitude0=lon0*pi/180
      if (abs(latitude0).lt.1e-6) latitude0=1e-6
      if (abs(latitude0-pi/2).lt.1e-6) latitude0=pi/2-1e-6
      if (abs(latitude0+pi/2).lt.1e-6) latitude0=-pi/2+1e-6
      n=sin(latitude0)
      r0=6371/tan(latitude0)
      L0=log(tan(pi/4+latitude0/2))
      initialized=.true.
      end
!
!  subroutine GeoToGlobalLambert(la,lo,X,Y)
!
      subroutine GeoToGlobalLambert(la,lo,X,Y)
      
      implicit none

      real X,Y,la,lo

      common /GlobalLambert/L0,latitude0,longitude0,n,pi,r0,
     &         initialized
      logical initialized
      real L0,latitude0,longitude0,n,pi,r0
      
      real L,gamma,r,lat,lon

      lat=la
      lon=lo
      if (.not.initialized) then
        print *,'GeoToGlobalLambert: InitGlobalLambert missing'
        stop
      endif
      if (lon.lt.0 ) lon=lon+360
      if (abs(lon).ge.360.or.abs(lat).ge.90) then
        print *,'GeoToGlobalLambert: parameter values out of range'
        stop
      endif
      lon=lon*pi/180
      lat=lat*pi/180
      L=log(tan(pi/4+lat/2))
      gamma=n*(lon-longitude0)
      r=r0*exp(n*(L0-L))
      X=r*sin(gamma)
      Y=R0-r*cos(gamma)
      end
c      { -------------------- }
      subroutine GlobalLambertToGeo(X,Y,lat,lon)

      implicit none

      real X,Y,lat,lon
      
      common /GlobalLambert/L0,latitude0,longitude0,n,pi,r0,
     &                        initialized
      logical initialized
      real L0,latitude0,longitude0,n,pi,r0

      real L,gamma,r

      if (.not.initialized) then
        print *,'GlobalLambertToGeo: InitGlobalLambert missing'
        stop
      endif
      if (abs(r0-Y).lt.1e-6) then
        gamma=pi/2
        if (X.lt.0) gamma=-gamma
      else
        gamma=atan(X/(r0-Y))
        if (abs(gamma).lt.1e-6) then 
          r=r0-Y 
        else 
          r=X/sin(gamma)
        endif
      endif
      lat=atan(exp(L0-log(r/r0)/n))*360/pi-90
      lon=(longitude0+gamma/n)*180/pi
      if (lon.lt.-180) then 
        lon=lon+360 
      else if (lon.gt.180) then
        lon=lon-360
      endif
      end
c      { -------------------- }
