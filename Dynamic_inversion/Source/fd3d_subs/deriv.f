
      subroutine dayz3d(nxb,nxe,nyb,nye,nzb,nze,nylim,nzlim,alp)

c     routine to perform exponential damping along yz intersect

c     nxb   starting point for damping in x dir    (integer)(sent)
c     nxe   ending point for damping in x dir      (integer)(sent)
c     nyb   starting point for damping in y dir    (integer)(sent)
c     nye   ending point for damping in y dir      (integer)(sent)
c     nzb   starting point for damping in z dir    (integer)(sent)
c     nze   ending point for damping in z dir      (integer)(sent)
c     nylim current limit in y-dir for damping     (integer)(sent)
c     nzlim current limit in y-dir for damping     (integer)(sent)
c     alp   parameter controlling damping          (real)   (sent)
 
      include 'parstat'
 
      do 60 k =  nzb,nze
      do 50 j =  nyb,nye
      ic = iabs(nzlim+1-k) + iabs(nylim+1-j)
      a2(j) = exp(-alp*ic)
   50 continue
      do 60 j =  nyb,nye
      do 60 i =  nxb,nxe
      u1(i,j,k) = u1(i,j,k)*a2(j)
      v1(i,j,k) = v1(i,j,k)*a2(j)
      w1(i,j,k) = w1(i,j,k)*a2(j)
   60 continue
      return
      end
c
      subroutine daxy3d(nxb,nxe,nyb,nye,nzb,nze,nxlim,nylim,alp)

c     routine to perform exponential damping along xy intersect

c     nxb   starting point for damping in x dir    (integer)(sent)
c     nxe   ending point for damping in x dir      (integer)(sent)
c     nyb   starting point for damping in y dir    (integer)(sent)
c     nye   ending point for damping in y dir      (integer)(sent)
c     nzb   starting point for damping in z dir    (integer)(sent)
c     nze   ending point for damping in z dir      (integer)(sent)
c     nylim current limit in y-dir for damping     (integer)(sent)
c     nxlim current limit in x-dir for damping     (integer)(sent)
c     alp   parameter controlling damping          (real)   (sent)

      include 'parstat'
 
      do 60 i =  nxb,nxe
      do 50 j =  nyb,nye
      ic = iabs(nxlim+1-i) + iabs(nylim+1-j)
      a2(j) = exp(-alp*ic)
   50 continue
      do 60 j =  nyb,nye
      do 60 k =  nzb,nze
      u1(i,j,k) = u1(i,j,k)*a2(j)
      v1(i,j,k) = v1(i,j,k)*a2(j)
      w1(i,j,k) = w1(i,j,k)*a2(j)
   60 continue
      return
      end
c
      subroutine daxz3d(nxb,nxe,nyb,nye,nzb,nze,nxlim,nzlim,alp)

c     routine to perform exponential damping along xz intersect

c     nxb   starting point for damping in x dir    (integer)(sent)
c     nxe   ending point for damping in x dir      (integer)(sent)
c     nyb   starting point for damping in y dir    (integer)(sent)
c     nye   ending point for damping in y dir      (integer)(sent)
c     nzb   starting point for damping in z dir    (integer)(sent)
c     nze   ending point for damping in z dir      (integer)(sent)
c     nzlim current limit in z-dir for damping     (integer)(sent)
c     nxlim current limit in x-dir for damping     (integer)(sent)
c     alp   parameter controlling damping          (real)   (sent)

      include 'parstat'
 
      do 60 k =  nzb,nze
      do 50 i =  nxb,nxe
      ic = iabs(nzlim+1-k) + iabs(nxlim+1-i)
      a1(i) = exp(-alp*ic)
   50 continue
      do 60 i =  nxb,nxe
      do 60 j =  nyb,nye
      u1(i,j,k) = u1(i,j,k)*a1(i)
      v1(i,j,k) = v1(i,j,k)*a1(i)
      w1(i,j,k) = w1(i,j,k)*a1(i)
   60 continue
      return
      end
c
      subroutine dxyz3d(nxb,nxe,nyb,nye,nzb,nze,nxlim,nylim,
     +nzlim,alp)

c     routine to perform exponential damping at xyz corner

c     nxb   starting point for damping in x dir    (integer)(sent)
c     nxe   ending point for damping in x dir      (integer)(sent)
c     nyb   starting point for damping in y dir    (integer)(sent)
c     nye   ending point for damping in y dir      (integer)(sent)
c     nzb   starting point for damping in z dir    (integer)(sent)
c     nze   ending point for damping in z dir      (integer)(sent)
c     nylim current limit in y-dir for damping     (integer)(sent)
c     nzlim current limit in z-dir for damping     (integer)(sent)
c     nxlim current limit in x-dir for damping     (integer)(sent)
c     alp   parameter controlling damping          (real)   (sent)

      include 'parstat'
 
      do 50 k =  nzb,nze
      do 50 j =  nyb,nye
      do 50 i =  nxb,nxe
      ic = iabs(nzlim+1-k) + iabs(nylim+1-j) + iabs(nxlim+1-i)
      a = exp(-alp*ic)
      u1(i,j,k) = u1(i,j,k)*a
      v1(i,j,k) = v1(i,j,k)*a
      w1(i,j,k) = w1(i,j,k)*a
   50 continue
      return
      end
c
      subroutine dazz3d(nxb,nxe,nyb,nye,nzb,nze,nzlim,alp)
      
c     routine to perform exponential damping along z = const

c     nxb   starting point for damping in x dir    (integer)(sent)
c     nxe   ending point for damping in x dir      (integer)(sent)
c     nyb   starting point for damping in y dir    (integer)(sent)
c     nye   ending point for damping in y dir      (integer)(sent)
c     nzb   starting point for damping in z dir    (integer)(sent)
c     nze   ending point for damping in z dir      (integer)(sent)
c     nzlim current limit in z-dir for damping     (integer)(sent)
c     alp   parameter controlling damping          (real)   (sent)

      include 'parstat'
 
      do 50 k =  nzb,nze
      ic = iabs(nzlim+1-k)
      a3(k) = exp(-alp*ic)
   50 continue

      do 60 k =  nzb,nze
      do 60 j =  nyb,nye
      do 60 i =  nxb,nxe
      u1(i,j,k) = u1(i,j,k)*a3(k)
      v1(i,j,k) = v1(i,j,k)*a3(k)
      w1(i,j,k) = w1(i,j,k)*a3(k)
   60 continue
      return
      end

      subroutine daxx3d(nxb,nxe,nyb,nye,nzb,nze,nxlim,alp)

c     routine to perform exponential damping along x = const

c     nxb   starting point for damping in x dir    (integer)(sent)
c     nxe   ending point for damping in x dir      (integer)(sent)
c     nyb   starting point for damping in y dir    (integer)(sent)
c     nye   ending point for damping in y dir      (integer)(sent)
c     nzb   starting point for damping in z dir    (integer)(sent)
c     nze   ending point for damping in z dir      (integer)(sent)
c     nxlim current limit in x-dir for damping     (integer)(sent)
c     alp   parameter controlling damping          (real)   (sent)

      include 'parstat'

      do 50 i =  nxb,nxe
      ic = iabs(nxlim+1-i)
      a1(i) = exp(-alp*ic)
   50 continue
      do 60 i =  nxb,nxe
      do 60 k =  nzb,nze
      do 60 j =  nyb,nye
      u1(i,j,k) = u1(i,j,k)*a1(i)
      v1(i,j,k) = v1(i,j,k)*a1(i)
      w1(i,j,k) = w1(i,j,k)*a1(i)
   60 continue
      return
      end
c
      subroutine dayy3d(nxb,nxe,nyb,nye,nzb,nze,nylim,alp)

c     routine to perform exponential damping along y = const

c     nxb   starting point for damping in x dir    (integer)(sent)
c     nxe   ending point for damping in x dir      (integer)(sent)
c     nyb   starting point for damping in y dir    (integer)(sent)
c     nye   ending point for damping in y dir      (integer)(sent)
c     nzb   starting point for damping in z dir    (integer)(sent)
c     nze   ending point for damping in z dir      (integer)(sent)
c     nylim current limit in y-dir for damping     (integer)(sent)
c     alp   parameter controlling damping          (real)   (sent)
     
      include 'parstat'

      do 50 j =  nyb,nye
      ic =  iabs(nylim+1-j)
      a2(j) = exp(-alp*ic)
   50 continue
      do 60 j =  nyb,nye
      do 60 k =  nzb,nze
      do 60 i =  nxb,nxe
      u1(i,j,k) = u1(i,j,k)*a2(j)
      v1(i,j,k) = v1(i,j,k)*a2(j)
      w1(i,j,k) = w1(i,j,k)*a2(j)
   60 continue
      return
      end

      subroutine dahs3d(nxt,nyt,nzt,mno)

c     routine to perform exponential damping at boundary
c     using a free surface

c     nxt   nodal points in x dir               (integer)(sent)
c     nyt   nodal points in y dir               (integer)(sent)
c     nzt   nodal points in z dir               (integer)(sent)
c     nmo   # of nodal points  to start damping (integer)(sent)

      include 'parstat'
 
      alp = .03/float(mno) 
      nxb = nxt-mno
      nye = nyt-mno
      nze = nzt-mno

c     boundary x=nx      

      call daxx3d(nxb+1,nxt,mno+1,nye,mno+1,nze,nxb-1,alp)

c     boundary x=1

      call daxx3d(1,mno,mno+1,nye,mno+1,nze,mno,alp)
 
c     boundary y=ny

      call dayy3d(mno+1,nxb,nye+1,nyt,mno+1,nze,nye-1,alp)

c     boundary y=1

      call dayy3d(mno+1,nxb,1,mno,mno+1,nze,mno,alp)

c     boundary z=nz

      call dazz3d(mno+1,nxb,mno+1,nye,nze+1,nzt,nze-1,alp)

c     boundary z=1

      call dazz3d(mno+1,nxb,mno+1,nye,1,mno,mno,alp)
 
c     line edge (x=1,y=1)
 
      call daxy3d(1,mno,1,mno,mno+1,nze,mno,mno,alp)

c     line edge (x=nx,y=1)

      call daxy3d(nxb+1,nxt,1,mno,mno+1,nze,nxb+1,mno,alp)
 
c     line edge (x=1,y=ny)

      call daxy3d(1,mno,nye+1,nyt,mno+1,nze,mno,nye-1,alp)

c     line edge (x=nx,y=ny)

      call daxy3d(nxb+1,nxt,nye+1,nyt,mno+1,nze,nxb-1,nye-1,alp)
 
c     line edge (y=1,z=1)

      call dayz3d(mno+1,nxb,1,mno,1,mno,mno,mno,alp)

c     line edge (y=ny,z=1)

      call dayz3d(mno+1,nxb,nye+1,nyt,1,mno,nye-1,mno,alp)
 
c     line edge (y=1,z=nz)

      call dayz3d(mno+1,nxb,1,mno,nze+1,nzt,mno,nze-1,alp)

c     line edge (y=ny,z=nz)

      call dayz3d(mno+1,nxb,nye+1,nyt,nze+1,nzt,nye-1,nze-1,alp)
 
c     line edge (x=1,z=1)

      call daxz3d(1,mno,mno+1,nye,1,mno,mno,mno,alp)

c     line edge (x=nx,z=1)

      call daxz3d(nxb+1,nxt,mno+1,nye,1,mno,nxb-1,mno,alp)
 
c     line edge (x=1,z=nz)

      call daxz3d(1,mno,mno+1,nye,nze+1,nzt,mno,nze-1,alp)

c     line edge (x=nx,z=nz)

      call daxz3d(nxb+1,nxt,mno+1,nye,nze+1,nzt,nxb-1,nze-1,alp)
 
c     corner (1,1,1)

      call dxyz3d(1,mno,1,mno,1,mno,mno,mno,mno,alp)

c     corner (nx,1,1)

      call dxyz3d(nxb+1,nxt,1,mno,1,mno,nxb-1,mno,mno,alp)
 
c     corner (1,ny,1)

      call dxyz3d(1,mno,nye+1,nyt,1,mno,mno,nye-1,mno,alp)

c     corner (nx,ny,1)

      call dxyz3d(nxb+1,nxt,nye+1,nyt,1,mno,nxb-1,nye-1,mno,alp)
 
c     corner (1,1,nz)

      call dxyz3d(1,mno,1,mno,nze+1,nzt,mno,mno,nze-1,alp)

c     corner (nx,1,nz)

      call dxyz3d(nxb+1,nxt,1,mno,nze+1,nzt,nxb-1,mno,nze-1,alp)
 
c     corner (1,ny,nz)

      call dxyz3d(1,mno,nye+1,nyt,nze+1,nzt,mno,nye-1,nze-1,alp)

c     corner (nx,ny,nz)

      call dxyz3d(nxb+1,nxt,nye+1,nyt,nze+1,nzt,nxb-1,nye-1,nze-1,alp)

      return
      end
     
      subroutine dvel(nxt,nyt,nzt,dh,dt)

c     4th order finite-difference of velocity components

c     nxt   nodal points in x dir  (integer)(sent)
c     nyt   nodal points in y dir  (integer)(sent)
c     nzt   nodal points in z dir  (integer)(sent)
c     dh    spatial discretization (real)   (sent)
c     dt    temporal discretization(real)   (sent)

      include 'parstat'
c
      dth =dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
c     Find displacement fields at time t+1/2
c
      call uxx1(3,nxt-2,3,nyt-2,3,nzt-2,dh,dt)
      call vyy1(3,nxt-3,3,nyt-3,3,nzt-2,dh,dt)
      call wzz1(3,nxt-3,3,nyt-2,3,nzt-3,dh,dt)
      return
      end
c
      subroutine uxx0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = dt/dh
c
c     Find u-displacement fields at time t+1/2
c
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe
c
      d=d1(i,j,k)
      u1(i,j,k)=u1(i,j,k)+(dth/d)*(
     +(xx(i,j,k)-xx(i-1,j,k))+ 
     +(xy(i,j,k)-xy(i,j-1,k))+
     +(xz(i,j,k)-xz(i,j,k-1)))

   50 continue
c
      return
      end
c
      subroutine uxx1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
c     Find u-displacement fields at time t+1/2
c
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe
c
      d=d1(i,j,k)
      u1(i,j,k)=u1(i,j,k)+(dth/d)*(
     +     c1*(xx(i,j,k)-xx(i-1,j,k))+ 
     +     c2*(xx(i+1,j,k)-xx(i-2,j,k))+
c
     +     c1*(xy(i,j,k)-xy(i,j-1,k))+
     +     c2*(xy(i,j+1,k)-xy(i,j-2,k))+
c
     +     c1*(xz(i,j,k)-xz(i,j,k-1))+
     +     c2*(xz(i,j,k+1)-xz(i,j,k-2)))

   50 continue
c
      return
      end
c
      subroutine vyy0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2nd order finite-difference of v1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = dt/dh
c
c     Find v-displacement fields at time t+1/2
c
      do 50 k= nzb,nze
         do 50 j= nyb,nye
            do 50 i= nxb,nxe

               d=d1(i,j,k)
               v1(i,j,k)=v1(i,j,k)+(dth/d)* (
     +              (xy(i+1,j,k)-xy(i,j,k))+
c     
     +              (yy(i,j+1,k)-yy(i,j,k))+
c     
     +              (yz(i,j,k)-yz(i,j,k-1)))
c
   50 continue

      return
      end
c
      subroutine vyy1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of v1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
c     Find v-displacement fields at time t+1/2
c
      do 50 k= nzb,nze
         do 50 j= nyb,nye
            do 50 i= nxb,nxe

               d=d1(i,j,k)
               v1(i,j,k)=v1(i,j,k)+(dth/d)* (
     +              c1*(xy(i+1,j,k)-xy(i,j,k))+
     +              c2*(xy(i+2,j,k)-xy(i-1,j,k))+
c
     +              c1*(yy(i,j+1,k)-yy(i,j,k))+
     +              c2*(yy(i,j+2,k)-yy(i,j-1,k))+  
c     
     +              c1*(yz(i,j,k)-yz(i,j,k-1))+
     +              c2*(yz(i,j,k+1)-yz(i,j,k-2)))
c
   50 continue

      return
      end
c
      subroutine wzz0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2nd order finite-difference of w1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = dt/dh
c
c     Find w-displacement fields at time t+1/2
c
      do 50 k= nzb,nze
         do 50 j= nyb,nye
            do 50 i= nxb,nxe

               d=d1(i,j,k)
               w1(i,j,k)=w1(i,j,k)+(dth/d)*(
     +              (xz(i+1,j,k)-xz(i,j,k))+
c     
     +              (yz(i,j,k)-yz(i,j-1,k))+
c     
     +              (zz(i,j,k+1)-zz(i,j,k)))
c     
   50 continue

      return
      end
c
      subroutine wzz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of w1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
c     Find w-displacement fields at time t+1/2
c
      do 50 k= nzb,nze
         do 50 j= nyb,nye
            do 50 i= nxb,nxe
               
               d=d1(i,j,k)
               w1(i,j,k)=w1(i,j,k)+(dth/d)*(
     +              c1*(xz(i+1,j,k)-xz(i,j,k))+
     +              c2*(xz(i+2,j,k)-xz(i-1,j,k))+
c
     +              c1*(yz(i,j,k)-yz(i,j-1,k))+
     +              c2*(yz(i,j+1,k)-yz(i,j-2,k))+
c
     +              c1*(zz(i,j,k+1)-zz(i,j,k))+
     +              c2*(zz(i,j,k+2)-zz(i,j,k-1)))
c     
   50 continue

      return
      end
c
      subroutine bnd2d(nxt,nyt,nzt,dt,dh)

c     bnd2d finds 2nd-order differencing of wave eq at bnd
c     to obtain the velocity values

c     nxt   nodal points in x dir          (integer)(sent)
c     nyt   nodal points in y dir          (integer)(sent)
c     nzt   nodal points in z dir          (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = dt/dh
      nzf = nzt -1
      nxf = nxt -1
      nyf = nyt -1
c
c     Find displacement fields at time t+1/2 at x=2 and x=nx-1 by 2nd
c     order differences
c
      call uxx0(2,2,2,nyt-1,2,nzt-1,dh,dt)
      call uxx0(nxt-1,nxt-1,2,nyt-1,2,nzt-1,dh,dt)
c
      call vyy0(2,2,2,nyt-2,2,nzt-1 ,dh,dt)
      call vyy0(nxt-2,nxt-2,2,nyt-2,2,nzt-1,dh,dt)
c
      call wzz0(2,2,2,nyt-1,2,nzt-2,dh,dt)
      call wzz0(nxt-2,nxt-2,2,nyt-1,2,nzt-2,dh,dt)
c
c     Find displacement fields at time t+1/2 at y=2 and y=ny-1
c
c
      call uxx0(3,nxt-2,2,2,2,nzt-1,dh,dt)
      call uxx0(3,nxt-2,nyt-1,nyt-1,2,nzt-1,dh,dt)
c
      call vyy0(3,nxt-3,2,2,2,nzt-1,dh,dt)
      call vyy0(3,nxt-3,nyt-2,nyt-2,2,nzt-1,dh,dt)
c
      call wzz0(3,nxt-3,2,2,2,nzt-2,dh,dt)
      call wzz0(3,nxt-3,nyt-1,nyt-1,2,nzt-2,dh,dt)
c
c     Find displacement fields at time t+1/2 at z=2 
c
      call uxx0(3,nxt-2,3,nyt-2,2,2,dh,dt)

      call vyy0(3,nxt-3,3,nyt-3,2,2,dh,dt)

      call wzz0(3,nxt-3,3,nyt-2,2,2,dh,dt)
c
      return
      end
c
c  ==============================================================
c  STRESS ADVANCE BY 4th Order FD 
c  ==============================================================
c
      subroutine dstres(nxt,nyt,nzt,dh,dt)

c     4th order finite-difference of stress components

c     nxt   nodal points in x dir          (integer)(sent)
c     nyt   nodal points in y dir          (integer)(sent)
c     nzt   nodal points in z dir          (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth =dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
c     Find displacement fields at time t+1/2
c
      nzf = nzt -2 
      nxf = nxt -2
      nyf = nyt -2
      do 60 k=3,nzf
         do 60 j=3,nyf
            do 60 i=2,nxf

               xl=lam1(i,j,k)
               xm=mu1(i,j,k)
               a = xl + 2.*xm
               b = xl
               
c
c     Find xx stress
c
               xx(i,j,k) = xx(i,j,k) + dth*a*               
     +              (c1*( u1(i+1,j,k) - u1(i,j,k)   )  +
     +              c2*( u1(i+2,j,k) - u1(i-1,j,k) )     )
c     
     +              +dth*b*( c1*( v1(i,j,k) - v1(i,j-1,k))+
     +              c2*( v1(i,j+1,k) - v1(i,j-2,k) )  +
c     
     +              c1*( w1(i,j,k)   - w1(i,j,k-1) )  +
     +              c2*( w1(i,j,k+1) - w1(i,j,k-2) ))
c     
c     Find yy stress
c
               yy(i,j,k) = yy(i,j,k) + dth*a*    
     +              (c1*( v1(i,j,k) - v1(i,j-1,k)   )  +
     +              c2*( v1(i,j+1,k) - v1(i,j-2,k) )     )
c     
     +              +dth*b*( c1*( u1(i+1,j,k)- u1(i,j,k) )+
     +              c2*( u1(i+2,j,k) - u1(i-1,j,k) )  +
c     
     +              c1*( w1(i,j,k)   - w1(i,j,k-1) )  +
     +              c2*( w1(i,j,k+1) - w1(i,j,k-2) ))
c     
c     Find zz stress
c
               zz(i,j,k) = zz(i,j,k) + dth*a*     
     +              (c1*( w1(i,j,k) - w1(i,j,k-1)   ) +
     +              c2*( w1(i,j,k+1) - w1(i,j,k-2) )  )
c     
     +              +dth*b*( c1*( u1(i+1,j,k)   - u1(i,j,k) )+
     +              c2*( u1(i+2,j,k) - u1(i-1,j,k) )  +
c     
     +              c1*( v1(i,j,k)   - v1(i,j-1,k) )  +
     +              c2*( v1(i,j+1,k) - v1(i,j-2,k) ))
 60         continue
c     
c     shear stresses
c

      dth=.5*dt/dh

      nzf = nzt -2
      nxf = nxt -2
      nyf = nyt -2
      do 70 k=2,nzf
         do 70 j=2,nyf
            do 70 i=3,nxf
c
c     Find xy stress
c
               xm1 = mu1(i,j,k)
               xm2 = mu1(i,j+1,k)
               xmu = xm1+xm2
               xy(i,j,k) = xy(i,j,k) + dth*xmu* 
     +              (c1*( u1(i,j+1,k)  - u1(i,j,k)   )  +
     +              c2*( u1(i,j+2,k)  - u1(i,j-1,k) )  +
c     
     +              c1*( v1(i,j,k)  - v1(i-1,j,k)   )  +
     +              c2*( v1(i+1,j,k)  - v1(i-2,j,k) ))
c
c     Find xz stress
c
               dth=.5*dt/dh
               
               xm1 = mu1(i,j,k)
               xm2 = mu1(i,j,k+1)
               xmu = xm1+xm2
               xz(i,j,k) = xz(i,j,k) + dth*xmu*
     +              (c1*( u1(i,j,k+1)  - u1(i,j,k)   )  +
     +              c2*( u1(i,j,k+2)  - u1(i,j,k-1) )  +
               
     +              c1*( w1(i,j,k)    - w1(i-1,j,k) )  +
     +              c2*( w1(i+1,j,k)  - w1(i-2,j,k) ))
c
c     Find yz stress
c
               dth=.5*dt/dh
               
               xm1 = mu1(i,j,k)
               xm2 = mu1(i+1,j+1,k+1)
               xmu = xm1+xm2
               yz(i,j,k) = yz(i,j,k) + dth*xmu* 
     +              (c1*( v1(i,j,k+1)  - v1(i,j,k)   )  +
     +              c2*( v1(i,j,k+2)  - v1(i,j,k-1) )  +
c     
     +              c1*( w1(i,j+1,k)    - w1(i,j,k) )  +
     +              c2*( w1(i,j+2,k)  - w1(i,j-1,k) ))
c     
   70 continue
c
      call sxy1(3,nxt-2,2,nyt-2,1,1,dh,dt)
      call sxy1(3,nxt-2,2,nyt-2,nzt-1,nzt,dh,dt)
c
      call sxz1(3,nxt-2,1,1,2,nzt-2,dh,dt)
      call sxz1(3,nxt-2,nyt-1,nyt,2,nzt-2,dh,dt)
c
      call syz1(1,2,2,nyt-2,2,nzt-2,dh,dt)
      call syz1(nxt-1,nxt-1,2,nyt-2,2,nzt-2,dh,dt)
c
      return
      end
c
      subroutine sxy0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2th order finite-difference of xy

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'

c
      dth = .5*dt/dh
c
      do 70 k=nzb,nze
         do 70 j=nyb,nye
            do 70 i=nxb,nxe
c
c     Find xy stress
c
               xm1 = mu1(i,j,k)
               xm2 = mu1(i,j+1,k)
               xmu = xm1+xm2
               
               xy(i,j,k) = xy(i,j,k) + dth*xmu*   
     +              (u1(i,j+1,k)-u1(i,j,k)+v1(i,j,k)-v1(i-1,j,k))
               
 70         continue

      return
      end
c
      subroutine sxy1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4th order finite-difference of xy

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = .5*dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
      do 70 k=nzb,nze
         do 70 j=nyb,nye
            do 70 i=nxb,nxe
c
c     Find xy stress
c
               xm1 = mu1(i,j,k)
               xm2 = mu1(i,j+1,k)
               xmu = xm1+xm2
               xy(i,j,k) = xy(i,j,k) + dth*xmu*   
     +              (c1*( u1(i,j+1,k)  - u1(i,j,k)   )  +
     +              c2*( u1(i,j+2,k)  - u1(i,j-1,k) )  +
c     
     +              c1*( v1(i,j,k)  - v1(i-1,j,k)   )  +
     +              c2*( v1(i+1,j,k)  - v1(i-2,j,k) ))
               
 70         continue

      return
      end
c
      subroutine syz0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2th order finite-difference of yz

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = .5*dt/dh
c
      do 70 k=nzb,nze
         do 70 j=nyb,nye
            do 70 i=nxb,nxe
c     
c     Find yz stress
c     
               xm1 = mu1(i,j,k)
               xm2 = mu1(i+1,j+1,k+1)
               xmu = xm1+xm2
               
               yz(i,j,k) = yz(i,j,k) + dth*xmu*  
     +              (v1(i,j,k+1)-v1(i,j,k)+w1(i,j+1,k)-w1(i,j,k))
c     
 70         continue

      return
      end
c
      subroutine syz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4th order finite-difference of yz

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = .5*dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
      do 70 k=nzb,nze
         do 70 j=nyb,nye
            do 70 i=nxb,nxe
c     
c     Find yz stress
c     
               xm1 = mu1(i,j,k)
               xm2 = mu1(i+1,j+1,k+1)
               xmu = xm1+xm2
               yz(i,j,k) = yz(i,j,k) + dth*xmu*  
     +              (c1*( v1(i,j,k+1)  - v1(i,j,k)   )  +
     +              c2*( v1(i,j,k+2)  - v1(i,j,k-1) )  +
               
     +              c1*( w1(i,j+1,k)    - w1(i,j,k) )  +
     +              c2*( w1(i,j+2,k)  - w1(i,j-1,k) ))
c     
   70 continue

      return
      end
c
      subroutine sxz0(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2th order finite-difference of xz

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = .5*dt/dh
c
      do 70 k=nzb,nze
         do 70 j=nyb,nye
            do 70 i=nxb,nxe
c     
c     Find xz stress
c     
               xm1 = mu1(i,j,k)
               xm2 = mu1(i,j,k+1)
               xmu = xm1+xm2
               xz(i,j,k) = xz(i,j,k) + dth*xmu*  
     +              (u1(i,j,k+1)-u1(i,j,k)+w1(i,j,k)-w1(i-1,j,k))
c     
 70         continue
            
      return
      end
c
      subroutine sxz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4th order finite-difference of xz

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'
c
      dth = .5*dt/dh
      c1 = 9./8.
      c2 = -1./24.
c
      do 70 k=nzb,nze
         do 70 j=nyb,nye
            do 70 i=nxb,nxe
c     
c     Find xz stress
c     
               xm1 = mu1(i,j,k)
               xm2 = mu1(i,j,k+1)
               xmu = xm1+xm2
               xz(i,j,k) = xz(i,j,k) + dth*xmu*  
     +              (c1*( u1(i,j,k+1)  - u1(i,j,k)   )  +
     +              c2*( u1(i,j,k+2)  - u1(i,j,k-1) )  +
c     
     +              c1*( w1(i,j,k)    - w1(i-1,j,k) )  +
     +              c2*( w1(i+1,j,k)  - w1(i-2,j,k) ))
c     
 70         continue
            
       return
       end
c
      subroutine strbnd(nxt,nyt,nzt,dh,dt)

c     2th order finite-difference of stresses at boundaries

c     nxt   nodal points in x dir  (integer)(sent)
c     nyt   nodal points in y dir  (integer)(sent)
c     nzt   nodal points in z dir  (integer)(sent)
c     dh    spatial discretization (real)   (sent)
c     dt    temporal discretization(real)   (sent)

      include 'parstat'
c
c     compute stresses at bnd with 2nd order difference operator
c
      call sxx(1,1,2,nyt,2,nzt,dh,dt)
      call sxx(nxt-1,nxt-1,2,nyt,2,nzt,dh,dt)

      call sxx(2,nxt-2,2,2,2,nzt,dh,dt)
      call sxx(2,nxt-2,nyt-1,nyt,2,nzt,dh,dt)

      call sxx(2,nxt-2,3,nyt-2,2,2,dh,dt)
      call sxx(2,nxt-2,3,nyt-2,nzt-1,nzt,dh,dt)
c
c     xy stress on x=2 plane and x=nxt-1
c
      call sxy0(2,2,1,nyt-1,1,nzt,dh,dt)
      call sxy0(nxt-1,nxt,1,nyt-1,1,nzt,dh,dt)
c
c     xy stress on y=2 plane and y=nyt-1
c
      call sxy0(3,nxt-2,1,1,1,nzt,dh,dt)
      call sxy0(3,nxt-2,nyt-1,nyt-1,1,nzt,dh,dt)
c
c     xz stress on z=1, x=2, x=nxt-1 planes
c
      call sxz0(2,nxt,1,nyt,1,1,dh,dt)
      call sxz0(2,nxt,1,nyt,nzt-1,nzt-1,dh,dt)
      call sxz0(2,2,1,nyt,2,nzt-2,dh,dt)
      call sxz0(nxt-1,nxt,1,nyt,2,nzt-2,dh,dt)
c
c     yz stress on z=1, y=1, y=nyt-1 planes
c
      call syz0(1,nxt-1,1,nyt-1,1,1,dh,dt)
      call syz0(1,nxt-1,1,nyt-1,nzt-1,nzt-1,dh,dt)
      call syz0(1,nxt-1,1,1,2,nzt-2,dh,dt)
      call syz0(1,nxt-1,nyt-1,nyt-1,2,nzt-2,dh,dt)
c
      return
      end
c
c  ==============================================================
c  STRESS FREE BOUNDARY CONDITION AT k=nzt
c  ==============================================================
c
      subroutine sxx(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2th order finite-difference of normal stresses 

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      include 'parstat'

c
      dth =dt/dh
c
c Find displacement fields at time t+1/2
c
      do 60 k=nzb,nze
         do 60 j=nyb,nye
            do 60 i=nxb,nxe

               xl=lam1(i,j,k)
               xm=mu1(i,j,k)
               a = xl + 2.*xm
               b = xl
c     
c     Find xx stress
c     
               xx(i,j,k) = xx(i,j,k) + 
     +              dth*a*( u1(i+1,j,k) - u1(i,j,k) ) +
     +              dth*b*( v1(i,j,k)  - v1(i,j-1,k)  +
     +              w1(i,j,k)   - w1(i,j,k-1) )  
c     
c     Find yy stress
c
               yy(i,j,k) = yy(i,j,k) + 
     +              dth*a*( v1(i,j,k) - v1(i,j-1,k)   )  +
     +              dth*b*(  u1(i+1,j,k)   - u1(i,j,k)   +
     +              w1(i,j,k)   - w1(i,j,k-1) )  
c     
c     Find zz stress
c
               zz(i,j,k) = zz(i,j,k) +
     +              dth*a*( w1(i,j,k) - w1(i,j,k-1) )  +
     +              dth*b*( u1(i+1,j,k) - u1(i,j,k)    +
     +              v1(i,j,k)   - v1(i,j-1,k) )  
c     
   60 continue
c
      return
      end
c
c ================================================================
c  ABSORBING B.C.
c ================================================================
c
      subroutine abc(nxt,nyt,nzt,dt,dh)
c
c     absorbing boundary condition
c
c     nxt   nodal points in x dir  (integer)(sent)
c     nyt   nodal points in y dir  (integer)(sent)
c     nzt   nodal points in z dir  (integer)(sent)
c     dh    spatial discretization (real)   (sent)
c     dt    temporal discretization(real)   (sent)

      include 'parstat'
c
      nxm  = nxt - 1
      nxm2 = nxt - 2
      nym  = nyt - 1
      nym2 = nyt - 2
      nzm  = nzt - 1
      nzm2 = nzt - 2
      dth=dt/dh
c
c    Apply abc at i=1 and i=nxt
c
c  composante u1
c
      do 50 k=1,nzt
         do 50 j=1,nyt
c
            d=d1(1,j,k)
            xl=lam1(1,j,k)
            xm=mu1(1,j,k)
            dp1=xl+2.*xm
            dp=dth*sqrt(dp1/d)
            u1(1,j,k)=u1(1,j,k)+dp*(u1(2,j,k)-u1(1,j,k))
            d=d1(nxt,j,k)
            xl=lam1(nxt,j,k)
            xm=mu1(nxt,j,k)
            dp1=xl+2.*xm
            dp=dth*sqrt(dp1/d)
            u1(nxt,j,k)=u1(nxt,j,k)-dp*(u1(nxt,j,k)-u1(nxm,j,k))
   50 continue
c
c  composante v1
c
      do 55 k=1,nzt
         do 55 j=1,nym
            d=d1(1,j,k)
            xm=mu1(1,j,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            v1(1,j,k)=v1(1,j,k)+ds*(v1(2,j,k)-v1(1,j,k))
            d=d1(nxm,j,k)
            xm=mu1(nxm,j,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            v1(nxm,j,k)=v1(nxm,j,k)-ds*(v1(nxm,j,k)-v1(nxm2,j,k))
   55 continue
c
c  composante w1
c
      do 57 k=1,nzt-1
         do 57 j=1,nyt
            d=d1(1,j,k)
            xm=mu1(1,j,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            w1(1,j,k)=w1(1,j,k)+ds*(w1(2,j,k)-w1(1,j,k))
            d=d1(nxm,j,k)
            xm=mu1(nxm,j,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            w1(nxm,j,k)=w1(nxm,j,k)-ds*(w1(nxm,j,k)-w1(nxm2,j,k))
   57 continue
c
c    Apply abcs at j=1 and j=nyt
c
c  composante u1
c
      do 60 k=1,nzt
         do 60 i=2,nxm
c
            d=d1(i,1,k)
            xm=mu1(i,1,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            u1(i,1,k)=u1(i,1,k)+ds*(u1(i,2,k)-u1(i,1,k))
            d=d1(i,nyt,k)
            xm=mu1(i,nyt,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            u1(i,nyt,k)=u1(i,nyt,k)-ds*(u1(i,nyt,k)-u1(i,nym,k))
   60 continue
c
c  composante v1
c
      do 65 k=1,nzt
         do 65 i=2,nxm2
            d=d1(i,1,k)
            xl=lam1(i,1,k)
            xm=mu1(i,1,k)
            dp1=xl+2.*xm
            dp=dth*sqrt(dp1/d)
            v1(i,1,k)=v1(i,1,k)+dp*(v1(i,2,k)-v1(i,1,k))
            d=d1(i,nym,k)
            xl=lam1(i,nym,k)
            xm=mu1(i,nym,k)
            dp1=xl+2.*xm
            dp=dth*sqrt(dp1/d)
            v1(i,nym,k)=v1(i,nym,k)-dp*(v1(i,nym,k)-v1(i,nym2,k))
 65      continue
c
c  composante w1
c
      do 66 k=1,nzt-1
         do 66 i=2,nxm2
            d=d1(i,1,k)
            xm=mu1(i,1,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            w1(i,1,k)=w1(i,1,k)+ds*(w1(i,2,k)-w1(i,1,k))
            d=d1(i,nyt,k)
            xm=mu1(i,nyt,k)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            w1(i,nyt,k)=w1(i,nyt,k)-ds*(w1(i,nyt,k)-w1(i,nym,k))
 66      continue
c
c    Apply abcs at k=0 and k=nzt
c
c  composante u1
c
      do 70 j=2,nym
         do 70 i=2,nxm
            d=d1(i,j,1)
            xm=mu1(i,j,1)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            u1(i,j,1)=u1(i,j,1)+ds*(u1(i,j,2)-u1(i,j,1))
            d=d1(i,j,nzt)
            xm=mu1(i,j,nzt)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            u1(i,j,nzt)=u1(i,j,nzt)-ds*(u1(i,j,nzt)-u1(i,j,nzm))
 70      continue
c
c composante v1
c     
      do 80 j=2,nym2
         do 80 i=2,nxm2
            d=d1(i,j,1)
            xm=mu1(i,j,1)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            v1(i,j,1)=v1(i,j,1)+ds*(v1(i,j,2)-v1(i,j,1))
            d=d1(i,j,nzt)
            xm=mu1(i,j,nzt)
            ds1=xm/d
            ds=dth*sqrt(ds1)
            v1(i,j,nzt)=v1(i,j,nzt)-ds*(v1(i,j,nzt)-v1(i,j,nzm))
 80   continue
c
c  composante w1
c
      do 90 j=2,nym
         do 90 i=2,nxm2
            d=d1(i,j,1)
            xl=lam1(i,j,1)
            xm=mu1(i,j,1)
            dp1=xl+2.*xm
            dp=dth*sqrt(dp1/d)
            w1(i,j,1)=w1(i,j,1)+dp*(w1(i,j,2)-w1(i,j,1))
            d=d1(i,j,nzt)
            xl=lam1(i,j,nzt)
            xm=mu1(i,j,nzt)
            dp1=xl+2.*xm
            dp=dth*sqrt(dp1/d)
            w1(i,j,nzm)=w1(i,j,nzm)-dp*(w1(i,j,nzm)-w1(i,j,nzm2))
   90 continue
      return
      end
c
c  ==============================================================
c  STRESS FREE BOUNDARY CONDITION AT k=nzt
c  Process stresses
c  ==============================================================
c
      subroutine fres(nxt,nyt,nzt,dh,dt)

c     free-surface B.C. for stresses

c     nxt   nodal points in x dir  (integer)(sent)
c     nyt   nodal points in y dir  (integer)(sent)
c     nzt   nodal points in z dir  (integer)(sent)
c     dh    spatial discretization (real)   (sent)
c     dt    temporal discretization(real)   (sent)

      include 'parstat'
c
      nyf = nyt - 2
      nxf = nxt - 2
      nzm = nzt-1

      nzt2=nzt+2
      nzt1=nzt+1
      nztm=nzt-1
      nztm2=nzt-2
c
      do 10 j=1,nyt
      do 10 i=1,nxt
c
      zz(i,j,nzt1) = -zz(i,j,nzt)
      zz(i,j,nzt2) = -zz(i,j,nztm)
c
      xz(i,j,nzt1) = -xz(i,j,nztm)
      xz(i,j,nzt2) = -xz(i,j,nztm2)
c
      yz(i,j,nzt1) = -yz(i,j,nztm)
      yz(i,j,nzt2) = -yz(i,j,nztm2)
c
c     zero yz and xz at free-surface.
c
      xz(i,j,nzt) = 0.
      yz(i,j,nzt) = 0.
c
   10 continue
      return
      end
c
c  ==============================================================
c  STRESS FREE BOUNDARY CONDITION AT k=nzt 
c  Process velocities
c  ==============================================================
c
      subroutine fuvw(nxt,nyt,nzt)

c     free-surface B.C. for velocities

c     nxt   nodal points in x dir (integer)(sent)
c     nyt   nodal points in y dir (integer)(sent)
c     nzt   nodal points in z dir (integer)(sent)

      include 'parstat'
c
      do 10 j=1,nyt
         do 10 i=2,nxt-1
            u1(i,j,nzt+1)=u1(i,j,nzt)-
     +           (w1(i,j,nzt)-w1(i-1,j,nzt))
   10 continue
c
      do 20 j=2,nyt-2
         do 20 i=2,nxt-1
            v1(i,j,nzt+1)=v1(i,j,nzt)-
     +           (w1(i,j+1,nzt)-w1(i,j,nzt))
   20 continue
c
      do 30 j=2,nyt-2
         do 30 i=2,nxt-2
            xl=lam1(i,j,nzt+1)
            xm=mu1(i,j,nzt+1)
            a=2.*xl
            b=xl+2.*xm
            w1(i,j,nzt+1)=w1(i,j,nzt)-(a/b)*(
     +           u1(i+1,j,nzt+1)-u1(i,j,nzt+1)+
     +           v1(i,j,nzt+1)-v1(i,j-1,nzt+1))
 30      continue
c
      return
      end




