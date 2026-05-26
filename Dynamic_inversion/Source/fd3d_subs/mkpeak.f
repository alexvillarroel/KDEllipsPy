      subroutine mkpeak (nxt,nyt,nzt,nell,xo,yo,a,b,phi,strin,
     $    strout,ioutput,r)
c      real peak_xz(120,120)
      real xo(2),yo(2),a(2),b(2),phi(2)
      real strin strbout    

      include 'parstat'
c      print *,'in mpeak nxt, nyt = ',nxt, nyt 
c      real a(10),b(10),phi(20),xo(10),yo(10)
c
c  read size of grid
c

      do i=1,nxt
         do j=1,nyt
            peak_xz(i,j)=0.
         enddo
      enddo
c
      strcen=0.5
      if(ioutput.eq.1)open(79,file='peakin.dat')
      do i = 1,nxt
            do j = 1,nyt
                  inside=0
                  do ie=1,nell
                      x=(i-xo(ie))*cos(phi(ie))+(j-yo(ie))*sin(phi(ie))
                      y=-(i-xo(ie))*sin(phi(ie))+(j-yo(ie))*cos(phi(ie))
                      rr=(x/a(ie))**2+(y/b(ie))**2
                      rr=sqrt(rr)
                      if(rr.lt.1.)then
                           inside=inside+1
                      endif
                  enddo
                  if(inside>0)then
                      peak_xz(i,j) = peak_xz(i,j)+strin
                  else
                      peak_xz(i,j) = peak_xz(i,j)+strin
                  endif
                  if(ioutput.eq.1)write(79,*)peak_xz(i,j)
           enddo
      enddo

c      print *,'in mkpeak nell ,xo, yo =  ',nell, xo(1),yo(1)
c      print *,'a, b, phi, = ',a(1),b(1),phi(1)
c      print *,'strin, strout = ', strin,strout
c      print *, 'peak_xz =', peak_xz(50,50), peak_xz(20,20)     
 
      if(ioutput.eq.1)close(79)
      return
      end


