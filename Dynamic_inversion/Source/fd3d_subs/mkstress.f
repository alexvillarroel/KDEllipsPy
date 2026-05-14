       subroutine mkstress(nxt,nyt,nzt,x0,y0,xo,yo,x1,y1,xa,ya,r,
     $        a,b,phi,strbin,strbout,strain,straout,nell,ioutput)

       include 'parstat'



c      Defining variables:      
       REAL :: x0,y0,xo(2),yo(2),phi(2),a(2),b(2),x1,y1,xa,ya,r
       REAL :: strbout,strain,straout,strbin
       REAL :: inside,inside1
        
c      Creating the elliptical patch and nucleation circle 
c      within the grid of the discretized fault: 
       if(ioutput.eq.1)open(69,file='stressin.dat')
       do j=1,nyt
          do i=1,nxt
             strinix(i,j)=0.
             inside=0
             inside1=0 
             do ie=1,nell
                ! Creating elliptical patch:
                x=(i-xo(ie))*cos(phi(ie))+(j-yo(ie))*sin(phi(ie))
                y=-(i-xo(ie))*sin(phi(ie))+(j-yo(ie))*cos(phi(ie))
                rr=(x/a(ie))**2+(y/b(ie))**2
                rr=sqrt(rr)
                  
                ! Check if points are inside the elliptical patch:
                if(rr .lt. 1)then
                   inside=inside+1
                endif

                ! Restrict the hypocenter inside a smaller ellipse 
                ! with semi-axes of (a-r) and (b-r):
                xh=(xa-xo(ie))*cos(phi(ie))+(ya-yo(ie))*sin(phi(ie))
                yh=-(xa-xo(ie))*sin(phi(ie))+(ya-yo(ie))*cos(phi(ie))
                rhypo=(xh/(a(ie)-r))**2+(yh/(b(ie)-r))**2
                rhypo=sqrt(rhypo)
                if(rhypo .lt. 1)then
                   inside1=inside1+1
                endif
             enddo

             ! Stress values for points inside and outside 
             ! of the elliptical patch:
             if(inside .gt. 0)then
                strinix(i,j) = strinix(i,j)+strbin
             else
                strinix(i,j) = strinix(i,j)+strbout
             endif

             ! Stress values for points inside and outside 
             ! of the nucleation circle:
             if(inside1 .gt. 0)then
                dr2=(i-xa)**2+(j-ya)**2-r**2
                if(dr2 .le. 0)then
                   strinix(i,j)=strinix(i,j)+strain
                else
                   strinix(i,j)=strinix(i,j)+straout
                endif
             endif

             if(ioutput .eq. 1)write (69,*)strinix(i,j)
          enddo
       enddo

c       print *,'x0, y0, xo, yo  = ',x0, y0, xo, yo 
c       print *,'xa, ya, r  = ',xa, ya, r 
c       print *,'a, b, phi  = ',a, b, phi 
c       print *,'strain, straout, =',strain, straout 
c       print *,'strbin, strbout,=', strbin, strbout  
c       print *,'strinix, =', strinix(50,50), strinix(55,55)
c       print *,'strinix, =', strinix(40,40), strinix(60,60)
c       print *,'strinix, =', strinix(10,10), strinix(90,90)  


       if(ioutput .eq. 1)close(69)
       return
       end
