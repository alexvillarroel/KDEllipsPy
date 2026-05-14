subroutine fft2cd(ntp,a,m,iwk)
  use spec_numbers
  implicit none
  include 'precision_subroutine.f90'
  
  integer :: ntp
  integer :: m
  integer,dimension(1:ntp) :: iwk
  complex(kind=pcal),dimension(1:ntp) :: a                                          
  
  
  integer :: i,isp,j,jj,jsp,k,k0,k1,k2,k3,kb,kn,mk,mm,mp,n,n4,n8,n2,lm,nn,jk
  real(kind=pcal) :: rad,c1,c2,c3,s1,s2,s3,ck,sk,sq,a0,a1,a2,a3,b0,b1,b2,b3,temp
  complex(kind=pcal) :: ak2
  logical :: done_it
                                                     
  sq = sqrt(two)/2
  sk = sin(pi/8)
  ck = cos(pi/8)
         
  mp = m+1                                                          
  n = 2**m                                                          
  iwk(1) = 1                                                        
  mm = (m/2)*2                                                      
  kn = n+1                                                          

  ! initialize work vector               

  do i=2,mp                                                      
     iwk(i) = iwk(i-1)+iwk(i-1)
  end do
  rad = twopi/n                                                     
  mk = m - 4                                                        
  kb = 1                                                            
  if (mm .ne. m) then                                           
     k2 = kn                                                           
     k0 = iwk(mm+1) + kb                                             
     done_it = .false.
     do while ((k0 .gt. kb) .or. (done_it .eqv. .false.))
        k2 = k2 - 1                                                       
        k0 = k0 - 1                                                       
        ak2 = a(k2)                                                       
        a(k2) = a(k0) - ak2                                               
        a(k0) = a(k0) + ak2                                               
        done_it = .true.
     end do
  end if
  c1 = one                                                          
  s1 = zero                                                         
  jj = 0                                                            
  k = mm - 1                                                        
  j = 4   
  
  if (k .ge. 1) go to 30                                            
  
  go to 70                                                          
  
20 continue
  
  do while (iwk(j) .le. jj)                                      
     jj = jj - iwk(j)                                                  
     j = j-1                                                           
     if (iwk(j) .gt. jj) go to 25                                      
     jj = jj - iwk(j)                                                  
     j = j - 1                                                         
     k = k + 2                                                         
     
  end do
  
25 continue
  
  jj = iwk(j) + jj                                                  
  j = 4                                                             
  
30 continue
  
  isp = iwk(k)           
  
  if (jj .ne. 0) then                                          
     ! reset trigonometric parameters
     c2 = jj * isp * rad                                               
     c1 = cos(c2)                                                      
     s1 = sin(c2)      
     
     c2 = c1 * c1 - s1 * s1                                            
     s2 = c1 * (s1 + s1)                                               
     c3 = c2 * c1 - s2 * s1                                            
     s3 = c2 * s1 + s2 * c1                                            
     
     goto 40
     
  end if
  
40 continue
  
  jsp = isp + kb                                                    
  ! determine fourier coefficients in groups of 4                     
  do i=1,isp                                                     
     k0 = jsp - i                                                   
     k1 = k0 + isp                                                  
     k2 = k1 + isp                                                  
     k3 = k2 + isp
     
     a0 = real(a(k0))                                               
     a1 = real(a(k1))                                                   
     a2 = real(a(k2))                                                    
     a3 = real(a(k3))
     
     b0 = aimag(a(k0))                                               
     b1 = aimag(a(k1))                                                   
     b2 = aimag(a(k2))                                                    
     b3 = aimag(a(k3))
     
     if (s1 .ne. zero) then                                     
        temp = a1                                                      
        a1 = a1 * c1 - b1 * s1                                         
        b1 = temp * s1 + b1 * c1                                       
        temp = a2                                                      
        a2 = a2 * c2 - b2 * s2                                         
        b2 = temp * s2 + b2 * c2                                       
        temp = a3                                                      
        a3 = a3 * c3 - b3 * s3                                         
        b3 = temp * s3 + b3 * c3 
     end if
     temp = a0 + a2                                                 
     a2 = a0 - a2                                                   
     a0 = temp                                                      
     temp = a1 + a3                                                 
     a3 = a1 - a3                                                   
     a1 = temp                                                      
     temp = b0 + b2                                                 
     b2 = b0 - b2                                                   
     b0 = temp                                                      
     temp = b1 + b3                                                 
     b3 = b1 - b3                                                   
     b1 = temp                                                      
     a(k0) = cmplx(a0+a1,b0+b1,pcal)                                     
     a(k1) = cmplx(a0-a1,b0-b1,pcal)                                     
     a(k2) = cmplx(a2-b3,b2+a3,pcal)                                     
     a(k3) = cmplx(a2+b3,b2-a3,pcal)                                     
  end do
  
  if (k .gt. 1) then                                            
     k = k - 2                                                         
     go to 30             
  end if

  kb = k3 + isp                                                     
  ! check for completion of final iteration                          
  if (kn .gt. kb) then                                          
     if (j .eq. 1) then                                           
        k = 3                                                             
        j = mk                                                            
        go to 20                
     end if
     j = j - 1                                                         
     c2 = c1                                                           
     if (j .eq. 2) then                                            
        c1 = c1 * ck + s1 * sk                                            
        s1 = s1 * ck - c2 * sk 
        c2 = c1 * c1 - s1 * s1                                            
        s2 = c1 * (s1 + s1)                                               
        c3 = c2 * c1 - s2 * s1                                            
        s3 = c2 * s1 + s2 * c1  
        go to 40                
     end if
     c1 = (c1 - s1) * sq                                               
     s1 = (c2 + s1) * sq     
     c2 = c1 * c1 - s1 * s1                                            
     s2 = c1 * (s1 + s1)                                               
     c3 = c2 * c1 - s2 * s1                                            
     s3 = c2 * s1 + s2 * c1     
     go to 40  
  end if
  
70 continue
                                                        
  ! permute the complex vector in reverse binary order to normal order                              
  if(m .gt. 1) then                                           
     mp = m+1                                                          
     jj = 1                                                            
     ! initialize work vector               
     iwk(1) = 1                                                        
     do i = 2,mp                                                   
        iwk(i) = iwk(i-1) * 2                                          
     end do
     n4 = iwk(mp-2)                                                    
     if (m .gt. 2) n8 = iwk(mp-3)                                      
     n2 = iwk(mp-1)                                                    
     lm = n2                                                           
     nn = iwk(mp)+1                                                    
     mp = mp-4                                                         
     ! determine indices and switch a       
     j = 2
     
     done_it = .false.
     
     do while ((j .le. lm) .or. (done_it .eqv. .false.))
        
        jk = jj + n2                                                      
        ak2 = a(j)                                                        
        a(j) = a(jk)                                                      
        a(jk) = ak2                                                       
        j = j+1                                                           
        if (jj .le. n4) then                                         
           jj = jj + n4                                                      
           go to 105 
        end if
        jj = jj - n4                                                      
        if (jj .le. n8) then                                        
           jj = jj + n8                                                      
           go to 105       
        end if
        jj = jj - n8                                                      
        k = mp
        do while (iwk(k) .lt. jj)                                 
           jj = jj - iwk(k)                                                  
           k = k - 1
        end do
        jj = iwk(k) + jj
        
105     continue
        
        if (jj .gt. j) then                                        
           k = nn - j                                                        
           jk = nn - jj                                                      
           ak2 = a(j)                                                        
           a(j) = a(jj)                                                      
           a(jj) = ak2                                                       
           ak2 = a(k)                                                        
           a(k) = a(jk)                                                      
           a(jk) = ak2
        end if
        j = j + 1                                                         
        ! cycle repeated until limiting number of changes is achieved             
        done_it = .true.
        
     end do
     
  end if
  
end subroutine fft2cd
