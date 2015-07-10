subroutine lennardjones(box,p,lx,ly,rmim,fxr,fyr,vr,R2)
      implicit none 
                integer, intent(in) :: box, lx, ly
                Real,intent(in) :: p(:,:)
                Real, intent(in) :: rmim
                real, intent(out) :: R2(box,box)
                real, intent(out) :: fxr(box),fyr(box),vr(box)
                real :: F(box,box),V(box,box),distx(box,box),disty(box,box),fx(box,box),fy(box,box)
                
                integer i,j
                                              
                do i=1,box 
                   distx(:,i) = p(:,1) - p(i,1) 
                   disty(:,i) = p(:,2) - p(i,2) 
                end do
                
                !Periodic condition of contour
                distx = distx -lx*nint(distx/lx)
                disty = disty -ly*nint(disty/ly) 
                
                R2 = distx*distx + disty*disty
                
                do i=1,box
                   do j=1,box
                      if(R2(i,j).gt.0.and.R2(i,j).lt.rmim) then
                         V(i,j) = 4*((1.0/(R2(i,j)**6))-(1/(R2(i,j)**3)))
                         F(i,j) = 48*((1.0/(R2(i,j)**13))-0.5*(1.0/(R2(i,j)**7)))
                         
                      else
                         V(i,j) = 0.0
                         F(i,j) = 0.0
                      end if
                   end do
                end do
      
                fx = F*distx
                fy = F*disty
                
                do i=1,(box-1)
                   fx(:,1) = fx(:,1) + fx(:,i+1)   
                   fy(:,1) = fy(:,1) + fy(:,i+1)   
                   V(:,1) = V(:,1) + V(:,i+1)   
                end do
                
                fxr = fx(:,1)
                fyr = fy(:,1)
                vr = V(:,1)
                
                return
end subroutine lennardjones


subroutine janusparticle(box,p,q,lx,ly,rmim,sig,a,C,C2,fxr,fyr,vr,R2,taur)
      implicit none 
                integer, intent(in) :: box, lx, ly
                Real,intent(in) :: p(:,:),q(:,:)
                Real, intent(in) :: rmim,sig,a,C,C2
                real, intent(out) :: R2(box,box)
                real, intent(out) :: fxr(box),fyr(box),vr(box),taur(box)
                real :: V(box,box),distx(box,box),disty(box,box),distqx(box,box)
                real :: distqy(box,box),fx(box,box),fy(box,box),tau(box,box),lua
                
                integer i,j
                
               
                               
                do i=1,box 
                   distx(:,i) = p(:,1) - p(i,1) 
                   disty(:,i) = p(:,2) - p(i,2) 
                end do

                do i=1,box 
                   distqx(:,i) = q(:,1) - q(i,1) 
                   distqy(:,i) = q(:,2) - q(i,2) 
                end do

                fx=0
                fy=0
                
                

                
                !Periodic condition of contour
                distx = distx -lx*nint(distx/lx)
                disty = disty -ly*nint(disty/ly) 
                
                R2 = distx*distx + disty*disty
                
                do i=1,box
                   do j=1,box
                      if(R2(i,j).gt.0.00001) then
                         lua = (distqx(i,j)*distx(i,j) + distqy(i,j)*disty(i,j))
                         
                         V(i,j) = lua*C*exp(-a*(sqrt(R2(i,j))-sig))/(R2(i,j)) + &
                              C2*4*((1.0/(R2(i,j)**6))-(1/(R2(i,j)**3)))

                         fx(i,j) = C*exp(-a*(sqrt(R2(i,j))-sig))*(lua &
                              *a*distx(i,j)/(R2(i,j))**(3./2) +2*lua* & 
                              distx(i,j)/(R2(i,j))**2 -distqx(i,j)/(R2(i,j))) &
                          + C2*48*((1.0/(R2(i,j)**(13./2.)))-0.5*(1.0/(R2(i,j)**(7./2.))))*distx(i,j)


                        
                         fy(i,j) = C*exp(-a*(sqrt(R2(i,j))-sig))*(lua &
                              *a*disty(i,j)/(R2(i,j))**(3./2) +2*lua* & 
                              disty(i,j)/(R2(i,j))**2 -distqy(i,j)/(R2(i,j))) &
                         + C2*48*((1.0/(R2(i,j)**(13./2.)))-0.5*(1.0/(R2(i,j)**(7./2.))))*disty(i,j)
                         


                         
                         
                      else
                         V(i,j) = 0.0
                         fx(i,j) = 0.0
                         fy(i,j) = 0.0
                      end if
                   end do
                end do

                
                tau = distx*fy - disty*fx                                


                do i=1,box
                   fxr(i)=SUM(fx(i,:))
                   fyr(i)=SUM(fy(i,:))
                   vr(i)=SUM(V(i,:))
                   taur(i)=SUM(tau(i,:))
                end do
                
               


                
                return
              end subroutine janusparticle



