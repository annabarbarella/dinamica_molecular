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
