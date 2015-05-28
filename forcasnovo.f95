subroutine lennardJones(x,y,lx,ly,xsize,ysize,rmim,fxr,fyr,vr,R2)
implicit none 
integer, intent(in) :: lx,ly,xsize,ysize
Real, intent(in) :: rmim
Real, intent(in) :: x(:)
Real, intent(in) :: y(:)
real, intent(out) :: R2(xsize,ysize)
real, intent(out) :: fxr(xsize),fyr(ysize),vr(xsize)
integer :: i,j
real :: F(xsize,ysize),V(xsize,ysize),distx(xsize,ysize),disty(xsize,ysize),fx(xsize,ysize),fy(xsize,ysize)

write(*,*)shape(x),shape(y)

do i=1,xsize 
   distx(i,:) = x - x(i) 
   disty(i,:) = y - y(i)
end do

R2 = distx*distx + disty*disty

do i=1,xsize
   do j=1,ysize
      if(R2(i,j).ge.0.and.R2(i,j).lt.rmim) then
         V(i,j) = 4*((R2(i,j)**(-6))-(R2(i,j)**(-3)))
         F(i,j) = 48*((R2(i,j)**(-13))-0.5*(R2(i,j)**(-7)))
      else
         V(i,j) = 0.0
         F(i,j) = 0.0
      end if
   end do
end do

fx = F*distx
fy = F*disty

do i=1,(xsize-1)
   fx(1,:) = fx(1,:) + fx(i+1,:)   
   fy(1,:) = fy(1,:) + fy(i+1,:)   
   V(1,:) = V(1,:) + V(i+1,:)   
end do

fxr = fx(1,:)
fyr = fy(1,:)
vr = V(1,:)

return
end subroutine lennardJones
