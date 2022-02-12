!program for simple harmonic oscillator using RK4
program s_h_m
implicit none
integer::i,itrn
real*8::KE,PE,TE,w,t,t1,v,v1,x,x1,h,k1,k2,s1,s2,l1,l2,p1,p2,f1,f2,init_x,init_v
print*,"type the initial position"
read*,init_x
print*,"type the initial velocity"
read*,init_v
!print*,"type the value of spring constant"
!read*,k
print*,"the value of step size"
read*,h
print*,"type of number of iterations"
read*,itrn
open(unit=1,file="time_velocity.dat")
open(unit=2,file="time_position.dat")
open(unit=3,file="position_velocity.dat")
open(unit=4,file="tote.dat")
v=init_v
x=init_x
t=3.1415926/2*w
!print*,"type the angular frequency"
!read*,w
w=1.0
do i=1,itrn
k1=h*f1(t,v,x);                          k2=h*f2(t,v,x)
s1=h*f1(t+h/2.0,v+(k1)/2.0,x+(k2)/2.0);    s2=h*f2(t+h/2.0,v+(k1)/2.0,x+(k2)/2.0)
l1=h*f1(t+(h)/2.0,v+(s1)/2.0,x+(s2)/2.0);    l2=h*f2(t+h/2.0,v+(s1)/2.0,x+(s2)/2.0)
p1=h*f1(t+h,v+l1,x+l2);                  p2=h*f2(t+h,v+l1,x+l2)
v1=v+(k1+2*s1+2*l1+p1)/6.0
v=v1
x1=x+(k2+2*s2+2*l2+p2)/6.0
x=x1
t1=t+h
t=t1
TE=0.5*(v1**2)+0.5*(x1**2)
KE=0.5*(v1**2)
PE = 0.5*x1**2
write(1,*) t,v1
write(2,*) t,x1
write(3,*) x1,v1
write(4,*) t,KE,PE,TE
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 close(1)
 close(2)
 close(3)
 close(4)
 end program s_h_m
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function f1(t,v,x)
implicit none
real*8,intent(in)::t,v,x


f1=-x
end function f1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function f2(t,v,x)
implicit none
real*8,intent(in)::t,v,x
f2=v
end function f2
