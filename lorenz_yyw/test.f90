!!!!!!!!!!!!!!LORENZ63!!!!!!!!!!!!!!
!        To average x purbs        !
!         Author: Yiwei Ye         !
!        Finished: 2023.3.6        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program LORENZ63
	implicit none
	integer, parameter :: purbnum = 80000
	integer, parameter :: nt = 3000
	real *8, parameter :: f = 0.d0
	real *8, parameter :: h = 0.01
	real *8, parameter :: delta_r = 0.000000001
	real *8 :: ave(0:nt),costf(0:nt)
	real *8 :: x0(3),purb(3),x(0:nt,3),x_nopurb(0:nt,3),x_purb(0:nt,3),init(3)
    real *8 :: x_init(0:5000,3)
	real *8 :: a(3),b(3),delta,norm_x0
	integer :: ipurb,i,j

!!!!!!!!!!!!!!INITIALIZATION!!!!!!!!!!!!
    open(33,file='initial63.txt')
    read(33,*) init
    
!!!!!!!!!!!!!!SPIN-UP!!!!!!!!!!!!!!!!!!!
    !call spinup(init,x0,h,4000,f)
    call model(5000,init,x_init,h,f)
    x0=x_init(5000,:)
    x(0,:)=x0

!!!!!!!!!!!!!!CONTROL_RUN!!!!!!!!!!!!!!
    !open(112,file='info.txt')
    call model(nt,x0,x,h,f)
    call Eu_norm(3,x0,norm_x0)
    delta=norm_x0*delta_r
    print*, "delta=", delta
    x_nopurb=x
    !write(112,*) 'x_nopurb='
    !write(112,*) x_nopurb

!!!!!!!!!!!!!!+PURB!!!!!!!!!!!!!!!!!!!!!
    open(12,file='error_output.txt')
    ave=0.0
    do ipurb = 1,purbnum
    	call normal_random(3, purb, 0.0D0, 0.2D0)
        call proj(3,purb, delta)
        x(0,:)=x0
        call model(nt,x0+purb,x,h,f)
        x_purb=x
        !write(112,*) 'i=',ipurb
        !write(112,*) x_purb
        costf=0.0
        do j=0,nt
        	do i=1,3
        		costf(j)=costf(j)+(x_purb(j,i)-x_nopurb(j,i))**2
        	enddo
        	ave(j)=ave(j)+costf(j)
        enddo
    enddo

    ave=ave/purbnum
    do j=0,nt
      write(12,*) ave(j)
    enddo
    close(12)
    !close(112)

end program LORENZ63
 

Subroutine normal_random(n, nrandom, mean, std)
    Implicit None
    Integer, Intent (In) :: n
    Real *8, Intent (In) :: mean, std
    Real *8, Intent (Out) :: nrandom(n)
    Integer :: i, j, k
    Real *8 :: u1, u2, pai
  
    pai = dasin(1.D0)*2.D0
  
    Do k = 1, n
      Call random_number(u1)
      Call random_number(u2)
      nrandom(k) = mean + std*dsqrt(-2.0D0*dlog(u1))*dcos(2.D0*pai*u2)
    End Do
  
    Return
End Subroutine normal_random
  
subroutine proj(n, x, delta)
    !USE commonpara
    IMPLICIT NONE
    integer n,i
    real*8  sn,delta
    real*8  x(n)

    call Eu_norm(n,x,sn)
    if (sn .gt. delta) then
      do i=1,n
        x(i)=(delta/sn)*x(i)
      enddo
    endif
    return
end subroutine proj


subroutine Eu_norm(m,X0,norm_X0)
    implicit none
    integer m
    real*8 X0(m),norm_X0
    integer i
    
    norm_X0=0.d0
    do i=1,m
       norm_X0=norm_X0+X0(i)**2.d0
    enddo
    norm_X0=dsqrt(norm_X0)
    
    return
end subroutine Eu_norm
