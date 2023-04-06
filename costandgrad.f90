!-------------cost function--------------------------------
subroutine evalf(m,X0,dX0,func)
    use commonpara
    implicit none
    integer m
    real*8 X0(m),dX0(m),func
    real*8 X1(0:nt,m),X2(0:nt,m)
    
    call model(m,nt,X0,X1,F,h)
    call model(m,nt,X0+dX0,X2,F,h)
    
    call Eu_norm(m,X2(nt,:)-X1(nt,:),func)
    
    func = -0.5d0*func**2.d0
    
    return
end subroutine evalf

!------------gradient----------------------------------------
subroutine evalg(m,X0,dX0,gfunc)
    use commonpara
    implicit none
    integer m
    real*8 X0(m),dX0(m),gfunc(m)
    real*8 X1(0:nt,m),X2(0:nt,m)
    real*8 ddX(0:nt,m)
    
    call model(m,nt,X0,X1,F,h)
    call model(m,nt,X0+dX0,X2,F,h)
    
    call adj_model(m,nt,X0+dX0,F,h,ddX,X2(nt,:)-X1(nt,:))
    
    gfunc=-1.d0*ddX(0,:)
    
    return
end subroutine evalg


subroutine proj(n, x)
    USE commonpara
    IMPLICIT NONE
    integer n,i
    real*8  sn
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
    !--------------------------------------------------------------------