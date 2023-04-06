!----the adjoint model ----------------------------------------
subroutine adj_model(m,nt,X0,F,h,ddX,ddXnt)
    implicit none
    integer m,nt
    real*8 X0(m),F,h,ddXnt(m)
    real*8 ddX(0:nt,m),matrix_J(m,m)
    real*8 J1(m,m),J2(m,m),J3(m,m),J4(m,m)
    real*8 X(0:nt,m),K1(m),K2(m),K3(m),K4(m)
    real*8 ddK1(m),ddK2(m),ddK3(m),ddK4(m)
    integer i,j
    
    X(0,:)=X0
    do i=0,nt-1
      call G(X(i,:),K1,m,F)
      call G(X(i,:)+h/2.d0*K1,K2,m,F)
      call G(X(i,:)+h/2.d0*K2,K3,m,F)
      call G(X(i,:)+h*K3,K4,m,F)
      X(i+1,:)=X(i,:)+h/6.d0*(K1+2.d0*K2+2.d0*K3+K4)
    enddo
    
    ddX(nt,:)=ddXnt
    
    do i=nt-1,0,-1
      call G(X(i,:),K1,m,F)
      call G(X(i,:)+h/2.d0*K1,K2,m,F)
      call G(X(i,:)+h/2.d0*K2,K3,m,F)
      call G(X(i,:)+h*K3,K4,m,F)
    
       call get_matrix(m,X(i,:),J1)
       call get_matrix(m,X(i,:)+h/2.d0*K1,J2)
       call get_matrix(m,X(i,:)+h/2.d0*K2,J3)
       call get_matrix(m,X(i,:)+h*K3,J4)
    
       call get_matrix_Jn(m,h,J1,J2,J3,J4,matrix_J)
    
      do j=1,m
       ddX(i,j)=dot_product(matrix_J(:,j),ddX(i+1,:))
      enddo
    
    enddo
    
    return
    end subroutine adj_model