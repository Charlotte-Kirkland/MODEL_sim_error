subroutine tlm_model(m,nt,X0,F,h,dX,dX0)
    implicit none
    integer m,nt
    real*8 X0(m),F,h,ddXnt(m),dX0(m)
    real*8 dX(0:nt,m),matrix_J(m,m)
    real*8 J1(m,m),J2(m,m),J3(m,m),J4(m,m)
    real*8 X(0:nt,m),K1(m),K2(m),K3(m),K4(m)
    
    integer i,j
    
    X(0,:)=X0
    do i=0,nt-1
      call G(X(i,:),K1,m,F)
      call G(X(i,:)+h/2.d0*K1,K2,m,F)
      call G(X(i,:)+h/2.d0*K2,K3,m,F)
      call G(X(i,:)+h*K3,K4,m,F)
      X(i+1,:)=X(i,:)+h/6.d0*(K1+2.d0*K2+2.d0*K3+K4)
    enddo
    
    dX(0,:)=dX0
    do i=0,nt-1
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
       dX(i+1,j)=dot_product(matrix_J(j,:),dX(i,:))
      enddo
    
    enddo
    
    return
    end subroutine tlm_model

    subroutine get_matrix_Jn(m,h,J1,J2,J3,J4,matrix_J)
        implicit none
        integer m
        real*8 h
        real*8 J1(m,m),J2(m,m),J3(m,m),J4(m,m),matrix_J(m,m)
        real*8 J5(m,m),J6(m,m),J7(m,m),J_trans(m,m),J_I(m,m)
        integer i,j
        
        J_I=0.d0
        do i=1,m
          J_I(i,i)=1.d0
        enddo
        
        J_trans=J_I+(h/2.d0)*J1
        do i=1,m
          do j=1,m
            J5(i,j)=dot_product(J2(i,:),J_trans(:,j))
          enddo
        enddo
        
        J_trans=J_I+(h/2.d0)*J5
        do i=1,m
          do j=1,m
            J6(i,j)=dot_product(J3(i,:),J_trans(:,j))
          enddo
        enddo
        
        J_trans=J_I+h*J6
        do i=1,m
          do j=1,m
            J7(i,j)=dot_product(J4(i,:),J_trans(:,j))
          enddo
        enddo
        
        matrix_J=J_I+h/6.d0*(J1+2.d0*J5+2.d0*J6+J7)
        
        return
    end subroutine get_matrix_Jn
        
    subroutine get_matrix(m,X0,matrix_J)
        implicit none
        integer m
        real*8 X0(m),matrix_J(m,m)
        integer j
        
        matrix_J=0.d0
        matrix_J(1,1)=-1.d0
        matrix_J(1,2)=X0(m)
        matrix_J(1,m-1)=-X0(m)
        matrix_J(1,m)=X0(2)-X0(m-1)
        matrix_J(2,1)=X0(3)-X0(m)
        matrix_J(2,2)=-1.d0
        matrix_J(2,3)=X0(1)
        matrix_J(2,m)=-X0(1)
        do j=3,m-1
        matrix_J(j,j-2)=-X0(j-1)
        matrix_J(j,j-1)=X0(j+1)-X0(j-2)
        matrix_J(j,j)=-1.d0
        matrix_J(j,j+1)=X0(j-1)
        enddo
        matrix_J(m,1)=X0(m-1)
        matrix_J(m,m-2)=-X0(m-1)
        matrix_J(m,m-1)=X0(1)-X0(m-2)
        matrix_J(m,m)=-1.d0
        
        return
    end subroutine get_matrix        



