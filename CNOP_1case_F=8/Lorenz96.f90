
    subroutine model(m,nt,X0,X,F,h)
        implicit none
        integer m,nt
        real*8 X0(m),F,h
        real*8 X(0:nt,m),K1(m),K2(m),K3(m),K4(m)
        integer i
  
        X(0,:)=X0
  
        do i=0,nt-1
            call G(X(i,:),K1,m,F)
            call G(X(i,:)+h/2.d0*K1,K2,m,F)
            call G(X(i,:)+h/2.d0*K2,K3,m,F)
            call G(X(i,:)+h*K3,K4,m,F)
            X(i+1,:)=X(i,:)+h/6.d0*(K1+2.d0*K2+2.d0*K3+K4)
        enddo
  
        return
    end subroutine model
  
    subroutine G(X,K,m,F)
        implicit none
        integer m
        real*8 X(m),K(m),F
        integer j
  
        K(1)=(X(2)-X(m-1))*X(m)-X(1)+F
        K(2)=(X(3)-X(m))*X(1)-X(2)+F
        do j=3,m-1
          K(j)=(X(j+1)-X(j-2))*X(j-1)-X(j)+F
        enddo
        K(m)=(X(1)-X(m-2))*X(m-1)-X(m)+F
        return
    end subroutine G
  