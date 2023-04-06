
    subroutine model(nt,X0,X,h,F)
        implicit none
        integer nt
        real*8 X0(3),F,h
        real*8 X(0:nt,3),K1(3),K2(3),K3(3),K4(3)
        integer i

        X(0,:)=X0
  
        do i=0,nt-1
            call G(X(i,:),K1,F)
            call G(X(i,:)+h/2.d0*K1,K2,F)
            call G(X(i,:)+h/2.d0*K2,K3,F)
            call G(X(i,:)+h*K3,K4,F)
            X(i+1,:)=X(i,:)+h/6.d0*(K1+2.d0*K2+2.d0*K3+K4)
            !print *,i
        enddo
    
        return
    end subroutine model
  
    subroutine G(X,DX,F)
        implicit none
        !integer m
        real*8 X(3),DX(3),F,s,b,r
  
        s = 10.0
        b = 8.0/3
        r = 28.0

        DX(1)=s*(X(2)-X(1))+F
        DX(2)=r*X(1)-X(1)*X(3)-X(2)+F
        DX(3)=X(1)*X(2)-b*X(3)+F
        return
    end subroutine G
  