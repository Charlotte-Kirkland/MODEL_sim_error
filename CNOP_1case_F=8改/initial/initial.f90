program main
    implicit none
    integer,parameter :: m=40

    real*8 :: x(m),y(m)
    !real*8 :: X(0:nt,m)
    integer:: i,j,k

    call random_seed()

    call random_number(x)
    call random_number(y)
    x=(x-y)*10
    print *,x

    open(11,file='initial3.txt')
    write(11,*) x

end program main
