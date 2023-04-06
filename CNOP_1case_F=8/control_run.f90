program main
    implicit none
    integer,parameter :: m=40
    real*8,parameter  :: F=8.d0
    real*8,parameter  :: h=0.05d0      ! 6h
    integer,parameter :: nt=4*365*10   ! ten years 
    real*8 :: X0(m)
    real*8 :: X(0:nt,m)
    integer:: i

    call read_initial(m,X0)

    X(0,:) = X0      

    call model(m,nt,X0,X,F,h)

    open(11,file='result/control_run.txt')

    do i=0,nt
        write(11,220) X(i,:)
    enddo

220 format (<m>(f24.12))

end program main
  
  subroutine read_initial(m,X0)
    implicit none
    integer :: m 
    real*8  :: X0(m)

    open(1,file='initial.txt')
    read(1,*) X0
  end subroutine read_initial