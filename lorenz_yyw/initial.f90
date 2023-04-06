!!!!!!!!!!!!!!LORENZ_INITIAL!!!!!!!!!!!!!!
!      To produce an initial field       !
!           Author: Yiwei Ye             !
!          Finished: 2023.3.15           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program initial
	implicit none
    integer, parameter :: m=3
	real *8 :: init(m)
	real *8 :: a(m),b(m)

    call random_seed()
    call random_number(a)
    call random_number(b)
    init=(a-b)*50
    open(11,file='initial63.txt')
    write(11,*) init
    close(11)

end program initial
 
