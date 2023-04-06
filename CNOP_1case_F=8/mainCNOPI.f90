module commonpara
  implicit none
  real*8,parameter  :: F=8.d0
  real*8,parameter  :: h=0.05d0      ! 6 hours
  integer,parameter :: nt=4*11        ! 7 days
  real*8,Parameter :: r_delta = 0.20
  REAL*8 :: delta      !the constraint
  
end module commonpara


Program maincnop
    use commonpara
    Implicit None
    Integer, Parameter :: m = 40
    Real *8 :: x0(m),norm_x0
    Real *8 :: x(0:nt,m)
    Real *8 :: mcnop(10,m),cnop(m)
    Real *8 :: mfunc(10)
    character(len=3):: ctemp1
    character(len=2):: ctemp2
    integer :: intemp
  
  !!parameter for spg iter!!
    Real*8, Parameter :: eps = 1.0D-10, eps2 = 1.0D-10
    Integer, Parameter :: maxit = 3000, maxfc = 3000
    Integer, Parameter :: mm = 10
    Logical, Parameter :: output = .FALSE.
    Real*8 :: pginfn, pgtwon, func
    Integer :: iter, fcnt, gcnt, flag, icnop, ptr, i
  
    Call readcase(m, x0)
    x(0,:) = x0
    Call model(m,nt,x0,x,F,h)
    call Eu_norm(m,x0,norm_x0)
    delta = norm_x0*r_delta
    print*, "delta=", delta
  
    Open (11, File='result/case.txt')
    Do i = 0, nt
      Write (11, 220) x(i,:)
    End Do
    Close (11)

  
    Print *, 'begin to search CNOP>>>>>>>>>>>>>>>'
  
    Do icnop = 1,10
      Call initialcnop(m,icnop,cnop)
      Call spg(m,cnop, mm, eps, eps2, maxit, maxfc, output,&
               func, pginfn, pgtwon, iter, fcnt, gcnt, flag,X0)
      mcnop(icnop, :) = cnop
      mfunc(icnop) = func
      Write (*, 5000) icnop, func
    End Do
  
    ptr = minloc(mfunc, 1)
    CNOP  = mcnop(ptr, :)
    func = mfunc(ptr)
  
    Print *, '********************'
    Print *, '****The best CNOP:'
    Write (*, 5000) ptr, func
    Print *, '********************'
    Print *, ''
  
    
    x(0, :) = x0
    Call model(m, nt, x0+cnop, x, f, h)



    intemp = r_delta*100
    write(ctemp1,'(i3.3)')  intemp
    write(ctemp2,'(i2.2)')  nt/4
    Open (12, File='result/CNOP_r'//ctemp1//'_d'//ctemp2//'.txt')
    Open (13, File='result/CNOPcase_r'//ctemp1//'_d'//ctemp2//'.txt')
    Write (12, 220) cnop
    Do i = 0, nt
      Write (13, 220) x(i,:)
    End Do
    Close (12)
    Close (13)
  
  
  
5000 Format (' Icnop = ', I3, ' obejct_function = ', F15.4)
220 Format (<m>(F24.12))
  
  End Program maincnop
  
  Subroutine readcase(m, x0)
    Implicit None
    Integer :: m, i, k
    Real *8 :: x0(m)
  
    k = 4*365 + 1
    Open (1, File='result/control_run.txt')
    Do i = 0, k
      Read (1, 230) x0
    End Do
  
    !Read (1, 230) x0
    Close (1)
230 Format (<m>(F24.12))
  End Subroutine readcase
  
  Subroutine initialcnop(m, icnop, cnop)
    Implicit None
    Integer :: m, icnop
    Real *8 ::  cnop(m), mcnop(10, m)
    Integer :: i
  
    Call normal_random(m, mcnop(1,:), 0.0D0, 0.2D0)
    Call normal_random(m, mcnop(2,:), 0.2D0, 0.3D0)
    Call normal_random(m, mcnop(3,:), -0.1D0, 0.4D0)
    mcnop(4, :) = 2.0D0
    mcnop(5, :20) = 2.0D0
    mcnop(5, 21:m) = -2.0D0
    mcnop(6:10, :) = -2.0D0*mcnop(1:5, :)
    cnop = mcnop(icnop, :)
    Return
  End Subroutine initialcnop
  
  
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
  