program cal_1
    implicit none
    real(8) :: dx, dt, rho1, rho2, c1, c2, kappa1, kappa2, pi=4*atan(1.0_8), s1(5), s2(5), vec(5), zp=25_8, zs=1400_8, f0
    integer :: ios, nt, nx, i, j, nb
    real(8), allocatable :: u(:, :), u_nxt(:), u_now(:), u_pre(:), force(:), zz(:)

    ! constants [km,s]
    c1 = 3.6_8
    c2 = 4.7_8
    rho1 = 2700000_8 ![kg/km]
    rho2 = 3400000_8 ![kg/km]
    kappa1 = rho1 * c1**2
    kappa2 = rho2 * c2**2

    dx = 1.0_8
    dt = 0.1_8
    nt = 11000
    nx = 2001
    nb = 1971

    allocate(u(nt+2, nx), u_nxt(nx), u_now(nx), u_pre(nx), force(nx), zz(nx))
    
    u(:,:) = 0.0_8

    ! stencil
    s1 = (/ -1.0_8, 2-2*(c1*dt/dx)**2, (c1*dt/dx)**2, (c1*dt/dx)**2, dt**2/rho1 /)
    s2 = (/ -1.0_8, 2-2*(c2*dt/dx)**2, (c2*dt/dx)**2, (c2*dt/dx)**2, dt**2/rho2 /)

    ! force
    f0 = - 1.0_8 ![N/m=kg/s^2]
    zz = [(i, i=1,nx)]
    force = f0 * (1-2*(pi*((zz-1)*dx - zs)/zp)**2) * exp(-(pi*((zz-1)*dx - zs)/zp)**2)

    ! putput force
    open(unit=11, iostat=ios, file='force.dat', action='write', &
     & form='formatted', status='replace')
     do i = 1,nx
         write(11,'(f16.8)') force(i)
     enddo
    close(11)

    ! calculate step by step
    do i = 1,nt
        u_pre = u(i,:)
        u_now = u(i+1,:)
        do j = 2,nb-1
            vec = (/ u_pre(j), u_now(j), u_now(j+1), u_now(j-1), force(j) /)
            u_nxt(j) = dot_product(s1, vec)
        enddo
        do j = nb+1,nx-1
            vec = (/ u_pre(j), u_now(j), u_now(j+1), u_now(j-1), force(j) /)
            u_nxt(j) = dot_product(s2, vec)
        enddo
        u_nxt(1) = u_nxt(2)
        u_nxt(nx) = u_nxt(nx-1)
        u_nxt(nb) = (kappa2 * u_nxt(nb+1) + kappa1 * u_nxt(nb-1)) / (kappa2 + kappa1)
        u(i+2,:) = u_nxt
    enddo

    open(unit=12, file='output_u.dat', status='replace', action='write')

    ! output
    do i = 2, nt+2
        write(12, '(F10.4, 1X, ' // trim(adjustl(itoa(nx))) // 'F16.8)') (i-2)*dt, (u(i,j)*10.0_8**(6), j=1,nx)
    end do

    close(12)

    stop
    contains
    function itoa(n)
        integer, intent(in) :: n
        character(len=10) :: itoa
        write(itoa, '(I0)') n
    end function itoa

end program cal_1