program cal_2
    implicit none
    integer :: nw, nx, nt, nb, i, j, info
    real(8) :: rho1, rho2, c1, c2, kappa1, kappa2, dw, dx, dt, T, w_i, time, x
    real(8) :: zp=25.0_8, zs=1400.0_8, f0=-1.0_8 ![N/km=kg/s^2]
    real(8), parameter :: pi=4*atan(1.0_8)
    complex(8) :: omega, weight
    real(8), allocatable :: mat_T(:,:), mat_H(:,:), u(:,:), A(:)
    complex(8), allocatable :: vec_c(:), u_f(:,:), g(:), omegas(:), ab(:,:), complex_phase(:)
    integer, allocatable :: ipiv(:)

    ! constants [km,s]
    c1 = 3.6_8
    c2 = 4.7_8
    rho1 = 2700000_8 ![kg/km]
    rho2 = 3400000_8 ![kg/km]
    kappa1 = rho1 * c1**2
    kappa2 = rho2 * c2**2

    w_i = 2.3_8 * 10.0_8**(-3)

    dx = 1.0_8
    nx = 2000
    nb = 1971
    dt = 0.1_8
    nt = 11000
    T = 2048
    dw = 2 * pi / T
    nw = 10240 ! =T/(2*dt)
    

    allocate(mat_T(nx+1, nx+1), mat_H(nx+1, nx+1), g(nx+1), vec_c(nx+1), u_f(nx+1, 2*nw), ipiv(nx+1), u(nx+1, nt+1), omegas(2*nw), ab(4,nx+1), complex_phase(2*nw), A(nx+1))

    ! matrix setting
    call TH_matrix_setting(nx, nb, rho1, rho2, kappa1, kappa2, mat_T, mat_H, dx)

    do i = 1, nx+1
        x = real(i-1, kind=8) * dx
        if(i == 1) then
            A(i) = (force(x, f0, zp, zs) / 3 + force(x+dx, f0, zp, zs) / 6) * dx
        elseif(i == nx+1) then
            A(i) = (force(x, f0, zp, zs) / 3 + force(x-dx, f0, zp, zs) / 6) * dx
        else
            A(i) = (2 * force(x, f0, zp, zs) / 3 + force(x-dx, f0, zp, zs) / 6 + force(x+dx, f0, zp, zs) / 6) * dx
        endif
    enddo


    write(*,*) "matrix set was done"

    
    omegas = [(cmplx(i, 0.0_8, kind=8), i=0,2*nw-1)] * dw - cmplx(0, w_i, kind=8) 
    
    ! solve w^2T - H = -g in each omega
    do i = 0, nw
        omega = omegas(i+1)
        call make_band_matrix(ab, mat_T, mat_H, omega, nx, 1, 1)
        g = cmplx(0.0_8, -A, kind=8) / omega
        vec_c = -g
        call zgbsv(nx+1, 1, 1, 1, ab, 4, ipiv, vec_c, nx+1, info)
        if(info /= 0) then
            write(*,*) "error at", i
        endif
        
        if(i == nw) then
            u_f(:, i+1) = cmplx(real(vec_c), 0, kind=8)
        elseif(i == 0) then
            u_f(:, i+1) = vec_c
        else
            u_f(:, i+1) = vec_c
            u_f(:, 2*nw-i+1) = conjg(vec_c)
        endif
        if(mod(i,1000) == 0) then
            write(*,*) "u_f:", i
        endif
    enddo

    omegas = omegas + cmplx(0, w_i, kind=8)

    ! fourier transform
    do j = 1, nt+1
        time = real(j-1, kind=8) * dt
        complex_phase = exp(omegas * cmplx(0.0_8, time, kind=8))
        weight = dw * exp(w_i * time) / (2.0_8 * pi) ! 
    
        u(:, j) = real(matmul(u_f, complex_phase), kind=8) * weight

        if(mod(j,1000) == 0) then
            write(*,*) "u:", j
        endif
    enddo

    ! output
    open(unit=12, file='output_u_fem.dat', status='replace', action='write')

    do i = 1, nt+1
        write(12, '(F10.4, 1X, ' // trim(adjustl(itoa(nx+1))) // 'F16.8)') (i-1)*dt, (u(j,i)*10.0_8**(6), j=1,nx+1)
    end do

    close(12)

    deallocate(mat_T, mat_H, vec_c, u_f, g, omegas, u, ipiv, ab, complex_phase, A)
    

    stop
    contains

    ! create matrix T, H
    subroutine TH_matrix_setting(nx, nb, rho1, rho2, kappa1, kappa2, mat_T, mat_H, dx)
        implicit none
        integer, intent(in) :: nx, nb
        real(8), intent(in) :: rho1, rho2, kappa1, kappa2, dx
        real(8), intent(out) :: mat_T(:, :), mat_H(:, :)
        real(8) :: rho, kappa
        integer :: i
    
        mat_T(:,:) = 0.0_8
        mat_H(:,:) = 0.0_8
    
        do i = 1, nx+1
            if (i == 1) then
                rho = rho1
                kappa = kappa1
                mat_T(i,i) = 5 * rho / 12
                mat_H(i,i) = kappa
                mat_T(i,i+1) = rho / 12
                mat_T(i+1,i) = rho / 12
                mat_H(i,i+1) = -kappa
                mat_H(i+1,i) = -kappa
    
            else if (i >= 2 .and. i <= nb - 1) then
                rho = rho1
                kappa = kappa1
                mat_T(i,i) = 10 * rho / 12
                mat_T(i,i+1) = rho / 12
                mat_T(i+1,i) = rho / 12
                mat_H(i,i+1) = -kappa
                mat_H(i+1,i) = -kappa
                mat_H(i,i) = 2 * kappa
    
            else if (i == nb) then
                mat_T(i,i) = 5 * (rho1 + rho2) / 12
                mat_H(i,i) = kappa1 + kappa2
    
            else if (i >= nb + 1 .and. i <= nx) then
                rho = rho2
                kappa = kappa2
                mat_T(i,i) = 10 * rho / 12
                mat_T(i,i-1) = rho / 12
                mat_T(i-1,i) = rho / 12
                mat_H(i,i) = 2 * kappa
                mat_H(i,i-1) = -kappa
                mat_H(i-1,i) = -kappa
    
            else if (i == nx + 1) then
                rho = rho2
                kappa = kappa2
                mat_T(i,i) = 5 * rho / 12
                mat_T(i,i-1) = rho / 12
                mat_T(i-1,i) = rho / 12
                mat_H(i,i) = kappa
                mat_H(i,i-1) = -kappa
                mat_H(i-1,i) = -kappa
            end if
        end do
    
        mat_T = mat_T * dx
        mat_H = mat_H / dx
    
    end subroutine TH_matrix_setting

    ! transform w^2T-H to band matrix for zgbsv
    subroutine make_band_matrix(A_band, T, H, omega, nx, kl, ku)
        implicit none
        integer, intent(in) :: nx, kl, ku
        complex(8), intent(in) :: omega
        real(8), intent(in) :: T(:,:), H(:,:)
        complex(8), intent(out) :: A_band(:,:)
        integer :: j, ldab
        ldab = 2*kl + ku + 1
    
        A_band(:,:) = (0.0_8, 0.0_8)
    
        do j = 1, nx+1
            if (j > 1) then
                A_band(2,j) = omega**2 * cmplx(T(j-1,j), 0.0_8, kind=8) - cmplx(H(j-1,j), 0.0_8, kind=8)
            endif
            A_band(3,j) = omega**2 * cmplx(T(j,j), 0.0_8, kind=8) - cmplx(H(j,j), 0.0_8, kind=8)
            if (j < nx+1) then
                A_band(4,j) = omega**2 * cmplx(T(j+1,j), 0.0_8, kind=8) - cmplx(H(j+1,j), 0.0_8, kind=8)
            endif
        end do
        
    end subroutine
    
    ! force
    real(8) function force(x, f0, zp, zs)
        implicit none
        real(8), intent(in) :: x, f0, zp, zs

        force = f0 * (1 - 2 * (pi * (x-zs) / zp)**2) * exp(-(pi * (x-zs) / zp)**2)
        return
    endfunction force

    function itoa(n)
        integer, intent(in) :: n
        character(len=10) :: itoa
        write(itoa, '(I0)') n
    end function itoa


endprogram cal_2