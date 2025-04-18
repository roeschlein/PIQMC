program phi_trim_fix
    implicit none

    ! Declare variables
    real*8 :: box_size_x, box_size_y, box_size_z
    integer :: N1, N2, N3, HO, KS_low, KS_high
    integer :: active_window_size
    integer :: x_start, x_stop, y_start, y_stop, z_start, z_stop
    integer :: x_size, y_size, z_size, i, j, k
    complex*16, allocatable :: phi_trim(:), phi(:)

    ! Open and read data from files
    open(unit=12, file="delta_xyz", status="old", action="read")
    read(12, *) box_size_x, box_size_y, box_size_z
    read(12, *) N1, N2, N3
    close(12)


    open(unit=13, file="active_window", status="old", action="read")
    read(13, *) HO, KS_low, KS_high
    close(13)

    active_window_size = KS_high - KS_low + 1

    open(unit=15, file="trim_xyz", status="old", action="read")
    read(15, *) x_start, x_stop
    read(15, *) y_start, y_stop
    read(15, *) z_start, z_stop
    close(15)

    ! Compute dimensions
    if (x_start > x_stop) x_stop = x_stop + N1
    if (y_start > y_stop) y_stop = y_stop + N2
    if (z_start > z_stop) z_stop = z_stop + N3

    x_size = x_stop - x_start + 1
    y_size = y_stop - y_start + 1
    z_size = z_stop - z_start + 1

    print*, "phi size: ", x_size * y_size * z_size
    print*, "x-dimensions: ", x_size
    print*, "y-dimensions: ", y_size
    print*, "z-dimensions: ", z_size

    ! Allocate memory
    allocate(phi_trim(x_size * y_size * z_size * active_window_size))
    allocate(phi(N1 * N2 * N3 * active_window_size))

    ! Read input phi data
    open(10, file="phi_temp.bin", form="unformatted", access="stream", status="old")
    read(10) phi
    close(10)
    print*, "Orbitals read"

    ! Call a subroutine for trimming phi
    call trim_phi(phi, phi_trim, x_start, x_stop, y_start, y_stop, z_start, z_stop, &
                  x_size, y_size, z_size, active_window_size, N1, N2, N3)

    ! Save output
    open(15, file="phi_trim.bin", form="unformatted", access="stream", status="replace")
    write(15) phi_trim
    close(15)
    print*, "phi_trim.bin written"

    open(17, file="trim_dimensions", status="unknown")
    write(17, *) x_size*y_size*z_size
    write(17, *) x_size
    write(17, *) y_size
    write(17, *) z_size
    close(17)


contains

    subroutine trim_phi(phi, phi_trim, x_start, x_stop, y_start, y_stop, z_start, z_stop, &
                    x_size, y_size, z_size, active_window_size, N1, N2, N3)
    implicit none
    integer, intent(in) :: x_start, x_stop, y_start, y_stop, z_start, z_stop
    integer, intent(in) :: x_size, y_size, z_size, active_window_size, N1, N2, N3
    complex*16, intent(inout) :: phi(:)
    complex*16, intent(out) :: phi_trim(:)
    integer :: i, j, k, l, counter
    integer :: i_mod, j_mod, k_mod

    counter = 1
    do l = 1, active_window_size
        do i = 0, x_size - 1
            do j = 0, y_size - 1
                do k = 0, z_size - 1
                    ! Apply periodic boundary conditions using modulo operations
                    i_mod = mod(x_start + i - 1, N1) + 1
                    j_mod = mod(y_start + j - 1, N2) + 1
                    k_mod = mod(z_start + k - 1, N3) + 1

                    ! Read from phi using wrapped indices
                    phi_trim(counter) = phi(k_mod + (j_mod-1)*N3 + (i_mod-1)*N2*N3 + (l-1)*N1*N2*N3)

                    counter = counter + 1
                end do
            end do
        end do
    end do
    end subroutine trim_phi


end program phi_trim_fix

