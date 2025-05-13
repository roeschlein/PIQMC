program laplacian_direct
    implicit none
    integer :: nx, ny, nz, trim_x_end, trim_x_start
    integer :: phi_size, trim_y_end, trim_y_start, trim_z_end, trim_z_start
    integer :: x, y, z, i, N1, N2, N3, num_fermi_freq, num_fermi_pairs
    real(8) :: box_size_x, box_size_y, box_size_z, T, pre
    integer :: xp1, xm1, yp1, ym1, zp1, zm1
    real(8), allocatable :: lap(:,:)
    real(8) :: ax, ay, az
    integer :: idx_center, idx_neighbor

    ! Read mat freq stuff
    open(unit=14, file="matsubara_params", status="old", action="read")
    read(14, *) num_fermi_freq
    read(14, *) num_fermi_pairs
    read(14, *) T
    close(14)

    ! Read the grid dimensions from delta_xyz
    open(unit=17, file="delta_xyz", status="old", action="read")
    read(17, *) box_size_x, box_size_y, box_size_z
    read(17, *) N1, N2, N3
    close(17)

    ! Read trim ranges
    open(unit=20, file="trim_xyz", status="old", action="read")
    read(20, *) trim_x_start, trim_x_end
    read(20, *) trim_y_start, trim_y_end
    read(20, *) trim_z_start, trim_z_end
    close(20)

    nx = trim_x_end - trim_x_start + 1
    ny = trim_y_end - trim_y_start + 1
    nz = trim_z_end - trim_z_start + 1

    ax = box_size_x/N1
    ay = box_size_y/N2
    az = box_size_z/N3

    pre = (ax * ay * az * T) / (8 * 3.141592653 * 51.4220674763)

    phi_size = nx * ny * nz
    allocate(lap(phi_size, phi_size), source=0.0d0)

    i = 1
    do x = 1, nx
        do y = 1, ny
            do z = 1, nz

                ! Compute wrapped neighbor indices
                xp1 = modulo(x, nx) + 1
                xm1 = modulo(x-2, nx) + 1
                yp1 = modulo(y, ny) + 1
                ym1 = modulo(y-2, ny) + 1
                zp1 = modulo(z, nz) + 1
                zm1 = modulo(z-2, nz) + 1

                ! Current center index
                idx_center = (x-1)*ny*nz + (y-1)*nz + z

                ! Neighbor in +x
                idx_neighbor = (xp1-1)*ny*nz + (y-1)*nz + z
                lap(i, idx_neighbor) = pre / ax**2

                ! Neighbor in -x
                idx_neighbor = (xm1-1)*ny*nz + (y-1)*nz + z
                lap(i, idx_neighbor) = pre / ax**2

                ! Neighbor in +y
                idx_neighbor = (x-1)*ny*nz + (yp1-1)*nz + z
                lap(i, idx_neighbor) = pre / ay**2

                ! Neighbor in -y
                idx_neighbor = (x-1)*ny*nz + (ym1-1)*nz + z
                lap(i, idx_neighbor) = pre / ay**2

                ! Neighbor in +z
                idx_neighbor = (x-1)*ny*nz + (y-1)*nz + zp1
                lap(i, idx_neighbor) = pre / az**2

                ! Neighbor in -z
                idx_neighbor = (x-1)*ny*nz + (y-1)*nz + zm1
                lap(i, idx_neighbor) = pre / az**2

                ! Center value
                lap(i, idx_center) = -2.0d0*pre*(1.0d0/ax**2 + 1.0d0/ay**2 + 1.0d0/az**2)

                ! Move to next row
                i = i + 1

            end do
        end do
    end do
    
    open(11, file="lap.bin", form="unformatted", access="stream", status="replace")
    write(11) lap
    close(11)

    print *, "Laplacian matrix written to 'lap.bin'"

    deallocate(lap)
end program laplacian_direct

