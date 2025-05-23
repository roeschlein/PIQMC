#!/bin/bash

read num_bose < matsubara_params

cat <<EOF > W_1.f90
program Wm_xy_optimized_dense
    implicit none
    integer :: N1, N2, N3, active_window_size, HO, KS_low, KS_high
    integer*8 :: phi_size, phi_x, phi_y, phi_z, idx1, idx2
    integer :: num_non_neg_bos_freq
    integer :: i, j, k, x, y, z, xx, yy, zz
    complex*16, dimension(:,:,:,:), allocatable :: phi
    complex*16, dimension(:,:,:), allocatable :: G
    complex*16 :: W_sum
    real*8 :: box_size_x, box_size_y, box_size_z, dV
    integer :: trim_x_start, trim_x_end, trim_y_start, trim_y_end, trim_z_start, trim_z_end
    real*8, dimension(:,:), allocatable :: W

    ! Read the grid dimensions from delta_xyz
    open(unit=17, file="delta_xyz", status="old", action="read")
    read(17, *) box_size_x, box_size_y, box_size_z
    read(17, *) N1, N2, N3
    close(17)

    ! Read the active window information
    open(unit=18, file="active_window", status="old", action="read")
    read(18, *) HO, KS_low, KS_high
    close(18)

    active_window_size = KS_high - KS_low + 1
    dV = (box_size_x/N1) * (box_size_y/N2) * (box_size_z/N3)

    ! Read trim ranges
    open(unit=20, file="trim_xyz", status="old", action="read")
    read(20, *) trim_x_start, trim_x_end
    read(20, *) trim_y_start, trim_y_end
    read(20, *) trim_z_start, trim_z_end
    close(20)

    phi_x = trim_x_end - trim_x_start + 1
    phi_y = trim_y_end - trim_y_start + 1
    phi_z = trim_z_end - trim_z_start + 1
    phi_size = phi_x * phi_y * phi_z

    allocate(phi(active_window_size, N1, N2, N3))
    open(10, file="phi.bin", form="unformatted", access="stream", status="old")
    read(10) phi
    close(10)
    print*, "Orbitals read"

    ! Read Matsubara parameters
    open(15, file="matsubara_params", status="old", action="read")
    read(15, *) num_non_neg_bos_freq
    close(15)

    allocate(G(num_non_neg_bos_freq, active_window_size, active_window_size))
    open(13, file="Gm_ij.bin", form="unformatted", access="stream", status="old")
    read(13) G
    close(13)
    print*, "Gm_ij read"

    k = 1  ! The first bosonic frequency needs its own script because of the factor of 2 (see notes)

    allocate(W(phi_size, phi_size))
    W = 0.0d0

    print*, "Filling full dense W matrix of size", phi_size, "x", phi_size

    do x = trim_x_start, trim_x_end
    do y = trim_y_start, trim_y_end
    do z = trim_z_start, trim_z_end
        idx1 = (x - trim_x_start) * phi_y * phi_z + (y - trim_y_start) * phi_z + z - trim_z_start + 1

        do xx = trim_x_start, trim_x_end
        do yy = trim_y_start, trim_y_end
        do zz = trim_z_start, trim_z_end
            idx2 = (xx - trim_x_start) * phi_y * phi_z + (yy - trim_y_start) * phi_z + zz - trim_z_start + 1

            if (idx2 < idx1) cycle  ! Exploit symmetry

            W_sum = 0
            do i = 1, active_window_size
                do j = 1, active_window_size
                    W_sum = W_sum + ((dV)**2) * G(k,i,j) * CONJG(phi(i,x,y,z)) * phi(j,x,y,z) * &
                            CONJG(phi(j,xx,yy,zz)) * phi(i,xx,yy,zz)
                end do
            end do

            W(idx1, idx2) = dreal(W_sum)
            W(idx2, idx1) = W(idx1, idx2)  ! Fill symmetric lower triangle
        end do
        end do
        end do
    end do
    end do
    end do

    print*, "Full dense Wm_xy matrix filled!"

    open(14, file="W_1.bin", form="unformatted", access="stream", status="replace")
    write(14) W
    close(14)

    print*, "Wm_xy written to W_1.bin"

end program Wm_xy_optimized_dense
EOF


for k_boson in $(seq 2 ${num_bose})
do
	cat <<EOF > W_${k_boson}.f90
program Wm_xy_optimized_dense
    implicit none
    integer :: N1, N2, N3, active_window_size, HO, KS_low, KS_high
    integer*8 :: phi_size, phi_x, phi_y, phi_z, idx1, idx2
    integer :: num_non_neg_bos_freq
    integer :: i, j, k, x, y, z, xx, yy, zz
    complex*16, dimension(:,:,:,:), allocatable :: phi
    complex*16, dimension(:,:,:), allocatable :: G
    complex*16 :: W_sum
    real*8 :: box_size_x, box_size_y, box_size_z, dV
    integer :: trim_x_start, trim_x_end, trim_y_start, trim_y_end, trim_z_start, trim_z_end
    real*8, dimension(:,:), allocatable :: W

    ! Read the grid dimensions from delta_xyz
    open(unit=17, file="delta_xyz", status="old", action="read")
    read(17, *) box_size_x, box_size_y, box_size_z
    read(17, *) N1, N2, N3
    close(17)

    ! Read the active window information
    open(unit=18, file="active_window", status="old", action="read")
    read(18, *) HO, KS_low, KS_high
    close(18)

    active_window_size = KS_high - KS_low + 1
    dV = (box_size_x/N1) * (box_size_y/N2) * (box_size_z/N3)

    ! Read trim ranges
    open(unit=20, file="trim_xyz", status="old", action="read")
    read(20, *) trim_x_start, trim_x_end
    read(20, *) trim_y_start, trim_y_end
    read(20, *) trim_z_start, trim_z_end
    close(20)

    phi_x = trim_x_end - trim_x_start + 1
    phi_y = trim_y_end - trim_y_start + 1
    phi_z = trim_z_end - trim_z_start + 1
    phi_size = phi_x * phi_y * phi_z

    allocate(phi(active_window_size, N1, N2, N3))
    open(10, file="phi.bin", form="unformatted", access="stream", status="old")
    read(10) phi
    close(10)
    print*, "Orbitals read"

    ! Read Matsubara parameters
    open(15, file="matsubara_params", status="old", action="read")
    read(15, *) num_non_neg_bos_freq
    close(15)

    allocate(G(num_non_neg_bos_freq, active_window_size, active_window_size))
    open(13, file="Gm_ij.bin", form="unformatted", access="stream", status="old")
    read(13) G
    close(13)
    print*, "Gm_ij read"

    k = ${k_boson}

    allocate(W(phi_size, phi_size))
    W = 0.0d0

    print*, "Filling full dense W matrix of size", phi_size, "x", phi_size

    do x = trim_x_start, trim_x_end
    do y = trim_y_start, trim_y_end
    do z = trim_z_start, trim_z_end
        idx1 = (x - trim_x_start) * phi_y * phi_z + (y - trim_y_start) * phi_z + z - trim_z_start + 1

        do xx = trim_x_start, trim_x_end
        do yy = trim_y_start, trim_y_end
        do zz = trim_z_start, trim_z_end
            idx2 = (xx - trim_x_start) * phi_y * phi_z + (yy - trim_y_start) * phi_z + zz - trim_z_start + 1

            if (idx2 < idx1) cycle  ! Exploit symmetry

            W_sum = 0
            do i = 1, active_window_size
                do j = 1, active_window_size
                    W_sum = W_sum + 2 * ((dV)**2) * G(k,i,j) * CONJG(phi(i,x,y,z)) * phi(j,x,y,z) * &
                            CONJG(phi(j,xx,yy,zz)) * phi(i,xx,yy,zz)
                end do
            end do

            W(idx1, idx2) = dreal(W_sum)
            W(idx2, idx1) = W(idx1, idx2)  ! Fill symmetric lower triangle
        end do
        end do
        end do
    end do
    end do
    end do

    print*, "Full dense Wm_xy matrix filled!"

    open(14, file="W_${k_boson}.bin", form="unformatted", access="stream", status="replace")
    write(14) W
    close(14)

    print*, "Wm_xy written to W_${k_boson}.bin"

end program Wm_xy_optimized_dense
EOF
done
