#!/bin/bash

read num_bose < ../matsubara_params

for k_boson in $(seq 1 ${num_bose})
do
	cat <<EOF> diag_${k_boson}.f90
program diagonalize_Sk
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none

    ! Variables
    integer :: phi_x, phi_y, phi_z, phi_size, i
    integer :: trim_x_start, trim_x_end
    integer :: trim_y_start, trim_y_end
    integer :: trim_z_start, trim_z_end
    integer :: info, lwork
    character(len=1) :: jobz, uplo
    real(wp), allocatable :: lap(:,:), W(:,:), Sk(:,:), eigvals(:), work(:)
    real(wp), allocatable :: wkopt_array(:)
    real(wp) :: wkopt

    ! Read trim ranges
    open(unit=20, file="../trim_xyz", status="old", action="read")
    read(20, *) trim_x_start, trim_x_end
    read(20, *) trim_y_start, trim_y_end
    read(20, *) trim_z_start, trim_z_end
    close(20)

    phi_x = trim_x_end - trim_x_start + 1
    phi_y = trim_y_end - trim_y_start + 1
    phi_z = trim_z_end - trim_z_start + 1
    phi_size = phi_x * phi_y * phi_z

    ! Allocate matrices
    allocate(lap(phi_size, phi_size))
    allocate(W(phi_size, phi_size))
    allocate(Sk(phi_size, phi_size))
    allocate(eigvals(phi_size))

    ! Read input matrices
    open(10, file="../lap.bin", form="unformatted", access="stream", status="old")
    read(10) lap
    close(10)

    open(11, file="../W_${k_boson}.bin", form="unformatted", access="stream", status="old")
    read(11) W
    close(11)

    Sk = -lap - W

    ! Set parameters for dsyev
    jobz = 'V'
    uplo = 'U'

    ! Workspace query
    allocate(wkopt_array(1))
    call dsyev(jobz, uplo, phi_size, Sk, phi_size, eigvals, wkopt_array, -1, info)
    if (info /= 0) then
        print*, "Error in DSYEV workspace query, info =", info
        stop
    end if

    wkopt = wkopt_array(1)
    lwork = int(wkopt)
    deallocate(wkopt_array)
    allocate(work(lwork))

    ! Actual diagonalization
    call dsyev(jobz, uplo, phi_size, Sk, phi_size, eigvals, work, lwork, info)
    if (info /= 0) then
        print*, "DSYEV failed, info =", info
        stop
    end if

    
    ! Write output
    open(20, file="eigvals_${k_boson}.dat", status="replace")
    do i = 1,phi_size
    write(20, *) eigvals(i)
    end do
    close(20)

    open(21, file="eigvecs_${k_boson}.bin", form="unformatted", access="stream", status="replace")
    write(21) Sk
    close(21)

    ! Cleanup
    deallocate(lap, W, Sk, eigvals, work)
end program diagonalize_Sk
EOF
done

