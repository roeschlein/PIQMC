#!/bin/bash


for A0_int in $(seq 1 100)
do
	cat <<EOF> A0_config_gen_${A0_int}.f90
program normal_random
    implicit none
    integer :: phi_x, phi_y, phi_z, phi_size, num_non_neg_bos_freq
    integer :: trim_x_start, trim_x_end, Nt, skip, ti, config
    integer :: trim_y_start, trim_y_end
    integer :: trim_z_start, trim_z_end
    real :: u1, u2
    real :: z0, z1
    real*8, allocatable :: eigval(:)
    real*8, allocatable :: non_neg_bos_freq(:)
    real*8, allocatable :: eigvec(:,:)
    complex*8, allocatable :: v(:)
    real*8 :: weight, delta, T
    complex*16, allocatable :: impt_smpl(:,:)
    complex*16, allocatable :: Aw(:,:)
    real*8, allocatable :: At(:,:)
    logical :: contains_nan
    integer :: seed
    integer :: i, j, n, w

    Nt = 120

    open(15, file="../matsubara_params", status="old", action="read")
    read(15, *) num_non_neg_bos_freq
    read(15, *) skip
    read(15, *) T
    close(15)

    delta = 1/(Nt*T)

    allocate(non_neg_bos_freq(num_non_neg_bos_freq))
    open(17, file='../non_neg_bosonic_frequencies', status="old", action="read")
    do i =1, num_non_neg_bos_freq
        read(17,*) non_neg_bos_freq(i)
    end do
    close(17)

    open(unit=20, file="../trim_xyz", status="old", action="read")
    read(20, *) trim_x_start, trim_x_end
    read(20, *) trim_y_start, trim_y_end
    read(20, *) trim_z_start, trim_z_end
    close(20)

    phi_x = trim_x_end - trim_x_start + 1
    phi_y = trim_y_end - trim_y_start + 1
    phi_z = trim_z_end - trim_z_start + 1
    phi_size = phi_x * phi_y * phi_z

    do config = 1+10*(${A0_int}-1), 10*${A0_int}
        print*, 'Generating configuration '//trim(str(config))//' '

        allocate(eigval(phi_size))
        allocate(impt_smpl(num_non_neg_bos_freq ,phi_size))
        n = phi_size

        call random_seed()
        do j = 1, num_non_neg_bos_freq
            open(10, file='eigvals_'//trim(str(j))//'.dat', status="old", action="read")
            do i =1, phi_size
                read(10,*) eigval(i)
            end do
            close(10)

            do i = 1, phi_size
                call random_number(u1)
                call random_number(u2)
                if (u1 == 0.0) u1 = 1.0e-10

                z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * acos(-1.0) * u2)
                z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * acos(-1.0) * u2)

                weight = 1/(2*sqrt(abs(eigval(i))))
                impt_smpl(j ,i) = weight * z0 + complex(0,1) * weight * z1
            end do
        end do

        deallocate(eigval)
        allocate(eigvec(phi_size, phi_size))
        allocate(Aw(num_non_neg_bos_freq, phi_size))
        allocate(v(phi_size))

        Aw = 0.0d0
        v = 0.0d0

        do w = 1, num_non_neg_bos_freq
            open(18, file='eigvecs_'//trim(str(w))//'.bin', form="unformatted", access="stream", status="old")
            read(18) eigvec
            close(18)

            do j = 2, phi_size
                do i = 1, phi_size
                    v(i) = v(i) + impt_smpl(w,j) * eigvec(i,j)
                end do
            end do

            do i = 1, phi_size
                Aw(w,i) = Aw(w,i) + v(i)
            end do
        end do

        deallocate(v)
        deallocate(impt_smpl)
        deallocate(eigvec)

        allocate(At(Nt,phi_size))
        At = 0.0d0

        do ti = 1, Nt
            do i = 1, phi_size
                At(ti,i) = T*real(Aw(1,i))
            end do

            do w = 2, num_non_neg_bos_freq
                do i = 1, phi_size
                    At(ti,i) = T * At(ti,i) + 2 *T* (real(Aw(w,i)) * cos(non_neg_bos_freq(w) * delta*(ti-1)) - &
                          aimag(Aw(w,i)) * sin(non_neg_bos_freq(w) * delta*(ti-1)))
                end do
            end do
        end do

        deallocate(Aw)

        ! Check for NaN
        contains_nan = .false.
        do ti = 1, Nt
            do i = 1, phi_size
                if (At(ti,i) .ne. At(ti,i)) then
                    contains_nan = .true.
                    exit
                end if
            end do
            if (contains_nan) exit
        end do

        if (contains_nan) then
            print*, 'NaN detected in At for config ', config, '- regenerating...'
            deallocate(At)
            cycle
        end if

        open(20, file='./A0/A0_'//trim(str(config))//'.dat', status="replace")
        do ti = 1,Nt
            do i = 1,phi_size
                write(20, *) At(ti,i)
            end do
        end do
        close(20)
        deallocate(At)
    end do

contains

    character(len=20) function str(k)
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end program normal_random
EOF
done

