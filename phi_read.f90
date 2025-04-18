program phi_read
        implicit none
    real*8 :: box_size_x, box_size_y, box_size_z
    integer :: N1, N2, N3, HO, KS_low, KS_high
    integer :: active_window_size, phi_size, l, x, y, z
    complex*16, allocatable :: phi_trim(:)
    complex*16, allocatable :: phi(:,:,:,:)

    ! Open and read data from files
    open(unit=12, file="delta_xyz", status="old", action="read")
    read(12, *) box_size_x, box_size_y, box_size_z
    read(12, *) N1, N2, N3
    close(12)


    open(unit=13, file="active_window", status="old", action="read")
    read(13, *) HO, KS_low, KS_high
    close(13)

    active_window_size = KS_high - KS_low + 1
    phi_size = N1*N2*N3*active_window_size

    allocate(phi_trim(phi_size))
    allocate(phi(active_window_size,N1,N2,N3))


    open(10, file="phi_trim.bin", form="unformatted", access="stream", status="old")
    read(10) phi_trim
    close(10)
    print*, "Orbitals read"

    do l = 1, active_window_size
        do x = 1, N1
            do y = 1, N2
                do z = 1, N3

                phi(l,x,y,z) = phi_trim((l-1)*N1*N2*N3 + (x-1)*N2*N3 + (y-1)*N3 + z)

                end do
            end do
        end do
    end do
         

    open(14, file="phi.bin", form="unformatted", access="stream", status="replace")
    write(14) phi
    close(14)
    print*, "phi(l,x,y,z) written as phi.bin"
     
end program phi_read
