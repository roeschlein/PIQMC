program density_check
        implicit none
        integer :: active_window_size, HO_index, low_index, high_index
        integer :: i, j, k, l
        integer :: N1, N2, N3
        real*8 :: box_size_x, box_size_y, box_size_z
        complex*16, dimension(:,:,:,:), allocatable :: phi
        real*8, dimension(:), allocatable :: rho_x, rho_y, rho_z, rho

        ! Read the grid dimensions from delta_xyz
        open(unit=12, file="delta_xyz", status="old", action="read")
        read(12, *) box_size_x, box_size_y, box_size_z
        read(12, *) N1, N2, N3
        close(12)

         ! Read the grid dimensions from delta_xyz
        open(unit=13, file="active_window", status="old", action="read")
        read(13, *) HO_index
        read(13, *) low_index
        read(13, *) high_index
        close(13)

        active_window_size = high_index - low_index + 1

        ! Read orbitals in active window
        allocate(phi(active_window_size,N1,N2,N3))
        open(10,file="phi.bin", form="unformatted", access="stream", status="old")
        read(10) phi
        close(10)

        allocate(rho_x(N1))
        allocate(rho_y(N2))
        allocate(rho_z(N3))


        do l = 1, HO_index-low_index 
        do i = 1, N1 
              do j = 1, N2 
                  do k = 1, N3 

                      rho_x(i) = rho_x(i) + phi(l,i,j,k) * CONJG(phi(l,i,j,k))
                  
                      rho_y(j) = rho_y(j) + phi(l,i,j,k) * CONJG(phi(l,i,j,k)) 
                  
                      rho_z(k) = rho_z(k) + phi(l,i,j,k) * CONJG(phi(l,i,j,k))
                  
                  end do 
              end do 
        end do
        end do


        do i = 1, N1
           write(*,'("x{",I0,"}: ",F10.6)') i, rho_x(i)/MAXVAL(rho_x)
        end do

        do j = 1, N2
           write(*,'("y{",I0,"}: ",F10.6)') j, rho_y(j)/MAXVAL(rho_y)
        end do

        do k = 1, N3
           write(*,'("z{",I0,"}: ",F10.6)') k, rho_z(k)/MAXVAL(rho_z)
        end do
end program density_check
