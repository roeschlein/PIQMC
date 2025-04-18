program inverse_propagator
        implicit none
        real*8, parameter :: PI = acos(-1.0d0)
        integer :: num_fermi_freq, active_window_size,i,j,HO_index, num_fermi_pairs
        integer :: HO, KS_low, KS_high
        complex*16, dimension(:,:), allocatable :: M
        real*8, dimension(:), allocatable :: fermi_freq
        real*8, dimension(:,:), allocatable :: KS_eigenvalues
        real*8 :: T, temp
       
        ! Read the grid dimensions from delta_xyz
        open(unit=13, file="active_window", status="old", action="read")
        read(13, *) HO, KS_low, KS_high
        close(13)

        active_window_size = KS_high - KS_low + 1


        ! Read mat freq stuff
        open(unit=14, file="matsubara_params", status="old", action="read")
        read(14, *) num_fermi_freq
        read(14, *) num_fermi_pairs
        read(14, *) T
        close(14)

        ! This allows us to be symmetric about 0
        allocate(fermi_freq(num_fermi_freq))
        do i = 1, num_fermi_freq/2
                fermi_freq(2*i-1) = (2*i-1)*PI*T !Fermionic Matsubara Frequencies
                fermi_freq(2*i) = -(2*i-1)*PI*T
        end do
        
            ! Bubble Sort Algorithm Sort the Fermi Mat Freq
    do i = 1, num_fermi_freq-1
        do j = 1, num_fermi_freq-i
            if (fermi_freq(j) > fermi_freq(j+1)) then
                temp = fermi_freq(j)
                fermi_freq(j) = fermi_freq(j+1)
                fermi_freq(j+1) = temp
            end if
        end do
    end do
    

    ! Read in Kohn Sham Eigenvalues in the active window
    allocate(KS_eigenvalues(3,active_window_size))
    open(10, file="energy_pop_shifted.dat")
    do i = 1, active_window_size
        read(10,*) KS_eigenvalues(1,i), KS_eigenvalues(2,i), KS_eigenvalues(3,i) ! index eigenvalue occupation
        KS_eigenvalues(2,i) = KS_eigenvalues(2,i)/T ! /tilde{/epsilon}_i = \beta(\epsilon_i - \mu) 
    end do
    close(10)


    
    
    ! Fill Inv Propagator Matrix
    allocate(M(num_fermi_freq, active_window_size))

        do i = 1, num_fermi_freq
                do j = 1, active_window_size
                        M(i,j) = 1/(-complex(0,1)*fermi_freq(i)+KS_eigenvalues(2,j))
                end do
        end do

   open(12, file="M.bin", form="unformatted", access="stream", status="replace")
   write(12) M
   close(11)

   print*, "Inverse propagators written to M(i,j);i=fermionic matsubara frequency index,j=orbital index"

end program inverse_propagator
        

