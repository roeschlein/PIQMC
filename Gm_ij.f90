program Gm_ij
        implicit none
        integer :: active_window_size, HO, KS_high, KS_low
        integer*8 :: phi_size
        integer :: num_non_neg_bos_freq, num_fermi_pairs, i, j, wm, wf
        real*8, dimension(:,:), allocatable :: fermi_pairs
        integer, dimension(:), allocatable :: fermi_pairs_len
        integer, dimension(:,:), allocatable :: fermi_pair_raw
        complex*16, dimension(:,:), allocatable :: M
        complex*16, dimension(:,:,:), allocatable :: G
        real*8 :: T
        


        open(unit=10, file="active_window", status="old", action="read")
        read(10, *) HO, KS_low, KS_high
        close(10)

        active_window_size = KS_high - KS_low + 1

        ! Read mat freq stuff
        open(unit=15, file="matsubara_params", status="old", action="read")
        read(15, *) num_non_neg_bos_freq
        read(15, *) num_fermi_pairs
        read(15, *) T
        close(15)


        ! Read in inverse propagators
        allocate(M(num_non_neg_bos_freq,active_window_size))
        open(13, file="M.bin", form="unformatted", access="stream", status="old")
        read(13) M
        close(13)
        print*, "Inverse propagators read"


        ! Read in file with how many pairs per bosonic frequency
        allocate(fermi_pairs_len(num_non_neg_bos_freq))
        open(12, file="fermi_freq_pairs_len")
        do i = 1, num_non_neg_bos_freq
                read(12,*) fermi_pairs_len(i)
        end do
        close(12)

        ! Read in fermi pairs (wp,wk) s.t.  wp-wk=wn for each bosonic frequency
        allocate(fermi_pair_raw(2,num_fermi_pairs))

        open(11, file="fermi_freq_pairs")
        do i = 1, num_fermi_pairs
                read(11,*) fermi_pair_raw(1,i), fermi_pair_raw(2,i)
        end do
        close(11)
        print*, "Fermionic Matsubara frequency pairs read [(wp,wk) s.t. wp-wk=wn for each bosonic frequency]"


        allocate(G(num_non_neg_bos_freq, active_window_size, active_window_size))
        do wm = 1, num_non_neg_bos_freq
            do wf = 1, fermi_pairs_len(wm)
                do i = 1, active_window_size
                do j = 1, active_window_size
                G(wm,i,j) = G(wm,i,j) + M(fermi_pair_raw(1,wf),i)*M(fermi_pair_raw(2,wf),j)
                end do
                end do
            end do
        end do
        print*,"Gm_ij filled"

        open(14, file="Gm_ij.bin", form="unformatted", access="stream", status="replace")
        write(14) G
        close(14)
        print*, "Gm_ij written as Gm_ij.bin"


end program Gm_ij
