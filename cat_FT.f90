program read_all_evcr
    implicit none

    integer :: N1, N2, N3, HO_index
    real*8 :: box_size_x, box_size_y, box_size_z
    integer :: starting_index, ending_index, i
    integer*8 :: num_elements, total_elements 
    integer :: num_files, io_status
    character(len=50) :: filename
    complex*16, dimension(:), allocatable :: evcr_all
    complex*16, dimension(:, :, :), allocatable :: evcr_data

    ! Read the grid dimensions from delta_xyz
    open(unit=12, file="delta_xyz", status="old", action="read")
    read(12, *) box_size_x, box_size_y, box_size_z
    read(12, *) N1, N2, N3
    close(12)

    allocate(evcr_data(N1, N2, N3))

    open(unit=17, file="active_window", status="old", action="read")
    read(17, *) HO_index
    read(17, *) starting_index
    read(17, *) ending_index
    close(17)

    num_files = ending_index - starting_index + 1

    ! Calculate total elements
    num_elements = N1 * N2 * N3
    total_elements = num_files * num_elements

    ! Allocate memory for the combined array
    allocate(evcr_all(total_elements))

    ! Loop over each file and read data into evcr_all
    do i = starting_index, ending_index
        ! Construct filename without zero-padding
        write(filename, '("evcr", I0, ".bin")') i

        ! Open the binary file for reading
        open(unit=13, file=filename, form="unformatted", access="stream", status="old", action="read", iostat=io_status)
        if (io_status /= 0) then
            print*, "Warning: Could not open file ", trim(filename)
            cycle
        end if

        ! Read the complex data
        read(13) evcr_data
        close(13)

        if ((i - starting_index) * num_elements + 1 > total_elements) then
        print*, "Error: Out-of-bounds access in evcr_all at index", i
        print*, (i - starting_index) * num_elements + 1
        stop
        end if
        
        ! Copy into evcr_all while maintaining memory layout
        evcr_all((i - starting_index) * num_elements + 1:(i - starting_index + 1) &
               * num_elements) = reshape(evcr_data, [num_elements])
    end do

    ! Save the combined data to phi.bin
    open(unit=14, file="phi_temp.bin", form="unformatted", access="stream", status="replace")
    write(14) evcr_all
    close(14)

    print*, "All the real-space KS orbitals in active window have been successfully combined into phi_temp.bin"

end program read_all_evcr

