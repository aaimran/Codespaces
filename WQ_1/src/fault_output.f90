module fault_output

  use common, only : wp
  use mpi
  use mpi3dbasic, only : pw,ps,MPI_REAL_PW,MPI_REAL_PS,error
  use mpi3dio
  use datatypes, only : fault_type

  implicit none

contains

  subroutine init_fault_output(w_fault, name, T, C, comm)

    use mpi3dcomm, only : cartesian3d_t, allocate_array_boundary

    implicit none

    character(*), intent(in) :: name

    type(fault_type), intent(inout) :: T
    type(cartesian3d_t), intent(in) :: C
    integer, intent(in) :: comm
    logical, intent(in) :: w_fault
    character(len=256),dimension(8) :: filename
    character(len=16) :: extension1,extension2,extension8
    character(len=12) :: extension3
    character(len=15) :: extension4,extension5,extension6,extension7
    character (len = 1) :: direction
    integer :: i
    extension1 = '_interface.Uface'
    extension2 = '_interface.Vface'
    extension3 = '_interface.S'
    extension4 = '_interface.Uhat'
    extension5 = '_interface.Vhat'
    extension6 = '_interface.Svel'
    extension7 = '_interface.trup'
    extension8 = '_interface.state'
    filename(1) = trim(adjustl(trim(name) // extension1))
    filename(2) = trim(adjustl(trim(name) // extension2))
    filename(3) = trim(adjustl(trim(name) // extension3))
    filename(4) = trim(adjustl(trim(name) // extension4))
    filename(5) = trim(adjustl(trim(name) // extension5))
    filename(6) = trim(adjustl(trim(name) // extension6))
    filename(7) = trim(adjustl(trim(name) // extension7))
    filename(8) = trim(adjustl(trim(name) // extension8))


    call subarray(C%nr,C%ns,C%mr,C%pr,C%ms,C%ps, MPI_REAL_PS, T%array_s)

    if(w_fault .eqv. .true.) then
       ! open files
       do i = 1,8
          call open_file_distributed(T%handles(i), filename(i), &
               'write', comm, T%array_s, ps)
       end do

    end if
    
    direction = 'q'

    call allocate_array_boundary(T%slip, C, 4, direction, ghost_nodes = .true.)
    call allocate_array_boundary(T%time_rup, C, 1, direction, ghost_nodes = .true.)
    call allocate_array_boundary(T%W, C, 1, direction, ghost_nodes = .true.)
    call allocate_array_boundary(T%Uhat_pluspres, C, 6, direction, ghost_nodes = .true.)
    call allocate_array_boundary(T%Vhat_pluspres, C, 6, direction, ghost_nodes = .true.)
    call allocate_array_boundary(T%U_pluspres, C, 9, direction, ghost_nodes = .true.)
    call allocate_array_boundary(T%V_pluspres, C, 9, direction, ghost_nodes = .true.)

  end subroutine init_fault_output


  subroutine write_fault(U,S,W,fault)

    use mpi
    use mpi3dbasic, only: rank

    implicit none

    real(kind = wp), dimension(:,:,:), intent(in) :: U
    real(kind = wp), dimension(:,:,:), intent(in) :: S
    real(kind = wp), dimension(:,:,:), intent(in) :: W
    type(fault_type), intent(in) :: fault
    integer :: usize

    ! write data to file Uface
    usize = size(U,1)
    call write_file_distributed(fault%handles(1), U)

    ! write data to file S
    call write_file_distributed(fault%handles(3), S)

    ! write data to file W
    call write_file_distributed(fault%handles(8), W)

  end subroutine write_fault


  subroutine write_hats(Uhat,Vhat,Svel,trup,fault)

    use mpi

    implicit none

    real(kind = wp), dimension(:,:,:), intent(in) :: Uhat,Vhat
    real(kind = wp), dimension(:,:,:), intent(in) :: Svel
    real(kind = wp), dimension(:,:), intent(in) :: trup
    type(fault_type), intent(in) :: fault

    ! write data to file Uhat
    call write_file_distributed(fault%handles(4), Uhat)

    ! write data to file 2
    call write_file_distributed(fault%handles(5), Vhat)

    ! write data to file S vel
    call write_file_distributed(fault%handles(6), Svel)

    ! write data to file trup
    call write_file_distributed(fault%handles(7), trup)

  end subroutine write_hats


  subroutine destroy_fault(fault)


    integer :: i
    type(fault_type), intent(inout) :: fault

    do i = 1,8
       call close_file_distributed(fault%handles(i))
    end do

  end subroutine destroy_fault


end module fault_output
