!> @brief domain module holds multiple blocks and interfaces, and is used by time-stepping routines

!> @todo Hari, I've never been careful about public vs private in modules, maybe you can help with that?

!> domain = set of blocks and interfaces, with fields defined at time t
!> blocks and interfaces are allocatable to allow for any number of them,
!> but initialization is specific to two blocks and one interface (which I think
!> should be the initial focus since that is sufficient to do many interesting problems),
!> but we should generalize this to an arbitrary number of blocks and interfaces
!> with an unstructured connectivity pattern (presumably then the domain_type would also
!> contain some type that would handle the block/interface layout)

module domain

  use mpi
  use common, only : wp
  use datatypes, only : domain_type, block_type, iface_type, mms_type,fault_type,moment_tensor
  use block, only : block_temp_parameters, eval_block_mms, norm_fields_block
  use slice_output, only : init_slice_output, end_slice_output
  use seismogram
  use plane_output, only : init_plane_output, write_plane_output, end_plane_output

  implicit none

  logical :: in_block_comm(2), in_fault_comm(2), in_hslice_comm(2), in_vslice_comm(2)

contains


  subroutine init_domain(D, infile, input_filename)

    !> initialize the domain at t=0

    use block, only : init_block, block_time_step
    use iface, only : init_iface
    use fault_output, only : init_fault_output
    use mms, only : init_mms
    use interface_condition, only : init_vel_state
    use mpi3dbasic, only: is_master, rank, nprocs, new_communicator, warning, print_run_info
    use mpi3d_interface
    use, intrinsic :: iso_c_binding, only: c_int

    implicit none

    interface
      subroutine c_exit(status) bind(C, name="exit")
        import :: c_int
        integer(c_int), value :: status
      end subroutine c_exit
    end interface

    type(domain_type),intent(out) :: D
    integer, intent(in) :: infile
    character(*), intent(in), optional :: input_filename


    integer :: i, cart_size(3)
    integer :: n, ierr, cart_rank, comm_cart, coord(3), normal(3), &
               block_comms(2), fault_comms(2), hslice_comms(2), &
               vslice_comms(2), my_block_comm
    logical :: periodic(3) = (/ .false., .false., .false. /), reorder=.true.
    type(interface3d) :: II

    ! Variables for parsing the input data file
    integer :: stat
    character(256) :: name, problem,response,plastic_model,fd_type
    character(256) :: response_norm
    logical :: invalid_response
    integer :: nt, nblocks
    logical :: output_fault_topo, w_fault, interpol, use_topography, mollify_source
    integer :: w_stride, ny, nz, order
    real(kind = wp) :: CFL, t_final, topo
    ! interface conditions
    character(64) :: coupling !< locked, slip-weakening_friction, linear_friction
    character(64) :: mesh_source, type_of_mesh, material_source !< cartesian or curvilinear

    real(kind = wp) :: dt, dtmin, dt1, dt2, spatialsetep1(3), spatialsetep2(3)
    real(kind = wp) :: dt_i, spat(3)
    
    real(kind = wp) :: meshvolume1, meshvolume2, totalvolume, ratio1, ratio2
    integer :: nprocs_1, nprocs_2

    type(block_temp_parameters) :: btp(2)


    ! namelists to read in parameters

    namelist /problem_list/ name, problem,response,plastic_model,nblocks, &
                            nt, CFL, coupling, fd_type, order, t_final, &
                            mesh_source, type_of_mesh,material_source, &
                            interpol, w_stride, w_fault, use_topography, topo, mollify_source

    namelist /block_list/ btp


    !---------------------------------------------------------------------------
    !                       END OF VARIABLES
    !---------------------------------------------------------------------------
    ! set default parameters and read problem list

    name = 'default'
    problem = 'TPV5'
    response = 'elastic'
    plastic_model = 'default'
    nblocks = 2          ! NEW: Default to 2 blocks for backward compatibility
    nt = 0
    CFL = 0.5_wp
    coupling = 'locked'
    t_final = 0.0_wp
    w_fault = .true.
    output_fault_topo = .true.
    mesh_source = 'compute'
    material_source = 'hardcode'
    w_stride = 1
    interpol = .false.
    use_topography = .false.
    mollify_source = .false.
    topo = 1.0_wp
    fd_type = 'traditional'
    order = 5
    

    read(infile,nml=problem_list,iostat=stat)
    if (stat>0) stop 'error reading namelist problem_list'

     ! Validate `response` only on the master rank, then broadcast to all ranks.
     invalid_response = .false.
     if (is_master()) then
       response_norm = trim(adjustl(response))
       select case (response_norm)
       case ('elastic','plastic','anelastic','low-pass','anelastic-Q','anelastic-Q8','anelastic-Qf','constant-Q-4M', &
         'constant-Q-8M','frequency-Q-4M','frequency-Q-8M')
         ! ok
       case default
         invalid_response = .true.
         write(*,*) ':: Error: invalid response option in &problem_list: "' // trim(response_norm) // '"'
         write(*,*) ':: Supported response options are: elastic, plastic, anelastic, low-pass, anelastic-Q, anelastic-Q8, anelastic-Qf, constant-Q-4M, constant-Q-8M, frequency-Q-4M, frequency-Q-8M'
       end select
     end if

     call MPI_BCAST(invalid_response, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(response_norm, len(response_norm), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     if (invalid_response) then
       call MPI_FINALIZE(ierr)
       call c_exit(1_c_int)
     end if

    if (is_master() .and. present(input_filename)) then
      call print_run_info(input_filename, name, problem)
      write(*,'(/,/)')
    end if

    D%name = name
    D%problem = problem
    D%response = response_norm
    D%plastic_model = plastic_model
    D%coupling = coupling
    D%w_fault = w_fault
    D%output_fault_topo = output_fault_topo
    D%t_final = t_final
    D%type_of_mesh = type_of_mesh
    D%mesh_source = mesh_source
    D%material_source = material_source
    D%CFL = CFL
    D%w_stride = w_stride
    D%interpol = interpol
    D%use_topography = use_topography
    D%topo = topo
    D%mollify_source = mollify_source 
    D%fd_type = fd_type
    D%order = order
    D%debug = .false.  ! Default off; enable for diagnostics when needed




    rewind(infile)
    read(infile,nml=block_list,iostat=stat)
    if (stat>0) stop 'error reading namelist block_list'


    
       

    !> set initial time

    D%t = 0.0_wp

    ! compute global time-step using only the blocks specified by `nblocks`
     dtmin = huge(1.0_wp)
     do i = 1, max(1,nblocks)
       spat(1) = (btp(i)%bqrs(1)-btp(i)%aqrs(1))/real(btp(i)%nqrs(1)-1)
       spat(2) = (btp(i)%bqrs(2)-btp(i)%aqrs(2))/real(btp(i)%nqrs(2)-1)
       spat(3) = (btp(i)%bqrs(3)-btp(i)%aqrs(3))/real(btp(i)%nqrs(3)-1)
       dt_i = block_time_step(spat, CFL, btp(i)%rho_s_p)
       dtmin = min(dtmin, dt_i)
     end do

     D%dt = dtmin

        if (is_master()) call warning('warning: current method for setting time step does not use material properties; ' // &
          'mesh information from files is ignored (only scalars from btp%... in input file are used)', 'init_domain')
    
    !> set number of blocks and interfaces
    !> (now supports 1 or 2 blocks; read from input file)
    
    ! Validate nblocks
    if (nblocks < 1 .or. nblocks > 2) then
      if (is_master()) call warning('Invalid nblocks value; must be 1 or 2. Using default nblocks = 2', 'init_domain')
      nblocks = 2
    end if

    D%nblocks = nblocks
    D%nifaces = max(0, nblocks - 1)  ! 0 for 1-block, 1 for 2-blocks

    ! MPI distribution: conditional on nblocks
    if (D%nblocks == 1) then
      ! 1-block mode: all ranks in block 1
      in_block_comm(1) = .true.
      in_block_comm(2) = .false.
    else
      ! 2-block mode: split ranks 50/50 (can be generalized with work ratio)
      ! compute the volume of work in each block
      meshvolume1 = real(btp(1)%nqrs(1)*btp(1)%nqrs(2)*btp(1)%nqrs(3))
      meshvolume2 = real(btp(2)%nqrs(1)*btp(2)%nqrs(2)*btp(2)%nqrs(3))

      ! compute the work ratio
      totalvolume =  meshvolume1 +  meshvolume2
      ratio1 =  meshvolume1/totalvolume
      ratio2 =  meshvolume2/totalvolume

      ! assign processes using the work ratio
      nprocs_1 =  int((nprocs + 1)*ratio1)
      nprocs_2 =  int((nprocs + 1)*ratio2)

      ! Currently using simple 50/50 split
      ! Can be generalized by using ratio1, ratio2 instead:
      ! in_block_comm(1) = (rank <= nprocs_1)
      ! in_block_comm(2) = (rank > nprocs_1)
      in_block_comm(1) = (rank < (nprocs + 1)/2)
      in_block_comm(2) = (rank >= nprocs/2)
    end if

    !in_block_comm(1) = (rank .le. nprocs_1)
    !in_block_comm(2) = (rank > nprocs_1)

    call new_communicator(in_block_comm(1), block_comms(1))
    call new_communicator(in_block_comm(2), block_comms(2))

    if (in_block_comm(1)) my_block_comm = block_comms(1)
    if (in_block_comm(2)) my_block_comm = block_comms(2)

    allocate(D%B(D%nblocks))
    allocate(D%I(D%nifaces))
    allocate(D%seismometers(D%nblocks))
    allocate(D%plane_outputs(D%nblocks))

    !> initialize blocks
    do i = 1,D%nblocks
       ny = btp(i)%nqrs(2)
       nz = btp(i)%nqrs(3)

    if((btp(i)%profile_type == 'read_from_memomry_fractal') .and. &
         (btp(i)%topography_type == 'read_topo_from_memory')) then

       ny = btp(i)%faultsize(1)
       nz = btp(i)%faultsize(2)

    end if
       
      if (.not.in_block_comm(i)) cycle
      call init_block(D%mesh_source, D%type_of_mesh, D%material_source,&
           D%response, D%fd_type,  D%order, D%interpol, D%use_topography, topo, D%B(i), &
         problem, btp(i), block_comms(i),infile,i, ny, nz, D%debug)      
  
      call MPI_Comm_size(block_comms(i),n,ierr)
      cart_size = (/0,0,0/)
      call MPI_Dims_create(n,3,cart_size,ierr)      

      call MPI_Cart_create(block_comms(i),3,cart_size,periodic,reorder,comm_cart,ierr)
      call MPI_Comm_rank(comm_cart,cart_rank,ierr)
      call MPI_Cart_coords(comm_cart,cart_rank,3,coord,ierr)

      
    end do


    D%nt = floor(D%t_final/dtmin)
   
    if (is_master()) then
      write (*,*) name  
      write (*,*) "Domain time parameters :"
      write (*,*) "        Number of time steps = ", D%nt
      write (*,*) "        Initial dt = ", D%dt
      write (*,*) "        Final time = ", D%t_final
    end if

    ! Setup communicator for interface partitioning.
    ! Only create interface communicators for 2-block mode (when nifaces > 0)

    if (D%nifaces > 0) then
      if (in_block_comm(1)) then
        normal = [1, 0, 0]
        call new_interface(coord, cart_size, normal, MPI_COMM_WORLD, D%B(1)%G%C, II)
      end if
      if (in_block_comm(2).and. .not.in_block_comm(1)) then
        normal = [-1, 0, 0]
        call new_interface(coord, cart_size, normal, MPI_COMM_WORLD, D%B(2)%G%C, II)
      end if
    end if

    !> initialize interfaces (only if nifaces > 0, i.e., for 2-block mode)

    do i = 1, D%nifaces
       !> create an interface that joins two blocks
       !> first define which blocks will be coupled: block 1 to block 2
       D%I(i)%im = 1
       D%I(i)%ip = 2
       !> and in which direction they will be coupled: q
       D%I(i)%direction = 'q'
       !> then initialize interface
       if (in_block_comm(1)) call init_iface(D%I(i),D%B(1), II)
       if (in_block_comm(2) .and. .not.in_block_comm(1)) call init_iface(D%I(i),D%B(2), II)
       if (in_block_comm(1)) call init_vel_state(D%problem, D%I(i), D%B(1), -1.0_wp, 2)
       if (in_block_comm(2) .and. .not.in_block_comm(1)) call init_vel_state(D%problem, D%I(i), D%B(2), 1.0_wp, 1)
    end do

    ! Exchange materials across interface (only for 2-block mode)
    if (D%nifaces > 0) call exchange_materials_interface(D)

    ! Fault output initialization (only for 2-block mode with interface)
    if (D%nifaces > 0) then
      in_fault_comm(1) = (in_block_comm(1) .and. D%I(1)%II%on_interface)
      in_fault_comm(2) = (in_block_comm(2) .and. D%I(1)%II%on_interface)
      call new_communicator(in_fault_comm(1), fault_comms(1))
      call new_communicator(in_fault_comm(2), fault_comms(2))

      if (in_fault_comm(1)) call init_fault_output(D%w_fault,D%name, D%fault, D%B(1)%G%C, fault_comms(1))
      if (in_fault_comm(2) .and. .not.in_fault_comm(1)) call init_fault_output(D%w_fault,D%name, D%fault, D%B(2)%G%C, fault_comms(2))
    else
      ! 1-block mode: no interface communicators
      in_fault_comm(1) = .false.
      in_fault_comm(2) = .false.
    end if

    call init_mms(infile, D%mms_vars)

    ! Initialize seismometers and plane outputs for each block (generalized for 1 or 2 blocks)
    do i = 1, D%nblocks
      D%seismometers(i)%block_num = i
      if (in_block_comm(i)) call init_seismogram(infile, D%seismometers(i), D%name, D%B(i)%G)
      if (in_block_comm(i)) call init_plane_output(infile, D%name, D%plane_outputs(i), D%B(i)%G, i)
    end do
    
    !> above assumes that blocks contain fields defined on a structured mesh,
    !> which is a cube in the computational domain; the computational coordinates
    !> are q,r,s, and there is a mapping to a curvilinear mesh in physical coordinates x,y,z
    !> I recommend that we use q,r,s in all references to directions, sides, etc., because
    !> it is possible that a boundary in q is mapped to a boundary in y or z rather than x
    !> (this will be a different naming convention than in original version of 3D code)

  end subroutine init_domain


  subroutine close_domain(D)

    use fault_output, only: destroy_fault

    type(domain_type),intent(inout) :: D

    if ( D%w_fault .eqv.  .true.) then
       if (in_fault_comm(1)) call destroy_fault(D%fault)
       if (in_fault_comm(2) .and. .not.in_fault_comm(1)) call destroy_fault(D%fault)
    end if
    !call end_slice_output(D%slicer)
    if (in_block_comm(1)) call destroy_seismogram(D%seismometers(1))
    if (in_block_comm(2)) call destroy_seismogram(D%seismometers(2))

    if (in_block_comm(1)) call end_plane_output(D%plane_outputs(1))
    if (in_block_comm(2)) call end_plane_output(D%plane_outputs(2))

  end subroutine close_domain

  subroutine enforce_bound_iface_conditions(D, stage)

    !> enforce boundary and interface conditions

    use block, only : enforce_bound_conditions

    implicit none

    type(domain_type),intent(inout) :: D
    integer, intent(in) :: stage

    !> enforce boundary conditions on external sides of each block
    !> enforce interface conditions to couple blocks (only for 2-block mode)

    if (in_block_comm(1)) then
      call enforce_bound_conditions(D%B(1), D%mms_vars, D%t)
      if (D%nifaces > 0) then
        call enforce_iface_conditions(D%problem, D%coupling, D%I(1), &
          D%B(1),2,D%t, stage, D%mms_vars, D%fault)
      end if
    end if
    if (in_block_comm(2)) then
      call enforce_bound_conditions(D%B(2), D%mms_vars, D%t)
      if (D%nifaces > 0) then
        call enforce_iface_conditions(D%problem, D%coupling, D%I(1), &
          D%B(2),1,D%t, stage, D%mms_vars, D%fault)
      end if
    end if

  end subroutine enforce_bound_iface_conditions

  subroutine enforce_iface_conditions(problem, coupling, I,B, ib, t, stage, mms_vars, handles)

    use datatypes, only : block_type, iface_type, mms_type, fault_type, block_boundary
    use RHS_Interior, only : Impose_Interface_Condition
    use mpi3dbasic, only : rank

    implicit none

    character(*), intent(in) :: problem, coupling
    type(iface_type),intent(inout) :: I
    type(block_type),intent(inout) :: B
    type(mms_type), intent(inout) :: mms_vars
    type(fault_type), intent(inout) :: handles
    real(kind = wp),intent(in) :: t
    integer, intent(in) :: stage, ib

    select case(I%direction)

    case('q')

!!! coupling: side 2 of block 1, B(2) <==> side 1 of block 2, Bp(1)

!!! this solves interface conditions for hat variables, constructs SAT forcing terms,
!!! and adds SAT forcing terms to rates

      if (B%boundary_vars%Rx == 0 .or. B%boundary_vars%Lx == 0) then
        call Impose_Interface_Condition(problem, coupling, I, B, ib, &
                              t, stage, mms_vars, handles)
      end if

    case('r','s')

       stop 'interfaces in r and s direction not implemented in enforce_iface_conditions'

    end select

  end subroutine enforce_iface_conditions

  subroutine exchange_fields(D)

    use block, only : exchange_fields_block

    implicit none

    type(domain_type),intent(inout) :: D
    integer :: i

    if (in_block_comm(1)) call exchange_fields_block(D%B(1))
    if (in_block_comm(2)) call exchange_fields_block(D%B(2))

  end subroutine exchange_fields

  subroutine exchange_fields_interface(D)

    use mpi3dbasic, only : nprocs
    use block, only : copy_fields_to_boundary
    use boundary, only : exchange_fields_across_interface

    implicit none
    type(domain_type),intent(inout) :: D

    ! Skip interface field exchange if no interfaces exist (1-block mode)
    if (D%nifaces == 0) return

    if (nprocs == 1) then
      call copy_fields_to_boundary(D%B(1))
      if (D%nblocks >= 2) then
        call copy_fields_to_boundary(D%B(2))
        D%B(1)%B(2)%Fopp(:,:,:) = D%B(2)%B(1)%F(:,:,:)
        D%B(2)%B(1)%Fopp(:,:,:) = D%B(1)%B(2)%F(:,:,:)
      end if
    else
      if (in_block_comm(1)) then
        call copy_fields_to_boundary(D%B(1))
        call exchange_fields_across_interface(D%B(1)%B(2), D%B(1)%G%C, D%I(1)%II)
      end if
      if (in_block_comm(2)) then
        call copy_fields_to_boundary(D%B(2))
        call exchange_fields_across_interface(D%B(2)%B(1), D%B(2)%G%C, D%I(1)%II)
      end if
    end if

  end subroutine exchange_fields_interface

  subroutine exchange_materials_interface(D)

    use mpi3dbasic, only : nprocs
    use block, only : copy_fields_to_boundary
    use boundary, only : exchange_materials_across_interface

    implicit none
    type(domain_type),intent(inout) :: D

    ! Skip interface material exchange if no interfaces exist (1-block mode)
    if (D%nifaces == 0) return

    if (nprocs == 1) then
      if (D%nblocks >= 2) then
        D%B(1)%B(2)%Mopp = D%B(2)%B(1)%M
        D%B(2)%B(1)%Mopp = D%B(1)%B(2)%M
      end if
    else
      if (in_block_comm(1)) then
        call exchange_materials_across_interface(D%B(1)%B(2), D%B(1)%G%C, D%I(1)%II)
      end if
      if (in_block_comm(2)) then
        call exchange_materials_across_interface(D%B(2)%B(1), D%B(2)%G%C, D%I(1)%II)
      end if
    end if

  end subroutine exchange_materials_interface

  subroutine write_output(D)

    !> @brief write fields (and, in some cases, rates) at time t
    !> note that the way it is written exposes the details of how fields are stored
    !> in blocks and on interfaces; a potentially better way is to use subroutines
    !> that retrieve certain fields and return them in an output array

    use,intrinsic :: iso_fortran_env, only : output_unit

    implicit none

    type(domain_type),intent(inout) :: D

    if ( D%w_fault .eqv.  .true.) then
       call write_fault_output(D)
       !call write_slice_output(D)
       call write_hat_output(D)
    end if

     if (in_block_comm(1)) call write_plane_output(D%plane_outputs(1), D%t, D%B(1)%F%F, D%B(1)%G)
     if (in_block_comm(2)) call write_plane_output(D%plane_outputs(2), D%t, D%B(2)%F%F, D%B(2)%G)
    
    call write_seismogram_output(D)

  end subroutine write_output

  subroutine write_fault_output(D)

    use fault_output
    use mpi3dbasic, only : rank
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1, mq2, mr2, ms2, pq2, pr2, ps2

    if ( D%w_fault .eqv.  .true.) then
       if (in_fault_comm(1)) then
          
          mq1 = D%B(1)%G%C%mq
          mr1 = D%B(1)%G%C%mr
          ms1 = D%B(1)%G%C%ms
          pq1 = D%B(1)%G%C%pq
          pr1 = D%B(1)%G%C%pr
          ps1 = D%B(1)%G%C%ps
          
          call write_fault(D%B(1)%F%F(pq1,mr1:pr1,ms1:ps1,:),  &
               D%I(1)%S(mr1:pr1,ms1:ps1,:), &
               D%I(1)%W(mr1:pr1,ms1:ps1,:), D%fault)
          
       end if
       if (in_fault_comm(2)) then
          
          mq2 = D%B(2)%G%C%mq
          mr2 = D%B(2)%G%C%mr
          ms2 = D%B(2)%G%C%ms
          pq2 = D%B(2)%G%C%pq
          pr2 = D%B(2)%G%C%pr
          ps2 = D%B(2)%G%C%ps
          
          call write_file_distributed(D%fault%handles(2), D%B(2)%F%F(mq2,mr2:pr2,ms2:ps2,:))
       end if
    end if


  end subroutine write_fault_output

  subroutine write_hat_output(D)
    use fault_output
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1

    if ( D%w_fault .eqv.  .true.) then

       if (.not.in_fault_comm(1)) return
       
       mq1 = D%B(1)%G%C%mq
       mr1 = D%B(1)%G%C%mr
       ms1 = D%B(1)%G%C%ms
       pq1 = D%B(1)%G%C%pq
       pr1 = D%B(1)%G%C%pr
       ps1 = D%B(1)%G%C%ps
       
       call write_hats(D%fault%Uhat_pluspres(mr1:pr1,ms1:ps1,:), &
            D%fault%Vhat_pluspres(mr1:pr1,ms1:ps1,:), &
            D%fault%Uhat_pluspres(mr1:pr1,ms1:ps1,1:3), &
            D%fault%time_rup(mr1:pr1,ms1:ps1,1),D%fault)
       
    end if
       
     end subroutine write_hat_output

  subroutine write_slice_output(D)

    use slice_output
    implicit none

    type(domain_type), intent(in) :: D

    if (in_block_comm(1)) call write_slice(D%B(1)%F%F,D%B(1)%G%C, D%slicer)

  end subroutine write_slice_output

  subroutine write_seismogram_output(D)

    use slice_output
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1, mq2, mr2, ms2, pq2, pr2, ps2

    if(in_block_comm(1)) then
        mq1 = D%B(1)%G%C%mq
        mr1 = D%B(1)%G%C%mr
        ms1 = D%B(1)%G%C%ms
        pq1 = D%B(1)%G%C%pq
        pr1 = D%B(1)%G%C%pr
        ps1 = D%B(1)%G%C%ps

        call write_seismogram(D%seismometers(1), D%t, D%B(1)%F%F)

    end if

    if(in_block_comm(2)) then
        mq2 = D%B(2)%G%C%mq
        mr2 = D%B(2)%G%C%mr
        ms2 = D%B(2)%G%C%ms
        pq2 = D%B(2)%G%C%pq
        pr2 = D%B(2)%G%C%pr
        ps2 = D%B(2)%G%C%ps

        call write_seismogram(D%seismometers(2), D%t, D%B(2)%F%F) 

    end if


  end subroutine write_seismogram_output

  subroutine eval_mms(D)

    implicit none
    type(domain_type),intent(inout) :: D

    integer :: i

    do i = 1, D%nblocks
      if (.not.in_block_comm(i)) cycle
       call eval_block_mms(D%B(i), D%t, D%mms_vars)
    end do

  end subroutine eval_mms

  subroutine norm_fields(D)

    implicit none
    type(domain_type),intent(inout) :: D

    integer :: i

     if (in_block_comm(1)) call norm_fields_block(D%B(1))
     if (in_block_comm(2)) call norm_fields_block(D%B(2))

  end subroutine norm_fields

  subroutine scale_rates(D,A)

    !> @brief multiply all rates by RK coefficient A

    use block, only : scale_rates_block
    use iface, only : scale_rates_iface

    implicit none

    type(domain_type),intent(inout) :: D
    real(kind = wp),intent(in) :: A

    integer :: i

    ! first within blocks and on their boundaries

      if (in_block_comm(1)) call scale_rates_block(D%B(1),A)
      if (in_block_comm(2)) call scale_rates_block(D%B(2),A)


    ! and then on interfaces

    do i = 1,D%nifaces
       call scale_rates_iface(D%I(i),A)
    end do

  end subroutine scale_rates

  subroutine set_rates(D)

    !> @brief set rates using the PDE
    !> @details set rates using the PDE (in a low storage RK method, the new rates
    !> at the current stage are added to the old rates, instead of overwriting
    !> the rates array)

    use block, only : set_rates_block
    use moment_tensor, only : set_moment_tensor, moment_tensor_body_force, set_moment_tensor_smooth

    implicit none

    type(domain_type),intent(inout) :: D

       if (in_block_comm(1)) call set_rates_block(D%B(1), D%type_of_mesh)
       if (in_block_comm(2)) call set_rates_block(D%B(2), D%type_of_mesh)


    ! Guard with in_block_comm: D%B(2) is out of bounds in 1-block mode, and
    ! D%B(i)%MT is uninitialized on ranks that do not own block i.
    if (in_block_comm(1)) then
    if (D%B(1)%MT%use_moment_tensor) then
       !call set_moment_tensor(D%B(1),D%t)
       if(D%mollify_source) then

          call set_moment_tensor_smooth(D%B(1),D%t)

       else

          call set_moment_tensor(D%B(1),D%t)
       end if

       !call moment_tensor_body_force(D%B(1),D%t)
    end if
    end if

    if (in_block_comm(2)) then
    if (D%B(2)%MT%use_moment_tensor) then
       !call set_moment_tensor(D%B(2),D%t)

       if(D%mollify_source) then

          call set_moment_tensor_smooth(D%B(2),D%t)

       else

          call set_moment_tensor(D%B(2),D%t)
       end if

       !call set_moment_tensor2(D%B(2),D%t)

       !call moment_tensor_body_force(D%B(2),D%t)
    end if
    end if
  end subroutine set_rates


  subroutine update_fields(D,dt,stage,RKstage)

    !> @brief use rates to update fields

    use block, only : update_fields_block
    use iface, only : update_fields_iface
    use plastic, only : update_fields_plastic
    implicit none

    type(domain_type),intent(inout) :: D
    real(kind = wp),intent(in) :: dt
    integer,intent(in) :: stage,RKstage

    integer :: i

    ! first within blocks and on their boundaries

      if (in_block_comm(1)) call update_fields_block(D%B(1),dt)
      if (in_block_comm(2)) call update_fields_block(D%B(2),dt)


    ! and then on interfaces

    do i = 1,D%nifaces
       call update_fields_iface(D%I(i),dt)
    end do

    if (stage==RKstage) then
       !if (in_block_comm(1)) call update_fields_plastic(D%B(1),D%B(1)%P,D%B(1)%G,D%B(1)%M,D%dt,D%t,D%problem,D%response)
       !if (in_block_comm(2)) call update_fields_plastic(D%B(2),D%B(2)%P,D%B(2)%G,D%B(2)%M,D%dt,D%t,D%problem,D%response)
       if (in_block_comm(1)) call update_fields_plastic(D%B(1),D%B(1)%P,D%B(1)%G,D%B(1)%M,D%dt,D%t,D%problem,D%response,&
            D%plastic_model)
       if (in_block_comm(2)) call update_fields_plastic(D%B(2),D%B(2)%P,D%B(2)%G,D%B(2)%M,D%dt,D%t,D%problem,D%response,&
            D%plastic_model)
    end if
  end subroutine update_fields

end module domain
