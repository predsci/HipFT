!#######################################################################
!
!      _    _ _       ______ _______
!     | |  | (_)     |  ____|__   __|
!     | |__| |_ _ __ | |__     | |
!     |  __  | | '_ \|  __|    | |
!     | |  | | | |_) | |       | |
!     |_|  |_|_| .__/|_|       |_|
!              | |
!              |_|
!
! ****** HipFT: High Performance Flux Transport.
!
!     Authors:  Ronald M. Caplan
!               Miko M. Stulajter
!
!     HipFT incorporates code from multiple tools developed by
!     Zoran Mikic and Jon A. Linker.
!
!     Predictive Science Inc.
!     www.predsci.com
!     San Diego, California, USA 92121
!
!#######################################################################
! Copyright 2022 Predictive Science Inc.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!#######################################################################
!
!#######################################################################
module ident
!
!-----------------------------------------------------------------------
! ****** Set the name, version, and date of code.
!-----------------------------------------------------------------------
!
      character(*), parameter :: cname='HipFT'
      character(*), parameter :: cvers='0.15.0'
      character(*), parameter :: cdate='11/22/2022'
!
end module
!#######################################################################
module number_types
!
!-----------------------------------------------------------------------
! ****** Set precisions for REALs.
!-----------------------------------------------------------------------
!
      use iso_fortran_env
!
! ****** Use double precision.
!
      integer, parameter :: r_typ = REAL64
!
end module
!#######################################################################
module constants
!
!-----------------------------------------------------------------------
! ****** Constants in r_typ precision for use throughout the code.
! ****** Used for simplicity and readability.
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero = 0.0_r_typ
      real(r_typ), parameter :: one = 1.0_r_typ
      integer*8,   parameter :: one_int = 1
      real(r_typ), parameter :: two = 2.0_r_typ
      integer*8,   parameter :: two_int = 2
      real(r_typ), parameter :: three = 3._r_typ
      integer*8,   parameter :: three_int = 3
      real(r_typ), parameter :: four = 4._r_typ
      integer*8,   parameter :: four_int = 4
      real(r_typ), parameter :: six = 6._r_typ
      real(r_typ), parameter :: nine = 9._r_typ
      real(r_typ), parameter :: ten = 10._r_typ
      real(r_typ), parameter :: fifteen = 15._r_typ
      real(r_typ), parameter :: sixteen = 16._r_typ
      real(r_typ), parameter :: half = 0.5_r_typ
      real(r_typ), parameter :: quarter = 0.25_r_typ
      real(r_typ), parameter :: twentyfour = 24.0_r_typ
      real(r_typ), parameter :: twentyfive = 25.0_r_typ
      real(r_typ), parameter :: fivehundred_i = 0.002_r_typ
      real(r_typ), parameter :: three_quarter = 0.75_r_typ
      real(r_typ), parameter :: two_third = 0.66666666666666666_r_typ
      real(r_typ), parameter :: third = 0.33333333333333333_r_typ
!
      real(r_typ), parameter :: pi = 3.1415926535897932_r_typ
      real(r_typ), parameter :: pi_two = 1.5707963267948966_r_typ
      real(r_typ), parameter :: pi_i = 0.3183098861837907_r_typ
      real(r_typ), parameter :: twopi = 6.2831853071795864_r_typ
      real(r_typ), parameter :: twopi_i = 0.15915494309189535_r_typ
      real(r_typ), parameter :: threepi_two = 4.71238898038469_r_typ
      real(r_typ), parameter :: threepi_four = 2.356194490192345_r_typ
!
      real(r_typ), parameter :: d2r = 0.017453292519943295_r_typ
      real(r_typ), parameter :: r2d = 57.29577951308232_r_typ
!
      real(r_typ), parameter :: rsun_cm = 6.96e10_r_typ
      real(r_typ), parameter :: rsun_cm2 = 4.84416e21_r_typ
!
      real(r_typ), parameter :: &
                         diff_km2_s_to_rs2_s = 2.0643413925221e-12_r_typ
      real(r_typ), parameter :: &
                         diff_km2_s_to_rs2_hr = 7.43162901307967e-09_r_typ
      real(r_typ), parameter :: &
                         km_s_to_rs_hr = 0.005172413793103448_r_typ
      real(r_typ), parameter :: &
                         m_s_to_rs_hr = 5.172413793103448e-06_r_typ
      real(r_typ), parameter :: &
                         km_s_to_rs_s = 1.4367816091954023e-06_r_typ
      real(r_typ), parameter :: output_flux_fac = 1.0e-21_r_typ
!
      real(r_typ), parameter :: small_value = tiny(one)
      real(r_typ), parameter :: large_value = huge(one)
      real(r_typ), parameter :: safety = 0.95_r_typ
!
end module
!#######################################################################
module input_parameters
!
!-----------------------------------------------------------------------
! ****** Input parameters.
!-----------------------------------------------------------------------
!
      use number_types
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical        :: verbose = .false.
!
      character(512) :: initial_map_filename = 'br_input.h5'
      character(512) :: output_map_root_filename = 'hipft_brmap'
      character(512) :: output_map_directory = 'output_maps'
!
! ****** Number of realizations ********
!
      integer :: n_realizations = 1
!
! ****** Validation Mode ********
!
      integer        :: validation_run = 0
!
! ****** Output map cadence and type ********
!
      integer        :: output_map_idx_cadence = 0
      real(r_typ)    :: output_map_time_cadence = zero
      logical        :: output_map_2d = .true.
!
! ****** Restarts ********
!
      logical        :: restart_run = .false.
      character(512) :: restart_file = ' '
!
! ****** Time ********
!
      real(r_typ)    :: time_start = zero
      real(r_typ)    :: time_end = one
      real(r_typ)    :: time_output_cadence = zero
      real(r_typ)    :: time_restart_cadence = zero
!
! ****** Timestep ********
!
      real(r_typ)    :: dt_max_increase_fac = zero
      real(r_typ)    :: dt_min = 1.0e-15_r_typ
      real(r_typ)    :: dt_max = huge(one)
!
! ****** Algorithm options.
!
      logical :: strang_splitting = .true.
!
! ****** Analysis options.
!
      real(r_typ) :: pole_flux_lat_limit = 30.0_r_typ
!
!-----------------------------------------------------------------------
!
! ****** FLOWS ********
!
! ****** Activate the flow advance.
!
      logical :: advance_flow = .false.
!
! ****** Add a rigid rotation vp velocity (km/s) of omega*sin(theta).
!
      real(r_typ)    :: flow_vp_rigid_omega = zero
!
! ****** Add a rigid rotation velocity (km/s) of omega for a rigid
! ****** rotation through the poles.
!
      real(r_typ)    :: flow_rigid_omega = zero
!
! ****** Add a constant vt velocity (km/s).
!
      real(r_typ)    :: flow_vt_const = zero
!
! ****** Attenuate the veolcity based on the value of Br.
! ****** This causes flow to be updated each step.
!
      logical    :: flow_attenuate = .false.
!
! ****** Built-in differential roation and meridianal flow models.
! ****** For each, setting "1" sets the model/params used in the AFT code.
!
      integer :: flow_dr_model = 0
      integer :: flow_mf_model = 0
!
! ****** Time dependent flows from file (used for conflow).
!
      logical :: use_flow_from_files = .false.
      character(512) :: flow_list_filename = ' '
      character(512) :: flow_root_dir = '.'
!
! ****** Algorithm options.
!        Can set upwind to central differencing by also setting UPWIND=0.
! ****** 1: FE+Upwind.
! ****** 2: RK3TVD+Upwind.
! ****** 3: RK3TVD+WENO3.
!
      integer :: flow_num_method = 3
!
! ****** Upwind coefficient.
!
      real(r_typ) :: upwind = one
!
!-----------------------------------------------------------------------
!
! ****** DIFFUSION ********
!
      logical :: advance_diffusion = .false.
!
      real(r_typ)    :: diffusion_coef_constant = zero
      character(512) :: diffusion_coef_filename = ' '
      logical        :: diffusion_coef_grid = .false.
      real(r_typ)    :: diffusion_coef_factor = diff_km2_s_to_rs2_hr
!
! ****** Algorithm options.
!
! ****** Select diffusion algorithm.
! ******   1: Explicit Euler (1st-order)
! ******   2: Explicit RKL2 Super Time-stepping (2nd-order)
! ******   3: Explicit RKG2 Super Time-stepping (2nd-order)
!
      integer :: diffusion_num_method = 3
!
! ****** Set number of diffusion subcycles per flow step.
! ****** For diffusion-only runs with RKG2(RKL2), set to ~30(60).
! ****** For flow+diffusion runs, this usually can be ~1.
!
      integer :: diffusion_subcycles = 30
!
!-----------------------------------------------------------------------
!
! ****** SOURCES ********
!
      logical :: advance_source = .false.
      character(512) :: source_filename = ' '
!
!-----------------------------------------------------------------------
!
! ****** DATA ASSIMILATION ********
!
      logical :: assimilate_data = .false.
      character(512) :: assimilate_data_map_list_filename = ' '
      character(512) :: assimilate_data_map_root_dir = '.'
!
      integer, parameter :: IO_DATA_IN = 8
!
end module
!#######################################################################
module output
!
!-----------------------------------------------------------------------
! ****** Output.
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical :: output_current_map = .false.
!
! ****** File sequence number.
!
      integer*8 :: idx_out = 0
!
      integer, parameter :: IO_HIST_NUM = 20
      integer, parameter :: IO_HIST_SOL = 21
      integer, parameter :: IO_MAP_OUT_LIST = 22
!
      character(15) :: io_hist_num_filename = 'history_num.dat'
      character(15) :: io_hist_sol_filename = 'history_sol.dat'
      character(19) :: io_map_output_list_filename = 'map_output_list.txt'
!
      real(r_typ), dimension(:,:,:), allocatable :: fout
      real(r_typ), dimension(:), allocatable :: pout
!
end module
!#######################################################################
module globals
!
!-----------------------------------------------------------------------
! ****** Internal global variables.
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Current time.
!
      real(r_typ) :: time = 0.
!
! ****** Current step.
!
      integer*8 :: ntime = 0
!
! ****** Current time-step.
!
      real(r_typ) :: dtime_global = 0.
      character(50) :: dtime_reason = ''
!
! ****** Explicit Euler diffusion stable time-step and number of cycles.
!
      real(r_typ) :: dtime_diffusion_euler = 0.
      integer*8 :: n_stable_diffusion_cycles = 1
      real(r_typ) :: dtime_diffusion_used = 0.
!
! ****** Explicit Euler advection stable time-step.
!
      real(r_typ) :: dtmax_flow = 0.
      real(r_typ) :: dtime_advection_used = 0.
!
! ****** Current output file sequence number.
!
      integer*8 :: output_seq = 0
!
! ****** Flag to indicate the flow needs updating.
!
      logical :: flow_needs_updating = .true.
!
! ****** Flag to indicate the time step needs updating.
!
      logical :: timestep_needs_updating = .true.
!
! ****** Flag to indicate the stable flow time step needs updating.
!
      logical :: timestep_flow_needs_updating = .true.
!
! ****** Flag to indicate the diffusion time step needs updating.
!
      logical :: timestep_diff_needs_updating = .true.
!
! ****** History analysis variables.
!
      real(r_typ) :: h_minbr, h_maxbr, h_minabsbr,                   &
                     h_fluxp, h_fluxm, h_valerr,                     &
                     h_fluxp_pn, h_fluxm_pn, h_fluxp_ps, h_fluxm_ps, &
                     h_area_pn, h_area_ps, h_eq_dipole, h_ax_dipole
!
      real(r_typ) :: u0max,u0min
!
end module
!#######################################################################
module mesh
!
      use number_types
!
      implicit none
!
      integer :: nt,ntm,ntm1,ntm2
      integer :: np,npm,npm1,npm2
      integer :: nr
!
      real(r_typ), dimension(:), allocatable :: t,p
      real(r_typ), dimension(:), allocatable :: th,ph
      real(r_typ), dimension(:), allocatable :: dt,dth,dp,dph
      real(r_typ), dimension(:), allocatable :: dt_i,dth_i,dp_i,dph_i
      real(r_typ), dimension(:), allocatable :: st,sth,ct,cth
      real(r_typ), dimension(:), allocatable :: st_i,sth_i
!
end module
!#######################################################################
module fields
!
      use number_types
!
      implicit none
!
      real(r_typ), dimension(:,:,:), allocatable :: f
      real(r_typ), dimension(:,:,:), allocatable :: fold
      real(r_typ), dimension(:,:,:), allocatable :: fval_u0
      real(r_typ), dimension(:,:,:), allocatable :: diffusion_coef
      real(r_typ), dimension(:,:,:), allocatable :: source
      real(r_typ), dimension(:,:,:), allocatable :: vt
      real(r_typ), dimension(:,:,:), allocatable :: vp
!
end module
!#######################################################################
module data_assimilation
!
      use number_types
!
      implicit none
!
      integer :: num_maps_in_list = 0
      integer :: current_map_input_idx = 1
      real(r_typ) :: time_of_next_input_map = 0.
!
      real(r_typ) :: map_time_initial_hr
!
      real(r_typ), dimension(:), allocatable :: map_times_actual_ut_jd
!
      character(19), dimension(:), allocatable :: &
                     map_times_requested_ut_str, map_times_actual_ut_str
!
      character(512), dimension(:), allocatable :: map_files_rel_path
!
end module
!#######################################################################
module flow_from_files
!
      use number_types
!
      implicit none
!
      integer :: num_flows_in_list = 0
      integer :: current_flow_input_idx = 1
      real(r_typ) :: time_of_next_input_flow = 0.
!
      real(r_typ) :: flow_time_initial_hr
!
      real(r_typ), dimension(:), allocatable :: flow_times_actual_ut_jd
      real(r_typ), dimension(:), allocatable :: flow_times_actual_ut_jd0
!
      character(512), dimension(:), allocatable :: flow_files_rel_path_t
      character(512), dimension(:), allocatable :: flow_files_rel_path_p
!
      real(r_typ), dimension(:,:,:),allocatable :: flow_from_file_vt_old
      real(r_typ), dimension(:,:,:),allocatable :: flow_from_file_vp_old
!
end module
!#######################################################################
module sts
!
      use number_types
!
      implicit none
!
      integer*8 :: sts_s
!
      logical :: need_to_load_sts = .true.
!
      real(r_typ), dimension(:), allocatable :: sts_uj
      real(r_typ), dimension(:), allocatable :: sts_vj
      real(r_typ), dimension(:), allocatable :: sts_ubj
      real(r_typ), dimension(:), allocatable :: sts_gj
      real(r_typ), dimension(:), allocatable :: sts_b
!
      real(r_typ), dimension(:,:,:), allocatable :: u0
      real(r_typ), dimension(:,:,:), allocatable :: dty0
      real(r_typ), dimension(:,:,:), allocatable :: ykm1
      real(r_typ), dimension(:,:,:), allocatable :: ukm1
      real(r_typ), dimension(:,:,:), allocatable :: ukm2
!
end module
!#######################################################################
module weno
!
      use number_types
!
      implicit none
!
      real(r_typ), dimension(:), allocatable :: D_C_CPt
      real(r_typ), dimension(:), allocatable :: D_C_MCt
      real(r_typ), dimension(:), allocatable :: D_M_Tt
      real(r_typ), dimension(:), allocatable :: D_CP_Tt
      real(r_typ), dimension(:), allocatable :: D_P_Tt
      real(r_typ), dimension(:), allocatable :: D_MC_Tt
!
      real(r_typ), dimension(:), allocatable :: D_C_CPp
      real(r_typ), dimension(:), allocatable :: D_C_MCp
      real(r_typ), dimension(:), allocatable :: D_M_Tp
      real(r_typ), dimension(:), allocatable :: D_CP_Tp
      real(r_typ), dimension(:), allocatable :: D_P_Tp
      real(r_typ), dimension(:), allocatable :: D_MC_Tp
!
      real(r_typ), dimension(:,:,:), allocatable :: alpha_t
      real(r_typ), dimension(:,:,:), allocatable :: alpha_p
!
      real(r_typ), parameter :: weno_eps = 1.0e-6_r_typ
!
      real(r_typ) :: p0p,p1p,p0m,p1m
      real(r_typ) :: B0p,B1p,B0m,B1m
      real(r_typ) :: w0p,w1p,w0m,w1m
      real(r_typ) :: wp_sum,wm_sum
      real(r_typ) :: OM0p,OM1p,OM0m,OM1m
      real(r_typ) :: up,um
!
end module
!#######################################################################
module matrix_storage
!
      use number_types
!
      implicit none
!
      real(r_typ), dimension(:,:,:,:), allocatable :: coef
!
end module
!#######################################################################
module timing
!
      use number_types
!
      implicit none
!
      real(r_typ) :: wtime_tmp = 0.
      real(r_typ) :: wtime_tmp_mpi = 0.
!
      real(r_typ) :: wtime_setup = 0.
      real(r_typ) :: wtime_update = 0.
      real(r_typ) :: wtime_flux_transport = 0.
      real(r_typ) :: wtime_flux_transport_advection = 0.
      real(r_typ) :: wtime_flux_transport_diffusion = 0.
      real(r_typ) :: wtime_source = 0.
      real(r_typ) :: wtime_analysis = 0.
      real(r_typ) :: wtime_io = 0.
      real(r_typ) :: wtime_mpi_overhead = 0.
!
      real(r_typ) :: wtime_total = 0.
!
end module
!#######################################################################
module rp1d_def
!
!-----------------------------------------------------------------------
! ****** Define a structure to hold a REAL 1D pointer.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      type :: rp1d
        real(r_typ), dimension(:), pointer, contiguous :: f
      end type
!
end module
!#######################################################################
module ds_def
!
!-----------------------------------------------------------------------
! ****** Definition of the IO data structure.
!-----------------------------------------------------------------------
!
      use number_types
      use rp1d_def
!
      implicit none
!
      integer, parameter, private :: mxdim = 3
!
      type :: ds
        integer :: ndim
        integer, dimension(mxdim) :: dims
        logical :: scale
        logical :: hdf32
        type(rp1d), dimension(mxdim) :: scales
        real(r_typ), dimension(:,:,:), pointer, contiguous :: f
      end type
!
end module
!#######################################################################
module read_2d_file_interface
      interface
        subroutine read_2d_file (fname,ln1,ln2,fin,s1,s2,ierr)
        use number_types
        use ds_def
        use timing
        implicit none
        character(*) :: fname
        real(r_typ), dimension(:,:), allocatable :: fin
        real(r_typ), dimension(:), allocatable :: s1,s2
        integer :: ln1,ln2,ierr
        end subroutine
      end interface
end module
!#######################################################################
module read_3d_file_interface
      interface
        subroutine read_3d_file (fname,ln1,ln2,ln3,fin,s1,s2,s3,ierr)
        use number_types
        use ds_def
        use timing
        implicit none
        character(*) :: fname
        real(r_typ), dimension(:,:,:), allocatable :: fin
        real(r_typ), dimension(:), allocatable :: s1,s2,s3
        integer :: ln1,ln2,ln3,ierr
        end subroutine
      end interface
end module
!#######################################################################
module assimilate_new_data_interface
      interface
        subroutine assimilate_new_data (new_data)
        use number_types
        use input_parameters
        use globals
        use constants
        use fields
        use mesh
        use data_assimilation
        implicit none
        real(r_typ), dimension(:,:,:,:), allocatable :: new_data
        end subroutine
      end interface
end module
!#######################################################################
module mpidefs
!
!-----------------------------------------------------------------------
! ****** MPI variables, processor topology, and processor information.
!-----------------------------------------------------------------------
!
      use number_types
      use mpi
!
      implicit none
!
! ****** Array of nr per rank.
!
      integer, dimension(:), allocatable :: nr_arr
!
! ****** Total number of processors.
!
      integer :: nproc
!
! ****** Total number of processors per node.
!
      integer :: nprocsh
!
! ****** Processor rank of this process in communicator
! ****** MPI_COMM_WORLD.
!
      integer :: iprocw
!
! ****** Processor rank of this process in communicator
! ****** comm_shared.
!
      integer :: iprocsh
!
! ****** Flag to designate that this is the processor with
! ****** rank 0 in communicator MPI_COMM_WORLD.
!
      logical :: iamp0
!
! ****** Processor rank of this process in communicator
! ****** comm_shared.
!
      integer :: iproc
!
! ****** Processor rank in communicator comm_shared for the
! ****** processor that has rank 0 in MPI_COMM_WORLD.
!
      integer :: iproc0
!
! ****** Communicator over all shared processors on the node.
!
      integer :: comm_shared
!
! ****** Number type for REALs to be used in MPI calls.
!
      integer :: ntype_real
!
! ****** GPU device number for current rank.
!
      integer :: igpu
!
! ****** Broadcast sizes
!
      integer :: flow_t_size,flow_p_size
!
      end module
!#######################################################################
program HIPFT
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use globals
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical :: the_run_is_done
      character(*), parameter :: &
                           FMT = '(a,i10,a,f12.6,a,f12.6,a,f12.6,a,a,a)'
!
!-----------------------------------------------------------------------
!
! ****** Initialize MPI.
!
      call initialize_mpi
!
! ****** Set up and initialize the run.
!
      call setup
!
! ****** START MAIN LOOP ********************************************
!
      if (iamp0) then
        write (*,*)
        write (*,*) '>running'
        write (*,*)
        flush(OUTPUT_UNIT)
      end if
!
      do
!
! ****** Update step number we are about to start. (ntime=0 initially)
!
        ntime = ntime + 1
!
! ****** Update prescribed quantities including data assimilation and/or
! ****** flows and set/check the timestep.
!
        call update_step
!
        if (verbose .and. iamp0) then
          write(*,FMT) 'Starting step: ',ntime,'  Time:',time,' / ', &
                        time_end,'  dt to use:',dtime_global, &
                        ' (',TRIM(dtime_reason),')'
          flush(OUTPUT_UNIT)
        end if
!
! ****** Evolve the field for one timestep.
!
        call flux_transport_step
!
! ****** Update the time now that we have taken the step.
!
        time = time + dtime_global
!
! ****** Perform analysis.
!
        call analysis_step
!
! ****** Output.
!
        call output_step
!
        if (iamp0) then
          write(*,FMT) 'Completed step: ',ntime,'  Time: ',time,' / ', &
                        time_end,'  dt: ',dtime_global, &
                        ' (',TRIM(dtime_reason),')'
          flush(OUTPUT_UNIT)
        end if
!
! ****** Check if the run is done.
!
        if (the_run_is_done()) exit
!
      enddo
!
! ****** END MAIN LOOP **********************************************
!
! ****** Perform end-of-run tasks.
!
      call wrap_it_up
!
! ****** Get wall-clock time.
!
      wtime_total = MPI_Wtime() - wtime_total
!
! ****** Write out time profiling.
!
      call write_timing
!
! ****** End the run.
!
      call endrun (.false.)
!
end program HIPFT
!#######################################################################
subroutine endrun (ifstop)
!
!-----------------------------------------------------------------------
!
! ****** End the run and exit the code.
!
!-----------------------------------------------------------------------
!
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical :: ifstop
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr = 0
!
!-----------------------------------------------------------------------
!
! ****** Exit MPI gracefully.
!
      call MPI_Finalize (ierr)
!
! ****** Call the STOP statement if requested.
!
      if (ifstop) stop
!
end subroutine
!#######################################################################
subroutine set_local_number_of_realizations
!
!
!-----------------------------------------------------------------------
! ****** Set number of realizations per MPI rank
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use input_parameters
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      nr = FLOOR(REAL(n_realizations,r_typ)/REAL(nproc,r_typ))
!
      if (nr*nproc .lt. n_realizations) then
        write(*,*)
        write(*,*) 'Warning!! Load imbalance detected.'
        write(*,*)
        if (iamp0) then
          nr = nr + 1
        elseif (mod(n_realizations,iprocsh) .gt. iprocsh) then
          nr = nr + 1
        end if
      end if
!
      if (iamp0) then
        allocate (nr_arr(0:nproc-1))
      end if
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Gather(nr,1,MPI_INT,nr_arr,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
      if (iamp0) then
!$acc enter data create(nr_arr)
      end if
!
end subroutine
!#######################################################################
subroutine initialize_mpi
!
!-----------------------------------------------------------------------
!
! ****** Initialize MPI.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr,tcheck
!
!-----------------------------------------------------------------------
!
      call MPI_Init_thread (MPI_THREAD_FUNNELED,tcheck,ierr)
!
! ****** Start wall clock timer.
!
      wtime_total = MPI_Wtime()
      wtime_tmp_mpi = wtime_total
!
! ****** Get the total number of processors.
!
      call MPI_Comm_size (MPI_COMM_WORLD,nproc,ierr)
!
! ****** Get the index (rank) of the local processor in
! ****** communicator MPI_COMM_WORLD in variable IPROCW.
!
      call MPI_Comm_rank (MPI_COMM_WORLD,iprocw,ierr)
!
! ****** Create a shared communicator for all ranks in the node.
!
      call MPI_Comm_split_type (MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED, &
                                0,MPI_INFO_NULL,comm_shared,ierr)
!
! ****** Get the total number of processors in node.
!
      call MPI_Comm_size (comm_shared,nprocsh,ierr)
!
! ****** Get the index (rank) of the local processor in the local node.
!
      call MPI_Comm_rank (comm_shared,iprocsh,ierr)
!
! ****** Set the flag to designate whether this processor
! ****** has rank 0 in communicator MPI_COMM_WORLD.
!
      if (iprocw .eq. 0) then
        iamp0 = .true.
      else
        iamp0 = .false.
      end if
!
! ****** Set the number type for communicating REAL
! ****** numbers in MPI calls.
!
      ntype_real = MPI_REAL8
!
! ****** Set the GPU device number based on the shared rank.
! ****** This assumes the code was launched with 1 MPI rank per GPU.
!
!$acc set device_num(iprocsh)
!
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
end subroutine
!#######################################################################
subroutine setup
!
!-----------------------------------------------------------------------
!
! ****** Setup the code run.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use input_parameters
      use timing
      use fields
      use mesh
      use output
      use globals, ONLY : time
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(512) :: cmd
!
!-----------------------------------------------------------------------
!
      wtime_tmp = MPI_Wtime()
!
      call read_input_file
!
! ****** Set up output directory.
!
      if (output_map_idx_cadence.gt.0 .or.     &
          output_map_time_cadence.gt.0.0) then
        if (iamp0) then
          cmd = 'mkdir -p '//trim(output_map_directory)
          call EXECUTE_COMMAND_LINE (cmd,exitstat=ierr)
          if (ierr.ne.0) then
            write (*,*) '### Could not make output subdirectory, using ./'
            output_map_directory = '.'
          end if
        end if
      end if
!
      if (restart_run) then
!        call load_restart
         print*,'RESTARTS ARE NOT IMPLEMENTED YET.'
         stop
      else
        call load_initial_condition
        time = time_start
      end if
!
! ****** Allocate flow arrays.
!
      if (advance_flow) then
        allocate (vt(nt,npm,nr))
        allocate (vp(ntm,np,nr))
        vt(:,:,:) = 0.
        vp(:,:,:) = 0.
!$acc enter data copyin(vt,vp)
      end if
!
      call set_unchanging_quantities
!
      if (assimilate_data) call load_data_assimilation
!
      if (use_flow_from_files) call load_flow_from_files
!      RMC:  For now, only have rank 0 do this.  Eventually, all ranks need to
!            create files for each realization they have...
      if (iamp0) call create_and_open_output_log_files
!
      call analysis_step
!
      if (output_map_idx_cadence .gt. 0 .or. &
          output_map_time_cadence .gt. 0.) then
        output_current_map = .true.
      end if
!
      call output_step
!
      if (iamp0) call write_welcome_message
!
      wtime_setup = MPI_Wtime() - wtime_tmp
!
end subroutine
!#######################################################################
function the_run_is_done ()
!
!-----------------------------------------------------------------------
!
! ****** Check if the run is done.
!
!-----------------------------------------------------------------------
!
      use globals, ONLY : time
      use input_parameters, ONLY : time_end
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical :: the_run_is_done
!
!-----------------------------------------------------------------------
!
      the_run_is_done = .false.
!
! ****** Check if the simulation time has reached the end time.
!
      if (time .ge. time_end) the_run_is_done = .true.
!
! ****** Check if a STOPRUN file has been created.
!
!        check_stoprun
!
! ****** Check if the maximum wall clock time has been reached.
!
!        check_wallclock
!
      return
end function
!#######################################################################
subroutine wrap_it_up
!
!-----------------------------------------------------------------------
!
! ****** Perform end-of-run tasks.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use timing
      use input_parameters
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Write the final Br map.
!
      call write_map (trim(output_map_root_filename)//'_final.h5')
!
      if (iamp0) then
        write(*,*) ' '
        write(*,*) 'Final 2D data written out to ', &
                     trim(output_map_root_filename)//'_final.h5'
!
!     output_final_restart?
!
        write(*,*)
        write(*,*) "Run completed at:"
        call timestamp
        write(*,*)
        flush(OUTPUT_UNIT)
      end if
!
!****** Synchronize all MPI ranks to get accurate wall clock time.
!
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
end subroutine
!#######################################################################
subroutine write_timing
!
!-----------------------------------------------------------------------
!
! ****** Write out timing.
!
!-----------------------------------------------------------------------
!
      use timing
      use number_types
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Timing buffers.
!
      integer, parameter :: lbuf=10
      real(r_typ), dimension(lbuf) :: sbuf
      real(r_typ), dimension(lbuf,0:nproc-1) :: tbuf
!
! ****** Timing statistics.
!
      real(r_typ), dimension(lbuf) :: tmin,tmax,tavg,tsdev
!
!-----------------------------------------------------------------------
!
      integer :: ierr,irank
      character(*), parameter :: FMT = '(a20,f20.2)'
!
!-----------------------------------------------------------------------
!
!
! ****** Gather the timing information for all processors into TBUF.
!
      sbuf(1) = wtime_total
      sbuf(2) = wtime_setup
      sbuf(3) = wtime_update
      sbuf(4) = wtime_flux_transport
      sbuf(5) = wtime_flux_transport_advection
      sbuf(6) = wtime_flux_transport_diffusion
      sbuf(7) = wtime_source
      sbuf(8) = wtime_analysis
      sbuf(9) = wtime_io
      sbuf(10) = wtime_mpi_overhead
!
      call MPI_Allgather (sbuf,lbuf,ntype_real,       &
                          tbuf,lbuf,ntype_real,MPI_COMM_WORLD,ierr)
!
! ****** Calculate the timing statistics.
!
      tavg=sum(tbuf,dim=2)/nproc
      tmin=minval(tbuf,dim=2)
      tmax=maxval(tbuf,dim=2)
!
      tsdev(:)=0.
!
      do irank=0,nproc-1
        tsdev(:)=tsdev(:)+(tbuf(:,irank)-tavg(:))**2
      enddo
!
      tsdev(:)=sqrt(tsdev(:)/nproc)
!
      if (iamp0) then
!
        do irank=0,nproc-1
          write(*,"(a40)") repeat("-", 40)
          write(*,*)
          write (*,*) 'Processor id = ',irank
          write(*,*)
          write(*,"(a40)") repeat("-", 40)
          write(*,FMT) "Wall clock time:   ",tbuf(1,irank)
          write(*,"(a40)") repeat("-", 40)
          write(*,FMT) "--> Setup:         ",tbuf(2,irank)
          write(*,FMT) "--> Update:        ",tbuf(3,irank)
          write(*,FMT) "--> Flux transport:",tbuf(4,irank)
          write(*,FMT) "    --> Advecton:  ",tbuf(5,irank)
          write(*,FMT) "    --> Diffusion  ",tbuf(6,irank)
          write(*,FMT) "    --> Source:    ",tbuf(7,irank)
          write(*,FMT) "--> Analysis:      ",tbuf(8,irank)
          write(*,FMT) "--> I/O:           ",tbuf(9,irank)
          write(*,FMT) "--> MPI overhead:  ",tbuf(10,irank)
          write(*,"(a40)") repeat("-", 40)
          write(*,*)
        enddo
!
        write (*,*)
        write (*,"(a40)") repeat("-", 40)
        write (*,*)
        write (*,*) 'Average times:'
        write (*,*) '-------------'
        write (*,*)
        write (*,300) 'Avg         Min         Max      S. Dev'
        write (*,300) '---         ---         ---      ------'
        write (*,400) 'Wall clock time:        ', &
                     tavg(1),tmin(1),tmax(1),tsdev(1)
        write (*,400) '--> Setup:              ', &
                     tavg(2),tmin(2),tmax(2),tsdev(2)
        write (*,400) '--> Update:             ', &
                     tavg(3),tmin(3),tmax(3),tsdev(3)
        write (*,400) '--> Flux transport:     ', &
                     tavg(4),tmin(4),tmax(4),tsdev(4)
        write (*,400) '    --> Advecton:       ', &
                     tavg(5),tmin(5),tmax(5),tsdev(5)
        write (*,400) '    --> Diffusion       ', &
                     tavg(6),tmin(6),tmax(6),tsdev(6)
        write (*,400) '    --> Source:         ', &
                     tavg(7),tmin(7),tmax(7),tsdev(7)
        write (*,400) '--> Analysis:           ', &
                     tavg(8),tmin(8),tmax(8),tsdev(8)
        write (*,400) '--> I/O:                ', &
                    tavg(9),tmin(9),tmax(9),tsdev(9)
        write (*,400) '--> MPI overhead:       ', &
                    tavg(10),tmin(10),tmax(10),tsdev(10)
!
  300   format (1x,33x,a)
  400   format (1x,a,4f12.3)
!
        write (*,*)
        write(*,"(a40)") repeat("-", 40)
        write (*,*)
        flush(OUTPUT_UNIT)
      end if
!
end subroutine
!#######################################################################
subroutine create_and_open_output_log_files
!
!-----------------------------------------------------------------------
!
! ****** Create and open output log files
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use output
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      call FFOPEN (IO_HIST_NUM,io_hist_num_filename,'rw',ierr)
!
      write (IO_HIST_NUM,'(a10,a,6(a22,a),a15)') &
      'STEP',' ',&
      'TIME',' ',&
      'DTIME',' ',&
      'DTIME_ADV_STB',' ',&
      'DTIME_ADV_USED',' ',&
      'DTIME_DIFF_STB',' ',&
      'DTIME_DIFF_USED',' ',&
      'N_DIFF_PER_STEP'
!
      close(IO_HIST_NUM)
!
      call FFOPEN (IO_HIST_SOL,io_hist_sol_filename,'rw',ierr)
!
      write (IO_HIST_SOL,'(a10,a,15(a22,a))') &
      'STEP',' ',&
      'TIME',' ',&
      'BR_MIN',' ',&
      'BR_MAX',' ',&
      'BR_ABS_MIN', ' ',&
      'FLUX_POSITIVE',' ',&
      'FLUX_NEGATIVE',' ',&
      'NPOLE_FLUX_POSITIVE',' ',&
      'NPOLE_FLUX_NEGATIVE',' ',&
      'SPOLE_FLUX_POSITIVE',' ',&
      'SPOLE_FLUX_NEGATIVE',' ',&
      'NPOLE_AREA',' ',&
      'SPOLE_AREA',' ',&
      'EQ_DIPOLE',' ',&
      'AX_DIPOLE',' ',&
      'VALIDATION_ERR_CVRMSD'
!
      close(IO_MAP_OUT_LIST)
!
      if (output_map_idx_cadence.gt.0 .or.     &
          output_map_time_cadence.gt.0.0) then
!
        call FFOPEN (IO_MAP_OUT_LIST,io_map_output_list_filename,&
                     'rw',ierr)
!
        write (IO_MAP_OUT_LIST,'(a10,a,a25,a,a)') &
                               'IDX',' ',&
                               'TIME',' ',&
                               'FILENAME'
!
        close(IO_MAP_OUT_LIST)
      end if
!
end subroutine
!#######################################################################
subroutine write_welcome_message
!
!-----------------------------------------------------------------------
!
! ****** Write welcome message
!
!-----------------------------------------------------------------------
!
      use ident
      use input_parameters
      use mesh
      use globals
      use fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      write (*,*) ''
      write (*,*) ''
      write (*,*) '      _    _ _       ______ _______ '
      write (*,*) '     | |  | (_)     |  ____|__   __|'
      write (*,*) '     | |__| |_ _ __ | |__     | |'
      write (*,*) '     |  __  | | ''_ \|  __|    | |'
      write (*,*) '     | |  | | | |_) | |       | |'
      write (*,*) '     |_|  |_|_| .__/|_|       |_|'
      write (*,*) '              | |'
      write (*,*) '              |_|'
      write (*,*) ''
      write (*,*) '     Version: ',cvers,' of ',cdate
      write (*,*) ''
      write (*,*) ' ****** HipFT: High Performance Flux Transport.'
      write (*,*) ''
      write (*,*) '        Authors:  Ronald M. Caplan'
      write (*,*) '                  Miko M. Stulajter'
      write (*,*) '                  Jon A. Linker'
      write (*,*) '                  Zoran Mikic'
      write (*,*) ''
      write (*,*) '        Predictive Science Inc.'
      write (*,*) '        www.predsci.com'
      write (*,*) '        San Diego, California, USA 92121'
      write (*,*)''
      write (*,*)''
      write (*,*) 'Start time = ',time
      write (*,*) 'End time = ',time_end
      write (*,*)
      write (*,*) 'Number of t mesh points = ',ntm
      write (*,*) 'Number of p mesh points = ',npm-1
      write (*,*)
      if (advance_diffusion) then
        write (*,'(a,f12.6,a)') ' Uniform diffusion = ', &
                                       diffusion_coef_constant,' km^2/s'
        write (*,*) ' Diffusion coefficient: '
        write (*,'(a,f12.6,a)') '   Minimum value = ', &
                  MINVAL(diffusion_coef/diffusion_coef_factor),' km^2/s'
        write (*,'(a,f12.6,a)') '   Maximum value = ', &
                  MAXVAL(diffusion_coef/diffusion_coef_factor),' km^2/s'
      end if
      write(*,*)
      write(*,*) "Run started at: "
      call timestamp
      flush(OUTPUT_UNIT)
!
end subroutine
!#######################################################################
subroutine load_initial_condition
!
!-----------------------------------------------------------------------
!
! ****** load_initial_condition
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use input_parameters
      use mesh
      use fields
      use globals
      use read_2d_file_interface
      use output
      use constants
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: val2_g_width = 0.03_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: n1,n2,n3,ierr,i,j,k,lsbuf, irank
      integer, dimension(:), allocatable :: displ,lbuf
      real(r_typ), dimension(:), allocatable :: fn1,fs1,fn2,fs2
      real(r_typ), dimension(:,:,:), allocatable :: f_local,gout
      real(r_typ), dimension(:,:), allocatable :: f_tmp2d
      real(r_typ), dimension(:), allocatable :: s1,s2,s3,s3t,tfout,gfout
!
!-----------------------------------------------------------------------
!
! ****** NOTE: At this time the user cannot specify NT,NP.
! ****** Instead, the resolution of the input map is used.
! ****** The input map is assumed to be in PT and with a
! ****** one-point overlap in phi.
!
! ****** Read the initial magnetic map.
!
      if (iamp0) then
        call read_2d_file (initial_map_filename,n1,n2,f_tmp2d,s1,s2,ierr)
      endif
!
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Bcast (n1,1,ntype_real,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast (n2,1,ntype_real,0,MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
      if (.not.iamp0) then
        allocate (s1(n1))
        allocate (s2(n2))
        allocate (f_tmp2d(n1,n2))
      endif
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Bcast (f_tmp2d,n1*n2,ntype_real,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast (s1,n1,ntype_real,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast (s2,n2,ntype_real,0,MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
      call set_local_number_of_realizations
      n3 = nr
      allocate (f_local(n1,n2,n3))

      do i=1,n3
        f_local(:,:,i) = f_tmp2d(:,:)
      enddo
!
      allocate (s3(n3))
      do i=1,n3
        s3(i)=real(i,r_typ)
      enddo
!
      if (validation_run .eq. 1) then
!
        allocate (fval_u0(n1,n2,n3))
!
! ****** Make initial solution f and output.
!
! ****** Have to put this on and off GPU since the routines
!        are GPU-only.
!
!$acc enter data create (f_local,fval_u0,f_tmp2d)
!$acc enter data copyin (s1,s2)
        do i=1,n3
          call tesseral_spherical_harmonic (0,6,s1,s2,n1,n2,f_tmp2d)
          do concurrent(k=1:n2,j=1:n1)
            f_local(j,k,i) = f_tmp2d(j,k)
          enddo
          call tesseral_spherical_harmonic (5,6,s1,s2,n1,n2,f_tmp2d)
          do concurrent(k=1:n2,j=1:n1)
            fval_u0(j,k,i) = f_tmp2d(j,k)
          enddo
        enddo
!$acc update self (f_local,fval_u0)
!$acc exit data delete (f_local,fval_u0,f_tmp2d,s1,s2)
!
        fval_u0(:,:,:) = 1000.0_r_typ*(f_local(:,:,:) +                  &
                              sqrt(14.0_r_typ/11.0_r_typ)*fval_u0(:,:,:))
!
        if (output_map_2d.and.n_realizations.eq.1) then
          call write_2d_file(                                        &
               (trim(output_map_root_filename)//'_initial_0.h5')      &
                         ,n1,n2,fval_u0(:,:,1),s1,s2,ierr)
        else
!
          wtime_tmp_mpi = MPI_Wtime()
!
          lsbuf=n1*n2*n3
          allocate (tfout(lsbuf))
          tfout(:)=reshape(fval_u0(:,:,:),(/lsbuf/))
          if (iamp0) then
            allocate (gfout(n1*n2*n_realizations))
            allocate (displ(0:nproc-1))
            allocate (lbuf(0:nproc-1))
            displ(0)=0
            do irank=0,nproc-1
              lbuf(irank)=n1*n2*nr_arr(irank)
            enddo
            do irank=1,nproc-1
              displ(irank)=displ(irank-1)+lbuf(irank-1)
            enddo
          end if
!
          call MPI_Gatherv (tfout,lsbuf,ntype_real,gfout,lbuf,       &
                            displ,ntype_real,0,MPI_COMM_WORLD,ierr)
          deallocate(tfout)
          wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
          if (iamp0) then
            allocate (s3t(n_realizations))
            allocate (gout(n1,n2,n_realizations))
            do i=1,n_realizations
              s3t(i)=real(i,r_typ)
            enddo
            gout=reshape(gfout(:),(/n1,n2,n_realizations/))
            call write_3d_file(                                        &
                (trim(output_map_root_filename)//'_initial_0.h5')      &
                               ,n1,n2,n_realizations,gout,s1,s2,s3t,ierr)
            deallocate (gout)
            deallocate (gfout)
            deallocate (displ)
            deallocate (lbuf)
            deallocate (s3t)
          end if
!
        endif
!
        f_local(:,:,:) = fval_u0(:,:,:)
!
! ****** Caculate final analytic solution and output.
!
        allocate (fout(n1,n2,n3))
        fout(:,:,:) = fval_u0(:,:,:)*exp(-42.0_r_typ*diffusion_coef_constant* &
                                     diffusion_coef_factor*time_end)
!
        if (output_map_2d.and.n_realizations.eq.1) then
          call write_2d_file(                                           &
               (trim(output_map_root_filename)//'_final_analytic.h5') &
                            ,n1,n2,fout(:,:,1),s1,s2,ierr)
        else
!
          wtime_tmp_mpi = MPI_Wtime()
          lsbuf=n1*n2*n3
          allocate (tfout(lsbuf))
          tfout(:)=reshape(fout(:,:,:),(/lsbuf/))
          if (iamp0) then
            allocate (gfout(n1*n2*n_realizations))
            allocate (displ(0:nproc-1))
            allocate (lbuf(0:nproc-1))
            displ(0)=0
            do irank=0,nproc-1
              lbuf(irank)=n1*n2*nr_arr(irank)
            enddo
            do irank=1,nproc-1
              displ(irank)=displ(irank-1)+lbuf(irank-1)
            enddo
          end if
!
          call MPI_Gatherv (tfout,lsbuf,ntype_real,gfout,lbuf,       &
                    displ,ntype_real,0,MPI_COMM_WORLD,ierr)
          deallocate(tfout)
          wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
          if (iamp0) then
            allocate (s3t(n_realizations))
            allocate (gout(n1,n2,n_realizations))
            do i=1,n_realizations
              s3t(i)=real(i,r_typ)
            enddo
            gout=reshape(gfout(:),(/n1,n2,n_realizations/))
            call write_3d_file(                                        &
                (trim(output_map_root_filename)//'_final_analytic.h5')      &
                               ,n1,n2,n_realizations,gout,s1,s2,s3t,ierr)
            deallocate (gout)
            deallocate (gfout)
            deallocate (displ)
            deallocate (lbuf)
            deallocate (s3t)
          end if
!
        endif
!
        deallocate (fout)
!
!$acc enter data copyin(fval_u0)
!
      elseif (validation_run .eq. 2) then
!
! ****** Make initial solution f and output.
!
        do k=1,n3
          do j=1,n2
            do i=1,n1
              f_local(i,j,k) = - EXP( (-(s2(j)-threepi_four)**2 - &
                                        (s1(i)-pi_two      )**2 )/val2_g_width) &
                               + EXP( (-(s2(j)-threepi_four)**2 - &
                                        (s1(i)-threepi_two )**2 )/val2_g_width)
            enddo
          enddo
        enddo
!
      elseif (validation_run .eq. 3) then
!
! ****** Make initial solution f and output.
!
        do k=1,n3
          do j=1,n2
            do i=1,n1
              f_local(i,j,k) = - EXP((-(s2(j)-pi_two)**2)/(val2_g_width))
            enddo
          enddo
        enddo
!
      end if
!
! ****** Set resolution values.
!
      ntm1 = n2
      npm1 = n1
!
! ****** Allocate main mesh grid arrays and data array with extended
! ****** phi grid.
!
      allocate (t(ntm1))
      allocate (p(npm1 + 1))
      allocate (f(ntm1,npm1 + 1,nr))
!
! ****** Transpose array since we assume PT but need TP for the code.
!
      do i=1,nr
        f(1:ntm1,1:npm1,i) = TRANSPOSE(f_local(1:npm1,1:ntm1,i))
      enddo
!
! ****** Impose perfect periodicity and set two-point overlap.
!
      f(:,1,:) = half*(f(:,1,:) + f(:,npm1,:))
      f(:,npm1,:) = f(:,1,:)
      f(:,npm1 + 1,:) = f(:,2,:)
!
! ****** Set main mesh grid arrays.
!
      t(:) = s2(:)
      p(1:npm1) = s1(:)
      p(npm1+1) = p(2) + twopi
!
! ****** Set additional meshes and precomputed mesh quantities.
!
      call set_mesh
!
! ****** Get the m=0 components near the poles.
!
      allocate (fn1(nr))
      allocate (fn2(nr))
      allocate (fs1(nr))
      allocate (fs2(nr))
      call get_m0 (f,fn1,fn2,fs1,fs2)
!
! ****** Set the pole values to have only an m=0 component.
!
      do i=1,nr
        f(1,:,i) = fn1(i)
        f(ntm,:,i) = fs1(i)
      enddo
!
!$acc enter data copyin(f)
!
! ****** Clean up memory.
!
      deallocate (s1)
      deallocate (s2)
      deallocate (s3)
      deallocate (fn1,fn2,fs1,fs2)
      deallocate (f_local)
      deallocate (f_tmp2d)
!
! ****** Write out initial condition.
!
      call write_map (trim(output_map_root_filename)//'_initial.h5')
!
end subroutine
!#######################################################################
subroutine tesseral_spherical_harmonic (m, l, pvec, tvec, np, nt, Y)
!
!-----------------------------------------------------------------------
!
! Subroutine for computing tesseral spherical harmonics using the
! First Modified Forward Column Recursion method from
! Holmes et.al. (2002) which reduces numerical errors near the poles.
!
! m is the order, while l is the degree (Y_lm)
!
! tvec is a nt length 1D array of theta values     [0,  pi]
! pvec is a np length 1D array of longitude values [0,2*pi]
! Y is a 2D array of size (np,nt)
!
! This code uses `do concurrent' and OpenACC for GPU-acceleration
! when using a compatible compiler (e.g. NVIDIA's nvfortran).
! For perfromance, it is best if lvec, tvec, and Y are already
! on the GPU.
!
! This code was adapted from the MATLAB code sphharm.m from the
! Chebfun package, Copyright 2018 by The University of Oxford
! and The Chebfun Developers.
! See http://www.chebfun.org/ for Chebfun information.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: pi = 3.141592653589793_r_typ
      real(r_typ), parameter :: one = 1.0_r_typ
      real(r_typ), parameter :: two = 2.0_r_typ
      real(r_typ), parameter :: three = 3.0_r_typ
      real(r_typ), parameter :: four = 4.0_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: l, m, nt, np
      real(r_typ), dimension(nt) :: tvec
      real(r_typ), dimension(np) :: pvec
      real(r_typ), dimension(np,nt) :: Y
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nt) :: cos_tvec, Pold, Poldold, Pl
      real(r_typ) :: anm, bnm, one_over_sqrt_4pi
      integer :: i, im, il, j, pos, abs_m
!
!-----------------------------------------------------------------------
!
! ****** Set |m|.
!
      abs_m = abs(m)
!
! ****** Create the temporary arrays on the GPU.
!
!$acc enter data create(cos_tvec, Pold, Poldold, Pl)
!
      do concurrent(i=1:nt)
!
! ****** Store cos(theta).
!
        cos_tvec(i) = cos(tvec(i))
!
! ****** Initialize P^0_0/u^0.
!
        Pold(i) = one
!
! ****** Initialize the recurrence.
!
        Poldold(i) = 0.
      enddo
!
! ****** Compute P^m_m/u^m.
!
      if (abs_m.gt.0) then
!
! ****** Compute P^1_1/u^1.
!
        do concurrent(i=1:nt)
          Pold(i) = sqrt(three)
        enddo
!
! ****** Compute P^m_m/u^m.
!
        do im=2,abs_m
          do concurrent(i=1:nt)
            Pold(i) = sqrt((two*im+one)/(two*im)) * Pold(i)
          enddo
        enddo
      end if
!
! ****** Compute P^m_l/u^m for m+1<=l.
!
      do il=abs_m+1,l
!
        anm = sqrt( ((two*il-one)*(two*il+one)) / &
                      ((il-abs_m)*(il+abs_m))    )
!
        bnm = sqrt( ((two*il+one)*(il+abs_m-one)*(il-abs_m-one)) / &
                    ((il-abs_m)*(il+abs_m)*(two*il-three))        )
!
        do concurrent(i=1:nt)
          Pl(i) = anm*cos_tvec(i)*Pold(i) - bnm*Poldold(i)
          Poldold(i) = Pold(i)
          Pold(i) = Pl(i)
        enddo
      enddo
!
! ****** Normalize and compute associated Legendre polynomials.
! ****** Note that there is no sqrt(2) term here as it is included
! ****** in the above recursion as it computed the fully
! ****** normalized P-bar.  The -1^m and 1/4pi terms
! ****** are for the spherical harmonic normalization.
!
      one_over_sqrt_4pi = one/sqrt(four*pi)
      do concurrent(i=1:nt)
        Pold(i) = ((-one)**abs_m) * (sin(tvec(i))**abs_m) * &
                  one_over_sqrt_4pi * Pold(i)
      enddo
!
! ****** Determine if the cos or sin term should be used.
!
      pos = abs( max(0,sign(1,m+1)) )
!
! ****** Compute the spherical harmonic.
!
      do concurrent(j=1:nt,i=1:np)
        Y(i,j) = Pold(j)*(                         &
                         (pos)*cos(m*pvec(i))      &
                   + (one-pos)*sin(abs_m*pvec(i))  &
                          )
      enddo
!
! ****** Remove the temporary arrays from the GPU:
!
!$acc exit data delete(cos_tvec, Pold, Poldold, Pl)
!
end subroutine
!#######################################################################
subroutine set_unchanging_quantities
!
!-----------------------------------------------------------------------
!
! ****** Set unchanging quantities.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Define the diffusion coefficent and matrix.
!
      if (advance_diffusion) call load_diffusion
!
      if (advance_flow.and.flow_num_method.eq.3) call load_weno
!
end subroutine
!#######################################################################
subroutine output_step
!
!-----------------------------------------------------------------------
!
! ****** Output.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use mpidefs
      use timing
      use globals
      use mesh
      use fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: t1
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
! RMC: This should eventually be for all ranks to write out all realizations.
      if (iamp0) call output_histories
!
      call output_map
!
!     call output_flow(?)
!
!     call output_restart
!
      wtime_io = wtime_io + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine output_histories
!
!-----------------------------------------------------------------------
!
! ****** Write out histories.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use output
      use globals
      use sts
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      integer*8 :: niters
      character(*), parameter :: FMT='(i10,a,6(1pe23.15e3,a),i15)'
      character(*), parameter :: FMT2='(i10,a,15(1pe23.15e3,a))'
!
!-----------------------------------------------------------------------
!
      call FFOPEN (IO_HIST_NUM,io_hist_num_filename,'a',ierr)
!
      if (diffusion_num_method.eq.1) then
        niters = n_stable_diffusion_cycles
      else
        niters = sts_s
      end if
!
      write(IO_HIST_NUM,FMT) ntime,' ',                               &
                             time,' ',dtime_global,' ',               &
                             dtmax_flow,' ',dtime_advection_used, ' ',&
                             dtime_diffusion_euler,' ',               &
                             dtime_diffusion_used, ' ',niters
!
      close(IO_HIST_NUM)
!
      call FFOPEN (IO_HIST_SOL,io_hist_sol_filename,'a',ierr)

      write(IO_HIST_SOL,FMT2) ntime,' ',time,' ',                     &
                              h_minbr,' ',h_maxbr,' ',h_minabsbr,' ', &
                              h_fluxp,    ' ', h_fluxm,    ' ',       &
                              h_fluxp_pn, ' ', h_fluxm_pn, ' ',       &
                              h_fluxp_ps, ' ', h_fluxm_ps, ' ',       &
                              h_area_pn,  ' ', h_area_ps,  ' ',       &
                              h_eq_dipole,' ', h_ax_dipole,' ',       &
                              h_valerr

!
      close(IO_HIST_SOL)
!
end subroutine
!#######################################################################
subroutine output_map
!
!-----------------------------------------------------------------------
!
! ****** Output map to disk if it is time to do so.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use output
      use globals
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr = 0
!
! ****** File name.
!
      character(6) :: idx_str
      character(256) :: full_filename = ' '
      character(256) :: base_filename = ' '
      character(*), parameter :: FMT='(i10,a,1pe25.15,a,a)'
!
!-----------------------------------------------------------------------
!
      if (output_map_idx_cadence .gt. 0) then
        if (MODULO(ntime,output_map_idx_cadence).eq.0) then
          output_current_map = .true.
        end if
      end if
!
      if (output_map_time_cadence .gt. 0.0) then
        if (time .ge. time_start+idx_out*output_map_time_cadence) then
          output_current_map = .true.
        end if
      end if
!
      if (output_current_map) then
!
        idx_out = idx_out + 1
!
! ****** Formulate file name.
!
        write (idx_str,'(i6.6)') idx_out
!
        full_filename = trim(output_map_directory)//'/'// &
                   trim(output_map_root_filename)// &
                   '_idx'//idx_str//'.h5'
!
        base_filename = trim(output_map_root_filename)// &
                   '_idx'//idx_str//'.h5'
!
! ****** Write out the map.
!
        call write_map (full_filename)
!
! ****** Reset output flag.
!
        output_current_map = .false.
!
! ****** Record output in output text file log.
!
        call FFOPEN (IO_MAP_OUT_LIST, &
                     io_map_output_list_filename,'a',ierr)
!
        write (IO_MAP_OUT_LIST,FMT) idx_out,'   ',time,'    ', &
                                    trim(base_filename)
!
        close(IO_MAP_OUT_LIST)
!
      end if
!
end subroutine
!#######################################################################
subroutine write_map (fname)
!
!-----------------------------------------------------------------------
!
! ****** Write out current map "f" to disk.
!
!-----------------------------------------------------------------------
!
      use input_parameters, ONLY : output_map_2d,n_realizations
      use number_types
      use output
      use mpidefs
      use fields, ONLY : f
      use mesh, ONLY : t,npm,ntm,nr
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      integer :: ierr = 0
      integer :: i, lsbuf, irank
      integer, dimension(:), allocatable :: displ,lbuf
      real(r_typ), dimension(:), allocatable :: gfout,tfout,rscale
      real(r_typ), dimension(:,:,:), allocatable :: gout
!
!-----------------------------------------------------------------------
!
!$acc update self(f)
!
! ****** Write out file in pt format with single point overlap in phi.
!
      allocate (fout(npm-1,ntm,nr))
!
!   WARNING: -stdpar may eventually GPU-ize TRANSPOSE.
!   For now, it probably in-lines so OK, and/or requires
!   the tensor module to be explicitly used.
!   For now, this should be computed on the CPU no matter what...
!
      do i=1,nr
        fout(:,:,i) = TRANSPOSE(f(:,1:npm-1,i))
      enddo
!
        if (output_map_2d.and.n_realizations.eq.1) then
          call write_2d_file (fname,npm-1,ntm,fout(:,:,1),pout,t,ierr)
        else
!
          wtime_tmp_mpi = MPI_Wtime()
          lsbuf=(npm-1)*ntm*nr
          allocate (tfout(lsbuf))
          tfout(:)=reshape(fout(:,:,:),(/lsbuf/))
          if (iamp0) then
            allocate (gfout((npm-1)*ntm*n_realizations))
            allocate (displ(0:nproc-1))
            allocate (lbuf(0:nproc-1))
            displ(0)=0
            do irank=0,nproc-1
              lbuf(irank)=(npm-1)*ntm*nr_arr(irank)
            enddo
            do irank=1,nproc-1
              displ(irank)=displ(irank-1)+lbuf(irank-1)
            enddo
          end if
!
          call MPI_Gatherv (tfout,lsbuf,ntype_real,gfout,lbuf,       &
                    displ,ntype_real,0,MPI_COMM_WORLD,ierr)
          deallocate(tfout)
          wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
          if (iamp0) then
            allocate (rscale(n_realizations))
            do i=1,n_realizations
              rscale(i) = real(i,r_typ)
            enddo
            allocate (gout((npm-1),ntm,n_realizations))
            gout=reshape(gfout(:),(/(npm-1),ntm,n_realizations/))
            call write_3d_file(fname,npm-1,ntm,n_realizations,gout,pout,t,   &
                               rscale,ierr)
            deallocate (gout)
            deallocate (gfout)
            deallocate (displ)
            deallocate (lbuf)
            deallocate (rscale)
          end if
!
        endif
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_MAP:'
        write (*,*) '### Could not write the output data set!'
        write (*,*) 'IERR: ',ierr
        write (*,*) 'File name: ', fname
        call endrun(.true.)
      end if
!
      deallocate (fout)
!
end subroutine
!#######################################################################
subroutine analysis_step
!
!-----------------------------------------------------------------------
!
! ****** Compute analysis of run.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use mpidefs
      use timing
      use globals
      use mesh
      use fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: t1
      real(r_typ) :: sumfval2, fv, fv2, eqd1, eqd2
      real(r_typ) :: tav, da_t, da_p, sn_t, d_t, cs_t, cs_p, sn_p
      real(r_typ), dimension(:,:,:), allocatable :: fval
      integer :: i, j, k
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
! ****** If running validation, compute current exact solution.
!
      if (validation_run .eq. 1) then
        allocate (fval(npm-1,ntm,nr))
!$acc enter data create (fval)
!
        do concurrent (k=1:nr,j=1:ntm,i=1:npm-1)
          fval(i,j,k) = fval_u0(i,j,k)*exp(-42.0_r_typ*diffusion_coef_constant* &
                                       diffusion_coef_factor*time)
        enddo
!
      end if
!
! ***** Get max and min metrics.
!
      h_minbr = large_value
      h_maxbr = small_value
      h_minabsbr = large_value
      h_valerr = 0.
      sumfval2 = 0.
!
      do k=1,1
!$omp parallel do collapse(2) default(shared) &
!$omp         reduction(+:sumfval2,h_valerr) &
!$omp         reduction(max:h_maxbr) reduction(min:h_minbr,h_minabsbr)
!$acc parallel loop collapse(2) default(present) &
!$acc&         reduction(+:sumfval2,h_valerr) &
!$acc&         reduction(max:h_maxbr) reduction(min:h_minbr,h_minabsbr)
        do j=1,npm-2
          do i=1,ntm
            h_minbr = min(f(i,j,k),h_minbr)
            h_maxbr = max(f(i,j,k),h_maxbr)
            h_minabsbr = min(abs(f(i,j,k)),h_minabsbr)
            if (validation_run .eq. 1) then
              fv = (f(i,j,k) - fval(j,i,k))**2
              fv2 = fval(j,i,k)**2
            else
              fv = 0.
              fv2 = 0.
            end if
            h_valerr = h_valerr + fv
            sumfval2 = sumfval2 + fv2
          enddo
        enddo
!$omp end parallel do
      enddo
!
      if (validation_run .eq. 1) then
        h_valerr = sqrt(h_valerr/sumfval2)
!$acc exit data delete (fval)
      end if
!
! ****** Get integrated metrics.
!
      h_fluxp = 0.
      h_fluxm = 0.
      h_fluxp_pn = 0.
      h_fluxm_pn = 0.
      h_fluxp_ps = 0.
      h_fluxm_ps = 0.
      h_area_pn = 0.
      h_area_ps = 0.
      eqd1 = 0.
      eqd2 = 0.
      h_eq_dipole = 0.
      h_ax_dipole = 0.
!
      do k=1,1
!$omp parallel do   collapse(2) default(shared)                 &
!$omp reduction(+:h_fluxp,h_fluxm,h_fluxp_pn,h_fluxm_pn,        &
!$omp             h_fluxp_ps,h_fluxm_ps,eqd1,eqd2,h_ax_dipole,h_area_ps)
!$acc parallel loop collapse(2) default(present) &
!$acc& reduction(+:h_fluxp,h_fluxm,h_fluxp_pn,h_fluxm_pn,h_area_pn, &
!$acc&             h_fluxp_ps,h_fluxm_ps,h_area_ps,eqd1,eqd2,h_ax_dipole)
        do i=1,npm-1
          do j=1,ntm
!
            if (j.eq.1) then
              tav=half*(t(1)+t(2))
              sn_t = sin(tav)
              cs_t = cos(tav)
              d_t = dth(2)
              da_t=quarter*sn_t*d_t
            else if (j.eq.ntm) then
              tav=half*(t(ntm)+t(ntm-1))
              sn_t = sin(tav)
              cs_t = cos(tav)
              d_t = dth(ntm)
              da_t=quarter*sn_t*d_t
            else
              sn_t = sin(t(j))
              cs_t = cos(t(j))
              d_t = dt(j)
              da_t = sn_t*d_t
            end if
!
            if (i.eq.1) then
              da_p=half*dph(1)
            else if (i.eq.npm-1) then
              da_p=half*dph(npm-1)
            else
              da_p=dp(i)
            end if
            cs_p = cos(p(i))
            sn_p = sin(p(i))
!
            fv = f(j,i,k)*da_t*da_p
!
! ****** Fluxes.
!
            if (fv.gt.0.) then
              h_fluxp = h_fluxp + fv
            else
              h_fluxm = h_fluxm + fv
            end if
!
! ****** Polar fluxes and areas.
!
            if (t(j).le.pole_flux_lat_limit*d2r) then
              if (fv.gt.0.) then
                h_fluxp_pn = h_fluxp_pn + fv
              else
                h_fluxm_pn = h_fluxm_pn + fv
              end if
              h_area_pn = h_area_pn + da_t*da_p
            end if
!
            if (t(j).ge.pi-pole_flux_lat_limit*d2r) then
              if (fv.gt.0.) then
                h_fluxp_ps = h_fluxp_ps + fv
              else
                h_fluxm_ps = h_fluxm_ps + fv
              end if
              h_area_ps = h_area_ps + da_t*da_p
            end if
!
! ****** Dipoles.
!
            eqd1 = eqd1 + f(j,i,k)*sn_t*cs_p*da_t*da_p
            eqd2 = eqd2 + f(j,i,k)*sn_t*sn_p*da_t*da_p
!
            h_ax_dipole = h_ax_dipole + f(j,i,k)*cs_t*da_t*da_p
!
          enddo
        enddo
!$omp end parallel do
      enddo
!
! ****** Set equatorial dipole strength.
!
      eqd1 = three_quarter*pi_i*eqd1
      eqd2 = three_quarter*pi_i*eqd2
      h_eq_dipole = SQRT(eqd1**2 + eqd2**2)
!
! ****** Set axial dipole strength.
!
      h_ax_dipole = three_quarter*pi_i*h_ax_dipole
!
! ****** Set fluxes to be in units of Mx
!
      h_fluxp    = rsun_cm2*h_fluxp
      h_fluxm    = rsun_cm2*h_fluxm
      h_fluxp_pn = rsun_cm2*h_fluxp_pn
      h_fluxm_pn = rsun_cm2*h_fluxm_pn
      h_fluxp_ps = rsun_cm2*h_fluxp_ps
      h_fluxm_ps = rsun_cm2*h_fluxm_ps
!
! ****** Set polar areas to be in units of cm.
!
      h_area_pn = rsun_cm2*h_area_pn
      h_area_ps = rsun_cm2*h_area_ps
!
      wtime_analysis = wtime_analysis + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine update_step
!
!-----------------------------------------------------------------------
!
! ****** Update prescribed quantities and time step.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use input_parameters
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: t1
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
      if (advance_flow)    call update_flow
      if (advance_source)  call update_source
      if (assimilate_data) call update_field
!
      call update_timestep
!
      wtime_update = wtime_update + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine update_flow
!
!-----------------------------------------------------------------------
!
! ****** Update flows.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use mesh
      use fields
      use globals
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: j,k,i
      real(r_typ) :: br_avg
      real(r_typ) :: vrigid_pole
      real(r_typ), dimension(:), allocatable :: vt_npole,vt_spole
      real(r_typ), dimension(:), allocatable :: vp_npole,vp_spole
      real(r_typ) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,cp,cph
!
!-----------------------------------------------------------------------
!
! ***** Reset the flow.
!
      if (flow_needs_updating) then
!
!       Assume no more updates needed.  If needed,
!       this will be set to true.
        flow_needs_updating = .false.
!       Since the flow is changing, we need to update time step.
        timestep_flow_needs_updating = .true.
        timestep_needs_updating = .true.
!
! ***** Reset flows.
!
        do concurrent (i=1:nr,k=1:npm,j=1:nt)
          vt(j,k,i) = 0.
        enddo
        do concurrent (i=1:nr,k=1:np,j=1:ntm)
          vp(j,k,i) = 0.
        enddo
!
! ***** Add in flow from files. Note this overwrites,
! ***** so needs to be called first!
!
        if (use_flow_from_files) then
          call add_flow_from_files
          flow_needs_updating = .true.
        end if
!
! ***** Add in flows from ConFlow algorithm.
!       TO BE DEV: if (flow_from_conflow) call add_flow_from_conflow (metadata_structure,flow_array)
!
! ***** Add in analytic meridianal flow model.
!
        if (flow_mf_model.eq.1) call add_flow_meridianal_aft
!
! ***** Add in analytic differential rotation model.
!
        if (flow_dr_model.eq.1) call add_flow_differential_rotation_aft
!
! ***** Add in constant angular phi velocity (rigid rotation)
!
        if (flow_vp_rigid_omega.ne.0.) then
          do concurrent (i=1:nr,k=1:np,j=1:ntm)
            vp(j,k,i) = vp(j,k,i) + flow_vp_rigid_omega*km_s_to_rs_hr*st(j)
          enddo
        end if
!
! ***** Add in constant theta velocity.
!
        if (flow_vt_const .ne. 0.) then
          do concurrent (i=1:nr,k=1:npm,j=1:nt)
            vt(j,k,i) = vt(j,k,i) + flow_vt_const*km_s_to_rs_hr
          enddo
        end if
!
! ***** Add in rotated rigid rotation velocity.
!
        if (flow_rigid_omega .ne. 0.) then
          vrigid_pole = flow_rigid_omega*km_s_to_rs_hr
          do concurrent (i=1:nr,k=1:npm,j=1:nt)
            cp = cos(p(k))
            p1 = -2*cth(j)*abs(sth(j)*cp)
            p2 = (cp*cth(j)**2-cp)*abs(cp)**3
            p3 = (-3*cp*cth(j)**2+2*cp)*abs(cp)
            p4 = sign(one,cp)*cth(j)**2
            p5 = 2*abs(cp)**3*cth(j)*abs(sth(j))
            p6 = (-cp**4+3*cp**2-1)*cth(j)**2
            p7 = (-2*cp**3*sth(j)+2*sth(j)*cp)*cth(j)
            p8 = cp**4-2*cp**2+1
            p9 = abs(sth(j)*cp)*sin(p(k))*vrigid_pole
            p10 = (cp**4-3*cp**2+1)*cth(j)**2
            p11 = (2*cp**3*sth(j)-2*sth(j)*cp)*cth(j)
            p12 = cp**4+2*cp**2
            vt(j,k,i) = vt(j,k,i)+(-((p1+(p2+p3+p4)*sign(one,sth(j))    &
                  +p5+p6+p7+p8)*p9)/SQRT(p10+p11-p12))
            if ((abs(vt(j,k,i))>10.0**5) .or. (abs(vt(j,k,i))<-10.0**(-5))) then
              vt(j,k,i)=0
            end if
          enddo
          do concurrent (i=1:nr,k=1:np,j=1:ntm)
            cph = cos(ph(k))
            p1 = vrigid_pole*(cph**2*ct(j)+st(j)*cph-ct(j))*sign(one,cph)*abs(st(j))
            p2 = (cph**4-3*cph**2+1)*ct(j)**2
            p3 = (2*cph**3*st(j)-2*st(j)*cph)*ct(j)
            p4 = cph**4+2*cph**2
            vp(j,k,i) = vp(j,k,i)- p1/SQRT(p2+p3-p4)
            if ((abs(vp(j,k,i))>10.0**5) .or. (abs(vp(j,k,i))<-10.0**(-5))) then
              vp(j,k,i)=0
            end if
         end do
        end if
!
! ***** Attenuate velocity based on Br.
!
        if (flow_attenuate) then
!
! ****** Need to interpolate f at half-mesh point.
!
          do concurrent (i=1:nr,k=2:npm-1,j=2:nt-1)
! RMC: Should maybe use max(|br|) here instead of avg?
            br_avg = half*(f(j-1,k,i) + f(j,k,i))
            vt(j,k,i) = vt(j,k,i)*(one - tanh(abs(br_avg)*fivehundred_i))
          enddo
!
          do concurrent (i=1:nr,k=2:np-1,j=2:ntm-1)
            br_avg = half*(f(j,k-1,i) + f(j,k,i))
            vp(j,k,i) = vp(j,k,i)*(one - tanh(abs(br_avg)*fivehundred_i))
          enddo
!
! ****** Set poles. Maybe not needed since field is weak at pole?
!
          allocate (vt_npole(nr))
          allocate (vt_spole(nr))
          allocate (vp_npole(nr))
          allocate (vp_spole(nr))
!$acc enter data create(vt_npole,vt_spole,vp_npole,vp_spole)
!
          do concurrent (i=1:nr)
            vt_npole(i) = 0.
            vt_spole(i) = 0.
            vp_npole(i) = 0.
            vp_spole(i) = 0.
          enddo
!
!$omp parallel do   collapse(2) default(shared)
!$acc parallel loop collapse(2) default(present)
          do i=1,nr
            do k=2,npm-1
!$omp atomic update
!$acc atomic update
              vt_npole(i) = vt_npole(i) + vt(   2,k,i)*dp(k)
!$omp atomic update
!$acc atomic update
              vt_spole(i) = vt_spole(i) + vt(ntm1,k,i)*dp(k)
            enddo
          enddo
!$omp end parallel do
!
!$omp parallel do   collapse(2) default(shared)
!$acc parallel loop collapse(2) default(present)
          do i=1,nr
            do k=2,npm1
!$omp atomic update
!$acc atomic update
              vp_npole(i) = vp_npole(i) + vp(    2,k,i)*dph(k)
!$omp atomic update
!$acc atomic update
              vp_spole(i) = vp_spole(i) + vp(ntm-1,k,i)*dph(k)
            enddo
          enddo
!$omp end parallel do
!
! ****** Attenuate pole value to set new v bc (f has same val for all pole pts, use 1)
!
          do concurrent (i=1:nr)
            vt_npole(i) = vt_npole(i)*twopi_i*(one - tanh(f(1,1,i)*fivehundred_i))
            vt_spole(i) = vt_spole(i)*twopi_i*(one - tanh(f(ntm,1,i)*fivehundred_i))
            vp_npole(i) = vp_npole(i)*twopi_i*(one - tanh(f(1,1,i)*fivehundred_i))
            vp_spole(i) = vp_spole(i)*twopi_i*(one - tanh(f(ntm,1,i)*fivehundred_i))
          enddo
!
! ****** Now set the pole:
!
          do concurrent (i=1:nr,k=1:npm)
            vt( 1,k,i) = two*vt_spole(i)-vt(   2,k,i)
            vt(nt,k,i) = two*vt_spole(i)-vt(ntm1,k,i)
          enddo
!
          do concurrent (i=1:nr,k=1:np)
            vp(  1,k,i) = two*vp_spole(i)-vp(    2,k,i)
            vp(ntm,k,i) = two*vp_spole(i)-vp(ntm-1,k,i)
          enddo
!
!$acc exit data delete(vt_npole,vt_spole,vp_npole,vp_spole)
          deallocate (vt_npole)
          deallocate (vt_spole)
          deallocate (vp_npole)
          deallocate (vp_spole)
!
! ****** Set periodicity
!
          call set_periodic_bc_3d (vt,nt,npm,nr)
          call set_periodic_bc_3d (vp,ntm,np,nr)
!
! ****** Since this needs to be done each timestep as Br changes,
! ****** the full velocity profile needs to be reloaded and
! ****** re-attenuated.
!
          flow_needs_updating = .true.
!
        end if
!
! ****** If using WENO3 for flow, set new alpha values for LF.
!
        if (flow_num_method.eq.3) call set_lf_alpha
!
      end if
!
end subroutine
!#######################################################################
subroutine add_flow_from_files
!
!-----------------------------------------------------------------------
!
! ****** Update flows from file(s).
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use flow_from_files
      use fields
      use globals
      use mpidefs
      use mesh
      use read_2d_file_interface
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr = 0
      real(r_typ), dimension(:,:), allocatable :: new_flow_t,new_flow_p
      real(r_typ), dimension(:), allocatable :: s1,s2
      integer :: ln1,ln2,i,j,k
      character(1024) :: flowfile_t = ' '
      character(1024) :: flowfile_p = ' '
!
!-----------------------------------------------------------------------
!
     if (time .ge. time_of_next_input_flow) then
!
! ****** Read the flow data.
!
       flowfile_t = TRIM(flow_root_dir)//"/"&
                  //TRIM(flow_files_rel_path_t(current_flow_input_idx))

       flowfile_p = TRIM(flow_root_dir)//"/"&
                  //TRIM(flow_files_rel_path_p(current_flow_input_idx))
!
! ****** VT (NT,NPM) ******
!
       if (iamp0) then
         call read_2d_file (flowfile_t,ln1,ln2,new_flow_t,s1,s2,ierr)
         deallocate(s1)
         deallocate(s2)
       endif
!
       wtime_tmp_mpi = MPI_Wtime()
       if (.not.iamp0) then
        allocate (new_flow_t(npm1,nt))
       endif
       call MPI_Bcast (new_flow_t,flow_t_size,ntype_real,0,MPI_COMM_WORLD,ierr)
       wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!$acc enter data copyin (new_flow_t)
!
! ****** copy, transpose, and convert units.
!
       do concurrent (k=1:nr,j=1:npm1,i=1:nt)
         vt(i,j,k) = m_s_to_rs_hr*new_flow_t(j,i)
       enddo
       ! Set phi ghost cell.
       do concurrent (k=1:nr,i=1:nt)
         vt(i,npm,k) = vt(i,2,k)
       enddo
!$acc exit data delete (new_flow_t)
       deallocate(new_flow_t)
!
! ****** VP (NTM,NP)
!
       if (iamp0) then
         call read_2d_file (flowfile_p,ln1,ln2,new_flow_p,s1,s2,ierr)
         deallocate(s1)
         deallocate(s2)
       endif
!
       wtime_tmp_mpi = MPI_Wtime()
       if (.not.iamp0) then
        allocate (new_flow_p(np,ntm1))
       endif
       call MPI_Bcast (new_flow_p,flow_p_size,ntype_real,0,MPI_COMM_WORLD,ierr)
       wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!$acc enter data copyin (new_flow_p)
!
! ****** copy, transpose, and convert units.
!.
       do concurrent (k=1:nr,j=1:np,i=1:ntm)
         vp(i,j,k) = m_s_to_rs_hr*new_flow_p(j,i)
       enddo
!$acc exit data delete (new_flow_p)
       deallocate(new_flow_p)
!
! ******  TODO:  Add error checking here (make sure scales match run
!                until interp is added)
!
! ****** Update position in map list and next map time.
!
       current_flow_input_idx = current_flow_input_idx + 1
       if (verbose) print*, 'FlowIdx:',current_flow_input_idx
!
! ****** If there are no more flow files to load, loop back to
!        the beginning (use index 2 to avoid dt=0 situation).
!
       if (current_flow_input_idx .gt. num_flows_in_list) then
         if (verbose) print*, 'Resetting FlowIdx to 2'
         current_flow_input_idx = 2
         flow_times_actual_ut_jd(:) = flow_times_actual_ut_jd0(:) &
                                      + time/twentyfour
       end if
!
       time_of_next_input_flow = twentyfour * &
                       flow_times_actual_ut_jd(current_flow_input_idx) &
                     - flow_time_initial_hr
!
! ****** Save read flows.  This is needed because one might be
! ****** attenuating velocity when its not yet time to read a new file.
! ****** In that case, v[tp] would be set to 0.  Instead, set to
! ****** old value below.
!
       do concurrent (k=1:nr,j=1:npm,i=1:nt)
         flow_from_file_vt_old(i,j,k) = vt(i,j,k)
       enddo
!
       do concurrent (k=1:nr,j=1:np,i=1:ntm)
         flow_from_file_vp_old(i,j,k) = vp(i,j,k)
       enddo
!
     else
!
! ****** Set to last read in flow since its not time to change.
!
       do concurrent (k=1:nr,j=1:npm,i=1:nt)
         vt(i,j,k) = flow_from_file_vt_old(i,j,k)
       enddo
!
       do concurrent (k=1:nr,j=1:np,i=1:ntm)
         vp(i,j,k) = flow_from_file_vp_old(i,j,k)
       enddo
!
     end if
!
end subroutine
!#######################################################################
subroutine load_flow_from_files
!
!-----------------------------------------------------------------------
!
! ****** Load flow_from_files
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use globals
      use output
      use mesh
      use flow_from_files
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: io = 0
      integer :: i,j,k
      character(*), parameter :: &
                    FMT = '(F11.5,1X,A11,1X,A11)'
!
!-----------------------------------------------------------------------
!
! ****** Open file list of input maps.
!
      OPEN (unit=IO_DATA_IN,file=flow_list_filename, &
            form="FORMATTED",status='old')
!
! ****** Get the number of entries in the list (skip first row)
!
      num_flows_in_list = -1
      do
        READ (IO_DATA_IN,'(a)',iostat=io)
        if (io/=0) exit
        num_flows_in_list = num_flows_in_list + 1
      enddo
      REWIND (IO_DATA_IN)
!
! ****** Allocate space to store the list.
!
      allocate(flow_times_actual_ut_jd     (num_flows_in_list))
      allocate(flow_times_actual_ut_jd0    (num_flows_in_list))
      allocate(flow_files_rel_path_t       (num_flows_in_list))
      allocate(flow_files_rel_path_p       (num_flows_in_list))
!
! ****** Read the list (skip first row).
!
      READ (IO_DATA_IN,*)
      do i=1,num_flows_in_list
        READ (IO_DATA_IN,FMT) flow_times_actual_ut_jd0(i),    &
                              flow_files_rel_path_t(i),       &
                              flow_files_rel_path_p(i)
      enddo
      flow_times_actual_ut_jd(:) = flow_times_actual_ut_jd0(:) &
                                   + time_start/twentyfour
      CLOSE (IO_DATA_IN)
!
! ******* Initialize flow time.
!
      flow_time_initial_hr = twentyfour*flow_times_actual_ut_jd0(1)
!
      current_flow_input_idx = 1
!
      time_of_next_input_flow = time_start
!
! ****** Allocate arrays to store old flows.
!
      allocate (flow_from_file_vt_old(nt,npm,nr))
      allocate (flow_from_file_vp_old(ntm,np,nr))
!$acc enter data create (flow_from_file_vt_old,flow_from_file_vp_old)
!
      do concurrent (i=1:nr,k=1:npm,j=1:nt)
        flow_from_file_vt_old(j,k,i) = 0.
      enddo
      do concurrent (i=1:nr,k=1:np,j=1:ntm)
        flow_from_file_vp_old(j,k,i) = 0.
      enddo
!
! ******* Print some diagnostics.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'FLOWS FROM FILE'
        write (*,*) 'Number of flows in flow list: ', &
                     num_flows_in_list
        write (*,*) 'File name of first theta flow file:     ', &
                     trim(flow_files_rel_path_t(1))
        write (*,*) 'File name of first phi flow file:     ', &
                     trim(flow_files_rel_path_p(1))
      end if
!
end subroutine
!#######################################################################
subroutine update_source
!
!-----------------------------------------------------------------------
!
! ****** Update source quantities.
!
!-----------------------------------------------------------------------
!
      use input_parameters
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     - FILE
!     - Random flux
!     - Artificial AR emerge
!
end subroutine
!#######################################################################
subroutine load_data_assimilation
!
!-----------------------------------------------------------------------
!
! ****** Load data assimilation.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use globals
      use output
      use data_assimilation
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: io = 0
      integer :: i
      character(*), parameter :: &
                             FMT = '(A19,1X,A19,1X,F13.5,1X,A)'
!
!-----------------------------------------------------------------------
!
! ****** Open file list of input maps.
!
      OPEN (unit=IO_DATA_IN,file=assimilate_data_map_list_filename, &
            form="FORMATTED",status='old')
!
! ****** Get the number of entries in the list (skip first row)
!
      num_maps_in_list = -1
      do
        READ (IO_DATA_IN,'(a)',iostat=io)
        if (io/=0) exit
        num_maps_in_list = num_maps_in_list + 1
      enddo
      REWIND (IO_DATA_IN)
!
! ****** Allocate space to store the list.
!
      allocate(map_times_actual_ut_jd     (num_maps_in_list))
      allocate(map_times_requested_ut_str (num_maps_in_list))
      allocate(map_times_actual_ut_str    (num_maps_in_list))
      allocate(map_files_rel_path         (num_maps_in_list))
!
! ****** Read the list (skip first row).
!
      READ (IO_DATA_IN,*)
      do i=1,num_maps_in_list
        READ (IO_DATA_IN,FMT) map_times_requested_ut_str(i),  &
                              map_times_actual_ut_str(i),     &
                              map_times_actual_ut_jd(i),      &
                              map_files_rel_path(i)
      enddo
      CLOSE (IO_DATA_IN)
!
! ******* Initialize assimilation time.
!  RMC: This only works for exact times!  If inbetween files,
!       need to interpolate two assimilation frames, or modify
!       start_time?
!
      map_time_initial_hr = twentyfour*map_times_actual_ut_jd(1)
!
      i = MINLOC(ABS( (twentyfour*map_times_actual_ut_jd(:) &
                       - map_time_initial_hr)               &
                      - time_start),1)
!
      current_map_input_idx = i
!
      time_of_next_input_map = time_start
!
! ******* Print some diagnostics.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'DATA ASSIMILATION LOADED'
        write (*,*) 'Number of maps in map list: ', &
                     num_maps_in_list
        write (*,*) 'Start date:     ', &
                     trim(map_times_actual_ut_str(1))
        write (*,*) 'End   date:     ', &
                     trim(map_times_actual_ut_str(num_maps_in_list))
        write (*,*) 'File name of first map:     ', &
                     trim(map_files_rel_path(1))
        write (*,*) 'File name of first used map:     ', &
                     trim(map_files_rel_path(i))
      end if
!
end subroutine
!#######################################################################
subroutine assimilate_new_data (new_data)
!
!-----------------------------------------------------------------------
!
! ****** Assimilate data.
! ****** This assumes the uncertainty is in the second slice
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use globals
      use constants, ONLY : one,half
      use fields, ONLY : f
      use mesh
      use data_assimilation
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: j,k,i
      real(r_typ), dimension(:,:,:,:), allocatable :: new_data
!
!-----------------------------------------------------------------------
!
      do concurrent (i=1:nr,k=1:npm1,j=1:ntm)
         f(j,k,i) = (one - new_data(k,j,2,i))*         f(j,k,i) +  &
                    (      new_data(k,j,2,i))*new_data(k,j,1,i)
      enddo
!
! ****** Set periodicity (since k=np (npm) not set above).
!
      call set_periodic_bc_3d (f,ntm,npm,nr)
!
end subroutine
!#######################################################################
subroutine update_field
!
!-----------------------------------------------------------------------
!
! ****** Update field through data assimilation.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use mpidefs
      use globals
      use output
      use mesh
      use read_3d_file_interface
      use assimilate_new_data_interface
      use data_assimilation
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr = 0
      real(r_typ), dimension(:,:,:,:), allocatable :: new_data
      real(r_typ), dimension(:,:,:), allocatable :: new_data2d
      real(r_typ), dimension(:), allocatable :: s1,s2,s3
      integer :: ln1,ln2,nslices,i
      character(1024) :: mapfile = ' '
!
!-----------------------------------------------------------------------
!
     if (time .ge. time_of_next_input_map) then
!
! ****** Read the map data.
!
       if (iamp0) then
         mapfile = TRIM(assimilate_data_map_root_dir)//"/"&
                   //TRIM(map_files_rel_path(current_map_input_idx))
!
         call read_3d_file (mapfile,ln1,ln2,nslices,new_data2d,s1,s2,s3,ierr)
         deallocate(s1)
         deallocate(s2)
         deallocate(s3)
       endif
!
       wtime_tmp_mpi = MPI_Wtime()
       call MPI_Bcast (ln1,1,ntype_real,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast (ln2,1,ntype_real,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast (nslices,1,ntype_real,0,MPI_COMM_WORLD,ierr)
!
       if (.not.iamp0) allocate (new_data2d(ln1,ln2,nslices))
!
       call MPI_Bcast (new_data2d,ln1*ln2*nslices,ntype_real,0,MPI_COMM_WORLD,ierr)
       wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi

!
! ******  TODO:  Add error checking here (make sure scales match run
!                until interp is added)
       allocate(new_data(ln1,ln2,nslices,nr))
       do i=1,nr
         new_data(:,:,:,i) = new_data2d(:,:,:)
       enddo
       deallocate(new_data2d)
!
!$acc enter data copyin(new_data)
!
! ****** Assimilate the data.
!
       call assimilate_new_data (new_data)
!
! ****** Clean up.
!
!$acc exit data delete(new_data)
       deallocate(new_data)
!
! ****** Update position in map list and next map time.
!
       current_map_input_idx = current_map_input_idx + 1
!
! ****** If there are no more data files to assimilate, don't
! ****** assimilate anymore.
!
       if (current_map_input_idx .gt. num_maps_in_list) then
         assimilate_data=.false.
       else
         time_of_next_input_map =                                     &
          twentyfour*map_times_actual_ut_jd(current_map_input_idx) &
          - map_time_initial_hr
       end if
!
     end if
!
end subroutine
!#######################################################################
subroutine update_timestep
!
!-----------------------------------------------------------------------
!
! ****** Update time step if needed.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use globals
      use output
      use data_assimilation
      use flow_from_files
      use sts
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dtime_old, dt_mxup
      logical :: time_step_changed = .false.
!
!-----------------------------------------------------------------------
!
! ****** Check:  Does the timestep need updating?
!
      if (timestep_needs_updating) then
!
! ****** Save previous timestep.
!
        dtime_old = dtime_global
!
! ****** Default timestep is one giant step for remaining time.
!
        dtime_global = time_end - time
        dtime_reason = 'end'
        time_step_changed = .true.
!
        if (advance_flow.and.timestep_flow_needs_updating) then
          call get_flow_dtmax (dtmax_flow)
          dtime_reason = 'flowcfl'
          timestep_flow_needs_updating = .false.
        else
          dtmax_flow = huge(one)
        end if
        dtime_global = MIN(dtime_global,dtmax_flow)
!
        ! Check for negative values with source?
!
! ****** Check to make sure dtime_global doesnt go up too fast.
!
        if (dt_max_increase_fac.gt.0.and.dtime_old.gt.0.) then
          dt_mxup = (one + dt_max_increase_fac)*dtime_old
        else
          dt_mxup = huge(one)
        end if
!
        dtime_global = MIN(dtime_global,dt_mxup,dt_max)
!
        if (dtime_global .lt. dt_min) then
          write (*,*) 'Warning! Time step is smaller than DTMIN!'
        end if
!
        timestep_needs_updating = .false.
!
      end if
!
! ****** Check for next output time, cut dt to match exactly.
!
      if (output_map_time_cadence .gt. 0.0) then
        if (time + dtime_global .ge. (time_start+idx_out*output_map_time_cadence)) then
          dtime_global = (time_start+idx_out*output_map_time_cadence) - time
          dtime_reason = 'output'
          timestep_needs_updating = .true.
          time_step_changed = .true.
        end if
      end if
!
! ****** Check for next map assimilation time, cut dt to match exactly.
!
      if (assimilate_data) then
        if (time + dtime_global .ge. time_of_next_input_map) then
          dtime_global = time_of_next_input_map - time
          dtime_reason = 'dataassim'
          timestep_needs_updating = .true.
          time_step_changed = .true.
        end if
      end if
!
! ****** Check for next flow input time, cut dt to match exactly.
!
      if (use_flow_from_files) then
        if (time + dtime_global .ge. time_of_next_input_flow) then
          dtime_global = time_of_next_input_flow - time
          dtime_reason = 'flowfile'
          timestep_needs_updating = .true.
          flow_needs_updating = .true.
          time_step_changed = .true.
        end if
      end if
!
! ****** Check for end time.
!
      if (time + dtime_global .ge. time_end) then
        dtime_global = time_end - time
        dtime_reason = 'end'
        timestep_needs_updating = .false.
        time_step_changed = .true.
      end if
!
! ****** Make sure STS is recomputed if timestep has changed.
!
      if (advance_diffusion .and. time_step_changed) then
        need_to_load_sts = .true.
      end if
!
end subroutine
!#######################################################################
subroutine flux_transport_step
!
!-----------------------------------------------------------------------
!
! ****** Evolve the field by one time step.
! ****** Strang splitting is in a ARDRA sequence.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use input_parameters
      use globals
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dtime_local, dtime_local_half, t1
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
      dtime_local = dtime_global
!
! ****** If only advancing flow or diffusion, no need to split.
!
      if (strang_splitting.and.advance_flow.and.advance_diffusion) then
        dtime_local_half = dtime_local*half
        call advection_step (dtime_local_half)
!        if (advance_source)    call source_step    (dtime_local_half)
        call diffusion_step (dtime_local)
!        if (advance_source)    call source_step    (dtime_local_half)
        call advection_step (dtime_local_half)
      else
        if (advance_flow)      call advection_step (dtime_local)
!        if (advance_source)    call source_step    (dtime_local)
        if (advance_diffusion) call diffusion_step (dtime_local)
      end if
!
      wtime_flux_transport = wtime_flux_transport + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine advection_step (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with advection by the time-step DT.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use input_parameters
      use timing
      use globals
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(IN) :: dtime_local
      real(r_typ) :: t1
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
      dtime_advection_used = dtime_local
!
      if (flow_num_method.eq.1) then
        call advection_step_fe_upwind (dtime_local)
      elseif (flow_num_method.eq.2) then
        call advection_step_rk3tvd_upwind (dtime_local)
      elseif (flow_num_method.eq.3) then
        call advection_step_rk3tvd_weno3 (dtime_local)
      end if
!
      wtime_flux_transport_advection = wtime_flux_transport_advection &
                                       + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine diffusion_step (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with diffusion by the time-step DTIME_LOCAL.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use input_parameters
      use timing
      use globals, ONLY : dtime_diffusion_used,time
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(IN) :: dtime_local
      real(r_typ) :: dtime_local2,t1
      real(r_typ) :: time_stepped,timetmp
      integer :: i
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
      time_stepped = 0.
!
      dtime_local2 = dtime_local/diffusion_subcycles
      dtime_diffusion_used = dtime_local2
!
      do i=1,diffusion_subcycles
!
        if (diffusion_num_method.eq.1) then
!
          call diffusion_step_euler_cd (dtime_local2)
!
        else if (diffusion_num_method.eq.2.or.   &
                 diffusion_num_method.eq.3) then
!
          call diffusion_step_sts_cd (dtime_local2)
!
        end if
!
        time_stepped = time_stepped + dtime_local2
!
        if (verbose.and.MOD(i,1).eq.0) then
          write(*,*) '-->Diff subcycle #',i,' of ', &
                     diffusion_subcycles, 'time:',time_stepped
          flush(OUTPUT_UNIT)
          timetmp = time
          time = time + time_stepped
          call analysis_step
! RMC: Eventaully, all ranks should write out histories for all realizations.
          if (iamp0) call output_histories
          time = timetmp
        end if
!
      enddo
!
      wtime_flux_transport_diffusion = wtime_flux_transport_diffusion &
                                       + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine source_step (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with diffusion by the time-step DT.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use input_parameters
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(IN) :: dtime_local
      real(r_typ) :: t1
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
! ---> Add random flux, ARs, etc.
!
      wtime_source = wtime_source + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine timestamp ()
!
! ****** TIMESTAMP prints out the current YMDHMS date as a timestamp.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    12 January 2007
!  Author:
!    John Burkardt
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer d,h,m,mm,n,s,y
      character*(8) ampm
      character*(8) date
      character*(9) month(12)
      character*(10) time
!
!-----------------------------------------------------------------------
!
      save month
!
!-----------------------------------------------------------------------
!
      data month / &
      'January  ', 'February ', 'March    ', 'April    ', &
      'May      ', 'June     ', 'July     ', 'August   ', &
      'September', 'October  ', 'November ', 'December ' /
!
!-----------------------------------------------------------------------
!
      call DATE_AND_TIME (date,time)
!
      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm
!
      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if
!
      write ( *, &
       '(a1,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
       ' ',d, month(m), y, h, ':', n, ':', s, '.', mm, ampm
!
end subroutine
!#######################################################################
function wtime ()
!
!*********************************************************************72
!
!  WTIME returns a reading of the wall clock time.
!
!  Discussion:
!    To get the elapsed wall clock time, call WTIME before and after
!    a given operation, and subtract the first reading from the second.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Author:
!    John Burkardt
!  Parameters:
!    Output, r_typ WTIME, the wall clock reading, in seconds.
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer*8 :: clock_max
      integer*8 :: clock_rate
      integer*8 :: clock_reading
!
      real(r_typ) :: wtime
!
!-----------------------------------------------------------------------
!
      call SYSTEM_CLOCK (clock_reading,clock_rate,clock_max)
!
      wtime = real(clock_reading,r_typ)/real(clock_rate,r_typ)
!
      return
end function
!#######################################################################
subroutine set_periodic_bc_1d (a,n1)
!
!-----------------------------------------------------------------------
!
! ****** Set the periodic phi direction boundary condition for
! ****** a 1D array assuming a two-point overlap.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(n1) :: a
      integer :: n1
!
!-----------------------------------------------------------------------
!
      a(1)  = a(n1-1)
      a(n1) = a(2)
!
end subroutine
!#######################################################################
subroutine set_periodic_bc_2d (a,n1,n2)
!
!-----------------------------------------------------------------------
!
! ****** Set the periodic phi direction boundary condition for
! ****** a 2D array assuming a two-point overlap.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(n1,n2) :: a
      integer :: n1,n2,j
!
!-----------------------------------------------------------------------
!
      do concurrent (j=1:n1)
        a(j, 1) = a(j,n2-1)
        a(j,n2) = a(j,2)
      enddo
!
end subroutine
!#######################################################################
subroutine set_periodic_bc_3d (a,n1,n2,n3)
!
!-----------------------------------------------------------------------
!
! ****** Set the periodic phi direction boundary condition for
! ****** a 3D array assuming a two-point overlap.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(n1,n2,n3) :: a
      integer :: n1,n2,n3,j,i
!
!-----------------------------------------------------------------------
!
      do concurrent (j=1:n1,i=1:n3)
        a(j, 1,i) = a(j,n2-1,i)
        a(j,n2,i) = a(j,2   ,i)
      enddo
!
end subroutine
!#######################################################################
subroutine load_diffusion_matrix
!
!-----------------------------------------------------------------------
!
! ****** Load diffusion matrix coefs.  Currently, these never change
! ****** over the run, so they are set once.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields
      use input_parameters
      use matrix_storage
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: j,k,i
!
!-----------------------------------------------------------------------
!
! ****** Allocate coef array.
!
      allocate (coef(2:ntm-1,2:npm-1,5,nr))
      coef(:,:,:,:)=0.
!$acc enter data copyin(coef)
!
! ****** Set coef for internal points and phi boundary point at k=1.
!
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1)
!
        coef(j,k,1,i) = half*(diffusion_coef(j  ,k  ,i)  &
                            + diffusion_coef(j+1,k  ,i)) &
                       *dph_i(k  )*dp_i(k)*st_i(j)*st_i(j)
!
        coef(j,k,2,i) = half*(diffusion_coef(j  ,k  ,i)  &
                            + diffusion_coef(j  ,k+1,i)) &
                       *dth_i(j  )*dt_i(j)*st_i(j)*sth(j  )
!
        coef(j,k,4,i) = half*(diffusion_coef(j+1,k  ,i)  &
                            + diffusion_coef(j+1,k+1,i)) &
                       *dth_i(j+1)*dt_i(j)*st_i(j)*sth(j+1)
!
        coef(j,k,5,i) = half*(diffusion_coef(j  ,k+1,i)  &
                            + diffusion_coef(j+1,k+1,i)) &
                       *dph_i(k+1)*dp_i(k)*st_i(j)*st_i(j)
!
        coef(j,k,3,i)=-(coef(j,k,1,i)+coef(j,k,2,i)+coef(j,k,4,i)+coef(j,k,5,i))
!
      enddo
!
end subroutine
!#######################################################################
subroutine get_dtime_diffusion_euler (dtime_exp)
!
!-----------------------------------------------------------------------
!
! ****** Get the explicit Euler time step limits for diffusion.
!
!-----------------------------------------------------------------------
!
      use matrix_storage
      use mpidefs
      use mesh
      use constants
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dtime_exp
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,d,ierr
      real(r_typ) :: max_eig,gersh_rad
!
!-----------------------------------------------------------------------
!
! *** Estimate maximum eigenvalue of all realizations using Gershgorin disks:
!
      max_eig = 0.
!
!$omp parallel do   collapse(3) default(shared)  reduction(max:max_eig)
!$acc parallel loop collapse(3) default(present) reduction(max:max_eig)
      do i=1,nr
        do k=2,npm-1
          do j=2,ntm-1
              gersh_rad = 0.
!$acc loop seq
              do d=1,5
                gersh_rad = gersh_rad+abs(coef(j,k,d,i))
              enddo
              max_eig = max(gersh_rad,max_eig)
          enddo
        enddo
      enddo
!$omp end parallel do
!
! *** Compute the Euler time-step bound.
!
      dtime_exp = two/max_eig
!
! *** Apply safety factor.
!
      dtime_exp = safety*dtime_exp
!
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Allreduce (MPI_IN_PLACE,dtime_exp,1,ntype_real,  &
                          MPI_MIN,MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
end subroutine
!#######################################################################
subroutine diffusion_step_sts_cd (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Diffuse the field by one time step using super time-stepping.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields
      use input_parameters
      use sts
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,el
      integer*8 :: k
      real(r_typ), INTENT(IN) :: dtime_local
!
!-----------------------------------------------------------------------
!
!     This only needs to happen more than once if dtime_local changes.
!
      if (need_to_load_sts) then
        if (diffusion_num_method.eq.2) then
          call load_sts_rkl2 (dtime_local)
        elseif(diffusion_num_method.eq.3) then
          call load_sts_rkg2 (dtime_local)
        end if
        need_to_load_sts = .false.
      end if
!
! ****** Allocate scratch arrays.
!
      allocate (u0(ntm,npm,nr))
      allocate (dty0(ntm,npm,nr))
      allocate (ykm1(ntm,npm,nr))
      allocate (ukm1(ntm,npm,nr))
      allocate (ukm2(ntm,npm,nr))
!$acc enter data create(u0,dty0,ykm1,ukm1,ukm2)
!
      call ax(f,dty0)
!
      do concurrent (el=1:nr,j=1:npm,i=1:ntm)
        u0(i,j,el) = f(i,j,el)
        ukm2(i,j,el) = f(i,j,el)
        dty0(i,j,el) = dtime_local*dty0(i,j,el)
        ukm1(i,j,el) = f(i,j,el) + sts_ubj(1)*dty0(i,j,el)
      enddo
!
! ****** Inner s-step loop
!
      do k=2,sts_s
!
        call ax(ukm1,ykm1)
!
        do concurrent (el=1:nr,j=1:npm,i=1:ntm)
          f(i,j,el) =             sts_uj(k)*ukm1(i,j,el) + &
                                  sts_vj(k)*ukm2(i,j,el) + &
                   (one-sts_uj(k)-sts_vj(k))*u0 (i,j,el) + &
                    sts_ubj(k)*dtime_local*ykm1 (i,j,el) + &
                                 sts_gj(k)*dty0 (i,j,el)
!
          ukm2(i,j,el) = ukm1(i,j,el)
          ukm1(i,j,el) = f(i,j,el)
        enddo
!
      enddo
!
!$acc exit data delete(u0,dty0,ykm1,ukm1,ukm2)
      deallocate (u0,dty0,ykm1,ukm1,ukm2)
!
end subroutine
!#######################################################################
subroutine load_sts_rkl2 (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Set up parameters and coefficient arrays for STS advance.
! ****** This uses the RKL2 2nd-order STS as given in
! ****** Meyer, et. al. J. Comp. Phys. 257 (2014) 594-626
!
! ****** The Euler stable dt (dtime_diffusion_euler) must be already be set.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use mesh
      use sts
      use constants
      use globals
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dtime_local
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: sts_s_real,bj_bjm2,bj_bjm1,w
      integer*8 :: j
      logical, save :: first = .true.
!
!-----------------------------------------------------------------------
!
! ****** Compute number of iterations per super-step.
!
      sts_s_real = &
        half*(sqrt(nine + sixteen*(dtime_local/dtime_diffusion_euler)) &
              - one)
!
      sts_s = CEILING(sts_s_real)
!
! ****** Make sure s is at least 3.
!
      if (sts_s.lt.3) then
        sts_s = 3
      end if
!
! ****** Make sure s is odd.
!
      if (MOD(sts_s,2).eq.0) then
        sts_s = sts_s + 1
      end if
!
! ****** Allocate super-time-step coefficent arrays.
!
      if (need_to_load_sts.and..not.first) then
!$acc exit data delete(sts_uj,sts_vj,sts_ubj,sts_gj)
        deallocate (sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
      end if
!
      if (need_to_load_sts) then
        allocate (sts_uj(sts_s))
        allocate (sts_vj(sts_s))
        allocate (sts_ubj(sts_s))
        allocate (sts_gj(sts_s))
        allocate (sts_b(sts_s))
!$acc enter data create(sts_uj,sts_vj,sts_ubj,sts_gj)
        if (first) first=.false.
      end if
!
! ****** Compute super-time-step coefficents.
!
      w = four/(sts_s*sts_s + sts_s - two_int)
!
      !b0 = one/three
      sts_b(1) = one/three
      sts_b(2) = one/three
!
      sts_uj(1) = -9999._r_typ
      sts_vj(1) = -9999._r_typ
      sts_ubj(1) = four/(three*(sts_s*sts_s + sts_s - two_int))
      sts_gj(1) = -9999._r_typ
!
      sts_uj(2) = three/two
      sts_vj(2) = -half
      sts_ubj(2) = sts_uj(2)*w
      sts_gj(2) = -w
!
      do j=3,sts_s
        sts_b(j) = (j*j + j - two)/(two*j*(j+1))
        bj_bjm1 = sts_b(j)/sts_b(j-1)
        bj_bjm2 = sts_b(j)/sts_b(j-2)
!
        sts_uj(j) = bj_bjm1*(two - one/j)
        sts_vj(j) = bj_bjm2*(one/j - one)
        sts_ubj(j) = sts_uj(j)*w
        sts_gj(j) = (sts_b(j-1) - one)*sts_ubj(j)
      enddo
!
!$acc update device(sts_uj,sts_vj,sts_ubj,sts_gj)
!
      if (verbose) then
        write (*,*) '*** Diffusion: 2nd-order RKL2  ***'
        write (*,*) 'Super-time-step used              = ',dtime_local
        write (*,*) 'Euler time-step                   = ',dtime_diffusion_euler
        write (*,*) 'Number of STS iterations needed   = ',sts_s
        write (*,*) 'Number of Euler iterations needed = ',dtime_local/dtime_diffusion_euler
        write (*,*) 'Potential max speedup of STS      = ', &
                          (dtime_local/dtime_diffusion_euler)/sts_s,'X'
      end if
!
end subroutine
!#######################################################################
subroutine load_sts_rkg2 (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Set up parameters and coefficient arrays for STS advance.
! ****** This uses the RKG2 2nd-order STS as given in
! ****** Skaras, et. al. J. Comp. Phys. 425 (2021) 109879
!
! ****** The Euler stable dt (dtime_diffusion_euler) must be already be set.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use mesh
      use sts
      use constants
      use globals
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dtime_local
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: sts_s_real,bj_bjm2,bj_bjm1,w
      integer*8 :: j
      logical, save :: first = .true.
!
!-----------------------------------------------------------------------
!
! ****** Compute number of iterations per super-step.
!
      sts_s_real = half*(sqrt(twentyfive + &
                         twentyfour*(dtime_local/dtime_diffusion_euler)) &
                         - three)
!
      sts_s = CEILING(sts_s_real)
!
! ****** Make sure s is at least 3.
!
      if (sts_s.lt.3) then
        sts_s = 3
      end if
!
! ****** Make sure s is odd.
!
      if (MOD(sts_s,2).eq.0) then
        sts_s = sts_s + 1
      end if
!
! ****** Allocate super-time-step coefficent arrays.
!
      if (need_to_load_sts.and..not.first) then
!$acc exit data delete(sts_uj,sts_vj,sts_ubj,sts_gj)
        deallocate (sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
      end if
!
      if (need_to_load_sts) then
        allocate (sts_uj(sts_s))
        allocate (sts_vj(sts_s))
        allocate (sts_ubj(sts_s))
        allocate (sts_gj(sts_s))
        allocate (sts_b(sts_s))
!$acc enter data create(sts_uj,sts_vj,sts_ubj,sts_gj)
        if (first) first=.false.
      end if
!
! ****** Compute super-time-step coefficents.
!
      w = six/((sts_s + four_int)*(sts_s - one_int))
!
      !b0 = one
      sts_b(1)=one/three
      sts_b(2)=one/fifteen
!
      sts_uj(1) = one
      sts_vj(1) = -9999._r_typ
      sts_ubj(1) = w
      sts_gj(1) = -9999._r_typ
!
      sts_uj(2) = one/two
      sts_vj(2) = -one/ten
      sts_ubj(2) = sts_uj(2)*w
      sts_gj(2) = zero
!
      do j=3,sts_s
        sts_b(j) = (four*(j-1)*(j+4))/(three*j*(j+1)*(j+2)*(j+3))
        bj_bjm1 = sts_b(j)/sts_b(j-1)
        bj_bjm2 = sts_b(j)/sts_b(j-2)
!
        sts_uj(j) = bj_bjm1*(two + one/j)
        sts_vj(j) = -bj_bjm2*(one/j + one)
        sts_ubj(j) = sts_uj(j)*w
        sts_gj(j) = (half*j*(j + one)*sts_b(j-1) - one)*sts_ubj(j)
      enddo
!$acc update device(sts_uj,sts_vj,sts_ubj,sts_gj)
!
      if (verbose) then
        write (*,*) '*** Diffusion: 2nd-order RKG2  ***'
        write (*,*) 'Super-time-step used              = ',dtime_local
        write (*,*) 'Euler time-step                   = ',dtime_diffusion_euler
        write (*,*) 'Number of STS iterations needed   = ',sts_s
        write (*,*) 'Number of Euler iterations needed = ',dtime_local/dtime_diffusion_euler
        write (*,*) 'Potential max speedup of STS      = ', &
                          (dtime_local/dtime_diffusion_euler)/sts_s,'X'
      end if
!
end subroutine
!#######################################################################
subroutine read_2d_file (fname,ln1,ln2,fin,s1,s2,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read 2D data from H5 file FNAME.
! ****** This allocates arrays "fin", "s1" and "s2"
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use ds_def
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      real(r_typ), dimension(:,:), allocatable :: fin
      real(r_typ), dimension(:), allocatable :: s1,s2
      integer :: ln1,ln2
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      integer :: ierr
      real(r_typ) :: t1
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
! ****** Read the data.
!
      call rdh5 (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_2D_FILE:'
        write (*,*) '### Could not read the requested data set.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(fname)
        call endrun(.true.)
      end if
!
! ****** Check that it is a 2D dataset.
!
      if (s%ndim.ne.2) then
        write (*,*)
        write (*,*) '### ERROR in READ_2D_FILE:'
        write (*,*) '### The file does not contain a 2D field.'
        write (*,*) '### Number of dimensions = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        call endrun(.true.)
      end if
!
! ****** Check that the their are scales.
!
      if (.not.s%scale) then
        write (*,*)
        write (*,*) '### ERROR in READ_2D_FILE:'
        write (*,*) '### The data set does not contain scales.'
        write (*,*) 'File name: ',trim(fname)
        call endrun(.true.)
      end if
!
! ****** Get the resolution of the input map and set resolution for
!        returned arrays.
!
      ln1 = s%dims(1)
      ln2 = s%dims(2)
!
! ****** Allocate and load the scales and data.
!
      allocate (s1(ln1))
      allocate (s2(ln2))
      allocate (fin(ln1,ln2))
!
      s1(:) = s%scales(1)%f(:)
      s2(:) = s%scales(2)%f(:)
!
! ****** Read the data.
!
      fin(:,:) = s%f(:,:,1)
!
! ****** Free up memory.
!
      deallocate (s%scales(1)%f)
      deallocate (s%scales(2)%f)
      deallocate (s%f)
!
      wtime_io = wtime_io + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine read_3d_file (fname,ln1,ln2,ln3,fin,s1,s2,s3,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read 3D data from H5 file FNAME.
! ****** This allocates arrays "fin", "s1", "s2", and "s3"
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use ds_def
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      real(r_typ), dimension(:,:,:), allocatable :: fin
      real(r_typ), dimension(:), allocatable :: s1,s2,s3
      integer :: ln1,ln2,ln3
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      integer :: ierr
      real(r_typ) :: t1
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
! ****** Read the data.
!
      call rdh5 (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_3D_FILE:'
        write (*,*) '### Could not read the requested data set.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(fname)
        call endrun(.true.)
      end if
!
! ****** Check that it is a 3D dataset.
!
      if (s%ndim.ne.3) then
        write (*,*)
        write (*,*) '### ERROR in READ_3D_FILE:'
        write (*,*) '### The file does not contain a 3D field.'
        write (*,*) '### Number of dimensions = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        call endrun(.true.)
      end if
!
! ****** Check that the their are scales.
!
      if (.not.s%scale) then
        write (*,*)
        write (*,*) '### ERROR in READ_3D_FILE:'
        write (*,*) '### The data set does not contain scales.'
        write (*,*) 'File name: ',trim(fname)
        call endrun(.true.)
      end if
!
! ****** Get the resolution of the input map and set resolution for
!        returned arrays.
!
      ln1 = s%dims(1)
      ln2 = s%dims(2)
      ln3 = s%dims(3)
!
! ****** Allocate and load the scales and data.
!
      allocate (s1(ln1))
      allocate (s2(ln2))
      allocate (s3(ln3))
      allocate (fin(ln1,ln2,ln3))
!
      s1(:) = s%scales(1)%f(:)
      s2(:) = s%scales(2)%f(:)
      s3(:) = s%scales(3)%f(:)
!
! ****** Read the data.
!
      fin(:,:,:) = s%f(:,:,:)
!
! ****** Free up memory.
!
      deallocate (s%scales(1)%f)
      deallocate (s%scales(2)%f)
      deallocate (s%scales(3)%f)
      deallocate (s%f)
!
      wtime_io = wtime_io + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine rdh5 (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read a 1D, 2D, or 3D scientific data set from an HDF5 file.
! ****** The HDF5 file is currently assumed to contain only one
! ****** dataset (1D,2D,or 3D), with or without scales, in group "/",
! ****** and has no other data members.
!
!-----------------------------------------------------------------------
!
! ****** This routine allocates the required memory and returns
! ****** pointers to the data and scale arrays.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    File name to read from.
!
! ****** Output arguments:
!
!          S       : [structure of type DS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was read
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
! ****** Components of structure S:
!
!          NDIM    : [integer]
!                    Number of dimensions found in the data set.
!
!          DIMS    : [integer, dimension(3)]
!                    Number of points in the data set dimensions.
!                    For a 1D data set, DIMS(2)=DIMS(3)=1.
!                    For a 2D data set, DIMS(3)=1.
!
!          SCALE   : [logical]
!                    Flag to indicate the presence of scales (axes)
!                    in the data set.  SCALE=.false. means that scales
!                    were not found; SCALE=.true. means that scales
!                    were found.
!
!          HDF32   : [logical]
!                    Flag to indicate the precision of the data set
!                    read in.  HDF32=.true. means that the data is
!                    32-bit; HDF32=.false. means that the data is
!                    64-bit.
!
!          SCALES  : [structure of type RP1D, dimension(3)]
!                    This array holds the pointers to the scales
!                    when SCALE=.true., and is undefined otherwise.
!
!          F       : [real, pointer to a rank-3 array]
!                    This array holds the data set values.
!
! ****** The storage for the arrays pointed to by F, and the
! ****** scales (if present) in structure SCALES, is allocated by
! ****** this routine.
!
!-----------------------------------------------------------------------
!
      use iso_fortran_env
      use number_types
      use ds_def
      use hdf5
      use h5ds
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      character(*) :: fname
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,n_members,obj_type
!
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id     ! Dataspace identifier
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HID_T) :: datatype_id   ! Datatype identifiers
!
      integer(SIZE_T) :: prec
!
      integer(HSIZE_T),dimension(:), allocatable :: s_dims,maxpts
      integer(HSIZE_T),dimension(1) :: s_dims_i
!
      real(REAL32), dimension(:,:,:), allocatable :: f4
      real(REAL32), dimension(:),     allocatable :: f4dim
      real(REAL64), dimension(:,:,:), allocatable :: f8
      real(REAL64), dimension(:),     allocatable :: f8dim
!
      character(512) :: obj_name
      character(4), parameter :: cname='RDH5'
!
      logical :: is_scale
!
!-----------------------------------------------------------------------
!
! ****** Initialize dimension count and arrays.
!
      s%ndim = 0
      s%dims(:) = 1
!
! ****** Initialize hdf5 interface.
!
      call h5open_f (ierr)
!
! ****** Open hdf5 file.
!
      call h5Fopen_f (trim(fname),H5F_ACC_RDONLY_F,file_id,ierr)
!
! ****** Get information about the hdf5 file.
!
      call h5Gn_members_f (file_id,"/",n_members,ierr)
!
! ****** Make sure there is (at maximum) one 3D dataset with scales.
!
      if (n_members.eq.0.or.n_members.gt.4) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file contains too few/many datasets.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
! ****** Assume the Dataset is in index 0 and get its name.
!
      call h5Gget_obj_info_idx_f (file_id,"/",0,obj_name,obj_type,ierr)
!
! ****** Open Dataset.
!
      call h5Dopen_f (file_id,trim(obj_name),dset_id,ierr)
!
! ****** Make sure the Dataset is not a scale.
!
      call h5DSis_scale_f(dset_id,is_scale,ierr)
      if (is_scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file Dataset at index 0 is a scale.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
! ****** Get dimensions (need s_dims array for format requirements).
!
      call h5Dget_space_f (dset_id,dspace_id,ierr)
      call h5Sget_simple_extent_ndims_f (dspace_id,s%ndim,ierr)
!
      allocate(s_dims(s%ndim))
!
      allocate(maxpts(s%ndim))
      call h5Sget_simple_extent_dims_f (dspace_id,s_dims,maxpts,ierr)
      deallocate(maxpts)
!
      do j=1,s%ndim
        s%dims(j) = INT(s_dims(j))
      enddo
!
! ****** Get the floating-point precision of the data and set flag.
!
      call h5Dget_type_f (dset_id,datatype_id,ierr)
      call h5Tget_precision_f (datatype_id,prec,ierr)
!
      if (prec.eq.32) then
        s%hdf32=.true.
      elseif (prec.eq.64) then
        s%hdf32=.false.
      end if
!
! ****** Allocate the memory for the Dataset array in s.
!
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
!
! ****** Need to read the file in its own datatype, and then convert
! ****** to datatype of s%f.
!
      if (s%hdf32) then
        allocate (f4(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f4,s_dims,ierr)
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              s%f(i,j,k) = REAL(f4(i,j,k),r_typ)
            enddo
          enddo
        enddo
        deallocate (f4)
      else
        allocate (f8(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f8,s_dims,ierr)
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              s%f(i,j,k) = REAL(f8(i,j,k),r_typ)
            enddo
          enddo
        enddo
        deallocate (f8)
      end if
!
      deallocate(s_dims)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDH5:'
        write (*,*) '### Error while reading the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from H5DREAD_F) = ',ierr,']'
        ierr=4
        return
      end if
!
! ****** Close the hdf5 type descriptor.
!
      call h5Tclose_f (datatype_id,ierr)
!
! ****** Check if there might be scales present, if so, read them.
!
      if (n_members.gt.1) then
!
! ***** First check that the number of scale datasets match the # dim.
!
        if (n_members-1.ne.s%ndim) then
          write (*,*)
          write (*,*) '### ERROR in RDH5:'
          write (*,*) '### # scales does not match # dims.'
          write (*,*) 'File name: ',trim(fname)
          return
        end if
!
        s%scale=.true.
!
! ****** Loop through scales, make sure each is a scale, and read them.
!
        do i=1,n_members-1
!
! ****** Get the name of scale dataset.
!
          call h5Gget_obj_info_idx_f (file_id,"/",i, &
                                     obj_name,obj_type,ierr)
!
! ****** Open scale dataset.
!
          call h5Dopen_f (file_id,trim(obj_name),dim_id,ierr)
!
! ****** Make sure the scale is a scale.
!
          call h5DSis_scale_f (dim_id,is_scale,ierr)
          if (.not.is_scale) then
            write (*,*)
            write (*,*) '### ERROR in RDH5:'
            write (*,*) '### Scale is not a scale.'
            write (*,*) 'File name: ',trim(fname)
            return
          end if
!
! ****** Get dimension of scale.
!
          s_dims_i = s%dims(i)
!
! ****** Allocate scale.
!
          allocate (s%scales(i)%f(s_dims_i(1)))
!
! ****** Get the floating-point precision of the scale.
!
          call h5Dget_type_f (dim_id,datatype_id,ierr)
          call h5Tget_precision_f (datatype_id,prec,ierr)
!
! ****** Read in the scale data.
!
          if (prec.eq.32) then
            allocate (f4dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f4dim,s_dims_i,ierr)
            do j=1,s%dims(i)
              s%scales(i)%f(j) = REAL(f4dim(j),r_typ)
            enddo
            deallocate (f4dim)
          elseif (prec.eq.64) then
            allocate (f8dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f8dim,s_dims_i,ierr)
            do j=1,s%dims(i)
              s%scales(i)%f(j) = REAL(f8dim(j),r_typ)
            enddo
            deallocate (f8dim)
          end if
!
! ****** Close the scale dataset.
!
          call h5Dclose_f (dim_id,ierr)
!
        enddo
!
! ****** Allocate dummy scales (of length 1) for empty dimensions.
!
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      else
!
! ****** If scales are not present, allocate dummy
! ****** scales (of length 1) so that the pointers to the scales
! ****** are valid.
!
        s%scale = .false.
!
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      end if
!
! ****** Close the dataset.
!
      call h5Dclose_f (dset_id,ierr)
!
! ****** Close the file.
!
      call h5Fclose_f (file_id,ierr)
!
! ****** Close FORTRAN interface.
!
      call h5close_f (ierr)
!
end subroutine
!#######################################################################
subroutine wrh5 (fname,s,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write a 1D, 2D, or 3D scientific data set to an HDF5 file.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          FNAME   : [character(*)]
!                    File name to write to.
!
!          S       : [structure of type DS]
!                    A structure that holds the field, its
!                    dimensions, and the scales, with the
!                    components described below.
!
! ****** Output arguments:
!
!          IERR    : [integer]
!                    IERR=0 is returned if the data set was written
!                    successfully.  Otherwise, IERR is set to a
!                    nonzero value.
!
!-----------------------------------------------------------------------
!
      use iso_fortran_env
      use number_types
      use ds_def
      use hdf5
      use h5ds
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      type(ds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      character(8) ::   dimname
      integer :: i,j,k
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id,dspacedim_id   ! Dataspace identifiers
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HSIZE_T),dimension(3) :: s_dims
      integer(HSIZE_T),dimension(1) :: s_dims_i
!
      real(REAL32), dimension(:,:,:), allocatable :: f4
      real(REAL32), dimension(:),     allocatable :: f4dim
      real(REAL64), dimension(:,:,:), allocatable :: f8
      real(REAL64), dimension(:),     allocatable :: f8dim
!
!-----------------------------------------------------------------------
!
! ****** HDF5 calls are picky about the integer format for the dims
! ****** so the s%dims need to be converted to HSIZE_T integers.
!
! ****** Also, sometimes calls to wrhdf() for 1D and 2D datasets
! ****** do not have the unused dims(i) set to 1 (unset).
! ****** To avoid needing a function call to implicitly reshape
! ****** f(n), set the dims here.
!
      do i=1,3
         if (i.le.s%ndim) then
           s_dims(i) = INT(s%dims(i),HSIZE_T)
         else
           s_dims(i) = 1
         end if
      enddo
!
! ****** Initialize hdf5 interface.
!
      call h5open_f (ierr)
!
! ****** Create the file.
!
      call h5Fcreate_f (trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
!
! ****** Create the dataspace.
!
      call h5Screate_simple_f (s%ndim,s_dims,dspace_id,ierr)
!
! ****** Create and write the dataset (convert s%f to proper type).
!
      if (s%hdf32) then
        allocate (f4(s_dims(1),s_dims(2),s_dims(3)))
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              f4(i,j,k) = REAL(s%f(i,j,k),REAL32)
            enddo
          enddo
        enddo
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_REAL, &
                          dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_REAL,f4,s_dims,ierr)
        deallocate (f4)
      else
        allocate (f8(s_dims(1),s_dims(2),s_dims(3)))
        do k=1,s%dims(3)
          do j=1,s%dims(2)
            do i=1,s%dims(1)
              f8(i,j,k) = REAL(s%f(i,j,k),REAL64)
            enddo
          enddo
        enddo
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_DOUBLE, &
                          dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_DOUBLE,f8,s_dims,ierr)
        deallocate (f8)
      end if
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRH5:'
        write (*,*) '### Could not write the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from h5Dwrite_f) = ',ierr,']'
        ierr=4
        return
      end if
!
! ****** Check for scales.  If present, add them to the hdf5 dataset.
!
      if (s%scale) then
        do i=1,s%ndim
          if (i.eq.1) then
            dimname='dim1'
          elseif (i.eq.2) then
            dimname='dim2'
          elseif (i.eq.3) then
            dimname='dim3'
          end if
          s_dims_i = s_dims(i)
          call h5Screate_simple_f(1,s_dims_i,dspacedim_id,ierr)
          if (s%hdf32) then
            allocate (f4dim(s_dims_i(1)))
            do j=1,s%dims(i)
              f4dim(j) = REAL(s%scales(i)%f(j),REAL32)
            enddo
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_REAL, &
                              dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_REAL, &
                             f4dim,s_dims_i,ierr)
            deallocate (f4dim)
          else
            allocate (f8dim(s_dims_i(1)))
            do j=1,s%dims(i)
              f8dim(j) = REAL(s%scales(i)%f(j),REAL64)
            enddo
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_DOUBLE, &
                             dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_DOUBLE, &
                             f8dim,s_dims_i,ierr)
            deallocate (f8dim)
          end if
          call h5DSset_scale_f (dim_id,ierr,dimname)
          call h5DSattach_scale_f (dset_id,dim_id,i,ierr)
          call h5DSset_label_f(dset_id, i, dimname, ierr)
          call h5Dclose_f (dim_id,ierr)
          call h5Sclose_f (dspacedim_id,ierr)
        enddo
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in WRH5:'
          write (*,*) '### Could not write the scales.'
          write (*,*) 'File name: ',trim(fname)
          ierr = 5
          return
        end if
      end if
!
! ****** Close the dataset.
!
      call h5Dclose_f (dset_id,ierr)
!
! ****** Close the dataspace.
!
      call h5Sclose_f (dspace_id,ierr)
!
! ****** Close the file.
!
      call h5Fclose_f (file_id,ierr)
!
! ****** Close the hdf5 interface.
!
      call h5close_f (ierr)
!
end subroutine
!#######################################################################
subroutine write_2d_file (fname,ln1,ln2,f,s1,s2,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Write 2D data to H5 file FNAME.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use ds_def
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      real(r_typ), dimension(ln1,ln2) :: f
      real(r_typ), dimension(ln1) :: s1
      real(r_typ), dimension(ln2) :: s2
      integer :: ln1,ln2
      real(r_typ) :: t1
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
! ****** Set the structure components.
!
      s%ndim = 2
      s%dims(1) = ln1
      s%dims(2) = ln2
      s%dims(3) = 1
      s%scale = .true.
      s%hdf32 = .true.
!
      allocate (s%scales(1)%f(ln1))
      allocate (s%scales(2)%f(ln2))
      allocate (s%f(ln1,ln2,1))
!
      s%scales(1)%f(:) = s1(:)
      s%scales(2)%f(:) = s2(:)
      s%f(:,:,1) = f(:,:)
!
! ****** Write the data set.
!
      call wrh5 (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_2D_FILE:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
! ****** Free up memory.
!
      deallocate (s%scales(1)%f)
      deallocate (s%scales(2)%f)
      deallocate (s%f)
!
      wtime_io = wtime_io + (MPI_Wtime() - t1)
!
end subroutine
!#######################################################################
subroutine write_3d_file (fname,ln1,ln2,ln3,f,s1,s2,s3,ierr)
!
!-----------------------------------------------------------------------
! ****** Write 2D data to H5 file FNAME.
!-----------------------------------------------------------------------
!
      use number_types
      use ds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      real(r_typ), dimension(ln1,ln2,ln3) :: f
      real(r_typ), dimension(ln1) :: s1
      real(r_typ), dimension(ln2) :: s2
      real(r_typ), dimension(ln3) :: s3
      integer :: ln1,ln2,ln3
      type(ds) :: s
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Set the structure components.
!
      s%ndim = 3
      s%dims(1) = ln1
      s%dims(2) = ln2
      s%dims(3) = ln3
      s%scale = .true.
      s%hdf32 = .true.
!
      allocate (s%scales(1)%f(ln1))
      allocate (s%scales(2)%f(ln2))
      allocate (s%scales(3)%f(ln3))
      allocate (s%f(ln1,ln2,ln3))
!
      s%scales(1)%f(:) = s1(:)
      s%scales(2)%f(:) = s2(:)
      s%scales(3)%f(:) = s3(:)
      s%f(:,:,:) = f(:,:,:)
!
! ****** Write the data set.
!
      call wrh5 (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_3D_FILE:'
        write (*,*) '### Could not write the 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
!
! ****** Free up memory.
!
      deallocate (s%scales(1)%f)
      deallocate (s%scales(2)%f)
      deallocate (s%scales(3)%f)
      deallocate (s%f)
!
end subroutine
!#######################################################################
subroutine set_mesh
!
!-----------------------------------------------------------------------
!
! ****** Compute the mesh quantities.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use constants
      use mpidefs
      use output, ONLY : pout
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Tolerance for precision of coordinates.
!
      real(r_typ), parameter :: eps=1.e-6_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
! ****** Set mesh resolutions.
!
      nt = ntm1 + 1
      np = npm1 + 1
      ntm2 = ntm1 - 1
      npm2 = npm1 - 1
      flow_t_size = npm1*nt
      flow_p_size = np*ntm1
!
! ****** Set "main mesh" compute resolution.
! ****** The phi direction has an extra point for a two-point overlap.
!
      ntm = ntm1
      npm = np
!
! ****** Check that the mesh covers a complete spherical surface.
!
      if (abs(t(1)).gt.eps.or. &
          abs(t(ntm)-pi).gt.eps.or. &
          abs(p(1)).gt.eps.or. &
          abs(p(npm-1)-twopi).gt.eps) then
        write (*,*)
        write (*,*) '### ERROR in SET_MESH:'
        write (*,*) '### Anomaly in data file coordinates:'
        write (*,*)
        write (*,*) 'Expected t range:'
        write (*,*) 'Min: ',0.
        write (*,*) 'Max: ',pi
        write (*,*) 'Actual t range:'
        write (*,*) 'Min: ',t(1)
        write (*,*) 'Max: ',t(ntm)
        write (*,*)
        write (*,*) 'Expected p range:'
        write (*,*) 'Min: ',0.
        write (*,*) 'Max: ',twopi
        write (*,*) 'Actual p range:'
        write (*,*) 'Min: ',p(1)
        write (*,*) 'Max: ',p(npm-1)
        call endrun(.true.)
      end if
!
! ****** Allocate main mesh grid quantities.
!
      allocate (dt  (ntm))
      allocate (dt_i(ntm))
      allocate (st  (ntm))
      allocate (st_i(ntm))
      allocate (ct  (ntm))
      allocate (dp  (npm))
      allocate (dp_i(npm))
!
! ****** Allocate half mesh grid quantities.
!
      allocate (th   (nt))
      allocate (ph   (np))
      allocate (dth  (nt))
      allocate (dth_i(nt))
      allocate (sth  (nt))
      allocate (sth_i(nt))
      allocate (cth  (nt))
      allocate (dph  (np))
      allocate (dph_i(np))
!
! ****** Set theta grids.
!
      do i=2,ntm1
        th(i) = half*(t(i) + t(i-1))
        dth(i) = t(i) - t(i-1)
      enddo
      th(1) = th(2) - dth(2)
      th(nt) = th(ntm1) + dth(ntm1)
      dth(1) = dth(2)
      dth(nt) = dth(ntm1)
      dth_i(:) = one/dth(:)
!
      do i=1,ntm
        dt(i) = th(i+1) - th(i)
        !   half*(t(i+1) + t(i)) - half*(t(i) + t(i-1))
        ! = half*(t(i+1) - t(i-1))
        dt_i(i) = one/dt(i)
      enddo
!
      st(:) = sin(t(:))
      st(1) = 0.
      st(ntm) = 0.
!
      ct(:) = cos(t(:))
      ct(1) = one
      ct(ntm) = one
!
      st_i(2:ntm-1) = one/st(2:ntm-1)
      st_i(1) = 0.
      st_i(ntm) = 0.
!
      sth(:) = sin(th(:))
      sth_i(2:ntm1) = one/sth(2:ntm1)
      sth_i(1) = 0.
      sth_i(nt) = 0.
!
      cth(:) = cos(th(:))
!
! ****** Set phi grids.
!
      do i=2,np
        ph(i) = half*(p(i) + p(i-1))
        dph(i) = p(i) - p(i-1)
      enddo
      ph(1) = ph(npm1) - twopi
      dph(1) = dph(npm1)
!
      dph_i(:)=one/dph(:)
!
      do i=1,npm-1
        dp(i) = ph(i+1) - ph(i)
      enddo
      dp(npm) = dp(2)
!
      dp_i(:) = one/dp(:)
!
  ! ****** Set up output scale for phi (single point overlap)
      allocate (pout(npm-1))
      pout(:) = p(1:npm-1)
!
!$acc enter data copyin (t,p,dt,dt_i,st,st_i,ct,dp,dp_i,th,ph, &
!$acc&                   dth,dth_i,sth,sth_i,cth,dph,dph_i)
end subroutine
!#######################################################################
subroutine load_diffusion
!
!-----------------------------------------------------------------------
!
! ****** Define the diffusion coef array.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use mpidefs
      use fields
      use input_parameters
      use constants
      use read_2d_file_interface
      use globals
      use sts
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ierr
      real(r_typ), dimension(:), allocatable :: fn1,fs1
!
      integer :: nft,nfp
      real(r_typ), dimension(:), allocatable :: tf,pf
      real(r_typ), dimension(:,:), allocatable :: f_tmp2d
      real(r_typ), dimension(:,:,:), allocatable :: vf,diffusion_coef_file
!
!-----------------------------------------------------------------------
!
! ****** Allocate memory for the total diffusion coef on the half mesh.
!
      allocate (diffusion_coef(nt,np,nr))
!
! ****** Set the initial diffusion coef to the uniform value.
!
      diffusion_coef(:,:,:) = diffusion_coef_constant
!
! ****** Read the diffusion coef file, if it was specified.
!
      if (diffusion_coef_filename.ne.' ') then
!
! ****** Load the diffusion coef into an array (assumes PT).
!
        if (iamp0) then
          call read_2d_file (diffusion_coef_filename,nfp,nft,f_tmp2d,pf,tf,ierr)
        endif
        wtime_tmp_mpi = MPI_Wtime()
        call MPI_Bcast (nfp,1,ntype_real,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (nft,1,ntype_real,0,MPI_COMM_WORLD,ierr)
        if (.NOT. iamp0) then
          allocate (pf(nfp))
          allocate (tf(nft))
          allocate (f_tmp2d(nfp,nft))
        endif
        call MPI_Bcast (f_tmp2d,nfp*nft,ntype_real,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (pf,nfp,ntype_real,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (tf,nft,ntype_real,0,MPI_COMM_WORLD,ierr)
        wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
        f_tmp2d(:,:) = TRANSPOSE(f_tmp2d(:,:))
        allocate(vf(nfp,nft,nr))
        do i=1,nr
          vf(:,:,i) = f_tmp2d(:,:)
        enddo
        deallocate(f_tmp2d)
!
! ****** Interpolate the diffusion coef onto the half mesh (th,ph).
!
        allocate (diffusion_coef_file(nt,np,nr))
!
        do i=1,nr
          call interp2d (nft,nfp,tf,pf,vf(:,:,i),nt,np,th,ph, &
                         diffusion_coef_file(:,:,i),ierr)
        enddo
!
! ****** Error exit: interpolation error.
!
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in LOAD_DIFFUSION_COEF:'
          write (*,*) '### The scales in the diffusion coef file'// &
                 ' are not monotonically increasing.'
          write (*,*) 'File name: ',trim(diffusion_coef_filename)
          call endrun(.true.)
        end if
!
! ****** Enforce periodicity.
!
        diffusion_coef_file(:, 1,:) = diffusion_coef_file(:,npm1,:)
        diffusion_coef_file(:,np,:) = diffusion_coef_file(:,   2,:)
!
! ****** Set the pole value to only have an m=0 component.
!
        allocate (fn1(nr))
        allocate (fs1(nr))
!
        fn1(:) = 0.
        fs1(:) = 0.
        do i=1,nr
          do k=2,npm1
            fn1(i) = fn1(i) + diffusion_coef_file(2,k,i)*dph(k)
            fs1(i) = fs1(i) + diffusion_coef_file(ntm1,k,i)*dph(k)
          enddo

          fn1(i) = fn1(i)*twopi_i
          fs1(i) = fs1(i)*twopi_i
!
          diffusion_coef_file(1 ,:,i) = two*fn1(i) &
                                      - diffusion_coef_file(   2,:,i)
          diffusion_coef_file(nt,:,i) = two*fs1(i) &
                                      - diffusion_coef_file(ntm1,:,i)
        enddo
!
        if (verbose) then
          write (*,*)
          write (*,*) 'Diffusion coef from file: ', &
                                           trim(diffusion_coef_filename)
          write (*,*) 'Minimum value = ',minval(diffusion_coef_file)
          write (*,*) 'Maximum value = ',maxval(diffusion_coef_file)
        end if
!
! ****** Add the file diffusion coef to the uniform value.
!
        diffusion_coef(:,:,:) = diffusion_coef(:,:,:) &
                              + diffusion_coef_file(:,:,:)
!
        deallocate (diffusion_coef_file)
        deallocate (vf)
        deallocate (tf)
        deallocate (pf)
        deallocate (fn1)
        deallocate (fs1)
!
      end if
!
! ****** Add grid-based diffusion coef if requested.
!
      if (diffusion_coef_grid) then
        do i=1,nr
          do k=1,np
            do j=1,nt
              diffusion_coef(j,k,i) = diffusion_coef(j,k,i) &
                                 + (dth(j)**2 + (dph(k)*sth(j))**2)
            enddo
          enddo
        enddo

        if (verbose) then
          write (*,*)
          write (*,*) 'Grid-based diffusion coef is activated.'
        end if

      end if
!
      diffusion_coef(:,:,:) = diffusion_coef_factor*diffusion_coef(:,:,:)
!$acc enter data copyin(diffusion_coef)
!
      call load_diffusion_matrix
!
! ****** Get stable explicit Euler timestep.
!
      call get_dtime_diffusion_euler (dtime_diffusion_euler)
!
end subroutine
!#######################################################################
subroutine load_weno
!
!-----------------------------------------------------------------------
!
! ****** Setup WENO3.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use weno
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: j,k
      real(r_typ) :: dt_total,dp_total
!
!-----------------------------------------------------------------------
!
! ****** Allocate arrays.
!
      allocate (alpha_t(ntm,npm,nr))
      alpha_t(:,:,:) = 0.
      allocate (alpha_p(ntm,npm,nr))
      alpha_p(:,:,:) = 0.
!
      allocate (D_C_CPt(nt))
      allocate (D_C_MCt(nt))
      allocate (D_M_Tt(nt))
      allocate (D_CP_Tt(nt))
      allocate (D_P_Tt(nt))
      allocate (D_MC_Tt(nt))
!
      allocate (D_C_CPp(np))
      allocate (D_C_MCp(np))
      allocate (D_M_Tp(np))
      allocate (D_CP_Tp(np))
      allocate (D_P_Tp(np))
      allocate (D_MC_Tp(np))
!
! ****** Set grid weights.
!
      do j=2,ntm1
        dt_total = dt(j-1) + dt(j) + dt(j+1)
        D_C_CPt(j) =           dt(j  )/(dt(j  )+dt(j+1))
        D_C_MCt(j) =           dt(j  )/(dt(j-1)+dt(j  ))
        D_M_Tt (j) =           dt(j-1)/dt_total
        D_CP_Tt(j) = (dt(j  )+dt(j+1))/dt_total
        D_P_Tt (j) =           dt(j+1)/dt_total
        D_MC_Tt(j) = (dt(j-1)+dt(j  ))/dt_total
      enddo
!
      do k=2,npm1
        dp_total = dp(k-1) + dp(k) + dp(k+1)
        D_C_CPp(k) =           dp(k  )/(dp(k  )+dp(k+1))
        D_C_MCp(k) =           dp(k  )/(dp(k-1)+dp(k  ))
        D_M_Tp (k) =           dp(k-1)/dp_total
        D_CP_Tp(k) = (dp(k  )+dp(k+1))/dp_total
        D_P_Tp (k) =           dp(k+1)/dp_total
        D_MC_Tp(k) = (dp(k-1)+dp(k  ))/dp_total
      enddo
!
! ****** Enforce periodicity.
!
      call set_periodic_bc_1d (D_C_CPp,np)
      call set_periodic_bc_1d (D_C_MCp,np)
      call set_periodic_bc_1d (D_M_Tp,np)
      call set_periodic_bc_1d (D_CP_Tp,np)
      call set_periodic_bc_1d (D_P_Tp,np)
      call set_periodic_bc_1d (D_MC_Tp,np)
!
!$acc enter data copyin(D_C_CPt,D_C_MCt, &
!$acc&                  D_M_Tt, D_CP_Tt,D_P_Tt, D_MC_Tt, &
!$acc&                  D_C_CPp,D_C_MCp, &
!$acc&                  D_M_Tp, D_CP_Tp,D_P_Tp, D_MC_Tp, &
!$acc&                  alpha_t,alpha_p)
!
end subroutine
!#######################################################################
subroutine set_lf_alpha
!
!-----------------------------------------------------------------------
!
! ****** Set Local Lax-Friedrichs alpha foruse with WENO3.
! ****** This should be called any time the velocity field changes.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields, ONLY : vt,vp
      use weno
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
! ****** Get alpha for theta:
!
      do concurrent (i=1:nr,k=1:npm,j=3:ntm-1)
        alpha_t(j,k,i) = MAX(ABS(vt(j-2,k,i)), &
                             ABS(vt(j-1,k,i)), &
                             ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)), &
                             ABS(vt(j+2,k,i)))
      enddo
      j = 2
      do concurrent (i=1:nr,k=1:npm)
        alpha_t(j,k,i) = MAX(ABS(vt(j-1,k,i)), &
                             ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)), &
                             ABS(vt(j+2,k,i)))
      enddo
      j = 1
      do concurrent (i=1:nr,k=1:npm)
        alpha_t(j,k,i) = MAX(ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)), &
                             ABS(vt(j+2,k,i)))
      enddo
      j = ntm
      do concurrent (i=1:nr,k=1:npm)
        alpha_t(j,k,i) = MAX(ABS(vt(j  ,k,i)), &
                             ABS(vt(j-1,k,i)), &
                             ABS(vt(j-2,k,i)))
      enddo
!
! ****** Get alpha for phi:
!
      do concurrent (i=1:nr,k=3:np-2,j=2:ntm-1)
        alpha_p(j,k,i) = MAX(ABS(vp(j,k-2,i)), &
                             ABS(vp(j,k-1,i)), &
                             ABS(vp(j,k  ,i)), &
                             ABS(vp(j,k+1,i)), &
                             ABS(vp(j,k+2,i)))
      enddo
!
      k = 2
      do concurrent (i=1:nr,j=2:ntm-1)
        alpha_p(j,k,i) = MAX(ABS(vp(j,np-2,i)), &
                             ABS(vp(j,k-1 ,i)), &
                             ABS(vp(j,k   ,i)), &
                             ABS(vp(j,k+1 ,i)), &
                             ABS(vp(j,k+2 ,i)))
      enddo
!
      k = npm-1
      do concurrent (i=1:nr,j=2:ntm-1)
        alpha_p(j,k,i) = MAX(ABS(vp(j,k-2,i)), &
                             ABS(vp(j,k-1,i)), &
                             ABS(vp(j,k  ,i)), &
                             ABS(vp(j,k+1,i)), &
                             ABS(vp(j,3  ,i)))
      enddo
!
      call set_periodic_bc_3d (alpha_p,ntm,npm,nr)
!
end subroutine
!#######################################################################
subroutine load_source
!
!-----------------------------------------------------------------------
!
! ****** Define the source term for the diffusion equation.
!
! ****** This is read in from the file SOURCEFILE if it is not blank.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use mpidefs
      use fields
      use input_parameters
      use constants
      use read_2d_file_interface
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: k,ierr,i
      real(r_typ), dimension(:), allocatable :: fn1,fs1
!
      integer :: nft,nfp
      real(r_typ), dimension(:), allocatable :: tf,pf
      real(r_typ), dimension(:,:), allocatable :: sf
      real(r_typ), dimension(:,:), allocatable :: f_tmp2d
!
!-----------------------------------------------------------------------
!
! ****** Read the source file if it was specified.
!
      if (source_filename.ne.' ') then
!
! ****** Allocate memory for the source term.
!
        allocate (source(ntm,npm,nr))
        allocate (f_tmp2d(ntm,npm))
!
        source(:,:,:) = 0.
!
        if (iamp0) then
          call read_2d_file (source_filename,nfp,nft,sf,pf,tf,ierr)
        endif
        wtime_tmp_mpi = MPI_Wtime()
        call MPI_Bcast (nfp,1,ntype_real,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (nft,1,ntype_real,0,MPI_COMM_WORLD,ierr)
        if (.NOT. iamp0) then
          allocate (pf(nfp))
          allocate (tf(nft))
          allocate (sf(nfp,nft))
        endif
        call MPI_Bcast (sf,nfp*nft,ntype_real,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (pf,nfp,ntype_real,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (tf,nft,ntype_real,0,MPI_COMM_WORLD,ierr)
        wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
        sf(:,:) = TRANSPOSE(sf(:,:))
!
! ****** Interpolate the source term onto the main mesh (t,p).
!
        call interp2d (nft,nfp,tf,pf,sf,ntm,npm,t,p,f_tmp2d,ierr)
        do i=1,nr
          source(:,:,i)=f_tmp2d(:,:)
        enddo
!
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in LOAD_SOURCE:'
          write (*,*) '### The scales in the source term file are'// &
                        ' not monotonically increasing.'
          write (*,*) 'File name: ',trim(source_filename)
          call endrun(.true.)
        end if
!
! ****** Enforce periodicity.
!
        allocate (fn1(nr))
        allocate (fs1(nr))
        !
        do i=1,nr
          source(:,   1,i)=half*(source(:,1,i)+source(:,npm,i))
          source(:,npm-1,i)=source(:,1,i)
!
! ****** Set the pole value to only have an m=0 component.
!
          fn1(:)=0.
          fs1(:)=0.
          do k=1,npm-2
            fn1(i)=fn1(i)+source(  1,k,i)*dp(k)
            fs1(i)=fs1(i)+source(ntm,k,i)*dp(k)
          enddo
          fn1(i)=fn1(i)*twopi_i
          fs1(i)=fs1(i)*twopi_i
!
          source(  1,:,i)=fn1(i)
          source(ntm,:,i)=fs1(i)
        enddo
!
        if (verbose) then
          write (*,*)
          write (*,*) 'A source term was read in from file: ', &
                     trim(source_filename)
          write (*,*) 'Minimum value = ',minval(source(:,:,1))
          write (*,*) 'Maximum value = ',maxval(source(:,:,1))
        end if
!
        deallocate (f_tmp2d)
        deallocate (fn1)
        deallocate (fs1)
        deallocate (sf)
        deallocate (tf)
        deallocate (pf)
!
      end if
!
end subroutine
!#######################################################################
subroutine add_flow_differential_rotation_aft
!
!-----------------------------------------------------------------------
!
! ****** Add the AFT differential rotation flow to vp.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Parameters in m/s
!
      real(r_typ), dimension(:), allocatable :: t0,t2,t4
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
      allocate (t0(nr))
      allocate (t2(nr))
      allocate (t4(nr))
!$acc enter data create(t0,t2,t4)
      do concurrent (i=1:nr)
        t0(i) = 46.0_r_typ
        t2(i) = -262.0_r_typ
        t4(i) = -379.0_r_typ
      enddo
!
      do concurrent (i=1:nr,k=1:np,j=1:ntm)
        vp(j,k,i) = vp(j,k,i) + &
                  m_s_to_rs_hr*st(j)*(t0(i) + t2(i)*ct(j)**2 + t4(i)*ct(j)**4)
      enddo
!
!$acc exit data delete(t0,t2,t4)
      deallocate (t0)
      deallocate (t2)
      deallocate (t4)
!
end subroutine
!#######################################################################
subroutine add_flow_meridianal_aft
!
!-----------------------------------------------------------------------
!
! ****** Add the AFT meridianal flow to vt.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Paramters in m/s
!
      real(r_typ), dimension(:), allocatable :: s1,s3,s5
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
      allocate (s1(nr))
      allocate (s3(nr))
      allocate (s5(nr))
!$acc enter data create(s1,s3,s5)
      do concurrent (i=1:nr)
        s1(i) = 22.0_r_typ
        s3(i) = 11.0_r_typ
        s5(i) = -28.0_r_typ
      enddo
!
      do concurrent (i=1:nr,k=1:npm,j=1:nt)
        vt(j,k,i) = vt(j,k,i) -            &
                  m_s_to_rs_hr*sth(j)* &
                  (s1(i)*cth(j) + s3(i)*cth(j)**3 + s5(i)*cth(j)**5)
      enddo
!
!$acc exit data delete(s1,s3,s5)
      deallocate (s1)
      deallocate (s3)
      deallocate (s5)
!
end subroutine
!#######################################################################
subroutine get_flow_dtmax (dtmaxflow)
!
!-----------------------------------------------------------------------
!
! ****** Get the maximum time step that can be used for stable
! ****** (explicit) advection for all realizations.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use mesh
      use fields, ONLY : vt,vp
      use constants
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(OUT) :: dtmaxflow
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ierr
      real(r_typ) :: kdotv,deltat,dtmax
!
!-----------------------------------------------------------------------
!
      dtmax = huge(one)
!
!$omp parallel do   collapse(3) default(shared)  reduction (min:dtmax)
!$acc parallel loop collapse(3) default(present) reduction (min:dtmax)
      do i=1,nr
        do k=2,npm-1
          do j=2,ntm-1
            kdotv = MAX(ABS(vt(j,k,i)),ABS(vt(j+1,k,i)))*dt_i(j) &
                  + MAX(ABS(vp(j,k,i)),ABS(vp(j,k+1,i)))*st_i(j)*dp_i(k)
            deltat = one/MAX(kdotv,small_value)
            dtmax = MIN(dtmax,deltat)
          enddo
        enddo
      enddo
!$omp end parallel do
!
      dtmaxflow = half*safety*dtmax
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Allreduce (MPI_IN_PLACE,dtmaxflow,1,ntype_real, &
                          MPI_MIN,MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
end subroutine
!#######################################################################
subroutine diffusion_step_euler_cd (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Diffuse the field by a time step using 1st order Euler in time
! ****** and 2nd-order central differencing in space.
! ****** This routine does not use the matrix coefs which it could.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields
      use input_parameters
      use globals
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dtime_local
!
!-----------------------------------------------------------------------
!
      integer :: j,k,el
      integer*8 :: i
      real(r_typ) :: dtime_euler_used
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary field.
!
      allocate (fold(ntm,npm,nr))
!$acc enter data create(fold)
!
! ****** Subcycle at a stable time-step.
!
      n_stable_diffusion_cycles = CEILING(dtime_local/dtime_diffusion_euler,8)
      dtime_euler_used = dtime_local/n_stable_diffusion_cycles
!
      do i=1,n_stable_diffusion_cycles
!
        call ax (f,fold)
!
        do concurrent (el=1:nr,k=1:npm,j=1:ntm)
          f(j,k,el) = f(j,k,el) + dtime_euler_used*fold(j,k,el)
        enddo
!
      enddo
!
!$acc exit data delete(fold)
      deallocate (fold)
!
end subroutine
!#######################################################################
subroutine advection_step_fe_upwind (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with advection by the time-step dtime_local
! ****** using the Foward-Euler+Upwind method.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use constants
      use mesh
      use fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(IN) :: dtime_local
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ) :: cct,ccp
      real(r_typ), dimension(:), allocatable :: fn,fs
      real(r_typ), dimension(:,:,:), allocatable :: flux_t,flux_p
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary flux arrays on the half-mesh.
!
      allocate (flux_t(nt,npm,nr))
      allocate (flux_p(ntm,np,nr))
!$acc enter data create(flux_t,flux_p)
!
      do concurrent (i=1:nr,k=1:npm,j=1:nt)
        flux_t(j,k,i) = zero
      enddo
!
      do concurrent (i=1:nr,k=1:np,j=1:ntm)
        flux_p(j,k,i) = zero
      enddo
!
! ****** Compute the fluxes at the cell faces.
!
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm1)
        cct = sign(upwind,vt(j,k,i))
        flux_t(j,k,i) = vt(j,k,i)*half*((one-cct)*f(j,k,i)+(one+cct)*f(j-1,k,i))
      enddo
!
      do concurrent (i=1:nr,k=2:npm1,j=2:ntm-1)
        ccp = sign(upwind,vp(j,k,i))
        flux_p(j,k,i) = vp(j,k,i)*half*((one-ccp)*f(j,k,i)+(one+ccp)*f(j,k-1,i))
      enddo
!
! ****** Set periodicity of the flux (seam).
!
      call set_periodic_bc_3d (flux_t,nt,npm,nr)
      call set_periodic_bc_3d (flux_p,ntm,np,nr)
!
! ****** Advect F by one time step.
!
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1)
        f(j,k,i) = f(j,k,i) - dtime_local*( (  sth(j+1)*flux_t(j+1,k,i) &
                                             - sth(j  )*flux_t(j  ,k,i) &
                                            )*st_i(j)*dt_i(j)           &
                                           +(  flux_p(j,k+1,i)          &
                                             - flux_p(j,k  ,i)          &
                                            )*st_i(j)*dp_i(k)           &
                                          )
      enddo
!
! ****** Advect the values at the poles.
!
      allocate (fn(nr))
      allocate (fs(nr))
!$acc enter data create(fn,fs)
!
      do concurrent (i=1:nr)
        fn(i) = zero
        fs(i) = zero
      enddo
!
!$omp parallel do   collapse(2) default(shared)
!$acc parallel loop collapse(2) default(present)
        do i=1,nr
          do k=2,npm-1
!$omp atomic update
!$acc atomic update
            fn(i) = fn(i) + flux_t(   2,k,i)*dp(k)
!$omp atomic update
!$acc atomic update
            fs(i) = fs(i) + flux_t(ntm1,k,i)*dp(k)
          enddo
        enddo
!$omp end parallel do
! ****** Note that the south pole needs a sign change since the
! ****** theta flux direction is reversed.
      do concurrent (i=1:nr,k=2:npm-1)
        f(  1,k,i) = f(  1,k,i) - fn(i)*dtime_local*two*pi_i*dt_i(  1)
        f(ntm,k,i) = f(ntm,k,i) + fs(i)*dtime_local*two*pi_i*dt_i(ntm)
      enddo
!
! ****** Set periodic boundary condition.
!
      call set_periodic_bc_3d (f,ntm,npm,nr)
!
!$acc exit data delete(flux_t,flux_p,fn,fs)
      deallocate (flux_t)
      deallocate (flux_p)
      deallocate (fn)
      deallocate (fs)
!
end subroutine
!#######################################################################
subroutine advection_step_rk3tvd_upwind (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with advection by the time-step dtime_local
! ****** using the RK3TVD+Upwind method.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use constants
      use mesh
      use fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(IN) :: dtime_local
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ), dimension(:,:,:), allocatable :: f1,f2,aop
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays for RK.
!
      allocate (f1(ntm,npm,nr))
      allocate (f2(ntm,npm,nr))
      allocate (aop(ntm,npm,nr))
!$acc enter data create(f1,f2,aop)
!
      call advection_operator_upwind (f,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f1(j,k,i) = f(j,k,i) - dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator_upwind (f1,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f2(j,k,i) = three_quarter*f(j,k,i) + quarter*f1(j,k,i) &
                                       - quarter*dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator_upwind (f2,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f(j,k,i) = third*f(j,k,i) + two_third*f2(j,k,i) &
                              - two_third*dtime_local*aop(j,k,i)
      enddo
!
! ****** Deallocate temporary arrays.
!
!$acc exit data delete(f1,f2,aop)
      deallocate (f1)
      deallocate (f2)
      deallocate (aop)
!
end subroutine
!#######################################################################
subroutine advection_operator_upwind (ftemp,aop)
!
!-----------------------------------------------------------------------
!
! ****** Compute DIV(v*f) using Upwinding
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use constants
      use mesh
      use fields, ONLY: vt,vp
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntm,npm,nr), INTENT(IN) :: ftemp
      real(r_typ), dimension(ntm,npm,nr), INTENT(OUT) :: aop
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ) :: cct,ccp
      real(r_typ), dimension(:), allocatable :: fn,fs
      real(r_typ), dimension(:,:,:), allocatable :: flux_t,flux_p
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays
! ****** (we set all points here so no initialization needed)
!
      allocate (flux_t(nt,npm,nr))
      allocate (flux_p(ntm,np,nr))
!$acc enter data create(flux_t,flux_p)
!
! ****** Compute the fluxes at the cell faces.
!
      do concurrent (i=1:nr,k=1:npm,j=2:ntm1)
        cct = sign(upwind,vt(j,k,i))
        flux_t(j,k,i) = vt(j,k,i)*half*((one-cct)*ftemp(j  ,k,i) &
                                      + (one+cct)*ftemp(j-1,k,i))
      enddo
!
      do concurrent (i=1:nr,k=2:npm1,j=2:ntm-1)
        ccp = sign(upwind,vp(j,k,i))
        flux_p(j,k,i) = vp(j,k,i)*half*((one-ccp)*ftemp(j,k  ,i) &
                                      + (one+ccp)*ftemp(j,k-1,i))
      enddo
!
! ****** Set periodicity of the fluxp (seam).
!
      call set_periodic_bc_3d (flux_p,ntm,np,nr)
!
! ****** Compute advection operator F.
!
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1)
        aop(j,k,i) = (  (  sth(j+1)*flux_t(j+1,k,i)  &
                         - sth(j  )*flux_t(j  ,k,i)  &
                        )*st_i(j)*dt_i(j)            &
                      + (           flux_p(j,k+1,i)  &
                         -          flux_p(j,k  ,i)  &
                        )*st_i(j)*dp_i(k)            &
                     )
      enddo
!
! ****** Get the advection operator at the poles.
!
      allocate (fn(nr))
      allocate (fs(nr))
!$acc enter data create(fn,fs)
!
      do concurrent (i=1:nr)
        fn(i) = zero
        fs(i) = zero
      enddo
!
!$omp parallel do   collapse(2) default(shared)
!$acc parallel loop collapse(2) default(present)
        do i=1,nr
          do k=2,npm-1
!$omp atomic update
!$acc atomic update
            fn(i) = fn(i) + flux_t(   2,k,i)*dp(k)
!$omp atomic update
!$acc atomic update
            fs(i) = fs(i) + flux_t(ntm1,k,i)*dp(k)
          enddo
        enddo
!$omp end parallel do
! ****** Note that the south pole needs a sign change since the
! ****** theta flux direction is reversed.
        do concurrent (i=1:nr,k=2:npm-1)
          aop(  1,k,i) =  fn(i)*two*pi_i*dt_i(  1)
          aop(ntm,k,i) = -fs(i)*two*pi_i*dt_i(ntm)
        enddo
!
! ****** Set periodic boundary condition.
!
      call set_periodic_bc_3d (aop,ntm,npm,nr)
!
!$acc exit data delete(flux_t,flux_p,fn,fs)
      deallocate (flux_t)
      deallocate (flux_p)
      deallocate (fn)
      deallocate (fs)
!
end subroutine
!#######################################################################
subroutine advection_operator_weno3 (ftemp,aop)
!
!-----------------------------------------------------------------------
!
! ****** Compute DIV(v*f) with WENO3 using
! ****** Local Lax-Friedrichs flux splitting and 2nd order
! ****** theta-direction approximation reconstruction.
! ****** The scheme should be 4th-order in phi and 2nd-order in theta
! ****** in smooth regions away from the pole.
! ****** Next to the pole, theta becomes 1st order (upwinding).
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use constants
      use mesh
      use weno
      use fields, ONLY: vt,vp
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntm,npm,nr), INTENT(IN) :: ftemp
      real(r_typ), dimension(ntm,npm,nr), INTENT(OUT) :: aop
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ) :: cct
      real(r_typ), dimension(:), allocatable :: fn,fs
      real(r_typ), dimension(:,:,:), allocatable :: flux_t,flux_p
      real(r_typ), dimension(:,:,:), allocatable :: LP,LN
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays
! ****** (we set all points here so no initialization needed)
!
      allocate (flux_t(nt,npm,nr))
      allocate (LP(ntm,npm,nr))
      allocate (LN(ntm,npm,nr))
!$acc enter data create(flux_t,LP,LN)
!
! ****** Compute Lax-Friedrichs fluxes.
!
      do concurrent (i=1:nr,k=1:npm,j=1:ntm1)
        LP(j,k,i) = half*ftemp(j,k,i)*(vt(j+1,k,i) - alpha_t(j,k,i))
        LN(j,k,i) = half*ftemp(j,k,i)*(vt(j  ,k,i) + alpha_t(j,k,i))
      enddo
!
! ***** No need to seam since k loop covers all points.
!
      do concurrent (i=1:nr,k=1:npm,j=3:nt-2)
!
        p0m = (one + D_C_MCt(j-1))*LN(j-1,k,i) - D_C_MCt(j-1)*LN(j-2,k,i)
        p1m =        D_C_MCt(j  ) *LN(j-1,k,i) + D_C_CPt(j-1)*LN(j  ,k,i)
        p0p = (one + D_C_CPt(j  ))*LP(j  ,k,i) - D_C_CPt(j  )*LP(j+1,k,i)
        p1p =        D_C_CPt(j-1) *LP(j  ,k,i) + D_C_MCt(j  )*LP(j-1,k,i)
!
        B0m = four*(D_C_MCt(j-1)*(LN(j-1,k,i) - LN(j-2,k,i)))**2
        B1m = four*(D_C_CPt(j-1)*(LN(j  ,k,i) - LN(j-1,k,i)))**2
        B0p = four*(D_C_CPt(j  )*(LP(j+1,k,i) - LP(j  ,k,i)))**2
        B1p = four*(D_C_MCt(j  )*(LP(j  ,k,i) - LP(j-1,k,i)))**2
!
        w0m = D_P_Tt(j-1) *(one/(weno_eps + B0m)**2)
        w1m = D_MC_Tt(j-1)*(one/(weno_eps + B1m)**2)
        w0p = D_M_Tt(j) *(one/(weno_eps + B0p)**2)
        w1p = D_CP_Tt(j)*(one/(weno_eps + B1p)**2)
!
        wm_sum = w0m + w1m
        wp_sum = w0p + w1p
!
        OM0m = w0m/wm_sum
        OM1m = w1m/wm_sum
        OM0p = w0p/wp_sum
        OM1p = w1p/wp_sum
!
        um = OM0m*p0m + OM1m*p1m
        up = OM0p*p0p + OM1p*p1p
!
        flux_t(j,k,i) = up + um
!
      enddo
!
! ****** Compute upwind theta flux for points near the pole.
!
      do concurrent (i=1:nr,k=1:npm)
!
        cct=sign(upwind,vt(2,k,i))
        flux_t(2,k,i) = vt(2,k,i)*half*((one-cct)*ftemp(2,k,i) &
                                  + (one+cct)*ftemp(1,k,i))
!
        cct=sign(upwind,vt(ntm1,k,i))
        flux_t(ntm1,k,i) = vt(ntm1,k,i)*half*((one-cct)*ftemp(ntm1,k,i) &
                                        + (one+cct)*ftemp(ntm2,k,i))
!
      enddo
!
!$acc exit data delete (LP,LN)
      deallocate (LP)
      deallocate (LN)
!
! ### PHI ######################################
!
      allocate (flux_p(ntm,np,nr))
      allocate (LP(ntm,0:npm,nr))
      allocate (LN(ntm,0:npm,nr))
!$acc enter data create(flux_p,LP,LN)
!
! ****** Compute flux splitting:
!
      ! Compute N+(i) and P-(i):
      do concurrent (i=1:nr,k=1:npm-1,j=2:ntm-1)
        LP(j,k,i) = half*ftemp(j,k,i)*(vp(j,k+1,i) - alpha_p(j,k,i))
        LN(j,k,i) = half*ftemp(j,k,i)*(vp(j,k  ,i) + alpha_p(j,k,i))
      enddo
!
      do concurrent (i=1:nr,j=2:ntm-1)
        LP(j,0,i) = half*ftemp(j,npm-2,i)*(vp(j,npm-1,i) - alpha_p(j,npm-2,i))
        LN(j,0,i) = half*ftemp(j,npm-2,i)*(vp(j,npm-2,i) + alpha_p(j,npm-2,i))
      enddo
!
      do concurrent (i=1:nr,j=2:ntm-1)
        LP(j,npm,i) = half*ftemp(j,npm,i)*(vp(j,3,i) - alpha_p(j,npm,i))
        LN(j,npm,i) = half*ftemp(j,npm,i)*(vp(j,2,i) + alpha_p(j,npm,i))
      enddo
!
! ***** No need to seam since k loop covers all.
!
! ****** Now compute f+(i-1/2) and f-(i-1/2)
! ****** Note that N-(i) = N+(i-1)
!
      do concurrent (i=1:nr,k=2:npm1,j=2:ntm-1)
!
        p0m = (one + D_C_MCp(k-1))*LN(j,k-1,i) - D_C_MCp(k-1)*LN(j,k-2,i)
        p1m =        D_C_MCp(k  ) *LN(j,k-1,i) + D_C_CPp(k-1)*LN(j,k  ,i)
        p0p = (one + D_C_CPp(k  ))*LP(j,k  ,i) - D_C_CPp(k  )*LP(j,k+1,i)
        p1p =        D_C_CPp(k-1) *LP(j,k  ,i) + D_C_MCp(k  )*LP(j,k-1,i)
!
        B0m = four*(D_C_MCp(k-1)*(LN(j,k-1,i) - LN(j,k-2,i)))**2
        B1m = four*(D_C_CPp(k-1)*(LN(j,k  ,i) - LN(j,k-1,i)))**2
        B0p = four*(D_C_CPp(k  )*(LP(j,k+1,i) - LP(j,k  ,i)))**2
        B1p = four*(D_C_MCp(k  )*(LP(j,k  ,i) - LP(j,k-1,i)))**2
!
        w0m = D_P_Tp (k-1)*(one/(weno_eps + B0m)**2)
        w1m = D_MC_Tp(k-1)*(one/(weno_eps + B1m)**2)
        w0p = D_M_Tp (k)*(one/(weno_eps + B0p)**2)
        w1p = D_CP_Tp(k)*(one/(weno_eps + B1p)**2)
!
        wm_sum = w0m + w1m
        wp_sum = w0p + w1p
!
        OM0m = w0m/wm_sum
        OM1m = w1m/wm_sum
        OM0p = w0p/wp_sum
        OM1p = w1p/wp_sum
!
        um = OM0m*p0m + OM1m*p1m
        up = OM0p*p0p + OM1p*p1p
!
        flux_p(j,k,i) = up + um
!
      enddo
!
      call set_periodic_bc_3d (flux_p,ntm,np,nr)
!
!$acc exit data delete (LP,LN)
      deallocate (LP)
      deallocate (LN)
!
! ****** Now evaluation the advection operator for internal points.
!
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1)
        aop(j,k,i) = (  (  sth(j+1)*flux_t(j+1,k,i)  &
                         - sth(j  )*flux_t(j  ,k,i)  &
                          )*st_i(j)*dt_i(j)          &
                      + (           flux_p(j,k+1,i)  &
                         -          flux_p(j,k  ,i)  &
                          )*st_i(j)*dp_i(k)          &
                     )
      enddo
!
! ****** Get the advection operator at the poles.
!
      allocate (fn(nr))
      allocate (fs(nr))
!$acc enter data create(fn,fs)
!
      do concurrent (i=1:nr)
        fn(i) = zero
        fs(i) = zero
      enddo
!
!$omp parallel do   collapse(2) default(shared)
!$acc parallel loop collapse(2) default(present)
      do i=1,nr
        do k=2,npm-1
!$omp atomic update
!$acc atomic update
          fn(i) = fn(i) + flux_t(   2,k,i)*dp(k)
!$omp atomic update
!$acc atomic update
          fs(i) = fs(i) + flux_t(ntm1,k,i)*dp(k)
        enddo
      enddo
!$omp end parallel do
! ****** Note that the south pole needs a sign change since the
! ****** theta flux direction is reversed.
      do concurrent (i=1:nr,k=2:npm-1)
        aop(  1,k,i) =  fn(i)*two*pi_i*dt_i(  1)
        aop(ntm,k,i) = -fs(i)*two*pi_i*dt_i(ntm)
      enddo
!
! ****** Set periodic phi boundary condition.
!
      call set_periodic_bc_3d (aop,ntm,npm,nr)
!
!$acc exit data delete(flux_t,flux_p,fn,fs)
      deallocate (flux_t)
      deallocate (flux_p)
      deallocate (fn)
      deallocate (fs)
!
end subroutine
!#######################################################################
subroutine advection_step_rk3tvd_weno3 (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with advection by the time-step dtime_local
! ****** using the RK3TVD+WENO3 method.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use constants
      use mesh
      use fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(IN) :: dtime_local
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ), dimension(:,:,:), allocatable :: f1,f2,aop
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays for RK.
!
      allocate (f1(ntm,npm,nr))
      allocate (f2(ntm,npm,nr))
      allocate (aop(ntm,npm,nr))
!$acc enter data create(f1,f2,aop)
!
      call advection_operator_weno3 (f,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f1(j,k,i) = f(j,k,i) - dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator_weno3 (f1,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f2(j,k,i) = three_quarter*f(j,k,i) + quarter*f1(j,k,i) &
                                           - quarter*dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator_weno3 (f2,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f(j,k,i) = third*f(j,k,i) + two_third*f2(j,k,i) &
                                  - two_third*dtime_local*aop(j,k,i)
      enddo
!
! ****** Deallocate temporary arrays.
!
!$acc exit data delete(f1,f2,aop)
      deallocate (f1)
      deallocate (f2)
      deallocate (aop)
!
end subroutine
!#######################################################################
subroutine ax (x,y)
!
!-----------------------------------------------------------------------
!
! ****** Apply the diffuse operator A as y=Ax.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
      use mesh
      use fields
      use input_parameters
      use matrix_storage
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ), dimension(:), allocatable :: fn2_fn1,fs2_fs1
      real(r_typ), dimension(ntm,npm,nr), INTENT(IN) :: x
      real(r_typ), dimension(ntm,npm,nr), INTENT(OUT) :: y
!
!-----------------------------------------------------------------------
!
! ****** Compute y=Ax.
!
! ****** Compute inner points.
!
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1)
        y(j,k,i) =  coef(j,k,1,i)*x(j,  k-1,i) &
                  + coef(j,k,2,i)*x(j-1,k  ,i) &
                  + coef(j,k,3,i)*x(j  ,k  ,i) &
                  + coef(j,k,4,i)*x(j+1,k  ,i) &
                  + coef(j,k,5,i)*x(j,  k+1,i)
      enddo
!
! ****** Compute boundary points.
!
! ****** Get the m=0 components near the poles.
!
      allocate (fn2_fn1(nr))
      allocate (fs2_fs1(nr))
!$acc enter data create(fn2_fn1,fs2_fs1)
!
      do concurrent (i=1:nr)
        fn2_fn1(i) = zero
        fs2_fs1(i) = zero
      enddo
!
!$omp parallel do   collapse(2) default(shared)
!$acc parallel loop collapse(2) default(present)
      do i=1,nr
        do k=2,npm-1
!$omp atomic update
!$acc atomic update
          fn2_fn1(i) = fn2_fn1(i) + (diffusion_coef(1    ,k,i)     &
                                  +  diffusion_coef(2    ,k,i))    &
                                  * (x(2    ,k,i) - x(1  ,k,i))*dp(k)
!$omp atomic update
!$acc atomic update
          fs2_fs1(i) = fs2_fs1(i) + (diffusion_coef(nt-1 ,k,i)     &
                                  +  diffusion_coef(nt   ,k,i))    &
                                  * (x(ntm-1,k,i) - x(ntm,k,i))*dp(k)
        enddo
      enddo
!$omp end parallel do
!
      do concurrent (i=1:nr,k=1:npm)
        y(  1,k,i) = fn2_fn1(i)*dt_i(  1)*dt_i(  1)*pi_i
        y(ntm,k,i) = fs2_fs1(i)*dt_i(ntm)*dt_i(ntm)*pi_i
      enddo
!
! ****** Set the periodic boundary conditions.
!
      call set_periodic_bc_3d (y,ntm,npm,nr)
!$acc exit data delete(fn2_fn1,fs2_fs1)
      deallocate (fn2_fn1)
      deallocate (fs2_fs1)
!
end subroutine
!#######################################################################
subroutine get_m0 (f,fn1,fn2,fs1,fs2)
!
!-----------------------------------------------------------------------
!
! ****** Get the m=0 component of the field F near the North and
! ****** South poles.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
      use mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntm,npm,nr) :: f
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nr) :: fn1,fn2,fs1,fs2
!
!-----------------------------------------------------------------------
!
      integer :: k,i
!
!-----------------------------------------------------------------------
!
      fn1(:) = 0.
      fn2(:) = 0.
      fs1(:) = 0.
      fs2(:) = 0.
!
      do i=1,nr
        do k=2,npm-1
          fn1(i) = fn1(i) + f(1,k,i)*dp(k)
          fn2(i) = fn2(i) + f(2,k,i)*dp(k)
          fs1(i) = fs1(i) + f(ntm,k,i)*dp(k)
          fs2(i) = fs2(i) + f(ntm-1,k,i)*dp(k)
        enddo
      enddo
!
      do i=1,nr
        fn1(i) = fn1(i)*twopi_i
        fn2(i) = fn2(i)*twopi_i
        fs1(i) = fs1(i)*twopi_i
        fs2(i) = fs2(i)*twopi_i
      enddo
!
end subroutine
!#######################################################################
subroutine interp2d (nxi,nyi,xi,yi,fi,nx,ny,x,y,f,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate a 2D field from array FI(NXI,NYI), defined
! ****** on the mesh XI(NXI) x YI(NYI), into the array F(NX,NY),
! ****** defined on the mesh X(NX) x Y(NY).
!
! ****** Zero values are returned at data points outside the
! ****** bounds of the XI x YI mesh.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nxi,nyi
      real(r_typ), dimension(nxi) :: xi
      real(r_typ), dimension(nyi) :: yi
      real(r_typ), dimension(nxi,nyi) :: fi
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: i,j,iip1,jjp1
      integer :: ii=0,jj=0
      real(r_typ) :: dummy,ax,ay,xv,yv
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: flint
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Check that the scales XI and YI are monotonic.
!
      dummy=flint(.true.,zero,nxi,xi,xi,ierr)
      dummy=flint(.true.,zero,nyi,yi,yi,ierr)
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in INTRP2D:'
        write (*,*) '### Scales are not monotonically increasing.'
        ierr=1
        return
      end if
!
! ****** Interpolate the data.
!
      do j=1,ny
        yv=y(j)
        if (yv.lt.yi(1).or.yv.gt.yi(nyi)) then
          f(:,j)=0.
          cycle
        else
          call interp (nyi,yi,yv,jj,jjp1,ay,ierr)
          if (ierr.ne.0) then
            f(:,j)=0.
            cycle
          end if
        end if
        do i=1,nx
          xv=x(i)
          if (xv.lt.xi(1).or.xv.gt.xi(nxi)) then
            f(i,j)=0.
            cycle
          else
            call interp (nxi,xi,xv,ii,iip1,ax,ierr)
            if (ierr.ne.0) then
              f(i,j)=0.
              cycle
            end if
          end if
          f(i,j)=(one-ax)*((one-ay)*fi(ii  ,jj  )+ay*fi(ii  ,jjp1)) &
                 +ax *((one-ay)*fi(iip1,jj  )+ay*fi(iip1,jjp1))
        enddo
      enddo
!
end subroutine
!#######################################################################
function flint (check,x,n,xn,fn,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate a function linearly.
!
!-----------------------------------------------------------------------
!
! ****** The funcion is defined at N nodes, with values given by
! ****** FN(N) at positions XN(N).  The function value returned is
! ****** the linear interpolant at X.
!
! ****** Note that if X.lt.XN(1), the function value returned
! ****** is FN(1), and if X.gt.XN(N), the function value returned
! ****** is FN(N).
!
! ****** Call once with CHECK=.true. to check that the values
! ****** in XN(N) are monotonically increasing.  In this mode
! ****** the array XN(N) is checked, and X and FN(N) are not
! ****** accessed.  If the check is passed, IERR=0 is returned.
! ****** Otherwise, IERR=1 is returned.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical :: check
      real(r_typ) :: x
      integer :: n
      real(r_typ), dimension(n) :: xn,fn
      integer :: ierr
      real(r_typ) :: flint
!
!-----------------------------------------------------------------------
!
      integer :: i
      real(r_typ) :: x1,x2,alpha
!
!-----------------------------------------------------------------------
!
      ierr=0
      flint=0.
!
! ****** If CHECK=.true., check the abscissa table.
!
      if (check) then
        if (n.le.0) then
          write (*,*)
          write (*,*) '### ERROR in FLINT:'
          write (*,*) '### Invalid abscissa table dimension.'
          write (*,*) 'N = ',n
          ierr=1
          return
        end if
        do i=1,n-1
          if (xn(i+1).le.xn(i)) then
            write (*,*)
            write (*,*) '### ERROR in FLINT:'
            write (*,*) '### Abscissa table values are not'// &
                       ' monotonically increasing.'
            write (*,*) 'N = ',n
            write (*,*) 'XN = ',xn
            ierr=1
            return
          end if
        enddo
        return
      end if
!
! ****** Get the interpolated value.
!
      if (x.le.xn(1)) then
        flint=fn(1)
      else if (x.gt.xn(n)) then
        flint=fn(n)
      else
        do i=1,n-1
          if (x.ge.xn(i).and.x.lt.xn(i+1)) exit
        enddo
        x1=xn(i)
        x2=xn(i+1)
        alpha=(x-x1)/(x2-x1)
        flint=fn(i)*(one-alpha)+fn(i+1)*alpha
      end if
!
      return
end function
!#######################################################################
subroutine interp (n,x,xv,i,ip1,a,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Get the interpolation factor at XV from the table X(N).
!
!-----------------------------------------------------------------------
!
! ****** This routine does not do the actual interpolation.  Use the
! ****** returned values of I, IP1, and A to get the interpolated
! ****** value.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i,ip1
      real(r_typ) :: a
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Check if the x-scale has only one point.
!
      if (n.eq.1.and.xv.eq.x(1)) then
        ip1=1
        a=0.
        return
      end if
!
! ****** Find the interval and compute the interpolation factor.
!
      do i=1,n-1
        if (xv.ge.x(i).and.xv.le.x(i+1)) then
          ip1=i+1
          if (x(i).eq.x(i+1)) then
            a=0.
          else
            a=(xv-x(i))/(x(i+1)-x(i))
          end if
          return
        end if
      enddo
!
! ****** ERROR: the value was not found.
!
      ierr=1
!
end subroutine
!#######################################################################
subroutine ffopen (iun,fname,mode,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Open file FNAME and link it to unit IUN.
!
! ****** If there is an error, this routine returns IERR.ne.0.
!
!-----------------------------------------------------------------------
!
! ****** When MODE='r', the file must exist.
! ****** When MODE='w', the file is created.
! ****** When MODE='rw', the file must exist, but can be overwritten.
! ****** When MODE='a', the file is created if it does not exist,
! ******                otherwise, it is appended.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iun
      character(*) :: fname
      character(*) :: mode
      integer :: ierr
      logical :: ex
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      if (mode.eq.'r') then
        open (iun,file=fname,form="FORMATTED",status='old',err=900)
      else if (mode.eq.'rw') then
        open (iun,file=fname,form="FORMATTED",status='replace',err=900)
      else if (mode.eq.'w') then
        open (iun,file=fname,form="FORMATTED",status='new',err=900)
      elseif (mode.eq.'a') then
        inquire(file=fname, exist=ex)
        if (ex) then
          open (iun,file=fname,form="FORMATTED",position='append',err=900)
        else
          open (iun,file=fname,form="FORMATTED",status='new',err=900)
        end if
      else
        write (*,*)
        write (*,*) '### ERROR in FFOPEN:'
        write (*,*) '### Invalid MODE requested.'
        write (*,*) 'MODE = ',mode
        write (*,*) 'File name: ',trim(fname)
        ierr=2
        return
      end if
!
      return
!
  900 continue
!
      write (*,*)
      write (*,*) '### ERROR in FFOPEN:'
      write (*,*) '### Error while opening the requested file.'
      write (*,*) 'File name: ',trim(fname)
      write (*,*) 'MODE = ',mode
      ierr=1
!
end subroutine
!#######################################################################
subroutine read_input_file
!
!-----------------------------------------------------------------------
!
! ****** Read the input file.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      namelist /hipft_input_parameters/                                &
               verbose,                                                &
               initial_map_filename,output_map_root_filename,          &
               output_map_directory,n_realizations,validation_run,     &
               output_map_idx_cadence,output_map_time_cadence,         &
               output_map_2d,restart_run,restart_file,time_end,        &
               dt_max_increase_fac,dt_min,dt_max,strang_splitting,     &
               pole_flux_lat_limit,advance_flow,flow_vp_rigid_omega,   &
               flow_rigid_omega,flow_vt_const,flow_attenuate,          &
               flow_dr_model,flow_mf_model,use_flow_from_files,        &
               flow_list_filename,flow_root_dir,flow_num_method,       &
               upwind,advance_diffusion,diffusion_coef_constant,       &
               diffusion_coef_filename,diffusion_coef_grid,            &
               diffusion_coef_factor,diffusion_num_method,             &
               diffusion_subcycles,advance_source,source_filename,     &
               assimilate_data,assimilate_data_map_list_filename,      &
               assimilate_data_map_root_dir,time_start
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(80) :: infile='hipft.dat'
!
!-----------------------------------------------------------------------
!
! ****** Read the input file.
!
! RMC: Could have just rank 0 read this and then broadcast to all other ranks.
      call ffopen (8,infile,'r',ierr)
!
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in READ_INPUT_FILE:'
          write (*,*) '### Could not open the input file.'
          write (*,*) 'File name: ',trim(infile)
          call endrun(.true.)
        end if
      end if
!
      read (8,hipft_input_parameters)
      close (8)
!
! ****** Check inputs for issues.
!
      if (output_map_time_cadence.gt.0 .and. output_map_idx_cadence.gt.0) then
        if (iamp0) then
          write(*,*) "ERROR!  Cannot use both output_map_time_cadence"
          write(*,*) "        and output_map_idx_cadence at the same time!"
        endif
        call endrun (.true.)
      end if
!
      if (nproc .gt. n_realizations) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in SETUP:'
          write (*,*) '### There are more MPI ranks than realizations.'
          write (*,*) 'Number of MPI ranks    = ',nproc
          write (*,*) 'Number of realizations = ',n_realizations
        end if
        call endrun (.true.)
      end if
!
end subroutine
!#######################################################################
!
!-----------------------------------------------------------------------
!
! HIPFT Update log:
!
! 05/24/2021, RC, Version 0.1.0:
!   - Original version of the program.
!     Derived from DIFFUSE_ADVECT version 2.0.0 of 05/24/2021.
!
! 11/01/2021, RC, Version 0.1.1:
!   - Replaced all non-reduction OpenACC/MP loops with
!     `do concurrent`.  This requires alternative compiler
!     flags to properly parallelize.  See build.sh for
!     examples.
!   - Merged sds and number_types modules/routines into main code.
!
! 11/05/2021, RC, Version 0.2.0:
!   - Added dtmax input parameter to specify maximum allowed dt.
!   - Added history files:
!     history_num.dat contains numerical diagnostics while
!     history_sol.dat contains diagnostics of the solution
!     including flux, and validation comparison.
!   - Updated explicit Euler diffusion to use matrix coefs.
!
! 11/24/2021, RC, Version 0.2.1:
!   - Some simplifications to STS method.
!   - Updated stable flow time-step to be more robust.
!
! 02/17/2022, RC, Version 0.3.0:
!   - Added velocity attenuation.  Activate with the flag "-va".
!
! 02/21/2022, RC, Version 0.4.0:
!   - Added ability to output maps at a specified iteration cadence.
!     Use the new flag "-omci" to set the Output Map Cadence Index.
!     A list of all output is stored in an output text file
!     called "map_output_list.txt".
!     The output maps are stored in a subdirectory named "output_maps".
!     If omci is set to 0 or not set, the subdirectory and text file are
!     not created.
!     The initial and final map files are still always output in the
!     run directory as before.
!   - Some internal updates and simplifications.
!
! 03/21/2022, RC, Version 0.5.0:
!   - Added data assimilation.  Activate with the flag "-assimilate_data"
!     and set the flag "-assimilate_data_map_list_filename" to the name
!     of the CSV file obtained from the data aquisition tool pyoft.
!     It is assumed that the code starts with the first entry of the
!     CSV file at time=0.
!   - Added time cadence output option.  Set with the flag "-omct <DT>"
!     where <DT> is the time cadense requested for output.
!     Note that one cannot use both -omct and -omci in the same run.
!   - Updated timestep setting so that outputs and data assimilation are
!     always at exact times.
!   - Updated I/O routines, including adding a 3D reader.
!   - DEFAULT CHANGE: Updated dt_max_increase_fac to 0.0 (disabled).
!   - Updated text output of the run, including more information when
!     using the verbose "-v" flag.
!
! 04/07/2022, RC, Version 0.6.0:
!   - Added option to use Strang splitting in ARDRA form.
!     This will keep the accuracy second-order in time once
!     the advection is updated to a higher-order scheme.
!     To use, set "-strang" on the command line.
!   - Added option to use simple pole averaging in advection.
!     set -flowpoleavg to use.
!   - Added -output_map_directory flag to set output directory.
!   - Added -assimilate_data_map_root_dir flag to set the base
!     directory for the input assimilation data.
!   - Fixed polar diffusion coef.
!   - Small robustness updates.
!
! 04/11/2022, RC, Version 0.7.0:
!   - Added the RK3TVD time-stepping scheme for advection.
!     This is now the default.
!   - DEFAULT CHANGE: Made the default splitting method to be Strang
!                     splitting. The "-strang" input flag is now
!                     "-nostrang" and disables strang splitting.
!   - Added "-fm" input flag to set advection numerical method.
!     1) FE+UW  2) RK3TVD+UW  3) RK3TVD+WENO3 (coming soon).
!   - Removed "-flowpoleavg".  No use case.
!
! 04/29/2022, RC, Version 0.7.1:
!   - BUG FIX: Changed staggering of the velocity.  Now, Vt is on the
!              theta half mesh and the phi main mesh, while Vp is on the
!              theta main mesh and the phi half mesh.  This avoids averaging
!              in the advection operator (which was not being done).
!   - BUG FIX: The diffusion_coef was not being averaged.
!   - BUG FIX: Fixed issue when not outputting time-dept I/O.
!
! 05/17/2022, RC+MS, Version 0.8.0:
!   - Added implementation of the WENO3 advection scheme (-fm 3).
!   - BUG FIX: time-step calculation.
!
! 05/20/2022, RC+MS, Version 0.9.0:
!   - vrun has been changed to an integer to select which validation
!     run initial condition to use.
!   - Added Guassian blob initial validation condition (-vrun 2).
!   - Added -vt_const to allow adding a constant velocity in theta
!     in km/s (used for validation runs).
!   - BUG FIX: GPU issue with -vrun 1 initial condition.
!
! 06/07/2022, RC+MS, Version 0.10.0:
!   - BUG FIX: WENO3 coef fix for non-uniform grids.
!   - Added new validation run IC.
!   - Added new rotated rigid velocity option for testing purposes.
!
! 06/20/2022, RC, Version 0.11.0:
!   - Added additional diagnostics in history output:
!     - Northern and southern polar fluxes (+ & -) and areas.
!     - Equatorial and axial dipole strengths.
!   - Made fluxes in history output in Mx.
!
! 07/15/2022, RC, Version 0.12.0:
!   - Added time dependent flows from file.
!     Activate with input flag "-use_flow_from_files".
!     Must set flag "-assimilate_data_map_root_dir" to the location
!     of the flow files, and "-flow_list_filename" to the location
!     of the csv file listing all files.  This file must be in a
!     specific format which is currently output by the ConFlow code
!     (see there for details). Note that "-flow" still needs to be active.
!     The flow files for vt and vp are assumed to be on the correct
!     staggered grids.
!   - Removed vtfile and vpfile input options.  To use a static
!     file-based flow, use the new time-dept flow_from_files with
!     2 identical files.
!   - DEFAULT CHANGE: Changed polar latitude limit for
!                     analysis from 25 to 30 degrees.
!   - BUG FIX: Data asssimilation was using "requested times" instead
!              of "actual times".
!   - BUG FIX: Fixed units of axial/equatorial dipole analysis output.
!
! 07/25/2022, RC, Version 0.12.1:
!   - BUG FIX:  Diffusion STS factors were not being recomputed when
!               the time step changed.
!
! 07/28/2022, RC, Version 0.13.0:
!   - Updated data assimilation loading to reflect changes in
!     OFTpy's file list format.
!   - BUG FIX: Fixed integral in axial and equatorial dipole strength.
!
! 11/01/2022, MS+RC, Version 0.14.0:
!   - Changed input to namelist instead of command line flags.
!     See the namelist in the code for parameter names.
!   - Added multiple realizations.  Use n_realizations to set.
!     Currently, each realization is an identical computation.
!     If using 1 realization, use output_map_2d to output 2D file
!     instead of 3d file (currently 2d is default).
!
! 11/22/2022, MS+RC, Version 0.15.0:
!   - Added MPI to parallelize across realizations.
!
!-----------------------------------------------------------------------
!
!#######################################################################
