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
! ****** HipFT: High-performance Flux Transport.
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
      character(*), parameter :: cvers='1.11.0'
      character(*), parameter :: cdate='07/05/2024'
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
      integer(8),  parameter :: one_int = 1
      real(r_typ), parameter :: two = 2.0_r_typ
      integer(8),  parameter :: two_int = 2
      real(r_typ), parameter :: three = 3._r_typ
      integer(8),  parameter :: three_int = 3
      real(r_typ), parameter :: four = 4._r_typ
      integer(8),  parameter :: four_int = 4
      real(r_typ), parameter :: six = 6._r_typ
      real(r_typ), parameter :: sixth = 0.16666666666666666_r_typ
      real(r_typ), parameter :: nine = 9._r_typ
      real(r_typ), parameter :: ten = 10._r_typ
      real(r_typ), parameter :: fifteen = 15._r_typ
      real(r_typ), parameter :: sixteen = 16._r_typ
      real(r_typ), parameter :: half = 0.5_r_typ
      real(r_typ), parameter :: quarter = 0.25_r_typ
      real(r_typ), parameter :: twentyfour = 24.0_r_typ
      real(r_typ), parameter :: twentyfive = 25.0_r_typ
      real(r_typ), parameter :: three_quarter = 0.75_r_typ
      real(r_typ), parameter :: two_third = 0.66666666666666666_r_typ
      real(r_typ), parameter :: third = 0.33333333333333333_r_typ
!
      real(r_typ), parameter :: pi = 3.1415926535897932_r_typ
      real(r_typ), parameter :: pi_two = 1.5707963267948966_r_typ
      real(r_typ), parameter :: pi_four = 0.7853981633974483_r_typ
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
      real(r_typ), parameter :: input_flux_fac = 1.0e21_r_typ
!
      real(r_typ), parameter :: small_value = tiny(one)
      real(r_typ), parameter :: large_value = huge(one)
      real(r_typ), parameter :: safety = 0.95_r_typ
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
      logical :: output_current_history = .false.
!
! ****** Current output file sequence number.
!
      integer(8) :: output_seq = 0
!
! ****** File sequence number.
!
      integer(8) :: idx_out = 0
      integer(8) :: idx_hist = 0
!
      integer, parameter :: IO_HIST_NUM = 20
      integer, parameter :: IO_HIST_SOL = 21
      integer, parameter :: IO_MAP_OUT_LIST = 22
      integer, parameter :: IO_TMP = 23
!
      integer, parameter :: IO_DATA_IN = 8
!
      character(29) :: io_hist_num_filename = 'hipft_history_num.out'
      character(29) :: io_hist_sol_filename = 'hipft_history_sol.out'
      character(25) :: io_output_map_list_filename = 'hipft_output_map_list.out'
      character(26) :: io_output_flows_list_filename = 'hipft_output_flow_list.out'
!
! ****** Temporary arrays used for output.
!
      real(r_typ), dimension(:,:,:), allocatable :: fout
      real(r_typ), dimension(:), allocatable :: pout
!
! ****** History analysis variables.
!
      real(r_typ), dimension(:), allocatable :: h_minbr,     h_maxbr,    &
                                                h_minabsbr,  h_fluxp,    &
                                                h_fluxm,     h_valerr,   &
                                                h_fluxp_pn,  h_fluxm_pn, &
                                                h_fluxp_ps,  h_fluxm_ps, &
                                                h_area_pn,   h_area_ps,  &
                                                h_eq_dipole, h_ax_dipole
end module
!#######################################################################
module globals
!
!-----------------------------------------------------------------------
! ****** Internal general global variables.
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
! ****** Main Br field arrays
!
      real(r_typ), dimension(:,:,:), allocatable :: f,fold,fval_u0
!
! ****** Main terminal output format string.
!
      character(*), parameter :: &
        MAINFMT = '(a,i10,a,f14.3,a,f14.3,a,f12.6,a,a,a)'
!
! ****** Current time.
!
      real(r_typ) :: time = 0.
!
! ****** Current step.
!
      integer(8) :: ntime = 0
!
! ****** Current time-step.
!
      real(r_typ) :: dtime_global = 0.
      character(50) :: dtime_reason = ''
!
! ****** Flag to indicate the time step needs updating.
!
      logical :: timestep_needs_updating = .true.
!
! ****** Internal realization variables.
!
      integer, parameter :: MAX_REALIZATIONS = 2000
!
      integer, dimension(:), allocatable :: local_realization_indices
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
module diffusion
!
      use number_types
!
      implicit none
!
      real(r_typ), dimension(:,:,:), allocatable :: diffusion_coef
!
! ****** Diffusion constant coef array.
!
      real(r_typ), dimension(:), allocatable :: diffusion_coef_constant_rvec
!
! ****** Explicit Euler diffusion stable time-step and number of cycles.
!
      real(r_typ) :: dtime_diffusion_euler = 0.
      integer(8) :: n_stable_diffusion_cycles = 1
      logical :: auto_sc=.false.
!
! ****** Flag to indicate the diffusion time step needs updating.
!
      logical :: timestep_diff_needs_updating = .true.
!
end module
!#######################################################################
module sources
!
      use number_types
!
      implicit none
!
      real(r_typ), dimension(:,:,:), allocatable :: source
      real(r_typ), dimension(:,:,:), allocatable :: source_from_file_data
      real(r_typ), dimension(:,:,:,:), allocatable :: source_rfe
!
      integer :: rfe1 = 1
      integer :: rfe2 = 2
!
      real(r_typ) :: time_of_next_source_rfe
      real(r_typ) :: time_of_prev_source_rfe
      real(r_typ), dimension(:), allocatable :: uf_per_cell_per_hour_rvec
      real(r_typ), dimension(:), allocatable :: uf_per_area_per_hour_rvec
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
      integer, dimension(:), allocatable :: assimilate_data_lat_limit_tidx0_rvec
      integer, dimension(:), allocatable :: assimilate_data_lat_limit_tidx1_rvec
!
      real(r_typ), dimension(:), allocatable :: assimilate_data_lat_limit_rvec
      real(r_typ), dimension(:), allocatable :: assimilate_data_mu_limit_rvec
      real(r_typ), dimension(:), allocatable :: assimilate_data_mu_power_rvec
!
      real(r_typ) :: map_time_initial_hr
      real(r_typ), dimension(:), allocatable :: map_times_actual_ut_jd
      character(19), dimension(:), allocatable :: &
                     map_times_requested_ut_str, map_times_actual_ut_str
      character(512), dimension(:), allocatable :: map_files_rel_path
!
end module
!#######################################################################
module flows
!
      use number_types
!
      implicit none
!
! ****** Velocity fields.
!
      real(r_typ), dimension(:,:,:), allocatable :: vt
      real(r_typ), dimension(:,:,:), allocatable :: vp
!
! ****** Explicit Euler advection stable time-step.
!
      real(r_typ) :: dtmax_flow = 0.
      real(r_typ) :: dtime_advection_used = 0.
!
! ****** Flag to indicate the flow needs updating.
!
      logical :: flow_needs_updating = .true.
!
! ****** Flow attenuation.
!
      real(r_typ), dimension(:), allocatable :: flow_attenuate_value_i_rvec
!
! ****** Flow differential rotational.
!
      real(r_typ), dimension(:), allocatable :: flow_dr_coef_p0_rvec
      real(r_typ), dimension(:), allocatable :: flow_dr_coef_p2_rvec
      real(r_typ), dimension(:), allocatable :: flow_dr_coef_p4_rvec
!
! ****** Flow meridianal.
!
      real(r_typ), dimension(:), allocatable :: flow_mf_coef_p1_rvec
      real(r_typ), dimension(:), allocatable :: flow_mf_coef_p3_rvec
      real(r_typ), dimension(:), allocatable :: flow_mf_coef_p5_rvec
!
! ****** Flow polar boundary condition factors.
!
      real(r_typ) :: bc_flow_npole_fac, bc_flow_spole_fac
!
! ****** Flow from files.
!
      integer :: flow_from_files_num_flows = 0
      integer :: flow_from_file_current_idx = 1
      real(r_typ) :: time_of_next_input_flow = 0.
      real(r_typ) :: time_of_current_input_flow = 0.
      real(r_typ) :: flow_from_files_dt_of_last_flows = 0.
!
      real(r_typ), dimension(:), allocatable :: flow_from_files_times
!
      character(512), dimension(:), allocatable :: flow_from_files_rel_path_t
      character(512), dimension(:), allocatable :: flow_from_files_rel_path_p
!
      real(r_typ), dimension(:,:,:),allocatable :: flow_from_files_current_vt
      real(r_typ), dimension(:,:,:),allocatable :: flow_from_files_current_vp
!
! ****** Numerical methods for flow advance.
!
      integer :: advection_num_method_space
      integer :: advection_num_method_time
!
      integer, parameter :: FE      = 1
      integer, parameter :: SSPRK33 = 2
      integer, parameter :: SSPRK43 = 3
!
      integer, parameter :: UW      = 1
      integer, parameter :: WENO3   = 2
!
! ****** Flag to indicate the stable flow time step needs updating.
!
      logical :: timestep_flow_needs_updating = .true.
!
end module
!#######################################################################
module sts
!
      use number_types
!
      implicit none
!
      integer(8) :: sts_s
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
      real(r_typ) :: wtime_data_assimilation = 0.
      real(r_typ) :: wtime_analysis = 0.
      real(r_typ) :: wtime_io = 0.
      real(r_typ) :: wtime_mpi_overhead = 0.
!
      real(r_typ) :: wtime_total = 0.
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
      use globals, ONLY : MAX_REALIZATIONS
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer, private :: i
!
      integer :: verbose = 0
!
! ****** Resolution ********
! ****** These are autoset when reading in an initial_map_filename *****
! ****** They are currently only used for validation runs.
!
      integer :: res_nt = 0
      integer :: res_np = 0
!
! ****** Number of realizations ********
!
      integer        :: n_realizations = 1
!
! ****** Validation Mode ********
!
      integer        :: validation_run = 0
      real(r_typ)    :: validation_run_width = 0.03_r_typ
!
! ****** Initial map ********
!
      character(512) :: initial_map_filename = ''
      logical :: initial_map_flux_balance = .false.
!
! ****** Output options ********
!
      character(512) :: output_map_root_filename = 'hipft_brmap'
      character(512) :: output_map_directory     = 'output_maps'
      logical        :: output_map_flux_balance = .false.
      integer(8)     :: output_map_idx_start = 1
      logical        :: output_flows = .false.
      character(512) :: output_flows_root_filename = 'hipft_flow'
      character(512) :: output_flows_directory     = 'output_flows'
      integer        :: output_map_idx_cadence = 0
      real(r_typ)    :: output_map_time_cadence = zero
      logical        :: output_map_2d = .true.
      real(r_typ)    :: output_history_time_cadence = zero
      logical        :: output_single_precision = .true.
!
! ****** Time ********
!
      real(r_typ)    :: time_start = zero
      real(r_typ)    :: time_end = one
!
! ****** Timestep ********
!
      real(r_typ)    :: dt_min = 1.0e-15_r_typ
      real(r_typ)    :: dt_max = huge(one)
!
! ****** General algorithm options.
!
      logical :: strang_splitting = .false.
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
      real(r_typ) :: flow_attenuate_value = 500.0_r_typ
      real(r_typ), dimension(MAX_REALIZATIONS) :: flow_attenuate_values
      data (flow_attenuate_values(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1./
!
! ****** Analytic differential roation and meridianal flow models.
!
      integer :: flow_dr_model = 0
      integer :: flow_mf_model = 0
!
      real(r_typ) :: flow_dr_coef_p0_value = 46.0_r_typ
      real(r_typ) :: flow_dr_coef_p2_value = -262.0_r_typ
      real(r_typ) :: flow_dr_coef_p4_value = -379.0_r_typ
      real(r_typ), dimension(MAX_REALIZATIONS) :: flow_dr_coef_p0_values
      real(r_typ), dimension(MAX_REALIZATIONS) :: flow_dr_coef_p2_values
      real(r_typ), dimension(MAX_REALIZATIONS) :: flow_dr_coef_p4_values
      data (flow_dr_coef_p0_values(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1000000./
      data (flow_dr_coef_p2_values(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1000000./
      data (flow_dr_coef_p4_values(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1000000./
!
      real(r_typ) :: flow_mf_coef_p1_value = 22.0_r_typ
      real(r_typ) :: flow_mf_coef_p3_value = 11.0_r_typ
      real(r_typ) :: flow_mf_coef_p5_value = -28.0_r_typ
      real(r_typ), dimension(MAX_REALIZATIONS) :: flow_mf_coef_p1_values
      real(r_typ), dimension(MAX_REALIZATIONS) :: flow_mf_coef_p3_values
      real(r_typ), dimension(MAX_REALIZATIONS) :: flow_mf_coef_p5_values
      data (flow_mf_coef_p1_values(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1000000./
      data (flow_mf_coef_p3_values(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1000000./
      data (flow_mf_coef_p5_values(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1000000./
!
! ****** Time dependent flows from file (used for conflow).
!
      logical :: use_flow_from_files = .false.
      character(512) :: flow_list_filename = ' '
      character(512) :: flow_root_dir = '.'
!
! ****** Algorithm options.
!        Can set upwind to central differencing by setting UPWIND=0.
! ****** 1: Forward Euler + Upwind.
! ****** 2: SSPRK(3,3) + Upwind.
! ****** 3: SSPRK(3,3) + WENO3.
! ****** 4: SSPRK(4,3) + WENO3
!
      integer :: flow_num_method = 4
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
      real(r_typ), dimension(MAX_REALIZATIONS) :: diffusion_coef_constants
      data (diffusion_coef_constants(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1./
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
! ****** For diffusion-only runs with RK[G|L]2, robust to set to ~30(60).
! ****** For flow+diffusion runs, this usually can be ~1.
! ****** Set this to 0 to "auto" set the subcycles.
!
      integer :: diffusion_subcycles = 0
!
!-----------------------------------------------------------------------
!
! ****** SOURCES ********
!
      logical :: advance_source = .false.
!
      integer     :: source_rfe_model = 0
      real(r_typ) :: source_rfe_total_unsigned_flux_per_hour = 0.
      real(r_typ), dimension(MAX_REALIZATIONS) :: source_rfe_total_unsigned_flux_per_hours
      data (source_rfe_total_unsigned_flux_per_hours(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1./
      integer     :: source_rfe_seed = -9999
      real(r_typ) :: source_rfe_lifetime = 0.3_r_typ
!
      logical :: source_from_file = .false.
      character(512) :: source_filename = ' '
!
!-----------------------------------------------------------------------
!
! ****** DATA ASSIMILATION ********
!
      logical :: assimilate_data = .false.
      logical :: assimilate_data_balance_flux = .false.
!
      character(512) :: assimilate_data_map_list_filename = ' '
      character(512) :: assimilate_data_map_root_dir = '.'
!
! ****** Custom assimilation options.
!
      logical :: assimilate_data_custom_from_mu = .false.
!
      real(r_typ) :: assimilate_data_lat_limit = 0.
      real(r_typ), dimension(MAX_REALIZATIONS) :: assimilate_data_lat_limits
      data (assimilate_data_lat_limits(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1./
!
! ****** Mu parameters.
!
      real(r_typ) :: assimilate_data_mu_power = 4.0_r_typ
      real(r_typ), dimension(MAX_REALIZATIONS) :: assimilate_data_mu_powers
      data (assimilate_data_mu_powers(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1./
!
      real(r_typ) :: assimilate_data_mu_limit = 0.1_r_typ
      real(r_typ), dimension(MAX_REALIZATIONS) :: assimilate_data_mu_limits
      data (assimilate_data_mu_limits(i),i=1,MAX_REALIZATIONS) /MAX_REALIZATIONS*-1./
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
      integer, dimension(:), allocatable :: nr_arr, nr_disp_arr
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
      integer :: iproc
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
      use iso_fortran_env,  ONLY : OUTPUT_UNIT
      use globals,          ONLY : MAINFMT, ntime, time,              &
                                   dtime_global, dtime_reason
      use input_parameters, ONLY : time_end, verbose
      use mpidefs,          ONLY : iamp0, MPI_Wtime
      use timing,           ONLY : wtime_total
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
! ****** Set up and initialize the run.
!
      call setup
!
      if (iamp0) then
        write (*,*)
        write (*,*) '>running'
        write (*,*)
        FLUSH(OUTPUT_UNIT)
      end if
!
! ****** START MAIN LOOP ********************************************
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
! ****** Evolve the field for one timestep.
!
        call flux_transport_step
!
! ****** Update the time now that we have taken the step.
!
        time = time + dtime_global
!
! ****** Apply data assimilation. Since this step replaces data directly,
! ****** we apply is at the end of the step advance.
! ****** This avoids double random flux (model RFE + data RFE).
!
        call apply_data_assimilation
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
          if (verbose.gt.0) write(*,*) ' '
          write(*,MAINFMT)                                            &
            'Completed step: ',ntime,'  Time: ',time,' / ',           &
            time_end,'  dt: ',dtime_global,' (',TRIM(dtime_reason),')'
          FLUSH(OUTPUT_UNIT)
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
! ****** Write out time profile.
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
      integer :: i, idisp, nravg, nrem
!
!-----------------------------------------------------------------------
!
! ****** Set up arrays of local number of realizations and displacements.
!
      allocate (nr_arr(0:nproc-1))
      allocate (nr_disp_arr(0:nproc-1))

      nravg = FLOOR(REAL(n_realizations,r_typ)/REAL(nproc,r_typ))
      nr_arr(:) = nravg
!
      if (nravg*nproc .lt. n_realizations) then
        write(*,*)
        write(*,*) '   WARNING: Load imbalance detected.'
        write(*,*)
        nrem = n_realizations - nravg*nproc
        do i=0,nrem-1
          nr_arr(i) = nr_arr(i) + 1
        enddo
      end if
!
      idisp = 0
      do i=0,nproc-1
        nr_disp_arr(i) = idisp
        idisp = idisp + nr_arr(i)
      enddo
!
      if (iamp0 .and. verbose.gt.0) then
        write(*,*) ' '
        write(*,*) '   SET_REALIZATIONS: Realization distribution:'
        write(*,*) '   SET_REALIZATIONS: nr_arr',nr_arr(:)
        write(*,*) '   SET_REALIZATIONS: nr_disp_arr',nr_disp_arr(:)
      end if
!
! ****** Now set local number of realizations.
!
      nr = nr_arr(iproc)
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
! ****** communicator MPI_COMM_WORLD in variable IPROC.
!
      call MPI_Comm_rank (MPI_COMM_WORLD,iproc,ierr)
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
      if (iproc .eq. 0) then
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
! ****** We have to put an OpenACC device selection directive here,
! ****** as the NVIDIA compiler needs it for the device selection
! ****** for DC loops, while the omp device selection is used for
! ****** OpenMP target data movement directives.
!
!$ call omp_set_default_device (iprocsh)
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
      use mesh
      use output
      use globals
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
! ****** Initialize MPI.
!
      call initialize_mpi
!
      wtime_tmp = MPI_Wtime()
!
      if (iamp0) then
        call write_welcome_message
      end if
!
      if (iamp0) write(*,*) '--> Reading input file...'
      call read_input_file
!
! ****** Initialize realization arrays.
!
      if (iamp0) write(*,*) '--> Initializing realization parameters...'
      call set_realization_parameters
!
! ****** Set up output directory.
!
      if (output_map_idx_cadence.gt.0 .or.     &
          output_map_time_cadence.gt.0.) then
        if (iamp0) write(*,*) '--> Setting up output directories...'
!
        if (iamp0) then
          cmd = 'mkdir -p '//trim(output_map_directory)
          call EXECUTE_COMMAND_LINE (cmd,exitstat=ierr)
          if (ierr.ne.0) then
            write (*,*) '### Could not make map output subdirectory, using ./'
            output_map_directory = '.'
          end if
          if (output_flows) then
            cmd = 'mkdir -p '//trim(output_flows_directory)
            call EXECUTE_COMMAND_LINE (cmd,exitstat=ierr)
            if (ierr.ne.0) then
              write (*,*) '### Could not make flow output subdirectory, using ./'
              output_flows_directory = '.'
            end if
          end if
        end if
!
! ****** Make sure initial conditon is output.
!
        output_current_map = .true.
!
      end if
!
      if (iamp0) write(*,*) '--> Loading initial condition...'
      call load_initial_condition
!
      time = time_start
!
! ****** Write out some basic run info.
!
      if (iamp0.and.verbose.gt.0) then
        write (*,*) ''
        write (*,*) '    Start time: ',time
        write (*,*) '    End time:   ',time_end
        write (*,*) ''
        write (*,*) '    Grid Nt: ',ntm
        write (*,*) '    Grid Np: ',npm-1
        write (*,*) ''
        write (*,*) '    Number of realizations: ',nr
        write (*,*) ''
      endif

      if (advance_flow) then
        if (iamp0) write(*,*) '--> Loading flows...'
        call load_flows
      end if
!
      if (advance_diffusion) then
        if (iamp0) write(*,*) '--> Loading diffusion...'
        call load_diffusion
      end if
!
      if (advance_source) then
        if (iamp0) write(*,*) '--> Loading sources...'
        call load_source
      end if
!
      if (assimilate_data)  then
        if (iamp0) write(*,*) '--> Loading data assimilation...'
        call load_data_assimilation
      end if
!
      if (iamp0) write(*,*) '--> Creating history and log files...'
      call create_and_open_output_log_files
!
      if (iamp0) write(*,*) '--> Analyzing initial condition...'
      call analysis_step
!
      if (iamp0) write(*,*) '--> Initial output...'
      call output_step
!
      if (iamp0) then
        write (*,*)
        write (*,*) '--> Setup complete!'
        write (*,*)
      end if
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
      logical :: the_run_is_done, check_if_stoprun
!
!-----------------------------------------------------------------------
!
      the_run_is_done = .false.
!
! ****** Check if the simulation time has reached the end time.
!
      if (time .ge. time_end) then
        the_run_is_done = .true.
        return
      end if
!
! ****** Check if a STOPRUN file has been created.
!
      if (check_if_stoprun()) then
        the_run_is_done = .true.
        return
      end if
!
      return
end function
!#######################################################################
function check_if_stoprun ()
!
!-----------------------------------------------------------------------
!
! ****** Stop the present run if a file named "STOPRUN" exists in
! ****** the run directory.
!
! ****** If it does, return true.
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
      logical :: check_if_stoprun
      logical :: exists = .false.
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      check_if_stoprun = .false.
!
! ****** Check if the user has requested to stop the run
! ****** prematurely.
!
      if (iamp0) then
        inquire (file='STOPRUN',exist=exists)
        if (exists) then
          write (*,*)
          write (*,*) '### NOTE: The run was stopped prematurely via STOPRUN!'
          check_if_stoprun = .true.
        end if
      end if
!
! ****** Broadcast check_if_stoprun to all the processors.
!
      call MPI_Bcast (check_if_stoprun,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
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
        write(*,*) 'Final data written out to ', &
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
      integer, parameter :: lbuf=11
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
      sbuf(11) = wtime_data_assimilation
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
      tsdev(:) = 0.
!
      do irank=0,nproc-1
        tsdev(:) = tsdev(:) + (tbuf(:,irank) - tavg(:))**2
      enddo
!
      tsdev(:) = sqrt(tsdev(:)/nproc)
!
      if (iamp0) then
!
        call ffopen (8,'hipft_timing.out','rw',ierr)
!
        do irank=0,nproc-1
          write(8,"(a40)") repeat("-", 40)
          write(8,*)
          write(8,*) 'Processor id = ',irank
          write(8,*)
          write(8,"(a40)") repeat("-", 40)
          write(8,FMT) "Wall clock time:   ",tbuf(1,irank)
          write(8,"(a40)") repeat("-", 40)
          write(8,FMT) "--> Setup:         ",tbuf(2,irank)
          write(8,FMT) "--> Update:        ",tbuf(3,irank)
          write(8,FMT) "--> Data Assim:    ",tbuf(11,irank)
          write(8,FMT) "--> Flux transport:",tbuf(4,irank)
          write(8,FMT) "    --> Advecton:  ",tbuf(5,irank)
          write(8,FMT) "    --> Diffusion  ",tbuf(6,irank)
          write(8,FMT) "    --> Source:    ",tbuf(7,irank)
          write(8,FMT) "--> Analysis:      ",tbuf(8,irank)
          write(8,FMT) "--> I/O:           ",tbuf(9,irank)
          write(8,FMT) "--> MPI overhead:  ",tbuf(10,irank)
          write(8,"(a40)") repeat("-", 40)
          write(8,*)
        enddo
!
        write (8,*) 'Run time summary:'
        write (8,*) '-----------------'
        write (8,*)
        write (8,300) 'Avg         Min         Max      S. Dev'
        write (8,300) '---         ---         ---      ------'
        write (8,400) 'Wall clock time:        ', &
                     tavg(1),tmin(1),tmax(1),tsdev(1)
        write (8,400) '--> Setup:              ', &
                     tavg(2),tmin(2),tmax(2),tsdev(2)
        write (8,400) '--> Update:             ', &
                     tavg(3),tmin(3),tmax(3),tsdev(3)
        write (8,400) '--> Data Assimilation:  ', &
                     tavg(11),tmin(11),tmax(11),tsdev(11)
        write (8,400) '--> Flux transport:     ', &
                     tavg(4),tmin(4),tmax(4),tsdev(4)
        write (8,400) '    --> Advecton:       ', &
                     tavg(5),tmin(5),tmax(5),tsdev(5)
        write (8,400) '    --> Diffusion       ', &
                     tavg(6),tmin(6),tmax(6),tsdev(6)
        write (8,400) '    --> Source:         ', &
                     tavg(7),tmin(7),tmax(7),tsdev(7)
        write (8,400) '--> Analysis:           ', &
                     tavg(8),tmin(8),tmax(8),tsdev(8)
        write (8,400) '--> I/O:                ', &
                    tavg(9),tmin(9),tmax(9),tsdev(9)
        write (8,400) '--> MPI overhead:       ', &
                    tavg(10),tmin(10),tmax(10),tsdev(10)
        write (8,*)
        write(8,"(a40)") repeat("-", 40)
        write (8,*)
!
        close (8)
!
        write (*,"(a40)") repeat("-", 40)
        write (*,*)
        write (*,*) 'Run time summary:'
        write (*,*) '-----------------'
        write (*,*)
        write (*,300) 'Avg         Min         Max      S. Dev'
        write (*,300) '---         ---         ---      ------'
        write (*,400) 'Wall clock time:        ', &
                     tavg(1),tmin(1),tmax(1),tsdev(1)
        write (*,400) '--> Setup:              ', &
                     tavg(2),tmin(2),tmax(2),tsdev(2)
        write (*,400) '--> Update:             ', &
                     tavg(3),tmin(3),tmax(3),tsdev(3)
        write (*,400) '--> Data Assimilation:  ', &
                     tavg(11),tmin(11),tmax(11),tsdev(11)
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
        write (*,*)
        write(*,"(a40)") repeat("-", 40)
        write (*,*)
!
        flush(OUTPUT_UNIT)
      end if
!
  300 format (1x,33x,a)
  400 format (1x,a,4f12.3)
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
      use constants
      use mpidefs
      use mesh
      use globals
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr, i
!
!-----------------------------------------------------------------------
!
      do i=1,nr
!
        write(io_hist_num_filename,'(A19,I6.6,A4)') &
              "hipft_history_num_r", local_realization_indices(i), ".out"
!
        write(io_hist_sol_filename,'(A19,I6.6,A4)') &
              "hipft_history_sol_r", local_realization_indices(i), ".out"
!
        call ffopen (IO_HIST_NUM,io_hist_num_filename,'rw',ierr)
!
        write (IO_HIST_NUM,'(a10,a,7(a22,a))') &
        'STEP',' ',&
        'TIME',' ',&
        'DTIME',' ',&
        'DTIME_ADV_STB',' ',&
        'DTIME_ADV_USED',' ',&
        'DTIME_DIFF_STB',' ',&
        'N_DIFF_PER_STEP',' ',&
        'N_DIFF_CYCLES',' '
!
        close(IO_HIST_NUM)
!
        call ffopen (IO_HIST_SOL,io_hist_sol_filename,'rw',ierr)
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
        'VALIDATION_ERR_CVRMSD',' '
!
        close(IO_HIST_SOL)
!
      enddo
!
! ****** Only rank 0 writes the IO logs.
!
      if (iamp0 .and. ( output_map_idx_cadence.gt.zero .or.     &
                       output_map_time_cadence.gt.zero)) then
!
        call ffopen (IO_MAP_OUT_LIST,io_output_map_list_filename,&
                     'rw',ierr)
!
        write (IO_MAP_OUT_LIST,'(a10,a,a25,a,a)') &
                               'IDX',' ',&
                               'TIME',' ',&
                               'FILENAME'
!
        close(IO_MAP_OUT_LIST)

        if (output_flows) then
!
          call ffopen (IO_MAP_OUT_LIST,io_output_flows_list_filename,&
                       'rw',ierr)
!
          write (IO_MAP_OUT_LIST,'(a10,a,a25,a,a)') &
                                 'IDX',' ',&
                                 'TIME',' ',&
                                 'FILENAME'
!
          close(IO_MAP_OUT_LIST)
!
        end if
!
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
      use output
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
      write (*,*) '     |  __  | | ''_ \\|  __|    | |'
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
      write (*,*) ''
      write (*,*) "Run started at: "
      write (*,*) ''
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
      use globals
      use read_2d_file_interface
      use output
      use diffusion
      use constants
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n1, n2, n3, ierr, i, j, k, lsbuf, irank
      integer, dimension(:), allocatable :: displ,lbuf
      real(r_typ) :: d1,d2, st_i_val
      real(r_typ) :: vt_end_loc
      real(r_typ), dimension(:), allocatable :: fn1,fs1,fn2,fs2
      real(r_typ), dimension(:,:,:), allocatable :: f_local,gout
      real(r_typ), dimension(:,:), allocatable :: f_tmp2d
      real(r_typ), dimension(:), allocatable :: s1,s2,s3,s3t,tfout,gfout
!
!-----------------------------------------------------------------------
!
! ****** NOTE: If the user does not specify NT,NP,
! ****** the resolution of the input map is used.
! ****** The input map is assumed to be in PT and with a
! ****** one-point overlap in phi.
! ******
!
! ****** If we are starting with an empty uniform grid map, allocate
! ****** and set the scales.
!
! [RMC] Eventaully, we will allow interpolation, so can set res and
!       have a filename.
!
      if (initial_map_filename .eq. '') then

        n1 = res_np
        d1 = twopi/(n1-1)
        n2 = res_nt
        d2 = pi/(n2-1)
        n3 = nr

        allocate (f_local(n1,n2,n3))
        allocate (f_tmp2d(n1,n2))
!
! ****** Can add other inital conditon options here.
!
        f_local(:,:,:) = 0.
        f_tmp2d(:,:) = 0.

        allocate (s1(n1))
        do i=1,n1
          s1(i) = (i-1)*d1
        enddo

        allocate (s2(n2))
        do i=1,n2
          s2(i) = (i-1)*d2
        enddo

        allocate (s3(n3))
        do i=1,n3
          s3(i) = real(i,r_typ)
        enddo

      elseif (initial_map_filename .ne. '') then
!
! ****** Read the initial magnetic map.
!
        if (iamp0) then
          call read_2d_file (initial_map_filename,n1,n2,f_tmp2d,s1,s2,ierr)
        endif
!
        wtime_tmp_mpi = MPI_Wtime()
        call MPI_Bcast (n1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (n2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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
        n3 = nr
        allocate (f_local(n1,n2,n3))

        do i=1,n3
          f_local(:,:,i) = f_tmp2d(:,:)
        enddo
!
        allocate (s3(n3))
        do i=1,n3
          s3(i) = real(i,r_typ)
        enddo
!
      end if
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
!$omp target enter data map(alloc:f_local,fval_u0,f_tmp2d)
!$omp target enter data map(to:s1,s2)
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
!$omp target update from(f_local,fval_u0)
!$omp target exit data map(delete:f_local,fval_u0,f_tmp2d,s1,s2)
!
        fval_u0(:,:,:) = f_local(:,:,:) +                             &
                         sqrt(14.0_r_typ/11.0_r_typ)*fval_u0(:,:,:)
!
        if (output_map_2d.and.n_realizations.eq.1) then
          call write_2d_file(                                         &
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
        end if
!
        f_local(:,:,:) = fval_u0(:,:,:)
!
! ****** Caculate final analytic solution and output.
!
        allocate (fout(n1,n2,n3))
        do k=1,n3
          do j=1,n2
            do i=1,n1
              fout(i,j,k) = fval_u0(i,j,k)*exp(-42.0_r_typ*diffusion_coef_constant_rvec(k)* &
                                           diffusion_coef_factor*time_end)
            enddo
          enddo
        enddo
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
!$omp target enter data map(to:fval_u0)
!
      elseif (validation_run .eq. 2) then
!
! ****** Make initial solution f and output.
!
        allocate (fout(n1,n2,1))
!
        if (iamp0) then
          write(*,*) ' '
          write(*,*) '   VALIDATION 2 TEST ENABLED.'
        end if
!
! ****** Use time to set constant theta velocity so final output is exact solution.
! ****** for a constant theta velocity.
!
        if (ABS(flow_vt_const).gt.zero) then
          if (iamp0) then
            write(*,*) ' '
            write(*,*) '   WARNING:  Overriding FLOW_VT_CONST'
            write(*,*) '             from',flow_vt_const,' to ', &
                        (pi_two/(time_end-time_start))/km_s_to_rs_hr
          end if
          flow_vt_const = (pi_two/(time_end-time_start))/km_s_to_rs_hr
          vt_end_loc = threepi_four
        else
          vt_end_loc = pi_four
        end if
!
! ****** Use time to set constant phi velocity so final output is exact solution.
! ****** for a constant phi rigid rotation velocity.
!
        if (ABS(flow_vp_rigid_omega).gt.zero) then
          if (iamp0) then
            write(*,*) ' '
            write(*,*) '   WARNING:  Overriding FLOW_VP_RIGID_OMEGA'
            write(*,*) '             from',flow_vp_rigid_omega,' to ', &
                              (twopi/(time_end-time_start))/km_s_to_rs_hr
          end if
          flow_vp_rigid_omega = (twopi/(time_end-time_start))/km_s_to_rs_hr
        end if
!
        do k=1,n3
          do j=1,n2
            if (j.gt.1 .and. j.lt.n2) then
                st_i_val = one/SIN(s2(j))
            else
                st_i_val = zero
            end if
            do i=1,n1
              f_local(i,j,k) = -st_i_val*EXP(-(s2(j) - pi_four     )**2/validation_run_width - &
                                              (s1(i) - pi_two      )**2/validation_run_width)  &
                               +st_i_val*EXP(-(s2(j) - pi_four     )**2/validation_run_width - &
                                              (s1(i) - threepi_two )**2/validation_run_width)
!
              fout(i,j,k)    = -st_i_val*EXP(-(s2(j) - vt_end_loc  )**2/validation_run_width - &
                                              (s1(i) - pi_two      )**2/validation_run_width)  &
                               +st_i_val*EXP(-(s2(j) - vt_end_loc  )**2/validation_run_width - &
                                              (s1(i) - threepi_two )**2/validation_run_width)
            enddo
          enddo
        enddo
!
        if (iamp0) then
          call write_2d_file((trim(output_map_root_filename)//'_final_analytic.h5') &
                              ,n1,n2,fout(:,:,1),s1,s2,ierr)
        end if
!
        deallocate (fout)
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
! ****** Set the poles.
!
      allocate (fn1(nr))
      allocate (fs1(nr))
      fn1(:) = 0.
      fs1(:) = 0.
!
      do i=1,nr
        do k=2,npm-1
          fn1(i) = fn1(i) + f(1,k,i)*dp(k)
          fs1(i) = fs1(i) + f(ntm,k,i)*dp(k)
        enddo
      enddo
!
! ****** Set the pole values to have only an m=0 component.
!
      do i=1,nr
        f(1,:,i)   = fn1(i)*twopi_i
        f(ntm,:,i) = fs1(i)*twopi_i
      enddo
!
!$omp target enter data map(to:f)
!
! ****** Balance the flux if desired (should not be used for validation runs!).
!
      if (initial_map_flux_balance) then
        do i=1,nr
          call balance_flux (f(:,:,i))
        enddo
      end if
!
! ****** Clean up memory.
!
      deallocate (s1)
      deallocate (s2)
      deallocate (s3)
      deallocate (fn1,fs1)
      deallocate (f_local)
      deallocate (f_tmp2d)
!
! ****** Write out initial condition.
!
      call write_map (trim(output_map_root_filename)//'_initial.h5')
!
end subroutine
!#######################################################################
subroutine set_realization_parameters
!
!-----------------------------------------------------------------------
!
! ****** set_realization_parameters
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use input_parameters
      use data_assimilation
      use flows
      use mesh
      use globals
      use output
      use diffusion
      use constants
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i, ierr
      character(*), parameter :: FMTI='(100(A29))'
      character(*), parameter :: FMTR='(I29,a,100(1pe16.8e3,a))'
!
!-----------------------------------------------------------------------
!
! ****** Get local number of realizations and set index arrays.
!
      call set_local_number_of_realizations
!
! ****** Set realization index arrays.
!
      allocate (local_realization_indices(nr))
!
      do i=1,nr
        local_realization_indices(i) = nr_disp_arr(iproc) + i
      enddo
!
! ****** Allocate history arrays.
!
      allocate (h_minbr(nr))
      allocate (h_maxbr(nr))
      allocate (h_minabsbr(nr))
      allocate (h_fluxp(nr))
      allocate (h_fluxm(nr))
      allocate (h_valerr(nr))
      allocate (h_fluxp_pn(nr))
      allocate (h_fluxm_pn(nr))
      allocate (h_fluxp_ps(nr))
      allocate (h_fluxm_ps(nr))
      allocate (h_area_pn(nr))
      allocate (h_area_ps(nr))
      allocate (h_eq_dipole(nr))
      allocate (h_ax_dipole(nr))
!
! ****** Set up diffusion arrays.
!
      if (advance_diffusion) then
        allocate (diffusion_coef_constant_rvec(nr))
        do i=1,nr
          diffusion_coef_constant_rvec(i) =                           &
            diffusion_coef_constants(local_realization_indices(i))
        enddo
!$omp target enter data map(to:diffusion_coef_constant_rvec)
      end if
!
! ****** Set up flow attenuation arrays.
!
      if (flow_attenuate) then
        allocate (flow_attenuate_value_i_rvec(nr))
        do i=1,nr
          if (flow_attenuate_values(local_realization_indices(i)).gt.zero) then
            flow_attenuate_value_i_rvec(i) =                            &
              one/flow_attenuate_values(local_realization_indices(i))
          else
            flow_attenuate_value_i_rvec(i) = HUGE(one)
          end if
        enddo
!$omp target enter data map(to:flow_attenuate_value_i_rvec)
      end if
!
! ****** Set up AFT differential rotation flow arrays.
!
      if (flow_dr_model.eq.1) then
        allocate (flow_dr_coef_p0_rvec(nr))
        allocate (flow_dr_coef_p2_rvec(nr))
        allocate (flow_dr_coef_p4_rvec(nr))
        do i=1,nr
          flow_dr_coef_p0_rvec(i) = flow_dr_coef_p0_values(local_realization_indices(i))
          flow_dr_coef_p2_rvec(i) = flow_dr_coef_p2_values(local_realization_indices(i))
          flow_dr_coef_p4_rvec(i) = flow_dr_coef_p4_values(local_realization_indices(i))
        enddo
!$omp target enter data map(to:flow_dr_coef_p0_rvec,&
!$omp flow_dr_coef_p2_rvec,flow_dr_coef_p4_rvec)
      end if
!
! ****** Set up AFT meridianal flow arrays.
!
      if (flow_mf_model.eq.1) then
        allocate (flow_mf_coef_p1_rvec(nr))
        allocate (flow_mf_coef_p3_rvec(nr))
        allocate (flow_mf_coef_p5_rvec(nr))
        do i=1,nr
          flow_mf_coef_p1_rvec(i) = flow_mf_coef_p1_values(local_realization_indices(i))
          flow_mf_coef_p3_rvec(i) = flow_mf_coef_p3_values(local_realization_indices(i))
          flow_mf_coef_p5_rvec(i) = flow_mf_coef_p5_values(local_realization_indices(i))
        enddo
!$omp target enter data map(to:flow_mf_coef_p1_rvec,&
!$omp flow_mf_coef_p3_rvec,flow_mf_coef_p5_rvec)
      end if
!
! ****** Set up data assimilation arrays.
!
      if (assimilate_data) then
!
        allocate (assimilate_data_lat_limit_rvec(nr))
        allocate (assimilate_data_lat_limit_tidx0_rvec(nr))
        allocate (assimilate_data_lat_limit_tidx1_rvec(nr))
!
!$omp target enter data map(alloc:assimilate_data_lat_limit_tidx0_rvec,&
!$omp assimilate_data_lat_limit_tidx1_rvec)
!
        do i=1,nr
          assimilate_data_lat_limit_rvec(i) =                         &
            assimilate_data_lat_limits(local_realization_indices(i))
        enddo
!$omp target enter data map(to:assimilate_data_lat_limit_rvec)
!
        do concurrent (i=1:nr)
          assimilate_data_lat_limit_tidx0_rvec(i) = -1
          assimilate_data_lat_limit_tidx1_rvec(i) = -1
        enddo
!
        if (assimilate_data_custom_from_mu) then
!
          allocate (assimilate_data_mu_power_rvec(nr))
          allocate (assimilate_data_mu_limit_rvec(nr))
!
          do i=1,nr
            assimilate_data_mu_power_rvec(i) =                        &
              assimilate_data_mu_powers(local_realization_indices(i))
            assimilate_data_mu_limit_rvec(i) =                        &
              assimilate_data_mu_limits(local_realization_indices(i))
          enddo
!$omp target enter data map(to:assimilate_data_mu_power_rvec,&
!$omp assimilate_data_mu_limit_rvec)
!
        end if
!
      end if
!
      if (iamp0 .and. n_realizations .gt. 1) then
!
! ****** Write realization meta data file.
!
        call ffopen (IO_TMP,'hipft_realization_parameters.out','rw',ierr)

        write (IO_TMP,FMTI) 'realization_number           ', &
                            'assimilate_data_lat_limits   ', &
                            'assimilate_data_mu_powers    ', &
                            'assimilate_data_mu_limits    ', &
                            'flow_dr_coef_p0_values       ', &
                            'flow_dr_coef_p2_values       ', &
                            'flow_dr_coef_p4_values       ', &
                            'flow_mf_coef_p1_values       ', &
                            'flow_mf_coef_p3_values       ', &
                            'flow_mf_coef_p5_values       ', &
                            'flow_attenuate_values        ', &
                            'diffusion_coef_constants     '
!
        do i=1,n_realizations
!
          write (IO_TMP,FMTR) i,' ',                                &
                                assimilate_data_lat_limits(i),' ',  &
                                assimilate_data_mu_powers(i),' ',   &
                                assimilate_data_mu_limits(i),' ',   &
                                flow_dr_coef_p0_values(i),' ',      &
                                flow_dr_coef_p2_values(i),' ',      &
                                flow_dr_coef_p4_values(i),' ',      &
                                flow_mf_coef_p1_values(i),' ',      &
                                flow_mf_coef_p3_values(i),' ',      &
                                flow_mf_coef_p5_values(i),' ',      &
                                flow_attenuate_values(i),' ',       &
                                diffusion_coef_constants(i),' '
!
        end do
!
        close (IO_TMP)
!
      end if
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
!$omp target enter data map(alloc:cos_tvec,pold,poldold,pl)
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
!$omp target exit data map(delete:cos_tvec,pold,poldold,pl)
!
end subroutine
!#######################################################################
subroutine load_flows
!
!-----------------------------------------------------------------------
!
! ****** Load flow advance.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use mesh
      use constants
      use flows
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Allocate flow arrays.
!
      allocate (vt(nt,npm,nr))
      allocate (vp(ntm,np,nr))
!
      vt(:,:,:) = 0.
      vp(:,:,:) = 0.
!$omp target enter data map(to:vt,vp)
!
      if (use_flow_from_files) call load_flow_from_files
!
      if (advection_num_method_space.eq.WENO3) call load_weno
!
      bc_flow_npole_fac = sin(dt(  1)*half)/(twopi*(one-cos(dt(  1)*half)))
      bc_flow_spole_fac = sin(dt(ntm)*half)/(twopi*(one-cos(dt(ntm)*half)))
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
      call output_histories
!
      call output_map
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
      use mpidefs
      use mesh
      use diffusion
      use flows
      use sts
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr, i
      integer*8 :: niters
      character(*), parameter :: FMT='(i10,a,5(1pe23.15e3,a),i15,a,i15)'
      character(*), parameter :: FMT2='(i10,a,15(1pe23.15e3,a))'
!
!-----------------------------------------------------------------------
!
      if (output_history_time_cadence .gt. 0.0) then
        if (time .ge. time_start+idx_hist*output_history_time_cadence) then
          output_current_history = .true.
        end if
      else
        output_current_history = .true.
      end if
!
      if (output_current_history) then
!
        idx_hist = idx_hist + 1

        do i=1,nr
!
          write(io_hist_num_filename,'(A19,I6.6,A4)') &
                "hipft_history_num_r", local_realization_indices(i), ".out"
!
          write(io_hist_sol_filename,'(A19,I6.6,A4)') &
                "hipft_history_sol_r", local_realization_indices(i), ".out"
!
          call ffopen (IO_HIST_NUM,io_hist_num_filename,'a',ierr)
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
                                 dtime_diffusion_euler, ' ',niters, ' ',  &
                                 diffusion_subcycles
!
          close(IO_HIST_NUM)
!
          call ffopen (IO_HIST_SOL,io_hist_sol_filename,'a',ierr)
!
          write(IO_HIST_SOL,FMT2) ntime,' ',time,' ',                     &
                                  h_minbr(i),' ',h_maxbr(i),' ',h_minabsbr(i),' ', &
                                  h_fluxp(i),    ' ', h_fluxm(i),    ' ',       &
                                  h_fluxp_pn(i), ' ', h_fluxm_pn(i), ' ',       &
                                  h_fluxp_ps(i), ' ', h_fluxm_ps(i), ' ',       &
                                  h_area_pn(i),  ' ', h_area_ps(i),  ' ',       &
                                  h_eq_dipole(i),' ', h_ax_dipole(i),' ',       &
                                  h_valerr(i)
!
          close(IO_HIST_SOL)
!
        enddo
!
        output_current_history = .false.
!
      end if
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
      use mpidefs
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
        if (time .ge. time_start+ &
            (idx_out-(output_map_idx_start-1))*output_map_time_cadence) then
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
        if (iamp0) then
          call ffopen (IO_MAP_OUT_LIST, &
                       io_output_map_list_filename,'a',ierr)
!
          write (IO_MAP_OUT_LIST,FMT) idx_out,'   ',time,'    ', &
                                      trim(base_filename)
!
          close(IO_MAP_OUT_LIST)
        end if
!
        if (output_flows) then
!
          base_filename = trim(output_flows_root_filename)// &
                               '_vX_idx'//idx_str//'.h5'
!
! ****** Write out the flows.
!
          call write_flows (idx_str)
!
! ****** Record output in output text file log.
!
          if (iamp0) then
            call ffopen (IO_MAP_OUT_LIST, &
                         io_output_flows_list_filename,'a',ierr)
!
            write (IO_MAP_OUT_LIST,FMT) idx_out,'   ',time,'    ', &
                                        trim(base_filename)
!
            close(IO_MAP_OUT_LIST)
          end if
!
        end if
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
      use input_parameters
      use number_types
      use output
      use mpidefs
      use globals, ONLY : f
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
      integer :: i,j,k,lsbuf,irank
      integer, dimension(:), allocatable :: displ,lbuf
      real(r_typ), dimension(:), allocatable :: gfout,tfout,rscale
      real(r_typ), dimension(:,:,:), allocatable :: gout,ftmp
!
!-----------------------------------------------------------------------
!
      allocate (ftmp(ntm,npm,nr))
!$omp target enter data map(alloc:ftmp)
!
      do concurrent (k=1:nr,j=1:npm,i=1:ntm)
        ftmp(i,j,k) = f(i,j,k)
      enddo
!
      if (output_map_flux_balance) then
        do i=1,nr
          call balance_flux (ftmp(:,:,i))
        enddo
      end if
!$omp target update from(ftmp)
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
        fout(:,:,i) = TRANSPOSE(ftmp(:,1:npm-1,i))
      enddo
!
!$omp target exit data map(delete:ftmp)
      deallocate (ftmp)
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
      end if
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
subroutine write_flows (idx_str)
!
!-----------------------------------------------------------------------
!
! ****** Write out current flow fields to disk.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use number_types
      use output
      use mpidefs
      use flows, ONLY : vt,vp
      use mesh
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: idx_str
      character(512) :: full_filename
      integer :: ierr = 0
      integer :: i
!
!-----------------------------------------------------------------------
!
!$omp target update from(vt,vp)
!
      full_filename = trim(output_flows_directory)//'/'// &
                      trim(output_flows_root_filename)// &
                      '_vt_idx'//idx_str//'.h5'
!
! ****** Write out file in pt format.
! ****** For Vt, need single point overlap in phi.
!
      allocate (fout(npm-1,nt,nr))
!
!   WARNING: -stdpar may eventually GPU-ize TRANSPOSE.
!   For now, it probably in-lines so OK, and/or requires
!   the tensor module to be explicitly used.
!   For now, this should be computed on the CPU no matter what...
!
      do i=1,nr
        fout(:,:,i) = TRANSPOSE(vt(:,1:npm1,i))/m_s_to_rs_hr
      enddo
!
      if (n_realizations.eq.1) then
        call write_2d_file (full_filename,npm1,nt,fout(:,:,1),pout,th,ierr)
      else
        write(*,*) 'SORRY!  Multiple realization flow output not yet implemented!'
!        wtime_tmp_mpi = MPI_Wtime()
!        lsbuf=(npm-1)*ntm*nr
!        allocate (tfout(lsbuf))
!        tfout(:)=reshape(fout(:,:,:),(/lsbuf/))
!        if (iamp0) then
!          allocate (gfout((npm-1)*ntm*n_realizations))
!          allocate (displ(0:nproc-1))
!          allocate (lbuf(0:nproc-1))
!          displ(0)=0
!          do irank=0,nproc-1
!            lbuf(irank)=(npm-1)*ntm*nr_arr(irank)
!          enddo
!          do irank=1,nproc-1
!            displ(irank)=displ(irank-1)+lbuf(irank-1)
!          enddo
!        end if
!
!        call MPI_Gatherv (tfout,lsbuf,ntype_real,gfout,lbuf,       &
!                  displ,ntype_real,0,MPI_COMM_WORLD,ierr)
!        deallocate(tfout)
!        wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
!        if (iamp0) then
!          allocate (rscale(n_realizations))
!          do i=1,n_realizations
!            rscale(i) = real(i,r_typ)
!          enddo
!          allocate (gout((npm-1),ntm,n_realizations))
!          gout=reshape(gfout(:),(/(npm-1),ntm,n_realizations/))
!          call write_3d_file(fname,npm-1,ntm,n_realizations,gout,pout,t,   &
!                             rscale,ierr)
!          deallocate (gout)
!          deallocate (gfout)
!          deallocate (displ)
!          deallocate (lbuf)
!          deallocate (rscale)
!        end if
!
      end if
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_FLOWS:'
        write (*,*) '### Could not write the output flow!'
        write (*,*) 'IERR: ',ierr
        write (*,*) 'File name: ', full_filename
        call endrun(.true.)
      end if
!
      deallocate (fout)

      full_filename = trim(output_flows_directory)//'/'// &
                      trim(output_flows_root_filename)// &
                      '_vp_idx'//idx_str//'.h5'
!
! ****** Write out file in pt format.
! ****** For Vp, write full 2-point overlap in phi.
!
      allocate (fout(np,ntm,nr))
!
!   WARNING: -stdpar may eventually GPU-ize TRANSPOSE.
!   For now, it probably in-lines so OK, and/or requires
!   the tensor module to be explicitly used.
!   For now, this should be computed on the CPU no matter what...
!
      do i=1,nr
        fout(:,:,i) = TRANSPOSE(vp(:,:,i))/m_s_to_rs_hr
      enddo
!
      if (n_realizations.eq.1) then
        call write_2d_file (full_filename,np,ntm,fout(:,:,1),p,t,ierr)
      else
        write(*,*) 'SORRY!  Multiple realization flow output not yet implmented!'
!        wtime_tmp_mpi = MPI_Wtime()
!        lsbuf=(npm-1)*ntm*nr
!        allocate (tfout(lsbuf))
!        tfout(:)=reshape(fout(:,:,:),(/lsbuf/))
!        if (iamp0) then
!          allocate (gfout((npm-1)*ntm*n_realizations))
!          allocate (displ(0:nproc-1))
!          allocate (lbuf(0:nproc-1))
!          displ(0)=0
!          do irank=0,nproc-1
!            lbuf(irank)=(npm-1)*ntm*nr_arr(irank)
!          enddo
!          do irank=1,nproc-1
!            displ(irank)=displ(irank-1)+lbuf(irank-1)
!          enddo
!        end if
!
!        call MPI_Gatherv (tfout,lsbuf,ntype_real,gfout,lbuf,       &
!                  displ,ntype_real,0,MPI_COMM_WORLD,ierr)
!        deallocate(tfout)
!        wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
!        if (iamp0) then
!          allocate (rscale(n_realizations))
!          do i=1,n_realizations
!            rscale(i) = real(i,r_typ)
!          enddo
!          allocate (gout((npm-1),ntm,n_realizations))
!          gout=reshape(gfout(:),(/(npm-1),ntm,n_realizations/))
!          call write_3d_file(fname,npm-1,ntm,n_realizations,gout,pout,t,   &
!                             rscale,ierr)
!          deallocate (gout)
!          deallocate (gfout)
!          deallocate (displ)
!          deallocate (lbuf)
!          deallocate (rscale)
!        end if
!
      end if
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_FLOWS:'
        write (*,*) '### Could not write the output flow!'
        write (*,*) 'IERR: ',ierr
        write (*,*) 'File name: ', full_filename
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
      use diffusion
      use output
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: t1,p1,p2
      real(r_typ) :: sumfval2, fv, fv2, eqd1, eqd2, vtt, vpt
      real(r_typ) :: tav, da_t, da_p, sn_t, d_t, cs_t, cs_p, sn_p
      real(r_typ), dimension(:,:,:), allocatable :: fval
      integer :: i, j, k
      real(r_typ) :: h_minbr_tmp, h_maxbr_tmp, h_minabsbr_tmp, h_fluxp_tmp
      real(r_typ) :: h_fluxm_tmp, h_valerr_tmp, h_fluxp_pn_tmp, h_fluxm_pn_tmp
      real(r_typ) :: h_fluxp_ps_tmp, h_fluxm_ps_tmp, h_area_pn_tmp, h_area_ps_tmp
      real(r_typ) :: h_ax_dipole_tmp
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
! ****** Only need to perform this analysis if we are writing out
! ****** the history on this step.
!
      if (output_history_time_cadence .gt. 0.0) then
        if (time .ge. time_start+idx_hist*output_history_time_cadence) then
          output_current_history = .true.
        end if
      else
        output_current_history = .true.
      end if
!
      if (output_current_history) then
!
! ****** If compute current exact solution for validation tests.
!
      allocate (fval(npm-1,ntm,nr))
!$omp target enter data map(alloc:fval)
!
! ****** [RC] NOTE!  The phi-shift validation error over time is not
!                    correct due to periodic domain. TBD later.
!
      if (validation_run .eq. 1) then
        do concurrent (k=1:nr,j=1:ntm,i=1:npm-1)
          fval(i,j,k) = fval_u0(i,j,k)*                                  &
                        exp(-42.0_r_typ*diffusion_coef_constant_rvec(k)* &
                            diffusion_coef_factor*(time-time_start))
        enddo
      elseif (validation_run .eq. 2) then
        do concurrent (k=1:nr,j=1:ntm,i=1:npm-1)
!
          vtt =       flow_vt_const*km_s_to_rs_hr*(time-time_start)
          vpt = flow_vp_rigid_omega*km_s_to_rs_hr*(time-time_start)
!
          p1 = MODULO(pi_two+vpt,      twopi)
          p2 = MODULO(threepi_two+vpt, twopi)
!
          fval(i,j,k) = -st_i(j)*EXP(-(t(j) - pi_four - vtt)**2/validation_run_width - &
                                                   (p(i)-p1)**2/validation_run_width)  &
                        +st_i(j)*EXP(-(t(j) - pi_four - vtt)**2/validation_run_width - &
                                                   (p(i)-p2)**2/validation_run_width)
        enddo
      else
        do concurrent (k=1:nr,j=1:ntm,i=1:npm-1)
          fval(i,j,k) = f(j,i,k)
        enddo
      end if
!
! ***** Get max and min metrics.
!
      do k=1,nr
        h_minbr_tmp    = large_value
        h_maxbr_tmp    = small_value
        h_minabsbr_tmp = large_value
        h_valerr_tmp   = zero
        sumfval2       = zero
!
        do concurrent (j=1:npm-2,i=1:ntm) reduce(+:sumfval2,h_valerr_tmp) &
                                          reduce(max:h_maxbr_tmp)         &
                                          reduce(min:h_minbr_tmp,h_minabsbr_tmp)
          h_minbr_tmp = min(f(i,j,k),h_minbr_tmp)
          h_maxbr_tmp = max(f(i,j,k),h_maxbr_tmp)
          h_minabsbr_tmp = min(abs(f(i,j,k)),h_minabsbr_tmp)
          fv = (f(i,j,k) - fval(j,i,k))**2
          fv2 = fval(j,i,k)**2
          h_valerr_tmp = h_valerr_tmp + fv
          sumfval2 = sumfval2 + fv2
        enddo
!
        h_minbr(k)    = h_minbr_tmp
        h_maxbr(k)    = h_maxbr_tmp
        h_minabsbr(k) = h_minabsbr_tmp
        h_valerr(k)   = h_valerr_tmp
        h_valerr(k) = sqrt(h_valerr_tmp/sumfval2)
      enddo
!
!$omp target exit data map(delete:fval)
      deallocate (fval)
!
! ****** Get integrated metrics.
!
      do k=1,nr
        h_fluxp_tmp = 0.
        h_fluxm_tmp = 0.
        h_fluxp_pn_tmp = 0.
        h_fluxm_pn_tmp = 0.
        h_fluxp_ps_tmp = 0.
        h_fluxm_ps_tmp = 0.
        h_area_pn_tmp = 0.
        h_area_ps_tmp = 0.
        eqd1 = 0.
        eqd2 = 0.
        h_ax_dipole_tmp = 0.
!
        do concurrent (i=1:npm-1,j=1:ntm) reduce(+:h_fluxp_tmp,h_fluxm_tmp,h_fluxp_pn_tmp)&
                                          reduce(+:h_fluxm_pn_tmp,h_area_pn_tmp)          &
                                          reduce(+:h_fluxp_ps_tmp,h_fluxm_ps_tmp)         &
                                          reduce(+:h_area_ps_tmp,eqd1,eqd2,h_ax_dipole_tmp)
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
            h_fluxp_tmp = h_fluxp_tmp + fv
          else
            h_fluxm_tmp = h_fluxm_tmp + fv
          end if
!
! ****** Polar fluxes and areas.
!
          if (t(j).le.pole_flux_lat_limit*d2r) then
            if (fv.gt.0.) then
              h_fluxp_pn_tmp = h_fluxp_pn_tmp + fv
            else
              h_fluxm_pn_tmp = h_fluxm_pn_tmp + fv
            end if
            h_area_pn_tmp = h_area_pn_tmp + da_t*da_p
          end if
!
          if (t(j).ge.pi-pole_flux_lat_limit*d2r) then
            if (fv.gt.0.) then
              h_fluxp_ps_tmp = h_fluxp_ps_tmp + fv
            else
              h_fluxm_ps_tmp = h_fluxm_ps_tmp + fv
            end if
            h_area_ps_tmp = h_area_ps_tmp + da_t*da_p
          end if
!
! ****** Dipoles.
!
          eqd1 = eqd1 + f(j,i,k)*sn_t*cs_p*da_t*da_p
          eqd2 = eqd2 + f(j,i,k)*sn_t*sn_p*da_t*da_p
!
          h_ax_dipole_tmp = h_ax_dipole_tmp + f(j,i,k)*cs_t*da_t*da_p
!
        enddo
!
! ****** Set fluxes to be in units of Mx and
! ****** polar areas to be in units of cm.
!
        h_fluxp(k)    = rsun_cm2*h_fluxp_tmp
        h_fluxm(k)    = rsun_cm2*h_fluxm_tmp
        h_fluxp_pn(k) = rsun_cm2*h_fluxp_pn_tmp
        h_fluxm_pn(k) = rsun_cm2*h_fluxm_pn_tmp
        h_fluxp_ps(k) = rsun_cm2*h_fluxp_ps_tmp
        h_fluxm_ps(k) = rsun_cm2*h_fluxm_ps_tmp
        h_area_pn(k)  = rsun_cm2*h_area_pn_tmp
        h_area_ps(k)  = rsun_cm2*h_area_ps_tmp
!
! ****** Set axial dipole strength.
!
        h_ax_dipole(k) = three_quarter*pi_i*h_ax_dipole_tmp
!
! ****** Set equatorial dipole strength.
!
        eqd1 = three_quarter*pi_i*eqd1
        eqd2 = three_quarter*pi_i*eqd2
        h_eq_dipole(k) = SQRT(eqd1**2 + eqd2**2)
!
      enddo
!
! ****** Reset this even though the same logic will set it in
! ****** the output_histories() subroutine.
!
      output_current_history = .false.
!
      end if
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
      use globals
      use constants
      use flows
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: j,k,i
      real(r_typ) :: br_abs_max
      real(r_typ) :: vrigid_pole
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
! ***** Add in analytic meridianal flow model.
!
        if (flow_mf_model.eq.1) call add_flow_meridianal_analytic
!
! ***** Add in analytic differential rotation model.
!
        if (flow_dr_model.eq.1) call add_flow_differential_rotation_analytic
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
! ***** Add in rotated rigid rotation velocity (might not be right).
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
! ****** Use MAX(|Br|) of two staggered Br cells around velocity vectors.
! ****** We do not update the polar boundary ghost cell for vt.
!
          do concurrent (i=1:nr,k=2:npm-1,j=2:nt-1)
            br_abs_max = MAX(ABS(f(j-1,k,i)),ABS(f(j,k,i)))
            vt(j,k,i) = vt(j,k,i)*(one - tanh(br_abs_max*flow_attenuate_value_i_rvec(i)))
          enddo
!
          do concurrent (i=1:nr,k=2:np-1,j=1:ntm)
            br_abs_max = MAX(ABS(f(j,k-1,i)),ABS(f(j,k,i)))
            vp(j,k,i) = vp(j,k,i)*(one - tanh(br_abs_max*flow_attenuate_value_i_rvec(i)))
          enddo
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
        if (advection_num_method_space .eq. WENO3) call set_lf_alpha
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
      use flows
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
      real(r_typ) :: delta_time_flow
      integer :: ln1,ln2,i,j,k
      character(1024) :: flowfile_t = ' '
      character(1024) :: flowfile_p = ' '
!
!-----------------------------------------------------------------------
!
     if (time .ge. time_of_next_input_flow) then
!
! ****** Set flow idx to flows about to be read.

       flow_from_file_current_idx = flow_from_file_current_idx + 1
!
! ****** Check if we have reached the end of the file list.
! ****** If so, reset index to the beginning.
!
       if (flow_from_file_current_idx .gt. flow_from_files_num_flows) then
         flow_from_file_current_idx = 1

         if (verbose.gt.0) then
           write(*,*) ' '
           write(*,*) '   ADD_FLOW_FROM_FILE: Resetting flow index to 1.'
         end if
       end if
!
! ****** Read the flow data.
!
! ****** VT (NT,NPM) ******
!
       if (iamp0) then
         flowfile_t = TRIM(flow_root_dir)//"/"&
                    //TRIM(flow_from_files_rel_path_t(flow_from_file_current_idx))
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
!$omp target enter data map(to:new_flow_t)
!
! ****** copy, transpose, and convert units.
!
       do concurrent (k=1:nr,j=1:npm1,i=1:nt)
         flow_from_files_current_vt(i,j,k) = m_s_to_rs_hr*new_flow_t(j,i)
       enddo
       ! Set phi ghost cell.
       do concurrent (k=1:nr,i=1:nt)
         flow_from_files_current_vt(i,npm,k) = flow_from_files_current_vt(i,2,k)
       enddo
!$omp target exit data map(delete:new_flow_t)
       deallocate(new_flow_t)
!
! ****** VP (NTM,NP)
!
       if (iamp0) then
         flowfile_p = TRIM(flow_root_dir)//"/"&
                    //TRIM(flow_from_files_rel_path_p(flow_from_file_current_idx))
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
!$omp target enter data map(to:new_flow_p)
!
! ****** copy, transpose, and convert units.
!.
       do concurrent (k=1:nr,j=1:np,i=1:ntm)
         flow_from_files_current_vp(i,j,k) = m_s_to_rs_hr*new_flow_p(j,i)
       enddo
!$omp target exit data map(delete:new_flow_p)
       deallocate(new_flow_p)
!
! ****** Update next map time.
!
!      Save current input flow time.
       time_of_current_input_flow = time_of_next_input_flow
!
       if (flow_from_file_current_idx+1 .gt. flow_from_files_num_flows) then
         delta_time_flow = flow_from_files_dt_of_last_flows
       else
         delta_time_flow = flow_from_files_times(flow_from_file_current_idx+1) - &
                           flow_from_files_times(flow_from_file_current_idx)
       end if
!
       time_of_next_input_flow = time_of_current_input_flow + delta_time_flow
!
       if (verbose.gt.0) then
         write(*,*) ' '
         write(*,*) '   ADD_FLOW_FROM_FILE: Current flow index',flow_from_file_current_idx
        write (*,*) '   ADD_FLOW_FROM_FILE: Time of current flow:  ', &
                     time_of_current_input_flow
        write (*,*) '   ADD_FLOW_FROM_FILE: Time of next flow:  ', &
                     time_of_next_input_flow
       end if
!
     end if
!
! ****** Add the last read flow.
!
     do concurrent (k=1:nr,j=1:npm,i=1:nt)
       vt(i,j,k) = vt(i,j,k) + flow_from_files_current_vt(i,j,k)
     enddo
!
     do concurrent (k=1:nr,j=1:np,i=1:ntm)
       vp(i,j,k) = vp(i,j,k) + flow_from_files_current_vp(i,j,k)
     enddo
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
      use flows
      use mpidefs
      use timing
      use read_2d_file_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: io = 0
      integer :: i,j,k,ln1,ln2
      integer :: ierr = 0
      real(r_typ) :: flow_from_files_jd_time, time_start_flows
      real(r_typ) :: flow_from_files_time_span, delta_time_flow
      character(*), parameter :: &
                    FMT = '(F11.5,1X,A11,1X,A11)'
      character(1024) :: flowfile_t = ' '
      character(1024) :: flowfile_p = ' '
      real(r_typ), dimension(:,:), allocatable :: flow_t,flow_p
      real(r_typ), dimension(:), allocatable :: s1,s2
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
      flow_from_files_num_flows = -1
      do
        READ (IO_DATA_IN,'(a)',iostat=io)
        if (io/=0) exit
        flow_from_files_num_flows = flow_from_files_num_flows + 1
      enddo
      REWIND (IO_DATA_IN)
!
! ****** Allocate space to store the list.
!
      allocate(flow_from_files_times      (flow_from_files_num_flows))
      allocate(flow_from_files_rel_path_t (flow_from_files_num_flows))
      allocate(flow_from_files_rel_path_p (flow_from_files_num_flows))
!
! ****** Read the list (skip first row).
!
      READ (IO_DATA_IN,*)
      do i=1,flow_from_files_num_flows
        READ (IO_DATA_IN,FMT) flow_from_files_jd_time,       &
                              flow_from_files_rel_path_t(i), &
                              flow_from_files_rel_path_p(i)
!
! ****** Make list of times.
!
        flow_from_files_times(i) = twentyfour*(flow_from_files_jd_time)
      enddo
      CLOSE (IO_DATA_IN)
!
! ****** Make list of times start at zero (usually the case anyways).
!
      do i=1,flow_from_files_num_flows
        flow_from_files_times(i) = flow_from_files_times(i) &
                                   - flow_from_files_times(1)
      enddo
!
! ****** Get time span of full list of flows.
!
      flow_from_files_time_span = flow_from_files_times(flow_from_files_num_flows) &
                                  - flow_from_files_times(1)
!
! ****** Get start time within flows.
!
      time_start_flows = MODULO(time_start,flow_from_files_time_span)
!
! ****** Get currently active flow file index.
!
      flow_from_file_current_idx = MINLOC(ABS(flow_from_files_times(:) - &
                                              time_start_flows),1)
      if (flow_from_files_times(flow_from_file_current_idx) &
          .gt. time_start_flows) then
        flow_from_file_current_idx = flow_from_file_current_idx-1
      end if
!
! ****** Allocate and read current flows.
!
      allocate (flow_from_files_current_vt(nt,npm,nr))
      allocate (flow_from_files_current_vp(ntm,np,nr))
!
! ****** VT (NT,NPM) ******
!
       if (iamp0) then
         flowfile_t = TRIM(flow_root_dir)//"/" &
                    //TRIM(flow_from_files_rel_path_t(flow_from_file_current_idx))
         call read_2d_file (flowfile_t,ln1,ln2,flow_t,s1,s2,ierr)
         deallocate(s1)
         deallocate(s2)
       endif
!
       wtime_tmp_mpi = MPI_Wtime()
       if (.not.iamp0) then
        allocate (flow_t(npm1,nt))
       endif
       call MPI_Bcast (flow_t,flow_t_size,ntype_real,0,MPI_COMM_WORLD,ierr)
       wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
       do j=1,npm1
         do i=1,nt
           flow_from_files_current_vt(i,j,:) = m_s_to_rs_hr*flow_t(j,i)
         enddo
       enddo
       ! Set phi ghost cell.
       flow_from_files_current_vt(:,npm,:) = flow_from_files_current_vt(:,2,:)
       deallocate(flow_t)
!$omp target enter data map(to:flow_from_files_current_vt)
!
!
! ****** VP (NTM,NP)
!
       if (iamp0) then
         flowfile_p = TRIM(flow_root_dir)//"/"&
                    //TRIM(flow_from_files_rel_path_p(flow_from_file_current_idx))
         call read_2d_file (flowfile_p,ln1,ln2,flow_p,s1,s2,ierr)
         deallocate(s1)
         deallocate(s2)
       endif
!
       wtime_tmp_mpi = MPI_Wtime()
       if (.not.iamp0) then
        allocate (flow_p(np,ntm1))
       endif
       call MPI_Bcast (flow_p,flow_p_size,ntype_real,0,MPI_COMM_WORLD,ierr)
       wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
       do j=1,np
         do i=1,ntm
           flow_from_files_current_vp(i,j,:) = m_s_to_rs_hr*flow_p(j,i)
         enddo
       enddo
       deallocate(flow_p)
!$omp target enter data map(to:flow_from_files_current_vp)
!
! ****** Assume delta time from last conflow to first is the same as the
! ****** delta time of the last two files.
!
      flow_from_files_dt_of_last_flows = flow_from_files_times(flow_from_files_num_flows) - &
                                         flow_from_files_times(flow_from_files_num_flows-1)
!
! ****** Get simulation time of currently loaded flow (not the same as time_start!)
!
      time_of_current_input_flow =  time_start - time_start_flows   &
                                  + flow_from_files_times(flow_from_file_current_idx)
!
      if (flow_from_file_current_idx+1 .gt. flow_from_files_num_flows) then
        delta_time_flow = flow_from_files_dt_of_last_flows
      else
        delta_time_flow = flow_from_files_times(flow_from_file_current_idx+1) - &
                          flow_from_files_times(flow_from_file_current_idx)
      end if
!
      time_of_next_input_flow = time_of_current_input_flow + delta_time_flow
!
! ******* Print some diagnostics.
!
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '   FLOWS_FROM_FILE: Number of flows in flow list: ', &
                     flow_from_files_num_flows
        write (*,*) '   FLOWS_FROM_FILE: File name of first theta flow file:     ', &
                     trim(flow_from_files_rel_path_t(1))
        write (*,*) '   FLOWS_FROM_FILE: File name of first phi flow file:     ', &
                     trim(flow_from_files_rel_path_p(1))
        write (*,*) '   FLOWS_FROM_FILE: Time span of all flows:  ', &
                     flow_from_files_time_span
        write (*,*) '   FLOWS_FROM_FILE: Start time within flows:  ', &
                     time_start_flows
        write (*,*) '   FLOWS_FROM_FILE: DT of last two flows:  ', &
                     flow_from_files_dt_of_last_flows
        write (*,*) '   FLOWS_FROM_FILE: Time of current flow:  ', &
                     time_of_current_input_flow
        write (*,*) '   FLOWS_FROM_FILE: Time of next flow:  ', &
                     time_of_next_input_flow
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
      use sources
      use mesh
      use globals, ONLY : time
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ) :: alpha,tmp1,tmp2
!
!-----------------------------------------------------------------------
!
! ****** Reset source array.
!
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        source(j,k,i) = 0.
      enddo
!
! ****** Add in source from file.
!
      if (source_from_file) then
        do concurrent (i=1:nr,k=1:npm,j=1:ntm)
          source(j,k,i) = source(j,k,i) + source_from_file_data(j,k,i)
        enddo
      end if
!
! ****** Random flux emergence.
!
      if (source_rfe_model.gt.0) then
!
        if (time .ge. time_of_next_source_rfe) then
!
          time_of_prev_source_rfe = time_of_next_source_rfe
          time_of_next_source_rfe = time_of_prev_source_rfe + source_rfe_lifetime
!
          tmp2 = rfe2
          tmp1 = rfe1
          rfe1 = tmp2
          rfe2 = tmp1

          call generate_rfe (source_rfe(:,:,:,rfe2))
!$omp target update to(source_rfe(:,:,:,rfe2))
!
! ****** Flux balance the source term.
!
          do i=1,nr
            call balance_flux (source_rfe(:,:,i,rfe2))
          end do
!
        end if
!
        alpha = (time - time_of_prev_source_rfe)/source_rfe_lifetime
!
        do concurrent (i=1:nr,k=1:npm,j=1:ntm)
          source(j,k,i) = source(j,k,i) +                       &
                          (one-alpha)*source_rfe(j,k,i,rfe1) +  &
                                alpha*source_rfe(j,k,i,rfe2)
        enddo
      end if
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
      use constants
      use mesh, ONLY : t,ntm,nr
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: io = 0
      integer :: i,j
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
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '   LOAD_DATA_ASSIMILATION: Number of maps in map list: ', &
                     num_maps_in_list
        write (*,*) '   LOAD_DATA_ASSIMILATION: Start date:     ', &
                     trim(map_times_actual_ut_str(1))
        write (*,*) '   LOAD_DATA_ASSIMILATION: End   date:     ', &
                     trim(map_times_actual_ut_str(num_maps_in_list))
        write (*,*) '   LOAD_DATA_ASSIMILATION: File name of first map:     ', &
                     trim(map_files_rel_path(1))
        write (*,*) '   LOAD_DATA_ASSIMILATION: File name of first used map:     ', &
                     trim(map_files_rel_path(i))
      end if
!
! ****** Find indices of latitude limit in main mesh tvec.
!
      do concurrent(i=1:nr)
        do j=1,ntm
          if (t(j).ge.assimilate_data_lat_limit_rvec(i)*d2r .and.           &
                          assimilate_data_lat_limit_tidx0_rvec(i).eq.-1) then
            assimilate_data_lat_limit_tidx0_rvec(i) = j
          end if
        enddo
        do j=ntm,1,-1
          if (t(j).le.(pi-assimilate_data_lat_limit_rvec(i)*d2r) .and.      &
                          assimilate_data_lat_limit_tidx1_rvec(i).eq.-1) then
            assimilate_data_lat_limit_tidx1_rvec(i) = j
          end if
        enddo
      enddo
!
end subroutine
!#######################################################################
subroutine assimilate_new_data (new_data)
!
!-----------------------------------------------------------------------
!
! ****** Assimilate data.
! ****** This assumes the data is the 1st slice and the
! ****** weight is in the second slice.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use globals
      use constants, ONLY : one,half
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
      real(r_typ), dimension(:,:,:), allocatable :: deltaf
!
!-----------------------------------------------------------------------
!
      allocate(deltaf(ntm,npm,nr))
!$omp target enter data map(alloc:deltaf)

      do concurrent (i=1:nr,k=1:npm1,j=1:ntm)
        deltaf(j,k,i) = new_data(k,j,2,i)*(new_data(k,j,1,i) - f(j,k,i))
      enddo
      call set_periodic_bc_3d (deltaf,ntm,npm,nr)
!
! ****** Balance the added flux if requested.
!
      if (assimilate_data_balance_flux) then
        if (verbose.gt.0) write(*,*) '   ASSIMILATE_NEW_DATA: About to balance DBr'
        call balance_flux (deltaf)
      end if
!
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f(j,k,i) = f(j,k,i) + deltaf(j,k,i)
      enddo
      call set_periodic_bc_3d (f,ntm,npm,nr)
!
!$omp target exit data map(delete:deltaf)
      deallocate(deltaf)
!
end subroutine
!#######################################################################
subroutine apply_data_assimilation
!
!-----------------------------------------------------------------------
!
! ****** Apply data assimilation.
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
      use sources
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
      integer :: npm_nd,ntm_nd,nslices,i,j,k,l
      character(1024) :: mapfile = ' '
!
!-----------------------------------------------------------------------
!
      wtime_tmp = MPI_Wtime()
!
! ****** If it is time to load new data, read it from file.
!
      if (assimilate_data.and.time.ge.time_of_next_input_map) then
!
        if (iamp0) then
          if (verbose.gt.0) then
            write(*,*)
            write(*,*) '   UPDATE_FIELD: Loading field index ',current_map_input_idx
            write(*,*) '   UPDATE_FIELD: Time of next input map: ',time_of_next_input_map
          end if
!
          mapfile = TRIM(assimilate_data_map_root_dir)//"/"&
                    //TRIM(map_files_rel_path(current_map_input_idx))
!
          call read_3d_file (mapfile,npm_nd,ntm_nd,nslices,new_data2d,s1,s2,s3,ierr)
! ******  TODO:  Add error checking here
!               (make sure scales match run until interp is added)
          deallocate(s1)
          deallocate(s2)
          deallocate(s3)
        end if
!
        wtime_tmp_mpi = MPI_Wtime()
        call MPI_Bcast (npm_nd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (ntm_nd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (nslices,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
        if (.not.iamp0) allocate (new_data2d(npm_nd,ntm_nd,nslices))
!
        call MPI_Bcast (new_data2d,npm_nd*ntm_nd*nslices,ntype_real,0,MPI_COMM_WORLD,ierr)
        wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
! ****** Allocate and set rank-local array across realizations
! ****** of assimilation data based on read-in 2D values.
!
!$omp target enter data map(to:new_data2d)
        allocate(new_data(npm_nd,ntm_nd,nslices,nr))
!$omp target enter data map(alloc:new_data)
        do concurrent(i=1:nr,j=1:nslices,k=1:ntm_nd,l=1:npm_nd)
          new_data(l,k,j,i) = new_data2d(l,k,j)
        enddo
!$omp target exit data map(delete:new_data2d)
        deallocate(new_data2d)
!
! ****** Modify the data based on user choices.
! ****** Note:  The data is stored in slices defined as:
! ******  [1] The data
! ******  [2] The default uncertainty/weights (default is mu^4 with mu limit of 0.1?)
! ******  [3] The "mu" value of the data (including negative values on back of Sun)
!
! ****** Make a custom uncertaintity map using mu.
!
        if (assimilate_data_custom_from_mu) then

          do concurrent(i=1:nr,k=1:ntm_nd,l=1:npm_nd)
            if (new_data(l,k,3,i).gt.assimilate_data_mu_limit_rvec(i)) then
              new_data(l,k,2,i) = new_data(l,k,3,i)**assimilate_data_mu_power_rvec(i)
            else
              new_data(l,k,2,i) = zero
            end if
          enddo

        end if
!
! ****** Apply a latitude limiter on data (solid cutoff for now).
!
        do concurrent(i=1:nr)
          if (assimilate_data_lat_limit_rvec(i).gt.0) then
            do concurrent(k=1:assimilate_data_lat_limit_tidx0_rvec(i),l=1:npm_nd)
              new_data(l,k,2,i) = zero
            enddo
            do concurrent(k=assimilate_data_lat_limit_tidx1_rvec(i):ntm_nd,l=1:npm_nd)
              new_data(l,k,2,i) = zero
            enddo
          end if
        enddo
!
! ****** Assimilate the data.
!
        call assimilate_new_data (new_data)
!
! ****** Clean up.
!
!$omp target exit data map(delete:new_data)
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

        if (iamp0.and.verbose.gt.0) then
          write(*,*) '   UPDATE_FIELD: Time of next input map: ', &
                      time_of_next_input_map
        end if
!
      end if
!
      wtime_data_assimilation = wtime_data_assimilation + MPI_Wtime() - wtime_tmp
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
      use sources
      use data_assimilation
      use flows
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
        if (verbose.gt.1) then
          write(*,*)
          write(*,*) '   UPDATE_TIMESTEP: dtime_old = ',dtime_old
        end if
!
! ****** Default timestep is one giant step for remaining time.
!
        dtime_global = time_end - time
        if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_endrun = ',dtime_global
        dtime_reason = 'endrun'
        time_step_changed = .true.
!
! ****** Flow CFL time step limit.
!
        if (advance_flow) then
          dtmax_flow = huge(one)
          if (timestep_flow_needs_updating) then
            call get_flow_dtmax (dtmax_flow)
            timestep_flow_needs_updating = .false.
            if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_flow = ',dtmax_flow
          end if
          if (dtmax_flow .lt. dtime_global) then
            dtime_global = dtmax_flow
            dtime_reason = 'flowcfl'
          end if
        end if
!
! ****** Maximum allowed timestep.
!
        if (dt_max .lt. dtime_global) then
          dtime_global = dt_max
          dtime_reason = 'dtmax'
          if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_dtmax = ',dtime_global
        end if
!
        if (dtime_global .lt. dt_min) then
          write (*,*) '   UPDATE_TIMESTEP: Warning! Time step is smaller than DTMIN!'
        end if
!
        timestep_needs_updating = .false.
!
      end if
!
! ****** Check for next output time, cut dt to match exactly.
!
      if (output_map_time_cadence .gt. 0.0) then
        if (time + dtime_global .ge. time_start+ &
            (idx_out-(output_map_idx_start-1))*output_map_time_cadence) then
          dtime_global = (time_start+ &
                         (idx_out-(output_map_idx_start-1))*output_map_time_cadence) &
                         - time
          dtime_reason = 'output'
          timestep_needs_updating = .true.
          timestep_flow_needs_updating = .true.
          time_step_changed = .true.
          if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_output = ',dtime_global
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
          timestep_flow_needs_updating = .true.
          time_step_changed = .true.
          if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_assim = ',dtime_global
        end if
      end if
!
! ****** Check for next flow input time, cut dt to match exactly.
!
      if (advance_flow .and. use_flow_from_files) then
        if (time + dtime_global .ge. time_of_next_input_flow) then
          dtime_global = time_of_next_input_flow - time
          dtime_reason = 'flowfile'
          timestep_needs_updating = .true.
          flow_needs_updating = .true.
          timestep_flow_needs_updating = .true.
          time_step_changed = .true.
          if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_nxtflow = ',dtime_global
        end if
      end if
!
! ****** Check for next source rfe switch.
!
      if (advance_source .and. source_rfe_model.gt.0) then
        if (time + dtime_global .ge. time_of_next_source_rfe) then
          dtime_global = MIN(time_of_next_source_rfe - time,   &
                             half*source_rfe_lifetime)
          dtime_reason = 'rfeswap'
          timestep_needs_updating = .true.
          timestep_flow_needs_updating = .true.
          time_step_changed = .true.
          if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_nxtrfe = ',dtime_global
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
        if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_end = ',dtime_global
      end if
!
! ****** Make sure STS is recomputed if timestep has changed.
!
      if (advance_diffusion .and. time_step_changed) then
        need_to_load_sts = .true.
      end if
!
      if (verbose.gt.1) write(*,*) '   UPDATE_TIMESTEP: dtime_global = ',dtime_global
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
      if (strang_splitting) then
        dtime_local_half = dtime_local*half
        if (advance_source)    call source_step    (dtime_local_half)
        if (advance_flow)      call advection_step (dtime_local_half)
        if (advance_diffusion) call diffusion_step (dtime_local)
        if (advance_flow)      call advection_step (dtime_local_half)
        if (advance_source)    call source_step    (dtime_local_half)
      else
        if (advance_flow)      call advection_step (dtime_local)
        if (advance_diffusion) call diffusion_step (dtime_local)
        if (advance_source)    call source_step    (dtime_local)
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
      use flows
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
      if (advection_num_method_time .eq. 1) then
        call advection_step_fe (dtime_local)
      elseif (advection_num_method_time .eq. 2) then
        call advection_step_ssprk33 (dtime_local)
      elseif (advection_num_method_time .eq. 3) then
        call advection_step_ssprk43 (dtime_local)
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
      use constants
      use globals
      use mesh
      use diffusion
      use sts, ONLY : need_to_load_sts
      use matrix_storage
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
      integer, parameter :: diffusion_subcycles_max = 60
      integer :: i
      logical :: we_are_done
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
! ****** Start diffusion subcycles.
!
      i = 1
      time_stepped = 0.
      we_are_done = .false.
!
! ****** Set the initial subcycle time step.
!
      if (auto_sc) then
!
! This will allow the first dtflux calculation to be unbound.
! After it is set, it is used as a lower bound.
!
        dtime_local2 = zero
      else
        dtime_local2 = dtime_local/REAL(diffusion_subcycles,r_typ)
      end if
!
! ****** Start the subcycle loop.
!
      do
!
! ****** Compute and set the flux time-step if requested.
!
        if (auto_sc) then
          call get_dtime_diffusion_ptl (dtime_local2)
          need_to_load_sts = .true.
        end if
!
! ****** Lower time step to reach the end time exactly if
! ****** it would otherwise go too far.
!
        if (auto_sc .and. time_stepped + dtime_local2 .ge. dtime_local) then
          dtime_local2 = dtime_local - time_stepped
          need_to_load_sts = .true.
          we_are_done = .true.
        end if
!
! ****** Take one last big step if hitting max # cycles.
! ****** Note: this is not the same effect of having that number from the start.
!
        if (auto_sc .and. i .eq. diffusion_subcycles_max-1) then
          dtime_local2 = dtime_local - time_stepped
          need_to_load_sts = .true.
          we_are_done = .true.
        end if
!
! ****** If using the flux time-step, compute the needed subcycles
! ****** if we were to equal-step from this point to the end.
! ****** This is a DIAGNOSTIC only!  This is not used!
!
        if (auto_sc .and. verbose.gt.0) then
          diffusion_subcycles = CEILING((dtime_local - time_stepped)/dtime_local2)
        end if
!
! ****** Take the diffusion step.
!
        if (diffusion_num_method.eq.1) then
!
          call diffusion_step_fe_cd (dtime_local2)
!
        else if (diffusion_num_method.eq.2.or.   &
                 diffusion_num_method.eq.3) then
!
          call diffusion_step_sts_cd (dtime_local2)
!
        end if
!
! ****** Update the amount of the large time step that has been stepped.
!
        time_stepped = time_stepped + dtime_local2
!
        if (verbose.gt.1) then
          timetmp = time
          time = time + time_stepped
          call analysis_step
          call output_histories
          write(*,*) ' '
          write(*,*) '   DIFFUSION_STEP: --> Cycle #',i,' SC: ',diffusion_subcycles, &
                     ' dtcycle:',dtime_local2,' time_stepped: ',time_stepped
          flush(OUTPUT_UNIT)
          time = timetmp
        end if
!
! ****** If we are using a set number of subcycles, check if done.
!
        if (.not. auto_sc .and. i .eq. diffusion_subcycles) then
          we_are_done = .true.
        end if
!
! ****** Exit the loop if the full time step is complete.
!
        if (we_are_done) exit
!
        i = i + 1
!
      enddo
!
! ****** If using auto-cycle, set this for history output.
!
      if (auto_sc) then
        diffusion_subcycles = i
        if (verbose.gt.0) write(*,*) '   DIFFUSION_STEP: Total number of cycles = ',diffusion_subcycles
      end if
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
! ****** Advance the field with the source term by the time-step DT.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use timing
      use globals, ONLY : f
      use sources
      use mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(IN) :: dtime_local
      real(r_typ) :: t1
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
      t1 = MPI_Wtime()
!
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f(j,k,i) = f(j,k,i) + dtime_local*source(j,k,i)
      enddo
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
      real(r_typ), INTENT(INOUT), dimension(n1,n2,n3) :: a
      integer, INTENT(IN) :: n1,n2,n3
      integer :: i,j
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
      use diffusion
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
!$omp target enter data map(to:coef)
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
        coef(j,k,3,i)=-(coef(j,k,1,i) + coef(j,k,2,i) &
                       +coef(j,k,4,i) + coef(j,k,5,i))
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
! *** Estimate maximum eigenvalue of all realizations
! *** using Gershgorin disks:
!
      max_eig = 0.
!
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1) reduce(max:max_eig)
        gersh_rad = zero
        do d=1,5
          gersh_rad = gersh_rad+ABS(coef(j,k,d,i))
        enddo
        max_eig = MAX(gersh_rad,max_eig)
      enddo
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
subroutine get_dtime_diffusion_ptl (dtime_ptl)
!
!-----------------------------------------------------------------------
!
! ****** Get the practical time step limit for diffusion.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use mesh
      use constants
      use timing
      use input_parameters, ONLY : verbose
      use globals, ONLY : f
      use diffusion
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), INTENT(INOUT) :: dtime_ptl
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ierr
      real(r_typ) :: axabsmax, dtime_in, deltau, deltaf
      real(r_typ), dimension(:,:,:), allocatable :: Af
      real(r_typ), parameter :: safe = 0.95_r_typ
      real(r_typ), parameter :: fmin = 1e-14_r_typ
!
!-----------------------------------------------------------------------
!
      dtime_in = dtime_ptl
      dtime_ptl = HUGE(one)
!
      if (verbose.gt.1) then
        write (*,*) ' '
        write (*,*) '   GET_DTIME_DIFFUSION_PTL: dtime_in          = ',dtime_in
      end if
!
      allocate (Af(ntm,npm,nr))
!$omp target enter data map(alloc:Af)
!
      call diffusion_operator_cd (f,Af)
!
! ****** Find maximum value of |F| (|Af|).
!
      axabsmax = -one
!
      do concurrent (i=1:nr,k=1:npm,j=1:ntm) reduce(max:axabsmax)
        if (ABS(f(j,k,i)) .gt. fmin) then
          axabsmax = MAX(ABS(Af(j,k,i)),axabsmax)
        end if
      enddo
!
! ****** Get maximum over all MPI ranks.
!
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Allreduce (MPI_IN_PLACE,axabsmax,1,ntype_real,     &
                          MPI_MAX,MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
      if (verbose.gt.1) write (*,*) '   GET_DTIME_DIFFUSION_PTL: MAX(|F|)          = ',axabsmax
!
      if (axabsmax .gt. zero) then
!
        do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1) reduce(min:dtime_ptl)
!
! ****** If we are at the location of max(|F|), calculate the PTL.
! ****** If there are multiple locations of |F|=max(|F|), take min(PTL).
!
          if (axabsmax .eq. ABS(Af(j,k,i)) .and. &
              ABS(f(j,k,i)) .gt. fmin) then
!
            deltau =  f(j,k,i) -  f(j-1,k,i)
            deltaf = Af(j,k,i) - Af(j-1,k,i)
            if (deltau*deltaf.lt.-fmin) dtime_ptl = MIN(-deltau/deltaf,dtime_ptl)
!
            deltau =  f(j,k,i) -  f(j+1,k,i)
            deltaf = Af(j,k,i) - Af(j+1,k,i)
            if (deltau*deltaf.lt.-fmin) dtime_ptl = MIN(-deltau/deltaf,dtime_ptl)
!
            deltau =  f(j,k,i) -  f(j,k-1,i)
            deltaf = Af(j,k,i) - Af(j,k-1,i)
            if (deltau*deltaf.lt.-fmin) dtime_ptl = MIN(-deltau/deltaf,dtime_ptl)
!
            deltau =  f(j,k,i) -  f(j,k+1,i)
            deltaf = Af(j,k,i) - Af(j,k+1,i)
            if (deltau*deltaf.lt.-fmin) dtime_ptl = MIN(-deltau/deltaf,dtime_ptl)
          end if
!
        enddo
!
! ****** Get global minimum across all MPI ranks.
!
        wtime_tmp_mpi = MPI_Wtime()
        call MPI_Allreduce (MPI_IN_PLACE,dtime_ptl,1,ntype_real,  &
                            MPI_MIN,MPI_COMM_WORLD,ierr)
        wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() - wtime_tmp_mpi
!
! ****** Add safety factor.
!
        dtime_ptl = safe*dtime_ptl
!
        if (verbose.gt.1) write (*,*) '   GET_DTIME_DIFFUSION_PTL: PTL timestep      = ',dtime_ptl
!
! ****** Don't let the new dt be smaller than the previous one.
!
        if (dtime_ptl .lt. dtime_in) dtime_ptl = dtime_in
!
        if (verbose.gt.1) write (*,*) '   GET_DTIME_DIFFUSION_PTL: PTL timestep used = ',dtime_ptl
      end if
!
!$omp target exit data map(delete:Af)
      deallocate (Af)
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
      use globals
      use diffusion
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
!$omp target enter data map(alloc:u0,dty0,ykm1,ukm1,ukm2)
!
      call diffusion_operator_cd (f,dty0)
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
        call diffusion_operator_cd (ukm1,ykm1)
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
!$omp target exit data map(delete:u0,dty0,ykm1,ukm1,ukm2)
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
      use diffusion
      use flows
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
!$omp target exit data map(delete:sts_uj,sts_vj,sts_ubj,sts_gj)
        deallocate (sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
      end if
!
      if (need_to_load_sts) then
        allocate (sts_uj(sts_s))
        allocate (sts_vj(sts_s))
        allocate (sts_ubj(sts_s))
        allocate (sts_gj(sts_s))
        allocate (sts_b(sts_s))
!$omp target enter data map(alloc:sts_uj,sts_vj,sts_ubj,sts_gj)
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
!$omp target update to(sts_uj,sts_vj,sts_ubj,sts_gj)
!
      if (verbose.gt.1) then
        write (*,*) ' '
        write (*,*) '   LOAD_STS:  *** Diffusion: 2nd-order RKL2  ***'
        write (*,*) '   LOAD_STS:  Super-time-step used              = ',dtime_local
        write (*,*) '   LOAD_STS:  Euler time-step                   = ',dtime_diffusion_euler
        write (*,*) '   LOAD_STS:  Number of STS iterations needed   = ',sts_s
        write (*,*) '   LOAD_STS:  Number of Euler iterations needed = ',CEILING(dtime_local/dtime_diffusion_euler)
        write (*,*) '   LOAD_STS:  Potential max speedup of STS      = ', &
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
      use diffusion
      use flows
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
!$omp target exit data map(delete:sts_uj,sts_vj,sts_ubj,sts_gj)
        deallocate (sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
      end if
!
      if (need_to_load_sts) then
        allocate (sts_uj(sts_s))
        allocate (sts_vj(sts_s))
        allocate (sts_ubj(sts_s))
        allocate (sts_gj(sts_s))
        allocate (sts_b(sts_s))
!$omp target enter data map(alloc:sts_uj,sts_vj,sts_ubj,sts_gj)
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
!$omp target update to(sts_uj,sts_vj,sts_ubj,sts_gj)
!
      if (verbose.gt.1) then
        write (*,*) ' '
        write (*,*) '   LOAD_STS:  *** Diffusion: 2nd-order RKG2  ***'
        write (*,*) '   LOAD_STS:  Super-time-step used              = ',dtime_local
        write (*,*) '   LOAD_STS:  Euler time-step                   = ',dtime_diffusion_euler
        write (*,*) '   LOAD_STS:  Number of STS iterations needed   = ',sts_s
        write (*,*) '   LOAD_STS:  Number of Euler iterations needed = ',CEILING(dtime_local/dtime_diffusion_euler)
        write (*,*) '   LOAD_STS:  Potential max speedup of STS      = ', &
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
      use input_parameters
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
      s%hdf32 = output_single_precision
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
      use input_parameters
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
      s%hdf32 = output_single_precision
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
!
      th(1) = th(2) - dth(2)
      th(nt) = th(ntm1) + dth(ntm1)
!
      dth(1) = dth(2)
      dth(nt) = dth(ntm1)
!
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
      sth(:) = sin(ABS(th(:)))
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
      allocate (pout(npm1))
      pout(:) = p(1:npm1)
!
!$omp target enter data map(to:t,p,dt,dt_i,st,st_i,ct,dp,dp_i,th,ph,&
!$omp dth,dth_i,sth,sth_i,cth,dph,dph_i)
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
      use input_parameters
      use constants
      use read_2d_file_interface
      use globals
      use sts
      use diffusion
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
! ****** Allocate memory for the total diffusion coef on the half-half mesh.
!
      allocate (diffusion_coef(nt,np,nr))
!
! ****** Set the initial diffusion coef to the uniform value.
!
      do i=1,nr
        do k=1,np
          do j=1,nt
            diffusion_coef(j,k,i) = diffusion_coef_constant_rvec(i)
          enddo
        enddo
      enddo
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
        call MPI_Bcast (nfp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (nft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (.not. iamp0) then
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
        if (verbose.gt.0) then
          write (*,*)
          write (*,*) '   LOAD_DIFFUSION: Diffusion coef from file: ', &
                                           trim(diffusion_coef_filename)
          write (*,*) '   LOAD_DIFFUSION: Minimum value = ',minval(diffusion_coef_file)
          write (*,*) '   LOAD_DIFFUSION: Maximum value = ',maxval(diffusion_coef_file)
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

        if (verbose.gt.0) then
          write (*,*)
          write (*,*) '   LOAD_DIFFUSION: Grid-based diffusion is activated.'
        end if

      end if
!
      diffusion_coef(:,:,:) = diffusion_coef_factor*diffusion_coef(:,:,:)
!$omp target enter data map(to:diffusion_coef)
!
      call load_diffusion_matrix
!
! ****** Get stable explicit Euler timestep.
!
      call get_dtime_diffusion_euler (dtime_diffusion_euler)
!
! ****** Set flag to auto-compute number of subcycles
! ****** using flux time-step.
!
      if (diffusion_subcycles.eq.0) then
        auto_sc = .true.
      else
        auto_sc = .false.
      end if
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
      do j=2,nt-2
        dt_total = dt(j-1) + dt(j) + dt(j+1)
        D_C_CPt(j) =           dt(j  )/(dt(j  )+dt(j+1))
        D_C_MCt(j) =           dt(j  )/(dt(j-1)+dt(j  ))
        D_M_Tt (j) =           dt(j-1)/dt_total
        D_CP_Tt(j) = (dt(j  )+dt(j+1))/dt_total
        D_P_Tt (j) =           dt(j+1)/dt_total
        D_MC_Tt(j) = (dt(j-1)+dt(j  ))/dt_total
      enddo
!
      do k=2,np-1
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
!$omp target enter data map(to:d_c_cpt,d_c_mct,d_m_tt,d_cp_tt,d_p_tt,&
!$omp d_mc_tt,d_c_cpp,d_c_mcp,d_m_tp,d_cp_tp,d_p_tp,d_mc_tp,alpha_t,&
!$omp alpha_p)
!
end subroutine
!#######################################################################
subroutine set_lf_alpha
!
!-----------------------------------------------------------------------
!
! ****** Set Local Lax-Friedrichs alpha for use with WENO3.
! ****** This should be called any time the velocity field changes.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use flows, ONLY : vt,vp
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
! ****** Get alpha for theta:  (ntm=ntm1=nt-1)
!
! ****** Get inner points.
!
      do concurrent (i=1:nr,k=1:npm,j=3:nt-3)
        alpha_t(j,k,i) = MAX(ABS(vt(j-2,k,i)), &
                             ABS(vt(j-1,k,i)), &
                             ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)), &
                             ABS(vt(j+2,k,i)), &
                             ABS(vt(j+3,k,i)))
      enddo
!
! ****** Edge cases.
!
      do concurrent (i=1:nr,k=1:npm)
        j = 1
        alpha_t(j,k,i) = MAX(ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)), &
                             ABS(vt(j+2,k,i)), &
                             ABS(vt(j+3,k,i)))
        j = 2
        alpha_t(j,k,i) = MAX(ABS(vt(j-1,k,i)), &
                             ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)), &
                             ABS(vt(j+2,k,i)), &
                             ABS(vt(j+3,k,i)))
        j = nt-2
        alpha_t(j,k,i) = MAX(ABS(vt(j-2,k,i)), &
                             ABS(vt(j-1,k,i)), &
                             ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)), &
                             ABS(vt(j+2,k,i)))
        j = nt-1
        alpha_t(j,k,i) = MAX(ABS(vt(j-2,k,i)), &
                             ABS(vt(j-1,k,i)), &
                             ABS(vt(j  ,k,i)), &
                             ABS(vt(j+1,k,i)))
      enddo
!
! ****** Get alpha for phi: (npm=np)
!
! ****** Get inner points.
!
      do concurrent (i=1:nr,k=3:np-3,j=2:ntm-1)
        alpha_p(j,k,i) = MAX(ABS(vp(j,k-2,i)), &
                             ABS(vp(j,k-1,i)), &
                             ABS(vp(j,k  ,i)), &
                             ABS(vp(j,k+1,i)), &
                             ABS(vp(j,k+2,i)), &
                             ABS(vp(j,k+3,i)))
      enddo
!
! ****** Edge cases (periodic).
!
      do concurrent (i=1:nr,j=2:ntm-1)
        k = 2
        alpha_p(j,k,i) = MAX(ABS(vp(j,np-2,i)), &
                             ABS(vp(j,k-1 ,i)), &
                             ABS(vp(j,k   ,i)), &
                             ABS(vp(j,k+1 ,i)), &
                             ABS(vp(j,k+2 ,i)), &
                             ABS(vp(j,k+3 ,i)))
        k = np-2
        alpha_p(j,k,i) = MAX(ABS(vp(j,k-2,i)), &
                             ABS(vp(j,k-1,i)), &
                             ABS(vp(j,k  ,i)), &
                             ABS(vp(j,k+1,i)), &
                             ABS(vp(j,k+2,i)), &
                             ABS(vp(j,3  ,i)))
        k = np-1
        alpha_p(j,k,i) = MAX(ABS(vp(j,k-2,i)), &
                             ABS(vp(j,k-1,i)), &
                             ABS(vp(j,k  ,i)), &
                             ABS(vp(j,k+1,i)), &
                             ABS(vp(j,3  ,i)), &
                             ABS(vp(j,4  ,i)))
      enddo
!
! ****** By seaming the periodic dimension,
! ****** the values for k=1 and k=npm (np) are set.
!
      call set_periodic_bc_3d (alpha_p,ntm,npm,nr)
!
end subroutine
!#######################################################################
subroutine load_source
!
!-----------------------------------------------------------------------
!
! ****** Load source term.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use mpidefs
      use sources
      use output
      use input_parameters
      use constants
      use read_2d_file_interface
      use globals, ONLY : time,local_realization_indices
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ierr
!
      integer :: nft,nfp,num_cells,seed_n
      integer,     dimension(:),   allocatable :: seed_new
      real(r_typ), dimension(:),   allocatable :: tf,pf,rscale
      real(r_typ), dimension(:,:), allocatable :: sf
      real(r_typ), dimension(:,:), allocatable :: f_tmp2d
      real(r_typ) :: get_gaussian
!
!-----------------------------------------------------------------------
!
! ****** Initialize source array.
!
      allocate (source(ntm,npm,nr))
      source(:,:,:) = 0.
!$omp target enter data map(to:source)
!
! ****** Add the source file if it was specified.
!
      if (source_from_file) then
!
! ****** Allocate memory for the source term.
!
        allocate (source_from_file_data(ntm,npm,nr))
        source_from_file_data(:,:,:) = 0.
!
        allocate (f_tmp2d(ntm,npm))
!
        if (iamp0) then
          call read_2d_file (source_filename,nfp,nft,sf,pf,tf,ierr)
        endif
        wtime_tmp_mpi = MPI_Wtime()
        call MPI_Bcast (nfp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast (nft,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (.not. iamp0) then
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
! ****** Interpolate (badly) the source term onto the main mesh (t,p).
!
        call interp2d (nft,nfp,tf,pf,sf,ntm,npm,t,p,f_tmp2d,ierr)
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
        do i=1,nr
          source_from_file_data(:,:,i) = f_tmp2d(:,:)
        enddo
!
        deallocate (sf)
        deallocate (tf)
        deallocate (pf)
        deallocate (f_tmp2d)
!
        if (verbose.gt.0) then
          write (*,*)
          write (*,*) '   LOAD_SOURCE: A source term was read in from file: ', &
                     trim(source_filename)
          write (*,*) '   LOAD_SOURCE: Minimum value = ',minval(source_from_file_data(:,:,:))
          write (*,*) '   LOAD_SOURCE: Maximum value = ',maxval(source_from_file_data(:,:,:))
        end if
!$omp target enter data map(to:source_from_file_data)
      end if
!
      if (source_rfe_model.gt.0) then
!
        allocate (source_rfe(ntm,npm,nr,2))
        source_rfe(:,:,:,:) = 0.
!
        if (source_rfe_seed .ne. -9999) then
          call RANDOM_SEED (size=seed_n)
          allocate (seed_new(seed_n))
          seed_new(:)=source_rfe_seed
          call RANDOM_SEED (put=seed_new)
          deallocate (seed_new)
        else
          call RANDOM_SEED ()
        end if
!
        num_cells = ntm*(npm-2)
!
! ****** Get unsigned flux per hour per cell/hour in code units
! ****** and set up realization arrays.
!
        if (source_rfe_model.eq.1) then
!
          allocate (uf_per_cell_per_hour_rvec(nr))
          do i=1,nr
            uf_per_cell_per_hour_rvec(i) = &
            source_rfe_total_unsigned_flux_per_hours(     &
                                 local_realization_indices(i))* &
                                 input_flux_fac/(rsun_cm2*num_cells)
          enddo
!$omp target enter data map(to:uf_per_cell_per_hour_rvec)
!
        elseif (source_rfe_model.eq.2) then
!
          allocate (uf_per_area_per_hour_rvec(nr))
          do i=1,nr
            uf_per_area_per_hour_rvec(i) = &
            source_rfe_total_unsigned_flux_per_hours(     &
                                 local_realization_indices(i))* &
                                   input_flux_fac/(rsun_cm2*four*pi)
          enddo
!$omp target enter data map(to:uf_per_area_per_hour_rvec)
!
        end if
!
! ****** Generate both starting frames of the rfe.
!
        call generate_rfe (source_rfe(:,:,:,rfe1))
        call generate_rfe (source_rfe(:,:,:,rfe2))
!$omp target enter data map(to:source_rfe)
!
! ****** Flux balance the source term.
!
        do i=1,nr
          call balance_flux (source_rfe(:,:,i,rfe1))
          call balance_flux (source_rfe(:,:,i,rfe2))
        end do
!
! ****** Set the time for the next rfe generation.
!
        time_of_next_source_rfe = time + source_rfe_lifetime
        time_of_prev_source_rfe = time
!
      end if
!
end subroutine
!#######################################################################
subroutine add_flow_differential_rotation_analytic
!
!-----------------------------------------------------------------------
!
! ****** Add analytical differential rotation flow to vp.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use constants
      use globals
      use flows
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
      do concurrent (i=1:nr,k=1:np,j=1:ntm)
        vp(j,k,i) = vp(j,k,i) +                            &
                  m_s_to_rs_hr*st(j)*(                     &
                        flow_dr_coef_p0_rvec(i) +          &
                        flow_dr_coef_p2_rvec(i)*ct(j)**2 + &
                        flow_dr_coef_p4_rvec(i)*ct(j)**4)
      enddo
!
end subroutine
!#######################################################################
subroutine add_flow_meridianal_analytic
!
!-----------------------------------------------------------------------
!
! ****** Add analytical meridianal flow to vt.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use constants
      use globals
      use flows
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
      do concurrent (i=1:nr,k=1:npm,j=1:nt)
        vt(j,k,i) = vt(j,k,i) -                          &
                  m_s_to_rs_hr*sth(j)*                   &
                  (flow_mf_coef_p1_rvec(i)*cth(j) +      &
                   flow_mf_coef_p3_rvec(i)*cth(j)**3 +   &
                   flow_mf_coef_p5_rvec(i)*cth(j)**5)
      enddo
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
      use flows
      use constants
      use input_parameters, ONLY : strang_splitting
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
      do concurrent (i=1:nr,k=2:npm-1,j=2:ntm-1) reduce(min:dtmax)
        kdotv = MAX(ABS(vt(j,k,i)),ABS(vt(j+1,k,i)))*dt_i(j) &
              + MAX(ABS(vp(j,k,i)),ABS(vp(j,k+1,i)))*st_i(j)*dp_i(k)
!        deltat = (one/sqrt(two))/MAX(kdotv,small_value)
        deltat = half/MAX(kdotv,small_value)
        dtmax = MIN(dtmax,deltat)
      enddo
!
      dtmaxflow = safety*dtmax
!
! ****** Since the flow time step will be cut in half when using
! ****** Strang splitting, we can double the time-step limit here.
!
      if (strang_splitting) then
        dtmaxflow = two*dtmaxflow
      end if
!
! ****** The CFL for the SSPRK(4,3) method is 2, so we can double dt!
!
      if (advection_num_method_time .eq. SSPRK43) then
        dtmaxflow = two*dtmaxflow
      end if
!
      wtime_tmp_mpi = MPI_Wtime()
      call MPI_Allreduce (MPI_IN_PLACE,dtmaxflow,1,ntype_real, &
                          MPI_MIN,MPI_COMM_WORLD,ierr)
      wtime_mpi_overhead = wtime_mpi_overhead + MPI_Wtime() &
                           - wtime_tmp_mpi
!
end subroutine
!#######################################################################
subroutine diffusion_step_fe_cd (dtime_local)
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
      use input_parameters
      use globals
      use diffusion
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
!$omp target enter data map(alloc:fold)
!
! ****** Subcycle at a stable time-step.
!
      n_stable_diffusion_cycles = CEILING(dtime_local/dtime_diffusion_euler,8)
      dtime_euler_used = dtime_local/n_stable_diffusion_cycles
!
      do i=1,n_stable_diffusion_cycles
!
        call diffusion_operator_cd (f,fold)
!
        do concurrent (el=1:nr,k=1:npm,j=1:ntm)
          f(j,k,el) = f(j,k,el) + dtime_euler_used*fold(j,k,el)
        enddo
!
      enddo
!
!$omp target exit data map(delete:fold)
      deallocate (fold)
!
end subroutine
!#######################################################################
subroutine advection_step_fe (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with advection by the time-step dtime_local
! ****** using the Foward-Euler method.
!
!-----------------------------------------------------------------------
!
      use number_types
      use globals, ONLY : f
      use mesh
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
      real(r_typ), dimension(:,:,:), allocatable :: aop
!
!-----------------------------------------------------------------------
!
      allocate (aop(ntm,npm,nr))
!$omp target enter data map(alloc:aop)
!
      call advection_operator (f,aop)
!
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f(j,k,i) = f(j,k,i) - dtime_local*aop(j,k,i)
      enddo
!
!$omp target exit data map(delete:aop)
      deallocate (aop)
!
end subroutine
!#######################################################################
subroutine advection_step_ssprk33 (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with advection by the time-step dtime_local
! ****** using the low storage SSPRK(3,3) (RK3TVD) method.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use constants
      use mesh
      use globals, ONLY : f
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
      real(r_typ), dimension(:,:,:), allocatable :: f1,aop
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays for RK.
!
      allocate (f1(ntm,npm,nr))
      allocate (aop(ntm,npm,nr))
!$omp target enter data map(alloc:f1,aop)
!
      call advection_operator (f,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f1(j,k,i) =               f(j,k,i)     &
                  - dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator (f1,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f1(j,k,i) = three_quarter*f(j,k,i)     &
                       + quarter*f1(j,k,i)     &
          - quarter*dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator (f1,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f(j,k,i) =          third*f(j,k,i)     &
                     + two_third*f1(j,k,i)     &
        - two_third*dtime_local*aop(j,k,i)
      enddo
!
! ****** Deallocate temporary arrays.
!
!$omp target exit data map(delete:f1,aop)
      deallocate (f1)
      deallocate (aop)
!
end subroutine
!#######################################################################
subroutine advection_step_ssprk43 (dtime_local)
!
!-----------------------------------------------------------------------
!
! ****** Advance the field with advection by the time-step dtime_local
! ****** using the low-storage SSPRK(4,3) method.
!
!-----------------------------------------------------------------------
!
      use number_types
      use input_parameters
      use constants
      use mesh
      use globals, ONLY : f
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
      real(r_typ), dimension(:,:,:), allocatable :: f1,aop
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays for RK.
!
      allocate (f1(ntm,npm,nr))
      allocate (aop(ntm,npm,nr))
!$omp target enter data map(alloc:f1,aop)
!
      call advection_operator (f,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f1(j,k,i) =                    f(j,k,i)         &
                  - half*dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator (f1,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f1(j,k,i) =                   f1(j,k,i)         &
                  - half*dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator (f1,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f1(j,k,i) =          two_third*f(j,k,i)         &
                              + third*f1(j,k,i)         &
                 - sixth*dtime_local*aop(j,k,i)
      enddo
!
      call advection_operator (f1,aop)
      do concurrent (i=1:nr,k=1:npm,j=1:ntm)
        f(j,k,i) =                    f1(j,k,i)         &
                  - half*dtime_local*aop(j,k,i)
      enddo
!
! ****** Deallocate temporary arrays.
!
!$omp target exit data map(delete:f1,aop)
      deallocate (f1)
      deallocate (aop)
!
end subroutine
!#######################################################################
subroutine advection_operator (ftemp,aop)
!
!-----------------------------------------------------------------------
!
! ****** Compute DIV(v*f)
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use flows, ONLY : advection_num_method_space,UW,WENO3
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
      if (advection_num_method_space .eq. UW) then
        call advection_operator_upwind (ftemp,aop)
      elseif (advection_num_method_space .eq. WENO3) then
        call advection_operator_weno3 (ftemp,aop)
      end if
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
      use flows
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
      real(r_typ) :: fn,fs
      real(r_typ), dimension(:,:,:), allocatable :: flux_t,flux_p
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays
! ****** (we set all points here so no initialization needed)
!
      allocate (flux_t(nt,npm,nr))
      allocate (flux_p(ntm,np,nr))
!$omp target enter data map(alloc:flux_t,flux_p)
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
      do concurrent (i=1:nr)
        fn = zero
        fs = zero
        do concurrent(k=2:npm-1) reduce(+:fn,fs)
          fn = fn + flux_t(   2,k,i)*dp(k)
          fs = fs + flux_t(ntm1,k,i)*dp(k)
        enddo
! ****** Note that the south pole needs a sign change since the
! ****** theta flux direction is reversed.
        do concurrent(k=2:npm-1)
          aop(  1,k,i) =  fn*bc_flow_npole_fac
          aop(ntm,k,i) = -fs*bc_flow_spole_fac
        enddo
      enddo
!
! ****** Set periodic boundary condition.
!
      call set_periodic_bc_3d (aop,ntm,npm,nr)
!
!$omp target exit data map(delete:flux_t,flux_p)
      deallocate (flux_t)
      deallocate (flux_p)
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
      use flows
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
      real(r_typ) :: p0p,p1p,p0m,p1m
      real(r_typ) :: B0p,B1p,B0m,B1m
      real(r_typ) :: w0p,w1p,w0m,w1m
      real(r_typ) :: wp_sum,wm_sum
      real(r_typ) :: OM0p,OM1p,OM0m,OM1m
      real(r_typ) :: up,um
      real(r_typ) :: cct
      real(r_typ) :: fn,fs
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
!$omp target enter data map(alloc:flux_t,lp,ln)
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
        w0m = D_P_Tt(j-1) /(dt(j-1) + B0m)**2
        w1m = D_MC_Tt(j-1)/(dt(j-1) + B1m)**2
        w0p = D_M_Tt(j)   /(dt(j  ) + B0p)**2
        w1p = D_CP_Tt(j)  /(dt(j  ) + B1p)**2
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
        flux_t(2,k,i) = vt(2,k,i)*half*((one - cct)*ftemp(2,k,i)      &
                                      + (one + cct)*ftemp(1,k,i))
!
        cct=sign(upwind,vt(ntm1,k,i))
        flux_t(ntm1,k,i) = vt(ntm1,k,i)*half*((one - cct)*ftemp(ntm1,k,i) &
                                            + (one + cct)*ftemp(ntm2,k,i))
!
      enddo
!
!$omp target exit data map(delete:lp,ln)
      deallocate (LP)
      deallocate (LN)
!
! ### PHI ######################################
!
      allocate (flux_p(ntm,np,nr))
      allocate (LP(ntm,0:npm,nr))
      allocate (LN(ntm,0:npm,nr))
!$omp target enter data map(alloc:flux_p,lp,ln)
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
        w0m = D_P_Tp (k-1)/(dp(k-1) + B0m)**2
        w1m = D_MC_Tp(k-1)/(dp(k-1) + B1m)**2
        w0p = D_M_Tp (k)  /(dp(k  ) + B0p)**2
        w1p = D_CP_Tp(k)  /(dp(k  ) + B1p)**2
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
!$omp target exit data map(delete:lp,ln)
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
      do concurrent (i=1:nr)
        fn = zero
        fs = zero
        do concurrent(k=2:npm-1) reduce(+:fn,fs)
          fn = fn + flux_t(   2,k,i)*dp(k)
          fs = fs + flux_t(ntm1,k,i)*dp(k)
        enddo
! ****** Note that the south pole needs a sign change since the
! ****** theta flux direction is reversed.
        do concurrent(k=2:npm-1)
          aop(  1,k,i) =  fn*bc_flow_npole_fac
          aop(ntm,k,i) = -fs*bc_flow_spole_fac
        enddo
      enddo
!
! ****** Set periodic phi boundary condition.
!
      call set_periodic_bc_3d (aop,ntm,npm,nr)
!
!$omp target exit data map(delete:flux_t,flux_p)
      deallocate (flux_t)
      deallocate (flux_p)
!
end subroutine
!#######################################################################
subroutine diffusion_operator_cd (x,y)
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
      use diffusion
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
      real(r_typ) :: fn2_fn1,fs2_fs1
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
        y(j,k,i) =  coef(j,k,1,i)*x(j,  k-1,i)  &
                  + coef(j,k,2,i)*x(j-1,k  ,i)  &
                  + coef(j,k,3,i)*x(j  ,k  ,i)  &
                  + coef(j,k,4,i)*x(j+1,k  ,i)  &
                  + coef(j,k,5,i)*x(j,  k+1,i)
      enddo
!
! ****** Compute boundary points.
!
      do concurrent(i=1:nr)
        fn2_fn1 = zero
        fs2_fs1 = zero
        do concurrent(k=2:npm-1) reduce(+:fn2_fn1,fs2_fs1)
          fn2_fn1 = fn2_fn1 + (diffusion_coef(1    ,k,i)        &
                            +  diffusion_coef(2    ,k,i))       &
                             * (x(2    ,k,i) - x(1  ,k,i))*dp(k)
          fs2_fs1 = fs2_fs1 + (diffusion_coef(nt-1 ,k,i)        &
                            +  diffusion_coef(nt   ,k,i))       &
                             * (x(ntm-1,k,i) - x(ntm,k,i))*dp(k)
        enddo
        do concurrent(k=1:npm)
          y(  1,k,i) = fn2_fn1*dt_i(  1)*dt_i(  1)*pi_i
          y(ntm,k,i) = fs2_fs1*dt_i(ntm)*dt_i(ntm)*pi_i
        enddo
      enddo
!
! ****** Set the periodic boundary conditions.
!
      call set_periodic_bc_3d (y,ntm,npm,nr)
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
               res_nt,                                                 &
               res_np,                                                 &
               n_realizations,                                         &
               initial_map_filename,                                   &
               initial_map_flux_balance,                               &
               validation_run,                                         &
               validation_run_width,                                   &
               time_start,                                             &
               time_end,                                               &
               output_map_idx_start,                                   &
               output_map_directory,                                   &
               output_map_root_filename,                               &
               output_map_idx_cadence,                                 &
               output_map_time_cadence,                                &
               output_map_2d,                                          &
               output_map_flux_balance,                                &
               output_flows,                                           &
               output_flows_directory,                                 &
               output_history_time_cadence,                            &
               output_single_precision,                                &
               dt_min,                                                 &
               dt_max,                                                 &
               strang_splitting,                                       &
               pole_flux_lat_limit,                                    &
               advance_flow,                                           &
               flow_vp_rigid_omega,                                    &
               flow_rigid_omega,                                       &
               flow_vt_const,                                          &
               flow_attenuate,                                         &
               flow_attenuate_value,                                   &
               flow_attenuate_values,                                  &
               flow_dr_model,                                          &
               flow_dr_coef_p0_value,                                  &
               flow_dr_coef_p2_value,                                  &
               flow_dr_coef_p4_value,                                  &
               flow_dr_coef_p0_values,                                 &
               flow_dr_coef_p2_values,                                 &
               flow_dr_coef_p4_values,                                 &
               flow_mf_model,                                          &
               flow_mf_coef_p1_value,                                  &
               flow_mf_coef_p3_value,                                  &
               flow_mf_coef_p5_value,                                  &
               flow_mf_coef_p1_values,                                 &
               flow_mf_coef_p3_values,                                 &
               flow_mf_coef_p5_values,                                 &
               use_flow_from_files,                                    &
               flow_list_filename,                                     &
               flow_root_dir,                                          &
               flow_num_method,                                        &
               upwind,                                                 &
               advance_diffusion,                                      &
               diffusion_coef_filename,                                &
               diffusion_coef_grid,                                    &
               diffusion_coef_factor,                                  &
               diffusion_num_method,                                   &
               diffusion_subcycles,                                    &
               diffusion_coef_constant,                                &
               diffusion_coef_constants,                               &
               assimilate_data,                                        &
               assimilate_data_balance_flux,                           &
               assimilate_data_map_list_filename,                      &
               assimilate_data_map_root_dir,                           &
               assimilate_data_custom_from_mu,                         &
               assimilate_data_mu_power,                               &
               assimilate_data_mu_powers,                              &
               assimilate_data_lat_limit,                              &
               assimilate_data_lat_limits,                             &
               assimilate_data_mu_limit,                               &
               assimilate_data_mu_limits,                              &
               advance_source,                                         &
               source_from_file,                                       &
               source_filename,                                        &
               source_rfe_model,                                       &
               source_rfe_total_unsigned_flux_per_hour,                &
               source_rfe_total_unsigned_flux_per_hours,               &
               source_rfe_lifetime,                                    &
               source_rfe_seed
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(len=:), allocatable :: infile
      integer arglen
!
!-----------------------------------------------------------------------
!
      if (command_argument_count().ge.1) then
        call GET_COMMAND_ARGUMENT(1,length=arglen)
        allocate(character(arglen) :: infile)
        call GET_COMMAND_ARGUMENT(1,value=infile)
      else
        allocate(character(9):: infile)
        infile='hipft.in'
      end if
!
!-----------------------------------------------------------------------
!
! ****** Read the input file (all MPI ranks read it).
!
      call ffopen (8,infile,'r',ierr)
!
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in READ_INPUT_FILE:'
          write (*,*) '### Could not open the input file.'
          write (*,*) 'File name: ',trim(infile)
        end if
        call endrun(.true.)
      end if
!
      read (8,hipft_input_parameters)
      close (8)
!
! ****** Write the NAMELIST parameter values to file.
!
      if (iamp0) then
        call ffopen (8,'hipft_run_parameters_used.out','rw',ierr)
        write (8,hipft_input_parameters)
        close (8)
      end if
!
! ****** Check & process input parameters.
!
      call process_input_parameters
!
end subroutine
!#######################################################################
subroutine process_input_parameters
!
!-----------------------------------------------------------------------
!
! ****** Process the input parameters.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use globals
      use mpidefs
      use flows
      use output
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      if (output_map_time_cadence.gt.0            &
          .and. output_map_idx_cadence.gt.0) then
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
      if (n_realizations .gt. MAX_REALIZATIONS) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in SETUP:'
          write (*,*) '### Requested realizations greater than maximum allowed.'
          write (*,*) 'Number of requested realizations = ',n_realizations
          write (*,*) 'Maximum number of realizations allowed:  ',MAX_REALIZATIONS
        end if
        call endrun (.true.)
      end if
!
      if (initial_map_filename .eq. '') then
        if (res_np .lt. 3 .or. res_nt .lt. 3) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in SETUP:'
            write (*,*) '### Must specify input map or proper res for initial map.'
            write (*,*) 'Input map name = ',initial_map_filename
            write (*,*) 'res_nt,res_np:  ',res_nt,res_np
          end if
          call endrun (.true.)
        end if
      end if
!
      if (initial_map_filename .ne. '' .and. res_np*res_nt .gt. 0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### WARNING in SETUP:'
          write (*,*) '### You have specified a resolution and provided an '
          write (*,*) '### input map file.  Interpolation is not yet'
          write (*,*) '### supported.  Using the native resolution of the map.'
          write (*,*) 'Input map name = ',initial_map_filename
          write (*,*) 'res_nt,res_np:  ',res_nt,res_np
        end if
        res_nt = 0
        res_np = 0
      end if
!
! ****** Set realizations to defaults if they have not been set.
!
      do i=1,MAX_REALIZATIONS
        if (assimilate_data_lat_limits(i) .lt. 0.) then
          assimilate_data_lat_limits(i) = assimilate_data_lat_limit
        end if
        if (assimilate_data_mu_limits(i) .lt. 0.) then
          assimilate_data_mu_limits(i) = assimilate_data_mu_limit
        end if
        if (assimilate_data_mu_powers(i) .lt. 0.) then
          assimilate_data_mu_powers(i) = assimilate_data_mu_power
        end if
        if (flow_dr_coef_p0_values(i) .lt. -100000.) then
          flow_dr_coef_p0_values(i) = flow_dr_coef_p0_value
        end if
        if (flow_dr_coef_p2_values(i) .lt. -100000.) then
          flow_dr_coef_p2_values(i) = flow_dr_coef_p2_value
        end if
        if (flow_dr_coef_p4_values(i) .lt. -100000.) then
          flow_dr_coef_p4_values(i) = flow_dr_coef_p4_value
        end if
        if (flow_mf_coef_p1_values(i) .lt. -100000.) then
          flow_mf_coef_p1_values(i) = flow_mf_coef_p1_value
        end if
        if (flow_mf_coef_p3_values(i) .lt. -100000.) then
          flow_mf_coef_p3_values(i) = flow_mf_coef_p3_value
        end if
        if (flow_mf_coef_p5_values(i) .lt. -100000.) then
          flow_mf_coef_p5_values(i) = flow_mf_coef_p5_value
        end if
        if (flow_attenuate_values(i) .lt. 0.) then
          flow_attenuate_values(i) = flow_attenuate_value
        end if
        if (diffusion_coef_constants(i) .lt. 0.) then
          diffusion_coef_constants(i) = diffusion_coef_constant
        end if
        if (source_rfe_total_unsigned_flux_per_hours(i) .lt. 0.) then
          source_rfe_total_unsigned_flux_per_hours(i) = source_rfe_total_unsigned_flux_per_hour
        end if
      enddo
!
! ****** Set numerical method presets.
!
      if (flow_num_method .eq. 1) then
        advection_num_method_time  = FE
        advection_num_method_space = UW
      elseif (flow_num_method .eq. 2) then
        advection_num_method_time  = SSPRK33
        advection_num_method_space = UW
      elseif (flow_num_method .eq. 3) then
        advection_num_method_time  = SSPRK33
        advection_num_method_space = WENO3
      elseif (flow_num_method .eq. 4) then
        advection_num_method_time  = SSPRK43
        advection_num_method_space = WENO3
      else
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in PROCESS_INPUT_PARAMETERS:'
          write (*,*) '### Invalid flow_num_method specified.'
          write (*,*) 'flow_num_method  = ',flow_num_method
          write (*,*) 'Valid values: {1,2,3,4}'
          write (*,*) '  1: FE+Upwind'
          write (*,*) '  2: SSPRK(3,3)+Upwind'
          write (*,*) '  3: SSPRK(3,3)+WENO3'
          write (*,*) '  4: SSPRK(4,3)+WENO3'
        end if
        call endrun (.true.)
      end if
!
! ****** Set up starting output index.
!
      idx_out = output_map_idx_start - 1
!
end subroutine
!#######################################################################
function get_gaussian() result(fn_val)
!
!-----------------------------------------------------------------------
!
! ****** Generate a random normal deviate using the polar method.
! ****** Reference: Marsaglia,G. & Bray,T.A.
! ****** 'A convenient method for generating normal variables',
! ****** Siam Rev., vol.6, 260-264, 1964.
! ****** Adapted from Alam Miller's implementation from:
! ****** https://jblevins.org/mirror/amiller
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
      real(r_typ) :: fn_val
!
!-----------------------------------------------------------------------
!
      real(r_typ)       :: u,s
      real(r_typ), save :: v,sln
      logical, save     :: second = .false.
!
      if (second) then
!
! ****** If second, use the second random number generated on last call.
!
        second = .false.
        fn_val = v*sln
      else
!
! ****** First call; generate a pair of random normals.
!
        second = .true.
        do
          call RANDOM_NUMBER(u)
          call RANDOM_NUMBER(v)
          u = SCALE(u,1) - one
          v = SCALE(v,1) - one
          s = u*u + v*v
          if (s.lt.one .and. s.gt.0.) exit
        enddo
        sln = SQRT(-SCALE(LOG(s),1)/s)
        fn_val = u*sln
      end if
!
      return
end function get_gaussian
!#######################################################################
subroutine get_flux (map,flux_p,flux_m)
!
!-----------------------------------------------------------------------
!
! ****** Get the flux of the input theta-phi Br map "map".
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
      use mesh
      use input_parameters, ONLY : verbose
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntm,npm) :: map
      real(r_typ)  :: flux_p, flux_m
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: tav,sn_t,d_t,da_t,da_p,sn_p,fv
      integer :: k,j
!
!-----------------------------------------------------------------------
!
! ******* Get the positive and negative flux of the map.
!
      flux_p = zero
      flux_m = zero
!
      do concurrent (k=1:npm-1,j=1:ntm) reduce(+:flux_p,flux_m)
!
        if (j.eq.1) then
          tav = half*(t(1) + t(2))
          sn_t = sin(tav)
          d_t = dth(2)
          da_t = quarter*sn_t*d_t
        else if (j.eq.ntm) then
          tav=half*(t(ntm)+t(ntm-1))
          sn_t = sin(tav)
          d_t = dth(ntm)
          da_t = quarter*sn_t*d_t
        else
          da_t = st(j)*dt(j)
        end if
!
        if (k.eq.1) then
          da_p = half*dph(1)
        else if (k.eq.npm-1) then
          da_p = half*dph(npm-1)
        else
          da_p = dp(k)
        end if
!
        fv = map(j,k)*da_t*da_p
!
        if (fv.gt.zero) then
          flux_p = flux_p + fv
        else
          flux_m = flux_m + fv
        end if
!
      enddo
!
! ****** Report flux if in verbose mode.
!
      if (verbose.gt.0) then
        write(*,*)
        write(*,*) '   GET_FLUX: Flux(+):       ',flux_p
        write(*,*) '   GET_FLUX: Flux(-):       ',flux_m
        write(*,*) '   GET_FLUX: Flux(signed):  ',flux_p +     flux_m
        write(*,*) '   GET_FLUX: Flux(unsigned):',flux_p + ABS(flux_m)
      end if
!
end subroutine get_flux
!#######################################################################
subroutine balance_flux (map)
!
!-----------------------------------------------------------------------
!
! ****** Flux balance the input theta-phi Br map "map" using
! ****** multiplicative flux-balancing.
! ****** The input map is overwriten with the flux-balanced version.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
      use mesh
      use input_parameters, ONLY : verbose
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntm,npm) :: map
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: flux_p,flux_m,ratio,br
      integer :: k,j
!
!-----------------------------------------------------------------------
!
! ******* First, get the total positive and negative flux of the map.
!
      call get_flux (map,flux_p,flux_m)
!
! ****** Check if any flux balancing is required (only happens in trivial cases).
!
      if (flux_m + flux_p.eq.zero) return
!
! ******* Flux balance the map.
!
      ratio = SQRT(ABS(flux_p/flux_m))
!
      do concurrent (k=1:npm,j=1:ntm)
        br = map(j,k)
        if (br.gt.0) then
          map(j,k) = br/ratio
        else
          map(j,k) = br*ratio
        end if
      enddo
!
! ******* Check new flux if in verbose mode.
!
      if (verbose.gt.0) call get_flux (map,flux_p,flux_m)
!
end subroutine balance_flux
!#######################################################################
subroutine generate_rfe (rfe)
!
!-----------------------------------------------------------------------
!
! ****** Generate random flux source term.
!
!-----------------------------------------------------------------------
!
      use number_types
      use constants
      use mesh
      use sources
      use mpidefs
      use input_parameters, ONLY: verbose,source_rfe_model
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntm,npm,nr) :: rfe
      real(r_typ), dimension(ntm,npm)    :: rnd_nums
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ierr
      real(r_typ) :: get_gaussian
      real(r_typ) :: rnd_num,tav,sn_t,d_t,da_ti,fn1,fs1,normfac
!
!-----------------------------------------------------------------------
!
! ****** Set random number map on rank 0.  We use the same random numbers
! ****** for all realizations since parallel random numbers are tricky to
! ****** set right.  We could use different numbers for each realization
! ****** across 1 rank, but that can be a later update if requested.
!
      if (iamp0) then
        normfac = pi/sqrt(twopi)
        do k=2,npm-1
          do j=1,ntm
            rnd_nums(j,k) = normfac*get_gaussian()
          enddo
        enddo
      end if
!
! ****** Broadcast random number map to all ranks.
!
      call MPI_Bcast (rnd_nums,ntm*npm,ntype_real,0,MPI_COMM_WORLD,ierr)
!
      do k=2,npm-1
        do j=1,ntm
!
          if (source_rfe_model.eq.1) then
            if (j.eq.1) then
              tav=half*(t(1)+t(2))
              sn_t = sin(tav)
              d_t = dth(2)
              da_ti = one/(quarter*sn_t*d_t)
            else if (j.eq.ntm) then
              tav=half*(t(ntm)+t(ntm-1))
              sn_t = sin(tav)
              d_t = dth(ntm)
              da_ti = one/(quarter*sn_t*d_t)
            else
              da_ti = st_i(j)*dt_i(j)
            end if
          end if
!
          do i=1,nr
            if (source_rfe_model.eq.1) then
              rfe(j,k,i) = rnd_nums(j,k)*uf_per_cell_per_hour_rvec(i)*da_ti*dp_i(k)
            elseif (source_rfe_model .eq. 2) then
              rfe(j,k,i) = rnd_nums(j,k)*uf_per_area_per_hour_rvec(i)
            end if
          enddo
!
        enddo
      enddo
!
! ****** Boundary conditions.
!
! ****** Enforce periodicity.
!
      do i=1,nr
        rfe(:,  1,i) = rfe(:,npm-1,i)
        rfe(:,npm,i) = rfe(:,2,    i)
      enddo
!
! ****** Poles.
!
      do i=1,nr
        fn1 = zero
        fs1 = zero
        do k=1,npm-2
          fn1 = fn1 + rfe(  1,k,i)*dp(k)
          fs1 = fs1 + rfe(ntm,k,i)*dp(k)
        enddo
        fn1 = fn1*twopi_i
        fs1 = fs1*twopi_i
        rfe(  1,:,i) = fn1
        rfe(ntm,:,i) = fs1
      enddo
!
      if (iamp0.and.verbose.gt.0) then
        write (*,*)
        write (*,*) '   GENERATE_RFE: A RFE source term was generated: '
        write (*,*) '   GENERATE_RFE: RFE Model = ',source_rfe_model
        write (*,*) '   GENERATE_RFE: Minimum value = ',MINVAL(rfe(:,:,:))
        write (*,*) '   GENERATE_RFE: Maximum value = ',MAXVAL(rfe(:,:,:))
      end if
!
end subroutine generate_rfe
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
!   - Added the SSPRK33 time-stepping scheme for advection.
!     This is now the default.
!   - DEFAULT CHANGE: Made the default splitting method to be Strang
!                     splitting. The "-strang" input flag is now
!                     "-nostrang" and disables strang splitting.
!   - Added "-fm" input flag to set advection numerical method.
!     1) FE+UW  2) SSPRK33+UW  3) SSPRK33+WENO3 (coming soon).
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
! 02/10/2023, RC, Version 0.16.0:
!   - BUG FIX: Advection time step was not being recalculated when
!              it should have been in some cases.
!
! 02/23/2023, RC, Version 0.17.0:
!   - Added new option OUTPUT_FLOWS which will output the flows
!     at the same cadense as the maps.  Set OUTPUT_FLOWS_DIRECTORY
!     to set the directory for the flows to go.
!   - Added FLOW_ATTENUATE_VALUE input parameter to set the
!     flow attenuation saturation level.  The default is 500 Gauss.
!
! 03/01/2023, RC, Version 0.17.1:
!   - BUG FIX:  The WENO3 "alpha" was not taking the max velocity over
!               a wide enough stencil.
!
! 03/14/2023, RC, Version 0.17.2:
!   - BUG FIX:  The WENO3 "weno_eps" was too large for some use cases
!               (specifically coronal hole advection) and could cause
!               ringing. Changed it to be much smaller.
!
! 03/15/2023, RC, Version 0.18.0:
!   - Added ability to redefine the data assimilation weights based on
!     a chosen power of mu.  Set assimilate_data_custom_from_mu=.true.
!     to activate, and set assimilate_data_mu_power=<REAL> to the
!     desired power of mu.
!     Also added assimilate_data_mu_limit input parameter to set
!     a cutoff for mu (default is 0.1).
!   - Added an optional latitude limit on the data assimilation.
!     To use, set assimilate_data_lat_limit=<REAL> to the
!     latitude number of degrees from the pole to zero-out the
!     assilation weights.
!
! 03/24/2023, RC, Version 0.18.1:
!   - Cleaned up some code.
!   - Added writing of namelist parameters to file.
!
! 04/28/2023, RC+MS, Version 0.19.0:
!   - Added ability to specify input file as command line argument.
!     If none specified, the default "hipft.in" is assumed.
!   - Added functionality to output history files for each realization.
!   - Added ability to set multiple values of data assimilation
!     latitude limit over realizations.
!   - Note that the maximum number of realizations is limited to 2000.
!   - Added STOPRUN feature.  Now, if a file is created called "STOPRUN"
!     in the working directory, the run will gracefully stop
!     (better than Ctrl-C!)
!   - Changed default names for output text files to be ".out" and
!     default input file to be hipft.in. This includes history files.
!   - Added new output text file for all timing data.  The summary is
!     output both to the terminal and the file.
!   - Added new output text file that contains realization parameters.
!   - Added experimental auto diffusion subcycle feature based on
!     flux-time-step.  To use, set diffusion_subcycles=0.
!
! 05/04/2023, RC, Version 0.20.0:
!   - Small update to avoid an MPI call.
!   - Added ability to set multiple values of data assimilation
!     mu power and mu limit over realizations
!     (assimilate_data_mu_limits, assimilate_data_mu_powers).
!   - Added ability to set multiple values of flow attenuation
!     over realizations (flow_attenuate_values).
!   - Added ability to set multiple values of diffusion coef constant
!     over realizations (diffusion_coef_constants).
!   - Changed flow attenuation to use maximum of absolute values of Br
!     around staggered velocity component.
!
! 05/24/2023, RC, Version 0.21.0:
!   - Updated auto diffusion subcycle feature based on
!     flux-time-step.  It should be more robust and faster.
!
! 05/31/2023, RC, Version 0.21.1:
!   - BUG FIX:  Automatic diffusion subcycle feature fixed.
!
! 06/02/2023, RC, Version 0.22.0:
!   - BUG FIX:  Fixed incorrect MPI type and duplicate output.
!   - Added output_history_time_cadence parameter to allow user to
!     set a cadence for the history output.
!
! 06/30/2023, RC, Version 0.23.0:
!   - Modified CFL condition in update_timestep to reflect
!     which method one is using to get maximum possible.
!   - Refactored advection routines to have less repeated code.
!
! 07/04/2023, RC, Version 0.24.0:
!   - Added SSPRK(4,3) time stepping scheme for advection.
!     To use, set flow_num_method to "4" (will use WENO3 as well).
!
! 07/10/2023, RC, Version 0.24.1:
!   - Changed the way realization parameters are written out.
!   - BUG FIX:  Flow attenuation values were not being put on GPU,
!               unclear how things were working...
!   - Modified flow output to be in units of m/s.
!
! 07/11/2023, RC, Version 0.25.0:
!   - Updated dtflux caluclation to not allow the time step
!     within the cycles to decrease.
!
! 08/18/2023, RC, Version 0.25.1:
!   - Updated PTL time step calculation to avoid
!     small steps near zero-valued field.
!
! 08/18/2023, RC, Version 0.25.2:
!   - BUGFIX: Seems the flow CFL was being set too high.
!             It needed a 1/2 factor.
!   - Yet another small update to the auto diffusion time step.
!
! 08/29/2023, RC, Version 0.26.0:
!   - Added ability to specify the resolution of the run if one is not
!     using an input map.  This will set the initial map to a uniform
!     grid with zero value.  This is useful for validation tests,
!     and building up a data assimilated map from scratch.
!     To use, set res_nt and res_np to desired resolution.
!     NOTE!  There is currently NO interpolation so you cannot specify
!     a resolution and an input map together.
!
! 09/07/2023, RC, Version 0.27.0:
!   - Updated flow boundary conditions at the poles to be more accurate.
!   - Updated outputs for verbose mode.
!   - Added logical input parameter OUTPUT_SINGLE_PRECISION.
!     Set to .FALSE. to output maps in double precision
!     (good for validation tests).  The default is .TRUE..
!   - Updated automatic diffusion time step (PTL) to newer algorithm.
!
! 10/18/2023, MS+RC, Version 0.28.0:
!   - Added realizations for analytic flow profile parameters for
!     differential rotation and meridianal flows.
!
! 11/22/2023, RC, Version 0.29.0:
!   - Refactored to avoid the need for atomics.
!   - Replaced reduction OpenACC/MP loops with do concurrent reduce.
!     This requires a Fortran 2023 compiler, or preprocessing to work.
!     The proper preprocessing is part of the GCC build scripts.
!   - Replaced OpenACC data movementa nd device selection with
!     OpenMP target data movements and device selection API call.
!     This allows code to compile to Intel GPUs with IFX.
!
! 11/28/2023, RC, Version 0.30.0:
!   - Cleaned up PTL routine to sync it with the version in MAS.
!
! 11/28/2023, RC, Version 0.30.1:
!   - Added back in OpenACC device selection.  The NVIDIA
!     compiler needs this for device selection for DC loops
!     even though it uses the omp version for omp data directives.
!     This should not effect other compilers.
!
! 12/05/2023, RC, Version 0.31.0:
!   - Modified polar boundary conditon loops to use sequential
!     inner loops.  The nested 'do concurrent' loops
!     were giving wrong answers on Intel GPUs with the IFX compiler.
!     Once that compiler is fixed, these loops should be profiled
!     for performance for various sizes of realizations on GPUs/CPUs
!     for nested DC, DC within do, and do within DC.
!     Note:  The NVIDIA compiler parallelizes the sequential do loops
!     without asking, so this change has no effect on NVIDIA GPU runs.
!
! 12/07/2023, RC, Version 1.0.0:
!   - Removed dt_max_increase_fac feature.
!   - Updated version number to reflect first public release.
!
! 12/20/2023, RC, Version 1.0.1:
!   - Added a conditional compilation sentinel in front of the
!     OpenMP device selection call (!$) for added portability.
!
! 01/12/2024, RC, Version 1.0.2:
!   - Bug fix in analytic flows.  The MF flow parameter loading
!     was under a conditional of activating DR flows, not MF.
!
! 01/24/2024, RC+MS, Version 1.1.0:
!   - Updated WENO3 scheme to use cell spacing for the smoother eps.
!     This increases the order of accuracy from 2nd to 3rd in
!     smooth regions near extrema
!     [ See Cravero & Semplice J Sci Comput (2016) 67:1219-1246 ]
!     [ DOI: 10.1007/s10915-015-0123-3                          ]
!
! 01/24/2024, RC+MS, Version 1.2.0:
!   - Made validation run 2 be 3, and set 2 to be blob condition
!     for initial and final output (for phi rigid rotation tests).
!     Validation 3 flips the initial solution for use with
!     symmetric theta blob tests.
!
! 02/02/2024, RC, Version 1.2.1:
!   - Changed defaults to disable Strang splitting.  In our tests,
!     the temporal error is so much smaller than the spatial errors,
!     that the reduced splitting error with Strang is not noticable,
!     and so not worth the extra computation.  It can be turned on
!     in the input file if desired.
!
! 02/08/2024, RC, Version 1.3.0:
!   - Refactored modules to make a bit more sense.
!   - Added ability to use source terms in the model.
!     Use ADVANCE_SOURCE, SOURCE_FROM_FILE, SOURCE_FILENAME,
!     SOURCE_RFE_MODEL, SOURCE_RFE_SIGMA, and
!     SOURCE_RFE_TOTAL_UNSIGNED_FLUX_PER_HOUR to use.
!     Currently, the RFE mode uses a fixed random number sequence
!     that is not updated each step (a planned future feature).
!     The input unit for SOURCE_RFE_TOTAL_UNSIGNED_FLUX_PER_HOUR
!     is 1e21 Mx/hour.
!
! 02/28/2024, RC, Version 1.4.0:
!   - Improved random flux source term implementation.
!   - Added routine to balance flux of a map and use it
!     on random flux source term.
!   - Added source_rfe_seed input parameter to set the random number
!     generator seed.
!   - RFE now generates two frames, interpolates in time between them,
!     and generates a new one when needed.
!   - The length of the time betwen frames is controlled by
!     SOURCE_RFE_LIFETIME (in hours).
!   - The timestep of the code is set to ensure rfe frame switching
!     and is set to half the lifespan.
!   - Note that the current RFE implementation may not be ideal
!     when using multiple realizations or MPI processes.
!     A future update will fix this issue.
!
! 03/04/2024, RC, Version 1.5.0:
!   - Added ability to flux balance the initial map, and/or
!     the output maps.  This is especilly useful when running
!     in a map smoothing/processing mode.
!     to activate set INITIAL_MAP_FLUX_BALANCE=.TRUE.
!     and/or          OUTPUT_MAP_FLUX_BALANCE=.TRUE.
!   - Updated validation test 3 function to exact solution to constant
!     theta velocity.
!   - Updated validation tests 2 and 3 to output errors in history
!     file when used with constant vt and vp velocities.
!   - Added VALIDATION_RUN_WIDTH to input parameters to specify
!     the width of the Gaussians in validation tests 2 and 3.
!   - Validation runs 2 and 3 now auto-set the velocities based on the
!     time-time_start.  For test 3, if a phi rotation is desired,
!     set FLOW_VP_RIGID_OMEGA to any value other than zero for the code
!     to autoset a full rotation.
!
! 03/04/2024, RC, Version 1.6.0:
!   - Fixed bugs preventing GCC compilation.
!   - Added OUTPUT_MAP_IDX_START input parameter to allow runs to
!     begin map output at arbitrary indices.  This helps when
!     restarting a pervious run in combinaton with TIME_START.
!
! 03/05/2024, RC, Version 1.6.1:
!   - Overhauled flows from file feature.  Now, it loops properly,
!     respects the time_start parameter, and is more flexible.
!
! 03/05/2024, RC, Version 1.7.0:
!   - COMPATIBILITY CHANGE: VERBOSE is now an integer.
!     The higher the number, the more info printed (2 is max for now).
!   - Added ability to flux balance the data assimilation window.
!     Use ASSIMILATE_DATA_BALANCE_FLUX=.TRUE. to activate.
!
! 03/06/2024, RC, Version 1.7.1:
!   - Fixed issue with RFE.  The first set of values during setup were
!     not being flux balanced.
!   - Added missing SOURCE_RFE_SEED to the namelist.
!
! 05/27/2024, RC, Version 1.8.0:
!   - Changed when data assimilation takes place and renamed routine.
!     Now, data assimilation takes place at the end of the current step.
!     This is more consistent with the way the time step is set based on
!     data assimilation, avoids double random flux (model+data), and
!     if the output map time is synchronized, makes the data
!     written to the map not have a step applied to it.
!   - Added new RFE model (source_rfe_model=2).  This is similar to
!     the first model, but uses unsigned flux per area per hour instead
!     of unsigned flux per cell per hour for the random sampling.
!     This causes it to have a more unified Br distribution across the
!     map, irrespective of cell size, and the Br is distributed across
!     the sigma much more closely.
!   - Added data assimilation timer.
!
! 05/29/2024, RC/MS/BKJ, Version 1.9.0:
!   - Added RFE unsigned flux as a realization parameter.
!     Use SOURCE_RFE_TOTAL_UNSIGNED_FLUX_PER_HOURS to set.
!   - SOURCE_RFE_SIGMA is deprecated since it was not being applied.
!   - Small OpenMP target bug fix (thanks BKJ!).
!   - Removed mu and sigma from get_gaussian().
!
! 05/30/2024, RC, Version 1.9.1:
!   - Small fix to setting the random number seed.
!
! 06/19/2024, MS/RC, Version 1.10.0:
!   - Removed validation test 2 and replaced it with a more general
!     version of test 3.  Now there are only 2 validation tests.
!     For validation test 2, one can set any nonzero value for
!     the constant Vt and/or constant rigid Vp velocities to auto set
!     the velocities for the preset distances over the run time.
!
! 07/05/2024, RC, Version 1.11.0:
!   - Added back DC to inner boundary condition loops
!     since the Intel compiler has been fixed
!     (see 0.31.0 of 12/05/2023 for details).
!   - Changed PTL numerical check value to 1e-14 instead of 1e-16.
!     Also added it for the delta check below zero.
!     The PTL was failing (giving tiny dt) for HipFT runs
!     at super high resolutions.
!
!-----------------------------------------------------------------------
!
!#######################################################################

