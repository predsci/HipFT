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
      character(*), parameter :: cvers='0.7.1'
      character(*), parameter :: cdate='04/29/2022'
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
      real(r_typ), parameter :: zero=0.0_r_typ
      real(r_typ), parameter :: one=1.0_r_typ
      integer*8,   parameter :: one_int=1
      real(r_typ), parameter :: two=2.0_r_typ
      integer*8,   parameter :: two_int=2
      real(r_typ), parameter :: three=3._r_typ
      integer*8,   parameter :: three_int=3
      real(r_typ), parameter :: four=4._r_typ
      integer*8,   parameter :: four_int=4
      real(r_typ), parameter :: six=6._r_typ
      real(r_typ), parameter :: nine=9._r_typ
      real(r_typ), parameter :: ten=10._r_typ
      real(r_typ), parameter :: fifteen=15._r_typ
      real(r_typ), parameter :: sixteen=16._r_typ
      real(r_typ), parameter :: half=0.5_r_typ
      real(r_typ), parameter :: quarter=0.25_r_typ
      real(r_typ), parameter :: twentyfour=24.0_r_typ
      real(r_typ), parameter :: twentyfive=25.0_r_typ
      real(r_typ), parameter :: fivehundred_i=0.002_r_typ
      real(r_typ), parameter :: three_quarter=0.75_r_typ
      real(r_typ), parameter :: two_third=0.66666666666666666_r_typ
      real(r_typ), parameter :: third=0.33333333333333333_r_typ
!
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
      real(r_typ), parameter :: pi_i=0.3183098861837907_r_typ
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
      real(r_typ), parameter :: twopi_i=0.15915494309189535_r_typ
!
      real(r_typ), parameter :: rsun_cm=6.96e10_r_typ
!
      real(r_typ), parameter :: &
                         diff_km2_s_to_rs2_s=2.0643413925221e-12_r_typ
      real(r_typ), parameter :: &
                         diff_km2_s_to_rs2_hr=7.43162901307967e-09_r_typ
      real(r_typ), parameter :: &
                         km_s_to_rs_hr=0.005172413793103448_r_typ
      real(r_typ), parameter :: &
                         m_s_to_rs_hr=5.172413793103448e-06_r_typ
      real(r_typ), parameter :: &
                         km_s_to_rs_s=1.4367816091954023e-06_r_typ
      real(r_typ), parameter :: output_flux_fac=1.0e-21_r_typ
!
      real(r_typ), parameter :: small_value=tiny(one)
      real(r_typ), parameter :: large_value=huge(one)
      real(r_typ), parameter :: safety=0.95_r_typ
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
! ****** Validation Mode ********
!
      logical        :: validation_run=.false.
!
! ****** Output map cadence ********
!
      integer        :: output_map_idx_cadence = 0
      real(r_typ)    :: output_map_time_cadence = 0.0
!
! ****** Restarts ********
!
      logical        :: restart_run = .false.
      character(512) :: restart_file = ' '
!
! ****** Time ********
!
      real(r_typ)    :: time_end = one
      real(r_typ)    :: time_output_cadence = zero
      real(r_typ)    :: time_restart_cadence = zero
!
! ****** Timestep ********
!
      real(r_typ)    :: dt_max_increase_fac = 0.0_r_typ
      real(r_typ)    :: dt_min = 1.0e-16_r_typ
      real(r_typ)    :: dt_max = huge(one)
!
! ****** Algorithm options.
!
      logical :: strang_splitting = .true.
!
!-----------------------------------------------------------------------
!
! ****** FLOWS ********
!
! ****** Activate the flow advance.
!
      logical :: advance_flow = .false.
!
! ****** Flow files (currently static).
!
      character(512) :: flow_vt_filename = ' '
      character(512) :: flow_vp_filename = ' '
!
! ****** Add a rigid rotation vp velocity (km/s) of omega*sin(theta).
!
      real(r_typ)    :: flow_vp_rigid_omega = 0.
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
! ****** Algorithm options.
!        Can set upwind to central differencing by also setting UPWIND=0.
! ****** 1: FE+Upwind.
! ****** 2: RK3TVD+Upwind.
! ****** 3: RK3TVD+WENO3.
!
      integer :: flow_num_method = 2
!
! ****** Upwind coefficient.
!
      real(r_typ) :: upwind = 1._r_typ
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
! ****** For diffusion-only runs with RKG2/(RKL2), set to ~30/(60).
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
      real(r_typ), dimension(:,:), allocatable :: fout
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
      real(r_typ) :: h_minbr, h_maxbr, h_minabsbr, &
                     h_fluxp, h_fluxm, h_valerr
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
!  [RMC] Make f->br??
      real(r_typ), dimension(:,:), allocatable :: f
      real(r_typ), dimension(:,:), allocatable :: fold
      real(r_typ), dimension(:,:), allocatable :: fval_u0
      real(r_typ), dimension(:,:), allocatable :: diffusion_coef
      real(r_typ), dimension(:,:), allocatable :: source
      real(r_typ), dimension(:,:), allocatable :: vt
      real(r_typ), dimension(:,:), allocatable :: vp
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
      real(r_typ), dimension(:), allocatable :: &
                     map_times_requested_ut_jd, map_times_actual_ut_jd
!
      character(19), dimension(:), allocatable :: &
                     map_times_requested_ut_str, map_times_actual_ut_str
!
      character(512), dimension(:), allocatable :: map_files_rel_path
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
      real(r_typ), dimension(:,:), allocatable :: u0
      real(r_typ), dimension(:,:), allocatable :: dty0
      real(r_typ), dimension(:,:), allocatable :: ykm1
      real(r_typ), dimension(:,:), allocatable :: ukm1
      real(r_typ), dimension(:,:), allocatable :: ukm2
!
end module
!#######################################################################
module matrix_storage
!
      use number_types
!
      implicit none
!
      real(r_typ), dimension(:,:,:), allocatable :: coef
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
!
      real(r_typ) :: wtime_setup = 0.
      real(r_typ) :: wtime_update = 0.
      real(r_typ) :: wtime_flux_transport = 0.
      real(r_typ) :: wtime_flux_transport_advection = 0.
      real(r_typ) :: wtime_flux_transport_diffusion = 0.
      real(r_typ) :: wtime_source = 0.
      real(r_typ) :: wtime_analysis = 0.
      real(r_typ) :: wtime_io = 0.
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
        real(r_typ), dimension(:,:,:), allocatable :: new_data
        end subroutine
      end interface
end module
!#######################################################################
program HIPFT
!
!-----------------------------------------------------------------------
!
      use number_types
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
      logical :: the_run_is_done
      real(r_typ) :: wtime
      character(*), parameter :: FMT = '(a,i8,a,f12.6,a,f12.6,a,f12.6)'
!
!-----------------------------------------------------------------------
!
! ****** Start wall clock timer.
!
      wtime_tmp = wtime()
!
! ****** Set up and initialize the run.
!
      call setup
!
! ****** START MAIN LOOP ********************************************
!
      write (*,*)
      write (*,*) '>running'
      write (*,*)
      flush(OUTPUT_UNIT)
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
        if (verbose) then
          write(*,FMT) 'Starting step: ',ntime,'  Time:',time,' /', &
                        time_end,'  dt to use:',dtime_global
          flush(OUTPUT_UNIT)
        endif
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
        write(*,FMT) 'Completed step: ',ntime,'  Time:',time,' /', &
                      time_end,'  dt used:',dtime_global
        flush(OUTPUT_UNIT)
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
      wtime_total = wtime() - wtime_tmp
!
! ****** Write out time profiling.
!
      call write_timing
!
end program HIPFT
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
      use input_parameters
      use timing
      use fields
      use mesh
      use output
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: wtime
      integer :: ierr
      character(512) :: cmd
!
!-----------------------------------------------------------------------
!
      wtime_tmp = wtime()
!
      call read_input
!
! ****** Set up output directory.
!
      if (output_map_idx_cadence.gt.0 .or.     &
          output_map_time_cadence.gt.0.0) then
        cmd = 'mkdir -p '//trim(output_map_directory)
        call EXECUTE_COMMAND_LINE (cmd,exitstat=ierr)
        if (ierr.ne.0) then
          write (*,*) '### Could not make output subdirectory, using ./'
          output_map_directory = '.'
        end if
      end if
!
      if (restart_run) then
!        call load_restart
         print*,'RESTARTS ARE NOT IMPLEMENTED YET.'
         stop
      else
        call load_initial_condition
      end if
!
! ****** Allocate flow arrays.
!
      if (advance_flow) then
        allocate (vt(nt,npm))
        allocate (vp(ntm,np))
        vt(:,:) = 0.
        vp(:,:) = 0.
!$acc enter data copyin(vt,vp)
      end if
!
      call set_unchanging_quantities
!
      if (assimilate_data) call load_data_assimilation

      call create_and_open_output_log_files
!
      call analysis_step
!
      if (output_map_idx_cadence .gt. 0 .or. &
         output_map_time_cadence .gt. 0.0) then
        output_current_map = .true.
      end if 
!
      call output_step
!
      call write_welcome_message
!
      wtime_setup = wtime() - wtime_tmp
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
      use input_parameters
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Write the final Br map.
!
      call write_map (trim(output_map_root_filename)//'_final.h5')
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
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*), parameter :: FMT = '(a20,f20.2)'
!
!-----------------------------------------------------------------------
!
      write(*,"(a40)") repeat("-", 40)
      write(*,FMT) "Wall clock time:   ",wtime_total
      write(*,"(a40)") repeat("-", 40)
      write(*,FMT) "--> Setup:         ",wtime_setup
      write(*,FMT) "--> Update:        ",wtime_update
      write(*,FMT) "--> Flux transport:",wtime_flux_transport
      write(*,FMT) "    --> Advecton:  ",wtime_flux_transport_advection
      write(*,FMT) "    --> Diffusion  ",wtime_flux_transport_diffusion
      write(*,FMT) "    --> Source:    ",wtime_source
      write(*,FMT) "--> Analysis:      ",wtime_analysis
      write(*,FMT) "--> I/O:           ",wtime_io
      write(*,"(a40)") repeat("-", 40)
!
      write(*,*)
      flush(OUTPUT_UNIT)
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
      write (IO_HIST_NUM,'(a10,a,6(a22,a),a15)')                      &
      'STEP',' ','TIME',' ','DTIME',' ','DTIME_ADV_STB',' ',          &
      'DTIME_ADV_USED',' ','DTIME_DIFF_STB',' ','DTIME_DIFF_USED',    &
      ' ','N_DIFF_PER_STEP'
!
      close(IO_HIST_NUM)
!
      call FFOPEN (IO_HIST_SOL,io_hist_sol_filename,'rw',ierr)
!
      write (IO_HIST_SOL,'(a10,a,7(a22,a))')                          &
      'STEP',' ','TIME',' ','BR_MIN',' ','BR_MAX',' ','BR_ABS_MIN',   &
      ' ','FLUX_POSITIVE',' ','FLUX_NEGATIVE',' ',                    &
      'VALIDATION_ERR_CVRMSD'
!
      close(IO_MAP_OUT_LIST)
!
      if (output_map_idx_cadence.gt.0 .or.     &
          output_map_time_cadence.gt.0.0) then
!
        call FFOPEN (IO_MAP_OUT_LIST,io_map_output_list_filename,'rw',ierr)
!
        write (IO_MAP_OUT_LIST,'(a10,a,a25,a,a)') 'IDX',' ','TIME',' ','FILENAME'
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
      write (*,*) '                  Jon A. Linker'
      write (*,*) '                  Zoran Mikic'
      write (*,*) ''
      write (*,*) '        Predictive Science Inc.'
      write (*,*) '        www.predsci.com'
      write (*,*) '        San Diego, California, USA 92121'
      write (*,*)''
      write (*,*)''
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
      use input_parameters
      use mesh
      use fields
      use globals
      use read_2d_file_interface
      use output
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: fn1,fs1,fn2,fs2
      real(r_typ), dimension(:,:), allocatable :: f_local
      real(r_typ), dimension(:), allocatable :: s1,s2
      integer :: n1,n2,ierr
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
      call read_2d_file (initial_map_filename,n1,n2,f_local,s1,s2,ierr)
!
      if (validation_run) then
!
        allocate (fval_u0(n1,n2))
!
! ****** Make initial solution f and output.
!
        call tesseral_spherical_harmonic (0,6,s1,s2,n1,n2,f_local)
        call tesseral_spherical_harmonic (5,6,s1,s2,n1,n2,fval_u0)
!
        fval_u0(:,:) = 1000.0_r_typ*(f_local(:,:) +                  &
                              sqrt(14.0_r_typ/11.0_r_typ)*fval_u0(:,:))
!
        call write_2d_file(                                           &
            (trim(output_map_root_filename)//'_initial_0.h5')       &
                           ,n1,n2,fval_u0,s1,s2,ierr)

        f_local(:,:) = fval_u0(:,:)
!
! ****** Caculate final analytic solution and output.
!
        allocate (fout(n1,n2))
        fout(:,:) = fval_u0(:,:)*exp(-42.0_r_typ*diffusion_coef_constant* &
                                     diffusion_coef_factor*time_end)
!
        call write_2d_file(                                           &
             (trim(output_map_root_filename)//'_final_analytic.h5') &
                          ,n1,n2,fout,s1,s2,ierr)
!
        deallocate (fout)
!
!$acc enter data copyin(fval_u0)
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
      allocate (f(ntm1,npm1 + 1))
!
! ****** Transpose array since we assume PT but need TP for the code.
!
      f(1:ntm1,1:npm1) = TRANSPOSE(f_local(1:npm1,1:ntm1))
!
! ****** Impose perfect periodicity and set two-point overlap.
!
      f(:,1) = half*(f(:,1) + f(:,npm1))
      f(:,npm1) = f(:,1)
      f(:,npm1 + 1) = f(:,2)
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
      call get_m0 (f,fn1,fn2,fs1,fs2)
!
! ****** Set the pole values to have only an m=0 component.
!
      f(1,:) = fn1
      f(ntm,:) = fs1
!$acc enter data copyin(f)
!
! ****** Clean up memory.
!
      deallocate (s1)
      deallocate (s2)
      deallocate (f_local)
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
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
!
      call output_histories
!
      call output_map
!
!     call output_flow(?)
!
!     call output_restart
!
      wtime_io = wtime_io + (wtime() - t1)
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
      character(*), parameter :: FMT='(i10,a,6(1pe22.15,a),i15)'
      character(*), parameter :: FMT2='(i10,a,7(1pe22.15,a))'
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

      write(IO_HIST_SOL,FMT2) ntime,' ',                              &
                             time,' ',h_minbr,' ',h_maxbr,' ',        &
                             h_minabsbr, ' ',h_fluxp,' ',h_fluxm,     &
                             ' ',h_valerr
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
! ****** File sequence number.
!
      integer*8, save :: idx = 0
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
        if (time .ge. idx_out*output_map_time_cadence) then
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
      use number_types
      use output
      use fields, ONLY : f
      use mesh, ONLY : t,npm,ntm
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      integer :: ierr = 0
!
!-----------------------------------------------------------------------
!
!$acc update self(f)
!
! ****** Write out file in pt format with single point overlap in phi.
!
      allocate (fout(npm-1,ntm))
!
!   WARNING: -stdpar may eventually GPU-ize TRANSPOSE.
!   For now, it probably in-lines so OK, and/or requires
!   the tensor module to be explicitly used.
!   For now, this should be computed on the CPU no matter what...
!
      fout(:,:) = TRANSPOSE(f(:,1:npm-1))
!
      call write_2d_file (fname,npm-1,ntm,fout,pout,t,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_MAP:'
        write (*,*) '### Could not write the output data set!'
        write (*,*) 'IERR: ',ierr
        write (*,*) 'File name: ', fname
        call exit (1)
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
      real(r_typ) :: t1,wtime
      real(r_typ) :: sumfval2,fv,fv2
      real(r_typ), dimension(:,:), allocatable :: fval
      integer :: i,j
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
!
! ****** If running validation, compute current exact solution.
!
      if (validation_run) then
        allocate (fval(npm-1,ntm))
!$acc enter data create (fval)
!
        do concurrent (i=1:npm-1,j=1:ntm)
          fval(i,j) = fval_u0(i,j)*exp(-42.0_r_typ*diffusion_coef_constant* &
                                       diffusion_coef_factor*time)
        enddo
!
      end if
!
! ***** Get metrics of the current field.
!
      h_fluxp = 0.
      h_fluxm = 0.
      h_minbr = large_value
      h_maxbr = small_value
      h_minabsbr = large_value
      h_valerr = 0.
      sumfval2 = 0.
!
!$omp parallel do collapse(2) default(shared) &
!$omp         reduction(+:sumfval2,h_valerr) &
!$omp         reduction(max:h_maxbr) reduction(min:h_minbr,h_minabsbr)
!$acc parallel loop collapse(2) default(present) &
!$acc         reduction(+:sumfval2,h_valerr) &
!$acc         reduction(max:h_maxbr) reduction(min:h_minbr,h_minabsbr)
      do j=1,npm-2
        do i=1,ntm
          h_minbr = min(f(i,j),h_minbr)
          h_maxbr = max(f(i,j),h_maxbr)
          h_minabsbr = min(abs(f(i,j)),h_minabsbr)
          if (validation_run) then
            fv = (f(i,j) - fval(j,i))**2
            fv2 = fval(j,i)**2
          else
            fv = 0.
            fv2 = 0.
          end if
          h_valerr = h_valerr + fv
          sumfval2 = sumfval2 + fv2
        enddo
      enddo
!$omp end parallel do
!
      if (validation_run) then
        h_valerr = sqrt(h_valerr/sumfval2)
!$acc exit data delete (fval)
      end if
!
! ****** Get flux in units of 1e21 Mx
!
      call get_flux (h_fluxp,h_fluxm)
      h_fluxp = output_flux_fac*h_fluxp*rsun_cm**2
      h_fluxm = output_flux_fac*h_fluxm*rsun_cm**2
!
      wtime_analysis = wtime_analysis + (wtime() - t1)
!
end subroutine
!#######################################################################
subroutine get_flux (fluxp,fluxm)
!
!-----------------------------------------------------------------------
!
! ****** Calculate the total positive and negative flux.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use globals
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
      real(r_typ) :: fluxp,fluxm
!
!-----------------------------------------------------------------------
!
      integer :: i,j
      real(r_typ) :: tav,da_t,da_p
      real(r_typ) :: fv,area_check
!
!-----------------------------------------------------------------------
!
      fluxp = 0.
      fluxm = 0.
      area_check = 0.
!
!$omp parallel do collapse(2) default(shared) &
!$omp         reduction(+:fluxp,fluxm,area_check)
!$acc parallel loop collapse(2) default(present) &
!$acc         reduction(+:fluxp,fluxm,area_check)
      do i=1,npm-1
        do j=1,ntm
          if (j.eq.1) then
            tav=half*(t(1)+t(2))
            da_t=quarter*sin(tav)*dth(2)
          else if (j.eq.ntm) then
            tav=half*(t(ntm)+t(ntm-1))
            da_t=quarter*sin(tav)*dth(ntm)
          else
            da_t=sin(t(j))*dt(j)
          end if
!
          if (i.eq.1) then
            da_p=half*dph(1)
          else if (i.eq.npm-1) then
            da_p=half*dph(npm-1)
          else
            da_p=dp(i)
          end if
!
          fv = f(j,i)*da_t*da_p
!
          if (fv.gt.0.) then
            fluxp = fluxp + fv
          else
            fluxm = fluxm + fv
          end if
!
          area_check = area_check + da_t*da_p
        enddo
      enddo
!$omp end parallel do
!
!      print*,'ZMFLUX:'
!      print*,'Area int:   ', area_check*1.0e-21_r_typ*rsun_cm**2
!      print*,'Area exact: ', 1.0e-21_r_typ*four*pi*rsun_cm**2
!      print*,'Flux+: ',fluxp
!      print*,'Flux-: ',fluxm
!
end subroutine
!#######################################################################
subroutine get_flux_trap (fluxp,fluxm)
!
!-----------------------------------------------------------------------
!
! ****** Calculate the total positive and negative flux.
!
!-----------------------------------------------------------------------
!
      use input_parameters
      use globals
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
      real(r_typ) :: fluxp,fluxm
!
!-----------------------------------------------------------------------
!
      integer :: j,k
      real(r_typ) :: fav,favp,favm,area_check
      real(r_typ), dimension(:), allocatable :: int_tmp_p, int_tmp_m, area_chk
!
!-----------------------------------------------------------------------
!
      fluxp = 0.
      fluxm = 0.
      area_check = 0.
!
!$acc update self(f)
!
      allocate (int_tmp_p(npm-1))
      allocate (int_tmp_m(npm-1))
      allocate (area_chk(npm-1))
      int_tmp_p(:)=0.
      int_tmp_m(:)=0
      area_chk(:)=0.
!
! ****** First, do integration along theta direction.
!
      do j=1,ntm-1
        do k=1,npm-1
          fav = half*(f(j+1,k) + f(j,k))
          if (fav.gt.0.) then
            int_tmp_p(k) = int_tmp_p(k) + fav*dth(j+1)*sth(j+1)
          else
            int_tmp_m(k) = int_tmp_m(k) + fav*dth(j+1)*sth(j+1)
          end if
          area_chk(k) = area_chk(k) + dth(j+1)*sth(j+1)
        enddo
      enddo
!
! ****** Now, do integration along phi direction.
!
      do k=1,npm-2
        favp = half*(int_tmp_p(k+1) + int_tmp_p(k))
        favm = half*(int_tmp_m(k+1) + int_tmp_m(k))
        fluxp = fluxp + favp*dph(k+1)
        fluxm = fluxm + favm*dph(k+1)
        area_check = area_check + area_chk(k)*dph(k+1)
      enddo
!
      fluxp = 1.0e-21_r_typ*fluxp*rsun_cm**2
      fluxm = 1.0e-21_r_typ*fluxm*rsun_cm**2
      area_check = area_check*1.0e-21_r_typ*rsun_cm**2
!
      deallocate (int_tmp_p)
      deallocate (int_tmp_m)
      deallocate (area_chk)

      print*,'TRAPZ:'
      print*,'Area int:   ', area_check
      print*,'Area exact: ', 1.0e-21_r_typ*four*pi*rsun_cm**2
      print*,'Flux+: ',fluxp
      print*,'Flux-: ',fluxm
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
      use input_parameters
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
!
      if (advance_flow)    call update_flow
      if (advance_source)  call update_source
      if (assimilate_data) call update_field
!
      call update_timestep
!
      wtime_update = wtime_update + (wtime() - t1)
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
      integer :: j,k
      real(r_typ) :: br_avg
      real(r_typ) :: vt_npole,vt_spole,vp_npole,vp_spole
!
!-----------------------------------------------------------------------
!
! ***** Reset the flow.
!
      if (flow_needs_updating) then
!
!       Assume no more updates needed.  If time-dept flow,
!       this will be set to true each time after this.
        flow_needs_updating = .false.
!       Since the flow is changing, we need to update time step.
        timestep_flow_needs_updating = .true.
        timestep_needs_updating = .true.
!
        do concurrent (k=1:npm,j=1:nt)
          vt(j,k) = 0.
        end do
!        
        do concurrent (k=1:np,j=1:ntm)
          vp(j,k) = 0.
        end do        
!
! ***** Add in flow from file (this should save current state and check
!       if a new one is needed/interp to avoid tons of I/O.
!
        call add_flow_from_file
!
! ***** Add in flow from ConFlow algorithm.
!
!       call add_flow_from_conflow (metadata_structure,flow_array)
!
! ***** Add in meridianal flow model.
!
        if (flow_mf_model.eq.1) call add_flow_meridianal_aft
!
! ***** Add in differential rotation model.
!
        if (flow_dr_model.eq.1) call add_flow_differential_rotation_aft
!
! ***** Add in constant angular velocity (rigid rotation)
!
        if (flow_vp_rigid_omega.gt.0.) then
          do concurrent (k=1:np,j=1:ntm)
            vp(j,k) = vp(j,k)+flow_vp_rigid_omega*km_s_to_rs_hr*st(j)
          end do
        endif
!
! ***** Attenuate velocity based on Br.
!
        if (flow_attenuate) then
!
!
! ****** Need to interpolate f at half-mesh point.
!
          do concurrent (k=2:npm-1,j=2:nt-1)
            br_avg = half*(f(j-1,k) + f(j,k))
            vt(j,k) = vt(j,k)*(one - tanh(abs(br_avg)*fivehundred_i))
          end do
          
          do concurrent (k=2:np-1,j=2:ntm-1)
            br_avg = half*(f(j,k-1) + f(j,k))
            vp(j,k) = vp(j,k)*(one - tanh(abs(br_avg)*fivehundred_i))
          end do
!
! ****** Set poles. Maybe not needed since field is weak at pole?
!
          vt_npole = 0.
          vt_spole = 0.
          vp_npole = 0.
          vp_spole = 0.
!
!$acc parallel loop present(vt) reduction(+:vt_npole,vt_spole)
          do k=2,npm-1
            vt_npole = vt_npole + vt(   2,k)*dp(k)
            vt_spole = vt_spole + vt(ntm1,k)*dp(k)
          enddo
!
!$acc parallel loop present(vp) reduction(+:vp_npole,vp_spole)
          do k=2,npm1
            vp_npole = vp_npole + vp(    2,k)*dph(k)
            vp_spole = vp_spole + vp(ntm-1,k)*dph(k)
          enddo  
!
          vt_npole = vt_npole*twopi_i
          vt_spole = vt_spole*twopi_i
          vp_npole = vp_npole*twopi_i
          vp_spole = vp_spole*twopi_i

! ****** Attenuate pole value to set new v bc
!        (f should have same value for all polar points so can use 1)

!$acc update host(f(1,1),f(ntm,1))
          vt_npole = vt_npole*(one - tanh(f(1,1)*fivehundred_i))
          vt_spole = vt_spole*(one - tanh(f(ntm,1)*fivehundred_i))
          vp_npole = vp_npole*(one - tanh(f(1,1)*fivehundred_i))
          vp_spole = vp_spole*(one - tanh(f(ntm,1)*fivehundred_i))
!
! ****** Now set the pole:
!
          do concurrent (k=1:npm)
            vt( 1,k) = two*vt_spole-vt(   2,k)
            vt(nt,k) = two*vt_spole-vt(ntm1,k)
          end do
!
          do concurrent (k=1:np)
            vp(  1,k) = two*vp_spole-vp(    2,k)
            vp(ntm,k) = two*vp_spole-vp(ntm-1,k)
          end do  
!
! ****** Set periodicity
!
          call set_periodic_bc_2d (vt,nt,npm)
          call set_periodic_bc_2d (vp,ntm,np)
!
! ****** Since this needs to be done each timestep as Br changes,
! ****** the full velocity profile needs to be reloaded and
! ****** re-attenuated.
!
          flow_needs_updating = .true.
!
        end if
!
      end if
!
end subroutine
!#######################################################################
subroutine add_flow_from_file
!
!-----------------------------------------------------------------------
!
! ****** Update flows from file(s).
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
  ! [] THESE DO NOT WORK WITH GPU AND SHOULD BE ADDING TO FLOW
  ! [] (NOT OVERWRITE)  NEED TO FIX THIS!
      if (flow_vt_filename.ne.' ') call load_vt
!
      if (flow_vp_filename.ne.' ') call load_vp
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
                             FMT = '(F11.5,1X,F11.5,1X,A19,1X,A19,1X,A)'
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
      end do
      REWIND (IO_DATA_IN)
!
! ****** Allocate space to store the list.
!
      allocate(map_times_requested_ut_jd  (num_maps_in_list))
      allocate(map_times_actual_ut_jd     (num_maps_in_list))
      allocate(map_times_requested_ut_str (num_maps_in_list))
      allocate(map_times_actual_ut_str    (num_maps_in_list))
      allocate(map_files_rel_path         (num_maps_in_list))
!
! ****** Read the list (skip first row).
!
      READ (IO_DATA_IN,*)
      do i=1,num_maps_in_list
        READ (IO_DATA_IN,FMT) map_times_requested_ut_jd(i),  &
                              map_times_actual_ut_jd(i),     &
                              map_times_requested_ut_str(i), &
                              map_times_actual_ut_str(i),    &
                              map_files_rel_path(i)
      enddo
      CLOSE (IO_DATA_IN)
!
! ******* Initialize assimilation time.
!
      map_time_initial_hr = twentyfour*map_times_requested_ut_jd(1)
!
      current_map_input_idx = 1
!
      time_of_next_input_map = 0.
!
! ******* Print some diagnostics.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'DATA ASSIMILATION'
        write (*,*) 'Number of maps in map list: ', &
                     num_maps_in_list
        write (*,*) 'Start date (actual):     ', &
                     trim(map_times_actual_ut_str(1))
        write (*,*) 'End   date (actual):     ', &
                     trim(map_times_actual_ut_str(num_maps_in_list))
        write (*,*) 'File name of first map:     ', &
                     trim(map_files_rel_path(1))
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
      integer :: j,k
      real(r_typ), dimension(:,:,:), allocatable :: new_data
!
!-----------------------------------------------------------------------
!
      do concurrent (k=1:npm1,j=1:ntm)
         f(j,k) = (one - new_data(k,j,2))*       f(j,k) +       &
                  (      new_data(k,j,2))*new_data(k,j,1)
      enddo
!
! ****** Set periodicity (since k=np (npm) not set above).
!
      call set_periodic_bc_2d (f,ntm,npm)
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
      use globals
      use output
      use read_3d_file_interface
      use assimilate_new_data_interface
      use data_assimilation
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr = 0
      real(r_typ), dimension(:,:,:), allocatable :: new_data
      real(r_typ), dimension(:), allocatable :: s1,s2,s3
      integer :: ln1,ln2,nslices
      character(1024) :: mapfile = ' '
!
!-----------------------------------------------------------------------
!
     if (time .ge. time_of_next_input_map) then
!
! ****** Read the map data.
!
       mapfile = TRIM(assimilate_data_map_root_dir)//"/"&
                 //TRIM(map_files_rel_path(current_map_input_idx))
!
       call read_3d_file (mapfile,ln1,ln2,nslices,new_data,s1,s2,s3,ierr)
!
! ******  TODO:  Add error checking here (make sure scales match run
!                until interp is added)
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
       deallocate(s1)
       deallocate(s2)
       deallocate(s3)
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
          twentyfour*map_times_requested_ut_jd(current_map_input_idx) &
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
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dtime_old, dt_mxup
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
!
        if (advance_flow.and.timestep_flow_needs_updating) then
          call get_flow_dtmax (dtmax_flow)
          dtime_global = MIN(dtime_global,dtmax_flow)
          timestep_flow_needs_updating = .false.
        end if
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
        if (verbose) write (*,*) 'Stable timestep: ',dtime_global
        timestep_needs_updating = .false.
!
      end if
!
! ****** Check for next output time, cut dt to match exactly.
!
      if (output_map_time_cadence .gt. 0.0) then
        if (time + dtime_global .ge. idx_out*output_map_time_cadence) then
          dtime_global = idx_out*output_map_time_cadence - time
          timestep_needs_updating = .true.
        end if
      end if
!
! ****** Check for next map assimilation time, cut dt to match exactly.
!
      if (assimilate_data) then
        if (time + dtime_global .ge. time_of_next_input_map) then
          dtime_global = time_of_next_input_map - time
          timestep_needs_updating = .true.
        end if
      end if
!
! ****** Check for end time.
!
      if (time + dtime_global .ge. time_end) then
        dtime_global = time_end - time
        timestep_needs_updating = .false.
      end if
!
! ****** Re-compute Euler diffusion time-step if needed.
!
      if (advance_diffusion.and.timestep_diff_needs_updating) then
        call get_dtime_diffusion_euler (dtime_diffusion_euler)
        timestep_diff_needs_updating = .false.
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
      real(r_typ) :: dtime_local, dtime_local_half, t1, wtime
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
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
      wtime_flux_transport = wtime_flux_transport + (wtime() - t1)
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
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
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
                                       + (wtime() - t1)
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
      real(r_typ) :: dtime_local2,t1,wtime
      real(r_typ) :: time_stepped,timetmp
      integer :: i
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
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
          call output_histories
          time = timetmp
        end if
!
      enddo
!
      wtime_flux_transport_diffusion = wtime_flux_transport_diffusion &
                                       + (wtime() - t1)
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
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
!
! ---> Add random flux, ARs, etc.
!
      wtime_source = wtime_source + (wtime() - t1)
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
      integer :: j,k
!
!-----------------------------------------------------------------------
!
! ****** Allocate coef array.
!
      allocate (coef(2:ntm-1,2:npm-1,5))
      coef(:,:,:)=0.
!$acc enter data copyin(coef)
!
! ****** Set coef for internal points and phi boundary point at k=1.
!
      do concurrent (j=2:ntm-1,k=2:npm-1)
!
        coef(j,k,1) = half*(diffusion_coef(j  ,k  )  &
                          + diffusion_coef(j+1,k  )) &
                     *dph_i(k  )*dp_i(k)*st_i(j)*st_i(j)
!
        coef(j,k,2) = half*(diffusion_coef(j  ,k  )  &
                          + diffusion_coef(j  ,k+1)) &
                     *dth_i(j  )*dt_i(j)*st_i(j)*sth(j  )
!
        coef(j,k,4) = half*(diffusion_coef(j+1,k  )  &
                          + diffusion_coef(j+1,k+1)) &
                     *dth_i(j+1)*dt_i(j)*st_i(j)*sth(j+1)
!
        coef(j,k,5) = half*(diffusion_coef(j  ,k+1)  &
                          + diffusion_coef(j+1,k+1)) &
                     *dph_i(k+1)*dp_i(k)*st_i(j)*st_i(j)
!
        coef(j,k,3)=-(coef(j,k,1)+coef(j,k,2)+coef(j,k,4)+coef(j,k,5))
!
      enddo
!
end subroutine
!#######################################################################
subroutine get_dtime_diffusion_euler (dtime_exp)
!
!-----------------------------------------------------------------------
!
! ****** Get the explicit Euler time step limit for thermal conduction.
!
!-----------------------------------------------------------------------
!
      use matrix_storage
      use mesh
      use constants
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
      integer :: j,k,d
      real(r_typ) :: max_eig,gersh_rad
!
!-----------------------------------------------------------------------
!
! *** Estimate maximum eigenvalue using Gershgorin disks:
!
      max_eig = 0.
!
!$omp parallel do default(shared) collapse(2) reduction(max:max_eig)
!$acc parallel loop default(present) collapse(2) reduction(max:max_eig)
      do k=2,npm-1
        do j=2,ntm-1
            gersh_rad = 0.
!$acc loop seq
            do d=1,5
              gersh_rad = gersh_rad+abs(coef(j,k,d))
            enddo
            max_eig = max(gersh_rad,max_eig)
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
      integer :: i,j
      integer*8 :: k
      real(r_typ), INTENT(IN) :: dtime_local
!
!-----------------------------------------------------------------------
!
!     This only needs to happen more than once if dtime_local changes...
!
      if (need_to_load_sts) then
        if (diffusion_num_method.eq.2) then
          call load_sts_rkl2 (dtime_local)
        elseif(diffusion_num_method.eq.3) then
          call load_sts_rkg2 (dtime_local)
        endif
        need_to_load_sts = .false.
      end if
!
! ****** Allocate scratch arrays.
!
      allocate (u0(ntm,npm))
      allocate (dty0(ntm,npm))
      allocate (ykm1(ntm,npm))
      allocate (ukm1(ntm,npm))
      allocate (ukm2(ntm,npm))
!$acc enter data create(u0,dty0,ykm1,ukm1,ukm2)
!
      call ax(f,dty0)
!
      do concurrent (j=1:npm,i=1:ntm)
        u0(i,j) = f(i,j)
        ukm2(i,j) = f(i,j)
        dty0(i,j) = dtime_local*dty0(i,j)
        ukm1(i,j) = f(i,j) + sts_ubj(1)*dty0(i,j)
      enddo
!
! ****** Inner s-step loop
!
      do k=2,sts_s
!
        call ax(ukm1,ykm1)
!
        do concurrent (j=1:npm,i=1:ntm)
          f(i,j) =               sts_uj(k)*ukm1(i,j) + &
                                 sts_vj(k)*ukm2(i,j) + &
                   (one-sts_uj(k)-sts_vj(k))*u0(i,j) + &
                    sts_ubj(k)*dtime_local*ykm1(i,j) + &
                                 sts_gj(k)*dty0(i,j)
!
          ukm2(i,j) = ukm1(i,j)
          ukm1(i,j) = f(i,j)
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
      endif
!
! ****** Make sure s is odd.
!
      if (MOD(sts_s,2).eq.0) then
        sts_s = sts_s + 1
      endif
!
! ****** Allocate super-time-step coefficent arrays.
!
      if (need_to_load_sts.and..not.first) then
!$acc exit data delete(sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
        deallocate (sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
      end if
!
      if (need_to_load_sts) then
        allocate (sts_uj(sts_s))
        allocate (sts_vj(sts_s))
        allocate (sts_ubj(sts_s))
        allocate (sts_gj(sts_s))
        allocate (sts_b(sts_s))
!$acc enter data create(sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
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
!$acc update device(sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
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
      endif
!
! ****** Make sure s is odd.
!
      if (MOD(sts_s,2).eq.0) then
        sts_s = sts_s + 1
      endif
!
! ****** Allocate super-time-step coefficent arrays.
!
      if (need_to_load_sts.and..not.first) then
!$acc exit data delete(sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
        deallocate (sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
      end if
!
      if (need_to_load_sts) then
        allocate (sts_uj(sts_s))
        allocate (sts_vj(sts_s))
        allocate (sts_ubj(sts_s))
        allocate (sts_gj(sts_s))
        allocate (sts_b(sts_s))
!$acc enter data create(sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
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
!$acc update device(sts_uj,sts_vj,sts_ubj,sts_gj,sts_b)
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
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
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
        call exit (1)
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
        call exit (1)
      end if
!
! ****** Check that the their are scales.
!
      if (.not.s%scale) then
        write (*,*)
        write (*,*) '### ERROR in READ_2D_FILE:'
        write (*,*) '### The data set does not contain scales.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
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
      wtime_io = wtime_io + (wtime() - t1)
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
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
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
        call exit (1)
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
        call exit (1)
      end if
!
! ****** Check that the their are scales.
!
      if (.not.s%scale) then
        write (*,*)
        write (*,*) '### ERROR in READ_3D_FILE:'
        write (*,*) '### The data set does not contain scales.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
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
      wtime_io = wtime_io + (wtime() - t1)
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
      endif
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
      endif
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
      end do
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
            end do
            deallocate (f4dim)
          elseif (prec.eq.64) then
            allocate (f8dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f8dim,s_dims_i,ierr)
            do j=1,s%dims(i)
              s%scales(i)%f(j) = REAL(f8dim(j),r_typ)
            end do
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
         endif
      end do
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
      endif
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
          endif
          s_dims_i = s_dims(i)
          call h5Screate_simple_f(1,s_dims_i,dspacedim_id,ierr)
          if (s%hdf32) then
            allocate (f4dim(s_dims_i(1)))
            do j=1,s%dims(i)
              f4dim(j) = REAL(s%scales(i)%f(j),REAL32)
            end do
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_REAL, &
                              dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_REAL, &
                             f4dim,s_dims_i,ierr)
            deallocate (f4dim)
          else
            allocate (f8dim(s_dims_i(1)))
            do j=1,s%dims(i)
              f8dim(j) = REAL(s%scales(i)%f(j),REAL64)
            end do
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_DOUBLE, &
                             dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_DOUBLE, &
                             f8dim,s_dims_i,ierr)
            deallocate (f8dim)
          endif
          call h5DSset_scale_f (dim_id,ierr,dimname)
          call h5DSattach_scale_f (dset_id,dim_id,i,ierr)
          call h5DSset_label_f(dset_id, i, dimname, ierr)
          call h5Dclose_f (dim_id,ierr)
          call h5Sclose_f (dspacedim_id,ierr)
        end do
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in WRH5:'
          write (*,*) '### Could not write the scales.'
          write (*,*) 'File name: ',trim(fname)
          ierr = 5
          return
        endif
      endif
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
      real(r_typ) :: t1,wtime
!
!-----------------------------------------------------------------------
!
      type(ds) :: s
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      t1 = wtime()
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
      wtime_io = wtime_io + (wtime() - t1)
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
        call exit (1)
      end if
!
! ****** Allocate main mesh grid quantities.
!
      allocate (dt(ntm))
      allocate (dt_i(ntm))
      allocate (st(ntm))
      allocate (st_i(ntm))
      allocate (ct(ntm))
      allocate (dp(npm))
      allocate (dp_i(npm))
!
! ****** Allocate half mesh grid quantities.
!
      allocate (th(nt))
      allocate (ph(np))
      allocate (dth(nt))
      allocate (dth_i(nt))
      allocate (sth(nt))
      allocate (sth_i(nt))
      allocate (cth(nt))
      allocate (dph(np))
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
        dt(i) = th(i+1) - th(i) !   half*(t(i+1) + t(i)) - half*(t(i) + t(i-1))
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
!$acc                    dth,dth_i,sth,sth_i,cth,dph,dph_i)
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
      use fields
      use input_parameters
      use constants
      use read_2d_file_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: j,k,ierr
      real(r_typ) :: fn1,fs1
!
      integer :: nft,nfp
      real(r_typ), dimension(:), allocatable :: tf,pf
      real(r_typ), dimension(:,:), allocatable :: vf,diffusion_coef_file
!
!-----------------------------------------------------------------------
!
! ****** Allocate memory for the total diffusion coef on the half mesh.
!
      allocate (diffusion_coef(nt,np))
!
! ****** Set the initial diffusion coef to the uniform value.
!
      diffusion_coef(:,:) = diffusion_coef_constant
!
! ****** Read the diffusion coef file, if it was specified.
!
      if (diffusion_coef_filename.ne.' ') then
!
! ****** Load the diffusion coef into an array (assumes PT).
!
        call read_2d_file (diffusion_coef_filename,nfp,nft,vf,pf,tf,ierr)
        vf(:,:) = TRANSPOSE(vf(:,:))
!
! ****** Interpolate the diffusion coef onto the half mesh (th,ph).
!
        allocate (diffusion_coef_file(nt,np))
!
        call interp2d (nft,nfp,tf,pf,vf,nt,np,th,ph, &
                       diffusion_coef_file,ierr)
!
! ****** Error exit: interpolation error.
!
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in LOAD_DIFFUSION_COEF:'
          write (*,*) '### The scales in the diffusion coef file'// &
                 ' are not monotonically increasing.'
          write (*,*) 'File name: ',trim(diffusion_coef_filename)
          call exit (1)
        end if
!
! ****** Enforce periodicity.
!
        diffusion_coef_file(:, 1) = diffusion_coef_file(:,npm1)
        diffusion_coef_file(:,np) = diffusion_coef_file(:,   2)
!
! ****** Set the pole value to only have an m=0 component.
!
        fn1 = 0.
        fs1 = 0.
        do k=2,npm1
          fn1 = fn1 + diffusion_coef_file(2,k)*dph(k)
          fs1 = fs1 + diffusion_coef_file(ntm1,k)*dph(k)
        enddo
        fn1 = fn1*twopi_i
        fs1 = fs1*twopi_i
!
        diffusion_coef_file( 1,:) = two*fn1 &
                                    - diffusion_coef_file(   2,:)
        diffusion_coef_file(nt,:) = two*fs1 &
                                    - diffusion_coef_file(ntm1,:)
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
        diffusion_coef(:,:) = diffusion_coef(:,:) &
                            + diffusion_coef_file(:,:)
!
        deallocate (diffusion_coef_file)
        deallocate (vf)
        deallocate (tf)
        deallocate (pf)
!
      end if
!
! ****** Add grid-based diffusion coef if requested.
!
      if (diffusion_coef_grid) then
        do k=1,np
          do j=1,nt
            diffusion_coef(j,k) = diffusion_coef(j,k) &
                             + (dth(j)**2 + (dph(k)*sth(j))**2)
          enddo
        enddo

        if (verbose) then
          write (*,*)
          write (*,*) 'Grid-based diffusion coef is activated.'
        end if

      end if
!
      diffusion_coef(:,:) = diffusion_coef_factor*diffusion_coef(:,:)
!$acc enter data copyin(diffusion_coef)
!
      call load_diffusion_matrix
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
      use fields
      use input_parameters
      use constants
      use read_2d_file_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: k,ierr
      real(r_typ) :: fn1,fs1
!
      integer :: nft,nfp
      real(r_typ), dimension(:), allocatable :: tf,pf
      real(r_typ), dimension(:,:), allocatable :: sf
!
!-----------------------------------------------------------------------
!
! ****** Allocate memory for the source term.
!
      allocate (source(ntm,npm))
!
      source=0.
!
! ****** Read the source file if it was specified.
!
      if (source_filename.ne.' ') then
!
        call read_2d_file (source_filename,nfp,nft,sf,pf,tf,ierr)
        sf(:,:) = TRANSPOSE(sf(:,:))
!
! ****** Interpolate the source term onto the main mesh (t,p).
!
        call interp2d (nft,nfp,tf,pf,sf,ntm,npm,t,p,source,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in LOAD_SOURCE:'
          write (*,*) '### The scales in the source term file are'// &
                        ' not monotonically increasing.'
          write (*,*) 'File name: ',trim(source_filename)
          call exit (1)
        end if
!
! ****** Enforce periodicity.
!
        source(:,   1)=half*(source(:,1)+source(:,npm))
        source(:,npm-1)=source(:,1)
!
! ****** Set the pole value to only have an m=0 component.
!
        fn1=0.
        fs1=0.
        do k=1,npm-2
          fn1=fn1+source(  1,k)*dp(k)
          fs1=fs1+source(ntm,k)*dp(k)
        enddo
        fn1=fn1*twopi_i
        fs1=fs1*twopi_i
!
        source(  1,:)=fn1
        source(ntm,:)=fs1
!
        if (verbose) then
          write (*,*)
          write (*,*) 'A source term was read in from file: ', &
                     trim(source_filename)
          write (*,*) 'Minimum value = ',minval(source)
          write (*,*) 'Maximum value = ',maxval(source)
        end if
!
        deallocate (sf)
        deallocate (tf)
        deallocate (pf)
!
      end if
!
end subroutine
!#######################################################################
subroutine load_vt
!
!-----------------------------------------------------------------------
!
! ****** Define the theta component of the flow.
!
! ****** This is read in from the file VT_FILE if it is not blank.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields
      use input_parameters
      use constants
      use read_2d_file_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: k,ierr
      real(r_typ) :: fn1,fs1
!
      integer :: nft,nfp
      real(r_typ), dimension(:), allocatable :: tf,pf
      real(r_typ), dimension(:,:), allocatable :: sf
!
!-----------------------------------------------------------------------
!
! ****** Read the vt file if it was specified.
!
        call read_2d_file (flow_vt_filename,nfp,nft,sf,pf,tf,ierr)
        sf(:,:) = TRANSPOSE(sf(:,:))
!
! ****** Interpolate vt onto the half mesh (th,p).
!
        call interp2d (nft,nfp,tf,pf,sf,nt,npm,th,p,vt,ierr)
! ****** Error exit: interpolation error.
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in LOAD_VT:'
          write (*,*) '### The scales in the vt file are'// &
                        ' not monotonically increasing.'
          write (*,*) 'File name: ',trim(flow_vt_filename)
          call exit (1)
        end if
!
! ****** Enforce periodicity.
!
        vt(:,  1)=vt(:,npm-1)
        vt(:,npm)=vt(:,    2)
!
! ****** Set the pole value to only have an m=0 component.
!
        fn1=0.
        fs1=0.
        do k=2,npm-1
          fn1=fn1+vt(   2,k)*dp(k)
          fs1=fs1+vt(ntm1,k)*dp(k)
        enddo
        fn1=fn1*twopi_i
        fs1=fs1*twopi_i
!
        vt( 1,:)=two*fn1-vt(   2,:)
        vt(nt,:)=two*fs1-vt(ntm1,:)
!
        if (verbose) then
          write (*,*)
          write (*,*) 'The flow component vt was read in'// &
                     ' from file: ',trim(flow_vt_filename)
          write (*,*) 'Minimum value = ',minval(vt)
          write (*,*) 'Maximum value = ',maxval(vt)
        end if
!
        deallocate (sf)
        deallocate (tf)
        deallocate (pf)
!
end subroutine
!#######################################################################
subroutine load_vp
!
!-----------------------------------------------------------------------
!
! ****** Define the phi component of the flow.
!
! ****** This is read in from the file VP_FILE if it is not blank.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields
      use input_parameters
      use constants
      use read_2d_file_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: k,ierr
      real(r_typ) :: fn1,fs1
!
      integer :: nft,nfp
      real(r_typ), dimension(:), allocatable :: tf,pf
      real(r_typ), dimension(:,:), allocatable :: sf
!
!-----------------------------------------------------------------------
!
! ****** Read the vp file if it was specified.
!
        call read_2d_file (flow_vp_filename,nfp,nft,sf,pf,tf,ierr)
        sf(:,:) = TRANSPOSE(sf(:,:))
!
! ****** Interpolate vp onto the correct mesh (t,ph).
!
        call interp2d (nft,nfp,tf,pf,sf,ntm,np,t,ph,vp,ierr)
!
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in LOAD_VP:'
          write (*,*) '### The scales in the vp file are'// &
                        ' not monotonically increasing.'
          write (*,*) 'File name: ',trim(flow_vp_filename)
          call exit (1)
        end if
!
! ****** Enforce periodicity.
!
        vp(:, 1)=vp(:,npm1)
        vp(:,np)=vp(:,   2)
!
! ****** Set the pole value to only have an m=0 component.
!
        fn1=0.
        fs1=0.
        do k=2,npm1
          fn1=fn1+vp(    2,k)*dph(k)
          fs1=fs1+vp(ntm-1,k)*dph(k)
        enddo
        fn1=fn1*twopi_i
        fs1=fs1*twopi_i
!
        vp(  1,:)=two*fn1-vp(     2,:)
        vp(ntm,:)=two*fs1-vp(ntm- 1,:)
!
        if (verbose) then
          write (*,*)
          write (*,*) 'The flow component vp was read in'// &
                     ' from file: ',trim(flow_vp_filename)
          write (*,*) 'Minimum value = ',minval(vp)
          write (*,*) 'Maximum value = ',maxval(vp)
        end if
!
        deallocate (sf)
        deallocate (tf)
        deallocate (pf)
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
      real(r_typ), parameter :: t0 = 46.0_r_typ
      real(r_typ), parameter :: t2 = -262.0_r_typ
      real(r_typ), parameter :: t4 = -379.0_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: j,k
!
!-----------------------------------------------------------------------
!
      do concurrent (k=1:np,j=1:ntm)
        vp(j,k) = vp(j,k) + &
                  m_s_to_rs_hr*st(j)*(t0 + t2*ct(j)**2 + t4*ct(j)**4)
      enddo
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
      real(r_typ), parameter :: s1 = 22.0_r_typ
      real(r_typ), parameter :: s3 = 11.0_r_typ
      real(r_typ), parameter :: s5 = -28.0_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: j,k
!
!-----------------------------------------------------------------------
!
      do concurrent (k=1:npm,j=1:nt)
        vt(j,k) = vt(j,k) -            &
                  m_s_to_rs_hr*sth(j)* &
                  (s1*cth(j) + s3*cth(j)**3 + s5*cth(j)**5)
      enddo
!
end subroutine
!#######################################################################
subroutine get_flow_dtmax (dtmaxflow)
!
!-----------------------------------------------------------------------
!
! ****** Get the maximum time step tha\t can be used for stable
! ****** (explicit) advection.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mesh
      use fields, ONLY : vt,vp
      use constants
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
      integer :: j,k
      real(r_typ) :: kdotv,deltat,dtmax
!
!-----------------------------------------------------------------------
!
      dtmax = huge(one)
!
!$omp parallel do default(shared) collapse(2) reduction(min:dtmax) &
!$omp private(kdotv,deltat)
!$acc parallel loop default(present) collapse(2) reduction(min:dtmax) &
!$acc private(kdotv,deltat)
      do k=2,npm-1
        do j=2,ntm-1
          kdotv = MAX(ABS(vt(j,k)),ABS(vt(j+1,k)))*dt_i(j) &
                + MAX(ABS(vp(j,k)),ABS(vp(j,k+1)))*st_i(j)*dp_i(k)
          deltat = one/MAX(kdotv,small_value)
          dtmax = MIN(dtmax,deltat)
        enddo
      enddo
!$omp end parallel do
!
      dtmaxflow = half*safety*dtmax
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
      integer :: j,k
      integer*8 :: i
      real(r_typ) :: dtime_euler_used
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary field.
!
      allocate (fold(ntm,npm))
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
        do concurrent (k=1:npm,j=1:ntm)
          f(j,k) = f(j,k) + dtime_euler_used*fold(j,k)
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
      integer :: j,k
      real(r_typ) :: cct,ccp,fn,fs
      real(r_typ), dimension(:,:), allocatable :: flux_t,flux_p
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary flux arrays on the half-mesh.
!
      allocate (flux_t(nt,npm))
      allocate (flux_p(ntm,np))
!$acc enter data create(flux_t,flux_p)
!
      do concurrent (k=1:npm,j=1:nt)
        flux_t(j,k) = zero
      enddo
!
      do concurrent (k=1:np,j=1:ntm)
        flux_p(j,k) = zero
      enddo      
!
! ****** Compute the fluxes at the cell faces.
!
      do concurrent (k=2:npm-1,j=2:ntm1)
        cct=sign(upwind,vt(j,k))
        flux_t(j,k)=vt(j,k)*half*((one-cct)*f(j,k)+(one+cct)*f(j-1,k))
      enddo
!
      do concurrent (k=2:npm1,j=2:ntm-1)
        ccp=sign(upwind,vp(j,k))
        flux_p(j,k)=vp(j,k)*half*((one-ccp)*f(j,k)+(one+ccp)*f(j,k-1))
      enddo
!
! ****** Set periodicity of the flux (seam).
!
      call set_periodic_bc_2d (flux_t,nt,npm)
      call set_periodic_bc_2d (flux_p,ntm,np)
!
! ****** Advect F by one time step.
!
      do concurrent (k=2:npm-1,j=2:ntm-1)
        f(j,k) = f(j,k) - dtime_local*(  (  sth(j+1)*flux_t(j+1,k) &
                                          - sth(j  )*flux_t(j  ,k) &
                                         )*st_i(j)*dt_i(j)         &
                                       + (  flux_p(j,k+1)          &
                                          - flux_p(j,k  )          &
                                         )*st_i(j)*dp_i(k)         &
                                      )
      enddo
!
! ****** Advect the values at the poles.
!
      fn = zero
      fs = zero
!
!$omp parallel do default(shared) reduction(+:fn,fs)
!$acc parallel loop default(present) reduction(+:fn,fs)
      do k=2,npm-1
        fn = fn + flux_t(   2,k)*dp(k)
        fs = fs + flux_t(ntm1,k)*dp(k)
      enddo
!$omp end parallel do
! ****** Note that the south pole needs a sign change since the
! ****** theta flux direction is reversed.
      do concurrent (k=2:npm-1)
        f(  1,k) = f(  1,k) - fn*dtime_local*two*pi_i*dt_i(  1)
        f(ntm,k) = f(ntm,k) + fs*dtime_local*two*pi_i*dt_i(ntm)
      end do
!
! ****** Set periodic boundary condition.
!
      call set_periodic_bc_2d (f,ntm,npm)
!
!$acc exit data delete(flux_t,flux_p)
      deallocate (flux_t)
      deallocate (flux_p)
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
      integer :: j,k
      real(r_typ), dimension(:,:), allocatable :: f1,f2,aop
!
!-----------------------------------------------------------------------
!
! ****** Allocate temporary arrays for RK.
!
      allocate (f1(ntm,npm))
      allocate (f2(ntm,npm))
      allocate (aop(ntm,npm))
!$acc enter data create(f1,f2,aop)
!
      call advection_operator_upwind (f,aop)
      do concurrent (k=1:npm,j=1:ntm)
        f1(j,k) = f(j,k) - dtime_local*aop(j,k)
      enddo
!
      call advection_operator_upwind (f1,aop)
      do concurrent (k=1:npm,j=1:ntm)
        f2(j,k) = three_quarter*f(j,k) + quarter*f1(j,k) &
                                       - quarter*dtime_local*aop(j,k)
      enddo
!
      call advection_operator_upwind (f2,aop)
      do concurrent (k=1:npm,j=1:ntm)
        f(j,k) = third*f(j,k) + two_third*f2(j,k) &
                              - two_third*dtime_local*aop(j,k)
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
      real(r_typ), dimension(ntm,npm), INTENT(IN) :: ftemp
      real(r_typ), dimension(ntm,npm), INTENT(OUT) :: aop
!
!-----------------------------------------------------------------------
!
      integer :: j,k
      real(r_typ) :: cct,ccp,fn,fs
      real(r_typ), dimension(:,:), allocatable :: flux_t,flux_p
!
!-----------------------------------------------------------------------
!
      allocate (flux_t(nt,npm))
      allocate (flux_p(ntm,np))
!$acc enter data create(flux_t,flux_p)
!
      do concurrent (k=1:npm,j=1:nt)
        flux_t(j,k) = zero
      enddo
!
      do concurrent (k=1:np,j=1:ntm)
        flux_p(j,k) = zero
      enddo    
!
! ****** Compute the fluxes at the cell faces.
!
      do concurrent (k=2:npm-1,j=2:ntm1)
        cct=sign(upwind,vt(j,k))
        flux_t(j,k)=vt(j,k)*half*((one-cct)*ftemp(j,k)+(one+cct)*ftemp(j-1,k))
      enddo
!
      do concurrent (k=2:npm1,j=2:ntm-1)
        ccp=sign(upwind,vp(j,k))
        flux_p(j,k)=vp(j,k)*half*((one-ccp)*ftemp(j,k)+(one+ccp)*ftemp(j,k-1))
      enddo
!
! ****** Set periodicity of the flux (seam).
!
      call set_periodic_bc_2d (flux_t,nt,npm)
      call set_periodic_bc_2d (flux_p,ntm,np)
!
! ****** Compute advection operator F.
!
      do concurrent (k=2:npm-1,j=2:ntm-1)
        aop(j,k) = (  (  sth(j+1)*flux_t(j+1,k)  &
                       - sth(j  )*flux_t(j  ,k)  &
                      )*st_i(j)*dt_i(j)          &
                    + (           flux_p(j,k+1)  &
                       -          flux_p(j,k  )  &
                      )*st_i(j)*dp_i(k)          &
                   )
      enddo
!
! ****** Get the advection operator at the poles.
!
      fn = zero
      fs = zero
!
!$omp parallel do default(shared) reduction(+:fn,fs)
!$acc parallel loop default(present) reduction(+:fn,fs)
        do k=2,npm-1
          fn = fn + flux_t(   2,k)*dp(k)
          fs = fs + flux_t(ntm1,k)*dp(k)
        enddo
!$omp end parallel do
! ****** Note that the south pole needs a sign change since the
! ****** theta flux direction is reversed.
        do concurrent (k=2:npm-1)
          aop(  1,k) = fn*two*pi_i*dt_i(  1)
          aop(ntm,k) = -fs*two*pi_i*dt_i(ntm)
        end do
!
! ****** Set periodic boundary condition.
!
      call set_periodic_bc_2d (aop,ntm,npm)
!
!$acc exit data delete(flux_t,flux_p)
      deallocate (flux_t)
      deallocate (flux_p)
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
      integer :: j,k
      real(r_typ) :: fn2_fn1,fs2_fs1
      real(r_typ), dimension(ntm,npm), INTENT(IN) :: x
      real(r_typ), dimension(ntm,npm), INTENT(OUT) :: y
!
!-----------------------------------------------------------------------
!
! ****** Compute y=Ax.
!
! ****** Compute inner points.
!
      do concurrent (k=2:npm-1,j=2:ntm-1)
        y(j,k) =  coef(j,k,1)*x(j,  k-1) &
                + coef(j,k,2)*x(j-1,k  ) &
                + coef(j,k,3)*x(j  ,k  ) &
                + coef(j,k,4)*x(j+1,k  ) &
                + coef(j,k,5)*x(j,  k+1)
      enddo
!
! ****** Compute boundary points.
!
! ****** Get the m=0 components near the poles.
!
      fn2_fn1 = zero
      fs2_fs1 = zero
!
!$omp parallel do default(shared) reduction(+:fn2_fn1,fs2_fs1)
!$acc parallel loop default(present) reduction(+:fn2_fn1,fs2_fs1)
      do k=2,npm-1
        fn2_fn1 = fn2_fn1 + (  diffusion_coef(1    ,k)     &
                             + diffusion_coef(2    ,k))    &
                           *(x(2    ,k) - x(1  ,k))*dp(k)
        fs2_fs1 = fs2_fs1 + (  diffusion_coef(nt-1,k)     &
                             + diffusion_coef(nt  ,k))    &
                           *(x(ntm-1,k) - x(ntm,k))*dp(k)
      enddo
!$omp end parallel do
!
      do concurrent (k=1:npm)
        y(  1,k) = fn2_fn1*dt_i(  1)*dt_i(  1)*pi_i
        y(ntm,k) = fs2_fs1*dt_i(ntm)*dt_i(ntm)*pi_i
      enddo
!
! ****** Set the periodic boundary conditions.
!
      call set_periodic_bc_2d (y,ntm,npm)
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
      real(r_typ), dimension(ntm,npm) :: f
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: fn1,fn2,fs1,fs2
!
!-----------------------------------------------------------------------
!
      integer :: k
!
!-----------------------------------------------------------------------
!
      fn1 = 0.
      fn2 = 0.
      fs1 = 0.
      fs2 = 0.
!
      do k=2,npm-1
        fn1 = fn1 + f(1,k)*dp(k)
        fn2 = fn2 + f(2,k)*dp(k)
        fs1 = fs1 + f(ntm,k)*dp(k)
        fs2 = fs2 + f(ntm-1,k)*dp(k)
      enddo
!
      fn1 = fn1*twopi_i
      fn2 = fn2*twopi_i
      fs1 = fs1*twopi_i
      fs2 = fs2*twopi_i
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
subroutine read_input
!
!-----------------------------------------------------------------------
!
! ****** Read parameters from the command line arguments.
!
!-----------------------------------------------------------------------
!
      use ident
      use number_types
      use constants
      use syntax
      use paragraph_def
      use get_usage_line_interface
      use print_par_interface
      use delete_par_interface
      use input_parameters
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Storage the for usage line.
!
      type(paragraph), pointer :: usage
!
! ****** Storage for the error message.
!
      character(72) :: errmsg
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(512) :: arg
      logical :: set
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fpval
      integer, external :: intval
!
!-----------------------------------------------------------------------
!
! ****** Define the syntax.
!
      call defarg (GROUP_K ,'-v',' ',' ')
      call defarg (GROUP_A ,'initial_map_filename',' ',' ')
      call defarg (GROUP_A ,'output_map_root_filename',' ',' ')
      call defarg (GROUP_KA,'-output_map_directory','output_maps','<dir>')
      call defarg (GROUP_KA,'-visc','0.','<val>')
      call defarg (GROUP_KA,'-diffusion_coef_filename','<none>','<file>')
      call defarg (GROUP_K,'-diffusion_coef_grid',' ',' ')
      call defarg (GROUP_KA,'-diffusion_coef_factor','1.','<val>')
      call defarg (GROUP_KA,'-s','<none>','<file>')
      call defarg (GROUP_KA,'-time','1.0','<val>')
      call defarg (GROUP_KA,'-dtmax','1.0e16','<val>')
      call defarg (GROUP_KA,'-dm','3','<val>')
      call defarg (GROUP_KA,'-diffusion_subcycles','30',' ')
      call defarg (GROUP_KA,'-omci','0','<val>')
      call defarg (GROUP_KA,'-omct','0.0','<val>')
      call defarg (GROUP_KA,'-uw','1.','<val>')
      call defarg (GROUP_KA,'-vtfile','<none>','<file>')
      call defarg (GROUP_KA,'-vpfile','<none>','<file>')
      call defarg (GROUP_KA,'-vpomega','0.','<val>')
      call defarg (GROUP_K,'-va',' ',' ')
      call defarg (GROUP_K,'-nostrang',' ',' ')
      call defarg (GROUP_KA,'-dr','0','<val>')
      call defarg (GROUP_KA,'-mf','0','<val>')
      call defarg (GROUP_K,'-diff',' ',' ')
      call defarg (GROUP_K,'-flow',' ',' ')
      call defarg (GROUP_KA,'-fm','2','<val>')
      call defarg (GROUP_KA,'-difftfac','1.','<val>')
      call defarg (GROUP_KA,'-diffpfac','1.','<val>')
      call defarg (GROUP_K,'-vrun',' ',' ')
      call defarg (GROUP_K,'-assimilate_data',' ',' ')
      call defarg (GROUP_KA ,'-assimilate_data_map_list_filename','<none>','<file>')
      call defarg (GROUP_KA ,'-assimilate_data_map_root_dir','.','<dir>')
!
! ****** Parse the command line.
!
      call parse (errmsg,ierr)
!
      if (ierr.ne.0) then
!
        write (*,*)
        write (*,*) '### ',cname,' Version ',cvers,' of ',cdate,'.'
!
! ****** Print the usage line.
!
        call get_usage_line (usage)
!
        write (*,*)
        write (*,*) 'Usage:'
        write (*,*)
!
        call print_par (usage)
!
        call delete_par (usage)
!
        write (*,*)
        write (*,*) '-time <#> Set the end time. (Default 1.0)'
        write (*,*)
       write (*,*) 'Viscosity is set as: diffusion_coef_factor*(visc + visc_file '// &
                    '+ visc_grid)'
        write (*,*) ' '
        write (*,*) ' -diffusion_coef_factor <#> (Overall multiplier, Default: 1.0)'
        write (*,*) ' -visc <#> (Uniform viscosity)'
        write (*,*) '           Default: If neither diffusion_coef_filename or '// &
                        'diffusion_coef_grid set, 1, otherwise, 0.'
        write (*,*) ' -diffusion_coef_filename <fname> (Read from file, Default: none)'
        write (*,*) '                    The file must contain a 2D'// &
                                      ' field with scales.'
        write (*,*) '                    Any region outside the '// &
                                        'domain is set to 0.'
        write (*,*) ' -diffusion_coef_grid (Grid viscosity, Default: disabled)'
        write (*,*) '           Defined by: visc_grid = dt^2 + '// &
                                '(sin(t)*dt)^2'
        write (*,*)
        write (*,*) 'Source term can be specified using: -s <file>.'
        write (*,*) '            (Default: none)'
        write (*,*) '            The file must contain a 2D'// &
                                ' field with scales.'
        write (*,*) '            Any region outside the '// &
                                'domain is set to 0.'
        write (*,*) ' '
        write (*,*)
        write (*,*) 'The flow field is specified using'// &
                    ' -vtfile and -vpfile.'// &
                    '  Flow components that are not'
        write (*,*) 'specified have zero values.'
        write (*,*)
        write (*,*) 'The upwind coefficient for the advective'// &
                    ' advance can be set via -uw'
        write (*,*) '(default=1.).'
!
        write (*,*) '-exp (Use old explicit algorithm for visc).'
        write (*,*) '-vrun (Activate validation run, more details to follow)'
        write (*,*) '-va (attenuate the velocity)'
        call exit (1)
!
      end if
!
! ****** Set the parameters.
!
! ****** Verbose flag.
!
      call fetcharg ('-v',set,arg)
      verbose = set
!
! ****** Input.
!
      call fetcharg ('initial_map_filename',set,arg)
      initial_map_filename = trim(arg)
!
! ****** Output.
!
      call fetcharg ('output_map_root_filename',set,arg)
      if (set) then
        output_map_root_filename = trim(arg)
      else
        output_map_root_filename = '.'
      end if      
!
      call fetcharg ('-output_map_directory',set,arg)
      if (set) then
        output_map_directory = trim(arg)
      else
        output_map_directory = 'output_maps'
      end if
!
      call fetcharg ('-omci',set,arg)
      output_map_idx_cadence = intval(arg,'-omci')
!
      call fetcharg ('-omct',set,arg)
      output_map_time_cadence = fpval(arg,'-omct')
!
      if (output_map_time_cadence.gt.0 .and. output_map_idx_cadence.gt.0) then
        write(*,*) "ERROR!  Cannot use both output_map_time_cadence"
        write(*,*) "        and output_map_idx_cadence at the same time!"
        call wrap_it_up
        STOP 1
      end if
!
! ****** Viscosity file.
!
      call fetcharg ('-diffusion_coef_filename',set,arg)
      if (set) then
        diffusion_coef_filename = trim(arg)
      else
        diffusion_coef_filename = ' '
      end if
!
! ****** Source term file.
!
      call fetcharg ('-s',set,arg)
      if (set) then
        source_filename = trim(arg)
      else
        source_filename = ' '
      end if
!
! ****** Uniform viscosity.
!
      call fetcharg ('-visc',set,arg)
      diffusion_coef_constant = fpval(arg,'-visc')
!
! ****** Viscosity multiplier.
!
      call fetcharg ('-diffusion_coef_factor',set,arg)
      if (set) diffusion_coef_factor = fpval(arg,'-diffusion_coef_factor')
!
! ****** Activate auto-viscosity.
!
      call fetcharg ('-diffusion_coef_grid',set,arg)
      diffusion_coef_grid = set
!
! ****** Upwind coefficient.
!
      call fetcharg ('-uw',set,arg)
      upwind=fpval(arg,'-uw')
!
! ****** File containing the theta component of the flow.
!
      call fetcharg ('-vtfile',set,arg)
      if (set) then
        flow_vt_filename = trim(arg)
      else
        flow_vt_filename = ' '
      end if
!
! ****** File containing the phi component of the flow.
!
      call fetcharg ('-vpfile',set,arg)
      if (set) then
        flow_vp_filename = trim(arg)
      else
        flow_vp_filename = ' '
      end if
!
      call fetcharg ('-vpomega',set,arg)
      flow_vp_rigid_omega = fpval(arg,'-vpomega')
!
      call fetcharg ('-va',set,arg)
      if (set) flow_attenuate = .true.
!
      call fetcharg ('-dr',set,arg)
      flow_dr_model = intval(arg,'-dr')
!
      call fetcharg ('-mf',set,arg)
      flow_mf_model = intval(arg,'-mf')
!
! ****** Time to diffuse for.
!
      call fetcharg ('-time',set,arg)
      time_end=fpval(arg,'-time')
!
! ****** Maximum allowed timestep.
!
      call fetcharg ('-dtmax',set,arg)
      dt_max=fpval(arg,'-dtmax')
!
! ****** Use original explicit algorithm.
!
      call fetcharg ('-dm',set,arg)
      if (set) diffusion_num_method = intval(arg,'-dm')
!
! ****** Set number of STS sub-cycles per flow step.
!
      call fetcharg ('-diffusion_subcycles',set,arg)
      diffusion_subcycles = intval(arg,'-diffusion_subcycles')
!
      call fetcharg ('-nostrang',set,arg)
      if (set) strang_splitting = .false.
!
      call fetcharg ('-diff',set,arg)
      if (set) advance_diffusion = .true.
!
      call fetcharg ('-flow',set,arg)
      if (set) advance_flow = .true.
!
      call fetcharg ('-fm',set,arg)
      if (set) flow_num_method = intval(arg,'-fm')
!
      call fetcharg ('-vrun',set,arg)
      validation_run = set
!
! ****** Data assimilation.
!
      call fetcharg ('-assimilate_data',set,arg)
      assimilate_data = set
!
      call fetcharg ('-assimilate_data_map_list_filename',set,arg)
      if (set) then
        assimilate_data_map_list_filename = trim(arg)
      else
        assimilate_data_map_list_filename = ' '
      end if
!
      call fetcharg ('-assimilate_data_map_root_dir',set,arg)
      if (set) then
        assimilate_data_map_root_dir = trim(arg)
      else
        assimilate_data_map_root_dir = '.'
      end if
!
end subroutine
!#######################################################################
!
!-----------------------------------------------------------------------
!
! ****** Update log:
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
!   - Updated default of dt_max_increase_fac to 0.0 (disabled).
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
!   - Made the default splitting method to be Strang splitting.
!     The "-strang" input flag is now "-nostrang" and disables
!     strang splitting.
!   - Added "-fm" input flag to set advection numerical method.
!     1) FE+UW  2) RK3TVD+UW  3) RK3TVD+WENO3 (coming soon).
!   - Removed "-flowpoleavg".  No use case.
!
! 04/29/2022, RC, Version 0.7.1:
!   - BUG FIX: Changed staggering of the velocity.  Now, Vt is on the 
!     theta half mesh and the phi main mesh, while Vp is on the 
!     theta main mesh and the phi half mesh.  This avoids averaging in 
!     the advection operator (which was not being done).
!   - BUG FIX: The diffusion_coef was not being averaged.
!   - BUG FIX: Fixed issue when not outputting time-dept I/O.
!-----------------------------------------------------------------------
!
!#######################################################################
