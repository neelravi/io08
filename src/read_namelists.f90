!
! copyright (c) 2002-2020 quantum espresso group
! this file is distributed under the terms of the
! gnu general public license. see the file `license'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------
! tb
! included gate related stuff, search 'tb'
!---------------------------------------------
!
!----------------------------------------------------------------------------
module read_namelists_module
  !----------------------------------------------------------------------------
  !
  !  ... this module handles the reading of input namelists
  !  ... written by carlo cavazzoni, with many additions
  !  --------------------------------------------------
  !
  use kinds,     only : dp
  use input_parameters
  !
  implicit none
  !
  save
  !
  private
  !
  real(dp), parameter :: sm_not_set = -20.0_dp
  !
  public :: read_namelists, sm_not_set
  public :: check_namelist_read ! made public upon request of a.jay
  ! fixme: should the following ones be public?
  public :: control_defaults, system_defaults, &
       electrons_defaults, wannier_ac_defaults, ions_defaults, &
       cell_defaults, press_ai_defaults, wannier_defaults, control_bcast,&
       system_bcast, electrons_bcast, ions_bcast, cell_bcast, &
       press_ai_bcast, wannier_bcast, wannier_ac_bcast, control_checkin, &
       system_checkin, electrons_checkin, ions_checkin, cell_checkin, &
       wannier_checkin, wannier_ac_checkin, fixval
  !
  !  ... end of module-scope declarations
  !
  !  ----------------------------------------------
  !
  contains
     !
     !=-----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist control
     !
     !=-----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine control_defaults( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       character(len=20) ::    temp_string 
       !
       !
       if ( prog == 'pw' ) then
          title = ' '
          calculation = 'scf'
       else
          title = 'md simulation'
          calculation = 'cp'
       end if

       verbosity = 'default'
       if( prog == 'pw' ) restart_mode = 'from_scratch'
       if( prog == 'cp' ) restart_mode = 'restart'
       nstep  = 50
       if( prog == 'pw' ) iprint = 100000
       if( prog == 'cp' ) iprint = 10
       if( prog == 'pw' ) isave = 0
       if( prog == 'cp' ) isave = 100
       !
       tstress = .false.
       tprnfor = .false.
       tabps = .false.
       !
       if( prog == 'pw' ) dt  = 20.0_dp
       if( prog == 'cp' ) dt  =  1.0_dp
       !
       ndr = 50
       ndw = 50
       !
       ! ... use the path specified as outdir and the filename prefix
       ! ... to store output data
       !
       call get_environment_variable( 'espresso_tmpdir', outdir )
       if ( trim( outdir ) == ' ' ) outdir = './'
       if( prog == 'pw' ) prefix = 'pwscf'
       if( prog == 'cp' ) prefix = 'cp'
       !
       ! ... directory containing the pseudopotentials
       !
       call get_environment_variable( 'espresso_pseudo', pseudo_dir )
       if ( trim( pseudo_dir ) == ' ') then
          call get_environment_variable( 'home', pseudo_dir )
          pseudo_dir = trim( pseudo_dir ) // '/espresso/pseudo/'
       end if
       !
       ! ... max number of md steps added to the xml file. needs to be limited for very long 
       !     md simulations 
       call get_environment_variable('max_xml_steps', temp_string) 
            if ( trim(temp_string) .ne.  ' ')  read(temp_string, *) max_xml_steps 
       refg          = 0.05_dp
       max_seconds   = 1.e+7_dp
       ekin_conv_thr = 1.e-6_dp
       etot_conv_thr = 1.e-4_dp
       forc_conv_thr = 1.e-3_dp
       disk_io  = 'default'
       dipfield = .false.
       gate     = .false. !tb
       lberry   = .false.
       gdir     = 0
       nppstr   = 0
       wf_collect = .true.
       lelfield = .false.
       lorbm = .false.
       nberrycyc  = 1
       lecrpa   = .false.   
       tqmmm = .false.
       !
       saverho = .true.
       memory = 'default'
       !
       lfcpopt = .false.
       lfcpdyn = .false.
       !
       call get_environment_variable( 'qexml', input_xml_schema_file )
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist system
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine system_defaults( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       !
       !
       ibrav  = -1
       celldm = (/ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
       a = 0.0_dp
       b = 0.0_dp
       c = 0.0_dp
       cosab = 0.0_dp
       cosac = 0.0_dp
       cosbc = 0.0_dp
       nat    = 0
       ntyp   = 0
       nbnd   = 0
       tot_charge = 0.0_dp
       tot_magnetization = -1
       ecutwfc = 0.0_dp
       ecutrho = 0.0_dp
       nr1  = 0
       nr2  = 0
       nr3  = 0
       nr1s = 0
       nr2s = 0
       nr3s = 0
       nr1b = 0
       nr2b = 0
       nr3b = 0
       occupations = 'fixed'
       smearing = 'gaussian'
       degauss = 0.0_dp
       nspin = 1
       nosym = .false.
       nosym_evc = .false.
       force_symmorphic = .false.
       use_all_frac = .false.
       noinv = .false.
       ecfixed = 0.0_dp
       qcutz   = 0.0_dp
       q2sigma = 0.01_dp
       input_dft = 'none'
       ecutfock  = -1.0_dp
       starting_charge = 0.0_dp
!
! ... set starting_magnetization to an invalid value:
! ... in pw starting_magnetization must be set for at least one atomic type
! ... (unless the magnetization is set in other ways)
! ... in cp starting_magnetization must remain unset
!
       starting_magnetization = sm_not_set

       if ( prog == 'pw' ) then
          starting_ns_eigenvalue = -1.0_dp
          u_projection_type = 'atomic'
       end if
       !
       ! .. dft + u and its extensions
       !
       lda_plus_u = .false.
       lda_plus_u_kind = 0
       hubbard_u = 0.0_dp
       hubbard_u_back = 0.0_dp
       hubbard_v = 0.0_dp
       hubbard_j0 = 0.0_dp
       hubbard_j = 0.0_dp
       hubbard_alpha = 0.0_dp
       hubbard_alpha_back = 0.0_dp
       hubbard_beta = 0.0_dp
       hubbard_parameters = 'input'
       reserv = .false.
       reserv_back = .false.
       backall = .false.
       lback = -1
       l1back = -1
       hub_pot_fix = .false.
       step_pen=.false.
       a_pen=0.0_dp
       sigma_pen=0.01_dp
       alpha_pen=0.0_dp
       !
       ! ... exx
       !
       ace=.true.
       n_proj = 0    
       localization_thr = 0.0_dp
       scdm=.false.
       scdmden=1.0d0
       scdmgrd=1.0d0
       nscdm=1
       !
       ! ... electric fields
       !
       edir = 1
       emaxpos = 0.5_dp
       eopreg = 0.1_dp
       eamp = 0.0_dp
       ! tb gate related variables
       zgate = 0.5
       relaxz = .false.
       block = .false.
       block_1 = 0.45
       block_2 = 0.55
       block_height = 0.0
       !
       !  ... postprocessing of dos & phonons & el-ph
       !
       la2f = .false.
       !
       ! ... non collinear program variables
       !
       lspinorb = .false.
       lforcet = .false.
       starting_spin_angle=.false.
       noncolin = .false.
       lambda = 1.0_dp
       constrained_magnetization= 'none'
       fixed_magnetization = 0.0_dp
       b_field = 0.0_dp
       angle1 = 0.0_dp
       angle2 = 0.0_dp
       report =-1
       !
       no_t_rev = .false.
       !
       assume_isolated = 'none'
       !
       one_atom_occupations=.false.
       !
       spline_ps = .false.
       ! 
       real_space = .false.
       !
       ! ... dft-d, tkatchenko-scheffler, xdm
       !
       vdw_corr    = 'none'
       london      = .false.
       london_s6   = 0.75_dp
       london_rcut = 200.00_dp
       london_c6   = -1.0_dp
       london_rvdw = -1.0_dp
       ts_vdw          = .false.
       ts_vdw_isolated = .false.
       ts_vdw_econv_thr = 1.e-6_dp
       xdm = .false.
       xdm_a1 = 0.0_dp
       xdm_a2 = 0.0_dp
       dftd3_version = 3
       dftd3_threebody = .true.
       !
       ! ... esm
       !
       esm_bc='pbc'
       esm_efield=0.0_dp
       esm_w=0.0_dp
       esm_a=0.0_dp
       esm_zb=-2.0_dp
       esm_nfit=4
       esm_debug=.false.
       esm_debug_gpmax=0
       !
       ! ... fcp
       !
       fcp_mu          = 0.0_dp
       fcp_mass        = 10000.0_dp
       fcp_tempw       = 0.0_dp
       fcp_relax       = 'lm'
       fcp_relax_step  = 0.5_dp
       fcp_relax_crit  = 0.001_dp
       fcp_mdiis_size  = 4
       fcp_mdiis_step  = 0.2_dp
       !
       ! ... wyckoff
       !
       space_group=0
       uniqueb = .false.
       origin_choice = 1
       rhombohedral = .true.
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist electrons
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine electrons_defaults( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       !
       !
       emass = 400.0_dp
       emass_cutoff = 2.5_dp
       orthogonalization = 'ortho'
       ortho_eps = 1.e-9_dp
       ortho_max = 300
       electron_maxstep = 100
       scf_must_converge = .true.
       !
       ! ... ( 'sd' | 'cg' | 'damp' | 'verlet' | 'none' | 'diis' | 'cp-bo' )
       !
       electron_dynamics = 'none'
       electron_damping = 0.1_dp
       !
       ! ... ( 'zero' | 'default' )
       !
       electron_velocities = 'default'
       !
       ! ... ( 'nose' | 'not_controlled' | 'rescaling')
       !
       electron_temperature = 'not_controlled'
       ekincw = 0.001_dp
       fnosee = 1.0_dp
       ampre  = 0.0_dp
       grease = 1.0_dp
       conv_thr = 1.e-6_dp
       diis_size = 4
       diis_nreset = 3
       diis_hcut = 1.0_dp
       diis_wthr = 0.0_dp
       diis_delt = 0.0_dp
       diis_maxstep = 100
       diis_rot = .false.
       diis_fthr = 0.0_dp
       diis_temp = 0.0_dp
       diis_achmix = 0.0_dp
       diis_g0chmix = 0.0_dp
       diis_g1chmix = 0.0_dp
       diis_nchmix = 3
       diis_nrot = 3
       diis_rothr  = 0.0_dp
       diis_ethr   = 0.0_dp
       diis_chguess = .false.
       mixing_mode = 'plain'
       mixing_fixed_ns = 0
       mixing_beta = 0.7_dp
       mixing_ndim = 8
       diagonalization = 'david'
       diago_thr_init = 0.0_dp
       diago_cg_maxiter = 20
       diago_ppcg_maxiter = 20
       diago_david_ndim = 2
       diago_full_acc = .false.
       !
       sic = 'none'
       sic_epsilon = 0.0_dp
       sic_alpha = 0.0_dp
       force_pairing = .false.
       !
       fermi_energy = 0.0_dp
       n_inner = 2
       niter_cold_restart=1
       lambda_cold=0.03_dp
       rotation_dynamics = "line-minimization"
       occupation_dynamics = "line-minimization"
       rotmass = 0.0_dp
       occmass = 0.0_dp
       rotation_damping = 0.0_dp
       occupation_damping = 0.0_dp
       !
       tcg     = .false.
       maxiter = 100
       passop  = 0.3_dp
       niter_cg_restart = 20
       etresh  = 1.e-6_dp
       !
       epol   = 3
       efield = 0.0_dp
       epol2  = 3
       efield2 = 0.0_dp
       efield_cart(1)=0.d0
       efield_cart(2)=0.d0
       efield_cart(3)=0.d0
       efield_phase='none'
       !
       occupation_constraints = .false.
       !
       adaptive_thr   =  .false.
       conv_thr_init  =  0.1e-2_dp
       conv_thr_multi =  0.1_dp
       !
       ! ... cp-bo ...
       tcpbo = .false.
       emass_emin = 200.0_dp
       emass_cutoff_emin = 6.0_dp
       electron_damping_emin = 0.35_dp
       dt_emin = 4.0_dp
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist wannier_ac
     !
     !----------------------------------------------------------------------
     subroutine wannier_ac_defaults( prog )
       !----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       !
       !
       plot_wannier = .false.
       use_energy_int = .false.
       print_wannier_coeff = .false.
       nwan = 0
       constrain_pot = 0.d0
       plot_wan_num = 0
       plot_wan_spin = 1
       !
       return
       !
     end subroutine

     !=----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist ions
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine ions_defaults( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       !
       ! ... ( 'sd' | 'cg' | 'damp' | 'verlet' | 'none' | 'bfgs' | 'beeman' )
       !
       ion_dynamics = 'none'
       ion_radius   = 0.5_dp
       ion_damping  = 0.1_dp
       !
       ! ... ( 'default' | 'from_input' )
       !
       ion_positions = 'default'
       !
       ! ... ( 'zero' | 'default' | 'from_input' )
       !
       ion_velocities = 'default'
       !
       ! ... ( 'nose' | 'not_controlled' | 'rescaling' | 'berendsen' |
       !       'andersen' | 'initial' )
       !
       ion_temperature = 'not_controlled'
       !
       tempw       = 300.0_dp
       fnosep      = -1.0_dp
       fnosep(1)   = 1.0_dp
       nhpcl       = 0
       nhptyp      = 0
       ndega       = 0
       tranp       = .false.
       amprp       = 0.0_dp
       greasp      = 1.0_dp
       tolp        = 100.0_dp
       ion_nstepe  = 1
       ion_maxstep = 100
       delta_t     = 1.0_dp
       nraise      = 1
       !
       refold_pos       = .false.
       remove_rigid_rot = .false.
       !
       upscale           = 100.0_dp
       pot_extrapolation = 'atomic'
       wfc_extrapolation = 'none'
       !
       ! ... bfgs defaults
       !
       bfgs_ndim        = 1
       trust_radius_max = 0.8_dp   ! bohr
       trust_radius_min = 1.e-4_dp ! bohr
       trust_radius_ini = 0.5_dp   ! bohr
       w_1              = 0.01_dp
       w_2              = 0.50_dp
       !
       l_mplathe=.false.
       n_muller=0
       np_muller=1
       l_exit_muller=.false.
       

       return
       !
     end subroutine
     !
     !
     !=----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist cell
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine cell_defaults( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       !
       !
       cell_parameters = 'default'
       !
       ! ... ( 'sd' | 'pr' | 'none' | 'w' | 'damp-pr' | 'damp-w' | 'bfgs' )
       !
       cell_dynamics = 'none'
       !
       ! ... ( 'zero' | 'default' )
       !
       cell_velocities = 'default'
       press = 0.0_dp
       wmass = 0.0_dp
       !
       ! ... ( 'nose' | 'not_controlled' | 'rescaling' )
       !
       cell_temperature = 'not_controlled'
       temph = 0.0_dp
       fnoseh = 1.0_dp
       greash = 1.0_dp
       !
       ! ... ('all'* | 'volume' | 'x' | 'y' | 'z' | 'xy' | 'xz' | 'yz' | 'xyz' )
       !
       cell_dofree = 'all'
       cell_factor = 0.0_dp
       cell_nstepe = 1
       cell_damping = 0.1_dp
       press_conv_thr = 0.5_dp
       treinit_gvecs = .false.
       !
       return
       !
     end subroutine
     !
     !
     !=----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist press_ai
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     subroutine press_ai_defaults( prog )
     !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       !
       abivol = .false.
       abisur = .false.
       pvar = .false.
       fill_vac = .false.
       cntr = .false.
       scale_at = .false.
       t_gauss = .false.
       jellium = .false.

       p_ext = 0.0_dp
       p_in = 0.0_dp
       p_fin = 0.0_dp
       surf_t = 0.0_dp
       rho_thr = 0.0_dp
       dthr = 0.0_dp
       step_rad = 0.0_dp
       delta_eps = 0.0_dp
       delta_sigma = 0.0_dp
       r_j = 0.0_dp
       h_j = 0.0_dp

       n_cntr = 0
       axis = 3
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  variables initialization for namelist wannier
     !
     !-----------------------------------------------------------------------
     subroutine wannier_defaults( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2) :: prog   ! ... specify the calling program
       !
       !
       wf_efield = .false.
       wf_switch = .false.
       !
       sw_len = 1
       !
       efx0 = 0.0_dp
       efy0 = 0.0_dp
       efz0 = 0.0_dp
       efx1 = 0.0_dp
       efy1 = 0.0_dp
       efz1 = 0.0_dp
       !
       wfsd = 1
       !
       wfdt        = 5.0_dp
       maxwfdt     = 0.30_dp
       wf_q        = 1500.0_dp
       wf_friction = 0.3_dp
!=======================================================================
!exx_wf related
       exx_neigh        =  60 
       vnbsp            =  0
       exx_poisson_eps  =  1.e-6_dp
       exx_dis_cutoff   =  8.0_dp
       exx_ps_rcut_self =  6.0_dp
       exx_ps_rcut_pair =  5.0_dp
       exx_me_rcut_self = 10.0_dp
       exx_me_rcut_pair =  7.0_dp
       exx_use_cube_domain = .false.
!=======================================================================
       !
       nit    = 10
       nsd    = 10
       nsteps = 20
       !
       tolw = 1.e-8_dp
       !
       adapt = .true.
       !
       calwf  = 3
       nwf    = 0
       wffort = 40
       !
       writev = .false.
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist control
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine control_bcast()
       !-----------------------------------------------------------------------
       !
       use io_global, only : ionode_id
       use mp,        only : mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       call mp_bcast( title,         ionode_id, intra_image_comm )
       call mp_bcast( calculation,   ionode_id, intra_image_comm )
       call mp_bcast( verbosity,     ionode_id, intra_image_comm )
       call mp_bcast( restart_mode,  ionode_id, intra_image_comm )
       call mp_bcast( nstep,         ionode_id, intra_image_comm )
       call mp_bcast( iprint,        ionode_id, intra_image_comm )
       call mp_bcast( isave,         ionode_id, intra_image_comm )
       call mp_bcast( tstress,       ionode_id, intra_image_comm )
       call mp_bcast( tprnfor,       ionode_id, intra_image_comm )
       call mp_bcast( tabps,         ionode_id, intra_image_comm )
       call mp_bcast( dt,            ionode_id, intra_image_comm )
       call mp_bcast( ndr,           ionode_id, intra_image_comm )
       call mp_bcast( ndw,           ionode_id, intra_image_comm )
       call mp_bcast( outdir,        ionode_id, intra_image_comm )
       call mp_bcast( wfcdir,        ionode_id, intra_image_comm )
       call mp_bcast( prefix,        ionode_id, intra_image_comm )
       call mp_bcast( max_seconds,   ionode_id, intra_image_comm )
       call mp_bcast( ekin_conv_thr, ionode_id, intra_image_comm )
       call mp_bcast( etot_conv_thr, ionode_id, intra_image_comm )
       call mp_bcast( forc_conv_thr, ionode_id, intra_image_comm )
       call mp_bcast( pseudo_dir,    ionode_id, intra_image_comm )
       call mp_bcast( refg,          ionode_id, intra_image_comm )
       call mp_bcast( disk_io,       ionode_id, intra_image_comm )
       call mp_bcast( tefield,       ionode_id, intra_image_comm )
       call mp_bcast( tefield2,      ionode_id, intra_image_comm )
       call mp_bcast( dipfield,      ionode_id, intra_image_comm )
       call mp_bcast( lberry,        ionode_id, intra_image_comm )
       call mp_bcast( gdir,          ionode_id, intra_image_comm )
       call mp_bcast( nppstr,        ionode_id, intra_image_comm )
       call mp_bcast( point_label_type,   ionode_id, intra_image_comm )
       call mp_bcast( wf_collect,    ionode_id, intra_image_comm )
       call mp_bcast( lelfield,      ionode_id, intra_image_comm )
       call mp_bcast( lorbm,         ionode_id, intra_image_comm )
       call mp_bcast( nberrycyc,     ionode_id, intra_image_comm )
       call mp_bcast( saverho,       ionode_id, intra_image_comm )
       call mp_bcast( lecrpa,        ionode_id, intra_image_comm )
       call mp_bcast( tqmmm,         ionode_id, intra_image_comm )
       call mp_bcast( vdw_table_name,ionode_id, intra_image_comm )
       call mp_bcast( memory,        ionode_id, intra_image_comm )
       call mp_bcast( lfcpopt,       ionode_id, intra_image_comm )
       call mp_bcast( lfcpdyn,       ionode_id, intra_image_comm )
       call mp_bcast( input_xml_schema_file, ionode_id, intra_image_comm )
       call mp_bcast( gate,          ionode_id, intra_image_comm ) !tb
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist system
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine system_bcast()
       !-----------------------------------------------------------------------
       !
       use io_global, only : ionode_id
       use mp,        only : mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       call mp_bcast( ibrav,             ionode_id, intra_image_comm )
       call mp_bcast( celldm,            ionode_id, intra_image_comm )
       call mp_bcast( a,                 ionode_id, intra_image_comm )
       call mp_bcast( b,                 ionode_id, intra_image_comm )
       call mp_bcast( c,                 ionode_id, intra_image_comm )
       call mp_bcast( cosab,             ionode_id, intra_image_comm )
       call mp_bcast( cosac,             ionode_id, intra_image_comm )
       call mp_bcast( cosbc,             ionode_id, intra_image_comm )
       call mp_bcast( nat,               ionode_id, intra_image_comm )
       call mp_bcast( ntyp,              ionode_id, intra_image_comm )
       call mp_bcast( nbnd,              ionode_id, intra_image_comm )
       call mp_bcast( tot_charge,        ionode_id, intra_image_comm )
       call mp_bcast( tot_magnetization, ionode_id, intra_image_comm )
       call mp_bcast( ecutwfc,           ionode_id, intra_image_comm )
       call mp_bcast( ecutrho,           ionode_id, intra_image_comm )
       call mp_bcast( nr1,               ionode_id, intra_image_comm )
       call mp_bcast( nr2,               ionode_id, intra_image_comm )
       call mp_bcast( nr3,               ionode_id, intra_image_comm )
       call mp_bcast( nr1s,              ionode_id, intra_image_comm )
       call mp_bcast( nr2s,              ionode_id, intra_image_comm )
       call mp_bcast( nr3s,              ionode_id, intra_image_comm )
       call mp_bcast( nr1b,              ionode_id, intra_image_comm )
       call mp_bcast( nr2b,              ionode_id, intra_image_comm )
       call mp_bcast( nr3b,              ionode_id, intra_image_comm )
       call mp_bcast( occupations,       ionode_id, intra_image_comm )
       call mp_bcast( smearing,          ionode_id, intra_image_comm )
       call mp_bcast( degauss,           ionode_id, intra_image_comm )
       call mp_bcast( nspin,             ionode_id, intra_image_comm )
       call mp_bcast( nosym,             ionode_id, intra_image_comm )
       call mp_bcast( nosym_evc,         ionode_id, intra_image_comm )
       call mp_bcast( noinv,             ionode_id, intra_image_comm )
       call mp_bcast( force_symmorphic,  ionode_id, intra_image_comm )
       call mp_bcast( use_all_frac,      ionode_id, intra_image_comm )
       call mp_bcast( ecfixed,           ionode_id, intra_image_comm )
       call mp_bcast( qcutz,             ionode_id, intra_image_comm )
       call mp_bcast( q2sigma,           ionode_id, intra_image_comm )
       call mp_bcast( input_dft,         ionode_id, intra_image_comm )

       ! ... exx

       call mp_bcast( ace,                 ionode_id, intra_image_comm )
       call mp_bcast( localization_thr,    ionode_id, intra_image_comm )
       call mp_bcast( scdm,                ionode_id, intra_image_comm )
       call mp_bcast( scdmden,             ionode_id, intra_image_comm )
       call mp_bcast( scdmgrd,             ionode_id, intra_image_comm )
       call mp_bcast( nscdm,               ionode_id, intra_image_comm )
       call mp_bcast( n_proj,              ionode_id, intra_image_comm )
       call mp_bcast( nqx1,                   ionode_id, intra_image_comm )
       call mp_bcast( nqx2,                   ionode_id, intra_image_comm )
       call mp_bcast( nqx3,                   ionode_id, intra_image_comm )
       call mp_bcast( exx_fraction,           ionode_id, intra_image_comm )
       call mp_bcast( screening_parameter,    ionode_id, intra_image_comm ) 
       call mp_bcast( gau_parameter,          ionode_id, intra_image_comm )
       call mp_bcast( exxdiv_treatment,       ionode_id, intra_image_comm )
       call mp_bcast( x_gamma_extrapolation,  ionode_id, intra_image_comm )
       call mp_bcast( yukawa,                 ionode_id, intra_image_comm )
       call mp_bcast( ecutvcut,               ionode_id, intra_image_comm )
       call mp_bcast( ecutfock,               ionode_id, intra_image_comm )
       !
       call mp_bcast( starting_charge,        ionode_id, intra_image_comm )
       call mp_bcast( starting_magnetization, ionode_id, intra_image_comm )
       call mp_bcast( starting_ns_eigenvalue, ionode_id, intra_image_comm )
       call mp_bcast( u_projection_type,      ionode_id, intra_image_comm )
       call mp_bcast( lda_plus_u,             ionode_id, intra_image_comm )
       call mp_bcast( lda_plus_u_kind,        ionode_id, intra_image_comm )
       call mp_bcast( hubbard_u,              ionode_id, intra_image_comm )
       call mp_bcast( hubbard_u_back,         ionode_id, intra_image_comm )
       call mp_bcast( hubbard_j0,             ionode_id, intra_image_comm )
       call mp_bcast( hubbard_j,              ionode_id, intra_image_comm )
       call mp_bcast( hubbard_v,              ionode_id, intra_image_comm )
       call mp_bcast( hubbard_alpha,          ionode_id, intra_image_comm )
       call mp_bcast( hubbard_alpha_back,     ionode_id, intra_image_comm )
       call mp_bcast( hubbard_beta,           ionode_id, intra_image_comm )
       call mp_bcast( hub_pot_fix,            ionode_id,intra_image_comm )
       call mp_bcast( hubbard_parameters,     ionode_id,intra_image_comm )
       call mp_bcast( reserv,                 ionode_id,intra_image_comm )
       call mp_bcast( reserv_back,            ionode_id,intra_image_comm )
       call mp_bcast( backall,                ionode_id,intra_image_comm )
       call mp_bcast( lback,                  ionode_id,intra_image_comm )
       call mp_bcast( l1back,                 ionode_id,intra_image_comm )
       call mp_bcast( step_pen,               ionode_id, intra_image_comm )
       call mp_bcast( a_pen,                  ionode_id, intra_image_comm )
       call mp_bcast( sigma_pen,              ionode_id, intra_image_comm )
       call mp_bcast( alpha_pen,              ionode_id, intra_image_comm )
       call mp_bcast( edir,                   ionode_id, intra_image_comm )
       call mp_bcast( emaxpos,                ionode_id, intra_image_comm )
       call mp_bcast( eopreg,                 ionode_id, intra_image_comm )
       call mp_bcast( eamp,                   ionode_id, intra_image_comm )
       call mp_bcast( la2f,                   ionode_id, intra_image_comm )
       !
       ! ... non collinear broadcast
       !
       call mp_bcast( lspinorb,                  ionode_id, intra_image_comm )
       call mp_bcast( lforcet,                   ionode_id, intra_image_comm )
       call mp_bcast( starting_spin_angle,       ionode_id, intra_image_comm )
       call mp_bcast( noncolin,                  ionode_id, intra_image_comm )
       call mp_bcast( angle1,                    ionode_id, intra_image_comm )
       call mp_bcast( angle2,                    ionode_id, intra_image_comm )
       call mp_bcast( report,                    ionode_id, intra_image_comm )
       call mp_bcast( constrained_magnetization, ionode_id, intra_image_comm )
       call mp_bcast( b_field,                   ionode_id, intra_image_comm )
       call mp_bcast( fixed_magnetization,       ionode_id, intra_image_comm )
       call mp_bcast( lambda,                    ionode_id, intra_image_comm )
       !
       call mp_bcast( assume_isolated,           ionode_id, intra_image_comm )
       call mp_bcast( one_atom_occupations,      ionode_id, intra_image_comm )
       call mp_bcast( spline_ps,                 ionode_id, intra_image_comm )
       !
       call mp_bcast( vdw_corr,                  ionode_id, intra_image_comm )
       call mp_bcast( ts_vdw,                    ionode_id, intra_image_comm )
       call mp_bcast( ts_vdw_isolated,           ionode_id, intra_image_comm )
       call mp_bcast( ts_vdw_econv_thr,          ionode_id, intra_image_comm )
       call mp_bcast( london,                    ionode_id, intra_image_comm )
       call mp_bcast( london_s6,                 ionode_id, intra_image_comm )
       call mp_bcast( london_rcut,               ionode_id, intra_image_comm )
       call mp_bcast( london_c6,                 ionode_id, intra_image_comm )
       call mp_bcast( london_rvdw,               ionode_id, intra_image_comm )
       call mp_bcast( xdm,                       ionode_id, intra_image_comm )
       call mp_bcast( xdm_a1,                    ionode_id, intra_image_comm )
       call mp_bcast( xdm_a2,                    ionode_id, intra_image_comm )
       !
       call mp_bcast( no_t_rev,                  ionode_id, intra_image_comm )
       !
       ! ... esm method broadcast
       !
       call mp_bcast( esm_bc,             ionode_id, intra_image_comm )
       call mp_bcast( esm_efield,         ionode_id, intra_image_comm )
       call mp_bcast( esm_w,              ionode_id, intra_image_comm )
       call mp_bcast( esm_a,              ionode_id, intra_image_comm )
       call mp_bcast( esm_zb,             ionode_id, intra_image_comm )
       call mp_bcast( esm_nfit,           ionode_id, intra_image_comm )
       call mp_bcast( esm_debug,          ionode_id, intra_image_comm )
       call mp_bcast( esm_debug_gpmax,    ionode_id, intra_image_comm )
       !
       ! ... fcp
       !
       call mp_bcast( fcp_mu,          ionode_id, intra_image_comm )
       call mp_bcast( fcp_mass,        ionode_id, intra_image_comm )
       call mp_bcast( fcp_tempw,       ionode_id, intra_image_comm )
       call mp_bcast( fcp_relax,       ionode_id, intra_image_comm )
       call mp_bcast( fcp_relax_step,  ionode_id, intra_image_comm )
       call mp_bcast( fcp_relax_crit,  ionode_id, intra_image_comm )
       call mp_bcast( fcp_mdiis_size,  ionode_id, intra_image_comm )
       call mp_bcast( fcp_mdiis_step,  ionode_id, intra_image_comm )
       !
       !
       ! ... space group information
       !
       call mp_bcast( space_group,        ionode_id, intra_image_comm )
       call mp_bcast( uniqueb,            ionode_id, intra_image_comm )
       call mp_bcast( origin_choice,      ionode_id, intra_image_comm )
       call mp_bcast( rhombohedral,       ionode_id, intra_image_comm )
       !
       ! tb - gate broadcast
       !
       call mp_bcast( zgate,              ionode_id, intra_image_comm )
       call mp_bcast( relaxz,             ionode_id, intra_image_comm )
       call mp_bcast( block,              ionode_id, intra_image_comm )
       call mp_bcast( block_1,            ionode_id, intra_image_comm )
       call mp_bcast( block_2,            ionode_id, intra_image_comm )
       call mp_bcast( block_height,       ionode_id, intra_image_comm )

       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist electrons
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine electrons_bcast()
       !-----------------------------------------------------------------------
       !
       use io_global, only : ionode_id
       use mp,        only : mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       call mp_bcast( emass,                ionode_id, intra_image_comm )
       call mp_bcast( emass_cutoff,         ionode_id, intra_image_comm )
       call mp_bcast( orthogonalization,    ionode_id, intra_image_comm )
       call mp_bcast( electron_maxstep,     ionode_id, intra_image_comm )
       call mp_bcast( scf_must_converge,    ionode_id, intra_image_comm )
       call mp_bcast( ortho_eps,            ionode_id, intra_image_comm )
       call mp_bcast( ortho_max,            ionode_id, intra_image_comm )
       call mp_bcast( electron_dynamics,    ionode_id, intra_image_comm )
       call mp_bcast( electron_damping,     ionode_id, intra_image_comm )
       call mp_bcast( electron_velocities,  ionode_id, intra_image_comm )
       call mp_bcast( electron_temperature, ionode_id, intra_image_comm )
       call mp_bcast( conv_thr,             ionode_id, intra_image_comm )
       call mp_bcast( ekincw,               ionode_id, intra_image_comm )
       call mp_bcast( fnosee,               ionode_id, intra_image_comm )
       call mp_bcast( startingwfc,          ionode_id, intra_image_comm )
       call mp_bcast( ampre,                ionode_id, intra_image_comm )
       call mp_bcast( grease,               ionode_id, intra_image_comm )
       call mp_bcast( startingpot,          ionode_id, intra_image_comm )
       call mp_bcast( diis_size,            ionode_id, intra_image_comm )
       call mp_bcast( diis_nreset,          ionode_id, intra_image_comm )
       call mp_bcast( diis_hcut,            ionode_id, intra_image_comm )
       call mp_bcast( diis_wthr,            ionode_id, intra_image_comm )
       call mp_bcast( diis_delt,            ionode_id, intra_image_comm )
       call mp_bcast( diis_maxstep,         ionode_id, intra_image_comm )
       call mp_bcast( diis_rot,             ionode_id, intra_image_comm )
       call mp_bcast( diis_fthr,            ionode_id, intra_image_comm )
       call mp_bcast( diis_temp,            ionode_id, intra_image_comm )
       call mp_bcast( diis_achmix,          ionode_id, intra_image_comm )
       call mp_bcast( diis_g0chmix,         ionode_id, intra_image_comm )
       call mp_bcast( diis_g1chmix,         ionode_id, intra_image_comm )
       call mp_bcast( diis_nchmix,          ionode_id, intra_image_comm )
       call mp_bcast( diis_nrot,            ionode_id, intra_image_comm )
       call mp_bcast( diis_rothr,           ionode_id, intra_image_comm )
       call mp_bcast( diis_ethr,            ionode_id, intra_image_comm )
       call mp_bcast( diis_chguess,         ionode_id, intra_image_comm )
       call mp_bcast( mixing_fixed_ns,      ionode_id, intra_image_comm )
       call mp_bcast( mixing_mode,          ionode_id, intra_image_comm )
       call mp_bcast( mixing_beta,          ionode_id, intra_image_comm )
       call mp_bcast( mixing_ndim,          ionode_id, intra_image_comm )
       call mp_bcast( tqr,                  ionode_id, intra_image_comm )
       call mp_bcast( tq_smoothing,         ionode_id, intra_image_comm )
       call mp_bcast( tbeta_smoothing,      ionode_id, intra_image_comm )
       call mp_bcast( diagonalization,      ionode_id, intra_image_comm )
       call mp_bcast( diago_thr_init,       ionode_id, intra_image_comm )
       call mp_bcast( diago_cg_maxiter,     ionode_id, intra_image_comm )
       call mp_bcast( diago_ppcg_maxiter,   ionode_id, intra_image_comm )
       call mp_bcast( diago_david_ndim,     ionode_id, intra_image_comm )
       call mp_bcast( diago_full_acc,       ionode_id, intra_image_comm )
       call mp_bcast( sic,                  ionode_id, intra_image_comm )
       call mp_bcast( sic_epsilon ,         ionode_id, intra_image_comm )
       call mp_bcast( sic_alpha   ,         ionode_id, intra_image_comm )
       call mp_bcast( force_pairing ,       ionode_id, intra_image_comm )
       !
       ! ... ensemble-dft
       !
       call mp_bcast( fermi_energy,       ionode_id, intra_image_comm )
       call mp_bcast( n_inner,            ionode_id, intra_image_comm )
       call mp_bcast( niter_cold_restart, ionode_id, intra_image_comm )
       call mp_bcast( lambda_cold,        ionode_id, intra_image_comm )
       call mp_bcast( rotation_dynamics,  ionode_id, intra_image_comm )
       call mp_bcast( occupation_dynamics,ionode_id, intra_image_comm )
       call mp_bcast( rotmass,            ionode_id, intra_image_comm )
       call mp_bcast( occmass,            ionode_id, intra_image_comm )
       call mp_bcast( rotation_damping,   ionode_id, intra_image_comm )
       call mp_bcast( occupation_damping, ionode_id, intra_image_comm )
       !
       ! ... conjugate gradient
       !
       call mp_bcast( tcg,     ionode_id, intra_image_comm )
       call mp_bcast( maxiter, ionode_id, intra_image_comm )
       call mp_bcast( etresh,  ionode_id, intra_image_comm )
       call mp_bcast( passop,  ionode_id, intra_image_comm )
       call mp_bcast( niter_cg_restart, ionode_id, intra_image_comm )
       !
       ! ... electric field
       !
       call mp_bcast( epol,   ionode_id, intra_image_comm )
       call mp_bcast( efield, ionode_id, intra_image_comm )
       !
       call mp_bcast( epol2,   ionode_id, intra_image_comm )
       call mp_bcast( efield2, ionode_id, intra_image_comm )
       call mp_bcast( efield_cart,   ionode_id, intra_image_comm )
       call mp_bcast( efield_phase,   ionode_id, intra_image_comm )
       !
       ! ... occupation constraints ...
       !
       call mp_bcast( occupation_constraints, ionode_id, intra_image_comm )
       !
       ! ... real space ...
       call mp_bcast( real_space,         ionode_id, intra_image_comm )
       call mp_bcast( adaptive_thr,       ionode_id, intra_image_comm )
       call mp_bcast( conv_thr_init,      ionode_id, intra_image_comm )
       call mp_bcast( conv_thr_multi,     ionode_id, intra_image_comm )
       !
       ! ... cp-bo ...
       call mp_bcast( tcpbo,                 ionode_id, intra_image_comm )
       call mp_bcast( emass_emin,            ionode_id, intra_image_comm )
       call mp_bcast( emass_cutoff_emin,     ionode_id, intra_image_comm )
       call mp_bcast( electron_damping_emin, ionode_id, intra_image_comm )
       call mp_bcast( dt_emin,               ionode_id, intra_image_comm )
       !
       return
       !
     end subroutine
     !
     !
     !=----------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist ions
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine ions_bcast()
       !-----------------------------------------------------------------------
       !
       use io_global, only: ionode_id
       use mp,        only: mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       call mp_bcast( ion_dynamics,      ionode_id, intra_image_comm )
       call mp_bcast( ion_radius,        ionode_id, intra_image_comm )
       call mp_bcast( ion_damping,       ionode_id, intra_image_comm )
       call mp_bcast( ion_positions,     ionode_id, intra_image_comm )
       call mp_bcast( ion_velocities,    ionode_id, intra_image_comm )
       call mp_bcast( ion_temperature,   ionode_id, intra_image_comm )
       call mp_bcast( tempw,             ionode_id, intra_image_comm )
       call mp_bcast( fnosep,            ionode_id, intra_image_comm )
       call mp_bcast( nhgrp,             ionode_id, intra_image_comm )
       call mp_bcast( fnhscl,            ionode_id, intra_image_comm )
       call mp_bcast( nhpcl,             ionode_id, intra_image_comm )
       call mp_bcast( nhptyp,            ionode_id, intra_image_comm )
       call mp_bcast( ndega,             ionode_id, intra_image_comm )
       call mp_bcast( tranp,             ionode_id, intra_image_comm )
       call mp_bcast( amprp,             ionode_id, intra_image_comm )
       call mp_bcast( greasp,            ionode_id, intra_image_comm )
       call mp_bcast( tolp,              ionode_id, intra_image_comm )
       call mp_bcast( ion_nstepe,        ionode_id, intra_image_comm )
       call mp_bcast( ion_maxstep,       ionode_id, intra_image_comm )
       call mp_bcast( delta_t,           ionode_id, intra_image_comm )
       call mp_bcast( nraise,            ionode_id, intra_image_comm )
       call mp_bcast( refold_pos,        ionode_id, intra_image_comm )
       call mp_bcast( remove_rigid_rot,  ionode_id, intra_image_comm )
       call mp_bcast( upscale,           ionode_id, intra_image_comm )
       call mp_bcast( pot_extrapolation, ionode_id, intra_image_comm )
       call mp_bcast( wfc_extrapolation, ionode_id, intra_image_comm )
       !
       ! ... bfgs
       !
       call mp_bcast( bfgs_ndim,        ionode_id, intra_image_comm )
       call mp_bcast( trust_radius_max, ionode_id, intra_image_comm )
       call mp_bcast( trust_radius_min, ionode_id, intra_image_comm )
       call mp_bcast( trust_radius_ini, ionode_id, intra_image_comm )
       call mp_bcast( w_1,              ionode_id, intra_image_comm )
       call mp_bcast( w_2,              ionode_id, intra_image_comm )
       !
       call mp_bcast(l_mplathe,         ionode_id, intra_image_comm )
       call mp_bcast(n_muller,          ionode_id, intra_image_comm ) 
       call mp_bcast(np_muller,         ionode_id, intra_image_comm )
       call mp_bcast(l_exit_muller,     ionode_id, intra_image_comm )


       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist cell
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine cell_bcast()
       !-----------------------------------------------------------------------
       !
       use io_global, only: ionode_id
       use mp, only: mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       call mp_bcast( cell_parameters,  ionode_id, intra_image_comm )
       call mp_bcast( cell_dynamics,    ionode_id, intra_image_comm )
       call mp_bcast( cell_velocities,  ionode_id, intra_image_comm )
       call mp_bcast( cell_dofree,      ionode_id, intra_image_comm )
       call mp_bcast( press,            ionode_id, intra_image_comm )
       call mp_bcast( wmass,            ionode_id, intra_image_comm )
       call mp_bcast( cell_temperature, ionode_id, intra_image_comm )
       call mp_bcast( temph,            ionode_id, intra_image_comm )
       call mp_bcast( fnoseh,           ionode_id, intra_image_comm )
       call mp_bcast( greash,           ionode_id, intra_image_comm )
       call mp_bcast( cell_factor,      ionode_id, intra_image_comm )
       call mp_bcast( cell_nstepe,      ionode_id, intra_image_comm )
       call mp_bcast( cell_damping,     ionode_id, intra_image_comm )
       call mp_bcast( press_conv_thr,   ionode_id, intra_image_comm )
       call mp_bcast( treinit_gvecs,    ionode_id, intra_image_comm )
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist press_ai
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     subroutine press_ai_bcast()
       !----------------------------------------------------------------------
       !
       use io_global, only: ionode_id
       use mp,        only: mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       !
       call mp_bcast( abivol, ionode_id, intra_image_comm )
       call mp_bcast( abisur, ionode_id, intra_image_comm )
       call mp_bcast( t_gauss, ionode_id, intra_image_comm )
       call mp_bcast( cntr, ionode_id, intra_image_comm )
       call mp_bcast( p_ext, ionode_id, intra_image_comm )
       call mp_bcast( surf_t, ionode_id, intra_image_comm )
       call mp_bcast( pvar, ionode_id, intra_image_comm )
       call mp_bcast( p_in, ionode_id, intra_image_comm )
       call mp_bcast( p_fin, ionode_id, intra_image_comm )
       call mp_bcast( delta_eps, ionode_id, intra_image_comm )
       call mp_bcast( delta_sigma, ionode_id, intra_image_comm )
       call mp_bcast( fill_vac, ionode_id, intra_image_comm )
       call mp_bcast( scale_at, ionode_id, intra_image_comm )
       call mp_bcast( n_cntr, ionode_id, intra_image_comm )
       call mp_bcast( axis, ionode_id, intra_image_comm )
       call mp_bcast( rho_thr, ionode_id, intra_image_comm )
       call mp_bcast( dthr, ionode_id, intra_image_comm )
       call mp_bcast( step_rad, ionode_id, intra_image_comm )
       call mp_bcast( jellium, ionode_id, intra_image_comm )
       call mp_bcast( r_j, ionode_id, intra_image_comm )
       call mp_bcast( h_j, ionode_id, intra_image_comm )
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist wannier
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine wannier_bcast()
       !-----------------------------------------------------------------------
       !
       use io_global, only: ionode_id
       use mp,        only: mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       call mp_bcast( wf_efield,   ionode_id, intra_image_comm )
       call mp_bcast( wf_switch,   ionode_id, intra_image_comm )
       call mp_bcast( sw_len,      ionode_id, intra_image_comm )
       call mp_bcast( efx0,        ionode_id, intra_image_comm )
       call mp_bcast( efy0,        ionode_id, intra_image_comm )
       call mp_bcast( efz0,        ionode_id, intra_image_comm )
       call mp_bcast( efx1,        ionode_id, intra_image_comm )
       call mp_bcast( efy1,        ionode_id, intra_image_comm )
       call mp_bcast( efz1,        ionode_id, intra_image_comm )
       call mp_bcast( wfsd,        ionode_id, intra_image_comm )
       call mp_bcast( wfdt,        ionode_id, intra_image_comm )
       call mp_bcast( maxwfdt,     ionode_id, intra_image_comm )
       call mp_bcast( wf_q,        ionode_id, intra_image_comm )
       call mp_bcast( wf_friction, ionode_id, intra_image_comm )
       call mp_bcast( nit,         ionode_id, intra_image_comm )
       call mp_bcast( nsd,         ionode_id, intra_image_comm )
       call mp_bcast( nsteps,      ionode_id, intra_image_comm )
       call mp_bcast( tolw,        ionode_id, intra_image_comm )
       call mp_bcast( adapt,       ionode_id, intra_image_comm )
       call mp_bcast( calwf,       ionode_id, intra_image_comm )
       call mp_bcast( nwf,         ionode_id, intra_image_comm )
       call mp_bcast( wffort,      ionode_id, intra_image_comm )
       call mp_bcast( writev,      ionode_id, intra_image_comm )
!=================================================================
!exx_wf related
       call mp_bcast( exx_neigh,       ionode_id, intra_image_comm )
       call mp_bcast( exx_poisson_eps, ionode_id, intra_image_comm )
       call mp_bcast( exx_dis_cutoff,  ionode_id, intra_image_comm )
       call mp_bcast( exx_ps_rcut_self, ionode_id, intra_image_comm )
       call mp_bcast( exx_ps_rcut_pair, ionode_id, intra_image_comm )
       call mp_bcast( exx_me_rcut_self, ionode_id, intra_image_comm )
       call mp_bcast( exx_me_rcut_pair, ionode_id, intra_image_comm )
       call mp_bcast( exx_use_cube_domain, ionode_id, intra_image_comm )
       call mp_bcast( vnbsp,       ionode_id, intra_image_comm )
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------------=!
     !
     !  broadcast variables values for namelist wannier_new
     !
     !=----------------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     subroutine wannier_ac_bcast()
       !----------------------------------------------------------------------
       !
       use io_global, only: ionode_id
       use mp,        only: mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       !
       call mp_bcast( plot_wannier,ionode_id, intra_image_comm )
       call mp_bcast( use_energy_int,ionode_id, intra_image_comm )
       call mp_bcast( print_wannier_coeff,ionode_id, intra_image_comm )
       call mp_bcast( nwan,        ionode_id, intra_image_comm )
       call mp_bcast( plot_wan_num,ionode_id, intra_image_comm )
       call mp_bcast( plot_wan_spin,ionode_id, intra_image_comm )
!       call mp_bcast( wan_data,ionode_id, intra_image_comm )
       call mp_bcast( constrain_pot,   ionode_id, intra_image_comm )
       return
       !
     end subroutine

     !
     !=----------------------------------------------------------------------=!
     !
     !  check input values for namelist control
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine control_checkin( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = ' control_checkin '
       integer           :: i
       logical           :: allowed = .false.
       !
       !
       do i = 1, size( calculation_allowed )
          if( trim(calculation) == calculation_allowed(i) ) allowed = .true.
       end do
       if( .not. allowed ) &
          call errore( sub_name, ' calculation "'// &
                       & trim(calculation)//'" not allowed ',1)
       if( ndr < 50 ) &
          call errore( sub_name,' ndr out of range ', 1 )
       if( ndw > 0 .and. ndw < 50 ) &
          call errore( sub_name,' ndw out of range ', 1 )
       if( nstep < 0 ) &
          call errore( sub_name,' nstep out of range ', 1 )
       if( iprint < 1 ) &
          call errore( sub_name,' iprint out of range ', 1 )

       if( prog == 'pw' ) then
         if( isave > 0 ) &
           call infomsg( sub_name,' isave not used in pw ' )
       else
         if( isave < 1 ) &
           call errore( sub_name,' isave out of range ', 1 )
       end if

       if( dt < 0.0_dp ) &
          call errore( sub_name,' dt out of range ', 1 )
       if( max_seconds < 0.0_dp ) &
          call errore( sub_name,' max_seconds out of range ', 1 )

       if( ekin_conv_thr < 0.0_dp ) then
          if( prog == 'pw' ) then
            call infomsg( sub_name,' ekin_conv_thr not used in pw ')
          else
            call errore( sub_name,' ekin_conv_thr out of range ', 1 )
          end if
       end if

       if( etot_conv_thr < 0.0_dp ) &
          call errore( sub_name,' etot_conv_thr out of range ', 1 )
       if( forc_conv_thr < 0.0_dp ) &
          call errore( sub_name,' forc_conv_thr out of range ', 1 )
       if( prog == 'cp' ) then
          if( dipfield ) &
             call infomsg( sub_name,' dipfield not yet implemented ')
          if( lberry ) &
             call infomsg( sub_name,' lberry not implemented yet ')
          if( gdir /= 0 ) &
             call infomsg( sub_name,' gdir not used ')
          if( nppstr /= 0 ) &
             call infomsg( sub_name,' nppstr not used ')
       end if
       !
       if( prog == 'pw' .and. trim( restart_mode ) == 'reset_counters' ) then
         call infomsg ( sub_name, ' restart_mode == reset_counters' // &
                    & ' not implemented in pw ' )
       end if
       !
       if( refg < 0 ) &
         call errore( sub_name, ' wrong table interval refg ', 1 )
       !
       if( ( prog == 'cp' ) .and. ( trim(memory) == 'small' ) .and. wf_collect ) &
         call errore( sub_name, ' wf_collect = .true. is not allowed with memory = small ', 1 )

       allowed = .false.
       do i = 1, size( memory_allowed )
          if( trim(memory) == memory_allowed(i) ) allowed = .true.
       end do
       if( .not. allowed ) &
          call errore(sub_name, ' memory "' // trim(memory)//'" not allowed',1)
       ! tb
       if ( gate .and. tefield .and. (.not. dipfield) ) &
          call errore(sub_name, ' gate cannot be used with tefield if dipole correction is not active', 1)
       if ( gate .and. dipfield .and. (.not. tefield) ) &
          call errore(sub_name, ' dipole correction is not active if tefield = .false.', 1)

       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  check input values for namelist system
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine system_checkin( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = ' system_checkin '
       integer           :: i
       logical           :: allowed
       !
       !
       if( ( ibrav /= 0 ) .and. (celldm(1) == 0.0_dp) .and. ( a == 0.0_dp ) ) &
           call errore( ' iosys ', &
                      & ' invalid lattice parameters ( celldm or a )', 1 )
       !
       if( nat < 0 ) &
          call errore( sub_name ,' nat less than zero ', max( nat, 1) )
       !
       if( ntyp < 0 ) &
          call errore( sub_name ,' ntyp less than zero ', max( ntyp, 1) )
       if( ntyp < 0 .or. ntyp > nsx ) &
          call errore( sub_name , &
                       & ' ntyp too large, increase nsx ', max( ntyp, 1) )
       !
       if( nspin < 1 .or. nspin > 4 .or. nspin == 3 ) &
          call errore( sub_name ,' nspin out of range ', max(nspin, 1 ) )
       !
       if( ecutwfc < 0.0_dp ) &
          call errore( sub_name ,' ecutwfc out of range ',1)
       if( ecutrho < 0.0_dp ) &
          call errore( sub_name ,' ecutrho out of range ',1)
       !
       if( prog == 'cp' ) then
          if( degauss /= 0.0_dp ) &
             call infomsg( sub_name ,' degauss is not used in cp ')
       end if
       !
       if( ecfixed < 0.0_dp ) &
          call errore( sub_name ,' ecfixed out of range ',1)
       if( qcutz < 0.0_dp ) &
          call errore( sub_name ,' qcutz out of range ',1)
       if( q2sigma < 0.0_dp ) &
          call errore( sub_name ,' q2sigma out of range ',1)
       if( prog == 'cp' ) then
          if( any(starting_magnetization /= sm_not_set ) ) &
             call infomsg( sub_name ,&
                          & ' starting_magnetization is not used in cp ')
          if( la2f ) &
             call infomsg( sub_name ,' la2f is not used in cp ')
          if( any(hubbard_alpha /= 0.0_dp) ) &
             call infomsg( sub_name ,' hubbard_alpha is not used in cp ')
          if( nosym ) &
             call infomsg( sub_name ,' nosym not implemented in cp ')
          if( nosym_evc ) &
             call infomsg( sub_name ,' nosym_evc not implemented in cp ')
          if( noinv ) &
             call infomsg( sub_name ,' noinv not implemented in cp ')
       end if
       !
       ! ... control on sic variables
       !
       if ( sic /= 'none' ) then
          !
          if (sic_epsilon > 1.0_dp )  &
             call errore( sub_name, &
                        & ' invalid sic_epsilon, greater than 1.',1 )
          if (sic_epsilon < 0.0_dp )  &
             call errore( sub_name, &
                        & ' invalid sic_epsilon, less than 0 ',1 )
          if (sic_alpha > 1.0_dp )  &
             call errore( sub_name, &
                        & ' invalid sic_alpha, greater than 1.',1 )
          if (sic_alpha < 0.0_dp )  &
             call errore( sub_name, &
                        & ' invalid sic_alpha, less than 0 ',1 )
          !
          if ( .not. force_pairing ) &
             call errore( sub_name, &
                        & ' invalid force_pairing with sic activated', 1 )
          if ( nspin /= 2 ) &
             call errore( sub_name, &
                        & ' invalid nspin with sic activated', 1 )
          if ( tot_magnetization /= 1._dp )  &
             call errore( sub_name, &
                  & ' invalid tot_magnetization_ with sic activated', 1 )
          !
       endif
       !
       ! ... control on exx variables
       !
       do i = 1, size( exxdiv_treatment_allowed )
          if( trim(exxdiv_treatment) == exxdiv_treatment_allowed(i) ) allowed = .true.
       end do
       if( .not. allowed ) call errore(sub_name, &
           ' invalid exxdiv_treatment: '//trim(exxdiv_treatment), 1 )
       !
       if ( trim(exxdiv_treatment) == "yukawa" .and. yukawa <= 0.0 ) &
          call errore(sub_name, ' invalid value for yukawa', 1 )
       !
       if ( trim(exxdiv_treatment) == "vcut_ws" .and. ecutvcut <= 0.0 ) &
          call errore(sub_name, ' invalid value for ecutvcut', 1 )
       !
       if ( x_gamma_extrapolation .and. ( trim(exxdiv_treatment) == "vcut_ws" .or. &
                                          trim(exxdiv_treatment) == "vcut_spherical" ) ) &
          call errore(sub_name, ' x_gamma_extrapolation cannot be used with vcut', 1 )
       !
       ! tb - gate check
       !
       if ( gate .and. tot_charge == 0 ) &
          call errore(sub_name, ' charged plane (gate) to compensate tot_charge of 0', 1)
       return
       !
       ! ... control on fcp variables
       !
       allowed = .false.
       do i = 1, size(fcp_relax_allowed)
          if( trim(fcp_relax) == fcp_relax_allowed(i) ) allowed = .true.
       end do
       if( .not. allowed ) &
          call errore(sub_name, ' fcp_relax '''//trim(fcp_relax)//''' not allowed ', 1)
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  check input values for namelist electrons
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine electrons_checkin( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = ' electrons_checkin '
       integer           :: i
       logical           :: allowed = .false.
       !
       !
       do i = 1, size(electron_dynamics_allowed)
          if( trim(electron_dynamics) == &
              electron_dynamics_allowed(i) ) allowed = .true.
       end do
       if( .not. allowed ) &
          call errore( sub_name, ' electron_dynamics "'//&
                       & trim(electron_dynamics)//'" not allowed ',1)
       if( emass <= 0.0_dp ) &
          call errore( sub_name, ' emass less or equal 0 ',1)
       if( emass_cutoff <= 0.0_dp ) &
          call errore( sub_name, ' emass_cutoff less or equal 0 ',1)
       if( ortho_eps <= 0.0_dp ) &
          call errore( sub_name, ' ortho_eps less or equal 0 ',1)
       if( ortho_max < 1 ) &
          call errore( sub_name, ' ortho_max less than 1 ',1)
       if( fnosee <= 0.0_dp ) &
          call errore( sub_name, ' fnosee less or equal 0 ',1)
       if( ekincw <= 0.0_dp ) &
          call errore( sub_name, ' ekincw less or equal 0 ',1)
       if( occupation_constraints ) &
          call errore( sub_name, ' occupation_constraints not yet implemented ',1)

!
       return
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  check input values for namelist ions
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine ions_checkin( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = ' ions_checkin '
       integer           :: i
       logical           :: allowed = .false.
       !
       !
       allowed = .false.
       do i = 1, size(ion_dynamics_allowed)
          if( trim(ion_dynamics) == ion_dynamics_allowed(i) ) allowed = .true.
       end do
       if( .not. allowed ) &
          call errore( sub_name, ' ion_dynamics "'// &
                       & trim(ion_dynamics)//'" not allowed ',1)
       if( tempw <= 0.0_dp ) &
          call errore( sub_name,' tempw out of range ',1)
       if( fnosep( 1 ) <= 0.0_dp ) &
          call errore( sub_name,' fnosep out of range ',1)
       if( nhpcl > nhclm ) &
          call infomsg ( sub_name,' nhpcl should be less than nhclm')
       if( nhpcl < 0 ) &
          call infomsg ( sub_name,' nhpcl out of range ')
       if( ion_nstepe <= 0 ) &
          call errore( sub_name,' ion_nstepe out of range ',1)
       if( ion_maxstep < 0 ) &
          call errore( sub_name,' ion_maxstep out of range ',1)
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  check input values for namelist cell
     !
     !=----------------------------------------------------------------------=!
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine cell_checkin( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = ' cell_checkin '
       integer           :: i
       logical           :: allowed = .false.
       !
       !
       do i = 1, size(cell_dynamics_allowed)
          if( trim(cell_dynamics) == &
              cell_dynamics_allowed(i) ) allowed = .true.
       end do
       if( .not. allowed ) &
          call errore( sub_name, ' cell_dynamics "'// &
                       trim(cell_dynamics)//'" not allowed ',1)
       if( wmass < 0.0_dp ) &
          call errore( sub_name,' wmass out of range ',1)
       if( prog == 'cp' ) then
          if( cell_factor /= 0.0_dp ) &
             call infomsg( sub_name,' cell_factor not used in cp ')
       end if
       if( cell_nstepe <= 0 ) &
          call errore( sub_name,' cell_nstepe out of range ',1)
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  check input values for namelist wannier
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine wannier_checkin( prog )
       !-----------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = 'wannier_checkin'
       !
       if ( calwf < 1 .or. calwf > 5 ) &
          call errore( sub_name, ' calwf out of range ', 1 )
       !
       if ( wfsd < 1 .or. wfsd > 3 ) &
          call errore( sub_name, ' wfsd out of range ', 1 )      !
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  check input values for namelist wannier_new
     !
     !=----------------------------------------------------------------------=!
     !
     !----------------------------------------------------------------------
     subroutine wannier_ac_checkin( prog )
       !--------------------------------------------------------------------
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = 'wannier_new_checkin'
       !
       !
       if ( nwan > nwanx ) &
          call errore( sub_name, ' nwan out of range ', 1 )

       if ( plot_wan_num < 0 .or. plot_wan_num > nwan ) &
          call errore( sub_name, ' plot_wan_num out of range ', 1 )

       if ( plot_wan_spin < 0 .or. plot_wan_spin > 2 ) &
          call errore( sub_name, ' plot_wan_spin out of range ', 1 )
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  set values according to the "calculation" variable
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine fixval( prog )
       !-----------------------------------------------------------------------
       !
       use constants, only : e2
       !
       implicit none
       !
       character(len=2)  :: prog   ! ... specify the calling program
       character(len=20) :: sub_name = ' fixval '
       !
       !
       select case( trim( calculation ) )
          case ('scf', 'ensemble')
             if( prog == 'cp' ) then
                 electron_dynamics = 'damp'
                 ion_dynamics      = 'none'
                 cell_dynamics     = 'none'
             end if
          case ('nscf', 'bands')
             if( prog == 'cp' ) occupations = 'bogus'
             if( prog == 'cp' ) electron_dynamics = 'damp'
          case ( 'cp-wf' )
             if( prog == 'cp' ) then
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
             end if
             if ( prog == 'pw' ) &
                call errore( sub_name, ' calculation ' // &
                           & trim( calculation ) // ' not implemented ', 1 )
          case ( 'vc-cp-wf' )
             if( prog == 'cp' ) then
                electron_dynamics = 'verlet'
                ion_dynamics      = 'verlet'
                cell_dynamics     = 'pr'
             else if( prog == 'pw' ) then
                call errore( sub_name, ' calculation ' // &
                           & trim( calculation ) // ' not implemented ', 1 )
             end if
             !
!=========================================================================
!lingzhu kong
          case ( 'cp-wf-nscf' )
             if( prog == 'cp' ) then
                occupations       = 'fixed'
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
             end if
             if ( prog == 'pw' ) &
                call errore( sub_name, ' calculation ' // &
                           & trim( calculation ) // ' not implemented ', 1 )
!=========================================================================
          case ('relax')
             if( prog == 'cp' ) then
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
             else if( prog == 'pw' ) then
                ion_dynamics = 'bfgs'
             end if
          case ( 'md', 'cp' )
             if( prog == 'cp' ) then
                electron_dynamics = 'verlet'
                ion_dynamics      = 'verlet'
             else if( prog == 'pw' ) then
                ion_dynamics = 'verlet'
             end if
          case ('vc-relax')
             if( prog == 'cp' ) then
                electron_dynamics = 'damp'
                ion_dynamics      = 'damp'
                cell_dynamics     = 'damp-pr'
             else if( prog == 'pw' ) then
                ion_dynamics = 'bfgs'
                cell_dynamics= 'bfgs'
             end if
          case ( 'vc-md', 'vc-cp' )
             if( prog == 'cp' ) then
                electron_dynamics = 'verlet'
                ion_dynamics      = 'verlet'
                cell_dynamics     = 'pr'
             else if( prog == 'pw' ) then
                ion_dynamics = 'beeman'
             end if
             !
          case default
             !
             call errore( sub_name,' calculation '// &
                        & trim(calculation)//' not implemented ', 1 )
             !
       end select
       !
       if ( prog == 'pw' ) then
          !
          if ( calculation == 'nscf' .or. calculation == 'bands'  ) then
             !
             startingpot = 'file'
             startingwfc = 'atomic+random'
             !
          else if ( restart_mode == "from_scratch" ) then
             !
             startingwfc = 'atomic+random'
             startingpot = 'atomic'
             !
          else
             !
             startingwfc = 'file'
             startingpot = 'file'
             !
          end if
          !
       else if ( prog == 'cp' ) then
          !
          startingwfc = 'random'
          startingpot = ' '
          !
       end if
       !
       if ( trim( sic ) /= 'none' ) then
         force_pairing = ( nspin == 2 .and. ( tot_magnetization==0._dp .or. &
                                              tot_magnetization==1._dp ) )
       end if
       !
       return
       !
     end subroutine
     !
     !=----------------------------------------------------------------------=!
     !
     !  namelist parsing main routine
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     subroutine read_namelists( prog_, unit )
       !-----------------------------------------------------------------------
       !
       !  this routine reads data from standard input and puts them into
       !  module-scope variables (accessible from other routines by including
       !  this module, or the one that contains them)
       !  ----------------------------------------------
       !
       ! ... declare modules
       !
       use io_global, only : ionode, ionode_id
       use mp,        only : mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       !
       ! ... declare variables
       !
       character(len=*) :: prog_  ! specifies the calling program, allowed:
                                  !     prog = 'pw'     pwscf
                                  !     prog = 'cp'     cp
                                  !     prog = 'pw+ipi' pwscf + i-pi
       !
       integer, intent(in), optional :: unit
       !
       ! ... declare other variables
       !
       character(len=2) :: prog
       integer :: ios
       integer :: unit_loc=5
       !
       ! ... end of declarations
       !
       !  ----------------------------------------------
       !
       if(present(unit)) unit_loc = unit
       !
       prog = prog_(1:2) ! allowed: 'pw' or 'cp'
       if( prog /= 'pw' .and. prog /= 'cp' ) &
          call errore( ' read_namelists ', ' unknown calling program ', 1 )
       !
       ! ... default settings for all namelists
       !
       call control_defaults( prog )
       call system_defaults( prog )
       call electrons_defaults( prog )
       call ions_defaults( prog )
       call cell_defaults( prog )
       !
       ! ... here start reading standard input file
       !
       !
       ! ... control namelist
       !
       ios = 0
       if( ionode ) then
          read( unit_loc, control, iostat = ios )
       end if
       call check_namelist_read(ios, unit_loc, "control")
       !
       call control_bcast( )
       call control_checkin( prog )
       !
       ! ... fixval changes some default values according to the value
       ! ... of "calculation" read in control namelist
       !
       call fixval( prog )
       !
       ! ... system namelist
       !
       ios = 0
       if( ionode ) then
          read( unit_loc, system, iostat = ios )
       end if
       call check_namelist_read(ios, unit_loc, "system")
       !
       call system_bcast( )
       !
       call system_checkin( prog )
       !
       ! ... electrons namelist
       !
       ios = 0
       if( ionode ) then
          read( unit_loc, electrons, iostat = ios )
       end if
       call check_namelist_read(ios, unit_loc, "electrons")
       !
       call electrons_bcast( )
       call electrons_checkin( prog )
       !
       ! ... ions namelist - must be read only if ionic motion is expected,
       ! ...                 or if code called by i-pi via run_driver
       !
       ios = 0
       if ( ionode ) then
          if ( ( trim( calculation ) /= 'nscf'  .and. &
                 trim( calculation ) /= 'bands' ) .or. &
               ( trim( prog_ ) == 'pw+ipi' ) ) then
             read( unit_loc, ions, iostat = ios )
          end if
          !
          ! scf might (optionally) have &ions :: ion_positions = 'from_file'
          !
          if ( (ios /= 0) .and. trim( calculation ) == 'scf' ) then
             ! presumably, not found: rewind the file pointer to the location
             ! of the previous present section, in this case electrons
             rewind( unit_loc )
             read( unit_loc, electrons, iostat = ios )
          end if
          !
       end if
       !
       call check_namelist_read(ios, unit_loc, "ions")
       !
       call ions_bcast( )
       call ions_checkin( prog )
       !
       ! ... cell namelist
       !
       ios = 0
       if( ionode ) then
          if( trim( calculation ) == 'vc-relax' .or. &
              trim( calculation ) == 'vc-cp'    .or. &
              trim( calculation ) == 'vc-md'    .or. &
              trim( calculation ) == 'vc-md'    .or. & 
              trim( calculation ) == 'vc-cp-wf') then
             read( unit_loc, cell, iostat = ios )
          end if
       end if
       call check_namelist_read(ios, unit_loc, "cell")
       !
       call cell_bcast()
       call cell_checkin( prog )
       !
       ios = 0
       if( ionode ) then
          if (tabps) then
             read( unit_loc, press_ai, iostat = ios )
          end if
       end if
       call check_namelist_read(ios, unit_loc, "press_ai")
       !
       call press_ai_bcast()
       !
       ! ... wannier namelist
       !
       call wannier_defaults( prog )
       ios = 0
       if( ionode ) then
          if( trim( calculation ) == 'cp-wf'       .or. &
              trim( calculation ) == 'vc-cp-wf'    .or. &
              trim( calculation ) == 'cp-wf-nscf') then
             read( unit_loc, wannier, iostat = ios )
          end if
       end if
       call check_namelist_read(ios, unit_loc, "wannier")
       !
       call wannier_bcast()
       call wannier_checkin( prog )
       !
       ! ... wannier_new namelist
       !
       call wannier_ac_defaults( prog )
       ios = 0
       if( ionode ) then
          if( use_wannier ) then
             read( unit_loc, wannier_ac, iostat = ios )
          end if
       end if
       call check_namelist_read(ios, unit_loc, "wannier_ac")
       !
       call wannier_ac_bcast()
       call wannier_ac_checkin( prog )
       !
       return
       !
     end subroutine read_namelists
     !
     subroutine check_namelist_read(ios, unit_loc, nl_name)
       use io_global, only : ionode, ionode_id
       use mp,        only : mp_bcast
       use mp_images, only : intra_image_comm
       !
       implicit none
       integer,intent(in) :: ios, unit_loc
       character(len=*) :: nl_name
       character(len=512) :: line
       integer :: ios2
       !
       if( ionode ) then
         ios2=0
         if (ios /=0) then
           backspace(unit_loc)
           read(unit_loc,'(a512)', iostat=ios2) line
         end if
       end if

       call mp_bcast( ios2, ionode_id, intra_image_comm )
       if( ios2 /= 0 ) then
          call errore( ' read_namelists ', ' could not find namelist &'//trim(nl_name), 2)
       endif
       !
       call mp_bcast( ios, ionode_id, intra_image_comm )
       call mp_bcast( line, ionode_id, intra_image_comm )
       if( ios /= 0 ) then
          call errore( ' read_namelists ', &
                       ' bad line in namelist &'//trim(nl_name)//&
                       ': "'//trim(line)//'" (error could be in the previous line)',&
                       1 )
       end if
       !
     end subroutine check_namelist_read
     !
end module read_namelists_module
