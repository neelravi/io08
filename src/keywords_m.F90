! Licence information goes here
!
!

#if defined HAVE_CONFIG_H
#  include "config.h"
#endif

#define THIS_FILE "keywords.F90"

!=====================================================================
!
! This file is a part of parser module of CHAMP
! It contains the variables, their default values, a short descriptions.
! The variables are initiazed with their default values here.
! These values will be changed by the parsed input
! 
! Ravindra Shinde
!
!=====================================================================

MODULE keywords


    use iso_fortran_env
    use periodic_table, only: atom_t, element    

    implicit none

    public  ::  title
    public  ::  path_pool
    public  ::  file_input, file_output
    public  ::  file_basis 
    public  ::  file_pseudo
    public  ::  file_orbitals
    public  ::  file_determinants
    public  ::  file_jastrow
    public  ::  file_jastrow_deriv

    public  ::  irn,        rand_seed 
    public  ::  ijas,       form_jastrow                            ! form of Jastrow. (between 1 and 6, mostly we use 4)
    public  ::  isc,        form_jastrow_scaling                    !  isc      form of scaling function for ri,rj,rij in Jastrow (between 1 and 7, mostly use 2,4,6,7)
                                                                    !           2  [1-exp(scalek*r)]/scalek
                                                                    !           3  [1-exp{-scalek*r-(scalek*r)**2/2}]/scalek
                                                                    !           4  r/(1+scalek*r)
                                                                    !           5  r/{1+(scalek*r)**2}**.5
                                                                    !           6  Short-range version of 2 (range given bu cutjas)
                                                                    !           7  Short-range version of 4 (range given bu cutjas)
    public  ::  iperiodic,  periodic
    public  ::  ibasis,     form_basis
    public  ::  hb,         hbar
    public  ::  etrial,     energy_trial
    public  ::  eunit,      energy_unit
    public  ::  nstep,      nsteps
    public  ::  nblk,       nblocks                       
    public  ::  nblkeq,     nequil_blocks                           ! number of equilibration blocks
    public  ::  nconf,      nconfigs                                ! target number of MC configurations in dmc
    public  ::  nconf_new,  nconfigs_new                            ! number of new MC configs. saved per processor.    
    public  ::  idump,      unit_restart_dump                       ! unit of a dump restart file
    public  ::  irstar,     restart                                 ! Currently integer :: restart from restart file; proposed logical


    public  ::  isite,      generate_config                         ! call sites to generate starting MC config. in vmc
    public  ::  ipr,        print_level                             ! print level
    public  ::  imetro,     form_metropolis                         ! form of Metropolis (6 is most efficient choice for most systems)
                                                                    !                     1 simple algorithm with force-bias
                                                                    !                     6 accelerated Metropolis algorithm from Cyrus' 1993 PRL
    public  ::  delta,      step_size                               ! step-size for simple algorithm
    public  ::  deltar,     step_size_radial                        ! radial step-size for accelerated algorithm
    public  ::  deltat,     step_size_angular                       ! angular step-size for accelerated algorithm
    public  ::  fbias,      force_bias                              ! force-bias.  (Use 1 always).
    public  ::  idmc,       form_dmc                                ! form of dmc algorithm
                                                                    !             1  simple dmc
                                                                    !             2  improved dmc from Umrigar, Nightingale, Runge 1993 JCP
    public  ::  nfprod                                              ! number of products to undo for estimating population control bias in dmc
    public  ::  tau                                                 ! time-step in dmc
    public  ::  nloc,       pseudo_format                           ! nonlocal pseudopotential
                                                                    !             0  local
                                                                    !             1  in Fahy format
                                                                    !             2  in Troullier-Martins format (unformatted)
                                                                    !             3  in Troullier-Martins format (formatted)
    public  ::  nquad,      nquadrature                             ! number of angular quadrature points for nonlocal psp.
    public  ::  nelec                                               ! number of electrons
    public  ::  nup,        nalpha                                  ! number of up-spin electrons
    public  ::  ndown,      nbeta                                   ! number of down-spin electrons    
    public  ::  nctype,     ntypes_atom                             ! number of atom/center types 
    public  ::  ncent,      natoms, ncenters, ncentres              ! number of atoms/centers
    public  ::  iwctype                                             ! specify atom-type for each atom
    public  ::  znuc,       znuclear, atomic_number                 ! nuclear charge
    public  ::  cent        !atom_coords                             ! atom positions
    public  ::  ndet,       ndeterminants                           ! number of determinants in wavefunction
    public  ::  nbasis                                              ! number of basis functions
    public  ::  norb,       norbitals                               ! number of orbitals
    public  ::  cdet        !det_coeffs                              ! coefficients of determinants
    public  ::  iworbd                                              ! which orbitals enter in which determinants
    public  ::  ianalyt_lap, analytic_laplacian                     ! analytic laplacian or not


    public  ::  nspin2                                              ! 1,2,3,-1,-2 -> nspin2b=abs(nspin2)
                                                                    !  nspin2   > 0  nspin2 sets of a, c parms, nspin2b sets of b parms
                                                                    !              nocuspb=0  parallel e-e cusp conditions satisfied (b=1/2,1/4)
                                                                    !  nspin2   < 0  -> nspin2=1
                                                                    !                nspin2=1 sets of a and c parms, nspin2b sets of b parms
                                                                    !                -1 nocuspb=1 parallel e-e cusp conditions not satisfied (1/2,1/2)
                                                                    !                -2 nocuspb=0 parallel e-e cusp conditions satisfied (1/2,1/4)
    public  ::  nord,       order_poly                              ! order of the polynmial
    public  ::  norda,      order_poly_en                           ! order of the e-n polynmial in Jastrow4
    public  ::  nordb,      order_poly_ee                           ! order of the e-e polynmial in Jastrow4
    public  ::  nordc,      order_poly_een                          ! order of the e-e-n polynmial in Jastrow4        

    public  ::  cjas1                                               ! simple jastrow1 (0.5 to satisfy cusps, parallel-spins automatically take half this value)
    public  ::  cjas2                                               ! simple jastrow1 parameter
    public  ::  scalek,     scaling_jastrow                         ! scale factor for Jastrow
    public  ::  a1,a2                                               ! Jastrow parameters for Jastrow2
    public  ::  a,b,c                                               ! Jastrow parameters for Jastrow3
    public  ::  a4                                                  ! Jastrow parameters for Jastrow4,5,6
    public  ::  cutjas,     cutoff_jastrow                          ! cutoff for Jastrow4,5,6 if cutjas=6,7
    public  ::  itau_eff,   itau_effective

!    Following not yet added 
!    rlobx(y) Lobachevsky parameters for Fock expansion
!    ipq,iacc_rej,icross,icuspg,idiv_v

    private :: sp, dp

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)    


!  declarations

    character(len=132)              ::  title
    character(len=132)              ::  path_pool 
    character(len=132)              ::  file_input, file_output
    character(len=132)              ::  file_basis 
    character(len=132)              ::  file_pseudo
    character(len=132)              ::  file_orbitals
    character(len=132)              ::  file_determinants
    character(len=132)              ::  file_jastrow
    character(len=132)              ::  file_jastrow_deriv

    integer, target                 ::  irn
    integer, pointer                ::  rand_seed => irn

    integer, target                 ::  ijas
    integer, pointer                ::  form_jastrow => ijas

    integer, target                 ::  isc
    integer, pointer                ::  form_jastrow_scaling => isc

    integer, target                 ::  iperiodic
    integer, pointer                ::  periodic => iperiodic

    integer, target                 ::  ibasis
    integer, pointer                ::  form_basis => ibasis

    real(dp), target                ::  hb
    real(dp), pointer               ::  hbar => hb
 
    real(dp), target                ::  etrial
    real(dp), pointer               ::  energy_trial => etrial

    integer, target                 ::  eunit
    integer, pointer                ::  energy_unit => eunit

    integer, target                 ::  nstep
    integer, pointer                ::  nsteps => nstep

    integer, target                 ::  nblk
    integer, pointer                ::  nblocks => nblk

    integer, target                 ::  nblkeq
    integer, pointer                ::  nequil_blocks => nblkeq

    integer, target                 ::  nconf
    integer, pointer                ::  nconfigs => nconf

    integer, target                 ::  nconf_new
    integer, pointer                ::  nconfigs_new => nconf_new

    integer, target                 ::  idump
    integer, pointer                ::  unit_restart_dump => idump

    integer, target                 ::  irstar
    integer, pointer                ::  restart => irstar

    integer, target                 ::  isite
    integer, pointer                ::  generate_config => isite

    integer, target                 ::  ipr
    integer, pointer                ::  print_level => ipr

    integer, target                 ::  imetro
    integer, pointer                ::  form_metropolis => imetro

    integer, target                 ::  delta
    integer, pointer                ::  step_size => delta

    integer, target                 ::  deltar
    integer, pointer                ::  step_size_radial => deltar

    integer, target                 ::  deltat
    integer, pointer                ::  step_size_angular => deltat

    integer, target                 ::  fbias
    integer, pointer                ::  force_bias => fbias

    integer, target                 ::  idmc
    integer, pointer                ::  form_dmc => idmc

    integer                         ::  nfprod                            
    real(sp)                        ::  tau                               

    integer, target                 ::  nloc
    integer, pointer                ::  pseudo_format => nloc

    integer, target                 ::  nquad
    integer, pointer                ::  nquadrature => nquad

    integer                         ::  nelec                                               

    integer, target                 ::  nup
    integer, pointer                ::  nalpha => nup                                 

    integer, target                 ::  ndown
    integer, pointer                ::  nbeta  => ndown

    integer, target                 ::  nctype
    integer, pointer                ::  ntypes_atom => nctype

    integer, target                 ::  ncent
    integer, pointer                ::  natoms => ncent, ncenters  => ncent, ncentres => ncent

    integer                         ::  iwctype                                             

    integer, target                 ::  znuc
    integer, pointer                ::  znuclear => znuc, atomic_number => znuc

    real(dp), allocatable           ::  cent(:,:)
!    real(dp), pointer               ::  atom_coords => cent  !! this might have issues in performance

    integer, target                 ::  ndet
    integer, pointer                ::  ndeterminants => ndet

    integer                         ::  nbasis                                              

    integer, target                 ::  norb
    integer, pointer                ::  norbitals => norb

    real(dp), allocatable           ::  cdet(:)
!    integer, pointer                ::  det_coeffs => cdet   ! issues in performance

    integer                         ::  iworbd                                              

    integer, target                 ::  ianalyt_lap
    integer, pointer                ::  analytic_laplacian => ianalyt_lap

    integer                         ::  nspin2                                              

    integer, target                 ::  nord
    integer, pointer                ::  order_poly => nord

    integer, target                 ::  norda
    integer, pointer                ::  order_poly_en => norda

    integer, target                 ::  nordb
    integer, pointer                ::  order_poly_ee => nordb

    integer, target                 ::  nordc
    integer, pointer                ::  order_poly_een => nordc

    integer, target                 ::  scalek
    integer, pointer                ::  scaling_jastrow => scalek

    integer                         ::  cjas1                                               
    integer                         ::  cjas2                                               

    real(dp)                        ::  a1,a2                                               
    real(dp)                        ::  a,b,c                                               
    real(dp)                        ::  a4


    real(dp), target                ::  cutjas
    real(dp), pointer               ::  cutoff_jastrow => cutjas

    real(dp), target                ::  itau_eff
    real(dp), pointer               ::  itau_effective => itau_eff

end module