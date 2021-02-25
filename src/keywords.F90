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

    implicit none

    public  ::  title 
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
    public  ::  cent,       atom_coords                             ! atom positions
    public  ::  ndet,       ndeterminants                           ! number of determinants in wavefunction
    public  ::  nbasis                                              ! number of basis functions
    public  ::  norb,       norbitals                               ! number of orbitals
    public  ::  cdet,       det_coeffs                              ! coefficients of determinants
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
    public  ::  a4,b,c                                              ! Jastrow parameters for Jastrow4,5,6
    public  ::  cutjas,     cutoff_jastrow                          ! cutoff for Jastrow4,5,6 if cutjas=6,7
    public  ::  itau_eff,   itau_effective

!    Following not yet added 
!    rlobx(y) Lobachevsky parameters for Fock expansion
!    ipq,iacc_rej,icross,icuspg,idiv_v

!    private variables

    private ::

    interface fdf_bnames
        module procedure names
    end interface

!   derived data types

    type, public :: atom_t
        character(len=:)        :: name
        character(len=:)        :: symbol        
        integer(int32)          :: znuclear
        real(real64)            :: atomic_mass
    end type atom_t


!  declarations

    public  ::  title 
    public  ::  irn,        rand_seed 
    public  ::  ijas,       form_jastrow        
    public  ::  isc,        form_jastrow_scaling
    public  ::  iperiodic,  periodic
    public  ::  ibasis,     form_basis
    public  ::  hb,         hbar
    public  ::  etrial,     energy_trial
    public  ::  eunit,      energy_unit
    public  ::  nstep,      nsteps
    public  ::  nblk,       nblocks                       
    public  ::  nblkeq,     nequil_blocks    
    public  ::  nconf,      nconfigs         
    public  ::  nconf_new,  nconfigs_new     
    public  ::  idump,      unit_restart_dump
    public  ::  irstar,     restart          
    public  ::  isite,      generate_config       
    public  ::  ipr,        print_level           
    public  ::  imetro,     form_metropolis       
    public  ::  delta,      step_size             
    public  ::  deltar,     step_size_radial      
    public  ::  deltat,     step_size_angular     
    public  ::  fbias,      force_bias            
    public  ::  idmc,       form_dmc              
    public  ::  nfprod                            
    public  ::  tau                               
    public  ::  nloc,       pseudo_format         
    public  ::  nquad,      nquadrature           
    public  ::  nelec                                               
    public  ::  nup,        nalpha                                  
    public  ::  ndown,      nbeta                                   
    public  ::  nctype,     ntypes_atom                             
    public  ::  ncent,      natoms, ncenters, ncentres              
    public  ::  iwctype                                             
    public  ::  znuc,       znuclear, atomic_number                 
    public  ::  cent,       atom_coords                             
    public  ::  ndet,       ndeterminants                           
    public  ::  nbasis                                              
    public  ::  norb,       norbitals                               
    public  ::  cdet,       det_coeffs                              
    public  ::  iworbd                                              
    public  ::  ianalyt_lap, analytic_laplacian                     
    public  ::  nspin2                                              
    public  ::  nord,       order_poly                              
    public  ::  norda,      order_poly_en                           
    public  ::  nordb,      order_poly_ee                           
    public  ::  nordc,      order_poly_een                          
    public  ::  cjas1                                               
    public  ::  cjas2                                               
    public  ::  scalek,     scaling_jastrow                         
    public  ::  a1,a2                                               
    public  ::  a,b,c                                               
    public  ::  a4,b,c                                              
    public  ::  cutjas,     cutoff_jastrow                          
    public  ::  itau_eff,   itau_effective    


   contains
