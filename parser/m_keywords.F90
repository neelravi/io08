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
    public  ::  file_molecule
    public  ::  file_determinants

    public  ::  etrial,     energy_trial

    public  ::  tau                                                 ! time-step in dmc
    public  ::  nelec                                               ! number of electrons
    public  ::  nup,        nalpha                                  ! number of up-spin electrons
    public  ::  ndown,      nbeta                                   ! number of down-spin electrons    
    public  ::  ncent,      natoms, ncenters, ncentres              ! number of atoms/centers
    public  ::  iwctype                                             ! specify atom-type for each atom

    public  ::  cent        !atom_coords                             ! atom positions
    public  ::  ndet,       ndeterminants                           ! number of determinants in wavefunction
    public  ::  cdet        !det_coeffs                              ! coefficients of determinants
    public  ::  iworbd                                              ! which orbitals enter in which determinants

    public  ::  nspin1


    public  :: optimize_wave
    public  :: optimize_ci 

    public  :: ncore
    public  :: nextorb

    public  :: sr_tau 
    public  :: sr_eps 

    public  :: energy_tol
    public  :: opt_method
    public  :: multiple_adiag    


    private :: sp, dp

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)    


!  declarations

    character(len=132)              ::  title
    character(len=132)              ::  path_pool 
    character(len=132)              ::  file_input, file_output
    character(len=132)              ::  file_basis 
    character(len=132)              ::  file_molecule    
    character(len=132)              ::  file_determinants

    real(dp), target                ::  etrial
    real(dp), pointer               ::  energy_trial => etrial


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
    integer, allocatable            ::  iworbd(:,:)   ! to store orbital mapping in determinants

    real(dp), allocatable           ::  cent(:,:)

    integer, target                 ::  ndet
    integer, pointer                ::  ndeterminants => ndet


    real(dp), allocatable           ::  cdet(:)
    integer                         ::  nspin1

    logical                         :: optimize_wave
    logical                         :: optimize_ci
    logical                         :: multiple_adiag

    integer                         :: ncore
    integer                         :: nextorb

    real(dp)                        :: sr_tau, tau
    real(dp)                        :: sr_eps 
    real(dp)                        :: energy_tol

    character(len=20)               :: opt_method

end module