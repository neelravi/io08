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

!    use iochamp 

    implicit none

    public  ::

    private ::

    interface fdf_bnames
        module procedure names
    end interface


! Dynamic list for parsed_line structures
    type, public :: line_dlist
        character(len=MAX_LENGTH)  :: str
        type(parsed_line), pointer :: pline => null()
    !
        type(line_dlist), pointer  :: next => null()
        type(line_dlist), pointer  :: prev => null()
    end type line_dlist

   contains

   

!  title      title
!  irn        random number seeds (four 4-digit integers)
!  ijas       form of Jastrow. (between 1 and 6, mostly we use 4)
!  isc        form of scaling function for ri,rj,rij in Jastrow (between 1 and 7, mostly use 2,4,6,7,16,17)
!  iperiodic  0  finite system
!             >0 periodic system
!  ibasis     form of basis
!  hb         hbar=0.5 for Hartree units
!  etrial     guess for energy
!  eunit      'Hartree'
!  nstep      number of steps per block
!  nblk       number of blocks
!  nblkeq     number of equilibration blocks
!  nconf      target number of MC configurations in dmc
!  nconf_new  number of new MC configs. saved per processor.
!  idump      dump restart file
!  irstar     restart from restart file
!  isite      call sites to generate starting MC config. in vmc
!  ipr        print level
!  imetro     form of Metropolis (6 is most efficient choice for most systems)
!             1 simple algorithm with force-bias
!             6 accelerated Metropolis algorithm from Cyrus' 1993 PRL
!  delta      step-size for simple algorithm
!  deltar     radial step-size for accelerated algorithm
!  deltat     angular step-size for accelerated algorithm
!  fbias      force-bias.  (Use 1 always).
!  idmc       form of dmc algorithm
!             1  simple dmc
!             2  improved dmc from Umrigar, Nightingale, Runge 1993 JCP
!  ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v
!  nfprod     number of products to undo for estimating population control bias in dmc
!  tau        time-step in dmc
!  nloc       nonlocal pseudopotential
!             0  local
!             1  in Fahy format
!             2  in Troullier-Martins format (unformatted)
!             3  in Troullier-Martins format (formatted)
!  nquad      number of angular quadrature points for nonlocal psp.
!  nelec      number of electrons
!  nup        number of up-spin electrons
!  nctype     number of atom/center types
!  ncent      number of atoms/centers
!  iwctype    specify atom-type for each atom
!  znuc       nuclear charge
!  cent       atom positions
!  ndet       number of determinants in wavefunction
!  nbasis     number of basis functions
!  norb       number of orbitals
!  cdet       coefficients of determinants
!  iworbd     which orbitals enter in which determinants
!  ianalyt_lap analytic laplacian or not
!  ijas     form of Jastrow. (between 1 and 6, mostly we use 4)
!  isc      form of scaling function for ri,rj,rij in Jastrow (between 1 and 7, mostly use 2,4,6,7)
!           2  [1-exp(scalek*r)]/scalek
!           3  [1-exp{-scalek*r-(scalek*r)**2/2}]/scalek
!           4  r/(1+scalek*r)
!           5  r/{1+(scalek*r)**2}**.5
!           6  Short-range version of 2 (range given bu cutjas)
!           7  Short-range version of 4 (range given bu cutjas)
!  nspin2   1,2,3,-1,-2 -> nspin2b=abs(nspin2)
!  nspin2   > 0  nspin2 sets of a, c parms, nspin2b sets of b parms
!              nocuspb=0  parallel e-e cusp conditions satisfied (b=1/2,1/4)
!  nspin2   < 0  -> nspin2=1
!                nspin2=1 sets of a and c parms, nspin2b sets of b parms
!                -1 nocuspb=1 parallel e-e cusp conditions not satisfied (1/2,1/2)
!                -2 nocuspb=0 parallel e-e cusp conditions satisfied (1/2,1/4)
!  nord     order of the polynmial
!  norda    order of the e-n polynmial in Jastrow4
!  nordb    order of the e-e polynmial in Jastrow4
!  nordc    order of the e-e-n polynmial in Jastrow4
!  cjas1    simple jastrow1 (0.5 to satisfy cusps, parallel-spins automatically take half this value)
!  cjas2    simple jastrow1 parameter
!  scalek   scale factor for Jastrow
!  a1,a2    Jastrow parameters for Jastrow2
!  a,b,c    Jastrow parameters for Jastrow3
!  a4,b,c   Jastrow parameters for Jastrow4,5,6
!  cutjas   cutoff for Jastrow4,5,6 if cutjas=6,7
!  rlobx(y) Lobachevsky parameters for Fock expansion