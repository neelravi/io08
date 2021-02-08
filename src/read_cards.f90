!
! copyright (c) 2002-2014 quantum espresso group
! this file is distributed under the terms of the
! gnu general public license. see the file `license'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
module read_cards_module
   !---------------------------------------------------------------------------
   !
   ! ...  this module handles the reading of cards from standard input
   ! ...  original version written by carlo cavazzoni
   !
   use kinds,     only : dp
   use io_global, only : stdout
   use wy_pos,    only : wypos
   use parser,    only : field_count, read_line, get_field, parse_unit
   use io_global, only : ionode, ionode_id
   !
   use input_parameters
   !
   !
   implicit none
   !
   save
   !
   private
   !
   public :: read_cards
   !
   ! ... end of module-scope declarations
   !
   !  ----------------------------------------------
   !
contains
   !
   ! ... read cards ....
   !
   ! ... subroutines
   !
   !----------------------------------------------------------------------
   subroutine card_default_values( )
      !----------------------------------------------------------------------
      !
      use autopilot, only : init_autopilot
      !
      implicit none
      !
      !
      ! ... mask that control the printing of selected kohn-sham occupied
      ! ... orbitals, default allocation
      !
      call allocate_input_iprnks( 0, nspin )
      nprnks  = 0
      !
      ! ... simulation cell from standard input
      !
      trd_ht = .false.
      rd_ht  = 0.0_dp
      !
      ! ... reference simulation cell from standard input
      !
      ref_cell = .false.
      rd_ref_ht  = 0.0_dp
      !
      ! ... constraints
      !
      nconstr_inp    = 0
      constr_tol_inp = 1.e-6_dp
      !
      ! ... ionic mass initialization
      !
      atom_mass = 0.0_dp
      !
      ! ... k-points
      !
      k_points = 'gamma'
      tk_inp   = .false.
      nkstot   = 1
      nk1      = 0
      nk2      = 0
      nk3      = 0
      k1       = 0
      k2       = 0
      k3       = 0
      !
      ! ... electronic states
      !
      tf_inp = .false.
      !
      ! ... ion_velocities
      !
      tavel = .false.
      !
      call init_autopilot()
      !
      return
      !
   end subroutine card_default_values
   !
   !
   !----------------------------------------------------------------------
   subroutine read_cards ( prog, unit )
      !----------------------------------------------------------------------
      !
      use autopilot, only : card_autopilot
      !
      implicit none
      !
      integer, intent(in), optional  :: unit
      !
      character(len=2)           :: prog   ! calling program ( pw, cp, wa )
      character(len=256)         :: input_line
      character(len=80)          :: card
      character(len=1), external :: capital
      logical                    :: tend
      integer                    :: i
      !
      ! read_line reads from unit parse_unit
      !
      if (present(unit)) then
         parse_unit =  unit
      else
         parse_unit =  5
      end if
      !
      call card_default_values( )
      !
100   call read_line( input_line, end_of_file=tend )
      !
      if( tend ) goto 120
      if( input_line == ' ' .or. input_line(1:1) == '#' .or. &
                                 input_line(1:1) == '!' ) goto 100
      !
      read (input_line, *) card
      !
      do i = 1, len_trim( input_line )
         input_line( i : i ) = capital( input_line( i : i ) )
      enddo
      !
      if ( trim(card) == 'autopilot' ) then
         !
         call card_autopilot( input_line )
         if ( prog == 'pw' .and. ionode ) &
            write( stdout,'(a)') 'warning: card '//trim(input_line)//' ignored'
         !
      elseif ( trim(card) == 'atomic_species' ) then
         !
         call card_atomic_species( input_line )
         !
      elseif ( trim(card) == 'atomic_positions' ) then
         !
         call card_atomic_positions( input_line, prog )
         !
      elseif ( trim(card) == 'atomic_forces' ) then
         !
         call card_atomic_forces( input_line )
         !
      elseif ( trim(card) == 'constraints' ) then
         !
         call card_constraints( input_line )
         !
      elseif ( trim(card) == 'dipole' ) then
         !
         call errore('read_cards','card dipole no longer existing',1)
         !
      elseif ( trim(card) == 'esr' ) then
         !
         call errore('read_cards','card esr no longer existing',1)
         !
      elseif ( trim(card) == 'k_points' ) then
         !
         if ( ( prog == 'cp' ) ) then
            if( ionode ) &
               write( stdout,'(a)') 'warning: card '//trim(input_line)//' ignored'
         else
            call card_kpoints( input_line )
         endif
         !
      elseif ( trim(card) == 'additional_k_points' ) then
         !
         call card_add_kpoints( input_line )
 
      elseif ( trim(card) == 'occupations' ) then
         !
         call card_occupations( input_line )
         !
      elseif ( trim(card) == 'cell_parameters' ) then
         !
         call card_cell_parameters( input_line )
         !
      elseif ( trim(card) == 'ref_cell_parameters' ) then
         !
         call card_ref_cell_parameters( input_line )
         !
      elseif ( trim(card) == 'atomic_velocities' ) then
         !
         call card_ion_velocities( input_line )
         !
      elseif ( trim(card) == 'ksout' ) then
         !
         call card_ksout( input_line )
         if ( ( prog == 'pw' ) .and. ionode ) &
            write( stdout,'(a)') 'warning: card '//trim(input_line)//' ignored'
         !
      elseif ( trim(card) == 'plot_wannier' ) then
         !
         call card_plot_wannier( input_line )

      elseif ( trim(card) == 'wannier_ac' .and. ( prog == 'wa' )) then
         !
         call card_wannier_ac( input_line )

      else
         !
         if ( ionode ) &
            write( stdout,'(a)') 'warning: card '//trim(input_line)//' ignored'
         !
      endif
      !
      ! ... end of loop ... !
      !
      goto 100
      !
120      continue
      !
      return
      !
   end subroutine read_cards

   !
   ! ... description of the allowed input cards
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! atomic_species
   !
   !   set the atomic species been read and their pseudopotential file
   !
   ! syntax:
   !
   !    atomic_specie
   !      label(1)    mass(1)    psfile(1)
   !       ...        ...        ...
   !      label(n)    mass(n)    psfile(n)
   !
   ! example:
   !
   ! atomic_species
   !  o 16.0 o.blyp.upf
   !  h 1.00 h.fpmd.upf
   !
   ! where:
   !
   !      label(i)  ( character(len=4) )  label of the atomic species
   !      mass(i)   ( real )              atomic mass
   !                                      ( in u.m.a, carbon mass is 12.0 )
   !      psfile(i) ( character(len=80) ) file name of the pseudopotential
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_atomic_species( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      integer            :: is, ip, ierr
      character(len=4)   :: lb_pos
      character(len=256) :: psfile
      !
      !
      if ( taspc ) then
         call errore( ' card_atomic_species  ', ' two occurrences', 2 )
      endif
      if ( ntyp > nsx ) then
         call errore( ' card_atomic_species ', ' nsp out of range ', ntyp )
      endif
      !
      do is = 1, ntyp
         !
         call read_line( input_line )
         read( input_line, *, iostat=ierr ) lb_pos, atom_mass(is), psfile
            call errore( ' card_atomic_species ', &
                'cannot read atomic specie from: '//trim(input_line), abs(ierr))
         atom_pfile(is) = trim( psfile )
         lb_pos         = adjustl( lb_pos )
         atom_label(is) = trim( lb_pos )
         !
         do ip = 1, is - 1
            if ( atom_label(ip) == atom_label(is) ) then
               call errore( ' card_atomic_species ', &
                           & ' two occurrences of the same atomic label ', is )
            endif
         enddo
         !
      enddo
      taspc = .true.
      !
      return
      !
   end subroutine card_atomic_species
   !
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! atomic_positions
   !
   !   set the atomic positions in the cell
   !
   ! syntax:
   !
   !   atomic_positions (units_option)
   !     label(1) tau(1,1) tau(2,1) tau(3,1) mbl(1,1) mbl(2,1) mbl(3,1)
   !     label(2) tau(1,2) tau(2,2) tau(3,2) mbl(1,2) mbl(2,2) mbl(3,2)
   !      ...              ...               ...               ... ...
   !     label(n) tau(1,n) tau(2,n) tau(3,n) mbl(1,3) mbl(2,3) mbl(3,3)
   !
   ! example:
   !
   ! atomic_positions (bohr)
   !    o     0.0099    0.0099    0.0000  0 0 0
   !    h     1.8325   -0.2243   -0.0001  1 1 1
   !    h    -0.2243    1.8325    0.0002  1 1 1
   !
   ! where:
   !
   !   units_option == crystal   position are given in scaled units
   !   units_option == bohr      position are given in bohr
   !   units_option == angstrom  position are given in angstrom
   !   units_option == alat      position are given in units of alat
   !
   !   label(k) ( character(len=4) )  atomic type
   !   tau(:,k) ( real )              coordinates  of the k-th atom
   !   mbl(:,k) ( integer )           mbl(i,k) > 0 the i-th coord. of the
   !                                  k-th atom is allowed to be moved
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_atomic_positions( input_line, prog )
      !
      use wrappers, only: feval_infix
      !
      implicit none
      !
      character(len=256) :: input_line
      character(len=2)   :: prog
      character(len=4)   :: lb_pos
      integer            :: ia, k, is, nfield, idx, rep_i
      logical, external  :: matches
      logical            :: tend
      real(dp)           :: inp(3)
      integer            :: fieldused
      !
      integer            :: ifield, ierr
      real(dp)           :: field_value
      character(len=256) :: field_str, error_msg, wp
      !
      !
      if ( tapos ) then
         call errore( 'card_atomic_positions', 'two occurrences', 2 )
      endif
      if ( .not. taspc ) then
         call errore( 'card_atomic_positions', &
                     & 'atomic_species must be present before', 2 )
      endif
      if ( ntyp > nsx ) then
         call errore( 'card_atomic_positions', 'nsp out of range', ntyp )
      endif
      if ( nat < 1 ) then
         call errore( 'card_atomic_positions', 'nat out of range', nat )
      endif
      !
      call allocate_input_ions(ntyp,nat)
      !
      rd_if_pos = 1
      !
      sp_pos = 0
      rd_pos = 0.0_dp
      na_inp = 0
      lsg=.false.
      !
      if ( matches( "crystal_sg", input_line ) ) then
         atomic_positions = 'crystal'
         lsg=.true.
      elseif ( matches( "crystal", input_line ) ) then
         atomic_positions = 'crystal'
      elseif ( matches( "bohr", input_line ) ) then
         atomic_positions = 'bohr'
      elseif ( matches( "angstrom", input_line ) ) then
         atomic_positions = 'angstrom'
      elseif ( matches( "alat", input_line ) ) then
         atomic_positions = 'alat'
      else
         if ( trim( adjustl( input_line ) ) /= 'atomic_positions' ) then
            call errore( 'read_cards ', &
                        & 'unknown option for atomic_position: '&
                        & // input_line, 1 )
         endif
         call infomsg( 'read_cards ', &
            & 'deprecated: no units specified in atomic_positions card' )
         if ( prog == 'cp' ) atomic_positions = 'bohr'
         if ( prog == 'pw' ) atomic_positions = 'alat'
         call infomsg( 'read_cards ', &
            & 'atomic_positions: units set to '//trim(atomic_positions) )
      endif
      !
      reader_loop : do ia = 1,nat
         !
         call read_line( input_line, end_of_file = tend )
         if ( tend ) call errore( 'read_cards', &
                           'end of file reading atomic positions', ia )
         !
         call field_count( nfield, input_line )
         !
         ! read atom symbol (column 1)
         !
         call get_field(1, lb_pos, input_line)
         lb_pos = trim(lb_pos)
         !
         error_msg = 'error while parsing atomic position card.'
         !
         ! read field 2 (atom x coordinate or wyckoff position symbol)
         !
         call get_field(2, field_str, input_line)
         !     
         ! check if position ia is expressed in wyckoff positions
         !
         idx = len_trim(field_str)
         if ( lsg .and. (idx < 4) .and. &
              ( iachar(field_str(idx:idx)) > 64 .and. &
                iachar(field_str(idx:idx)) < 123 ) ) then
            !
            ! wyckoff positions
            !
            if ( nfield < 3 .and. nfield > 8 ) &
            call errore( 'read_cards', 'wrong number of columns ' // &
                           & 'in atomic_positions', ia )
            wp=field_str
            inp(:)=1.d5
            !
            do k = 3,min(nfield,5)
               ! read k-th field (coordinate k-2)
               call get_field(k, field_str, input_line)
               inp(k-2) = feval_infix(ierr, field_str )
               call errore('card_atomic_positions', error_msg, ierr)
            enddo
            !
            call wypos(rd_pos(1,ia),wp,inp,space_group, &
                 uniqueb,rhombohedral,origin_choice)
            !
            ! count how many fields were used to find wyckoff positions
            !
            fieldused=2
            if ( any (rd_pos(1:3,ia)==inp(1)) ) fieldused=fieldused+1
            if ( any (rd_pos(2:3,ia)==inp(2)) ) fieldused=fieldused+1
            if (      rd_pos(  3,ia)==inp(3)  ) fieldused=fieldused+1
            !
         else
            !
            ! no wyckoff positions 
            !
            if ( nfield /= 4 .and. nfield /= 7 ) &
            call errore( 'read_cards', 'wrong number of columns ' // &
                           & 'in atomic_positions', ia )
            !
            ! field just read is coordinate x
            !
            rd_pos(1,ia) = feval_infix(ierr, field_str )
            call errore('card_atomic_positions', error_msg, ierr)
            do k = 3,4
               ! read fields 3 and 4 (atom y and z coordinate)
               call get_field(k, field_str, input_line)
               rd_pos(k-1,ia) = feval_infix(ierr, field_str )
               call errore('card_atomic_positions', error_msg, ierr)
            end do
            !
            fieldused=4
            !
         endif
         ! read constraints if present (last 3 fields)
         if ( nfield-fieldused > 0 .and. nfield-fieldused /= 3 ) &
            call errore( 'read_cards', 'unexpected number of columns ' // &
                           & 'in atomic_positions', ia )
         do k = fieldused+1, nfield
            call get_field(k, field_str, input_line)
            read(field_str, *) rd_if_pos(k-fieldused,ia)
         enddo
         !
         match_label: do is = 1, ntyp
            !
            if ( trim(lb_pos) == trim( atom_label(is) ) ) then
               !
               sp_pos(ia) = is
               exit match_label
               !
            endif
            !
         enddo match_label
         !
         if( ( sp_pos(ia) < 1 ) .or. ( sp_pos(ia) > ntyp ) ) then
            !
            call errore( 'read_cards', 'species '//trim(lb_pos)// &
                           & ' in atomic_positions is nonexistent', ia )
            !
         endif
         !
         is = sp_pos(ia)
         !
         na_inp(is) = na_inp(is) + 1
         !
      enddo reader_loop
      !
      tapos = .true.
      !

      return
      !
   end subroutine card_atomic_positions
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! atomic_forces
   !
   !   read external forces (in atomic units) from standard input
   !
   ! syntax:
   !
   !   atomic_forces
   !     label fx(1) fy(1) fz(1)
   !     .....
   !     label fx(n) fy(n) fz(n)
   !
   ! example:
   !
   !   ???
   !
   ! where:
   !
   !   label (character(len=4))       atomic label
   !   fx(:), fy(:) and fz(:) (real)  x, y and z component of the external force
   !                                  acting on the ions whose coordinate are given
   !                                  in the same line in card atomic_position
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_atomic_forces( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      integer            :: ia, k, nfield
      character(len=4)   :: lb
      !
      !
      if( tforces ) then
         call errore( ' card_atomic_forces ', ' two occurrences ', 2 )
      endif
      !
      if( .not. tapos ) then
         call errore( ' card_atomic_forces ', &
                     & ' atomic_species must be present before ', 2 )
      endif
      !
      rd_for = 0.0_dp
      !
      do ia = 1, nat
         !
         call read_line( input_line )
         call field_count( nfield, input_line )
         if ( nfield == 4 ) then
            read(input_line,*) lb, ( rd_for(k,ia), k = 1, 3 )
         elseif( nfield == 3 ) then
            read(input_line,*) ( rd_for(k,ia), k = 1, 3 )
         else
            call errore( ' iosys ', ' wrong entries in atomic_forces ', ia )
         endif
         !
      enddo
      !
      tforces = .true.
      !
      return
      !
   end subroutine card_atomic_forces
   !
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! k_points
   !
   !   use the specified set of k points
   !
   ! syntax:
   !
   !   k_points (mesh_option)
   !     n
   !     xk(1,1) xk(2,1) xk(3,1) wk(1)
   !     ...     ...     ...     ...
   !     xk(1,n) xk(2,n) xk(3,n) wk(n)
   !
   ! example:
   !
   ! k_points
   !   10
   !    0.1250000  0.1250000  0.1250000   1.00
   !    0.1250000  0.1250000  0.3750000   3.00
   !    0.1250000  0.1250000  0.6250000   3.00
   !    0.1250000  0.1250000  0.8750000   3.00
   !    0.1250000  0.3750000  0.3750000   3.00
   !    0.1250000  0.3750000  0.6250000   6.00
   !    0.1250000  0.3750000  0.8750000   6.00
   !    0.1250000  0.6250000  0.6250000   3.00
   !    0.3750000  0.3750000  0.3750000   1.00
   !    0.3750000  0.3750000  0.6250000   3.00
   !
   ! where:
   !
   !   mesh_option == automatic  k points mesh is generated automatically
   !                             with monkhorst-pack algorithm
   !   mesh_option == crystal    k points mesh is given in stdin in scaled
   !                             units
   !   mesh_option == tpiba      k points mesh is given in stdin in units
   !                             of ( 2 pi / alat )
   !   mesh_option == gamma      only gamma point is used ( default in
   !                             cpmd simulation )
   !   mesh_option == tpiba_b    as tpiba but the weights gives the
   !                             number of points between this point
   !                             and the next
   !   mesh_option == crystal_b  as crystal but the weights gives the
   !                             number of points between this point and
   !                             the next
   !   mesh_option == tpiba_c    the code expects three k points 
   !                             k_0, k_1, k_2 in tpiba units.
   !                             these points define a rectangle
   !                             in reciprocal space with vertices k_0, k_1,
   !                             k_2, k_1+k_2-k_0:  k_0 + \alpha (k_1-k_0)+
   !                             \beta (k_2-k_0) with 0<\alpha,\beta < 1. 
   !                             the code produces a uniform mesh n1 x n2 
   !                             k points in this rectangle. n1 and n2 are 
   !                             the weights of k_1 and k_2. the weight of k_0
   !                             is not used. useful for contour plots of the 
   !                             bands.
   !   mesh_option == crystal_c  as tpiba_c but the k points are given
   !                             in crystal coordinates.
   ! 
   !
   !   n       ( integer )  number of k points
   !   xk(:,i) ( real )     coordinates of i-th k point
   !   wk(i)   ( real )     weights of i-th k point
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_kpoints( input_line )
      !
      use bz_form, only : transform_label_coord
      use cell_base, only : cell_base_init, celldm_cb => celldm
      implicit none
      !
      character(len=256) :: input_line, buffer
      integer            :: i, j
      integer            :: nkaux
      integer, allocatable :: wkaux(:)
      real(dp), allocatable :: xkaux(:,:)
      integer            :: npk_label, nch
      character(len=3), allocatable :: letter(:)
      integer, allocatable :: label_list(:)
      real(dp) :: delta, wk0
      real(dp) :: dkx(3), dky(3)
      logical, external  :: matches
      logical            :: tend,terr
      logical            :: kband = .false.
      logical            :: kband_plane = .false.
      !
      !
      if ( tkpoints ) then
         call errore( ' card_kpoints ', ' two occurrences', 2 )
      endif
      !
      if ( matches( "automatic", input_line ) ) then
         !  automatic generation of k-points
         k_points = 'automatic'
      elseif ( matches( "crystal", input_line ) ) then
         !  input k-points are in crystal (reciprocal lattice) axis
         k_points = 'crystal'
         if ( matches( "_b", input_line ) ) kband=.true.
         if ( matches( "_c", input_line ) ) kband_plane=.true.
      elseif ( matches( "tpiba", input_line ) ) then
         !  input k-points are in 2pi/a units
         k_points = 'tpiba'
         if ( matches( "_b", input_line ) ) kband=.true.
         if ( matches( "_c", input_line ) ) kband_plane=.true.
      elseif ( matches( "gamma", input_line ) ) then
         !  only gamma (k=0) is used
         k_points = 'gamma'
      else
         !  by default, input k-points are in 2pi/a units
         k_points = 'tpiba'
      endif
      !
      if ( k_points == 'automatic' ) then
         !
         ! ... automatic generation of k-points
         !
         nkstot = 0
         call read_line( input_line, end_of_file = tend, error = terr )
         if (tend) goto 10
         if (terr) goto 20
         read(input_line, *, end=10, err=20) nk1, nk2, nk3, k1, k2 ,k3
         if ( k1 < 0 .or. k1 > 1 .or. &
               k2 < 0 .or. k2 > 1 .or. &
               k3 < 0 .or. k3 > 1 ) call errore &
                  ('card_kpoints', 'invalid offsets: must be 0 or 1', 1)
         if ( nk1 <= 0 .or. nk2 <= 0 .or. nk3 <= 0 ) call errore &
                  ('card_kpoints', 'invalid values for nk1, nk2, nk3', 1)
         allocate ( xk(3,1), wk(1) ) ! prevents problems with debug flags
         !                           ! when init_startk is called in iosys
      elseif ( ( k_points == 'tpiba' ) .or. ( k_points == 'crystal' ) ) then
         !
         ! ... input k-points 
         !
         call read_line( input_line, end_of_file = tend, error = terr )
         if (tend) goto 10
         if (terr) goto 20
         read(input_line, *, end=10, err=20) nkstot
         !
         if (kband) then
!
!        only the initial and final k points of the lines are given
!
            nkaux=nkstot
            allocate(xkaux(3,nkstot), wkaux(nkstot))
            allocate ( letter(nkstot) )
            allocate ( label_list(nkstot) )
            npk_label=0
            do i = 1, nkstot
               call read_line( input_line, end_of_file = tend, error = terr )
               if (tend) goto 10
               if (terr) goto 20
               do j=1,256   ! loop over all characters of input_line
                  if ((ichar(input_line(j:j)) < 58 .and. &   ! a digit
                       ichar(input_line(j:j)) > 47) &
                 .or. ichar(input_line(j:j)) == 43 .or. &    ! the + sign
                      ichar(input_line(j:j))== 45 .or. &     ! the - sign
                      ichar(input_line(j:j))== 46 ) then     ! a dot .
!
!   this is a digit, therefore this line contains the coordinates of the
!   k point. we read it and exit from the loop on the characters
!
                     read(input_line,*, end=10, err=20) xkaux(1,i), &
                                           xkaux(2,i), xkaux(3,i), wk0
                     wkaux(i) = nint ( wk0 ) ! beware: wkaux is integer
                     exit
                  elseif ((ichar(input_line(j:j)) < 123 .and. &
                           ichar(input_line(j:j)) > 64))  then
!
!   this is a letter, not a space character. we read the next three 
!   characters and save them in the letter array, save also which k point
!   it is
!
                     npk_label=npk_label+1
                     read(input_line(j:),'(a3)') letter(npk_label)
                     label_list(npk_label)=i
!
!  now we remove the letters from input_line and read the number of points
!  of the line. the next two line should account for the case in which
!  there is only one space between the letter and the number of points.
!
                     nch=3
                     if ( ichar(input_line(j+1:j+1))==32 .or. &
                          ichar(input_line(j+2:j+2))==32 ) nch=2
                     buffer=input_line(j+nch:)
                     read(buffer,*,err=20) wkaux(i)
                     exit
                  endif
               enddo
            enddo
            if ( npk_label > 0 ) then
               call cell_base_init ( ibrav, celldm, a, b, c, cosab, &
                              cosac, cosbc, trd_ht, rd_ht, cell_units )
               call transform_label_coord(ibrav, celldm_cb, xkaux, letter, &
                    label_list, npk_label, nkstot, k_points, point_label_type )
            end if

            deallocate(letter)
            deallocate(label_list)
            ! count k-points first
            nkstot=sum(wkaux(1:nkaux-1))+1
            do i=1,nkaux-1
              if (wkaux(i)==0) nkstot=nkstot+1
            enddo
            allocate ( xk(3,nkstot), wk(nkstot) )
            !
            !  generate the points along the lines
            !
            call generate_k_along_lines(nkaux, xkaux, wkaux, xk, wk, nkstot)
            !
            !  workaround: discard current wk (contains the length of k-path, 
            !  never used), replace with wk=1 so that band occupations (wg)
            !  are correctly written to file - needed by berkeleygw interface
            !
            wk(:) = 1.0_dp
            deallocate(xkaux)
            deallocate(wkaux)
            !
         elseif (kband_plane) then
!
!        generate a uniform mesh of k points on the plane defined by
!        the origin k_0, and two vectors applied in k_0, k_1 and k_2.
!
            if (nkstot /= 3) call errore ('card_kpoints', &
                                'option _c requires 3 k points',i)
            nkaux=nkstot
            allocate(xkaux(3,nkstot), wkaux(nkstot))
            do i = 1, nkstot
               call read_line( input_line, end_of_file = tend, error = terr )
               if (tend) goto 10
               if (terr) goto 20
               read(input_line,*, end=10, err=20) xkaux(1,i), xkaux(2,i), &
                                                  xkaux(3,i), wk0
               wkaux(i) = nint ( wk0 ) ! beware: wkaux is integer
            enddo
            ! count k-points first
            nkstot = wkaux(2) * wkaux(3)
            allocate ( xk(3,nkstot), wk(nkstot) )
            call generate_k_in_plane(nkaux, xkaux, wkaux, xk, wk, nkstot)
            deallocate(xkaux)
            deallocate(wkaux)
         else
!
!    reads on input the k points
!
            allocate ( xk(3, nkstot), wk(nkstot) )
            do i = 1, nkstot
               call read_line( input_line, end_of_file = tend, error = terr )
               if (tend) goto 10
               if (terr) goto 20
               read(input_line,*, end=10, err=20) xk(1,i),xk(2,i),xk(3,i),wk(i)
            enddo
         endif
         !
      elseif ( k_points == 'gamma' ) then
         !
         nkstot = 1
         allocate ( xk(3,1), wk(1) )
         xk(:,1) = 0.0_dp
         wk(1) = 1.0_dp
         !
      endif
      !
      tkpoints  = .true.
      tk_inp = .true.
      !
      return
10     call errore ('card_kpoints', ' end of file while reading ' &
            & // trim(k_points) // ' k points', 1)
20     call errore ('card_kpoints', ' error while reading ' &
            & // trim(k_points) // ' k points', 1)
      !
   end subroutine card_kpoints

   subroutine card_add_kpoints( input_line )
     use additional_kpoints, only : nkstot_add, xk_add
     implicit none
     character(len=*),intent(in) :: input_line
     character(len=256) :: input_line_aux
     real(dp),allocatable :: xk_old(:,:), wk_old(:)
     integer :: nk1_old, nk2_old, nk3_old, nkstot_old
     integer :: k1_old,  k2_old,  k3_old
     logical, external  :: matches
     !
     if(.not.allocated(xk) .or. .not.allocated(wk))&
       call errore("add_kpoints", "additional_k_points must appear after k_points",1)
     if(.not.tkpoints) &
       call errore("add_kpoints", "additional_k_points must appear after k_points",2)
     if(matches( "automatic", input_line )) &
       call errore("add_kpoints", "additional_k_points cannot be 'automatic'", 3)

     ! back-up existing points
     nkstot_old = nkstot
     allocate(xk_old(3,nkstot_old))
     allocate(wk_old(nkstot_old))
     xk_old  = xk
     wk_old  = wk
     nk1_old = nk1
     nk2_old = nk2
     nk3_old = nk3
     k1_old  = k1
     k2_old  = k2
     k3_old  = k3
     deallocate(xk,wk)

     ! prepare to read k-points again
     nkstot = 0
     input_line_aux = trim(adjustl(input_line))
     input_line_aux = input_line_aux(12:)
     tkpoints = .false.
     call card_kpoints(input_line_aux)
     !
     ! backup new points to module
     nkstot_add = nkstot
     if(nkstot_add==0) call errore("add_kpoints", "no new k_points?",1)
     allocate(xk_add(3,nkstot_add))
     xk_add = xk

     ! put back previous stuff
     deallocate(xk, wk)
     nkstot = nkstot_old
     allocate(xk(3,nkstot))
     allocate(wk(nkstot))
     xk  = xk_old
     wk  = wk_old
     nk1 = nk1_old
     nk2 = nk2_old
     nk3 = nk3_old
     k1  = k1_old
     k2  = k2_old
     k3  = k3_old
     deallocate(xk_old,wk_old)

     return 
   end subroutine card_add_kpoints
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! occupations
   !
   !   use the specified occupation numbers for electronic states.
   !   note that you should specify 10 values per line maximum!
   !
   ! syntax (nspin == 1):
   !
   !   occupations
   !      f(1)  ....   ....  f(10)
   !      f(11) .... f(nbnd)
   !
   ! syntax (nspin == 2):
   !
   !   occupations
   !      u(1)  ....   ....  u(10)
   !      u(11) .... u(nbnd)
   !      d(1)  ....   ....  d(10)
   !      d(11) .... d(nbnd)
   !
   ! example:
   !
   ! occupations
   !  2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
   !  2.0 2.0 2.0 2.0 2.0 1.0 1.0
   !
   ! where:
   !
   !      f(:) (real)  these are the occupation numbers
   !                   for lda electronic states.
   !
   !      u(:) (real)  these are the occupation numbers
   !                   for lsd spin == 1 electronic states
   !      d(:) (real)  these are the occupation numbers
   !                   for lsd spin == 2 electronic states
   !
   !      note, maximum 10 values per line!
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_occupations( input_line )
      !
      use wrappers, only: feval_infix
      !
      implicit none
      !
      character(len=256) :: input_line, field_str
      integer            :: is, nx10, i, j, nspin0
      integer            :: nfield, nbnd_read, nf, ierr
      logical :: tef
      !
      !
      if ( tocc ) then
         call errore( ' card_occupations ', ' two occurrences', 2 )
      endif
      nspin0=nspin
      if (nspin == 4) nspin0=1
      !
      allocate ( f_inp ( nbnd, nspin0 ) )
      do is = 1, nspin0
         !
         nbnd_read = 0
         do while ( nbnd_read < nbnd)
            call read_line( input_line, end_of_file=tef )
            if (tef) call errore('card_occupations',&
                        'missing occupations, end of file reached',1)
            call field_count( nfield, input_line )
            !
            do nf = 1,nfield
               nbnd_read = nbnd_read+1
               if (nbnd_read > nbnd ) exit
               call get_field(nf, field_str, input_line)
               !
               f_inp(nbnd_read,is) = feval_infix(ierr, field_str )
               call errore('card_occupations',&
                  'error parsing occupation: '//trim(field_str), nbnd_read*ierr)
            enddo
         enddo
         !
      enddo
      !
      tf_inp = .true.
      tocc = .true.
      !
      return
      !
   end subroutine card_occupations
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! cell_parameters
   !
   !   use the specified cell dimensions
   !
   ! syntax:
   !
   !    cell_parameters (cell_option)
   !      ht(1,1) ht(1,2) ht(1,3)
   !      ht(2,1) ht(2,2) ht(2,3)
   !      ht(3,1) ht(3,2) ht(3,3)
   !
   !   cell_option == alat      lattice vectors in units of alat
   !   cell_option == bohr      lattice vectors in bohr
   !   cell_option == angstrom  lattice vectors in angstrom
   !
   ! example:
   !
   ! cell_parameters
   !    24.50644311    0.00004215   -0.14717844
   !    -0.00211522    8.12850030    1.70624903
   !     0.16447787    0.74511792   23.07395418
   !
   ! where:
   !
   !      ht(i,j) (real)  cell dimensions ( in a.u. ),
   !                      note the relation with lattice vectors:
   !                      ht(1,:) = a1, ht(2,:) = a2, ht(3,:) = a3
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_cell_parameters( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      integer            :: i, j
      logical, external  :: matches
      !
      !
      if ( tcell ) then
         call errore( ' card_cell_parameters ', ' two occurrences', 2 )
      endif
      !
      if ( matches( "bohr", input_line ) ) then
         cell_units = 'bohr'
      elseif ( matches( "angstrom", input_line ) ) then
         cell_units = 'angstrom'
      elseif ( matches( "alat", input_line ) ) then
         cell_units = 'alat'
      else
         cell_units = 'none'
         call infomsg( 'read_cards ', &
            & 'deprecated: no units specified in cell_parameters card' )
         ! cell parameters are set in cell_base_init
      endif
      !
      do i = 1, 3
         call read_line( input_line )
         read(input_line,*) ( rd_ht( i, j ), j = 1, 3 )
      enddo
      !
      trd_ht = .true.
      tcell  = .true.
      !
      return
      !
   end subroutine card_cell_parameters
   !
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! ref_cell_parameters
   !
   !   use the specified cell dimensions
   !
   ! syntax:
   !
   !    ref_cell_parameters (cell_option)
   !      rd_ref_ht(1,1) rd_ref_ht(1,2) rd_ref_ht(1,3)
   !      rd_ref_ht(2,1) rd_ref_ht(2,2) rd_ref_ht(2,3)
   !      rd_ref_ht(3,1) rd_ref_ht(3,2) rd_ref_ht(3,3)
   !
   !   cell_option == alat      lattice vectors in units of alat set by ref_alat keyword (default)
   !   cell_option == bohr      lattice vectors in bohr
   !   cell_option == angstrom  lattice vectors in angstrom
   !
   ! example:
   !
   ! ref_cell_parameters
   !    24.50644311    0.00004215   -0.14717844
   !    -0.00211522    8.12850030    1.70624903
   !     0.16447787    0.74511792   23.07395418
   !
   ! where:
   !
   !      rd_ref_ht(i,j) (real)  cell dimensions ( in a.u. ),
   !        note the relation with reference lattice vectors:
   !        rd_ref_ht(1,:) = ref_a1, rd_ref_ht(2,:) = ref_a2, rd_ref_ht(3,:) = re_a3
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_ref_cell_parameters( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      integer            :: i, j
      logical, external  :: matches
      !
      !
      if ( ref_cell ) then
         call errore( ' card_reference_cell_parameters ', ' two occurrences', 2 )
      endif
      !
      if ( matches( "bohr", input_line ) ) then
         ref_cell_units = 'bohr'
      elseif ( matches( "angstrom", input_line ) ) then
         ref_cell_units = 'angstrom'
      else
         ref_cell_units = 'alat'
      endif
      !
      do i = 1, 3
         call read_line( input_line )
         read(input_line,*) ( rd_ref_ht( i, j ), j = 1, 3 )
      enddo
      !
      ref_cell = .true.
      !
      return
      !
   end subroutine card_ref_cell_parameters
   !
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! atomic_velocities
   !
   !   read velocities (in atomic units) from standard input
   !
   ! syntax:
   !
   !   atomic_velocities
   !     label(1)  vx(1) vy(1) vz(1)
   !     ....
   !     label(n)  vx(n) vy(n) vz(n)
   !
   ! example:
   !
   !   ???
   !
   ! where:
   !
   !   label (character(len=4))       atomic label
   !   vx(:), vy(:) and vz(:) (real)  x, y and z velocity components of
   !                                  the ions
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_ion_velocities( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      integer            :: ia, k, is, nfield
      character(len=4)   :: lb_vel
      !
      !
      if( tionvel ) then
         call errore( ' card_ion_velocities ', ' two occurrences', 2 )
      endif
      !
      if( .not. tapos ) then
         call errore( ' card_ion_velocities ', &
                     & ' atomic_species must be present before ', 2 )
      endif
      !
      rd_vel = 0.0_dp
      sp_vel = 0
      !
      if ( ion_velocities == 'from_input' ) then
         !
         tavel = .true.
         !
         do ia = 1, nat
            !
            call read_line( input_line )
            call field_count( nfield, input_line )
            if ( nfield == 4 ) then
               read(input_line,*) lb_vel, ( rd_vel(k,ia), k = 1, 3 )
            else
               call errore( ' iosys ', &
                           & ' wrong entries in atomic_velocities ', ia )
            endif
            !
            match_label: do is = 1, ntyp
               if ( trim( lb_vel ) == atom_label(is) ) then
                  sp_vel(ia) = is
                  exit match_label
               endif
            enddo match_label
            !
            if ( sp_vel(ia) < 1 .or. sp_vel(ia) > ntyp ) then
               call errore( ' iosys ', ' wrong label in ion_velocities ', ia )
            endif
            !
         enddo
         !
      endif
      !
      tionvel = .true.
      !
      return
      !
   end subroutine
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! constraints
   !
   !   ionic constraints
   !
   ! syntax:
   !
   !    constraints
   !      nconstr constr_tol
   !      constr_type(.) constr(1,.) constr(2,.) ... { constr_target(.) }
   !
   ! where:
   !
   !      nconstr(integer)    number of constraints
   !
   !      constr_tol          tolerance for keeping the constraints
   !                          satisfied
   !
   !      constr_type(.)      type of constrain:
   !                          1: for fixed distances ( two atom indexes must
   !                             be specified )
   !                          2: for fixed planar angles ( three atom indexes
   !                             must be specified )
   !
   !      constr(1,.) constr(2,.) ...
   !
   !                          indices object of the constraint, as
   !                          they appear in the 'position' card
   !
   !      constr_target       target for the constrain ( in the case of
   !                          planar angles it is the cos of the angle ).
   !                          this variable is optional.
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_constraints( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      integer            :: i, nfield
      !
      !
      if ( tconstr ) call errore( 'card_constraints', 'two occurrences', 2 )
      !
      call read_line( input_line )
      !
      call field_count( nfield, input_line )
      !
      if ( nfield == 1 ) then
         !
         read( input_line, * ) nconstr_inp
         !
      elseif ( nfield == 2 ) then
         !
         read( input_line, * ) nconstr_inp, constr_tol_inp
         !
      else
         !
         call errore( 'card_constraints', 'too many fields', nfield )
         !
      endif
      write(stdout,'(5x,a,i4,a,f12.6)') &
         'reading',nconstr_inp,' constraints; tolerance:', constr_tol_inp
      !
      call allocate_input_constr()
      !
      do i = 1, nconstr_inp
         !
         call read_line( input_line )
         !
         read( input_line, * ) constr_type_inp(i)
         !
         call field_count( nfield, input_line )
         !
         if ( nfield > nc_fields + 2 ) &
            call errore( 'card_constraints', &
                        'too many fields for this constraint', i )
         !
         select case( constr_type_inp(i) )
         case( 'type_coord', 'atom_coord' )
            !
            if ( nfield == 5 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i)
               !
               write(stdout,'(7x,i3,a,i3,a,i2,a,2f12.6)') i, &
                  ') '//constr_type_inp(i)(1:4),int(constr_inp(1,i)) ,&
                  ' coordination wrt type:', int(constr_inp(2,i)), &
                  ' cutoff distance and smoothing:',  constr_inp(3:4,i)
            elseif ( nfield == 6 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               write(stdout,'(7x,i3,a,i3,a,i2,a,2f12.6,a,f12.6)') i, &
                  ') '//constr_type_inp(i)(1:4),int(constr_inp(1,i)) , &
                  ' coordination wrt type:', int(constr_inp(2,i)), &
                  ' cutoff distance and smoothing:',  constr_inp(3:4,i), &
                  '; target:', constr_target_inp(i)
            else
               !
               call errore( 'card_constraints', 'type_coord, ' // &
                           & 'atom_coord: wrong number of fields', nfield )
               !
            endif
            !
         case( 'distance' )
            !
            if ( nfield == 3 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i)
               !
               write(stdout,'(7x,i3,a,2i3)') &
                  i,') distance between atoms: ', int(constr_inp(1:2,i))
            elseif ( nfield == 4 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               write(stdout,'(7x,i3,a,2i3,a,f12.6)') i, &
                  ') distance between atoms: ', int(constr_inp(1:2,i)), &
                  '; target:',  constr_target_inp(i)
            else
               !
               call errore( 'card_constraints', &
                           & 'distance: wrong number of fields', nfield )
               !
            endif
            !
         case( 'planar_angle' )
            !
            if ( nfield == 4 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i)
               !
               write(stdout, '(7x,i3,a,3i3)') &
                  i,') planar angle between atoms: ', int(constr_inp(1:3,i))
            elseif ( nfield == 5 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               write(stdout, '(7x,i3,a,3i3,a,f12.6)') i, &
                  ') planar angle between atoms: ', int(constr_inp(1:3,i)), &
                  '; target:', constr_target_inp(i)
            else
               !
               call errore( 'card_constraints', &
                           & 'planar_angle: wrong number of fields', nfield )
               !
            endif
            !
         case( 'torsional_angle' )
            !
            if ( nfield == 5 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i)
               !
               write(stdout, '(7x,i3,a,4i3)') &
                  i,') torsional angle between atoms: ', int(constr_inp(1:4,i))
            elseif ( nfield == 6 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               write(stdout, '(7x,i3,a,4i3,a,f12.6)') i, &
                  ') torsional angle between atoms: ', int(constr_inp(1:4,i)),&
                  '; target:', constr_target_inp(i)
            else
               !
               call errore( 'card_constraints', &
                           & 'torsional_angle: wrong number of fields', nfield )
               !
            endif
            !
         case( 'bennett_proj' )
            !
            if ( nfield == 5 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i)
               !
               write(stdout, '(7x,i3,a,i3,a,3f12.6)') i, &
                  ') bennet projection of atom ', int(constr_inp(1,i)), &
                  ' along vector:', constr_inp(2:4,i)
            elseif ( nfield == 6 ) then
               !
               read( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               write(stdout, '(7x,i3,a,i3,a,3f12.6,a,f12.6)') i, &
                  ') bennet projection of atom ', int(constr_inp(1,i)), &
                  ' along vector:', constr_inp(2:4,i), &
                  '; target:', constr_target_inp(i)
            else
               !
               call errore( 'card_constraints', &
                           & 'bennett_proj: wrong number of fields', nfield )
               !
            endif
            !
         case default
            !
            call errore( 'card_constraints', 'unknown constraint ' // &
                        & 'type: ' // trim( constr_type_inp(i) ), 1 )
            !
         end select
         !
      enddo
      !
      tconstr = .true.
      !
      return
      !
   end subroutine card_constraints
   !
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! ksout
   !
   !   enable the printing of kohn sham states
   !
   ! syntax ( nspin == 2 ):
   !
   !   ksout
   !     nu
   !     iu(1) iu(2) iu(3) .. iu(nu)
   !     nd
   !     id(1) id(2) id(3) .. id(nd)
   !
   ! syntax ( nspin == 1 ):
   !
   !   ksout
   !     ns
   !     is(1) is(2) is(3) .. is(ns)
   !
   ! example:
   !
   !   ???
   !
   ! where:
   !
   !   nu (integer)     number of spin=1 states to be printed
   !   iu(:) (integer)  indexes of spin=1 states, the state iu(k)
   !                    is saved to file ks_up.iu(k)
   !
   !   nd (integer)     number of spin=2 states to be printed
   !   id(:) (integer)  indexes of spin=2 states, the state id(k)
   !                    is saved to file ks_dw.id(k)
   !
   !   ns (integer)     number of lda states to be printed
   !   is(:) (integer)  indexes of lda states, the state is(k)
   !                    is saved to file ks.is(k)
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_ksout( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      integer            :: i, s, nksx
      type occupancy_type
         integer, pointer :: occs(:)
      end type occupancy_type
      type(occupancy_type), allocatable :: is(:)
      !
      if ( tksout ) then
         call errore( ' card_ksout ', ' two occurrences', 2 )
      endif
      !
      nprnks = 0
      nksx   = 0
      !
      allocate ( is (nspin) )
      !
      do s = 1, nspin
         !
         call read_line( input_line )
         read(input_line, *) nprnks( s )
         !
         if ( nprnks( s ) < 1 ) then
            call errore( ' card_ksout ', ' wrong number of states ', 2 )
         endif
         !
         allocate( is(s)%occs( 1:nprnks(s) ) )
         !
         call read_line( input_line )
         read(input_line, *) ( is(s)%occs(i), i = 1, nprnks( s ) )
         !
         nksx = max( nksx, nprnks( s ) )
         !
      enddo
      !
      call allocate_input_iprnks( nksx, nspin )
      !
      do s = 1, nspin
         !
         do i = 1, nprnks( s )
            !
            iprnks( i, s ) = is(s)%occs(i)
            !
         enddo
         !
         deallocate( is(s)%occs )
         !
      enddo
      !
      deallocate( is )
      !
      tksout = .true.
      !
      return
      !
   end subroutine
   !
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   ! plot wannier
   !
   !   needed to specify the indices of the wannier functions that
   !   have to be plotted
   !
   ! syntax:
   !
   !   plot_wannier
   !     index1, ..., indexn
   !
   ! where:
   !
   !   index1, ..., indexn are indices of the wannier functions
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_plot_wannier( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      logical, external  :: matches
      !
      integer                    :: i, ib
      character(len=6)           :: i_char
      character(len=6), external :: int_to_char
      !
      !
      if ( twannier ) &
         call errore( 'card_plot_wannier', 'two occurrences', 2 )
      !
      if ( nwf > 0 ) then
         !
         if ( nwf > nwf_max ) &
            call errore( 'card_plot_wannier', 'too many wannier functions', 1 )
         !
         call read_line( input_line )
         !
         ib = 0
         !
         do i = 1, nwf_max
            !
            i_char = int_to_char( i )
            !
            if ( matches( ' ' // trim( i_char ) // ',', &
                           ' ' // trim( input_line ) // ',' ) ) then
               !
               ib = ib + 1
               !
               if ( ib > nwf ) &
                  call errore( 'card_plot_wannier', 'too many indices', 1 )
               !
               wannier_index(ib) = i
               !
            endif
            !
         enddo
         !
      endif
      !
      twannier = .true.
      !
      return
      !
   end subroutine card_plot_wannier
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !
   !
   ! template
   !
   !      this is a template card info section
   !
   ! syntax:
   !
   !    template
   !     rvalue ivalue
   !
   ! example:
   !
   !    ???
   !
   ! where:
   !
   !      rvalue (real)     this is a real value
   !      ivalue (integer)  this is an integer value
   !
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_template( input_line )
      !
      implicit none
      !
      character(len=256) :: input_line
      !
      !
      if ( ttemplate ) then
         call errore( ' card_template ', ' two occurrences', 2 )
      endif
      !
      ! ....  code here
      !
      ttemplate = .true.
      !
      return
      !
   end subroutine
   !
   !
   !------------------------------------------------------------------------
   !    begin manual
   !----------------------------------------------------------------------
   !wannier_ac
   !wannier# 1 10.5 15.7 2
   !atom 1
   !d 1 0.45
   !p 3 0.55
   !wannier# 2 10.5 15.7 1
   !atom 3
   !p 1 0.8
   !spin#2:
   !wannier# 1 10.5 15.7 2
   !atom 1
   !d 1 0.45
   !p 3 0.55
   !wannier# 2 10.5 15.7 1
   !atom 3
   !p 1 0.8
   !----------------------------------------------------------------------
   !    end manual
   !------------------------------------------------------------------------
   !
   subroutine card_wannier_ac( input_line )
      !
      use wannier_new, only: nwan

      implicit none
      !
      character(len=256) :: input_line
      integer :: i,j,k, nfield, iwan, ning, iatom,il,im,ispin
      logical :: tend
      real :: c, b_from, b_to
      character(len=10) :: text, lo

      ispin = 1
      !
      do i = 1, nwan
         !
         call read_line( input_line, end_of_file = tend )
         !
         if ( tend ) &
            call errore( 'read_cards', &
                        'end of file reading trial wfc composition', i )
         !
         call field_count( nfield, input_line )
         !
         if ( nfield == 4 ) then
            read(input_line,*) text, iwan, b_from, b_to
            ning = 1
         elseif ( nfield == 5 ) then
            read(input_line,*) text, iwan, b_from, b_to, ning
         else
            call errore( 'read_cards', &
                        'wrong format', nfield )
         endif
         if(iwan/=i) call errore( 'read_cards', 'wrong wannier order', iwan)

         ! read atom number
         call read_line( input_line, end_of_file = tend )
         read(input_line,*) text, iatom
         !
         wan_data(iwan,ispin)%iatom = iatom
         wan_data(iwan,ispin)%ning = ning
         wan_data(iwan,ispin)%bands_from = b_from
         wan_data(iwan,ispin)%bands_to = b_to
         !
         do j=1, ning
            call read_line( input_line, end_of_file = tend )
            !
            if ( tend ) &
               call errore( 'read_cards', &
                           'not enough wavefunctions', j )
            if (ning==1) then
               read(input_line,*) lo,im
               c = 1.d0
            else
               read(input_line,*) lo,im,c
            endif

            select case(trim(lo))
            case('s')
               il = 0
            case('p')
               il = 1
            case('d')
               il = 2
            case('f')
               il = 3
            case default
               call errore( 'read_cards', &
                           'wrong l-label', 1 )
            end select

            wan_data(iwan,ispin)%ing(j)%l = il
            wan_data(iwan,ispin)%ing(j)%m = im
            wan_data(iwan,ispin)%ing(j)%c = c
         enddo
      enddo

      !is there spin 2 information?
      call read_line( input_line, end_of_file = tend )
      !
      if ( .not. tend ) then
         read(input_line,*) text
         if ( trim(text) == 'spin#2:') then ! ok, there is spin 2 data
            ispin = 2
            !
            do i = 1, nwan
               !
               call read_line( input_line, end_of_file = tend )
               !
               if ( tend ) &
                  call errore( 'read_cards', &
                              'end of file reading trial wfc composition', i )
               !
               call field_count( nfield, input_line )
               !
               if ( nfield == 4 ) then
                  read(input_line,*) text, iwan, b_from, b_to
                  ning = 1
               elseif ( nfield == 5 ) then
                  read(input_line,*) text, iwan, b_from, b_to, ning
               else
                  call errore( 'read_cards', &
                              'wrong format', nfield )
               endif
               if(iwan/=i) call errore( 'read_cards', 'wrong wannier order', iwan)

               ! read atom number
               call read_line( input_line, end_of_file = tend )
               read(input_line,*) text, iatom
               !
               wan_data(iwan,ispin)%iatom = iatom
               wan_data(iwan,ispin)%ning = ning
               wan_data(iwan,ispin)%bands_from = b_from
               wan_data(iwan,ispin)%bands_to = b_to
               !
               do j=1, ning
                  call read_line( input_line, end_of_file = tend )
                  !
                  if ( tend ) &
                     call errore( 'read_cards', &
                                 'not enough wavefunctions', j )
                  if (ning==1) then
                     read(input_line,*) lo,im
                     c = 1.d0
                  else
                     read(input_line,*) lo,im,c
                  endif

                  select case(trim(lo))
                  case('s')
                     il = 0
                  case('p')
                     il = 1
                  case('d')
                     il = 2
                  case('f')
                     il = 3
                  case default
                     call errore( 'read_cards', &
                                 'wrong l-label', 1 )
                  end select

                  wan_data(iwan,ispin)%ing(j)%l = il
                  wan_data(iwan,ispin)%ing(j)%m = im
                  wan_data(iwan,ispin)%ing(j)%c = c
               enddo
            enddo
         else
         ! oups - that is not our data - lets's move one line up in input file
         ! not sure that a direct access to the parce_unit is safe enougth
         backspace(parse_unit)
         endif
      else
         ! ok, that's the end of file. but i will move one line up
         ! for a correct handling of eof in the parent read_cards subroutine
         ! otherwise (at least with gfortran on mac) there will be the read error
         backspace(parse_unit)
      endif
      !
      return
      !
   end subroutine card_wannier_ac
end module read_cards_module
