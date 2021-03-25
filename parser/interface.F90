!
!     Shows FDF capabilities..
!
PROGRAM iochamp
  USE fdf
  USE prec


! Note the following two modules are being used to store and process the parsed data  
  use keywords
  use periodic_table
!  
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug

  character(len=72)          :: fname, filename, fmt
  character(len=72)          :: molecule_name, key, comment
  character(2)               :: symbol(maxa)
  character(len=20)          :: chunks(10), subblock(10)
  character(len=30)          :: keyword(5) 
  integer(sp)                :: i, j, ia, na, number_of_atoms
  integer(sp)                :: isa(maxa)
  real(dp)                   :: coeff(maxa)
  real(sp)                   :: wmix
  real(dp)                   :: phonon_energy
  real(dp)                   :: xa(3, maxa)
  real(dp)                   :: listr(maxa)
  type(block_fdf)            :: bfdf, bfdf2
  type(parsed_line), pointer :: pline, pline2
  real(dp)                   :: float_value
  character(len=20)          :: real_format    = '(A, T20, F14.8)'
  character(len=20)          :: int_format     = '(A, T20, I8)'
  character(len=80)          :: string_format  = '(A, T40, A)'  
  character(len=80)          :: logical_format = '(A, T40, L)'    

! for determinants sections
  integer                    :: nelectrons, nexcitation, iostat
  integer, allocatable       :: det_alpha(:), det_beta(:)
  real(selected_real_kind(6,15)), allocatable :: det_coeff(:)
  character(len=20)          :: temp1, temp2, temp3, temp4, temp5
!------------------------------------------------------------------------- BEGIN

! Initialize
  call fdf_init('test.inp', 'test.out')


  write(6,*) '------------------------------------------------------'  
! strings/characters
  title = fdf_string('title', 'Default title')
  write(6,'(A)') 'Title of the calculation :: ', title

! Get the directory where the pooled data is kept
  path_pool = fdf_string('pool', './')
  write(6,fmt=string_format) 'pool directory location :: ', path_pool

  write(6,*) '------------------------------------------------------'  


! Get all the filenames from which the data is to be read
  file_basis = fdf_load_filename('basis', 'default.gbs')
  write(6,fmt=string_format) 'filename basis :: ', trim(file_basis)

  file_determinants = fdf_load_filename('determinants', 'default.det')
  write(6,fmt=string_format) 'filename determinants :: ', trim(file_determinants)


  write(6,*) '------------------------------------------------------'  

! Logical variables
  optimize_wave = fdf_boolean("optimize_wave", .false.)
  write(6,*) ' optimize_wavefunction = ', optimize_wave

! Integer numbers (keyword, default_value). The variable is assigned default_value when keyword is not present

  ncore = fdf_integer('ncore', 0)
  write(6,fmt=int_format) 'NCore =', ncore

! floats (keyword, default_value) variable is assigned default_value when keyword is not present

  sr_eps = fdf_double('sr_eps', 0.025d0)
  write(6,fmt=real_format) 'sr_eps:', sr_eps

  ! logical :: true, .true., yes, T, 1, and TRUE are equivalent
  debug = fdf_boolean('Debug', .TRUE.)
  write(6,'(A, L2)') 'Debug:', debug

! floats/integers/strings/boolean can be parsed generically using fdf_get

  sr_tau = fdf_get('sr_tau', 0.025d0)
  write(6,fmt=real_format) 'sr_tau:', sr_tau

  nspin1 = fdf_get('nspin1', 1)
  write(6,fmt=int_format) 'nspin1 from global.fdf ', nspin1

  energy_tol = fdf_get('energy_tol', 0.00001d0)
  write(6,fmt=real_format) 'energy_tol:', energy_tol

  opt_method = fdf_get('opt_method', "sr_n")
  write(6,fmt=string_format) 'Optimization method ', opt_method

  multiple_adiag = fdf_get('multiple_adiag', .false.)
  write(6,fmt=logical_format) 'multiple_adiag:', multiple_adiag


  ! mixed types in one line (for example, reading a number with units)
  tau = fdf_get('tau', 0.05)
  write(6,fmt=real_format) 'DMC tau = ', tau

  etrial = fdf_physical('etrial', -20.d0, 'eV')
  write(6,fmt=real_format) 'Energy CutOff in eV :: ', energy_trial

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'



!  write(6,'(A,4X)') 'optimize_wavefunction using bline', (subblock(i), i = 1, 4)

  if (fdf_block('general', bfdf)) then
    write(*,*) "inside general block"
    i = 1
    do while(fdf_bline(bfdf, pline))    
      doit = fdf_bsearch(pline, "pool")    
      write(*,*) "pool found", doit      
      i = i + 1
    enddo
  endif

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'
  
  
  

  if (.not. fdf_block('molecule', bfdf)) then
      !   External file reading
          write(6,*) 'Reading coordinates of the molecule from an external file'
          ia = 1

          open (unit=12,file=file_molecule, iostat=iostat, action='read' )
          if (iostat .ne. 0) stop "Problem in opening the molecule file"
          read(12,*) natoms
          print*, "natoms ", natoms
          if (.not. allocated(cent)) allocate(cent(3,natoms))
          
          read(12,'(A)')  key
          print*, "Comment :: ", trim(key)
          do i = 1, natoms
            read(12,*) symbol(i), cent(1,i), cent(2,i), cent(3,i)
          enddo
          close(12)

          write(6,*) 'Coordinates from Molecule load construct: '
          do ia= 1, natoms
            write(6,'(A4,3F10.6)') symbol(ia), (cent(i,ia),i=1,3)
          enddo
  endif
 
  write(6,'(A)')  
  write(6,*) '------------------------------------------------------'




  if (fdf_block('molecule', bfdf)) then
    !   External file reading
        write(6,*) 'Beginning of external file coordinates block  '
        ia = 1
    
        do while((fdf_bline(bfdf, pline)))
!         get the integer from the first line 
          if ((pline%id(1) .eq. "i") .and. (pline%ntokens .eq. 1)) then        ! check if it is the only integer present in a line
            natoms = fdf_bintegers(pline, 1)
            write(*,*) "Number of atoms = ", natoms
          endif

          if (.not. allocated(cent)) allocate(cent(3,natoms))
        
          if (pline%ntokens == 4) then
            symbol(ia) = fdf_bnames(pline, 1)
            do i= 1, 3
              cent(i,ia) = fdf_bvalues(pline, i)
            enddo
            ia = ia + 1
          endif
        enddo

        write(6,*) 'Coordinates from single line Molecule block: '
        do ia= 1, natoms
          write(6,'(A4,3F10.6)') symbol(ia), (cent(i,ia),i=1,3)
        enddo
      endif

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'


!  Molecule coordinate block begins here  for demonstration

  if (fdf_block('Coordinates', bfdf)) then
    ia = 1
    do while(fdf_bline(bfdf, pline))
      symbol(ia) = fdf_bnames(pline, 1)
      do i= 1, 3
        xa(i,ia) = fdf_bvalues(pline, i)
      enddo
      ia = ia + 1
    enddo
    write(6,*) 'Coordinates from explicit data block:'
    do j = 1, ia
      write(6,'(A, 4x, 3F10.6)') symbol(j), (xa(i,j),i=1,3) 
    enddo
  endif

  write(6,*) '------------------------------------------------------'







  !  Determinants as a block. read directly from the input file
!    under construction
    if (fdf_block('determinants', bfdf)) then
      ia = 1
      do while(fdf_bline(bfdf, pline))
        symbol(ia) = fdf_bnames(pline, 1)
        do i= 1, 3
          xa(i,ia) = fdf_bvalues(pline, i)
        enddo
        ia = ia + 1
      enddo
      na = ia - 1 

    endif

    ! if (fdf_block('Coordinates', bfdf)) then
    !   write(6,*) 'Coordinates:'
    !   do ia = 1, na
    !     write(6,'(A, 4x, 3F10.6)') symbol(ia), (xa(i,ia),i=1,3) 
    !   enddo
    ! endif


    write(6,*) '------------------------------------------------------'



  if (.not. fdf_block('determinants', bfdf)) then
    if ( fdf_load_defined('determinants') ) then
      !   External file reading
        write(6,'(A)')  " Determinants Block"

        write(6,*) '------------------------------------------------------'      

        write(6,*) 'Reading the determinants block from an external file '

        open (unit=11,file=file_determinants, iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the determinant file"
        read(11,*) temp1, temp2, nelectrons, temp3, nalpha

        read(11,*)  temp1, ndeterminants, nexcitation
        if (.not. allocated(det_coeff)) allocate(det_coeff(ndeterminants))           

        read(11,*) (det_coeff(i), i=1,ndeterminants)
!        write(*,'(<ndeterminants>(f11.8, 1x))') (det_coeff(i), i=1,ndeterminants)    ! for Intel Fortran  

        nbeta       = nelectrons - nalpha        

        write(*,*) "total number of       electrons ", nelectrons
        write(*,*) "      number of alpha electrons ", nalpha        
        write(*,*) "      number of beta  electrons ", nbeta

        write(fmt,*)  '(', ndeterminants, '(f11.8,1x))'
        write(*,fmt) (det_coeff(i), i=1,ndeterminants)

!       allocate the orbital mapping array        
        if (.not. allocated(iworbd)) allocate(iworbd(nelectrons, ndeterminants))
        
        do i = 1, ndeterminants
          read(11,*) (iworbd(j,i), j=1,nelectrons)
        enddo

        write(fmt,*)  '(i4,1x)'        
        do i = 1, ndeterminants
          write(*,'(<nelectrons>(i4, 1x))') (iworbd(j,i), j=1,nelectrons)
        enddo
     
        read(11,*) temp1
        if (temp1 == "end" ) write(*,*) "Determinant File read successfully "
        close(11)

    endif ! condition if load determinant is present

  endif ! condition determinant block not present

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'


  call fdf_shutdown()

!----------------------------------------------------------------------------END
END PROGRAM iochamp
