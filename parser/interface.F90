!
!
PROGRAM interface
  USE fdf
  USE prec


! Note the following two modules are being used to store and process the parsed data  
  use keywords
  use periodic_table
!  
  implicit none
!--------------------------------------------------------------- Local Variables
  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline

  logical                    :: debug

  integer(sp)                :: i, j, ia
  
  character(len=72)          :: fmt, key
  character(2), allocatable  :: symbol(:)
  
  character(len=20)          :: real_format    = '(A, T28, F14.8)'
  character(len=20)          :: int_format     = '(A, T34, I8)'
  character(len=80)          :: string_format  = '(A, T40, A)'  
  character(len=80)          :: logical_format = '(A, T40, L)'    

! for determinants sections
  integer                    :: nelectrons, iostat
  real(kind=8), allocatable  :: det_coeff(:)
  character(len=20)          :: temp1, temp2, temp3
!------------------------------------------------------------------------- BEGIN

! Initialize
  call fdf_init('test.inp', 'test.out')


  write(6,*) '------------------------------------------------------'  
! strings/characters
  title = fdf_string('title', 'Default title')
  write(6,fmt=string_format) ' Title of the calculation :: ', title

! Get the directory where the pooled data is kept
  path_pool = fdf_string('pool', './')
  write(6,fmt=string_format) ' pool directory location :: ', path_pool

  write(6,*) '------------------------------------------------------'  




! Get all the filenames from which the data is to be read
  file_basis = fdf_load_filename('basis', 'default.gbs')
  write(6,fmt=string_format) ' filename basis :: ', trim(file_basis)

  file_molecule = fdf_load_filename('molecule', '')
  write(6,fmt=string_format) ' filename molecule :: ', trim(file_molecule)

  file_determinants = fdf_load_filename('determinants', 'default.det')
  write(6,fmt=string_format) ' filename determinants :: ', trim(file_determinants)


  write(6,*) '------------------------------------------------------'  

! Logical variables
  optimize_wave = fdf_boolean("optimize_wave", .false.)
!  write(6,fmt=logical_format) 'optimize_wavefunction = ', optimize_wave


! Integer numbers (keyword, default_value). The variable is assigned default_value when keyword is not present
  nextorb  = fdf_integer('nextorb', 0)
!  write(6,fmt=int_format) ' NExtOrb =', nextorb

! floats (keyword, default_value) variable is assigned default_value when keyword is not present
  sr_eps = fdf_double('sr_eps', 0.025d0)
!  write(6,fmt=real_format) ' sr_eps:', sr_eps

! logical :: true, .true., yes, T, 1, and TRUE are equivalent
  debug = fdf_boolean('Debug', .TRUE.)
!  write(6,'(A, L2)') ' Debug:', debug



! floats/integers/strings/boolean can be parsed generically using fdf_get
  sr_tau = fdf_get('sr_tau', 0.025d0)
!  write(6,fmt=real_format) ' sr_tau:', sr_tau

  nspin1 = fdf_get('nspin1', 1)
!  write(6,fmt=int_format) ' nspin1 from global ', nspin1

  energy_tol = fdf_get('energy_tol', 0.00001d0)
!  write(6,fmt=real_format) ' energy_tol:', energy_tol

  opt_method = fdf_get('opt_method', "sr_n")
!  write(6,fmt=string_format) ' Optimization method ', opt_method

  multiple_adiag = fdf_get('multiple_adiag', .false.)
!  write(6,fmt=logical_format) ' multiple_adiag:', multiple_adiag


! mixed types in one line (for example, reading a number with units)
  tau = fdf_get('tau', 0.05)
!  write(6,fmt=real_format) ' DMC tau = ', tau

  etrial = fdf_physical('etrial', -20.d0, 'eV')
!  write(6,fmt=real_format) ' Energy CutOff in eV :: ', energy_trial


! Pretty printing of above-mentioned keywords
  write(6,'(A)')  
  write(6,*) '------------------------------------------------------'

  write(6,fmt=string_format) ' Optimization method ', opt_method  

  write(6,fmt=logical_format) ' Optimize wavefunctions :: ', optimize_wave
  write(6,fmt=logical_format) ' multiple_adiag :: ', multiple_adiag
  write(6,fmt=logical_format) ' Debug :: ', debug

  write(6,*) '-------------------------'  

  write(6,fmt=int_format) ' NExtOrb :: ', nextorb
  write(6,fmt=int_format) ' Nspin1 from global :: ', nspin1

  write(6,*) '-------------------------'  

  write(6,fmt=real_format) ' sr_tau :: ', sr_tau
  write(6,fmt=real_format) ' energy_tol :: ', energy_tol

  write(6,*) '-------------------------'  

  write(6,fmt=real_format) ' Trial Energy in eV :: ', energy_trial

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
          if (.not. allocated(symbol)) allocate(symbol(natoms))                    
          
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
    write(6,'(A)')  
    write(6,*) '------------------------------------------------------'        
  endif
 




  if (fdf_block('molecule', bfdf)) then
    !   External file reading
        write(6,*) 'Beginning of molecular coordinates block  '
        ia = 1
    
        do while((fdf_bline(bfdf, pline)))
!         get the integer from the first line 
          if ((pline%id(1) .eq. "i") .and. (pline%ntokens .eq. 1)) then        ! check if it is the only integer present in a line
            natoms = fdf_bintegers(pline, 1)
            write(*,*) "Number of atoms = ", natoms
          endif

          if (.not. allocated(cent)) allocate(cent(3,natoms))
          if (.not. allocated(symbol)) allocate(symbol(natoms))          
        
          if (pline%ntokens == 4) then
            symbol(ia) = fdf_bnames(pline, 1)
            do i= 1, 3
              cent(i,ia) = fdf_bvalues(pline, i)
            enddo
            ia = ia + 1
          endif
        enddo

        write(6,*) 'Coordinates from Molecule block: '
        do ia= 1, natoms
          write(6,'(A4,3F10.6)') symbol(ia), (cent(i,ia),i=1,3)
        enddo

        write(6,'(A)')  
        write(6,*) '------------------------------------------------------'     
  endif




  if (.not. fdf_block('determinants', bfdf)) then
    if ( fdf_load_defined('determinants') ) then
      !   External file reading
        write(6,'(A)')  " Determinants Block"

        write(6,*) '------------------------------------------------------'      

        write(6,*) 'Reading the determinants block from an external file '

        open (unit=11,file=file_determinants, iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the determinant file"
        read(11,*) temp1, temp2, nelectrons, temp3, nalpha

        read(11,*)  temp1, ndeterminants, iwctype
        if (.not. allocated(det_coeff)) allocate(det_coeff(ndeterminants))           

        read(11,*) (det_coeff(i), i=1,ndeterminants)
!        write(*,'(<ndeterminants>(f11.8, 1x))') (det_coeff(i), i=1,ndeterminants)    ! for Intel Fortran  

        nbeta       = nelectrons - nalpha        

        write(*,*) "total number of       electrons ", nelectrons
        write(*,*) "      number of alpha electrons ", nalpha        
        write(*,*) "      number of beta  electrons ", nbeta

        write(6,'(A)')  
        write(*,*) "Determinant Coefficients" 
        write(fmt,*)  '(', ndeterminants, '(f11.8,1x))'
        write(*,fmt) (det_coeff(i), i=1,ndeterminants)

!       allocate the orbital mapping array        
        if (.not. allocated(iworbd)) allocate(iworbd(nelectrons, ndeterminants))
        
        do i = 1, ndeterminants
          read(11,*) (iworbd(j,i), j=1,nelectrons)
        enddo

        write(6,'(A)')  
        write(*,*) "Spin-alpha and Spin-beta determinants" 
        write(fmt,*) '(', nelectrons, '(i4,1x))'        
        do i = 1, ndeterminants
!          write(*,'(<nelectrons>(i4, 1x))') (iworbd(j,i), j=1,nelectrons)       ! For Intel Fortran
          write(*,fmt) (iworbd(j,i), j=1,nelectrons)          
        enddo
     
        read(11,*) temp1
        if (temp1 == "end" ) write(*,*) "Determinant File read successfully "
        close(11)

    endif ! condition if load determinant is present
    write(6,'(A)')  
    write(6,*) '------------------------------------------------------'
  endif ! condition determinant block not present



  call fdf_shutdown()

!----------------------------------------------------------------------------END
END PROGRAM interface
