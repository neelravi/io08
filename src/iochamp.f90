!
!     Shows FDF capabilities..
!
PROGRAM iochamp
  USE fdf
  USE prec
  USE parse
  use io_fdf

! Note the following two modules are being used to store and process the parsed data  
  use keywords
  use periodic_table
!  
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug, check, val, logic(10)
  logical                    :: beginning, ending
  character(len=72)          :: fname, axis, status, filename, fmt
  character(len=72)          :: molecule_name, key, comment
  character(2)               :: symbol(maxa)
  character(len=20)          :: chunks(10), subblock(10)
  character(len=30)          :: keyword(5) 
  integer(sp)                :: i, j, ia, na, external_entry, number_of_atoms, ind
  integer(sp)                :: isa(maxa)
  real(dp)                   :: coeff(maxa)
  real(sp)                   :: wmix
  real(dp)                   :: cutoff, phonon_energy, factor
  real(dp)                   :: xa(3, maxa)
  real(dp)                   :: listr(maxa)
  type(block_fdf)            :: bfdf, bfdf2
  type(parsed_line), pointer :: pline, pline2
  !type(fdf_file)             :: fdffile 
  integer                    ::  max_iteration, max_iter, linecount, argument(5)
  real(dp)                   :: float_value
  character(len=20)          :: real_format = '(A, T20, F14.8)'
  character(len=20)          :: int_format = '(A, T20, I8)'

! for determinants sections
  integer                    :: nelectrons, nexcitation, iostat
  integer, allocatable       :: det_alpha(:), det_beta(:)
  real(selected_real_kind(6,15)), allocatable :: det_coeff(:)
  character(len=20)          :: temp1, temp2, temp3, temp4, temp5
!------------------------------------------------------------------------- BEGIN

! Initialize
  call fdf_init('test-champ.inp', 'test-champ.out')

! strings/characters
  title = fdf_string('title', 'Default title')
  write(6,'(A)') 'Title of the calculation :: ', title

  ! Get the directory where the pooled data is kept
  path_pool = fdf_string('pool', './')
  write(6,'(A)') 'pool directory location :: ', path_pool

  ! Get all the filenames from which the data is to be read
  file_pseudo = fdf_load_filename('pseudopot', 'default.psp')
  write(6,'(A)') 'filename pseuodpotential :: ', file_pseudo

  file_basis = fdf_load_filename('basis', 'default.bas')
  write(6,'(A)') 'filename basis :: ', file_basis

  file_determinants = fdf_load_filename('determinants', 'default.det')
  write(6,'(A)') 'filename determinants :: ', file_determinants

  file_orbitals = fdf_load_filename('orbitals', 'default.orb')
  write(6,'(A)') 'filename orbitals :: ', file_orbitals

  file_jastrow = fdf_load_filename('jastrow', 'default.jas')
  write(6,'(A)') 'filename jastrow :: ', file_jastrow

  file_jastrow_deriv = fdf_load_filename('jastrow_deriv', 'default.jasder')
  write(6,'(A)') 'filename jastrow derivatives :: ', file_jastrow_deriv


! &optwf ioptwf 1 ioptci 1 ioptjas 1 ioptorb 1
  optimize_wavefunction = fdf_boolean("optimize_wavefunction", .false.)
  write(6,*) ' optimize_wavefunction = ', optimize_wavefunction

  optimize_ci = fdf_boolean('optimize_ci', .false.)
  write(6,*) ' optimize_ci = ', optimize_ci

  optimize_jastrow = fdf_boolean('optimize_jastrow', .false.)
  write(6,*) ' optimize_jastrow = ', optimize_jastrow

  optimize_orbitals = fdf_boolean('optimize_orbitals', .false.)
  write(6,*) ' optimize_orbitals = ', optimize_orbitals

  write(6,'(A)')  
  write(6,*) '------------------------------------------------------'


!Integer numbers (keyword, default_value). The variable is assigned default_value when keyword is not present
  ! &optwf ncore 0 nextorb 280 no_active 0
  ! &optwf nblk_max 200 nopt_iter 2
  ncore = fdf_integer('ncore', 0)
  write(6,fmt=int_format) 'NCore =', ncore

  nextorb = fdf_integer('nextorb', 0)
  write(6,fmt=int_format) 'Next Orb =', nextorb

  no_active = fdf_integer('no_active', 0)
  write(6,fmt=int_format) 'no_active =', no_active

  nblk_max = fdf_integer('nblk_max', 0)
  write(6,fmt=int_format) 'nblk max =', nblk_max

  nopt_iter = fdf_integer('nopt_iter', 0)
  write(6,fmt=int_format) 'nopt_iter =', nopt_iter


! floats (keyword, default_value) variable is assigned default_value when keyword is not present

  ! &optwf sr_tau 0.025 sr_eps 0.001 sr_adiag 0.01
  ! &optwf isample_cmat 0 energy_tol 0.0
  
  sr_tau = fdf_get('sr_tau', 0.025d0)
  write(6,fmt=real_format) 'sr_tau:', sr_tau

  sr_eps = fdf_get('sr_eps', 0.001d0)
  write(6,fmt=real_format) 'sr_eps:', sr_eps

  sr_adiag = fdf_get('sr_adiag', 0.01d0)
  write(6,fmt=real_format) 'sr_adiag:', sr_adiag

  energy_tol = fdf_get('energy_tol', 0.00001d0)
  write(6,fmt=real_format) 'energy_tol:', energy_tol

  ! &optwf method sr_n multiple_adiag 0

  opt_method = fdf_get('opt_method', "sr_n")
  write(6,*) 'Optimization method ', opt_method

  multiple_adiag = fdf_get('multiple_adiag', .false.)
  write(6,*) 'multiple_adiag:', multiple_adiag





  ! logical :: true, .true., yes, T, and TRUE are equivalent
  debug = fdf_boolean('Debug', .TRUE.)
  write(6,'(A, L2)') 'Debug:', debug




! ianalyt_lap 1 isc 2 nspin1 1 nspin2 1 ifock 0
  analytic_laplacian = fdf_get('ianalyt_lap', 1)
  write(6,*) 'analytic laplacian from global.fdf pointer explained ', ianalyt_lap

  nspin1 = fdf_get('nspin1', 1)
  write(6,*) 'nspin1 from global.fdf ', nspin1

  nspin2 = fdf_get('nspin2', 1)
  write(6,*) 'nspin2 from global.fdf ', nspin2

  ifock = fdf_get('ifock', 1)
  write(6,*) 'ifock from global.fdf ', ifock


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
  



  if (fdf_block('molecule', bfdf)) then
    !   External file reading
        write(6,*) 'beginning of external file coordinates block  '
        ia = 1
!        write(*,*) "linecount", fdf_block_linecount("molecule")
    
        do while((fdf_bline(bfdf, pline)))

          if (pline%ntokens == 1) then      
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
    endif
  write(6,*) 'Coordinates from Molecule block: External file'
  do ia= 1, natoms
    write(6,'(A4,3F10.6)') symbol(ia), (cent(i,ia),i=1,3)
  enddo
  
 

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'


!  Molecule coordinate block begins here  

  if (fdf_block('Coordinates', bfdf)) then
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

  if (fdf_block('Coordinates', bfdf)) then
    write(6,*) 'Coordinates:'
    do ia = 1, na
      write(6,'(A, 4x, 3F10.6)') symbol(ia), (xa(i,ia),i=1,3) 
    enddo
  endif


  write(6,*) '------------------------------------------------------'


  if (fdf_block('inline_xyz', bfdf)) then
!   Forward reading 
    write(6,*) 'Reading an inline_xyz block  '
    ia = 1

    do while((fdf_bline(bfdf, pline)))

      if (pline%ntokens == 1) then
        number_of_atoms = fdf_bintegers(pline, 1)
        write(*,*) "Number of atoms", number_of_atoms
      endif
      na = number_of_atoms

      if (pline%ntokens == 4) then
        symbol(ia) = fdf_bnames(pline, 1)
        do i= 1, 3
          xa(i,ia) = fdf_bvalues(pline, i)
        enddo
        ia = ia + 1
      endif
    enddo

    write(6,*) 'Inline XYZ Coordinates block:'
    do ia= 1, na
      write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    enddo
  endif

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'



  ! if (fdf_block('inline_xyz2', bfdf)) then
  !   !   Forward reading 
  !       write(6,*) 'Reading an inline_xyz2 block  '
  !       ia = 1
    
  !       do while(fdf_bline(bfdf, pline))
    
  !         if (pline%ntokens == 1) then      
  !           number_of_atoms = fdf_bintegers(pline, 1)
  !           write(*,*) "Number of atoms", number_of_atoms
  !         endif
  !         na = number_of_atoms
        
  !         if (pline%ntokens == 4) then
  !           symbol(ia) = fdf_bnames(pline, 1)
  !           do i= 1, 3
  !             xa(i,ia) = fdf_bvalues(pline, i)
  !           enddo
  !           ia = ia + 1
  !         endif
  !       enddo
    
  !       write(6,*) 'Inline XYZ2 Coordinates block:'
  !       do ia= 1, na
  !         write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
  !       enddo
  !     endif
    
  !     write(6,'(A)')  

  !     write(6,*) '------------------------------------------------------'
    
        
    !   if (fdf_block('molecule2', bfdf)) then
    !     !   External file reading
    !         write(6,*) 'beginning of external file coordinates block  '
    !         ia = 1
    ! !        write(*,*) "linecount", fdf_block_linecount("molecule")
        
    !         do while((fdf_bline(bfdf, pline)))
    
    !           if (pline%ntokens == 1) then      
    !             number_of_atoms = fdf_bintegers(pline, 1)
    !             write(*,*) "number of atoms", number_of_atoms
    !           endif
    !           na = number_of_atoms
            
    !           if (pline%ntokens == 4) then
    !             symbol(ia) = fdf_bnames(pline, 1)
    !             do i= 1, 3
    !               xa(i,ia) = fdf_bvalues(pline, i)
    !             enddo
    !             ia = ia + 1
    !           endif
    !         enddo
    !     endif
    !   write(6,*) 'Coordinates from Molecule2 block: External file'
    !   do ia= 1, na
    !     write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
    !   enddo
      
     
    
    !   write(6,'(A)')  
    
    !   write(6,*) '------------------------------------------------------'
    
    





  ! if ( fdf_block('ListBlock',bfdf) ) then
  !    i = 0
  !    do while ( fdf_bline(bfdf,pline) )
  !       i = i + 1
  !       na = fdf_bnlists(pline)
  !       write(*,'(2(a,i0),a)') 'Listblock line: ',i,' has ',na,' lists'
  !       do ia = 1 , na
  !          j = -1

  !          call fdf_bilists(pline,ia,j,isa)
  !          write(*,'(tr5,2(a,i0),a)') 'list ',ia,' has ',j,' entries' 
  !          call fdf_bilists(pline,ia,j,isa)
  !          write(*,'(tr5,a,1000(tr1,i0))') 'list: ',isa(1:j)

  !       end do
  !    end do
  ! end if


  ! if ( fdf_islreal('list_floats') .and. fdf_islist('list_floats') &
  !     .and. (.not. fdf_islinteger('list_floats')) ) then
  !   na = -1
  !   call fdf_list('list_floats',na,listr)
  !   write(*,'(tr1,a,i0,a)') 'list_floats has ',na,' entries'
  !   if ( na < 2 ) stop 1
  !   call fdf_list('list_floats',na,listr)
  !   write(*,'(tr5,a,1000(tr1,f12.8))') 'list_floats: ',listr(1:na)
  ! else
  !   write(*,*)'list_floats was not recognized'
  !   stop 1
  ! end if

  ! write(6,'(A)')  
    
  ! write(6,*) '------------------------------------------------------'

  write(6,'(A)')  " Determinants Block"

  write(6,*) '------------------------------------------------------'


  if (fdf_block('determinants', bfdf)) then
    !   External file reading
        write(6,*) 'Beginning of external file determinant block  '
        ia = 1

!        call io_status()
!        call fdf_printfdf()
        print*, "printing label ", bfdf%label , trim(bfdf%mark%pline%line)


        print*, "pline obtained",  (fdf_bline(bfdf, pline))

        open (unit=11,file='TZ_1M_500.det', iostat=iostat, action='read' )
        read(11,*) temp1, temp2, nelectrons, temp3, nalpha
        write(*,'(a,1x,i3,1x,i3)') "write after read1", nelectrons, nalpha        
        read(11,*)  temp1, ndeterminants, nexcitation   
        allocate(det_coeff(ndeterminants))             
        write(*,'(a,1x,i3, 1x, i3)') "write after read2", ndeterminants, nexcitation
        read(11,*) (det_coeff(i), i=1,ndeterminants)
        write(fmt,*)  '(', ndeterminants, '(f11.8,1x))'
        write(*,fmt) (det_coeff(i), i=1,ndeterminants)
!        write(*,'(<ndeterminants>(f11.8, 1x))') (det_coeff(i), i=1,ndeterminants)    ! for Intel Fortran        
        close(11)

        
          if(fdf_bsearch(pline, "&electrons")) then
            nelectrons  =  integers(pline, 1)
            nalpha      =  integers(pline, 2)
            nbeta       = nelectrons - nalpha

            write(*,*) "total number of       electrons ", nelectrons
            write(*,*) "      number of alpha electrons ", nalpha        
            write(*,*) "      number of beta  electrons ", nalpha
            if (.not. allocated(det_alpha)) allocate(det_alpha(nalpha))
            if (.not. allocated(det_beta)) allocate(det_beta(nbeta))
          endif

          if(fdf_bsearch(pline, "determinants")) then
            ndeterminants   =  fdf_bintegers(pline, 1)          
            nexcitation     =  fdf_bintegers(pline, 2)
            write(*,*) "total number of determinants ", ndeterminants
            write(*,*) "      number of excitations  ", nexcitation 
            if (.not. allocated(det_coeff)) allocate(det_coeff(ndeterminants)) 
          endif


          na = nintegers(pline)
          write(*,'(2(a,i0),a)') 'number of integers: ', na, ' integers'
          write(*,'(tr5,a,<nalpha>(tr1,i0))') 'list: ', det_alpha(1:nalpha)


!          endif
!        enddo
  endif

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'







  call fdf_shutdown()

!----------------------------------------------------------------------------END
END PROGRAM iochamp
