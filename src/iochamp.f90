!
!     Shows FDF capabilities..
!
PROGRAM iochamp
  USE fdf
  USE prec
  USE parse
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug, check, val, logic(10)
  logical                    :: beginning, ending
  character(len=72)          :: fname, axis, status, filename, title
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
  integer                    :: nextorb, nblk_max, nopt_iter, max_iteration, max_iter, linecount, argument(5)
  real(dp)                   :: energy_tol, float_value
  real(dp)                   :: sr_tau, sr_eps, sr_adiag 
  character(len=20)          :: real_format = '(A, T20, F14.5)'
  character(len=20)          :: int_format = '(A, T20, I8)'

!------------------------------------------------------------------------- BEGIN

! Initialize
  call fdf_init('test-champ.inp', 'test-champ.out')

! Handle/Use fdf structure
  if (fdf_defined('new-style')) write(6,*) 'New-style stuff'

! strings/characters
  fname = fdf_string('title', 'Default title')
  write(6,'(A)') 'title of the calculation :: ', fname

!Integer numbers (keyword, default_value). The variable is assigned default_value when keyword is not present
  nextorb = fdf_integer('nextorb', 0)
  write(6,fmt=int_format) 'Next Orb =', nextorb

  nblk_max = fdf_integer('nblk_max', 0)
  write(6,fmt=int_format) 'nblk max =', nblk_max

  nopt_iter = fdf_integer('nopt_iter', 0)
  write(6,fmt=int_format) 'nopt_iter =', nopt_iter


! floats (keyword, default_value) variable is assigned default_value when keyword is not present
  sr_tau = fdf_get('sr_tau', 0.025d0)
  write(6,fmt=real_format) 'sr_tau:', sr_tau

  sr_eps = fdf_get('sr_eps', 0.001d0)
  write(6,fmt=real_format) 'sr_eps:', sr_eps

  sr_adiag = fdf_get('sr_adiag', 0.01d0)
  write(6,fmt=real_format) 'sr_adiag:', sr_adiag

  energy_tol = fdf_get('energy_tol', 0.00001d0)
  write(6,fmt=real_format) 'energy_tol:', energy_tol

  nopt_iter = fdf_integer('a', 0)
  write(6,fmt=int_format) 'a =', nopt_iter

  nopt_iter = fdf_integer('b', 0)
  write(6,fmt=int_format) 'b =', nopt_iter



  ! logical :: true, .true., yes, T, and TRUE are equivalent
  debug = fdf_boolean('Debug', .TRUE.)
  write(6,'(A, L2)') 'Debug:', debug


! mixed types in one line (for example, reading a number with units)

  max_iter = fdf_integer('max_iteration', 100)
  write(6,*) 'Examples: maximum_iter =', max_iter


  float_value = fdf_get('float_value', 0.00001d0)
  write(6,*) 'float_value :: ', float_value


  cutoff = fdf_physical('Energy_Cutoff', 8.d0, 'Ry')
  write(6,fmt=real_format) 'Energy CutOff in Rydberg :: ', cutoff

  phonon_energy = fdf_physical('phonon-energy', 0.01d0, 'eV')
  write(6,fmt=real_format) 'Phonon Energy in eV :: ', phonon_energy

  write(6,'(A)')  

  write(6,*) '------------------------------------------------------'


  ! block containing logical key-value pairs
  doit = fdf_boolean("optimize_wavefunction", .True.)
  write(6,*) ' optimize_wavefunction = ', doit

  doit = fdf_boolean('optimize_ci', .True.)
  write(6,*) ' optimize_ci = ', doit

  doit = fdf_boolean('optimize_jastrow', .True.)
  write(6,*) ' optimize_jastrow = ', doit

  doit = fdf_boolean('optimize_orbitals', .True.)
  write(6,*) ' optimize_orbitals = ', doit



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
            number_of_atoms = fdf_bintegers(pline, 1)
            write(*,*) "number of atoms", number_of_atoms
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
    endif
  write(6,*) 'Coordinates from Molecule block: External file'
  do ia= 1, na
    write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
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

  write(6,*) 'Coordinates:'
  do ia = 1, na
    write(6,'(A, 4x, 3F10.6)') symbol(ia), (xa(i,ia),i=1,3) 
  enddo

!  Molecule coordinate block ends here
  write(6,*) '------------------------------------------------------'

  write(6,'(A)')  

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



  if (fdf_block('inline_xyz2', bfdf)) then
    !   Forward reading 
        write(6,*) 'Reading an inline_xyz2 block  '
        ia = 1
    
        do while(fdf_bline(bfdf, pline))
    
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
    
        write(6,*) 'Inline XYZ2 Coordinates block:'
        do ia= 1, na
          write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
        enddo
      endif
    
      write(6,'(A)')  

      write(6,*) '------------------------------------------------------'
    
        



  if ( fdf_block('ListBlock',bfdf) ) then
     i = 0
     do while ( fdf_bline(bfdf,pline) )
        i = i + 1
        na = fdf_bnlists(pline)
        write(*,'(2(a,i0),a)') 'Listblock line: ',i,' has ',na,' lists'
        do ia = 1 , na
           j = -1

           call fdf_bilists(pline,ia,j,isa)
           write(*,'(tr5,2(a,i0),a)') 'list ',ia,' has ',j,' entries' 
           call fdf_bilists(pline,ia,j,isa)
           write(*,'(tr5,a,1000(tr1,i0))') 'list: ',isa(1:j)

        end do
     end do
  end if


  if ( fdf_islreal('list_floats') .and. fdf_islist('list_floats') &
      .and. (.not. fdf_islinteger('list_floats')) ) then
    na = -1
    call fdf_list('list_floats',na,listr)
    write(*,'(tr1,a,i0,a)') 'list_floats has ',na,' entries'
    if ( na < 2 ) stop 1
    call fdf_list('list_floats',na,listr)
    write(*,'(tr5,a,1000(tr1,f12.8))') 'list_floats: ',listr(1:na)
  else
    write(*,*)'list_floats was not recognized'
    stop 1
  end if


  call fdf_shutdown()

!----------------------------------------------------------------------------END
END PROGRAM iochamp
