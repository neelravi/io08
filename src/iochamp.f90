!
!     Shows FDF capabilities..
!
PROGRAM iochamp
  USE fdf
  USE prec
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug, check, val, logic(10)
  logical                    :: beginning, ending
  character(len=72)          :: fname, axis, status, filename, title
  character(len=72)          :: molecule_name, key, comment
  character(2)               :: symbol(maxa)
  character(len=20)          :: chunks(10), subblock(10)
  character(len=30)          :: keyword(5), argument(5)
  integer(sp)                :: i, j, ia, na, external_entry, number_of_atoms
  integer(sp)                :: isa(maxa)
  real(sp)                   :: wmix
  real(dp)                   :: cutoff, phonon_energy, factor
  real(dp)                   :: xa(3, maxa)
  real(dp)                   :: listr(maxa)
  type(block_fdf)            :: bfdf, bfdf2
  type(parsed_line), pointer :: pline, pline2
  !type(fdf_file)             :: fdffile 
  integer                    :: nextorb, nblk_max, nopt_iter, max_iteration, max_iter, linecount
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

! Integer numbers (keyword, default_value). The variable is assigned default_value when keyword is not present
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


! block containing other key-value pairs :: A temporary workaround  
  if (fdf_block('optimization_flags', bfdf)) then
    write(*,*) "inside optimization_flags block"
    linecount = fdf_block_linecount("optimization_flags")    
    i = 1
    do while(fdf_bline(bfdf, pline))   
!      write(*,*) "some debug info pline" , pline%ntokens, pline%line
      keyword(i) = fdf_bnames(pline, 1)
      argument(i) = fdf_bnames(pline, 2)
      i = i + 1
    enddo
  endif

  do i = 1, linecount
    write(6,'(*(A,4X,A))') keyword(i), argument(i)
  enddo


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
  

!  This block currently fails as it is not possible to parse within the scope of a block
  doit = fdf_get("optimize_wavefunction", .false.)
  write(6,*) 'outside  optimize_wavefunction', doit

  doit = fdf_get('opt.optimize_ci', .false.)
  write(6,*) 'outside  optimize_ci', doit

  doit = fdf_get('opt.optimize_jastrow', .false.)
  write(6,*) 'outside  optimize_jastrow:', doit

  doit = fdf_get('opt.optimize_orbitals', .false.)
  write(6,*) 'outside  optimize_orbitals:', doit




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
    
          if (pline%ntokens == 1) then      
!            molecule_name =  fdf_bnames(pline, 1) 
            molecule_name =  fdf_get('molecule', 'Unknown molecule') 
          endif
    
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

    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))

      if (pline%ntokens == 1) then      
        number_of_atoms = fdf_bintegers(pline, 1)
        write(*,*) "Number of atoms", number_of_atoms
      endif
      na = number_of_atoms

      if (pline%ntokens == 1) then      
        molecule_name =  fdf_string('', 'Unknown molecule') 
      endif

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


  if (fdf_block('inline_xyz2', bfdf)) then
    !   Forward reading 
        write(6,*) 'Reading an inline_xyz block  '
        ia = 1
    
        do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
    
          if (pline%ntokens == 1) then      
            number_of_atoms = fdf_bintegers(pline, 1)
            write(*,*) "Number of atoms", number_of_atoms
          endif
          na = number_of_atoms
    
          if (pline%ntokens == 1) then      
            molecule_name =  fdf_string('', 'Unknown molecule') 
          endif
    
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
    
    



!   if ( fdf_block('ListBlock',bfdf) ) then
!      i = 0
!      do while ( fdf_bline(bfdf,pline) )
!         i = i + 1
!         na = fdf_bnlists(pline)
!         write(*,'(2(a,i0),a)') 'Listblock line: ',i,' has ',na,' lists'
!         do ia = 1 , na
!            j = -1
!            call fdf_bilists(pline,ia,j,isa)
!            write(*,'(tr5,2(a,i0),a)') 'list ',ia,' has ',j,' entries'
!            call fdf_bilists(pline,ia,j,isa)
!            write(*,'(tr5,a,1000(tr1,i0))') 'list: ',isa(1:j)
!         end do
!      end do
!   end if

  ! Check lists
  if ( fdf_islinteger('MyList') .and. fdf_islist('MyList') &
      .and. (.not. fdf_islreal('MyList')) ) then
     na = -1
     call fdf_list('MyList',na,isa)
     if ( na < 2 ) stop 1
     write(*,'(tr1,a,i0,a)') 'MyList has ',na,' entries'
     call fdf_list('MyList',na,isa)
     write(*,'(tr5,a,1000(tr1,i0))') 'MyList: ',isa(1:na)
   else
     write(*,*)'MyList was not recognized'
     stop 1
   end if

  if ( fdf_islinteger('MyListOne') .and. fdf_islist('MyListOne') &
      .and. (.not. fdf_islreal('MyListOne')) ) then
     na = -1
     call fdf_list('MyListOne',na,isa)
     if ( na /= 1 ) stop 1
     write(*,'(tr1,a,i0,a)') 'MyListOne has ',na,' entries'
     call fdf_list('MyListOne',na,isa)
     write(*,'(tr5,a,1000(tr1,i0))') 'MyListOne: ',isa(1:na)
  else
     write(*,*)'MyListOne was not recognized'
     stop 1
  end if

  if ( fdf_islreal('MyListR') .and. fdf_islist('MyListR') &
      .and. (.not. fdf_islinteger('MyListR')) ) then
    na = -1
    call fdf_list('MyListR',na,listr)
    write(*,'(tr1,a,i0,a)') 'MyListR has ',na,' entries'
    if ( na < 2 ) stop 1
    call fdf_list('MyListR',na,listr)
    write(*,'(tr5,a,1000(tr1,f4.1))') 'MyListR: ',listr(1:na)
  else
    write(*,*)'MyListR was not recognized'
    stop 1
  end if

  if ( fdf_islreal('MyListROne') .and. fdf_islist('MyListROne') &
      .and. (.not. fdf_islinteger('MyListROne')) ) then
    na = -1
    call fdf_list('MyListROne',na,listr)
    if ( na /= 1 ) stop 1
    write(*,'(tr1,a,i0,a)') 'MyListROne has ',na,' entries'
    call fdf_list('MyListROne',na,listr)
    write(*,'(tr5,a,1000(tr1,f4.1))') 'MyListROne: ',listr(1:na)
  else
    write(*,*)'MyListROne was not recognized'
    stop 1
  end if

  if ( fdf_islist('externalentry') ) then
     write(*,*) 'externalentry is a list'
  else
     write(*,*) 'externalentry is not a list'
  end if

  external_entry = fdf_integer('externalentry', 60)
  write(6,*) 'ExternalEntry:', external_entry

  axis   = fdf_string('AxisXY', 'Cartesian')
  status = fdf_string('StatusXY', 'Enabled')
  write(6,*) 'Axis: ', TRIM(axis), ' | ', TRIM(status)

! Shutdown and deallocates fdf structure
  call fdf_shutdown()

!----------------------------------------------------------------------------END
END PROGRAM iochamp
