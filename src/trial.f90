program trial_reading
    use iso_fortran_env
    use periodic_table, only: atom_t, element


    integer ::  nalpha, nbeta, nelectrons, ndeterminants, nexcitation
    character(40)   :: temp1, temp2, temp3, fmt
    real(selected_real_kind(6,15)), allocatable :: det_coeff(:)

        type(atom_t) :: atom1

        atom1 = element("Hydrogen")

        print*, "atom info name ", atom1%name
        print*, "atom info symbol ", atom1%symbol
        print*, "atom info mass ", atom1%atomic_mass
        print*, "atom info charge ", atom1%znuclear                        

        print*, atom1

        open (unit=11,file='TZ_1M_500.det', iostat=iostat, action='read' )
        read(11,*) temp1, temp2, nelectrons, temp3, nalpha
        write(*,'(a,1x,i3,1x,i3)') "write after read1", nelectrons, nalpha        
        read(11,*)  temp1, ndeterminants, nexcitation        
        allocate(det_coeff(ndeterminants))
        write(*,'(a,1x,i3, 1x, i3)') "write after read2", ndeterminants, nexcitation
        read(11,*) (det_coeff(i), i=1,500)
        write(fmt,*)  '(', ndeterminants, '(f11.8,1x))'
        write(*,fmt) (det_coeff(i), i=1,500)
!        write(*,'(<ndeterminants>(f10.8, 1x))') (det_coeff(i), i=1,ndeterminants)
        close(11)

end program        