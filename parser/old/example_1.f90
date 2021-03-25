program test_m_config
     use m_config

     integer, parameter    :: dp = kind(0.0d0)
     type(CFG_t)           :: my_cfg

     ! Some dummy variables
     real(dp), allocatable :: trial_energy(:)
     integer               :: n_reals
     character(len=20)     :: fmt_string

     character(len=20)     :: sections  
     logical               :: optimize_wavefunction, optimize_ci
     logical               :: optimize_jastrow, optimize_orbitals

     ! general block
     character(len=100)    ::  title, filename, molecule
     character(len=50)     ::  output_directory 
     character(len=50)     ::  pool
     character(len=50)     ::  basis
     character(len=50)     ::  pseudo

     ! mixed block

     character(len=20)     :: unit 
     integer               :: maximum_iterations 
     logical               :: restart_vmc 



     ! title and external files
     call CFG_add(my_cfg, "title", "this/is/a/filename", &
          "A string containing a filename")

     call CFG_add(my_cfg, "filename", "this/is/a/filename", &
          "A string containing a filename")

     call CFG_add(my_cfg, "molecule", "h2o.xyz", &
          "Molecule's coordinates in xyz file format")



     ! General block
     call CFG_add(my_cfg, "general%output_directory", "./", &
          "output directory")

     call CFG_add(my_cfg, "general%pool", "./pool", &
          "a pool directory containing required files")

     call CFG_add(my_cfg, "general%basis", "./pool/basis", &
          "a basis file with its location")

     call CFG_add(my_cfg, "general%pseudo", "./pool/pseudo", &
          "a pseudopotential file with its location")


     ! a block containing mixed data
     call CFG_add(my_cfg, "mixed%unit", "eV", &
          "Energy unit")

     call CFG_add(my_cfg, "mixed%maximum_iterations", 250, &
          "Maximum iterations")       

     call CFG_add(my_cfg, "mixed%trial_energy", (/13.37_dp, 13.40_dp, 13.80_dp , 14.00_dp /), &
          "Trial energies", dynamic_size=.true.)

     call CFG_add(my_cfg, "mixed%restart_vmc", .true., &
          "Restart VMC ? ")


     ! optimization block logical
     call CFG_add(my_cfg, "optimization_flags%optimize_wavefunction", .false., &
          "optimize wavefunctions")

     call CFG_add(my_cfg, "optimization_flags%optimize_ci", .false., &
          "optimize ci")

     call CFG_add(my_cfg, "optimization_flags%optimize_orbitals", .false., &
          "optimize orbitals")

     call CFG_add(my_cfg, "optimization_flags%optimize_jastrow", .false., &
          "optimize jastrow")


     ! Sort the configuration (this can speed up looking for variables, but only if
     ! you have a sufficiently large number of them)
     call CFG_sort(my_cfg)


     print *, "Reading in example_1_input.cfg"
     call CFG_read_file(my_cfg, "example_1_input.cfg") ! Update values with file

     print *, "----------------------------------------"

     print *, "----------------------------------------"
     print *, "The code below demonstrates how to get values: "
     print *, "----------------------------------------"
     print *, ""
     !  Ravindra added stuff

     ! title and external files
     call CFG_get(my_cfg, "title", title)
     call CFG_get(my_cfg, "filename", filename)
     call CFG_get(my_cfg, "molecule", molecule)  


     call CFG_get(my_cfg, "general%output_directory", output_directory)
     call CFG_get(my_cfg, "general%pool", pool)
     call CFG_get(my_cfg, "general%basis", basis)
     call CFG_get(my_cfg, "general%pseudo", pseudo)


     call CFG_get(my_cfg, "mixed%unit", unit)
     call CFG_get(my_cfg, "mixed%maximum_iterations", maximum_iterations)       
     call CFG_get(my_cfg, "mixed%restart_vmc", restart_vmc)

     call CFG_get_size(my_cfg, "mixed%trial_energy", n_reals)
     ! Generate format string for trial energy values
     write(fmt_string, "(A,I0,A)") "(A25,", n_reals, "F10.5)"

     allocate(trial_energy(n_reals))
     call CFG_get(my_cfg, "mixed%trial_energy", trial_energy)
!     write(*, fmt_string) "Trial Energies ", trial_energy
     deallocate(trial_energy)


     call CFG_get(my_cfg, "optimization_flags%optimize_wavefunction", optimize_wavefunction)
     call CFG_get(my_cfg, "optimization_flags%optimize_ci", optimize_ci)
     call CFG_get(my_cfg, "optimization_flags%optimize_orbitals", optimize_orbitals)
     call CFG_get(my_cfg, "optimization_flags%optimize_jastrow", optimize_jastrow)


     ! final printing part
     call CFG_write(my_cfg, "stdout") ! Write to stdout
     call CFG_write(my_cfg, "example_1_output.cfg") ! Write to file
     call CFG_write_markdown(my_cfg, "example_1_output.md") ! Write markdown file

end program test_m_config
