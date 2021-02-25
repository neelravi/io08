module periodic_table

    public :: element, atom_t

    type :: atom_t
        character(len=20)        :: name
        character(len=3)        :: symbol        
        integer                 :: znuclear
        real                    :: atomic_mass
    end type atom_t

    contains 

    type(atom_t) function element(sym) result(atom)

    character(len=*) :: sym

    select case(sym)
        case("hydrogen", "Hydrogen", "H", "1")
                atom = atom_t(name="hydrogen", symbol="H", atomic_mass=1.00794, znuclear=1)
        case default
                error stop "Unknown element or symbol"
    end select
    return 
    end function      

end module