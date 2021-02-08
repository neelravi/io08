!
! copyright (c) 2011 quantum espresso group
! this file is distributed under the terms of the
! gnu general public license. see the file `license'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! read input data in pw and cp from text file (xml file to be implemented)
! note aug 2018 (pg): reading of old xml input file using iotk deleted
! 
!----------------------------------------------------------------------------
module read_input
   !---------------------------------------------------------------------------
   !
   use kinds,     only: dp
   !
   implicit none
   save
   !
   private
   public :: read_input_file, has_been_read
   !
   logical :: has_been_read = .false.
   !
   contains
   !
   !-------------------------------------------------------------------------
   subroutine read_input_file ( prog, input_file_ )
     !-------------------------------------------------------------------------
     !
     use input_parameters,      only : reset_input_checks
     use read_namelists_module, only : read_namelists
     use read_cards_module,     only : read_cards
     use io_global,             only : ionode, ionode_id, qestdin
     use mp,                    only : mp_bcast
     use mp_images,             only : intra_image_comm
     use open_close_input_file, only : open_input_file, close_input_file
     !
     implicit none
     !
     character(len=*), intent (in) :: prog
     character(len=*), intent (in) :: input_file_
     !
     logical :: xmlinput
     integer :: ierr
     !
     if ( ionode ) ierr = open_input_file( input_file_, xmlinput )
     call mp_bcast( ierr, ionode_id, intra_image_comm )
     if ( ierr > 0 ) call errore('read_input', 'opening input file',ierr)
     call mp_bcast( xmlinput, ionode_id, intra_image_comm )
     !
     call reset_input_checks () 
     !
     if ( xmlinput ) then
        !
        call errore('read_input', 'xml input disabled',1)
        !
     else
        !
        ! ... read namelists 
        !
        call read_namelists( prog, qestdin )
        !
        ! ... read cards (requires in input only first two letters of prog)
        !
        call read_cards ( prog(1:2), qestdin )
        !
     end if
     if ( ionode) ierr = close_input_file( )
     !
     has_been_read = .true.
     !
     return
     !
   end subroutine read_input_file
  !
end module read_input
