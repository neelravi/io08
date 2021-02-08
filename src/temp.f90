!
!
! ... subroutine field_count:   accepts two string (one of them is optional) 
!                               and one integer and count the number of fields
!                               in the string separated by a blank or a tab 
!                               character. if the optional string is specified
!                               (it has anyway len=1) it is assumed as the 
!                               separator character.
!                               ignores any character following the exclamation
!                               mark (fortran comment)
!
! ... subroutine con_cam:       counts the number of fields in a string 
!                               separated by the optional character
!
! ... subroutine field_compare: accepts two strings and one integer. counts the
!                               fields contained in the first string and 
!                               compares it with the integer. 
!                               if they are less than the integer calls the 
!                               routine error and show by the second string the
!                               name of the field where read-error occurred.
!
!
!----------------------------------------------------------------------------
module parser
  !----------------------------------------------------------------------------
  !
  integer, parameter :: stdin  = 5    ! unit connected to standard input
  integer :: stdout = 6    ! unit connected to standard output

  integer, parameter :: dp = selected_real_kind(16,200)
  !
  private
  !
  public :: parse_unit, field_count, read_line, get_field
  public :: version_parse, version_compare
  !
  integer :: parse_unit = 5 ! normally 5, but can be set otherwise
  !
  contains
  !
  !
  !--------------------------------------------------------------------------
  subroutine field_count( num, line, car )
    !--------------------------------------------------------------------------
    !
    implicit none
    !
    integer,                    intent(out) :: num
    character(len=*),           intent(in)  :: line
    character(len=1), optional, intent(in)  :: car
    character(len=1)                        :: sep1, sep2    
    integer                                 :: j
    !
    !
    num = 0
    !
    if ( .not. present(car) ) then
       !
       sep1 = char(32)  ! ... blank character
       sep2 = char(9)   ! ... tab character
       !
       do j = 2, max( len( line ), 256 )
          !
          if ( line(j:j) == '!' .or. line(j:j) == char(0) ) then
             !
             if ( line(j-1:j-1) /= sep1 .and. line(j-1:j-1) /= sep2 ) then
                !
                num = num + 1
                !
             end if   
             !
             exit
             !
          end if
          !
          if ( ( line(j:j) == sep1 .or. line(j:j) == sep2 ) .and. &
               ( line(j-1:j-1) /= sep1 .and. line(j-1:j-1) /= sep2 ) ) then
             !
             num = num + 1
             !
          end if
          !
       end do
       !
    else
       !
       sep1 = car
       !
       do j = 2, max( len( line ), 256 )
          ! 
          if ( line(j:j) == '!' .or. &
               line(j:j) == char(0) .or. line(j:j) == char(32) ) then
             !
             if ( line(j-1:j-1) /= sep1 ) num = num + 1
             !
             exit
             !
          end if
          !
          if ( line(j:j) == sep1 .and. line(j-1:j-1) /= sep1 ) num = num + 1
          !
       end do
       !
    end if
    !
    return
    !
  end subroutine field_count
  !
  !
  !--------------------------------------------------------------------------
  subroutine read_line( line, nfield, field, end_of_file, error )
    !--------------------------------------------------------------------------
    !
    use mp,        only : mp_bcast
    use mp_images, only : intra_image_comm
    use io_global, only : ionode, ionode_id
    !
    implicit none
    !
    character(len=*),           intent(out) :: line
    character(len=*), optional, intent(in)  :: field
    integer,          optional, intent(in)  :: nfield
    logical,          optional, intent(out) :: end_of_file, error
    logical                                 :: tend, terr
    !
    !
    if( len( line ) < 256 ) then
       call errore(' read_line ', ' input line too short ', max(len(line),1) )
    end if
    !
    tend = .false.
    terr = .false.
    if ( ionode ) then
30     read (parse_unit, fmt='(a256)', err=15, end=10) line
       if( line == ' ' .or. line(1:1) == '#' ) go to 30
       go to 20
10     tend = .true.
       go to 20
15     terr = .true.
20     continue
    end if
    !
    call mp_bcast( tend, ionode_id, intra_image_comm )
    call mp_bcast( terr, ionode_id, intra_image_comm )
    call mp_bcast( line, ionode_id, intra_image_comm )
    !
    if( present(end_of_file) ) then
       end_of_file = tend
    else if( tend ) then
       call infomsg(' read_line ', ' end of file ' )
    end if
    if( present(error) ) then
       error = terr
    else if( terr ) then
       call infomsg(' read_line ', ' read error ' )
    end if
    if( present(field) .and. .not.(tend.or.terr) ) &
     &call field_compare( line, nfield, field )
    !
  end subroutine read_line
  !
  !
  !--------------------------------------------------------------------------
  subroutine field_compare( str, nf, var )
    !--------------------------------------------------------------------------
    !
    implicit none
    !
    character(len=*), intent(in) :: var
    integer,          intent(in) :: nf
    character(len=*), intent(in) :: str
    integer                      :: nc
    !
    call field_count( nc, str )
    !
    if( nc < nf ) &
      call errore( ' field_compare ', &
                 & ' wrong number of fields: ' // trim( var ), 1 )
    !
    return
    !
  end subroutine field_compare
  !
  !
  !--------------------------------------------------------------------------
  subroutine con_cam(num, line, car)
    !--------------------------------------------------------------------------
    character(len=*) :: line
    character(len=1) :: sep
    character(len=1), optional :: car
    integer :: num, j

    num = 0
    if (len(line) .gt. 256 ) then
       write( stdout,*) 'riga ', line
       write( stdout,*) 'lunga ', len(line)
       num = -1
       return
    end if

    write( stdout,*) '1riga ', line
    write( stdout,*) '1lunga ', len(line)
    if ( .not. present(car) ) then
       sep=char(32)             !char(32) is the blank character
    else
       sep=car
    end if

    do j=2, max(len(line),256)
       if ( line(j:j) == '!' .or. line(j:j) == char(0)) then
          return
       end if
       if ( (line(j:j) .eq. sep) .and. &
            (line(j-1:j-1) .ne. sep) )  then
          num = num + 1
       end if
    end do
    return
  end subroutine con_cam
  !
  !--------------------------------------------------------------------------
  subroutine get_field(n, field, str, sep)
    !--------------------------------------------------------------------------
    ! extract whitespace-separated nth block from string
    implicit none
    integer,intent(in) :: n
    character(len=*),intent(out) :: field
    character(len=*),intent(in)  :: str
    character(len=1),optional,intent(in) :: sep
    integer :: i,j,z ! block start and end
    integer :: k     ! block counter
    character(len=1) :: sep1, sep2
    !print*, "------------- parser start -------------"
    !print '(3a)', "string: -->", str,"<--"
    if(present(sep)) then
      sep1 = sep
      sep2 = sep ! redundant, but easy
    else
      sep1 = char(32)  ! ... blank character
      sep2 = char(9)   ! ... tab char
    endif
    !
    k = 1 ! counter for the required block
    !
    do i = 1,len(str)
    ! look for the beginning of the required block
      z = max(i-1,1)
      !print '(2a1,3i4,2l)', str(i:i), str(z:z), i,z,k,n,&
      !       (str(i:i) == sep1 .or. str(i:i) == sep2), (str(z:z) /= sep1 .and. str(z:z) /= sep2)
      if( k == n) exit
      if( (str(i:i) == sep1 .or. str(i:i) == sep2) &
           .and. &
          (str(z:z) /= sep1 .and. str(z:z) /= sep2) &
        ) &
        k = k+1
    enddo
    !
    !print*, "i found: ",i
    do j = i,len(str)
    ! look for the beginning of the next block
      z = max(j-1,1)
      if( (str(j:j) == sep1 .or. str(j:j) == sep2) &
           .and. &
          (str(z:z) /= sep1 .and. str(z:z) /= sep2) &
        ) &
        k = k+1
      if( k >n) exit
    enddo
    !print*, "j found: ",j
    !
    if (j <= len(str)) then
      ! if we are here, the reqired block was followed by a separator
      ! and another field, we have to trash one char (a separator)
      field = trim(adjustl(str(i:j-1)))
      !print*, "taking: ",i,j-2
    else
      ! if we are here, it was the last block in str, we have to take
      ! all the remaining chars
      field = trim(adjustl(str(i:len(str))))
      !print*, "taking from ",i
    endif
    !print*, "------------- parser end -------------"

  end subroutine get_field

end module parser
