! -*-f90-*-
MODULE fileio

	USE constants
	IMPLICIT NONE
  ! IO buffer
  CHARACTER(LEN=1024)::buff

CONTAINS

  !-----------------------------------------------------------------

  SUBROUTINE GetLine( unit,buff,verbose,ierr )

    ! Read from a unit until a non-blank line is reached.
    ! The read line is returned in buff.

    IMPLICIT NONE
    ! Subroutine args
    INTEGER(I4)::unit,ierr
    LOGICAL::verbose
    CHARACTER(LEN=*)::buff

    READ(unit,'(a)',iostat=ierr) buff
    IF (ierr.ne.0) RETURN
    buff=ADJUSTL(buff)
    IF (verbose) WRITE(stdout,'(a)') TRIM(buff)
    DO WHILE (buff.eq.' ' .and. ierr.eq.0)
      READ(unit,'(a)',iostat=ierr) buff
      IF (ierr.ne.0) RETURN
      buff=ADJUSTL(buff)
      IF (verbose) WRITE(stdout,'(a)') TRIM(buff)
    ENDDO

  END SUBROUTINE GetLine

  !-----------------------------------------------------------------

  INTEGER FUNCTION StrTokCnt( buff )

    ! Find the number of whitespace separated tokens in the string
    ! "buff"

    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*)::buff
    ! Local variables
    INTEGER(I4)::i,numtok
    LOGICAL::outside

    outside=.true.
    numtok=0

    DO i=1,LEN(buff)
      IF (buff(i:i).eq.' ' .or. buff(i:i).eq.TAB) THEN
        outside=.true.
      ELSE
        IF (outside) THEN
          outside=.false.
          numtok=numtok+1
        ENDIF
      ENDIF
    ENDDO

    StrTokCnt=numtok
    RETURN
  END FUNCTION StrTokCnt

  !-----------------------------------------------------------------

  SUBROUTINE CleanFileName( FullFileName,FileName,Extension )

    ! Clean extensions on a filename

    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*)::FullFileName
    CHARACTER(LEN=*)::FileName
    CHARACTER(LEN=*)::Extension
    ! Local variables
    INTEGER(I4)::file_len,ex_len

    FileName=ADJUSTL(FileName)

    ex_len=len_trim(Extension)
    file_len=len_trim(FileName)

    IF (FileName(file_len-ex_len:file_len).eq.Extension) THEN
      FullFilename=TRIM(FileName)
      FileName=FileName(1:file_len-(ex_len+1))
    ELSE
      FullFileName=TRIM(FileName)//TRIM(Extension)
    ENDIF

    RETURN
  END SUBROUTINE CleanFileName

  !-----------------------------------------------------------------

  SUBROUTINE Usage()

    ! Produce usage information for the user

    IMPLICIT none

    WRITE(stderr,'(a)') 'Usage: create_thickness < args >'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') ' args: -file <filename>   default file name'
    WRITE(stderr,'(a)') '       -slice <filename>  slice file'
    WRITE(stderr,'(a)') '       -thick <filename>  thickness file'
    WRITE(stderr,'(a)') '       -node <filename>   input ipnode file'
    WRITE(stderr,'(a)') '       -elem <filename>   input ipelem file'
    WRITE(stderr,'(a)') '       -output <filename> output file name'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -prompt_slice      prompt for slice data'
    WRITE(stderr,'(a)') '       -read_slice        read slice data from file'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -prompt_thick      prompt for thickness data'
    WRITE(stderr,'(a)') '       -read_thick        read thickness data from file'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -new_versions      versions added to new points'
    WRITE(stderr,'(a)') '       -no_new_versions   no versions on new nodes'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -one_node          add one node per slice'
    WRITE(stderr,'(a)') '       -one_node...not    add thickness per slice'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -ratio_thick       specify thickness by ratio'
    WRITE(stderr,'(a)') '       -ratio_thick...not specify thickness by dimension'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -add_basis <#>     add basis number to the element file'
    WRITE(stderr,'(a)') '       -offset <#>        node number offset'
    WRITE(stderr,'(a)') '       -fudge <#>         scaling fudge factor'
    WRITE(stderr,'(a)') '       -ratio <#>         value for projecting inwards'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -redist <a,b,c,d,e:a,c,e> redistribute the nodes'
    WRITE(stderr,'(a)') '                       a,b,c,d,e along a spline defined'
    WRITE(stderr,'(a)') '                       by the nodes a,c,e. The spline '
    WRITE(stderr,'(a)') '                       definition syntax requires that'
    WRITE(stderr,'(a)') '                       the start and end points of the'
    WRITE(stderr,'(a)') '                       spline (a and e of the first set)'
    WRITE(stderr,'(a)') '                       be the first and last points of'
    WRITE(stderr,'(a)') '                       the points to move (a and e of the '
    WRITE(stderr,'(a)') '                       second set)'
    WRITE(stderr,*)
    WRITE(stderr,'(a)') '       -verb              verbose output'
    WRITE(stderr,'(a)') '       -help              display this usage info'
    STOP
  END SUBROUTINE Usage

  !-----------------------------------------------------------------

END MODULE fileio
