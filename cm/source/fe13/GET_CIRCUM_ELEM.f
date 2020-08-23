      SUBROUTINE GET_CIRCUM_ELEM(CIRCUM_XI,FIRST_ELEMENT,NXI,
     '     CIRCUM_ELEM,ERROR, *)

C#### Subroutine: GET_CIRCUM_ELEM
C###  Description:
C###    Store the list of elements located on the circumferential
C###    line around a surface element mesh
C###    values are stored in CIRCUM_ELEM
C###    CIRCUM_ELEM(0) contains the number of elements around the 
C###    circumference
C*** Created by Peter Bier, Mar 2003.

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'

!     Parameter list
      INTEGER CIRCUM_XI,FIRST_ELEMENT,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '     CIRCUM_ELEM(0:NEM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ne,ni
      LOGICAL FIN

      CALL ENTERS('GET_CIRCUM_ELEM',*9999)
C     build a list of circumference elements to process

      ne = FIRST_ELEMENT
      ni = 1                    ! CIRCUM_ELEM array index
      CIRCUM_ELEM(ni) = ne      ! store first element
      FIN = .FALSE.
      DO WHILE(.NOT.FIN)            
         ne = NXI(CIRCUM_XI,1,ne) ! fetch next elem in circum_xi dir
         
         CALL ASSERT(ne.NE.0,'>>did not find element while ' //
     '        'looping around the circumference' //
     '        ': check circumference xi direction is correct',
     '        ERROR,*9999)

         IF( ne .NE. FIRST_ELEMENT ) THEN
            ni = ni + 1
            CIRCUM_ELEM(ni) = ne
         ELSE                   ! back to first element
            FIN = .TRUE.
         ENDIF            
      ENDDO
      
      CIRCUM_ELEM(0) = ni       ! store total number of circumferential elems
      
C         WRITE(*,*) 'List of circumferential elements'
C         DO ni=1,CIRCUM_ELEM(0)
C            WRITE(*,*) CIRCUM_ELEM(ni)
C         ENDDO

      CALL EXITS('GET_CIRCUM_ELEM')
      RETURN
 9999 CALL ERRORS('GET_CIRCUM_ELEM',ERROR)
      CALL EXITS('GET_CIRCUM_ELEM')
      RETURN 1
      END


