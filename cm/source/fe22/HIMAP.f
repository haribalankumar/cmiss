      SUBROUTINE HIMAP(ISEG,ISMAP,STRING,ERROR,*)

C#### Subroutine: HIMAP
C###  Description:
C###    HIMAP hides increment segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER ISEG(*),ISMAP(NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,nomap

      CALL ENTERS('HIMAP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide map
C###  Description:
C###    Hide map projection.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIMAP',ERROR,*9999)
      ELSE
        DO nomap=1,NTMAP
          IF(ISEG(ISMAP(nomap)).EQ.2) THEN

C LKC 2-MAY-1998  use iw
C            CALL ACWK(4,1,ERROR,*9999)
C            CALL VISIB(iw,ISEG,ISMAP(nomap),'INVISIBLE',ERROR,*9999)
C            CALL DAWK(4,1,ERROR,*9999)
            iw=4
            CALL ACWK(iw,1,ERROR,*9999)
            CALL VISIB(iw,ISEG,ISMAP(nomap),'INVISIBLE',ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIMAP')
      RETURN
 9999 CALL ERRORS('HIMAP',ERROR)
      CALL EXITS('HIMAP')
      RETURN 1
      END


