      SUBROUTINE CAMAP(ISEG,ISMAP,STRING,ERROR,*)

C#### Subroutine: CAMAP
C###  Description:
C###    CAMAP cancels map segments.

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

      CALL ENTERS('CAMAP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel map
C###  Description:
C###    Cancel map window.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAMAP',ERROR,*9999)
      ELSE
        IMAP=0
c       CALL COWK(4,ERROR,*9999)

C LKC 2-MAY-1998 Use iw
        iw=4
C        CALL ACWK(4,1,ERROR,*9999)
        CALL ACWK(iw,1,ERROR,*9999)
        DO nomap=1,NTMAP
          IF(ISMAP(nomap).GT.0) THEN
            CALL DELETE_SEGMENT(ISMAP(nomap),ISEG,iw,ERROR,*9999)
          ENDIF
        ENDDO
C LKC 2-MAY-1998 Use iw
C        CALL DAWK(4,1,ERROR,*9999)
        CALL DAWK(iw,1,ERROR,*9999)
      ENDIF

      CALL EXITS('CAMAP')
      RETURN
 9999 CALL ERRORS('CAMAP',ERROR)
      CALL EXITS('CAMAP')
      RETURN 1
      END


