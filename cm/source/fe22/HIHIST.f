      SUBROUTINE HIHIST(ISEG,ISHIST,NPNODE,STRING,ERROR,*)

C#### Subroutine: HIHIST
C###  Description:
C###    HIHIST hides history segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISHIST(0:NPM),NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,nonode,np,nr

      CALL ENTERS('HIHIST',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide history
C###  Description:
C###    Hide history plot.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIHIST',ERROR,*9999)
      ELSE

C 2-MAY-1998 use iw
C        IF(IWKS(10).GT.0) THEN
C          CALL ACWK(10,1,ERROR,*9999)
        iw=10
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nr=1,NRT
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(ISEG(ISHIST(np)).EQ.2) THEN
                CALL VISIB(iw,ISEG,ISHIST(np),'INVISIBLE',ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
C 2-MAY-1998 use iw
C          CALL DAWK(10,1,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('HIHIST')
      RETURN
 9999 CALL ERRORS('HIHIST',ERROR)
      CALL EXITS('HIHIST')
      RETURN 1
      END


