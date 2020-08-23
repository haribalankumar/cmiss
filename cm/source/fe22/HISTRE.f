      SUBROUTINE HISTRE(ISEG,ISSTRE,NEELEM,STRING,ERROR,*)

C#### Subroutine: HISTRE
C###  Description:
C###    HISTRE hides element stress segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSTRE(NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,ne,noelem,nostre,nr

      CALL ENTERS('HISTRE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide stress
C###  Description:
C###    Hide principal stress vectors.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HISTRE',ERROR,*9999)
      ELSE
        iw=2*NJT-3
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nostre=1,NTSTRE
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISEG(ISSTRE(ne,nostre)).EQ.2) THEN
                  CALL VISIB(iw,ISEG,ISSTRE(ne,nostre),'INVISIBLE',
     '              ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('HISTRE')
      RETURN
 9999 CALL ERRORS('HISTRE',ERROR)
      CALL EXITS('HISTRE')
      RETURN 1
      END


