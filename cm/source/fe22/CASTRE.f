      SUBROUTINE CASTRE(ISEG,ISSTRE,NEELEM,STRING,ERROR,*)

C#### Subroutine: CASTRE
C###  Description:
C###    CASTRE cancels element stress segments.

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

      CALL ENTERS('CASTRE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel stress;s
C###  Description:
C###    Cancel stress segment.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CASTRE',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        iw=2*NJT-3
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nostre=1,NTSTRE
                IF(ISSTRE(ne,nostre).GT.0) THEN
                  CALL DELETE_SEGMENT(ISSTRE(ne,nostre),ISEG,iw,ERROR,
     '              *9999)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
          NTSTRE=0
        ENDIF
      ENDIF

      CALL EXITS('CASTRE')
      RETURN
 9999 CALL ERRORS('CASTRE',ERROR)
      CALL EXITS('CASTRE')
      RETURN 1
      END


