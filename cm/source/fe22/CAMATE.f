      SUBROUTINE CAMATE(ISEG,ISMATE,NEELEM,STRING,ERROR,*)

C#### Subroutine: CAMATE
C###  Description:
C###    CAMATE cancels material segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISMATE(NWM,NEM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,noelem,noiw,nr,NTIW

      CALL ENTERS('CAMATE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel materials;s
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to cancel the
C###    materials on.
C###  Description:
C###    Cancel material segment on specified workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAMAIL',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISMATE(iw,ne).GT.0) THEN
                  CALL DELETE_SEGMENT(ISMATE(iw,ne),ISEG,iw,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CAMATE')
      RETURN
 9999 CALL ERRORS('CAMATE',ERROR)
      CALL EXITS('CAMATE')
      RETURN 1
      END


