      SUBROUTINE SHMATE(ISEG,ISMATE,NEELEM,STRING,ERROR,*)

C#### Subroutine: SHMATE
C###  Description:
C###    SHMATE shows material segments.

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

      CALL ENTERS('SHMATE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show materials
C###  Parameter:   <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the
C###    materials on.
C###  Description:
C###    Make the material segments visible on the specified workstation.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHMATE',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISEG(ISMATE(iw,ne)).EQ.1) THEN
                  CALL VISIB(iw,ISEG,ISMATE(iw,ne),'VISIBLE',ERROR,
     '              *9999)
                ELSE IF(ISEG(ISMATE(iw,ne)).EQ.0) THEN
                  WRITE(OP_STRING,'('' >>Material in element '',I4,'
     '              //''' is not defined on '',I1)') ne,iw
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHMATE')
      RETURN
 9999 CALL ERRORS('SHMATE',ERROR)
      CALL EXITS('SHMATE')
      RETURN 1
      END


