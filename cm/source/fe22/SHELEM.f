      SUBROUTINE SHELEM(ISEG,ISELNO,NEELEM,STRING,ERROR,*)

C#### Subroutine: SHELEM
C###  Description:
C###    SHELEM shows element segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISELNO(NWM,NEM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,noelem,noiw,nr,NTIW

      CALL ENTERS('SHELEM',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show elements
C###  Description:
C###    Make the element number segments visible on the specified
C###    workstation.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which window to show the element number on. The
C###    default is to show the element numbers on all windows.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHELEM',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISEG(ISELNO(iw,ne)).EQ.1) THEN
                  CALL VISIB(iw,ISEG,ISELNO(iw,ne),'VISIBLE',ERROR,
     '              *9999)
                ELSE IF(ISEG(ISELNO(iw,ne)).EQ.0) THEN
                  WRITE(OP_STRING,'('' >>Element number '',I4,'
     '              //''' is not defined on '',I1)') ne,iw
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHELEM')
      RETURN
 9999 CALL ERRORS('SHELEM',ERROR)
      CALL EXITS('SHELEM')
      RETURN 1
      END


