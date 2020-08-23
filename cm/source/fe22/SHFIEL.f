      SUBROUTINE SHFIEL(ISEG,ISFIEL,NEELEM,NELIST,STRING,ERROR,*)

C#### Subroutine: SHFIEL
C###  Description:
C###    SHFIEL shows field segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISFIEL(NWM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),N3CO,ne,noelem,noiw,nolist,nr,NTIW
      LOGICAL CBBREV

      CALL ENTERS('SHFIEL',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show field
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specify which elements to show field in.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the
C###    field on.
C###  Description:
C###    Make the field visible in the specified elements on the
C###    specified workstation.

        OP_STRING(1)=STRING(1:IEND) //' <in (all/ELEMENT#s)[all]>'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHFIEL',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO nr=1,NRT
            DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
              NELIST(noelem)=NEELEM(noelem,nr)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              IF(ISEG(ISFIEL(iw,ne)).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISFIEL(iw,ne),'VISIBLE',ERROR,*9999)
              ELSE IF(ISEG(ISFIEL(iw,ne)).EQ.0) THEN
                WRITE(OP_STRING,
     '            '('' >>Field is not defined in element '',I1)') ne
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHFIEL')
      RETURN
 9999 CALL ERRORS('SHFIEL',ERROR)
      CALL EXITS('SHFIEL')
      RETURN 1
      END


