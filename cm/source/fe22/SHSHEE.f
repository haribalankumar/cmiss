      SUBROUTINE SHSHEE(ISEG,ISSHEE,NEELEM,NELIST,STRING,ERROR,*)

C#### Subroutine: SHSHEE
C###  Description:
C###    SHSHEE shows sheet segments.
C**** Valid commands:
C**** Show sheets in list1 at list2
C****   where list1 are element numbers (defaults to all elements)
C****      "  list2 are sheet   indices (defaults to noshee=1)

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSHEE(NWM,NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IL2(100),iw,IWK(6),N3CO,ne,noelem,noil2,noiw,
     '  nolist,noshee,nr,NTIL2,NTIW
      LOGICAL CBBREV

      CALL ENTERS('SHSHEE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show sheets
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specifies the element(s) in which to show the sheets. The
C###    all keyword specifies all currently defined elements.
C###  Parameter:    <at SHEET_INDICES[1]>
C###    Specify the sheet indices
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###    Make the sheet segments visible in the specified elements on
C###    the specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <in (all/ELEMENT#s)[all]>'
        OP_STRING(2)=BLANK(1:15) //'<at SHEET_INDICES[1]>'
        OP_STRING(3)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHSHEE',ERROR,*9999)
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
        IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),100,NTIL2,IL2,ERROR,*9999)
        ELSE
          NTIL2=1
          IL2(1)=NTSHEE
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
C ***       loop over sheet indices
            DO noil2=1,NTIL2
              noshee=IL2(noil2)
C ***         loop over elements
              DO nolist=1,NELIST(0)
                ne=NELIST(nolist)
                IF(ISEG(ISSHEE(iw,ne,noshee)).EQ.1) THEN
                  CALL VISIB(iw,ISEG,ISSHEE(iw,ne,noshee),'VISIBLE',
     '              ERROR,*9999)
                ELSE IF(ISEG(ISSHEE(iw,ne,noshee)).EQ.0) THEN
                  WRITE(OP_STRING,
     '              '('' >>Sheet in element '',I4,'' at '',I4,'
     '              //''' is not defined on '',I1)') ne,noshee,iw
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHSHEE')
      RETURN
 9999 CALL ERRORS('SHSHEE',ERROR)
      CALL EXITS('SHSHEE')
      RETURN 1
      END


