      SUBROUTINE SHFIBR(ISEG,ISFIBR,NEELEM,NELIST,STRING,ERROR,*)

C#### Subroutine: SHFIBR
C###  Description:
C###    SHFIBR shows fibre segments.
C**** Valid commands:
C**** Show fibres in list1 at list2
C****   where list1 are element numbers (defaults to all elements)
C****      "  list2 are fibre   indices (defaults to nofibr=1)

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISFIBR(NWM,NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IL2(100),iw,IWK(6),N3CO,ne,noelem,nofibr,
     '  nolist,noil2,noiw,nr,NTIL2,NTIW
      LOGICAL CBBREV

      CALL ENTERS('SHFIBR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show fibres
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specify which elements to show fibres in.
C###  Parameter:    <at FIBRE_INDICES[1]>
C###    Specify the fibre indices
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the
C###    fibres on.
C###  Description:
C###    Make the specified fibre segments in the specified elements
C###    visible on the specified worksation.

        OP_STRING(1)=STRING(1:IEND) //' <in (all/ELEMENT#s)[all]>'
        OP_STRING(2)=BLANK(1:15) //'<at FIBRE_INDICES[1]>'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHFIBR',ERROR,*9999)
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
          IL2(1)=NTFIBR
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
C ***       loop over fibre indices
            DO noil2=1,NTIL2
              nofibr=IL2(noil2)
C ***         loop over elements
              DO nolist=1,NELIST(0)
                ne=NELIST(nolist)
                IF(ISEG(ISFIBR(iw,ne,nofibr)).EQ.1) THEN
                  CALL VISIB(iw,ISEG,ISFIBR(iw,ne,nofibr),'VISIBLE',
     '              ERROR,*9999)
                ELSE IF(ISEG(ISFIBR(iw,ne,nofibr)).EQ.0) THEN
                  WRITE(OP_STRING,'('' >>Fibre in element '',I4,'
     '              //''' at '',I4,'' is not defined on '',I1)')
     '              ne,nofibr,iw
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHFIBR')
      RETURN
 9999 CALL ERRORS('SHFIBR',ERROR)
      CALL EXITS('SHFIBR')
      RETURN 1
      END


