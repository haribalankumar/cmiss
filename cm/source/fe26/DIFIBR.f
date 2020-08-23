      SUBROUTINE DIFIBR(IBT,IDO,INP,ISEG,ISFIPR,
     '  NAN,NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,
     '  NRE,NVJE,SE,XA,XE,XG,XP,CSEG,STRING,ERROR,*)

C#### Subroutine: DIFIBR
C###  Description:
C###    DIFIBR displays profiles of fibre angles in 3D geometries.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISEG(*),ISFIPR,NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER CSEG(*)*(*),ERROR*(*),
     '  STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,IW,IXI,N3CO,NOELEM,NOLIST,NR
      REAL*8 RFROMC,XIPOS(3),YMAGN
      LOGICAL CBBREV

      CALL ENTERS('DIFIBR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM display fibre
C###  Parameter:       <in (ELEMENTS/all)[all]>
C###    Specifies the elements into which the fibres are displayed
C###  Parameter:       <at XI_DIRECTION)[1]>
C###    Specifies the alignment of the fibres.
C###  Parameter:       <xi_1 VALUE[0.5]>
C###    Specifies where the fibres are to be displayed in terms of
C###    material coordinates
C###  Parameter:       <xi_2 VALUE[0.5]>
C###    Specifies where the fibres are to be displayed in terms of
C###    material coordinates
C###  Parameter:       <xi_3 VALUE[0.0]>
C###    Specifies where the fibres are to be displayed in terms of
C###    material coordinates
C###  Parameter:       <maximum MAGNITUDE>
C###    Specifies the maximum magnitude of the fibres to be displayed
C###  Description:
C###    Displays the fibres which have been predefined.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<in (ELEMENTS/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<at XI_DIRECTION[1]>'
        OP_STRING(4)=BLANK(1:15)//'<xi_1 VALUE[0.5]>'
        OP_STRING(5)=BLANK(1:15)//'<xi_2 VALUE[0.5]>'
        OP_STRING(6)=BLANK(1:15)//'<xi_3 VALUE[0.0]>'
        OP_STRING(7)=BLANK(1:15)//'<maximum MAGNITUDE>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','DIFIBR',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        IF(CBBREV(CO,'IN',1,noco,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO nr=1,NRT
            DO NOELEM=NELIST(0)+1,NELIST(0)+NEELEM(0,NR)
              NELIST(NOELEM)=NEELEM(NOELEM,NR)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,NR)
          ENDDO
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,*)
     '      ' NELIST=',(NELIST(NOLIST),NOLIST=0,NELIST(0))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
          IXI=IFROMC(CO(N3CO+1))
        ELSE
          IXI=1
        ENDIF
        IF(IXI.EQ.1) THEN
          IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(2)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(2)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(3)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(3)=0.0D0
          ENDIF
        ELSE IF(IXI.EQ.2) THEN
          IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(1)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(1)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(3)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(3)=0.0D0
          ENDIF
        ELSE IF(IXI.EQ.3) THEN
          IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(1)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(1)=0.5D0
          ENDIF
          IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
            XIPOS(2)=RFROMC(CO(N3CO+1))
          ELSE
            XIPOS(2)=0.5D0
          ENDIF
        ENDIF
        IF(CBBREV(CO,'MAXIMUM',1,noco+1,NTCO,N3CO)) THEN
          YMAGN=RFROMC(CO(N3CO+1))
        ELSE
          YMAGN=0.0D0
        ENDIF
        IW=33
        CALL ACWK(IW,0,ERROR,*9999)
        CALL SGFIPR(1,IBT,IDO,INP,ISEG,ISFIPR,IW,IXI,NAN,NBJ,
     '    NELIST,NKJE,NPF,NPNE,NRE,NVJE,
     '    SE,XA,XE,XG,XIPOS,XP,YMAGN,CSEG,ERROR,*9999)
        CALL DAWK(IW,0,ERROR,*9999)
      ENDIF

      CALL EXITS('DIFIBR')
      RETURN
 9999 CALL ERRORS('DIFIBR',ERROR)
      CALL EXITS('DIFIBR')
      RETURN 1
      END


