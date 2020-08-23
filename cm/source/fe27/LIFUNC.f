      SUBROUTINE LIFUNC(IBT,IDO,INP,LD,NBH,NBJ,NEELEM,NHE,NHP,
     '  NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NVHE,NVHP,NVJE,NW,NYNE,NYNP,
     '  CURVCORRECT,SE,XA,XE,XID,XP,YP,ZA,ZE,ZP,STRING,ERROR,*)

C#### Subroutine: LIFUNC
C###  Description:
C###    LIFUNC lists objective function.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,ID_TYPE,IEND,IEND1,IEND2,IFROMC,N3CO,
     '  noobje,NO_OBJECT,nr,nx
      CHARACTER FILE*100,OBJECT*20,TYPE*9
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LIFUNC',*9999)

      nx=1 !temporary
      nr=1 !update later

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        IF(NJ_LOC(njl_fibr,0,nr).GT.0) THEN !field variables defined
          TYPE='field'
        ELSE IF(ITYP2(nr,nx).GT.0) THEN !dependent variables defined
          TYPE='dependent'
        ELSE
C KAT 29May99: objective does nothing
C          TYPE='objective'
C LKC 6-NOV-2000 Zero-length string
C          TYPE=''
          TYPE='-'
        ENDIF
        CALL STRING_TRIM(TYPE,IBEG1,IEND1)
        IF(NTOBJE.GT.0) THEN
          OBJECT=OBJECT_NAME(NTOBJE)
        ELSE
          OBJECT=' '
        ENDIF
        CALL STRING_TRIM(OBJECT,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM list function<;FILENAME>
C###  Description:
C###    Lists the geometric postion and a function value with their first
C###    derivatives at each data point in a graphics object.
C###  Parameter:    <(field/dependent) VARIABLE#[1]>
C###    Use `field' to specify the function as the fibre angle
C###    VARIABLE#.  Use dependent to specify the function as the
C###    dependent variable VARIABLE#.  If fibre angles are defined the
C###    default is `field 1'.   Otherwise the default is `dependent 1'
C###  Parameter:    <with OBJECT_NAME>
C###    Specify the graphics object name.  The default is the last
C###    object drawn.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)
     '    //'<(field/dependent)['//TYPE(IBEG1:IEND1)
     '    //'] ID#[1]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<with OBJECT_NAME['//OBJECT(1:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIFUNC',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opfunc','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'FIELD',1,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
          ID_TYPE=IFROMC(CO(N3CO+1))
        ELSE IF(CBBREV(CO,'DEPENDENT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='DEPENDENT'
          ID_TYPE=IFROMC(CO(N3CO+1))
C KAT 29May99: objective does nothing
C        ELSE IF(CBBREV(CO,'OBJECTIVE',1,noco+1,NTCO,N3CO)) THEN
C          TYPE='OBJECTIVE'
C          ID_TYPE=1
        ELSE
          IF(NJ_LOC(njl_fibr,0,nr).GT.0) THEN !field variables defined
            TYPE='FIELD'
          ELSE IF(ITYP2(nr,nx).GT.0) THEN !dependent vars defined
            TYPE='DEPENDENT'
          ELSE
C KAT 29May99: objective does nothing
C            TYPE='OBJECTIVE'
            ERROR='>> Function not defined'
            GO TO 9999
          ENDIF
          ID_TYPE=1
        ENDIF

        IF(TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:9).EQ.'DEPENDENT') THEN
          CALL ASSERT(NTOBJE.GT.0,'>>no object defined',ERROR,*9999)
          IF(CBBREV(CO,'OBJECT',1,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            OBJECT=CO(N3CO+1)(IBEG:IEND)
          ELSE
            CALL STRING_TRIM(OBJECT_NAME(NTOBJE),IBEG,IEND)
            OBJECT=OBJECT_NAME(NTOBJE)(IBEG:IEND)
          ENDIF
          CALL STRING_TRIM(OBJECT,IBEG,IEND)
          DO noobje=1,NTOBJE
            IF(OBJECT_NAME(noobje)(IBEG:IEND).EQ.OBJECT(IBEG:IEND)) THEN
              NO_OBJECT=noobje
            ENDIF
          ENDDO
        ENDIF

        CALL OPFUNC(IBT,IDO,ID_TYPE,INP,LD,NBH,NBJ,NEELEM,NHE(1,nx),
     '    NHP(1,nr,nx),NKH(1,1,1,nr),NKHE,NKJE,NO_OBJECT,NPF,NPNE,
     '    NPNODE,nr,NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),NYNE,NYNP,
     '    CURVCORRECT,SE,XA,XE,XID,XP,YP,ZA,ZE,ZP,TYPE,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIFUNC')
      RETURN
 9999 CALL ERRORS('LIFUNC',ERROR)
      CALL EXITS('LIFUNC')
      RETURN 1
      END


