      SUBROUTINE LIMOTI(IBT,NBH,NEELEM,NHP,NKH,NPNODE,YP,
     '  STRING,FIX,ERROR,*)

C#### Subroutine: LIMOTI
C###  Description:
C###    LIMOTI lists motion parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NPNODE(0:NP_R_M,0:NRM)
      REAL*8 YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nr,nx
      CHARACTER FILE*100
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LIMOTI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list motion<;FILENAME>
C###  Description:
C###    List motion parameters to sreen or FILENAME.opmoti.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<for REGION[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIMOTI',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opmoti','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        nx=1 !temporary
        IF(CBBREV(CO,'FOR',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        CALL ASSERT(KTYP58(nr).GT.0,'>>no motion parameters defined',
     '    ERROR,*9999)

        CALL OPMOTI(IBT,NBH,NEELEM,NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,nr,
     '    FIX(1,1,nx),YP(1,1,nx),ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIMOTI')
      RETURN
 9999 CALL ERRORS('LIMOTI',ERROR)
      CALL EXITS('LIMOTI')
      RETURN 1
      END


