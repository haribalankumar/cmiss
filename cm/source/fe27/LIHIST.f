      SUBROUTINE LIHIST(NBH,NEELEM,NHE,NHP,NKH,NPLIST,NPNODE,NVHP,
     '  NYNE,NYNP,YP,ZA,ZP,STRING,ERROR,*)

C#### Subroutine: LIHIST
C###  Description:
C###    LIHIST lists history of solution at a node.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,N3CO,nolist,
     '  nonode,np,nr,nrc,nx
      CHARACTER CHAR*4,FILE*100,TYPE*12
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LIHIST',*9999)
      nrc=2 !temporary
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list history<;FILENAME>
C###  Parameter:    <at (NODE#s/all)[all]>
C###    Specifies the nodes to be included in the history file. The
C###    'all' command prompts all the currently defined nodes to
C###    to be defined.
C###  Parameter:    <(value/reaction)[value]>
C###    Gives the option to either list the value of the
C###    dependent variable at the specified node or the reaction at
C###    the specified node.
C###  Description:
C###    Lists history of dependent variables at specified nodes to
C###    screen or file FILENAME.history if qualifier present.
C###    Alternatively reactions are listed.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<at (NODE#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(value/reaction)[value]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIHIST',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG1,IEND1)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'AT',1,noco+1,noco+1,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NPM,NPLIST(0),NPLIST(1),ERROR,*9999)
        ELSE
          NPLIST(0)=0
          DO nr=1,NRT
            DO nonode=NPLIST(0)+1,NPLIST(0)+NPNODE(0,nr)
              NPLIST(nonode)=NPNODE(nonode,nr)
            ENDDO
            NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
          ENDDO
        ENDIF

        IF(CBBREV(CO,'VALUE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='VALUE'
        ELSE IF(CBBREV(CO,'REACTION',1,noco+1,NTCO,N3CO)) THEN
          TYPE='REACTION'
        ELSE
          TYPE='VALUE'
        ENDIF
        nx=1 !Temporary
        nr=1 !temporarily

        DO nolist=1,NPLIST(0)
          np=NPLIST(nolist)
          IF(OPFILE) THEN
C            CHAR=CFROMI(np,'(I4)')
            WRITE(CHAR,'(I4)') np
            CALL STRING_TRIM(CHAR,IBEG2,IEND2)
            IF(nolist.EQ.1) THEN
              CALL OPENF(3,'DISK',FILE(IBEG1:IEND1)//'.iotime','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
            ENDIF
            CALL OPENF(4,'DISK',FILE(IBEG1:IEND1)//'.ophist_'
     '        //CHAR(IBEG2:IEND2),'NEW','SEQUEN','FORMATTED',
     '        132,ERROR,*9999)
          ENDIF

          IF(NYT(nrc,1,nx).GT.0.AND.ITYP5(nr,nx).EQ.2) THEN !time-dependent
            CALL STRING_TRIM(FILE02,IBEG,IEND)
C!!! same unit as IOOUT for `set output'!
            CALL OPENF(9,'DISK',FILE02(IBEG:IEND)//'.history','OLD',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
            CALL OPHIST(NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '        nolist,np,NPNODE,nr,NVHP(1,1,1,nr),NYNE,NYNP,YP,ZA,ZP,
     '        OPFILE,TYPE,ERROR,*9999)
          ENDIF
        ENDDO

        IF(OPFILE) THEN
          CALL CLOSEF(3,ERROR,*9999)
          CALL CLOSEF(4,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('LIHIST')
      RETURN
 9999 CALL ERRORS('LIHIST',ERROR)
      CALL EXITS('LIHIST')
      RETURN 1
      END


