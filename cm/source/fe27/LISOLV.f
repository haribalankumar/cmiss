      SUBROUTINE LISOLV(NPNY,NRLIST,NXLIST,STRING,ERROR,*)

C#### Subroutine: LISOLV
C###  Description:
C###    LISOLV lists solution parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NPNY(0:6,NYM,0:NRCM,NXM),
     '  NRLIST(0:NRM),NXLIST(0:NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nolist,nr,nx,nxc
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,OPFILE


C LKC 25-NOV-97 INTEGER NPNODE(0:NP_R_M,0:NRM),

      CALL ENTERS('LISOLV',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list solve<;FILENAME>
C###  Description:
C###    List the details of the solution procedure and the level of
C###    solver output to the screen or to file FILENAME.opsolv if
C###    qualifier present.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LISOLV',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opsolv','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL ASSERT(CALL_SOLV,'>>Define solve first',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        DO nolist=1,NRLIST(0)
          nr=NRLIST(nolist)
          IF(ITYP2(nr,nx).GT.0) THEN
            CALL OPSOLV(NPNY(0,1,0,nx),nr,nx,nxc,
     '        ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' >>Problem type not defined for '
     '        //'region '',I2)') nr
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LISOLV')
      RETURN
 9999 CALL ERRORS('LISOLV',ERROR)
      CALL EXITS('LISOLV')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


