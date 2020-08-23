      SUBROUTINE DEBOUN(NBJ,NEELEM,NENP,NNB,NPNE,NRLIST,NXI,NYNP,
     '  YP,FIX,STRING,ERROR,*)

C#### Subroutine: DEBOUN
C###  Description:
C###    DEBOUN defines boundary conditions with prompted input or from
C###    filename.ipboun (or filename.irboun for coupled problems).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),NNB(4,4,4,NBFM),
     '  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,
     '  IEND2,IPFILE,no_nrlist,nr,nx
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEBOUN',*9999)

      nx=1 ! temporary cpb 22/11/94

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define boundary;d/g/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define boundary properties. Used for setting up simple boundary
C###    conditions on blocks of nodes. Boundary properties are
C###    read from or written to the file FILENAME (with extension
C###    .ipinit) in the directory specified by PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';d/g/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEINIT',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 6-Jun-1993
        CALL PARSE_QUALIFIERS('DGLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(CALL_EQUA,'>>Problem type not defined',ERROR,*9999)
        CALL ASSERT(CALL_MATE,'>>Material parameters not defined',
     '    ERROR,*9999)
        CALL ASSERT(CALL_INIT,'>>Initial conditions not defined',
     '    ERROR,*9999)

       IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'init',
     '        STATUS,ERR,ERROR,*9999)

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

              IF(ITYP1(nr,nx).EQ.3) THEN
C
C ??? cpb 22/1/94 Why are ny's being used here ?
C
                CALL IPBOUN(NBJ,NEELEM,NENP,NNB,NPNE,nr,NXI,NYNP,
     '            YP(1,1,nx),ALL_REGIONS,FIX(1,1,nx),ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'('' DEFINE BOUNDARY only '
     '            //'implemented for FE30 problems'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF

            ENDDO
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO

        ENDIF
      ENDIF

      CALL EXITS('DEBOUN')
      RETURN
 9999 CALL ERRORS('DEBOUN',ERROR)
      CALL EXITS('DEBOUN')
      RETURN 1
      END


