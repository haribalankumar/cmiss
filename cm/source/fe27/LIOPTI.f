      SUBROUTINE LIOPTI(LDR,NLNO,NMNO,NP1OPT,NP2OPT,NP3OPT,
     '  NRLIST,NXLIST,NYNO,PAOPTY,CM,CONTR,PAOPTI,PMIN,PMAX,
     '  RESID,RESIDM,STRING,ERROR,*)

C#### Subroutine: LIOPTI
C###  Description:
C###    LIOPTI lists optimisation parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER LDR(0:NDM),NLNO(NOPM,NXM),NMNO(1:2,0:NOPM,NXM),
     '  NP1OPT(NOPM),NP2OPT(NOPM),NP3OPT(NOPM),NRLIST(0:NRM),
     '  NXLIST(0:NXM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),PAOPTY(NOPM)
      REAL*8 CM(*),CONTR(*),PAOPTI(*),PMIN(*),PMAX(*),
     '  RESID(*),RESIDM(*)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,no_nrlist,nr,nx,nxc
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,CBBREV,OPFILE,SUMMARY

      CALL ENTERS('LIOPTI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list optimise<;FILENAME>
C###  Parameter:    <summary>
C###    Only output a summary
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to list the optimisation parameters
C###    for. The all value specifies all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number to list the optimisation parameters
C###    for.
C###  Description:
C###    Lists optimization parameters to screen or file FILENAME.opopti
C###    if the qualifier is present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<summary>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIOPTI',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opopti','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL ASSERT(KTYP26.GT.0.AND.KTYP27.GT.0,
     '    '>>No optimisation problem defined',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
        CALL ASSERT(nx.GT.0,
     '    '>>No nx defined for this optimisation class',ERROR,*9999)

        IF(CBBREV(CO,'SUMMARY',1,noco+1,NTCO,N3CO)) THEN
          SUMMARY=.TRUE.
        ELSE
          SUMMARY=.FALSE.
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL OPOPTI(LDR,NLNO(1,nx),NMNO(1,0,nx),NP1OPT,NP2OPT,NP3OPT,
     '      NYNO(0,1,1,nr,nx),PAOPTY,CM,CONTR,PAOPTI,PMIN,PMAX,
     '      RESID,RESIDM,SUMMARY,ERROR,*9999)
        ENDDO !no_nrlist

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIOPTI')
      RETURN
 9999 CALL ERRORS('LIOPTI',ERROR)
      CALL EXITS('LIOPTI')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


