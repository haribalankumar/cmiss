      SUBROUTINE OPINIT(ITHRES,NBH,NBJ,NEELEM,NHE,NHP,
     '  NKH,NODENVCB,NPNODE,nr,NVCB,NVHP,NW,NWQ,nx,nxc,NYNE,
     '  NYNP,AQ,THRES,YG,YP,YQ,FIX,FULL,ERROR,*)
C SMAR009 19/01/99 removed NVCNODE,
C#### Subroutine: OPINIT
C###  Description:
C###    OPINIT outputs initial and boundary conditions for class
C###    nx, region nr.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      INTEGER ITHRES(3,NGM,NEM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M),NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NODENVCB(NVCBM),NPNODE(0:NP_R_M),nr,NVCB(-1:3,NVCBM),
     '  NVHP(NHM,NPM,NCM),
     '  NW(NEM,3),NWQ(8,0:NQM,NAM),nx,nxc,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
C SMAR009 19/01/99 removed NVCNODE(2,NP_R_M),
      REAL*8 AQ(NMAQM,NQM),THRES(3,NGM,NEM),YG(NIYGM,NGM,NEM),
     '  YP(NYM,NIYM),YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),FULL
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('OPINIT',*9999)

      WRITE(OP_STRING,'(/'' Class '',I1,'' (nx='',I1,'
     '  //''') Region '',I1,'':'')') nxc,nx,nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

! Time-varying boundary conditions
      IF(ITYP5(nr,nx).EQ.2) THEN !Time integration problem
        IF(KTYP3_init(nx).EQ.1) THEN
          WRITE(OP_STRING,'('' Bdry conditions are constant in time'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_init(nx).EQ.2) THEN
          WRITE(OP_STRING,'('' Time-varying bdry conditions are defined'
     '      //' in USER_IPINIT subroutine'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_init(nx).EQ.3) THEN
          CALL STRING_TRIM(FILE03,IBEG,IEND)
          WRITE(OP_STRING,'('' Time-varying bdry conditions are defined'
     '      //' in '//FILE03(IBEG:IEND)//'.IPINIT_time'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_init(nx).EQ.4) THEN
          !using time variables
          WRITE(OP_STRING,'('' Time-varying bdry conditions are defined'
     '      //' by time variables'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_init(nx).EQ.5) THEN
          !using time variables
          WRITE(OP_STRING,'('' Time-varying bdry conditions are defined'
     '      //' by cellular process models'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF !KTYP3_init
      ENDIF !ityp5(nr)

      IF(ITYP1(nr,nx).EQ.3) THEN      !FE30 problems
        CALL OPINI3(ITHRES,NBH,NBJ,NEELEM,NHE,NHP,
     '    NKH,NPNODE,nr,NVHP,NWQ,nx,NYNP,
     '    FIX,FULL,AQ,THRES,YG,YP,YQ,ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.4) THEN !FE40 problems
        CALL OPINI4(NHP,NKH,NPNODE,nr,NVHP,nx,NYNP,YP,FIX,ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.5) THEN !FE50 problems
        CALL OPINI5(NBH,NEELEM,NHP,NKH,NPNODE,nr,
     '    NVHP,NW,nx,NYNE,NYNP,FIX,FULL,YP,ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.6) THEN !FE60 problems
        CALL OPINI6(NODENVCB,NPNODE,nr,NVCB,nx,NYNP,YP,
     '    ERROR,*9999)
C  SMAR009 18/01/99 removed NHP,NVCNODE,
      ELSE IF(ITYP1(nr,nx).EQ.9) THEN !FE90 problems
!news AJP 12/4/95
        IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.1) THEN
          !Static linear elasticity
          CALL OPINI4(NHP,NKH,NPNODE,nr,NVHP,nx,NYNP,YP,FIX,ERROR,*9999)
        ELSE !Rest of BEM
          CALL OPINI9(NHP,NKH,NPNODE,nr,NVHP,nx,NYNP,FIX,FULL,YP,
     '      ERROR,*9999)
        ENDIF
!newe
      ENDIF

      CALL EXITS('OPINIT')
      RETURN
 9999 CALL ERRORS('OPINIT',ERROR)
      CALL EXITS('OPINIT')
      RETURN 1
      END


