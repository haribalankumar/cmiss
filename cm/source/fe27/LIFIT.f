      SUBROUTINE LIFIT(NBJ,NEELEM,NKH,NMNO,NPNODE,NRLIST,NVHP,NXLIST,
     '  NYNP,WU,STRING,FIX,ERROR,*)

C#### Subroutine: LIFIT
C###  Description:
C###    LIFIT lists fitting information.

C**** KTYP8=1 ITYP6(nr,nx) =1 is linear    geometric fitting problem
C**** KTYP8=1 ITYP6(nr,nx) =2 is nonlinear geometric fitting problem
C**** KTYP8=2 NJ_LOC(njl_fibr,0,nr) =1,2 or 3 is
C****          fibre or sheet fitting problem
C****         (ITYP6(nr,nx)=1)
C**** KTYP8=3 NJ_LOC(njl_fibr,0,nr) > 0    is field fitting problem
C****         (ITYP6(nr,nx)=1)
C**** KTYP8=4 is potential field fitting problem
C****         (ITYP6(nr,nx)=1)
C**** KTYP8=5 is motion fitting problem with Fourier basis
C**** KTYP8=6 is data fitting by optimisation
C**** KTYP6=1 if fitting uses Gauss points
C**** KTYP12=1 for smoothing constraints (otherwise 0) defined in WU

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKH(NHM,NPM,NCM,0:NRM),NMNO(1:2,0:NOPM,NXM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 WU(0:NUM+1,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,no_nrlist,nr,nx,nxc
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,OPFILE

      CALL ENTERS('LIFIT',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list fit<;FILENAME>
C###  Description:
C###    Lists fitting parameters. The information is written to the
C###    file FILENAME (with extension .opfit)
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIFIT',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opfit','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL ASSERT(KTYP8.GT.0,'>>Define fit first',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
C CPB 8/6/94 Adding NX_LOC
        IF(KTYP8.EQ.6) THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx.GT.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)
        ELSE
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '      ERROR,*9999)
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL OPFIT(NBJ,NEELEM,NKH(1,1,1,nr),NMNO(1,0,nx),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNP,WU,
     '      FIX(1,1,nx),ERROR,*9999)
        ENDDO !no_nrlist

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIFIT')
      RETURN
 9999 CALL ERRORS('LIFIT',ERROR)
      CALL EXITS('LIFIT')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END


