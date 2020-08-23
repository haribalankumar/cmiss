      SUBROUTINE OPMESH(NBJ,NEELEM,NELIST,NENP,NORD,NPLIST,NPNE,nr,NVJE,
     &  nx,NXI,NYNP,CE,CP,XP,YP,LISOLN,ERROR,*)

C#### Subroutine: OPMESH
C###  Description:
C###    OPMESh outputs specialized mesh data.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),
     &  NENP(NPM,0:NEPM),NORD(5,NE_R_M),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),nr,NVJE(NNM,NBFM,NJM,NEM),nx,
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 CE(NMM,NEM),CP(NMM,NPM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      LOGICAL LISOLN
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER N3CO
      INTEGER*4 BRANCHES_PTR,diameters_PTR,lengths_PTR,NBRANCHES_PTR,
     &  STATISTICS_PTR
      CHARACTER OPT1(10,2)*25
      LOGICAL CBBREV

      DATA OPT1(1,1)/'Regular'/
     '     OPT1(1,2)/'Regular'/,
     '     OPT1(2,1)/'Circular'/,
     '     OPT1(2,2)/'Eccentric spheres mesh'/,
     '     OPT1(3,1)/'Regular fractal tree'/,
     '     OPT1(3,2)/'Regular fractal tree'/,
     '     OPT1(4,1)/'Stochastic fractal tree'/,
     '     OPT1(4,2)/'Stochastic fractal tree'/,
     '     OPT1(5,1)/' '/,
     '     OPT1(5,2)/'Open cylindrical mesh'/,
     '     OPT1(6,1)/'Coronary mesh '/,
     '     OPT1(6,2)/'Coronary mesh'/,
     '     OPT1(7,1)/' '/,
     '     OPT1(7,2)/' '/,
     '     OPT1(8,1)/'Purkinje fibre mesh'/,
     '     OPT1(8,2)/'Purkinje fibre mesh'/,
     '     OPT1(9,1)/'Regular mesh with incision '/,
     '     OPT1(9,2)/'Regular mesh with incision '/,
     '     OPT1(10,1)/' '/,
     '     OPT1(10,2)/' '/

      CALL ENTERS('OPMESH',*9999)
      CALL ASSERT(JTYP14.GE.1.AND.JTYP14.LE.5,'>>Invalid JTYP14',
     '  ERROR,*9999)
      IF(CBBREV(CO,'COORDINATES',4,noco+1,NTCO,N3CO).OR.LISOLN) THEN
        CALL ASSERT(JTYP14.EQ.3.OR.JTYP14.EQ.4,'>>Incorrect mesh type',
     '    ERROR,*9999)
        NBRANCHES_PTR=0
        BRANCHES_PTR=0
        STATISTICS_PTR=0
        diameters_PTR=0
        lengths_PTR=0
        CALL ALLOCATE_MEMORY(5*NE_R_M,1,INTTYPE,NBRANCHES_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(10*NE_R_M,1,DPTYPE,BRANCHES_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(21*NE_R_M,1,DPTYPE,STATISTICS_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NEM,1,DPTYPE,diameters_PTR,.TRUE.,ERROR,
     &    *9999)
        CALL ALLOCATE_MEMORY(NEM,1,DPTYPE,lengths_PTR,.TRUE.,ERROR,
     &    *9999)
        CALL OPMESH2(NBJ,%VAL(NBRANCHES_PTR),NEELEM,NELIST,NENP,NORD,
     &    NPLIST,NPNE,nr,NVJE,nx,NXI,NYNP,%VAL(BRANCHES_PTR),CP,
     &    %VAL(diameters_PTR),%VAL(lengths_PTR),%VAL(STATISTICS_PTR),XP,
     &    YP,LISOLN,.TRUE.,ERROR,*9999)
        CALL FREE_MEMORY(NBRANCHES_PTR,ERROR,*9999)
        CALL FREE_MEMORY(BRANCHES_PTR,ERROR,*9999)
        CALL FREE_MEMORY(STATISTICS_PTR,ERROR,*9999)
        CALL FREE_MEMORY(diameters_PTR,ERROR,*9999)
        CALL FREE_MEMORY(lengths_PTR,ERROR,*9999)
      ELSE
        IF(.NOT.LISOLN)THEN
          WRITE(OP_STRING,'(''  Specialized mesh is '',A)')
     '      OPT1(JTYP14,NJT-1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(JTYP14.EQ.1) THEN !Regular mesh
          CALL OPMESH1(ERROR,*9999)
        ELSE IF(JTYP14.EQ.3.OR.JTYP14.EQ.4) THEN !Fractal tree
          NBRANCHES_PTR=0
          BRANCHES_PTR=0
          STATISTICS_PTR=0
          diameters_PTR=0
          lengths_PTR=0
          CALL ALLOCATE_MEMORY(5*NE_R_M,1,INTTYPE,NBRANCHES_PTR,.TRUE.,
     &      ERROR,*9999)
          CALL ALLOCATE_MEMORY(10*NE_R_M,1,DPTYPE,BRANCHES_PTR,.TRUE.,
     &      ERROR,*9999)
          CALL ALLOCATE_MEMORY(21*NE_R_M,1,DPTYPE,STATISTICS_PTR,.TRUE.,
     &      ERROR,*9999)
          CALL ALLOCATE_MEMORY(NEM,1,DPTYPE,diameters_PTR,.TRUE.,ERROR,
     &      *9999)
          CALL ALLOCATE_MEMORY(NEM,1,DPTYPE,lengths_PTR,.TRUE.,ERROR,
     &      *9999)
          CALL OPMESH2(NBJ,%VAL(NBRANCHES_PTR),NEELEM,NELIST,NENP,NORD,
     &      NPLIST,NPNE,nr,NVJE,nx,NXI,NYNP,%VAL(BRANCHES_PTR),CP,
     &      %VAL(diameters_PTR),%VAL(lengths_PTR),%VAL(STATISTICS_PTR),
     &      XP,YP,LISOLN,.FALSE.,ERROR,*9999)
          CALL FREE_MEMORY(NBRANCHES_PTR,ERROR,*9999)
          CALL FREE_MEMORY(BRANCHES_PTR,ERROR,*9999)
          CALL FREE_MEMORY(STATISTICS_PTR,ERROR,*9999)
          CALL FREE_MEMORY(diameters_PTR,ERROR,*9999)
          CALL FREE_MEMORY(lengths_PTR,ERROR,*9999)
        ELSE IF(JTYP14.EQ.2) THEN !Eccentric spheres mesh
          CALL OPMESH3(ERROR,*9999)
        ELSE IF(JTYP14.EQ.5) THEN !Eccentric spheres mesh
          CALL OPMESH5(ERROR,*9999)
        ELSE IF(JTYP14.EQ.6) THEN !Coronary mesh
          CALL OPMESH6(ERROR,*9999)
        ELSE IF(JTYP14.EQ.8) THEN !Coronary mesh
          CALL OPMESH8(ERROR,*9999)
        ELSE IF(JTYP14.EQ.9) THEN !Coronary mesh
          CALL OPMESH9(ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('OPMESH')
      RETURN
 9999 CALL ERRORS('OPMESH',ERROR)
      CALL EXITS('OPMESH')
      RETURN 1
      END


