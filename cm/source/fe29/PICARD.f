      SUBROUTINE PICARD(CONSTR,M,N,NRHS,U,LDU,SM,LDSM,B,LDB,
     '  REG_PARAMETER,WORK,LWORK,ISEG,ISPLOTXY,CSEG,ERROR,*)

C#### Subroutine: PICARD
C###  Description:
C###    Visual inspection of the Picard condition.
C###  Reference:
C###   P. C. Hansen, "Regularization tools, a Matlab package for
C###  analysis and solution of discrete ill-posed problems," UNI.C, 1998.
C###  Note:
C###    LWORK >= 4*min(M,N)
CC JMB 13-OCT-2000

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'plot02.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER ISEG(*), ISPLOTXY(2), LDB, LDSM, LDU, LWORK, M, N, NRHS
      REAL*8 B(LDB,*), REG_PARAMETER(*), SM(LDSM,*), U(LDU,*), WORK(*)
      CHARACTER CONSTR, CSEG(*)*(*), ERROR*(*)
!     Local Variables
      INTEGER i, IBEG, ICHAR, IEND, INDEX, INFO, iw, j, MN, NOQUES
      INTEGER*4 D_PTS_PTR
      CHARACTER CHAR*5
      LOGICAL FILEIP
      EXTERNAL POLYLINE_DYNAM
      PARAMETER (iw = 10)
!     Functions
      INTEGER INDEX_POLYLINE

      CALL ENTERS('PICARD',*9999)

      ! Initialisation
      MN = MIN(M,N)
      IF( LWORK.LT.4*MN ) GOTO 9999
      CALL GSVALUES(CONSTR, MN, SM, LDSM, WORK(1), ERROR, *9999)

      DO i = 1,MN
          WORK(2*MN + i) = DBLE(i)
      ENDDO
      D_PTS_PTR = 0
      CALL ALLOCATE_MEMORY(3*MN, 6, DPTYPE, D_PTS_PTR, MEM_INIT,
     '  ERROR,*9999)
      DO j = 1,NRHS
        CALL FOURIERCOEFFS(CONSTR, M, N, U, LDU, B(1,j),
     '    WORK(MN + 1), ERROR, *9999)
        DO i = 1,MN
          WORK(3*MN + i) = DABS(WORK(MN + i))
        ENDDO
        ! Draw results
        CALL ACWK(iw, 0, ERROR, *9999)
        XSCALE = 'LINEAR'
        YSCALE = 'LOG'
        NEXTPLOT = 'REPLACE'
        ISHOLD = .TRUE.
        XMAX = MN
        XMIN = 1
        YMAX = WORK(1)
        YMIN = WORK(MN)
        DO i = 1,MN
          YMAX = MAX(YMAX,WORK(3*MN + i))
          YMIN = MIN(YMIN,WORK(3*MN + i))
        ENDDO
        INDEX = INDEX_POLYLINE(0, 'DOTTED', 'WIDTH1', 'BLUE')
        CALL ACWK(iw, 0, ERROR, *9999)
        CALL SGPLOTXY(POLYLINE_DYNAM, INDEX, ISEG, 1, ISPLOTXY(1), iw,
     '    MN, %VAL(D_PTS_PTR), WORK(2*MN + 1), WORK(1), CSEG,
     '    ERROR, *9999)
        NEXTPLOT = 'ADD'
        ISHOLD = .FALSE.
        INDEX = INDEX_POLYLINE(0, 'SOLID', 'WIDTH1', 'RED')
        CALL SGPLOTXY(POLYLINE_DYNAM, INDEX, ISEG, 1, ISPLOTXY(1), iw,
     '    MN, %VAL(D_PTS_PTR), WORK(2*MN + 1),WORK(3*MN + 1), CSEG,
     '    ERROR, *9999)
        TITLE = 'Picard Condition'
        XLABEL = 'n'
C LKC 6-NOV-2000 Zero length string
C        YLABEL = ''
        YLABEL = '-'
        XGRID = .TRUE.
        YGRID = .TRUE.
        CALL SGPLOTXY(%VAL(0), INDEX, ISEG, 2, ISPLOTXY(2), iw, %VAL(0),
     '    %VAL(D_PTS_PTR), %VAL(0), %VAL(0), CSEG, ERROR, *9999)
        CALL DAWK(iw, 0, ERROR, *9999)
        CALL REFRESH_GRAPHICS(0.0, ERROR, *9999)

        ! Determine truncation rank
        WRITE(CHAR,'(I5)') MN
        CALL STRING_TRIM(CHAR, IBEG, IEND)
        FORMAT = '('' Enter the truncation rank ['//CHAR(IBEG:IEND)//']:
     '    '',I3)'
        NOQUES = 0
        FILEIP = .FALSE.
        ICHAR = 999
        IDEFLT(1) = MN
        IOTYPE = 0
        CALL GINOUT(IOTYPE, 3, IVDU, IFILE, 0, 0, NOQUES, FILEIP,
     '    FORMAT, 1, ADATA, ADEFLT, CDATA, CDEFLT, ICHAR, IDATA, IDEFLT,
     '    0, MN, LDATA, LDEFLT, RDATA, RDEFLT, -RMAX, RMAX, INFO,
     '    ERROR, *9999)
        IF( ISTABILISE.EQ.2 ) THEN
          ! TGSVD
          REG_PARAMETER(i) = DBLE(IDATA(1))
        ELSEIF( ISTABILISE.EQ.3 ) THEN
          ! Tikhonov
          REG_PARAMETER(i) = WORK(IDATA(1))
        ENDIF
      ENDDO
      CALL FREE_MEMORY(D_PTS_PTR, ERROR, *9999)

      CALL EXITS('PICARD')
      RETURN
 9999 CALL ERRORS('PICARD',ERROR)
      CALL EXITS('PICARD')
      RETURN 1
      END

C---------------------------------------------------------------------
