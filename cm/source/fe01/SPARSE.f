      SUBROUTINE SPARSE(I,J,N,nz,NZMAX,NZTOT,ISC,ISR,SPARSENESS,
     '  ERROR,*)

C#### Subroutine: SPARSE
C###  Description:
C###    SPARSE gives the position (nz) in the sparse array that
C###    corresponds to position (i,j) in the non-sparse array at that
C###    position. If position (i,j) cannot be found in the sparsity
C###    structure arrays the routine gives nz=0.


C#### Comment: SPARSITY STRUCTURES
C###  Description:
C###    <HTML>
C###    The sparsity structures used are governed by the SPARSENESS
C###    parameter associated with the array. If SPARSENESS=0 the
C###    matrix is not sparse and the the non-sparse matrix leading
C###    dimension NMAX is used to calculate matrix positions. If
C###    SPARSENESS=1 the matrix has compressed row sparsity (see
C###    below) and the sparsity structure arrays ISR and ISC are used
C###    for the calculation. If SPARSENESS=2 the matrix has row-column
C###    sparsity (see below) and the sparsity structure arrays ISR and
C###    ISC are used for the calculation. If SPARSENESS=3 the matrix has
C###    compressed-column sparsity (see below) and the sparsity
C###    structure arrays ISR and ISC are used for the calculation.  If
C###    SPARSENESS=4 the matrix has sorted row-column sparsity (see
C###    below) and the sparsity structure arrays ISR and
C###    ISC are used for the calculation. IF SPARSENESS=5 then
C###    the matrix has row-column sparsity but with the column
C###    information appended to the end of the row information in a
C###    single array. Only ISC is used for the positions and it must
C###    be at least twice NZA long. If SPARSENESS=6 then the matrix has
C###    row-column sparsity, but bisection method is not suitable.
C###    <P>
C###    Sparsity storage structures:
C###    </P> <P>
C###    COMPRESSED-ROW FORMAT:
C###    </P>
C###      The sparsity scheme is based on storing a NxN matrix, GX, as a
C###      one dimensional array of length NZT (where nz=sxNxN, s is the
C###      sparsity of the array) that stores only the non-zero elements
C###      of GX. Two additional arrays ISR and ISC store the positions
C###      of the non-zero elements. ISR is of length N+1 and ISC is of
C###      length NZT. ISR(i) stores the position in ISC of the start of
C###      row i. The N+1 position of ISR stores the size of ISC+1 ie.
C###      NZT+1. The number of non-zero elements in row i can be found
C###      from ISR(i+1)-ISR(i). ISC(nz) gives the column number for
C###      non-zero element nz. See also COMPRESSED-COLUMN format.
C###
C###      Example of sparsity storage scheme on a NxN matrix (N=6). Here
C###      the sparsity is 8/36 or 22%
C###      <PRE>
C###
C###      GX  1 2 3 4 5 6
C###         ____________        GX(nz)
C###       1| 0 A 0 B 0 0          A B C D E F G H
C###       2| 0 0 C 0 0 0
C###       3| 0 0 0 0 D E        ISR(i)
C###       4| F 0 0 0 0 0          1 3 4 6 7 8 9
C###       5| 0 0 G 0 0 0        ISC(i)
C###       6| 0 0 0 0 0 H          2 4 3 5 6 1 3 6
C###
C###      </PRE> </P> <P>
C###    ROW-COLUMN FORMAT:
C###    <P>
C###      The sparsity scheme is based on storing a NxN matrix, GX, as a
C###      one dimensional array of length NZT (where nz=sxNxN, s is the
C###      sparsity of the array) that stores only the non-zero elements
C###      of GX. Two additional arrays ISC and ISR store the positions
C###      of the non-zero elements. Both ISR and ISC are of length NZT.
C###      ISR(nz) stores the row number of the non-zero element nz and
C###      ISC(nz) stores the column number of the element.
C###
C###      The elements are assumed to be arranged in GX with all the
C###      elements in each row grouped together.
C###
C###      Example of sparsity storage scheme on a NxN matrix (N=6). Here
C###      the sparsity is 8/36 or 22%
C###      <PRE>
C###
C###      GX  1 2 3 4 5 6
C###         ____________        GX(nz)
C###       1| 0 A 0 B 0 0          H A C D G B E F
C###       2| 0 0 C 0 0 0
C###       3| 0 0 0 0 D E        ISR(i)
C###       4| F 0 0 0 0 0          6 1 2 3 5 1 3 4
C###       5| 0 0 G 0 0 0        ISC(i)
C###       6| 0 0 0 0 0 H          6 2 3 5 3 4 6 1
C###
C###      </PRE> </P> <P>
C###    COMPRESSED-COLUMN FORMAT:
C###    </P>
C###      The sparsity scheme is based on storing a NxN matrix, GX, as a
C###      one dimensional array of length NZT (where nz=sxNxN, s is the
C###      sparsity of the array) that stores only the non-zero elements
C###      of GX. Two additional arrays ISR and ISC store the positions
C###      of the non-zero elements. ISR is of length NZT and ISC is of
C###      length N+1. ISC(i) stores the position in ISR of the start of
C###      column i. The N+1 position of ISC stores the size of ISR+1 ie.
C###      NZT+1. The number of non-zero elements in column i can be
C###      found from ISC(i+1)-ISC(i). ISR(nz) gives the row number
C###      fornon-zero element nz. See also COMPRESSED-ROW format.
C###
C###      Example of sparsity storage scheme on a NxN matrix (N=6). Here
C###      the sparsity is 8/36 or 22%
C###      <PRE>
C###
C###      GX  1 2 3 4 5 6
C###         ____________        GX(nz)
C###       1| 0 A 0 B 0 0          F A C G B D E H
C###       2| 0 0 C 0 0 0
C###       3| 0 0 0 0 D E        ISR(i)
C###       4| F 0 0 0 0 0          4 1 2 5 1 3 3 6
C###       5| 0 0 G 0 0 0        ISC(i)
C###       6| 0 0 0 0 0 H          1 2 3 5 6 7 9
C###
C###      </PRE> </P> <P>
C###    SORTED ROW-COLUMN FORMAT:
C###    <P>
C###      This sparsity scheme is the same as ROW-COLUMN format except
C###      that the row numbers ISR(i) are assumed to be in
C###      non-decreasing order and the column numbers ISC(i) within each
C###      row are assumed to be in increasing order.  This enables
C###      faster searching of the nz for given row and column numbers.
C###
C###      Example of sparsity storage scheme on a NxN matrix (N=6). Here
C###      the sparsity is 8/36 or 22%
C###      <PRE>
C###
C###      GX  1 2 3 4 5 6
C###         ____________        GX(nz)
C###       1| 0 A 0 B 0 0          A B C D E F G H
C###       2| 0 0 C 0 0 0
C###       3| 0 0 0 0 D E        ISR(i)
C###       4| F 0 0 0 0 0          1 1 2 3 3 4 5 6
C###       5| 0 0 G 0 0 0        ISC(i)
C###       6| 0 0 0 0 0 H          2 4 3 5 6 1 3 6
C###
C###    </PRE> </P> <P>
C###    UMFPACK ROW-COLUMN FORMAT:
C###    <P>
C###      This sparsity scheme is the same as ROW-COLUMN format except
C###      that the column indices ISC(NZA+i) are stored after the
C###      row indices ISC(i)
C###
C###      Example of sparsity storage scheme on a NxN matrix (N=6). Here
C###      the sparsity is 8/36 or 22%
C###      <PRE>
C###
C###      GX  1 2 3 4 5 6
C###         ____________        GX(nz)
C###       1| 0 A 0 B 0 0          H A C D G B E F
C###       2| 0 0 C 0 0 0
C###       3| 0 0 0 0 D E        ISC(i)
C###       4| F 0 0 0 0 0          6 1 2 3 5 1 3 4 6 2 3 5 3 4 6 1
C###       5| 0 0 G 0 0 0          (rows         )(columns        )
C###       6| 0 0 0 0 0 H
C###
C###    </PRE> </P> </HTML>

      IMPLICIT NONE
!     Parameter List
      INTEGER I,J,N,nz,NZMAX,NZTOT,ISC(*),ISR(*),SPARSENESS
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LOWLIMIT,MIDPOINT,UPLIMIT
      LOGICAL FOUNDCOLUMN,FOUNDROW

C      CALL ENTERS('SPARSE',*9999)

      IF(SPARSENESS.EQ.0) THEN !No sparsity
        nz=I+(J-1)*N
        IF(nz.GT.NZMAX.OR.nz.LT.1) THEN
          nz=0
          ERROR='>>Array coordinates outside range[0]'
          GOTO 9999
        ENDIF
      ELSE IF(SPARSENESS.EQ.1) THEN !Compressed Row format
        IF(I.LT.1.OR.I.GT.N) THEN
          ERROR='>>Array coordinates outside range[1]'
          GOTO 9999
        ENDIF
C***    Search for the column number in the sparsity list using the
C***    bisection (binary search) algorithm
        LOWLIMIT=ISR(I)
        UPLIMIT=ISR(I+1)
        MIDPOINT=(UPLIMIT+LOWLIMIT)/2
        DO WHILE(MIDPOINT.GT.LOWLIMIT)
          IF(ISC(MIDPOINT).GT.J) THEN
            UPLIMIT=MIDPOINT
          ELSE
            LOWLIMIT=MIDPOINT
          ENDIF
          MIDPOINT=(UPLIMIT+LOWLIMIT)/2
        ENDDO
        IF(ISC(LOWLIMIT).EQ.J) THEN
          nz=LOWLIMIT
        ELSE
          nz=0
        ENDIF
      ELSE IF(SPARSENESS.EQ.2) THEN !Row-Column format
        FOUNDROW=.FALSE.
        nz=1
        DO WHILE(.NOT.FOUNDROW.AND.nz.LE.NZTOT)
          IF(ISR(nz).EQ.I) THEN
            FOUNDROW=.TRUE.
            FOUNDCOLUMN=.FALSE.
            DO WHILE(.NOT.FOUNDCOLUMN.AND.nz.LE.NZTOT)
              IF(ISC(nz).EQ.J.AND.ISR(nz).EQ.I) THEN
                FOUNDCOLUMN=.TRUE.
              ELSE IF(ISR(nz).NE.I) THEN
                nz=NZTOT+1
              ELSE
                nz=nz+1
              ENDIF
            ENDDO
          ELSE
            nz=nz+1
          ENDIF
        ENDDO
        IF(.NOT.(FOUNDROW.AND.FOUNDCOLUMN)) THEN
          nz=0
        ENDIF
      ELSE IF(SPARSENESS.EQ.3) THEN !Compressed Column format
        IF(J.LT.1.OR.J.GT.N) THEN
          ERROR='>>Array coordinates outside range[3]'
          GOTO 9999
        ENDIF
C***    Search for the row number in the sparsity list using the
C***    bisection (binary search) algorithm
        LOWLIMIT=ISC(J)
        UPLIMIT=ISC(J+1)
        MIDPOINT=(UPLIMIT+LOWLIMIT)/2
        DO WHILE(MIDPOINT.GT.LOWLIMIT)
          IF(ISR(MIDPOINT).GT.I) THEN
            UPLIMIT=MIDPOINT
          ELSE
            LOWLIMIT=MIDPOINT
          ENDIF
          MIDPOINT=(UPLIMIT+LOWLIMIT)/2
        ENDDO
        IF(ISR(LOWLIMIT).EQ.I) THEN
          nz=LOWLIMIT
        ELSE
          nz=0
        ENDIF
      ELSE IF(SPARSENESS.EQ.4) THEN !Sorted Row-Column format
C***    Search for the coordinates in the sparsity list using the
C***    bisection (binary search) algorithm
        LOWLIMIT=1
        UPLIMIT=NZTOT
        DO WHILE(UPLIMIT.GT.LOWLIMIT)
          MIDPOINT=(UPLIMIT+LOWLIMIT)/2
          IF(ISR(MIDPOINT).EQ.I) THEN
            IF(ISC(MIDPOINT).LT.J) THEN
              LOWLIMIT=MIDPOINT+1
            ELSE
              UPLIMIT=MIDPOINT
            ENDIF
          ELSE IF(ISR(MIDPOINT).LT.I) THEN
            LOWLIMIT=MIDPOINT+1
          ELSE
            UPLIMIT=MIDPOINT
          ENDIF
        ENDDO
        IF(ISR(LOWLIMIT).EQ.I.AND.ISC(LOWLIMIT).EQ.J) THEN
          nz=LOWLIMIT
        ELSE
          nz=0
        ENDIF
      ELSE IF(SPARSENESS.EQ.5) THEN !Sorted UMFPACK Row-Column format
C***    Search for the coordinates in the sparsity list using the
C***    bisection (binary search) algorithm
        LOWLIMIT=1
        UPLIMIT=NZTOT+1
        MIDPOINT=(UPLIMIT+LOWLIMIT)/2
        DO WHILE(MIDPOINT.GT.LOWLIMIT)
          IF(ISC(MIDPOINT).EQ.I) THEN
            IF(ISC(NZTOT+MIDPOINT).GT.J) THEN
              UPLIMIT=MIDPOINT
            ELSE
              LOWLIMIT=MIDPOINT
            ENDIF
          ELSE IF(ISC(MIDPOINT).GT.I) THEN
            UPLIMIT=MIDPOINT
          ELSE
            LOWLIMIT=MIDPOINT
          ENDIF
          MIDPOINT=(UPLIMIT+LOWLIMIT)/2
        ENDDO
        IF(ISC(LOWLIMIT).EQ.I.AND.ISC(NZTOT+LOWLIMIT).EQ.J) THEN
          nz=LOWLIMIT
        ELSE
          nz=0
        ENDIF
C KAT 5Nov98: way too slow
C        FOUNDROW=.FALSE.
C        nz=1
C        DO WHILE(.NOT.FOUNDROW.AND.nz.LE.NZTOT)
C          IF(ISC(nz).EQ.I) THEN
C            FOUNDROW=.TRUE.
C            FOUNDCOLUMN=.FALSE.
C            DO WHILE(.NOT.FOUNDCOLUMN.AND.nz.LE.NZTOT)
C              IF(ISC(nz+NZTOT).EQ.J.AND.ISC(nz).EQ.I) THEN
C                FOUNDCOLUMN=.TRUE.
C              ELSE IF(ISC(nz).NE.I) THEN
C                nz=NZTOT+1
C              ELSE
C                nz=nz+1
C              ENDIF
C            ENDDO
C          ELSE
C            nz=nz+1
C          ENDIF
C        ENDDO
C        IF(.NOT.(FOUNDROW.AND.FOUNDCOLUMN)) THEN
C          nz=0
C        ENDIF
      ELSE IF(SPARSENESS.EQ.6)THEN
C MHT 02-02-99 Bisection method doesn't work for lung sparsity.
C The following is quick for lung sparsity, but not for general
C   row-column sparsity.
        IF(I.LT.1.OR.I.GT.N) THEN
          ERROR='>>Array coordinates outside range[6]'
          GOTO 9999
        ENDIF
        LOWLIMIT=ISR(I)
        UPLIMIT=ISR(I+1)
        nz=0
        DO WHILE(LOWLIMIT.LE.UPLIMIT)
          IF(ISC(LOWLIMIT).EQ.J)THEN
            nz=LOWLIMIT
            LOWLIMIT=UPLIMIT+1
          ELSE
            LOWLIMIT=LOWLIMIT+1
          ENDIF
        ENDDO
      ELSE
        ERROR='>>Unknown sparsity format'
        GOTO 9999
      ENDIF

      RETURN
C 9998 CALL EXITS('SPARSE')
C      RETURN
 9999 CALL ERRORS('SPARSE',ERROR)
C      CALL EXITS('SPARSE')
      RETURN 1
      END


