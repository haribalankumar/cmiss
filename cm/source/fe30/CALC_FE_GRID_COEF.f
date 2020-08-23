      SUBROUTINE CALC_FE_GRID_COEF(nb,NITB,nq,NQGP,NQGP_PIVOT,
     &  nr,nx_ext,NXQ,COEFFSEXT,CQ,GM,NQGW,PG,WG,XQ,BIDOMAIN,
     &  FIXQ,SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: CALC_FE_GRID_COEF
C###  Description:
C###    CALC_FE_GRID_COEF is used for the Grid-based Finite Element
C###    scheme. It finds the coeffs of the nonzeros in stiffness
C###    matrix (NQGW) and adds coefficients to global mass matrix (GM)
C###    for a single nq. The number of nonzeros depends on the number
C###    of local xi coords.
C***  Created by Scott Marsden 18 October 2000

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc' !NJL_FIBR

!     Parameter list
      INTEGER nb,NITB,nq,NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),nr,
     &  nx_ext,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 COEFFSEXT(NQGM),CQ(NMM,NQM),GM(NZ_GM_M),NQGW(NQGM),
     &  PG(NSM,NUM,NGM,NBM),WG(NGM,NBM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,FIXQ(NYQM,NIYFIXM,NXM),SOLVEEIGHTPROBLEM
!     Local variables
      INTEGER GB_ELEM(8),i,m,NGTB,ni1,ni2,ni2to,ni3,ni3to,NNTB
      LOGICAL INT_ELEM,EXT_ELEM

      CALL ENTERS('CALC_FE_GRID_COEF',*9999)

      NNTB = NNT(nb)
      NGTB = NGT(nb)

C initialise NQGP, NQGW a row of GM
      NQGP(0,nq)=0
      DO i=1,NQGM ! initialise variables
        NQGP(i,nq)=0
        COEFFSEXT(i)=0.0d0
        NQGW(i)=0.0d0
        GM(i+NQGM*(nq-1))=0.0d0
      ENDDO !i

C check that there are no nodal derivatives stored at nodes as only
C     linear basis functions can be interpolated at the moment
      CALL ASSERT(NKT(0,nb).EQ.1.AND.NNTB.EQ.2**NITB,
     &  '>>Must use a linear basis to interpolate'
     &  //' dependent variable for grid-based FE',ERROR,*9999)

C check that the number of Xi directions for the element is equal to
C     the number of global coordinates
      CALL ASSERT(NITB.EQ.NJ_LOC(NJL_GEOM,0,nr),
     &  '>>#Xi directions for element must equal #global'
     &  //' coordinates (for grid-based FE method)',ERROR,*9999)

C evaluate Xi loop variables
      ni2to = -1
      ni3to = -1
      IF(NITB.GE.2) THEN
        ni2to = 1
        IF(NITB.GE.3) THEN
          ni3to = 1
        ENDIF
      ENDIF

      m = NNTB+1 !this variable loops over element grid points
      DO ni3=-1,ni3to,2
        DO ni2=-1,ni2to,2
          DO ni1=-1,1,2
            m=m-1 ! m becomes NNTB for first time through loop

            INT_ELEM = .FALSE.
            EXT_ELEM = .FALSE.

C If problem is bidomain and the grid point is not an essential or
C           analytic bc, check that element can be created in the
C           extracellular domain.
C NOTE: absolute values are used on grid references because for
C           intracellular discontinuities, a negative sign is put
C           before the grid number (done in UPGRID_CONN).
            IF (BIDOMAIN) THEN
C!!! KAT 21Feb01: Need NIM>=3 even for 1d
              CALL ASSERT(NIM.GE.3,'Increase NIM to 3',ERROR,*9999)
              IF(.NOT.(FIXQ(nq,1,nx_ext).OR.FIXQ(nq,3,nx_ext))) THEN
                IF ( NXQ(ni1,0,nq).GT.0                   !gp in 1 dirn
     &            .AND.(NXQ(ni2*2,0,nq).GT.0.OR.NITB.LT.2) !gp in 2 dirn
     &            .AND.(NXQ(ni3*3,0,nq).GT.0.OR.NITB.LT.3) !gp in 3 dirn
     &            )THEN
                  IF ( ((NXQ(ni1,0,ABS(NXQ(ni2*2,1,nq))).GT.0
     &              .AND.NXQ(ni2*2,0,ABS(NXQ(ni1,1,nq))).GT.0)
     &              .OR.NITB.LT.2)   !conn in 2D (or 1D problem)
     &              .AND.((NXQ(ni1,0,ABS(NXQ(ni3*3,1,nq))).GT.0
     &              .AND.NXQ(ni2*2,0,ABS(NXQ(ni3*3,1,nq))).GT.0
     &              .AND.NXQ(ni3*3,0,ABS(NXQ(ni1,1,nq))).GT.0
     &              .AND.NXQ(ni3*3,0,ABS(NXQ(ni2*2,1,nq))).GT.0)
     &              .OR.NITB.LT.3))THEN
                    IF ((
     &                NXQ(ni1,0,ABS(NXQ(ni2*2,1,ABS(NXQ(ni3*3,1,
     &                nq))))).GT.0.AND.
     &                NXQ(ni2*2,0,ABS(NXQ(ni3*3,1,ABS(NXQ(ni1,1,
     &                nq))))).GT.0.AND.
     &                NXQ(ni3*3,0,ABS(NXQ(ni1,1,ABS(NXQ(ni2*2,1,
     &                nq))))).GT.0).OR.NITB.LT.3)THEN !can make extracellular elem
                      EXT_ELEM = .TRUE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF !bidomain
            ENDIF !not essential or analytic bc

C!!! KAT 22Feb00: Is NXQ(ni,1,nq) set if NXQ(ni,0,nq)=0 ?

C check that element can be created in the intracellular domain
            IF(.NOT.SOLVEEIGHTPROBLEM) THEN
C!!! DPN 10 September 2001 - KAT 21Feb01: Need NIM>=3 even for 1d
              CALL ASSERT(NIM.GE.3,'Increase NIM to 3',ERROR,*9999)
              IF ( NXQ(ni1,1,nq).GT.0                   !gp in 1 dirn
     &          .AND.(NXQ(ni2*2,1,nq).GT.0.OR.NITB.LT.2) !gp in 2 dirn
     &          .AND.(NXQ(ni3*3,1,nq).GT.0.OR.NITB.LT.3) !gp in 3 dirn
     &          )THEN
                IF ( ((NXQ(ni1,1,NXQ(ni2*2,1,nq)).GT.0
     &            .AND.NXQ(ni2*2,1,NXQ(ni1,1,nq)).GT.0)
     &            .OR.NITB.LT.2)   !conn in 2D or 1D problem
     &            .AND.((NXQ(ni1,1,NXQ(ni3*3,1,nq)).GT.0
     &            .AND.NXQ(ni2*2,1,NXQ(ni3*3,1,nq)).GT.0
     &            .AND.NXQ(ni3*3,1,NXQ(ni1,1,nq)).GT.0
     &            .AND.NXQ(ni3*3,1,NXQ(ni2*2,1,nq)).GT.0)
     &            .OR.NITB.LT.3))THEN
                  IF ((NXQ(ni1,1,NXQ(ni2*2,1,NXQ(ni3*3,1,nq)))
     &              .GT.0.AND.
     &              NXQ(ni2*2,1,NXQ(ni3*3,1,NXQ(ni1,1,nq)))
     &              .GT.0.AND.
     &              NXQ(ni3*3,1,NXQ(ni1,1,NXQ(ni2*2,1,nq)))
     &              .GT.0).OR.NITB.LT.3)THEN !can make intracellular elem

                    INT_ELEM = .TRUE.
                  ENDIF
                ENDIF
              ENDIF
            ENDIF

            IF(.NOT.INT_ELEM.AND.DOP)THEN
              WRITE(OP_STRING,
     &          '('' cannot make intracellular element: '',4I5)')
     &          nq,ni1,ni2,ni3
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

            IF(INT_ELEM.OR.EXT_ELEM) THEN
C create grid based element with grid points as nodes
C global grid numbers for the grid-based element are stored in GB_ELEM
              GB_ELEM(m) = nq
              GB_ELEM(m+ni1) = ABS(NXQ(ni1,1,nq))
              IF (NITB.GT.1)  THEN ! 2D
                GB_ELEM(m+ni2*2) = ABS(NXQ(ni2*2,1,nq))
                GB_ELEM(m+ni1+ni2*2) =
     &            ABS(NXQ(ni1,1,ABS(NXQ(ni2*2,1,nq))))
                IF (NITB.GT.2)  THEN ! 3D
                  GB_ELEM(m+ni3*4) = ABS(NXQ(ni3*3,1,nq))
                  GB_ELEM(m+ni1+ni3*4) =
     &              ABS(NXQ(ni1,1,ABS(NXQ(ni3*3,1,nq))))
                  GB_ELEM(m+ni2*2+ni3*4) =
     &              ABS(NXQ(ni2*2,1,ABS(NXQ(ni3*3,1,nq))))
                  GB_ELEM(m+ni1+ni2*2+ni3*4) =
     &              ABS(NXQ(ni1,1,ABS(NXQ(ni2*2,1,ABS(NXQ(ni3*3,1,
     &              nq))))))
                ENDIF
              ENDIF

              CALL CALC_FE_COEF(GB_ELEM,m,NGTB,NITB,NNTB,nq,
     &          NQGP(0,nq),nr,nx_ext,COEFFSEXT,CQ,GM,NQGW,
     &          PG(1,1,1,nb),WG(1,nb),XQ,BIDOMAIN,EXT_ELEM,FIXQ,
     &          INT_ELEM,SOLVEEIGHTPROBLEM,ERROR,*9999)

C KAT 21Feb01: Should move outside loop after check no surrounding
C             elements set COEFFSEXT diagonal to 1 for analytic or
C             essential bc
              IF (BIDOMAIN) THEN
                IF (FIXQ(nq,1,nx_ext).OR.FIXQ(nq,3,nx_ext))THEN
                  DO i=1,NQGP(0,nq)
                    IF (NQGP(i,nq).EQ.nq) COEFFSEXT(i) = 1.0d0
                  ENDDO !i
                ENDIF !analytic or essential bc
              ENDIF !bidomain

            ENDIF !can make internal or external element
          ENDDO !ni1
        ENDDO !ni2
      ENDDO !ni3

C essential bc for solve8 problem or no surrounding elements (intracellular)

C!!! KAT 21Feb01: Need to set diagonal of GM if there are no surrounding
C!!! intra elements even if there are surrounding extra elements.

      IF (NQGP(0,nq).EQ.0) THEN
        NQGP(0,nq) = 1
        NQGP(1,nq) = nq
        IF(SOLVEEIGHTPROBLEM) THEN !essential bc
C fix grid point to rhs value
          COEFFSEXT(1)= 1.0d0
        ELSE !no surrounding elements
          NQGW(1)= 0.0d0
          COEFFSEXT(1)= 0.0d0 !??? IS THIS RIGHT ???
          GM(1+NQGM*(nq-1)) = 1.0d0
        ENDIF
      ENDIF

C sort NQGP
      CALL ISORTP(NQGP(0,nq),NQGP(1,nq),NQGP_PIVOT(1,nq))

      CALL EXITS('CALC_FE_GRID_COEF')
      RETURN
 9999 CALL ERRORS('CALC_FE_GRID_COEF',ERROR)
      CALL EXITS('CALC_FE_GRID_COEF')
      RETURN 1
      END
