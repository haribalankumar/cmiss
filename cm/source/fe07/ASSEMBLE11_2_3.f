      SUBROUTINE ASSEMBLE11_2_3(ISC_GKK,ISR_GKK,NEELEM,NENQ,NPNODE,
     '  NQGP,NQGP_PIVOT,NQLIST,NQLIST2,NQLIST3,NQNP,NQNP_LOCAL,
     '  NQNP_LOCAL_PIV,NQS,NQXI,nr,nr2,NWQ,nx,nx2,NXQ,CQ,
     '  D_MATRIX,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GK,GKK,GQ,
     '  GQ_D_MATRIX,GKGK2,GQGQ2,GUQ,NQGW,PROPQ,FIXQ,
     '  SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: ASSEMBLE11_2_3
C###  Description:
C###    <HTML> <PRE>
C###    ASSEMBLE11_2_3 creates matrices for coupled BEM-FD
C###    problems where both are assembled into 1 matrix.
C###
C###    NOTE: GK,GQ must be fully populated
C###    </PRE> </HTML>
C***  Created by Martin Buist, May 2000

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),
     '  NPNODE(0:NP_R_M,0:NRM),NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),
     '  NQLIST(0:NQM),NQLIST2(0:NQM),NQLIST3(0:NQM),NQNP(NPM),
     '  NQNP_LOCAL(NPM),NQNP_LOCAL_PIV(NPM),NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  nr,nr2,NWQ(8,0:NQM),nx,nx2,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 CQ(NMM,NQM),D_MATRIX(NPM,*),DNUDXQ(3,3,NQM),
     '  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GK(NZ_GK_M),
     '  GKK(NZ_GKK_M,NXM),GQ(NZ_GQ_M),GQ_D_MATRIX(NPM,*),
     '  GUQ(3,3,NQM),NQGW(NQGM,NQM),GKGK2(NPM,NPM),GQGQ2(NPM,NPM),
     '  PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM),SOLVEEIGHTPROBLEM
!     Local variables
      INTEGER BCINNBND,DUMMY_LIST(0:1),IDUMMY,ncol,NITB,np,nq,nq2,nqq,
     '  nrow,ntemp,nzero,nzz,PLACEEXT
      REAL*8 ANS,COEFFSEXT(NQGM),PVTTMP,RDUMMY(NQGM)
      REAL TIME_START1(1),TIME_START2(1),TIME_STOP1(1),TIME_STOP2(1),
     '  ELAPSED_TIME
      CHARACTER LOCFILENAME*(MXCH)
      LOGICAL ERROR_FLAG,INLIST

      CALL ENTERS('ASSEMBLE11_2_3',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      NITB=NQXI(0,NQS(NEELEM(1,nr)))

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Boundary Element GK,GQ '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Region '',I6,'' Class '',I6)') nr2,nx2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Rows '',I6,'' Columns '',I6)')
     '    NYT(1,1,nx2),NYT(2,1,nx2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nzz=1,NYT(1,1,nx2)*NYT(2,1,nx2)
          WRITE(OP_STRING,'(I8,'' GK '',F12.6,'' GQ '',F12.6)')
     '      nzz,GK(nzz),GQ(nzz)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nzz
      ENDIF

!     Calc system sizes
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C$OMP PARALLEL DO
C$OMP&PRIVATE(nq),
C$OMP&SHARED(NWQ)
      DO nq=1,NQT
        NWQ(3,nq)=0
      ENDDO !nq
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
C$OMP&PRIVATE(nq,nqq,nq2,ntemp),
C$OMP&SHARED(NQGP,NQLIST,NWQ)
      DO nqq=1,NQLIST(0)
        nq=NQLIST(nqq)
        DO ntemp=1,NQGP(0,nq)
          nq2=NQGP(ntemp,nq)
          NWQ(3,nq2)=2
        ENDDO
      ENDDO !nq
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
C$OMP&PRIVATE(nq,nqq),
C$OMP&SHARED(NQLIST,NWQ)
      DO nqq=1,NQLIST(0)
        nq=NQLIST(nqq)
        NWQ(3,nq)=0
      ENDDO !nq
C$OMP END PARALLEL DO

      BCINNBND=0
      NQLIST2(0)=0
      DO nq=1,NQT
        IF(NWQ(3,nq).EQ.2) THEN
          NQLIST2(0)=NQLIST2(0)+1
          NQLIST2(NQLIST2(0))=nq
          BCINNBND=BCINNBND+1
        ENDIF
      ENDDO !nq
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'(/'' Time to calc system sizes : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      DO nq=1,NQLIST(0)
        NQLIST2(0)=NQLIST2(0)+1
        NQLIST2(NQLIST2(0))=NQLIST(nq)
      ENDDO
      CALL ISORT(NQLIST2(0),NQLIST2(1))

C$OMP PARALLEL DO
C$OMP&PRIVATE(nq),
C$OMP&SHARED(NQLIST3)
      DO nq=0,NQM
        NQLIST3(nq)=0
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
C$OMP&PRIVATE(nq),
C$OMP&SHARED(NQLIST2,NQLIST3)
      DO nq=1,NQLIST2(0)
        NQLIST3(NQLIST2(nq))=nq
      ENDDO
C$OMP END PARALLEL DO

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Inner Bound conn nq = '',I8)') BCINNBND
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Region = '',I8)') nr2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nqq=1,NPNODE(0,nr2)
          np=NPNODE(nqq,nr2)
          WRITE(OP_STRING,'(2I8)') np,NQNP(np)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

!     Create the D matrix
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      nrow=NQLIST(0)
      ncol=NQLIST2(0)
C$OMP PARALLEL DO
C$OMP&PRIVATE(nq,nq2),
C$OMP&SHARED(ncol,nrow,D_MATRIX)
      DO nq=1,nrow
        DO nq2=1,ncol
          D_MATRIX(nq,nq2)=0.0d0
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

      ERROR_FLAG=.FALSE.
      IF(SOLVEEIGHTPROBLEM) THEN
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(nq,nqq,nzero),
C$OMP&  SHARED(CQ,D_MATRIX,DNUDXQ,DXDXIQ,DXDXIQ2,NENQ,NQLIST,NQLIST3,
C$OMP&  NQS,NQXI,NXQ,NQGP,NQGP_PIVOT,NQGW,ERROR_FLAG)
        DO nqq=1,NQLIST(0)
          IF(.NOT.ERROR_FLAG) THEN
            nq=NQLIST(nqq)
            CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,NXQ,NQGW(1,nq),
     '        CQ(1,nq),DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,*100)
            DO nzero=1,NQGP(0,nq)
              D_MATRIX(nqq,NQLIST3(NQGP(nzero,nq)))=
     '          NQGW(NQGP_PIVOT(nzero,nq),nq)
            ENDDO
            GO TO 102
C           This statement is designed to be skipped if no error
C           occurs. However if a error occurs within a subroutine
C           the alternate return points to line 100 to set the flag
 100        CONTINUE
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*101)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '        //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
 102        CONTINUE
          ENDIF !.NOT.ERROR_FLAG
        ENDDO
C$OMP   END PARALLEL DO
      ELSE
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(COEFFSEXT,nq,nqq,nzero),
C$OMP&  SHARED(CQ,D_MATRIX,DNUDXQ,DXDXIQ,DXDXIQ2,NENQ,NQLIST,NQLIST3,
C$OMP&  NQS,NQXI,NXQ,NQGP,NQGP_PIVOT,ERROR_FLAG)
        DO nqq=1,NQLIST(0)
          IF(.NOT.ERROR_FLAG) THEN
            nq=NQLIST(nqq)
            CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,NXQ,COEFFSEXT,
     '        CQ(6,nq),DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,*200)
            DO nzero=1,NQGP(0,nq)
              D_MATRIX(nqq,NQLIST3(NQGP(nzero,nq)))=
     '          COEFFSEXT(NQGP_PIVOT(nzero,nq))
            ENDDO
            GO TO 202
C           This statement is designed to be skipped if no error
C           occurs. However if a error occurs within a subroutine
C           the alternate return points to line 200 to set the flag
 200        CONTINUE
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*201)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '        //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201        CONTINUE
 202        CONTINUE
          ENDIF !.NOT.ERROR_FLAG
        ENDDO
C$OMP   END PARALLEL DO
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time create D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Make a local copy of GK,GQ
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C$OMP PARALLEL DO
C$OMP&PRIVATE(nq,nzero),
C$OMP&SHARED(GK,GKGK2,GQ,GQGQ2,nr2,nx2)
      DO nq=1,NYT(1,1,nx2)
        DO nzero=1,NYT(2,1,nx2)
          GQGQ2(nq,nzero)=
     '      GQ(((nzero-1)*NYT(1,1,nx2))+nq)
          GKGK2(nq,nzero)=
     '      GK(((nzero-1)*NYT(1,1,nx2))+nq)
        ENDDO
      ENDDO
C$OMP END PARALLEL DO
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time duplicate GK,GQ : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START2)
!     Need to adjust GK,GQ to reflect grid numbering
!     swap 2 rows at once, element by element so only
!     one temporary memory space is needed. Both the rows
!     and the columns need to be pivoted.

C$OMP PARALLEL DO
C$OMP&PRIVATE(np,nqq),
C$OMP&SHARED(NQNP,NQNP_LOCAL,NPNODE,nr2)
      DO nqq=1,NPNODE(0,nr2)
        np=NPNODE(nqq,nr2)
        NQNP_LOCAL(nqq)=NQNP(np)
      ENDDO
C$OMP END PARALLEL DO
      CALL ISORTP(NPNODE(0,nr2),NQNP_LOCAL,NQNP_LOCAL_PIV)

      IF(DOP) THEN
        DO nqq=1,NPNODE(0,nr2)
          WRITE(OP_STRING,'(3I8)') nqq,NQNP_LOCAL(nqq),
     '      NQNP_LOCAL_PIV(nqq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

!     Pivot rows
      DO nqq=1,NPNODE(0,nr2)
        IF(NQNP_LOCAL_PIV(nqq).NE.nqq) THEN !must pivot
          DO nq2=nqq+1,NPNODE(0,nr2)
            IF(NQNP_LOCAL_PIV(nq2).EQ.nqq) THEN !found pivot
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Row pivot '',I8,'' with '',I8)')
     '            nqq,NQNP_LOCAL_PIV(nqq)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO nzero=1,NPNODE(0,nr2)
                PVTTMP=GKGK2(nqq,nzero)
                GKGK2(nqq,nzero)=
     '            GKGK2(NQNP_LOCAL_PIV(nqq),nzero)
                GKGK2(NQNP_LOCAL_PIV(nqq),nzero)=PVTTMP

                PVTTMP=GQGQ2(nqq,nzero)
                GQGQ2(nqq,nzero)=
     '            GQGQ2(NQNP_LOCAL_PIV(nqq),nzero)
                GQGQ2(NQNP_LOCAL_PIV(nqq),nzero)=PVTTMP
              ENDDO
              NQNP_LOCAL_PIV(nq2)=NQNP_LOCAL_PIV(nqq)
              NQNP_LOCAL_PIV(nqq)=nqq
              IF(DOP) THEN
                DO nzero=1,NPNODE(0,nr2)
                  WRITE(OP_STRING,'(2I8)') nzero,NQNP_LOCAL_PIV(nzero)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO

!     Pivot columns
C$OMP PARALLEL DO
C$OMP&PRIVATE(np,nqq),
C$OMP&SHARED(NQNP,NQNP_LOCAL,NPNODE,nr2)
      DO nqq=1,NPNODE(0,nr2)
        np=NPNODE(nqq,nr2)
        NQNP_LOCAL(nqq)=NQNP(np)
      ENDDO
C$OMP END PARALLEL DO
      CALL ISORTP(NPNODE(0,nr2),NQNP_LOCAL,NQNP_LOCAL_PIV)

      DO nqq=1,NPNODE(0,nr2)
        IF(NQNP_LOCAL_PIV(nqq).NE.nqq) THEN !must pivot
          DO nq2=nqq+1,NPNODE(0,nr2)
            IF(NQNP_LOCAL_PIV(nq2).EQ.nqq) THEN !found pivot
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Col pivot '',I8,'' with '',I8)')
     '            nqq,NQNP_LOCAL_PIV(nqq)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO nzero=1,NPNODE(0,nr2)
                PVTTMP=GKGK2(nzero,nqq)
                GKGK2(nzero,nqq)=
     '            GKGK2(nzero,NQNP_LOCAL_PIV(nqq))
                GKGK2(nzero,NQNP_LOCAL_PIV(nqq))=PVTTMP

                PVTTMP=GQGQ2(nzero,nqq)
                GQGQ2(nzero,nqq)=
     '            GQGQ2(nzero,NQNP_LOCAL_PIV(nqq))
                GQGQ2(nzero,NQNP_LOCAL_PIV(nqq))=PVTTMP
              ENDDO
              NQNP_LOCAL_PIV(nq2)=NQNP_LOCAL_PIV(nqq)
              NQNP_LOCAL_PIV(nqq)=nqq
              IF(DOP) THEN
                DO nzero=1,NPNODE(0,nr2)
                  WRITE(OP_STRING,'(2I8)') nzero,NQNP_LOCAL_PIV(nzero)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to pivot : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Create GQ*D
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C$OMP PARALLEL DO
C$OMP&PRIVATE(ANS,nq,nq2,nzero),
C$OMP&SHARED(D_MATRIX,GQ_D_MATRIX,GQGQ2,NQLIST,NQLIST2)
      DO nq=1,NQLIST(0)
        DO nq2=1,NQLIST2(0)
          ANS=0
          DO nzero=1,NQLIST(0)
            ANS=ANS+GQGQ2(nq,nzero)*D_MATRIX(nzero,nq2)
          ENDDO
          GQ_D_MATRIX(nq,nq2)=ANS
        ENDDO
      ENDDO
C$OMP END PARALLEL DO
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to create GQ*D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Create GK-(GQ*D)
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C$OMP PARALLEL DO
C$OMP&PRIVATE(nq,nq2),
C$OMP&SHARED(D_MATRIX,GQ_D_MATRIX,NQLIST,NQLIST2)
      DO nq=1,NQLIST(0)
        DO nq2=1,NQLIST2(0)
          D_MATRIX(nq,nq2)=GQ_D_MATRIX(nq,nq2)
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
C$OMP&PRIVATE(np,nq,nq2,nqq),
C$OMP&SHARED(D_MATRIX,GKGK2,NQLIST,NQLIST3)
      DO nq=1,NQLIST(0)
        DO nq2=1,NQLIST(0)
          np=NQLIST(nq2)
          nqq=NQLIST3(np)
          D_MATRIX(nq,nqq)=D_MATRIX(nq,nqq)+
     '      GKGK2(nq,nq2)
        ENDDO
      ENDDO
C$OMP END PARALLEL DO
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to create GK-GQ*D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Initialising
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      IF(SPARSEGKK(nx).EQ.0) THEN
        nzz=NQT*NQT
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(GKK,nx,nzz)
        DO nq=1,nzz
          GKK(nq,nx)=0.0d0
        ENDDO !nq
C$OMP   END PARALLEL DO
        NOT(1,1,nr,nx)=NQT
        NOT(2,1,nr,nx)=NQT
        NZZT(1,nr,nx)=nzz
        IF(NZ_GKK_M.LT.nzz) THEN
          WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
          GOTO 9999
        ENDIF
      ELSE IF(SPARSEGKK(nx).EQ.1 .OR. SPARSEGKK(nx).EQ.5) THEN !row/col #2
        nzz=0
        DO nq=NQR(1,nr),NQR(2,nr)
          DO nzero=1,NQGP(0,nq)
            nzz=nzz+1
          ENDDO !nzero
        ENDDO !nq
        NOT(1,1,nr,nx)=NQT
        NOT(2,1,nr,nx)=NQT
        NZZT(1,nr,nx)=nzz
        IF(NZ_GKK_M.LT.nzz) THEN
          WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
          GOTO 9999
        ELSE IF(NISC_GKKM.LT.2*nzz) THEN
          WRITE(ERROR,'(''>>Increase ISC_GKK_M to >= '',I12)') 2*nzz
          GOTO 9999
        ELSE IF(NISR_GKKM.LT.NQT+1) THEN
          WRITE(ERROR,'(''>>Increase ISR_GKK_M to >= '',I12)') NQT+1
          GOTO 9999
        ENDIF
      ENDIF
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to init matrix : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Assemble grid coeffs
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      IF(SOLVEEIGHTPROBLEM) THEN
        ntemp=NX_LIST(0)
        NX_LIST(0)=1
        IF(SPARSEGKK(nx).EQ.0) THEN
          DO nq=NQR(1,nr),NQR(2,nr)
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx)=
     '            COEFFSEXT(NQGP_PIVOT(nzero,nq))
              ENDDO !nz
            ELSE
              DO nzero=1,NQLIST(0)
                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                nq2=NQLIST2(nzero)
                GKK(((nq2-1)*NQT)+nq,nx)=D_MATRIX(nqq,nzero)
              ENDDO !nz
            ENDIF
          ENDDO !nq
        ELSE IF(SPARSEGKK(nx).EQ.1 .OR. SPARSEGKK(nx).EQ.5) THEN
          PLACEEXT=0
          DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
C                IF(DABS(COEFFSEXT(NQGP_PIVOT(nzero,nq))).GT.ZERO_TOL)
C     '            THEN
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx)=nq
                  ISR_GKK(PLACEEXT,nx)=NQGP(nzero,nq)
                  GKK(PLACEEXT,nx)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
C                ENDIF
              ENDDO !nzero
            ELSE
              DO nzero=1,NQLIST(0)
                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
C                IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                  nq2=NQLIST2(nzero)
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx)=nq
                  ISR_GKK(PLACEEXT,nx)=nq2
                  GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
C                ENDIF
              ENDDO !nzero
            ENDIF
          ENDDO !nq

          NZZT(1,nr,nx)=PLACEEXT
C$OMP     PARALLEL DO
C$OMP&    PRIVATE(nq),
C$OMP&    SHARED(ISC_GKK,ISR_GKK,nx,PLACEEXT)
          DO nq=1,PLACEEXT
            ISC_GKK(PLACEEXT+nq,nx)=ISR_GKK(nq,nx)
          ENDDO
C$OMP     END PARALLEL DO
        ELSE
          ERROR='>>Invalid sparsity pattern'
          GOTO 9999
        ENDIF
        NX_LIST(0)=ntemp
      ELSE
        IF(SPARSEGKK(nx).EQ.0) THEN
          DO nq=NQR(1,nr),NQR(2,nr)
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),RDUMMY,PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                GKK(((NQGP(nzero,nq)-1)*NQT)+nq,nx)=
     '            COEFFSEXT(NQGP_PIVOT(nzero,nq))
              ENDDO !nz
            ELSE
              DO nzero=1,NQLIST(0)
                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                nq2=NQLIST2(nzero)
                GKK(((nq2-1)*NQT)+nq,nx)=D_MATRIX(nqq,nzero)
              ENDDO !nz
            ENDIF
          ENDDO !nq
        ELSE IF(SPARSEGKK(nx).EQ.1 .OR. SPARSEGKK(nx).EQ.5) THEN
          PLACEEXT=0
          DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),RDUMMY,PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                PLACEEXT=PLACEEXT+1
                ISC_GKK(PLACEEXT,nx)=nq
                ISR_GKK(PLACEEXT,nx)=NQGP(nzero,nq)
                GKK(PLACEEXT,nx)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
              ENDDO !nzero
            ELSE
              DO nzero=1,NQLIST(0)
                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                nq2=NQLIST2(nzero)
                PLACEEXT=PLACEEXT+1
                ISC_GKK(PLACEEXT,nx)=nq
                ISR_GKK(PLACEEXT,nx)=nq2
                GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
              ENDDO !nzero
            ENDIF
          ENDDO !nq
          CALL ASSERT(PLACEEXT.LE.NISR_GKKM,'>>Increase NISR_GKKM',
     '      ERROR,*9999)

          NZZT(1,nr,nx)=PLACEEXT
C$OMP     PARALLEL DO
C$OMP&    PRIVATE(nq),
C$OMP&    SHARED(ISC_GKK,ISR_GKK,nx,PLACEEXT)
          DO nq=1,PLACEEXT
            ISC_GKK(PLACEEXT+nq,nx)=ISR_GKK(nq,nx)
          ENDDO
C$OMP     END PARALLEL DO
        ELSE
          ERROR='>>Invalid sparsity pattern'
          GOTO 9999
        ENDIF
      ENDIF

C SEN Convert from sparse 5 to sparse 1 -- hack until the code is
C     re-written.
      IF(SPARSEGKK(nx).EQ.1 .OR. SPARSEGKK(nx).EQ.5) THEN
        CALL CONVERT_SPARSE_5TO1(NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '    NZZT(1,nr,nx),ISC_GKK(1,nx),ISR_GKK(1,nx),ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time assemble matrix : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '(/'' Global stiffness matrix GKK - external:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '    //''', NOT(2,1,nr,nx)='',I5)') NOT(1,1,nr,nx),
     '    NOT(2,1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL OPSTFMAT(DUMMY_LIST,ISC_GKK(1,nx),ISR_GKK(1,nx),IOOP,
     '    NOT(1,1,nr,nx),NOT(2,1,nr,nx),
     '    NZZT(1,nr,nx),DUMMY_LIST,SPARSEGKK(nx),
     '    GKK(1,nx),GKK(1,nx),'GKK','GKK',.TRUE.,.TRUE.,.FALSE.,
     '    ERROR,*9999)
      ENDIF

      IF(IWRIT4(nr,nx).EQ.-3) THEN !Output global matrices
        IDUMMY=KTYP4
        KTYP4=1
        WRITE(LOCFILENAME,'(A)') FILE00
        WRITE(FILE00,'(''extcGKK'')')
        CALL WRITE_SOL_MATRIX(ISC_GKK(1,nx),ISR_GKK(1,nx),nr,nx,
     '    GKK(1,nx),RDUMMY,ERROR,*9999)
        WRITE(FILE00,'(A)') LOCFILENAME
        KTYP4=IDUMMY
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START2)
      nzero=0
      DO nq=1,NZZT(1,nr,nx)
        IF(DABS(GKK(nq,nx)).GT.ZERO_TOL) nzero=nzero+1
      ENDDO
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)

      WRITE(OP_STRING,'('' Time to count nzs : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      WRITE(OP_STRING,'(/'' Number of entries : '',I10)') NZZT(1,nr,nx)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Number of nonzeros: '',I10)') nzero
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_STOP1)
      ELAPSED_TIME=TIME_STOP1(1)-TIME_START1(1)
      WRITE(OP_STRING,'(/'' Time for ASSEMBLE11_2_3 : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('ASSEMBLE11_2_3')
      RETURN
 9999 CALL ERRORS('ASSEMBLE11_2_3',ERROR)
      CALL EXITS('ASSEMBLE11_2_3')
      RETURN 1
      END


