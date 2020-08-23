      SUBROUTINE ASSEMBLE11_4_5(ISC_GKK,ISR_GKK,NEELEM,NENQ,NPNODE,
     '  NQGP,NQGP_PIVOT,NQLIST,NQLIST2,NQLIST3,NQLIST4,NQNP,NQNP_LOCAL,
     '  NQNP_LOCAL_PIV,NQS,NQXI,nr,nr2,nr3,NWQ,nx,nx2,nx3,NXQ,CQ,
     '  D_MATRIX,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GK,GK2,GKK,GQ,GQ2,
     '  GQ_D_MATRIX,GKGK2,GQGQ2,GUQ,NQGW,PROPQ,FIXQ,
     '  SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: ASSEMBLE11_4_5
C###  Description:
C###    <HTML> <PRE>
C###    ASSEMBLE11_4_5 creates matrices for coupled BEM-FD
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
     '  NQLIST(0:NQM),NQLIST2(0:NQM),NQLIST3(0:NQM),NQLIST4(0:NQM),
     '  NQNP(NPM),
     '  NQNP_LOCAL(NPM),NQNP_LOCAL_PIV(NPM),NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  nr,nr2,nr3,NWQ(8,0:NQM),nx,nx2,nx3,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 CQ(NMM,NQM),D_MATRIX(NPM,*),DNUDXQ(3,3,NQM),
     '  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GK(NZ_GK_M),
     '  GK2(*),GKK(NZ_GKK_M,NXM),GQ(NZ_GQ_M),GQ2(*),GQ_D_MATRIX(NPM,*),
     '  GUQ(3,3,NQM),NQGW(NQGM,NQM),GKGK2(NPM,NPM),GQGQ2(NPM,NPM),
     '  PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM),SOLVEEIGHTPROBLEM
!     Local variables
      INTEGER BCINNBND,DUMMY_LIST(0:2),IDUMMY,ncol,NITB,np,nq,nq2,nqq,
     '  nrow,ntemp,nzero,nzz,PLACEEXT
      REAL*8 ANS,COEFFSEXT(NQGM),PVTTMP,RDUMMY(NQGM)
      REAL TIME_START1(1),TIME_START2(1),TIME_STOP1(1),TIME_STOP2(1),
     '  ELAPSED_TIME
      CHARACTER LOCFILENAME*(MXCH)
      LOGICAL ERROR_FLAG,INLIST

      CALL ENTERS('ASSEMBLE11_4_5',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      NITB=NQXI(0,NQS(NEELEM(1,nr)))

!     Diagnostics
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
        WRITE(OP_STRING,'('' Boundary Element GK,GQ '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Region '',I6,'' Class '',I6)') nr3,nx3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Rows '',I6,'' Columns '',I6)')
     '    NYT(1,1,nx3),NYT(2,1,nx3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nzz=1,NYT(1,1,nx3)*NYT(2,1,nx3)
          WRITE(OP_STRING,'(I8,'' GK '',F12.6,'' GQ '',F12.6)')
     '      nzz,GK2(nzz),GQ2(nzz)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nzz
      ENDIF

!     Calc system sizes
      CALL CPU_TIMER(CPU_USER,TIME_START2)

!     nqlist4 is just nqlist with the cusps added
      DO nq=0,NQLIST(0)
        NQLIST4(nq)=NQLIST(nq)
      ENDDO
      DO nq=1,CPLST(0,1)
        NQLIST4(0)=NQLIST4(0)+1
        NQLIST4(NQLIST4(0))=NQNP(CPLST(nq,1))
      ENDDO
      CALL ISORT(NQLIST4(0),NQLIST4(1))

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Connected grid points '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NQLIST4(0)
          WRITE(OP_STRING,'('' Index '',I8,'' Grid '',I8)') nq,
     '      NQLIST4(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      DO nq=1,NQT
        NWQ(3,nq)=0
      ENDDO !nq

!     In this loop, I don't think NQLIST4 should be used
!     because the points stored are not flux points.
      DO nqq=1,NQLIST(0)
        nq=NQLIST(nqq)
        DO ntemp=1,NQGP(0,nq)
          nq2=NQGP(ntemp,nq)
          NWQ(3,nq2)=2
        ENDDO
      ENDDO !nq

!     changed to nqlist4
      DO nqq=1,NQLIST4(0)
        nq=NQLIST4(nqq)
        NWQ(3,nq)=0
      ENDDO !nq

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

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Additional internal connected points '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NQLIST2(0)
          WRITE(OP_STRING,'('' Index '',I8,'' Grid '',I8)') nq,
     '      NQLIST2(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

!     Changed to nqlist4
      DO nq=1,NQLIST4(0)
        NQLIST2(0)=NQLIST2(0)+1
        NQLIST2(NQLIST2(0))=NQLIST4(nq)
      ENDDO
      CALL ISORT(NQLIST2(0),NQLIST2(1))

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' NQLIST2 - D columns '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NQLIST2(0)
          WRITE(OP_STRING,'('' Index '',I8,'' Grid '',I8)') nq,
     '      NQLIST2(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      DO nq=0,NQM
        NQLIST3(nq)=0
      ENDDO
      DO nq=1,NQLIST2(0)
        NQLIST3(NQLIST2(nq))=nq
      ENDDO

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Inner Bound conn nq = '',I8)') BCINNBND
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Region = '',I8)') nr2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Nodes =      Grid ='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nqq=1,NPNODE(0,nr2)
          np=NPNODE(nqq,nr2)
          WRITE(OP_STRING,'(2I8)') np,NQNP(np)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'('' Region = '',I8)') nr3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Nodes =      Grid ='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nqq=1,NPNODE(0,nr3)
          np=NPNODE(nqq,nr3)
          WRITE(OP_STRING,'(2I8)') np,NQNP(np)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

!     Create the D matrix
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      nrow=NQLIST4(0)
      ncol=NQLIST2(0)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Number of rows in D '',I8)') nrow
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Number of cols in D '',I8)') ncol
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO nq=1,nrow
        DO nq2=1,ncol
          D_MATRIX(nq,nq2)=0.0d0
        ENDDO
      ENDDO

      DUMMY_LIST(0)=2
      DO nq=1,CPLST(0,1)
        DUMMY_LIST(nq)=NQNP(CPLST(nq,1))
      ENDDO

      ERROR_FLAG=.FALSE.
      IF(SOLVEEIGHTPROBLEM) THEN
        DO nqq=1,nrow
          IF(.NOT.ERROR_FLAG) THEN
            nq=NQLIST4(nqq)
            IF(.NOT.INLIST(nq,DUMMY_LIST(1),DUMMY_LIST(0),IDUMMY)) THEN
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,NXQ,NQGW(1,nq),
     '          CQ(1,nq),DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,*100)
              DO nzero=1,NQGP(0,nq)
                D_MATRIX(nqq,NQLIST3(NQGP(nzero,nq)))=
     '            NQGW(NQGP_PIVOT(nzero,nq),nq)
C                D_MATRIX(nqq,NQLIST3(NQGP(nzero,nq)))=DBLE(nzero)
              ENDDO
            ENDIF
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
      ELSE
        DO nqq=1,nrow
          IF(.NOT.ERROR_FLAG) THEN
            nq=NQLIST4(nqq)
            IF(.NOT.INLIST(nq,DUMMY_LIST(1),DUMMY_LIST(0),IDUMMY)) THEN
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,NXQ,COEFFSEXT,
     '          CQ(6,nq),DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,*200)
              DO nzero=1,NQGP(0,nq)
                D_MATRIX(nqq,NQLIST3(NQGP(nzero,nq)))=
     '            COEFFSEXT(NQGP_PIVOT(nzero,nq))
              ENDDO
            ENDIF
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
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' D matrix values '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,nrow
          DO nq2=1,ncol
            WRITE(OP_STRING,'('' Row '',I6,'' Col '',I6,'' Value'
     '        //' '',F12.6)') nq,nq2,D_MATRIX(nq,nq2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time create D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Make a local copy of GK,GQ
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      DO nq=1,NPM
        DO nq2=1,NPM
          GKGK2(nq,nq2)=0.0d0
          GQGQ2(nq,nq2)=0.0d0
        ENDDO
      ENDDO
      DO nq=1,NYT(1,1,nx2)
        DO nzero=1,NYT(2,1,nx2)
          GQGQ2(nq,nzero)=
     '      GQ(((nzero-1)*NYT(1,1,nx2))+nq)
          GKGK2(nq,nzero)=
     '      GK(((nzero-1)*NYT(1,1,nx2))+nq)
C          GKGK2(nq,nzero)=nq
C          GQGQ2(nq,nzero)=nzero
        ENDDO
      ENDDO
      DO nq=NYT(1,1,nx2)+1,NYT(1,1,nx2)+NYT(1,1,nx3)
        DO nzero=NYT(2,1,nx2)+1,NYT(2,1,nx2)+NYT(2,1,nx3)
          GQGQ2(nq,nzero)=
     '      GQ2((((nzero-NYT(2,1,nx2))-1)*NYT(1,1,nx3))+
     '      (nq-NYT(1,1,nx2)))
          GKGK2(nq,nzero)=
     '      GK2((((nzero-NYT(2,1,nx2))-1)*NYT(1,1,nx3))+
     '      (nq-NYT(1,1,nx2)))
C          GKGK2(nq,nzero)=nq
C          GQGQ2(nq,nzero)=nzero
        ENDDO
      ENDDO
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time duplicate GK,GQ : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' GKGK2,GQGQ2 matrices '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NYT(1,1,nx2)+NYT(1,1,nx3)
          DO nzero=1,NYT(2,1,nx2)+NYT(2,1,nx3)
            IF((DABS(GKGK2(nq,nzero)).GT.ZERO_TOL).OR.
     '        (DABS(GQGQ2(nq,nzero)).GT.ZERO_TOL)) THEN
              WRITE(OP_STRING,'('' row,col,GK,GQ '',2I6,2F12.6)')
     '          nq,nzero,GKGK2(nq,nzero),GQGQ2(nq,nzero)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START2)
!     Need to adjust GKGK2,GQGQ2 to reflect grid numbering
!     swap 2 rows at once, element by element so only
!     one temporary memory space is needed. Both the rows
!     and the columns need to be pivoted.
!     Create 1 matrix which can be multiplied by the D matrix

      DO nqq=1,NPNODE(0,nr2)
        np=NPNODE(nqq,nr2)
        NQNP_LOCAL(nqq)=NQNP(np)
      ENDDO
      DO nqq=(NPNODE(0,nr2)+1),(NPNODE(0,nr2)+NPNODE(0,nr3))
        np=NPNODE((nqq-NPNODE(0,nr2)),nr3)
        NQNP_LOCAL(nqq)=NQNP(np)
      ENDDO
      CALL ISORTP(NQLIST4(0),NQNP_LOCAL,NQNP_LOCAL_PIV)

!     Pivot rows
      DO nqq=1,NQLIST4(0)
        IF(NQNP_LOCAL_PIV(nqq).NE.nqq) THEN !must pivot
          DO nq2=nqq+1,NQLIST4(0)
            IF(NQNP_LOCAL_PIV(nq2).EQ.nqq) THEN !found pivot
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Row pivot '',I8,'' with '',I8)')
     '            nqq,NQNP_LOCAL_PIV(nqq)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO nzero=1,NQLIST4(0)
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
                DO nzero=1,NQLIST4(0)
                  WRITE(OP_STRING,'(2I8)') nzero,NQNP_LOCAL_PIV(nzero)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO

!     Pivot Columns
      DO nqq=1,NPNODE(0,nr2)
        np=NPNODE(nqq,nr2)
        NQNP_LOCAL(nqq)=NQNP(np)
      ENDDO
      DO nqq=(NPNODE(0,nr2)+1),(NPNODE(0,nr2)+NPNODE(0,nr3))
        np=NPNODE((nqq-NPNODE(0,nr2)),nr3)
        NQNP_LOCAL(nqq)=NQNP(np)
      ENDDO
      CALL ISORTP(NQLIST4(0),NQNP_LOCAL,NQNP_LOCAL_PIV)

      DO nqq=1,NQLIST4(0)
        IF(NQNP_LOCAL_PIV(nqq).NE.nqq) THEN !must pivot
          DO nq2=nqq+1,NQLIST4(0)
            IF(NQNP_LOCAL_PIV(nq2).EQ.nqq) THEN !found pivot
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Col pivot '',I8,'' with '',I8)')
     '            nqq,NQNP_LOCAL_PIV(nqq)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO nzero=1,NQLIST4(0)
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
                DO nzero=1,NQLIST4(0)
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

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' GKGK2,GQGQ2 pivoted matrices '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NYT(1,1,nx2)+NYT(1,1,nx3)
          DO nzero=1,NYT(2,1,nx2)+NYT(2,1,nx3)
            IF((DABS(GKGK2(nq,nzero)).GT.ZERO_TOL).OR.
     '        (DABS(GQGQ2(nq,nzero)).GT.ZERO_TOL)) THEN
              WRITE(OP_STRING,'('' row,col,GK,GQ '',2I6,2F12.6)')
     '          nq,nzero,GKGK2(nq,nzero),GQGQ2(nq,nzero)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!     Diagnostics
      IF(DOP) THEN
        DO nq2=1,NQLIST4(0)
          nq=NQLIST4(nq2)
          DO nqq=1,NPNODE(0,nr2)
            np=NPNODE(nqq,nr2)
            IF(NQNP(np).EQ.nq) nzero=np
          ENDDO
          DO nqq=(NPNODE(0,nr2)+1),(NPNODE(0,nr2)+NPNODE(0,nr3))
            np=NPNODE((nqq-NPNODE(0,nr2)),nr3)
            IF(NQNP(np).EQ.nq) nzero=np
          ENDDO
          WRITE(OP_STRING,'('' Grid '',I8,'' Node '',I8)') nq,nzero
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

!     Create GQ*D
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C$OMP PARALLEL DO
C$OMP&PRIVATE(ANS,nq,nq2,nzero),
C$OMP&SHARED(D_MATRIX,GQ_D_MATRIX,GQGQ2,ncol,nrow)
      DO nq=1,nrow
        DO nq2=1,ncol
          ANS=0
          DO nzero=1,nrow
            ANS=ANS+GQGQ2(nq,nzero)*D_MATRIX(nzero,nq2)
          ENDDO
          GQ_D_MATRIX(nq,nq2)=ANS
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Rows in GQD '',I8)') nrow
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Cols in GQD '',I8)') ncol
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to create GQ*D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Create GK-(GQ*D)
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      DO nq=1,nrow
        DO nq2=1,ncol
          D_MATRIX(nq,nq2)=GQ_D_MATRIX(nq,nq2)
        ENDDO
      ENDDO
      DO nq=1,nrow
        DO nq2=1,nrow
          np=NQLIST4(nq2)
          nqq=NQLIST3(np)
          D_MATRIX(nq,nqq)=D_MATRIX(nq,nqq)+
     '      GKGK2(nq,nq2)
        ENDDO
      ENDDO

!     Diagnostics
      IF(DOP) THEN
        DO nq=1,nrow
          DO nq2=1,nrow
            np=NQLIST4(nq2)
            nqq=NQLIST3(np)
            WRITE(OP_STRING,'('' GK_GQD indices '',2I8,'' GK '
     '        //'indices '',2I8)') nq,nqq,nq,nq2
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to create GK-GQ*D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Initialising
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      IF(SPARSEGKK(nx).EQ.0) THEN
        nzz=NQT*NQT
        DO nq=1,nzz
          GKK(nq,nx)=0.0d0
        ENDDO !nq
        NOT(1,1,nr,nx)=NQT
        NOT(2,1,nr,nx)=NQT
        NZZT(1,nr,nx)=nzz
        IF(NZ_GKK_M.LT.nzz) THEN
          WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
          GOTO 9999
        ENDIF
      ELSE IF(SPARSEGKK(nx).EQ.1) THEN !row/col #2
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
          DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
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
C             changed from nqlist to nqlist4
C              DO nzero=1,NQLIST(0)
C                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
C              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO

C              WRITE(*,*) 'GKK row,D row',nq,nqq
              DO nzero=1,NQLIST2(0)
                nq2=NQLIST2(nzero)
                GKK(((nq2-1)*NQT)+nq,nx)=D_MATRIX(nqq,nzero)
              ENDDO !nz
            ENDIF
          ENDDO !nq
        ELSE IF(SPARSEGKK(nx).EQ.1) THEN
          PLACEEXT=0
          DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                IF(DABS(COEFFSEXT(NQGP_PIVOT(nzero,nq))).GT.ZERO_TOL)
     '            THEN
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx)=nq
                  ISR_GKK(PLACEEXT,nx)=NQGP(nzero,nq)
                  GKK(PLACEEXT,nx)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDIF
              ENDDO !nzero
            ELSE
C             changed from nqlist to nqlist4
C              DO nzero=1,NQLIST(0)
C                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
C              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                  nq2=NQLIST2(nzero)
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx)=nq
                  ISR_GKK(PLACEEXT,nx)=nq2
                  GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                ENDIF
              ENDDO !nzero
            ENDIF
          ENDDO !nq

          NZZT(1,nr,nx)=PLACEEXT
          DO nq=1,PLACEEXT
            ISC_GKK(PLACEEXT+nq,nx)=ISR_GKK(nq,nx)
          ENDDO
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
C             changed from nqlist to nqlist4
C              DO nzero=1,NQLIST(0)
C                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
C              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                nq2=NQLIST2(nzero)
                GKK(((nq2-1)*NQT)+nq,nx)=D_MATRIX(nqq,nzero)
              ENDDO !nz
            ENDIF
          ENDDO !nq
        ELSE IF(SPARSEGKK(nx).EQ.1) THEN
          PLACEEXT=0
          DO nq=NQR(1,nr),NQR(2,nr) !create line for each nq
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),RDUMMY,PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                IF(DABS(COEFFSEXT(NQGP_PIVOT(nzero,nq))).GT.ZERO_TOL)
     '            THEN
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx)=nq
                  ISR_GKK(PLACEEXT,nx)=NQGP(nzero,nq)
                  GKK(PLACEEXT,nx)=COEFFSEXT(NQGP_PIVOT(nzero,nq))
                ENDIF
              ENDDO !nzero
            ELSE
C             changed from nqlist to nqlist4
C              DO nzero=1,NQLIST(0)
C                IF(nq.EQ.NQLIST(nzero)) nqq=nzero
C              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                  nq2=NQLIST2(nzero)
                  PLACEEXT=PLACEEXT+1
                  ISC_GKK(PLACEEXT,nx)=nq
                  ISR_GKK(PLACEEXT,nx)=nq2
                  GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                ENDIF
              ENDDO !nzero
            ENDIF
          ENDDO !nq

          CALL ASSERT(PLACEEXT.LE.NISR_GKKM,'>>Increase NISR_GKKM',
     '      ERROR,*9999)

          NZZT(1,nr,nx)=PLACEEXT
          DO nq=1,PLACEEXT
            ISC_GKK(PLACEEXT+nq,nx)=ISR_GKK(nq,nx)
          ENDDO
        ELSE
          ERROR='>>Invalid sparsity pattern'
          GOTO 9999
        ENDIF
      ENDIF

C SEN Convert from sparse 5 to sparse 1 -- hack until the code is
C     re-written.
      IF(SPARSEGKK(nx).EQ.1) THEN
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
      WRITE(OP_STRING,'(/'' Time for ASSEMBLE11_4_5 : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('ASSEMBLE11_4_5')
      RETURN
 9999 CALL ERRORS('ASSEMBLE11_4_5',ERROR)
      CALL EXITS('ASSEMBLE11_4_5')
      RETURN 1
      END


