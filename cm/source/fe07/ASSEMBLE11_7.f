      SUBROUTINE ASSEMBLE11_7(ISC_GKK,ISR_GKK,NEELEM,NENQ,
     '  NP_INTERFACE,NPLIST,NPNODE,NQGP,NQGP_PIVOT,NQLIST,NQLIST2,
     '  NQLIST3,NQLIST4,NQNP,NQNP_LOCAL,NQNP_LOCAL_PIV,NQS,NQXI,nr,nr2,
     '  nr3,nr4,SYSSIZE,NWQ,nx,nx2,nx3,nx4,NXQ,CQ,D_MATRIX,DNUDXQ,
     '  DXDXIQ,DXDXIQ2,GCHQ,GK,GK2,GK3,GKK,GQ,GQ2,GQ3,GQ_D_MATRIX,
     '  GKGK2,GQGQ2,GUQ,NQGW,PROPQ,YQ,FIXQ,SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: ASSEMBLE11_7
C###  Description:
C###    <HTML> <PRE>
C###    ASSEMBLE11_7 creates matrices for coupled BEM-FD
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
     '  NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQLIST(0:NQM),
     '  NQLIST2(0:NQM),
     '  NQLIST3(0:NQM),NQLIST4(0:NQM),NQNP(NPM),NQNP_LOCAL(NPM),
     '  NQNP_LOCAL_PIV(NPM),NQS(NEQM),NQXI(0:NIM,NQSCM),nr,nr2,nr3,nr4,
     '  NWQ(8,0:NQM),nx,nx2,nx3,nx4,NXQ(-NIM:NIM,0:4,0:NQM),SYSSIZE
      REAL*8 CQ(NMM,NQM),D_MATRIX(NPM,*),DNUDXQ(3,3,NQM),
     '  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GK(NZ_GK_M),
     '  GK2(*),GK3(*),GKK(NZ_GKK_M,NXM),GQ(NZ_GQ_M),GQ2(*),GQ3(*),
     '  GQ_D_MATRIX(NPM,*),GUQ(3,3,NQM),NQGW(NQGM,NQM),GKGK2(NPM,NPM),
     '  GQGQ2(NPM,NPM),PROPQ(3,3,4,2,NQM),YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM),SOLVEEIGHTPROBLEM
!     Local variables
      INTEGER BCINNBND,DUMMY_LIST(0:2),F1,F2,IDUMMY,ncol,NCOUPNODES,
     '  NITB,np,npp,nq,nq2,nqq,nrow,ntemp,nzero,nzz,PLACEEXT,RMVDCOL,
     '  S1,S2
      REAL*8 ANS,COEFFSEXT(NQGM),PVTTMP,RDUMMY(NQGM)
      REAL TIME_START1(1),TIME_START2(1),TIME_STOP1(1),TIME_STOP2(1),
     '  ELAPSED_TIME
      CHARACTER LOCFILENAME*(MXCH)
      LOGICAL ERROR_FLAG,INLIST

      CALL ENTERS('ASSEMBLE11_7',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START1)
      NITB=NQXI(0,NQS(NEELEM(1,nr)))

C      IF(SALU_CONSISTENCY(nx2).AND.
C     '  (NOT(1,1,nr2,nx2).LT.NPNODE(0,nr2))) THEN
C        NOT(1,1,nr2,nx2)=NOT(1,1,nr2,nx2)+1
C        NOT(2,1,nr2,nx2)=NOT(2,1,nr2,nx2)+1
C      ENDIF
C      IF(SALU_CONSISTENCY(nx3).AND.
C     '  (NOT(1,1,nr3,nx3).LT.NPNODE(0,nr3))) THEN
C        NOT(1,1,nr3,nx3)=NOT(1,1,nr3,nx3)+1
C        NOT(2,1,nr3,nx3)=NOT(2,1,nr3,nx3)+1
C      ENDIF
C      IF(SALU_CONSISTENCY(nx4).AND.
C     '  (NOT(1,1,nr4,nx4).LT.NPNODE(0,nr4))) THEN
C        NOT(1,1,nr4,nx4)=NOT(1,1,nr4,nx4)+1
C        NOT(2,1,nr4,nx4)=NOT(2,1,nr4,nx4)+1
C      ENDIF

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' Regions: '',4I4)') nr,nr2,nr3,nr4
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Classes: '',4I4)') nx,nx2,nx3,nx4
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' GKs: '',3F12.6)') GK(1),GK2(1),GK3(1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' GQs: '',3F12.6)') GQ(1),GQ2(1),GQ3(1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

!     Get a list of torso surface node numbers
      NPLIST(0)=0
      NCOUPNODES=NQLIST(0)
      DO npp=1,NPNODE(0,nr4)
        np=NPNODE(npp,nr4)
        IF(NP_INTERFACE(np,0).EQ.1) THEN
          NPLIST(0)=NPLIST(0)+1
          NPLIST(NPLIST(0))=np
          NQNP(np)=NQT+NPLIST(0)
          NQLIST(0)=NQLIST(0)+1
          NQLIST(NQLIST(0))=NQNP(np)
        ENDIF
      ENDDO

      DO nq=0,NQLIST(0)
        NQLIST4(nq)=NQLIST(nq)
      ENDDO
      DO nq=1,CPLST(0,1)
        NQLIST4(0)=NQLIST4(0)+1
        NQLIST4(NQLIST4(0))=NQNP(CPLST(nq,1))
      ENDDO
      CALL ISORT(NQLIST4(0),NQLIST4(1))

      SYSSIZE=NQT+NPLIST(0)
      CALL ASSERT(NPM.GE.NQLIST(0),'>> Increase NPM',ERROR,*9999)
      CALL ASSERT(NQM.GE.SYSSIZE,'>> Increase NQM',ERROR,*9999)

!     Diagnostics
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' NQT: '',I8)') NQT
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Coupled nodes: '',I8)') NCOUPNODES
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO npp=1,NPLIST(0)
          np=NPLIST(npp)
          WRITE(OP_STRING,'('' Torso Node: '',I8,'' Grid Num: '',I8)')
     '      np,NQNP(np)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NQLIST(0)
          WRITE(OP_STRING,'('' Grid Num: '',I8,'' NQLIST:'',I8)')
     '      nq,NQLIST(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(/'' System Size: '',I8)') SYSSIZE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Reference node: '',I8)') SOL_ACT_FIX_NODE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

!     Calc system sizes
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      DO nq=1,NQT
        NWQ(3,nq)=0
      ENDDO !nq

      DO nqq=1,NCOUPNODES
        nq=NQLIST(nqq)
        DO ntemp=1,NQGP(0,nq)
          nq2=NQGP(ntemp,nq)
          NWQ(3,nq2)=2
        ENDDO
      ENDDO !nq

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

      DO nq=1,NQLIST4(0)
        NQLIST2(0)=NQLIST2(0)+1
        NQLIST2(NQLIST2(0))=NQLIST4(nq)
      ENDDO
      CALL ISORT(NQLIST2(0),NQLIST2(1))

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
        DO nqq=1,NPNODE(0,nr2)
          np=NPNODE(nqq,nr2)
          WRITE(OP_STRING,'(2I8)') np,NQNP(np)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'('' Region = '',I8)') nr3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nqq=1,NPNODE(0,nr3)
          np=NPNODE(nqq,nr3)
          WRITE(OP_STRING,'(2I8)') np,NQNP(np)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'('' Region = '',I8)') nr3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nqq=1,NPNODE(0,nr4)
          np=NPNODE(nqq,nr4)
          WRITE(OP_STRING,'(2I8)') np,NQNP(np)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

!     Create the D matrix
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      nrow=NPNODE(0,nr2)+NPNODE(0,nr3)+NPNODE(0,nr4)
      ncol=NPNODE(0,nr2)+NPNODE(0,nr3)+NPNODE(0,nr4)+BCINNBND
      CALL ASSERT(nrow.EQ.NQLIST4(0),' >>Matrix row size mismatch',
     '  ERROR,*9999)
      CALL ASSERT(ncol.EQ.NQLIST2(0),' >>Matrix col size mismatch',
     '  ERROR,*9999)
      DO nq=1,nrow
        DO nq2=1,ncol
          D_MATRIX(nq,nq2)=0.0d0
        ENDDO
      ENDDO

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Nrow,Ncol '',2I8)') nrow,ncol
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Nrow='',I6,''+'',I6,''+'',I6)')
     '    NPNODE(0,nr2),NPNODE(0,nr3),NPNODE(0,nr4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Ncol='',I6,''+'',I6,''+'',I6,''+'',I6)')
     '    NPNODE(0,nr2),NPNODE(0,nr3),NPNODE(0,nr4),BCINNBND
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DUMMY_LIST(0)=2
      DO nq=1,CPLST(0,1)
        DUMMY_LIST(nq)=NQNP(CPLST(nq,1))
      ENDDO

      ERROR_FLAG=.FALSE.
      NCOUPNODES=NCOUPNODES+2
      IF(SOLVEEIGHTPROBLEM) THEN
        DO nqq=1,NCOUPNODES
          IF(.NOT.ERROR_FLAG) THEN
            nq=NQLIST4(nqq)
            IF(.NOT.INLIST(nq,DUMMY_LIST(1),DUMMY_LIST(0),IDUMMY)) THEN
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,NXQ,NQGW(1,nq),
     '          CQ(1,nq),DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,*100)
              DO nzero=1,NQGP(0,nq)
                D_MATRIX(nqq,NQLIST3(NQGP(nzero,nq)))=
     '            NQGW(NQGP_PIVOT(nzero,nq),nq)
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
        DO nqq=1,NCOUPNODES
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
      NCOUPNODES=NCOUPNODES-2

!     Diagnostics
      IF(DOP) THEN
        DO nq=1,nrow
          DO nq2=1,ncol
            IF(DABS(D_MATRIX(nq,nq2)).GT.ZERO_TOL) THEN
              WRITE(OP_STRING,'(''row,col,D '',2I8,F12.6)')
     '          nq,nq2,D_MATRIX(nq,nq2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
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
        ENDDO
      ENDDO
      S1=NYT(1,1,nx2)+NYT(1,1,nx3)
      F1=NYT(1,1,nx2)+NYT(1,1,nx3)+NYT(1,1,nx4)
      S2=NYT(2,1,nx2)+NYT(2,1,nx3)
      F2=NYT(2,1,nx2)+NYT(2,1,nx3)+NYT(2,1,nx4)
      DO nq=S1+1,F1
        DO nzero=S2+1,F2
          GQGQ2(nq,nzero)=
     '      GQ3((((nzero-S2)-1)*NYT(1,1,nx4))+
     '      (nq-S1))
          GKGK2(nq,nzero)=
     '      GK3((((nzero-S2)-1)*NYT(1,1,nx4))+
     '      (nq-S1))
        ENDDO
      ENDDO

!     Diagnostics
      IF(DOP) THEN
        F1=NYT(1,1,nx2)+NYT(1,1,nx3)+NYT(1,1,nx4)
        F2=NYT(2,1,nx2)+NYT(2,1,nx3)+NYT(2,1,nx4)
        DO nq=1,F1
          DO nq2=1,F2
            WRITE(OP_STRING,'('' row,col,GKGK2,GQGQ2 '',2I8,2F12.6)')
     '        nq,nq2,GKGK2(nq,nq2),GQGQ2(nq,nq2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time duplicate GK,GQ : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

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
      S1=NPNODE(0,nr2)+NPNODE(0,nr3)
      F1=NPNODE(0,nr2)+NPNODE(0,nr3)+NPNODE(0,nr4)
      DO nqq=S1+1,F1
        np=NPNODE((nqq-S1),nr4)
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
      S1=NPNODE(0,nr2)+NPNODE(0,nr3)
      F1=NPNODE(0,nr2)+NPNODE(0,nr3)+NPNODE(0,nr4)
      DO nqq=S1+1,F1
        np=NPNODE((nqq-S1),nr4)
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

!     Create GQ*D
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C$OMP PARALLEL DO
C$OMP&PRIVATE(ANS,nq,nq2,nzero),
C$OMP&SHARED(D_MATRIX,GQ_D_MATRIX,GQGQ2,NPLIST,nrow,NQLIST2,NQLIST4)
      DO nq=1,nrow
C        DO nq2=1,ncol
        DO nq2=1,(NQLIST2(0)-NPLIST(0))
          ANS=0
C          DO nzero=1,nrow
          DO nzero=1,(NQLIST4(0)-NPLIST(0))
            ANS=ANS+GQGQ2(nq,nzero)*D_MATRIX(nzero,nq2)
          ENDDO
          GQ_D_MATRIX(nq,nq2)=ANS
        ENDDO
        DO nq2=(NQLIST2(0)-NPLIST(0))+1,NQLIST2(0)
          GQ_D_MATRIX(nq,nq2)=0.0d0
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to create GQ*D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!     Diagnostics
      IF(DOP) THEN
        DO nq=1,NQLIST4(0)
          DO nq2=1,NQLIST2(0)
            WRITE(OP_STRING,'('' row,col,GQD '',2I8,F12.6)')
     '        nq,nq2,GQ_D_MATRIX(nq,nq2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

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
      CALL CPU_TIMER(CPU_USER,TIME_STOP2)
      ELAPSED_TIME=TIME_STOP2(1)-TIME_START2(1)
      WRITE(OP_STRING,'('' Time to create GK-GQ*D : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      IF(DOP) THEN
        DO nq=1,NQLIST4(0)
          DO nq2=1,NQLIST2(0)
            WRITE(OP_STRING,'('' row,col,GKGQD '',2I8,F12.6)')
     '        nq,nq2,D_MATRIX(nq,nq2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

!     Initialising
      CALL CPU_TIMER(CPU_USER,TIME_START2)
      RMVDCOL=NQNP(SOL_ACT_FIX_NODE)
      IF(SPARSEGKK(nx).EQ.0) THEN
        nzz=(SYSSIZE-1)*(SYSSIZE-1)
        DO nq=1,nzz
          GKK(nq,nx)=0.0d0
        ENDDO !nq
        NOT(1,1,nr,nx)=SYSSIZE-1
        NOT(2,1,nr,nx)=SYSSIZE-1
        NZZT(1,nr,nx)=nzz
        IF(NZ_GKK_M.LT.nzz) THEN
          WRITE(ERROR,'(''>>Increase NZ_GKK_M to >= '',I12)') nzz
          GOTO 9999
        ENDIF
      ELSE IF(SPARSEGKK(nx).EQ.1) THEN !row/col #2
        nzz=0
        DO nq=1,SYSSIZE
          DO nzero=1,NQGP(0,nq)
            nzz=nzz+1
          ENDDO !nzero
        ENDDO !nq
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
        WRITE(OP_STRING,'('' Solve8 for testing only'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ntemp=NX_LIST(0)
        NX_LIST(0)=1
        IF(SPARSEGKK(nx).EQ.0) THEN
          DO nq=1,SYSSIZE
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),NQGW(1,nq),PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                GKK(((NQGP(nzero,nq)-1)*(SYSSIZE-1))+nq,nx)=
     '            COEFFSEXT(NQGP_PIVOT(nzero,nq))
              ENDDO !nz
            ELSE IF(nq.EQ.RMVDCOL) THEN
              !do nothing
              DO nzero=1,NQM
                YQ(nzero,1,1,nx3)=0.0d0
              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                YQ(NQLIST2(nzero),1,1,nx3)=D_MATRIX(nqq,nzero)
              ENDDO
            ELSE
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              IF(nq.LT.RMVDCOL) THEN
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      GKK(((nq2-1)*(SYSSIZE-1))+nq,nx)=
     '                  D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      GKK(((nq2-2)*(SYSSIZE-1))+nq,nx)=
     '                  D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ELSE !nq > rmvdcol
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      GKK(((nq2-1)*(SYSSIZE-1))+(nq-1),nx)=
     '                  D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      GKK(((nq2-2)*(SYSSIZE-1))+(nq-1),nx)=
     '                  D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ENDIF
            ENDIF
          ENDDO !nq
          SYSSIZE=SYSSIZE-1
        ELSE IF(SPARSEGKK(nx).EQ.1) THEN
          PLACEEXT=0
          DO nq=1,SYSSIZE!create line for each nq
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
            ELSE IF(nq.EQ.RMVDCOL) THEN
              !do nothing
              DO nzero=1,NQM
                YQ(nzero,1,1,nx3)=0.0d0
              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                YQ(NQLIST2(nzero),1,1,nx3)=D_MATRIX(nqq,nzero)
              ENDDO
            ELSE
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              IF(nq.LT.RMVDCOL) THEN
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq
                      ISR_GKK(PLACEEXT,nx)=nq2
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq
                      ISR_GKK(PLACEEXT,nx)=nq2-1
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ELSE !nq > rmvdcol
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq-1
                      ISR_GKK(PLACEEXT,nx)=nq2
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq-1
                      ISR_GKK(PLACEEXT,nx)=nq2-1
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ENDIF
            ENDIF
          ENDDO !nq

          NZZT(1,nr,nx)=PLACEEXT
          DO nq=1,PLACEEXT
            ISC_GKK(PLACEEXT+nq,nx)=ISR_GKK(nq,nx)
          ENDDO
          SYSSIZE=SYSSIZE-1
          NOT(1,1,nr,nx)=SYSSIZE
          NOT(2,1,nr,nx)=SYSSIZE
        ELSE
          ERROR='>>Invalid sparsity pattern'
          GOTO 9999
        ENDIF
        NX_LIST(0)=ntemp
      ELSE
        IF(SPARSEGKK(nx).EQ.0) THEN
          DO nq=1,SYSSIZE
            IF(.NOT.INLIST(nq,NQLIST(1),NQLIST(0),IDUMMY)) THEN
              CALL CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,NQXI,NWQ(1,nq),
     '          NXQ,nx,nx,COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          GCHQ(1,nq),GUQ(1,1,nq),RDUMMY,PROPQ(1,1,1,1,nq),
     '          .TRUE.,FIXQ,.FALSE.,.TRUE.,ERROR,*9999)
              DO nzero=1,NQGP(0,nq)
                GKK(((NQGP(nzero,nq)-1)*(SYSSIZE-1))+nq,nx)=
     '            COEFFSEXT(NQGP_PIVOT(nzero,nq))
              ENDDO !nz
            ELSE IF(nq.EQ.RMVDCOL) THEN
              !do nothing
              DO nzero=1,NQM
                YQ(nzero,1,1,nx3)=0.0d0
              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                YQ(NQLIST2(nzero),1,1,nx3)=D_MATRIX(nqq,nzero)
              ENDDO
            ELSE
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              IF(nq.LT.RMVDCOL) THEN
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      GKK(((nq2-1)*(SYSSIZE-1))+nq,nx)=
     '                  D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      GKK(((nq2-2)*(SYSSIZE-1))+nq,nx)=
     '                  D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ELSE !nq > rmvdcol
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      GKK(((nq2-1)*(SYSSIZE-1))+(nq-1),nx)=
     '                  D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      GKK(((nq2-2)*(SYSSIZE-1))+(nq-1),nx)=
     '                  D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ENDIF
            ENDIF
          ENDDO !nq
          SYSSIZE=SYSSIZE-1
        ELSE IF(SPARSEGKK(nx).EQ.1) THEN
          PLACEEXT=0
          DO nq=1,SYSSIZE !create line for each nq
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
            ELSE IF(nq.EQ.RMVDCOL) THEN
              !do nothing
              DO nzero=1,NQM
                YQ(nzero,1,1,nx3)=0.0d0
              ENDDO
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              DO nzero=1,NQLIST2(0)
                YQ(NQLIST2(nzero),1,1,nx3)=D_MATRIX(nqq,nzero)
              ENDDO
            ELSE
              DO nzero=1,NQLIST4(0)
                IF(nq.EQ.NQLIST4(nzero)) nqq=nzero
              ENDDO
              IF(nq.LT.RMVDCOL) THEN
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq
                      ISR_GKK(PLACEEXT,nx)=nq2
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq
                      ISR_GKK(PLACEEXT,nx)=nq2-1
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ELSE !nq > rmvdcol
                DO nzero=1,NQLIST2(0)
                  IF(DABS(D_MATRIX(nqq,nzero)).GT.ZERO_TOL) THEN
                    nq2=NQLIST2(nzero)
                    IF(nq2.EQ.RMVDCOL) THEN
                      !do nothing
                    ELSE IF(nq2.LT.RMVDCOL) THEN
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq-1
                      ISR_GKK(PLACEEXT,nx)=nq2
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ELSE !nq2 > rmvdcol
                      PLACEEXT=PLACEEXT+1
                      ISC_GKK(PLACEEXT,nx)=nq-1
                      ISR_GKK(PLACEEXT,nx)=nq2-1
                      GKK(PLACEEXT,nx)=D_MATRIX(nqq,nzero)
                    ENDIF
                  ENDIF
                ENDDO !nzero
              ENDIF
            ENDIF
          ENDDO !nq

          CALL ASSERT(PLACEEXT.LE.NISR_GKKM,'>>Increase NISR_GKKM',
     '      ERROR,*9999)

          NZZT(1,nr,nx)=PLACEEXT
          DO nq=1,PLACEEXT
            ISC_GKK(PLACEEXT+nq,nx)=ISR_GKK(nq,nx)
          ENDDO
          SYSSIZE=SYSSIZE-1
          NOT(1,1,nr,nx)=SYSSIZE
          NOT(2,1,nr,nx)=SYSSIZE
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

C      IF(SALU_CONSISTENCY(nx2).AND.
C     '  (NOT(1,1,nr2,nx2).LT.NPNODE(0,nr2))) THEN
C        NOT(1,1,nr2,nx2)=NOT(1,1,nr2,nx2)-1
C        NOT(2,1,nr2,nx2)=NOT(2,1,nr2,nx2)-1
C      ENDIF
C      IF(SALU_CONSISTENCY(nx3).AND.
C     '  (NOT(1,1,nr3,nx3).LT.NPNODE(0,nr3))) THEN
C        NOT(1,1,nr3,nx3)=NOT(1,1,nr3,nx3)-1
C        NOT(2,1,nr3,nx3)=NOT(2,1,nr3,nx3)-1
C      ENDIF
C      IF(SALU_CONSISTENCY(nx4).AND.
C     '  (NOT(1,1,nr4,nx4).LT.NPNODE(0,nr4))) THEN
C        NOT(1,1,nr4,nx4)=NOT(1,1,nr4,nx4)-1
C        NOT(2,1,nr4,nx4)=NOT(2,1,nr4,nx4)-1
C      ENDIF

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
      WRITE(OP_STRING,'(/'' Time for ASSEMBLE11_7 : '',F12.6)')
     '  ELAPSED_TIME
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('ASSEMBLE11_7')
      RETURN
 9999 CALL ERRORS('ASSEMBLE11_7',ERROR)
      CALL EXITS('ASSEMBLE11_7')
      RETURN 1
      END


