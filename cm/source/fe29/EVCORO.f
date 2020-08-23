      SUBROUTINE EVCORO(BC_POINTS,BRANCH,CALCULATED,CONECT,NBJ,NENQ,NEP,
     '  NPNE,NQET,NQNE,NQS,NXQ,NYNQ,XIP,XQ,YQ,STRING,ERROR,*)
C#### Subroutine: EVCORO
C###  Description:
C###    EVCORO evaluates coronary flow wash out curves,
C###    and regional hetergeneity.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER BC_POINTS(NQM),CONECT(-1:1,0:2,NQM),
     '  NBJ(NJM,NEM),NENQ(0:8,NQM),NEP(NPM),NPNE(NNM,NBFM,NEM),
     '  NQET(NQSCM),NQNE(NEQM,NQEM),NQS(NEQM),NYNQ(NHM,NQM,0:NRCM,NXM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8  XIP(NIM,NPM),XQ(NJM,NQM),YQ(NYQM,NIQM,NAM,NXM)
      LOGICAL CALCULATED(NQM),BRANCH(NQM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ADJACENT,COUNT,CURRENT,i,ii,IBEG,IEND,
     '  IFROMC,j,LABLED,ne,nj,NJTOT,no_bc_points,np,
     '  nq,nq_next,nq1,nq2,nqq,nr,nx,nxc,ny_va,ny_vv,N3CO,
     '  NUM_NODES,POINTS(3)
      REAL*8 DELTA_X,FLOW,TRANS_TIME_ART,TRANS_TIME_VIEN
      LOGICAL CBBREV,FOUND

      CALL ENTERS('EVCORO',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate coronary
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    evaluates washout curves for coronary out flow

        OP_STRING(1)=BLANK(1:15)//'<region>'
        OP_STRING(2)=BLANK(1:15)//'<class>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
      ELSE
        IF(CBBREV(CO,'REGION',2,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=2
        ENDIF

        IF(CBBREV(CO,'CLASS',2,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
        ELSE
          nxc=2
        ENDIF

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)

        DO nq=1,NQT
          CALCULATED(nq)=.FALSE.
        ENDDO
        CALCULATED(NQ_START(nr))=.TRUE. !ipinit NPS 4/2/97
        CURRENT=NXQ(1,1,NQ_START(nr),1)
        IF(CURRENT.EQ.0) THEN
          CURRENT=NXQ(-1,1,NQ_START(nr),1)
        ENDIF
        CONECT(1,0,NQ_START(nr))=1
        CONECT(-1,1,CURRENT)=NQ_START(nr)
        CONECT(-1,0,CURRENT)=1
        CONECT(1,1,NQ_START(nr))=CURRENT
        NUM_NODES=1
        DO WHILE (CURRENT.NE.NQ_START(nr))
          ADJACENT=0
          DO i=-1,1,2
            DO j=1,NXQ(i,0,CURRENT,1)
              ADJACENT=ADJACENT+1
              POINTS(ADJACENT)=NXQ(i,j,CURRENT,1)
              IF (CALCULATED(NXQ(i,j,CURRENT,1))) THEN
                CONECT(-1,1,CURRENT)=POINTS(ADJACENT)
                LABLED=ADJACENT
              ENDIF
            ENDDO
          ENDDO
          CONECT(0,1,CURRENT)=CURRENT
          CONECT(-1,0,CURRENT)=1
          COUNT=0
          DO i=1,ADJACENT
            IF (i.NE.LABLED) THEN
              COUNT=COUNT+1
              CONECT(1,COUNT,CURRENT)=POINTS(I)
            ENDIF
          ENDDO
          IF(ADJACENT.EQ.2) CONECT(1,0,CURRENT)=1
          IF(ADJACENT.EQ.3) THEN !recalculates NXQ to remove
            CONECT(1,0,CURRENT)=2 ! conectivity at distale
            DO i=1,ADJACENT !grid points at a bifurcation
              IF(i.NE.LABLED) THEN
                nq1=POINTS(i)
                DO ii=1,ADJACENT
                  IF((ii.NE.i).AND.(ii.NE.LABLED)) THEN
                    nq2=POINTS(ii)
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
            i=-3
            FOUND=.FALSE.
            DO WHILE((.NOT.FOUND).AND.(I.LE.1))
              I=I+2
              COUNT=0
              DO WHILE((COUNT.LT.NXQ(i,0,nq1,1)).AND.(.NOT.FOUND))
                COUNT=COUNT+1
                IF(NXQ(i,count,nq1,1).EQ.nq2) THEN
                  FOUND=.TRUE.
                ENDIF
              ENDDO
            ENDDO
            IF(FOUND) THEN
              IF(COUNT.GE.2) THEN
                NXQ(i,0,nq1,1)=NXQ(i,0,nq1,1)-1
              ELSE IF(COUNT.EQ.1) THEN
                NXQ(i,1,nq1,1)=NXQ(i,NXQ(i,0,nq1,1),nq1,1)
                NXQ(i,0,nq1,1)=NXQ(i,0,nq1,1)-1
              ENDIF
            ENDIF
            i=-3
            FOUND=.FALSE.
            DO WHILE((.NOT.FOUND) .AND.(I.LT.1))
              I=I+2
              COUNT=0
              DO WHILE((COUNT.LT.NXQ(i,0,nq2,1)).AND.(.NOT.FOUND))
                COUNT=COUNT+1
                IF(NXQ(i,count,nq2,1).EQ.nq1) THEN
                  FOUND=.TRUE.
                ENDIF
              ENDDO
            ENDDO
            IF(FOUND) THEN
              IF(COUNT.GE.2) THEN
                NXQ(i,0,nq2,1)=NXQ(i,0,nq2,1)-1
              ELSE IF(COUNT.EQ.1) THEN
                NXQ(i,1,nq2,1)=NXQ(i,NXQ(i,0,nq2,1),nq2,1)
                NXQ(i,0,nq2,1)=NXQ(i,0,nq2,1)-1
              ENDIF
            ENDIF
          ENDIF
          CALCULATED(CURRENT)=.TRUE.
          NUM_NODES=NUM_NODES+1
          IF(ADJACENT.EQ.1) THEN
            FOUND=.FALSE.
            DO WHILE ((.NOT.FOUND).AND.(CURRENT.NE.NQ_START(nr)))
              CURRENT=CONECT(-1,1,CURRENT)
              IF(CONECT(1,0,CURRENT).GT.1) THEN
                IF (.NOT.CALCULATED(CONECT(1,2,CURRENT))) THEN
                  FOUND=.TRUE.
                  CURRENT=CONECT(1,2,CURRENT)
                ENDIF
              ENDIF
            ENDDO
          ELSE
            CURRENT=CONECT(1,1,CURRENT)
          ENDIF
        ENDDO
C calculating the grid points at bifurcations start and end points

        no_bc_points=0
        DO nq=NQR(1,nr),NQR(2,nr) !ipinit NPS 4/2/97
          ADJACENT=0
          DO ii=-1,1,2
            DO nqq=1,NXQ(ii,0,nq,1)
              ADJACENT=ADJACENT+1
            ENDDO
          ENDDO
          CALL ASSERT(ADJACENT.LE.3,'trifurcation',ERROR,*9999)
          IF((ADJACENT.NE.2)) THEN
            no_bc_points=no_bc_points+1
            BC_POINTS(no_bc_points)=nq
            IF (ADJACENT.NE.3) THEN
              BRANCH(no_bc_points)=.FALSE.
            ELSE
              BRANCH(no_bc_points)=.TRUE.
C This is the place to put the conectivity for the branch points
C need to record the colection of three branch points around
C a bifurcation NPS 17/11/96
            ENDIF
          ENDIF
        ENDDO

        NJTOT=NJ_LOC(NJL_GEOM,0,nr)

        CALL OPENF(IFILE,'DISK','outflow.dat',
     '    'UNKNOWN','SEQUEN','FORMATTED',160,ERROR,*9999)

        WRITE(IFILE,*) 'number of ends',NO_BC_POINTS
        DO I=1,NO_BC_POINTS
          nq=BC_POINTS(I)
          FLOW=(YQ(NYNQ(2,nq,0,nx),1,1,nx)**2.0d0)*PI*
     '      YQ(NYNQ(3,nq,0,nx),1,1,nx)
          IF(CONECT(1,0,nq).eq.0) THEN !end point
            TRANS_TIME_ART=0.0d0
            TRANS_TIME_VIEN=0.0d0
            CURRENT=nq
            DO WHILE(CURRENT.NE.NQ_START(nr))
              nq_next=CONECT(-1,1,CURRENT)
              DELTA_X=0.0d0
              DO nj=1,NJTOT !calculate delta x for half time step
                DELTA_X=DELTA_X+((XQ(nj,CURRENT)-
     '            XQ(nj,nq_next))**2.0d0)
              ENDDO
              DELTA_X=DMAX1((DELTA_X**0.5d0),0.60D0)
              ny_va=NYNQ(3,CURRENT,0,nx)
              ny_vv=NYNQ(6,CURRENT,0,nx)
              IF (ABS(YQ(ny_va,1,1,nx)).GT.LOOSE_TOL) THEN
                TRANS_TIME_ART=TRANS_TIME_ART+(DELTA_X/
     '            YQ(ny_va,1,1,nx))
              ELSE
                WRITE(*,*) 'ARTER PORB',CURRENT,ny_va
              ENDIF

              IF (ABS(YQ(ny_vV,1,1,nx)).GT.LOOSE_TOL) THEN
                TRANS_TIME_VIEN=TRANS_TIME_VIEN+(DELTA_X/
     '            YQ(ny_vv,1,1,nx))
              ELSE
                WRITE(*,*) 'vien PORB',CURRENT,ny_vv
              ENDIF

              CURRENT=nq_next
            ENDDO
            WRITE(IFILE,*) trans_time_ART,trans_time_VIEN,flow
          ENDIF
        ENDDO
        CALL CLOSEF(IFILE,ERROR,*9999)

C out put regional flows
        CALL OPENF(IFILE,'DISK','regflow.dat',
     '    'UNKNOWN','SEQUEN','FORMATTED',160,ERROR,*9999)

        DO I=1,NO_BC_POINTS
          nq=BC_POINTS(I)
          FLOW=(YQ(NYNQ(2,nq,0,nx),1,1,nx)**2.0d0)*PI*
     '      YQ(NYNQ(3,nq,0,nx),1,1,nx)
          IF(CONECT(1,0,nq).eq.0) THEN !end point
            CALL ASSERT(NENQ(0,nq).EQ.1,'>>Not an end point',
     '        ERROR,*9999)
            ne=NENQ(1,nq)
            IF(NQNE(ne,1).eq.nq) THEN
              np=NPNE(1,NBJ(1,ne),ne)
            ELSE IF(NQNE(ne,NQET(NQS(ne))).eq.nq) THEN
              np=NPNE(NNT(NBJ(1,ne)),NBJ(1,ne),ne)
            ELSE
              WRITE(*,*) 'error in evcoro'
            ENDIF
            WRITE(IFILE,*) NEP(np),XIP(1,np),XIP(2,np),XIP(3,np),flow
          ENDIF
        ENDDO
        CALL CLOSEF(IFILE,ERROR,*9999)
      ENDIF
      CALL EXITS('EVCORO')
      RETURN
 9999 CALL ERRORS('EVCORO',ERROR)
      CALL EXITS('EVCORO')
      RETURN 1
      END


