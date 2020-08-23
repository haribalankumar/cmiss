      SUBROUTINE MARCH4(BC_POINTS,BRANCH,CALCULATED,CONECT,CQ,ER,ES,
     '  GKK,GRR,ISC_GKK,ISR_GKK,LGE,NBJ,NEELEM,NENQ,NHQ,NPNE,NQNE,NQS,
     '  nr,nx,NXQ,NYNQ,NYQNR,XP,XIP,XQ,YQ,TIME_VALUES,
     '  NTIME_POINTS,NTIME_NR,ERROR,*)

C#### Subroutine: MARCH4
C###  Description:
C###    MARCH4 performs time integration of linear or nonlinear
C###    equations using explicit or implicit finite differences.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER BC_POINTS(3,3,0:NQM),ISC_GKK(NISC_GKKM),
     '  ISR_GKK(NISR_GKKM),LGE(NHM*NSM,NRCM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NHQ(NRM),
     '  NPNE(NNM,NBFM,NEM),NQNE(NEQM,NQEM),NQS(NEQM),nr,nx,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNQ(NHM,NQM,0:NRCM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM),
     '  NTIME_POINTS(NTIMEVARSM),NTIME_NR(0:NTIMEVARSM,NRM)
      REAL*8 CQ(NMM,NQM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),GKK(NZ_GKK_M),
     '  GRR(NOM),TIME,XQ(NJM,NQM),XP(NKM,NVM,NJM,NPM),XIP(NIM,NPM),
     '  YQ(NYQM,NIQM),TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      CHARACTER ERROR*(*)
      LOGICAL CALCULATED(NQM)

!     Local Variables
      INTEGER ADJACENT,BIFUR_COUNT,CONECT(-1:1,0:2,NQM),COUNT,CURRENT,
     '  I,ii,j,jj,LABLED,DATA_COUNT,FILE_COUNT,OLD_CURENT,nb1,nb2,ne1,
     '  ne2,ne3,ne_current,nj,nj1,nj2,njj,np1,np2,nq,nqq,nq1,nq1_pre,
     '  nq2,nq2_pre,nq3,nq_old,no_bc_points,no_ne,no_nq,
     '  no_nynr,NUM_NODES,ny,ny_p,ny_r,ny_v,POINTS(3),TIME_STEPS,
     '  PATH_ARRAY(60),IBEG,IEND,n_term,N_TERM_Q(-1:1,0:NEM),np
      REAL*8 FA,FB,GRAD_DIF,LAMBDA,LAMBDA1,LEN_COUNT,STEP,VEL,TRACE1,
     '  TRACE2,TRACE3,XI_DIST1,XI_DIST2, XI_SUM1,XI_SUM2,
     '  xc_nq,yc_nq,zc_nq,xc_np,yc_np,zc_np
      REAL TIME_START2(1),TIME_STOP(1),TOT_BITIME,TOT_BOTIME,TOT_GRTIME
      CHARACTER FMT*100
      INTEGER*4 WORK_PTR
      LOGICAL BRANCH(NQM),DISCONT,DONE,FOUND,UPDATE_MATRIX,
     '  UPDATE_VECTOR,HALF_TIME_STEP,PRINT_FILE

      CALL ENTERS('MARCH4',*9999)

      FMT='(4E18.10)'
      IF(ITYP16(nr,nx).EQ.4) THEN !Lax-Wendroff
        IF(ITYP3(nr,nx).EQ.1) THEN !flow in elastic tube
          DO nq=1,NQT
            CALCULATED(nq)=.FALSE.
          ENDDO !nq
          CALCULATED(NQ_START(nr))=.TRUE.!ipinit NPS 4/2/97
          CURRENT=NXQ(1,1,NQ_START(nr),1)
          IF(CURRENT.EQ.0) CURRENT=NXQ(-1,1,NQ_START(nr),1)
          CONECT(1,0,NQ_START(nr))=1
          CONECT(-1,1,CURRENT)=NQ_START(nr)
          CONECT(-1,0,CURRENT)=1
          CONECT(1,1,NQ_START(nr))=CURRENT
          OLD_CURENT=NQ_START(nr)
          NUM_NODES=1
          DO WHILE (CURRENT.NE.NQ_START(nr))
            ADJACENT=0
            DO i=-1,1,2
              DO j=1,NXQ(i,0,CURRENT,1)
                ADJACENT=ADJACENT+1
                POINTS(ADJACENT)=NXQ(i,j,CURRENT,1)
                IF(CALCULATED(NXQ(i,j,CURRENT,1))) THEN
 !                  CONECT(-1,1,CURRENT)=POINTS(ADJACENT)
                  LABLED=ADJACENT
                ENDIF !CALCULATED
              ENDDO !j
            ENDDO !i
            CONECT(-1,1,CURRENT)=OLD_CURENT !NEW 28/5/99
            CONECT(0,1,CURRENT)=CURRENT
            CONECT(-1,0,CURRENT)=1
            COUNT=0
            DO i=1,ADJACENT
              IF(i.NE.LABLED) THEN
                IF(.NOT.CALCULATED(POINTS(I))) THEN
                  COUNT=COUNT+1
                  CONECT(1,COUNT,CURRENT)=POINTS(I)
                ENDIF !.not.CALCULATED
              ENDIF !i
            ENDDO !i
            IF(ADJACENT.EQ.2) THEN
              IF(.NOT.(CALCULATED(POINTS(1))
     '          .AND.CALCULATED(POINTS(2)))) THEN !a normal position
                CONECT(1,0,CURRENT)=1
              ENDIF
            ENDIF !ADJACENT=2
            IF(ADJACENT.EQ.3) THEN  !recalculates NXQ to remove
              CONECT(1,0,CURRENT)=2  !conectivity at diastole
              DO i=1,ADJACENT        !grid points at a bifurcation
                IF(i.NE.LABLED) THEN
                  nq1=POINTS(i)
                  DO ii=1,ADJACENT
                    IF((ii.NE.i).AND.(ii.NE.LABLED)) THEN
                      nq2=POINTS(ii)
                    ENDIF !ii.ne.i
                  ENDDO !ii
                ENDIF !i
              ENDDO !i
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
                ENDDO !while
              ENDDO !while
              IF(FOUND) THEN
                IF(COUNT.GE.2) THEN
                  NXQ(i,0,nq1,1)=NXQ(i,0,nq1,1)-1
                ELSE IF(COUNT.EQ.1) THEN
                  NXQ(i,1,nq1,1)=NXQ(i,NXQ(i,0,nq1,1),nq1,1)
                  NXQ(i,0,nq1,1)=NXQ(i,0,nq1,1)-1
                ENDIF !COUNT
              ENDIF !FOUND
              i=-3
              FOUND=.FALSE.
              DO WHILE((.NOT.FOUND) .AND.(I.LT.1))
                I=I+2
                COUNT=0
                DO WHILE((COUNT.LT.NXQ(i,0,nq2,1)).AND.(.NOT.FOUND))
                  COUNT=COUNT+1
                  IF(NXQ(i,count,nq2,1).EQ.nq1) THEN
                    FOUND=.TRUE.
                  ENDIF !NXQ
                ENDDO !WHILE
              ENDDO !WHILE
              IF(FOUND) THEN
                IF(COUNT.GE.2) THEN
                  NXQ(i,0,nq2,1)=NXQ(i,0,nq2,1)-1
                ELSE IF(COUNT.EQ.1) THEN
                  NXQ(i,1,nq2,1)=NXQ(i,NXQ(i,0,nq2,1),nq2,1)
                  NXQ(i,0,nq2,1)=NXQ(i,0,nq2,1)-1
                ENDIF !COUNT
              ENDIF !FOUND
            ENDIF !ADJACENT=3
            CALCULATED(CURRENT)=.TRUE.
            NUM_NODES=NUM_NODES+1
            IF((ADJACENT.EQ.2).AND.CALCULATED(POINTS(1))
     '        .AND.CALCULATED(POINTS(2))) THEN
              FOUND=.FALSE.
              DO WHILE ((.NOT.FOUND).AND.(CURRENT.NE.NQ_START(nr)))
                CURRENT=CONECT(-1,1,CURRENT)
                IF(CONECT(1,0,CURRENT).GT.1) THEN
                  IF(.NOT.CALCULATED(CONECT(1,2,CURRENT))) THEN
                    FOUND=.TRUE.
                    OLD_CURENT=CURRENT !NEW 28/5/99
                    CURRENT=CONECT(1,2,CURRENT)
                  ENDIF
                ENDIF
              ENDDO !WHILE
            ELSE IF(ADJACENT.EQ.1) THEN
              FOUND=.FALSE.
              DO WHILE ((.NOT.FOUND).AND.(CURRENT.NE.NQ_START(nr)))
                CURRENT=CONECT(-1,1,CURRENT)
                IF(CONECT(1,0,CURRENT).GT.1) THEN
                  IF (.NOT.CALCULATED(CONECT(1,2,CURRENT))) THEN
                    FOUND=.TRUE.
                    OLD_CURENT=CURRENT !NEW 28/5/99
                    CURRENT=CONECT(1,2,CURRENT)
                  ENDIF
                ENDIF
              ENDDO !WHILE
            ELSE
              OLD_CURENT=CURRENT           !NEW 28/5/99
              CURRENT=CONECT(1,1,CURRENT)
            ENDIF !ADJACENT=1
          ENDDO !ADJACENT=2
C***      Calculating the grid points at bifurcations and the start and
C***      end points where boundary condition need to be applied
          no_bc_points=0
          DO nq=NQR(1,nr),NQR(2,nr) !New BC_POINTS array NPS 29/5/99
            ADJACENT=0
            DO ii=-1,1,2
              DO nqq=1,NXQ(ii,0,nq,1)
                ADJACENT=ADJACENT+1
              ENDDO !nqq
            ENDDO !ii
            CALL ASSERT(ADJACENT.LE.3,'trifurcation',ERROR,*9999)
            IF((ADJACENT.NE.2)) THEN !normal pt or terminal pt
              no_bc_points=no_bc_points+1
              BC_POINTS(1,1,no_bc_points)=nq
              IF(ADJACENT.NE.3) THEN !not bifurcation
                BRANCH(no_bc_points)=.FALSE.
              ELSE !is bifurcation
                BRANCH(no_bc_points)=.TRUE.
                BC_POINTS(2,1,no_bc_points)=CONECT(1,1,nq)
                BC_POINTS(3,1,no_bc_points)=CONECT(1,2,nq)
C Determine the a1,b1,c1 for the bifurcation
                DO jj=1,3 !segments
                  nq1=NXQ(1,1,BC_POINTS(jj,1,no_bc_points),1)
                  nq2=NXQ(-1,1,BC_POINTS(jj,1,no_bc_points),1)
                  IF (nq1.NE.0) THEN
                    ne1=NENQ(1,nq1)
                  ELSE
                    ne1=0
                  ENDIF
                  IF (nq2.NE.0) THEN
                    ne2=NENQ(1,nq2)
                  ELSE
                    ne2=0
                  ENDIF
                  ne_current=NENQ(1,BC_POINTS(jj,1,no_bc_points))
                  IF(ne1.EQ.ne_current) THEN
                    BC_POINTS(jj,2,no_bc_points)=
     '                NXQ(1,1,BC_POINTS(jj,1,no_bc_points),1)
                  ELSE
                    BC_POINTS(jj,2,no_bc_points)=
     '                NXQ(-1,1,BC_POINTS(jj,1,no_bc_points),1)
                  ENDIF
                ENDDO !jj
C***            This determines the a2,b2,c2 for the bifurcation
                DO jj=1,3 !segments
                  nq1=CONECT(-1,1,BC_POINTS(jj,2,no_bc_points))
                  nq2=BC_POINTS(jj,1,no_bc_points)
                  IF (nq1.eq.nq2) THEN
                    BC_POINTS(jj,3,no_bc_points)=nq2
                  ELSE
                    BC_POINTS(jj,3,no_bc_points)=
     '                BC_POINTS(jj,2,no_bc_points)
                  ENDIF
                ENDDO !jj
C***            This determines the half space step grid points
C***            a1-a2,b1-b2,c1-c2 for the bifurcation

C***            This is the place to put the connectivity for the
C***            branch points.  Need to record the colection of three
C***            branch points around a bifurcation NPS 17/11/96
              ENDIF !ADJACENT.NE.3
            ENDIF !ADJACENT.NE.2
          ENDDO !nq
          BC_POINTS(1,1,0)=no_bc_points
        ENDIF !ITYP3(nr,nx)=1
      ENDIF !ITYP16(nr,nx)=4
C***  This section of code calculates the conectivity of a 1D branching
C***  network in a 3d host mesh, NXQ can not be used directly because
C***  adjacent elemnts may not have consistant Xi directions.

      DO no_nq=1,BC_POINTS(1,1,0) !now loops have been added this should
        nq=BC_POINTS(1,1,no_nq) ! be rechecked NPS 29/5/99
        IF(BRANCH(no_nq)) THEN
          DO i=1,2
            DISCONT=.FALSE.
            nq1=conect(1,I,nq)
            trace1=YQ(NYNQ(1,nq1,0),2)
            ne1=INT(CQ(8,nq1))
            nq2=conect(1,1,nq1)
            IF(nq2.NE.0) THEN
              trace2=YQ(NYNQ(1,nq2,0),2)
              ne2=INT(cq(8,nq2))
              nq3=CONECT(1,1,nq2)
              trace3=YQ(NYNQ(1,nq3,0),2)
              ne3=INT(cq(8,nq3))
              DO WHILE((CONECT(1,0,nq3).EQ.1).AND.(.NOT.DISCONT))
                IF((ne2.NE.ne1).OR.(ne2.NE.ne3).OR.DISCONT) THEN
                  GRAD_DIF=TRACE1-2.0d0*TRACE2+TRACE3
                  IF((DABS(GRAD_DIF).GT.0.5d0).OR.DISCONT) THEN
                    DISCONT=.TRUE.
                  ENDIF
                ENDIF !ne2
C!!! Something uninitialized here
                IF (((DABS(TRACE1-TRACE2)).GT.2.50D0).OR.
     '            ((DABS(TRACE2-TRACE3)).GT.2.50D0)) THEN
                   DISCONT=.TRUE.
                ENDIF
C the above code removes discontinuities in the trace do to linear
C interpolation of the host mesh. Once a full C1 conintous host
C mechanics mesh is implimented it can be removed. NPS 12/11/98
                nq1=nq2
                trace1=trace2
                ne1=ne2
                nq2=nq3
                trace2=trace3
                ne2=ne3
                nq3=conect(1,1,nq2)
                trace3=YQ(NYNQ(1,nq3,0),2)
                ne3=INT(cq(8,nq3))
              ENDDO !while
              IF(DISCONT) THEN
                nq1=conect(1,I,nq)
                TRACE1=YQ(NYNQ(1,nq1,0),2)
                DO WHILE(CONECT(1,0,nq1).EQ.1)
                  nq1=conect(1,1,nq1)
                  YQ(NYNQ(1,nq1,0),2)=TRACE1
                  YQ(NYNQ(4,nq1,0),2)=TRACE1
                ENDDO !WHILE
              ENDIF ! DISCONT
            ENDIF ! (nq2.NE.0)
          ENDDO !i
        ENDIF !BRANCH(no_nq)
      ENDDO !no_nq

      DO no_nq=1,BC_POINTS(1,1,0) ! This code removes dicontinuities
        nq=BC_POINTS(1,1,no_nq)   ! in the lambda vessel stretch values
        IF(BRANCH(no_nq)) THEN ! between vessel segments NPS 26/11/99
          DO I=1,2             ! checking both branches
            nq1=conect(1,I,nq)
            LAMBDA1=YQ(NYNQ(3,nq1,0),2)
            DO WHILE(CONECT(1,0,nq1).EQ.1) !while not end of element
              nq1=conect(1,1,nq1)
              YQ(NYNQ(3,nq1,0),2)=LAMBDA1
              YQ(NYNQ(6,nq1,0),2)=LAMBDA1
            ENDDO
          ENDDO !I=1,2
        ENDIF !BRANCH(no_nq)
      ENDDO !no_nq

      CALL CLOSEF(IFILE,ERROR,*9999)
      TIME=T_RESTART(nx)
      WORK_PTR=0
      PRINT_FILE=.TRUE.
      TOT_BITIME=0.0
      TOT_BOTIME=0.0
      TOT_GRTIME=0.0
      TIME_STEPS=INT((TFINISH-TSTART+1.0d-6)/DT)
      DO I=1,50
        PATH_ARRAY(I)=1
      ENDDO
      UPDATE_VECTOR=.TRUE.
      IF(ITYP16(nr,nx).NE.4) THEN ! not an explicit scheme
        CALL ALLOCATE_MEMORY(NOQT(1,1,nr,nx)*NOQT(2,1,nr,nx),0,CHARTYPE,
     '    WORK_PTR,MEM_INIT,ERROR,*9999)
      ELSE
        CALL ALLOCATE_MEMORY(1,0,CHARTYPE,
     '    WORK_PTR,MEM_INIT,ERROR,*9999)
      ENDIF
      UPDATE_MATRIX=.TRUE.
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C PM 26-JUL-01 : to define boundary conditions from a file
      IF(KTYP3_INIT(nx).EQ.3) THEN
        CALL STRING_TRIM(FILE03,IBEG,IEND)
        CALL OPENF(IOFILE1,'TERM',FILE03(IBEG:IEND),'UNKNOWN',
     '    'SEQUEN','UNFORMATTED',132,ERROR,*9999)
      ENDIF

C PM 29-NOV-01 : Determine grid points where transition occurs from arteries
C                to veins from the corrosponding nodes
      IF((VENOUS_NETWORK.EQ.'Y').OR.(VENOUS_NETWORK.EQ.'y')) THEN
        IF(N_VENOUS_GEOM.EQ.2) THEN

          DO nq=NQR(1,nr),NQR(2,nr)
            xc_nq=XQ(1,nq)
            yc_nq=XQ(2,nq)
            zc_nq=XQ(3,nq)

            DO n_term=1,N_TERM_P(0)
              np=N_TERM_P(n_term)
              xc_np=XP(1,1,1,np)
              yc_np=XP(1,1,2,np)
              zc_np=XP(1,1,3,np)

              IF((ABS(xc_np-xc_nq).LE.1.0D-07).AND.
     '          (ABS(yc_np-yc_nq).LE.1.0D-07).AND.
     '          (ABS(zc_np-zc_nq).LE.1.0D-07)) THEN
                N_TERM_Q(0,n_term)=nq
              ENDIF
            ENDDO
          ENDDO
          N_TERM_Q(0,0)=N_TERM_P(0)

        ENDIF
      ENDIF

      DO I=1,TIME_STEPS
        IF(ITYP16(nr,nx).eq.4) THEN
          HALF_TIME_STEP=.TRUE.
          UPDATE_VECTOR=.TRUE.
          TIME=TIME+(DT/2.0d0) !half time step

           CALL ASSEMBLE9(BC_POINTS,BRANCH,CQ,CONECT,ER,ES,GKK,GRR,
     '      HALF_TIME_STEP,ISC_GKK,ISR_GKK,LGE,1,1,nr,nx,NHQ,
     '      NQ_START(nr),N_TERM_Q,NXQ,NYNQ,NYQNR,UPDATE_MATRIX,
     '      UPDATE_VECTOR,%VAL(WORK_PTR),TIME,TOT_BITIME,TOT_BOTIME,
     '      TOT_GRTIME,XQ,YQ,TIME_VALUES,NTIME_POINTS,NTIME_NR,
     &      N_VENOUS_GEOM,VENOUS_NETWORK,ERROR,*9999)

C no call to solve is needed for an explicite scheme, in this
C case the right hand side is written directly into YQ
          HALF_TIME_STEP=.FALSE.
          TIME=TIME+(DT/2.0d0)
        ENDIF
        DO no_nynr=1,NYQNR(0,0,1,nr)
          ny=NYQNR(no_nynr,0,1,nr)
          YQ(ny,8)=YQ(ny,1) !temporary storage for previous time step
        ENDDO
        UPDATE_VECTOR=.TRUE.
          CALL ASSEMBLE9(BC_POINTS,BRANCH,CQ,CONECT,ER,ES,GKK,GRR,
     '      HALF_TIME_STEP,ISC_GKK,ISR_GKK,LGE,1,1,nr,nx,NHQ,
     '      NQ_START(nr),N_TERM_Q,NXQ,NYNQ,NYQNR,UPDATE_MATRIX,
     '      UPDATE_VECTOR,%VAL(WORK_PTR),TIME,TOT_BITIME,TOT_BOTIME,
     '    TOT_GRTIME,XQ,YQ,TIME_VALUES,NTIME_POINTS,NTIME_NR,
     &    N_VENOUS_GEOM,VENOUS_NETWORK,ERROR,*9999)
C no call to solve is needed for an explicite scheme, in this
C case the right hand side is written directly into YQ
        HALF_TIME_STEP=.FALSE.
        IF(ITYP16(nr,nx).eq.4) THEN !Lax-Wendroff
          IF(ITYP3(nr,nx).eq.1) THEN !flow in elastic tube
C PM 03-OCT-01
C            DATA_COUNT=DATA_COUNT+1
C            IF(DATA_COUNT.GE.((T1-T0)/(20*TINCR))) THEN
C              DATA_COUNT=0
C              IF (DATA_FILE_COUNT.EQ.0) THEN
C                FILE_COUNT=FILE_COUNT+1
C                PRINT_FILE=.TRUE.
C              ELSE IF(PRINT_FILE) THEN !only one file per time step
C              ENDIF
C              IF(PRINT_FILE) THEN
C                PRINT_FILE=.FALSE.
                nq=NQ_START(nr)
                LAMBDA=(YQ(NYNQ(3,nq,0),9)*(T1-TIME)/(T1-T0))
     '                +(((TIME-T0)/(T1-T0))*
     '                YQ(NYNQ(3,nq,0),2))
                COUNT=1
                LEN_COUNT=0.0d0
                nj=NJ_LOC(NJL_FIEL,1,nr)
                DONE=.FALSE.
                BIFUR_COUNT=1
                DO while(.NOT.DONE)
                  ny=NYNQ(3,nq,0)
C                  WRITE(IOFILE1,FMT) LEN_COUNT,YQ(ny,1),
C     '              (-1.0d0*YQ((ny+3),1)),LAMBDA
                  IF(conect(1,1,nq).EQ.0) THEN
                    DONE=.TRUE.
                  ELSE
                    nq_old=nq
                    COUNT=COUNT+1
                    IF (CONECT(1,0,nq).GT.1)THEN
                      COUNT=COUNT-1
                      ne1=NENQ(1,CONECT(1,1,nq))
                      ne2=NENQ(1,CONECT(1,2,nq))
                      nb1=NBJ(1,ne1)
                      nb2=NBJ(1,ne2)
                      IF(XP(1,1,nj,NPNE(1,nb1,ne1)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb1),nb1,ne1))) THEN
                        np1=NPNE(NNT(nb1),nb1,ne1)
                      ELSE
                        np1=NPNE(1,nb1,ne1)
                      ENDIF
                      IF(XP(1,1,nj,NPNE(1,nb2,ne2)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb2),nb2,ne2))) THEN
                        np2=NPNE(NNT(nb2),nb2,ne2)
                      ELSE
                        np2=NPNE(1,nb2,ne2)
                      ENDIF
                      IF(PATH_ARRAY(BIFUR_COUNT).EQ.1) THEN
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ENDIF
                      ELSE
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ENDIF
                      ENDIF
                      BIFUR_COUNT=BIFUR_COUNT+1
                    ELSE
                      nq=CONECT(1,1,nq)
                    ENDIF
                    STEP=0.0d0
                    LAMBDA=(YQ(NYNQ(3,nq,0),9)*
     '                (TFINISH-TIME)/(TFINISH-TSTART))
     '                +(((TIME-TSTART)/(TFINISH-TSTART))*
     '                YQ(NYNQ(3,nq,0),2))
                    DO njj=1,3
                      STEP=STEP+((XQ(njj,nq)-XQ(njj,nq_old))**2.0d0)
                    ENDDO
                    STEP=STEP*LAMBDA
                    LEN_COUNT=LEN_COUNT+(STEP**0.5d0)
                  ENDIF
                ENDDO
                nq=NQ_START(nr)
                COUNT=1
                LEN_COUNT=0.0d0
                nj=NJ_LOC(NJL_FIEL,1,nr)
                DONE=.FALSE.
                BIFUR_COUNT=1
                DO while(.NOT.DONE)
                  ny=NYNQ(1,nq,0)
C                  WRITE(IOFILE2,FMT) LEN_COUNT,YQ(ny,1),YQ(ny,2),
C     '              YQ(ny+3,1)
                  FA=PI*(YQ(NYNQ(2,nq,0),1)**2.0d0)*
     '              YQ(NYNQ(3,nq,0),1)
                  FB=-PI*(YQ(NYNQ(5,nq,0),1)**2.0d0)*
     '              YQ(NYNQ(6,nq,0),1)
C                  WRITE(IOFILE3,FMT) LEN_COUNT,FA,FB,0.0D0
                  IF(conect(1,1,nq).EQ.0) THEN
                    FA=PI*(YQ(NYNQ(2,nq,0),1)**2.0d0)*
     '                YQ(NYNQ(3,nq,0),1)
                    FB=-PI*(YQ(NYNQ(5,nq,0),1)**2.0d0)*
     '                YQ(NYNQ(6,nq,0),1)
                    DONE=.TRUE.
                  ELSE
                    nq_old=nq
                    COUNT=COUNT+1
                    IF (CONECT(1,0,nq).GT.1)THEN
                      COUNT=COUNT-1
                      ne1=NENQ(1,CONECT(1,1,nq))
                      ne2=NENQ(1,CONECT(1,2,nq))
                      nb1=NBJ(1,ne1)
                      nb2=NBJ(1,ne2)
                      IF(XP(1,1,nj,NPNE(1,nb1,ne1)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb1),nb1,ne1))) THEN
                        np1=NPNE(NNT(nb1),nb1,ne1)
                      ELSE
                        np1=NPNE(1,nb1,ne1)
                      ENDIF
                      IF(XP(1,1,nj,NPNE(1,nb2,ne2)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb2),nb2,ne2))) THEN
                        np2=NPNE(NNT(nb2),nb2,ne2)
                      ELSE
                        np2=NPNE(1,nb2,ne2)
                      ENDIF
                      IF(PATH_ARRAY(BIFUR_COUNT).EQ.1) THEN
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ENDIF
                      ELSE
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ENDIF
                      ENDIF
                      BIFUR_COUNT=BIFUR_COUNT+1
                    ELSE
                      nq=CONECT(1,1,nq)
                    ENDIF
                    STEP=0.0d0
                    LAMBDA=(YQ(NYNQ(3,nq,0),9)*
     '                (TFINISH-TIME)/(TFINISH-TSTART))
     '                +(((TIME-TSTART)/(TFINISH-TSTART))*
     '                YQ(NYNQ(3,nq,0),2))
                    DO njj=1,3
                      STEP=STEP+((XQ(njj,nq)-XQ(njj,nq_old))**2.0d0)
                    ENDDO
                    STEP=STEP*LAMBDA
                    LEN_COUNT=LEN_COUNT+(STEP**0.5d0)
                  ENDIF
                ENDDO
                nq=NQ_START(nr)
                COUNT=1
                LEN_COUNT=0.0d0
                nj=NJ_LOC(NJL_FIEL,1,nr)
                DONE=.FALSE.
                BIFUR_COUNT=1
                DO while(.NOT.DONE)
                  ny=NYNQ(2,nq,0)
                  YQ(ny,6)=1.0d0 !labels grid points along a path
                  IF(conect(1,1,nq).EQ.0) THEN
                    DONE=.TRUE.
                  ELSE
                    nq_old=nq
                    COUNT=COUNT+1
                    IF (CONECT(1,0,nq).GT.1)THEN
                      COUNT=COUNT-1
                      ne1=NENQ(1,CONECT(1,1,nq))
                      ne2=NENQ(1,CONECT(1,2,nq))
                      nb1=NBJ(1,ne1)
                      nb2=NBJ(1,ne2)
                      IF(XP(1,1,nj,NPNE(1,nb1,ne1)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb1),nb1,ne1))) THEN
                        np1=NPNE(NNT(nb1),nb1,ne1)
                      ELSE
                        np1=NPNE(1,nb1,ne1)
                      ENDIF
                      IF(XP(1,1,nj,NPNE(1,nb2,ne2)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb2),nb2,ne2))) THEN
                        np2=NPNE(NNT(nb2),nb2,ne2)
                      ELSE
                        np2=NPNE(1,nb2,ne2)
                      ENDIF
                      IF(PATH_ARRAY(BIFUR_COUNT).EQ.1) THEN
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ENDIF
                      ELSE
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ENDIF
                      ENDIF
                      BIFUR_COUNT=BIFUR_COUNT+1
                    ELSE
                      nq=CONECT(1,1,nq)
                    ENDIF
                    STEP=0.0d0
                    LAMBDA=(YQ(NYNQ(3,nq,0),9)*(TFINISH-TIME)/
     '                (TFINISH-TSTART))
     '                +(((TIME-TSTART)/(TFINISH-TSTART))*
     '                YQ(NYNQ(3,nq,0),2))
                    DO njj=1,3
                      STEP=STEP+((XQ(njj,nq)-XQ(njj,nq_old))**2.0d0)
                    ENDDO
                    STEP=STEP*LAMBDA
                    LEN_COUNT=LEN_COUNT+(STEP**0.5d0)
                  ENDIF
                ENDDO
                DO nq=NQR(1,nr),NQR(2,nr)
                  ny_p=NYNQ(1,nq,0)
                  ny_r=NYNQ(2,nq,0)
                  ny_v=NYNQ(3,nq,0)
                  VEL=YQ(ny_v,1)
                  IF(DABS(VEL).LT.LOOSE_TOL) THEN
                    VEL=LOOSE_TOL
                  ENDIF
                ENDDO
                DO nq=NQR(1,nr),NQR(2,nr)
                  ny_p=NYNQ(4,nq,0)
                  ny_r=NYNQ(5,nq,0)
                  ny_v=NYNQ(6,nq,0)
                  VEL=YQ(ny_v,1)
                  IF(DABS(VEL).LT.LOOSE_TOL) THEN
                    VEL=LOOSE_TOL
                  ENDIF
                  LAMBDA=YQ(NYNQ(3,nq,0),2)
                ENDDO
C              ENDIF
C            ENDIF
          ENDIF
        ENDIF
      ENDDO !time
      T_RESTART(nx)=TIME !PM 26-JUL-01
      WRITE(OP_STRING,'(''TIME'',F12.8,'' seconds'')') TIME
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      CALL CPU_TIMER(CPU_USER,TIME_STOP)

C PM 03-OCT-01 : added condition for coupled flow/solid mech
C PM 29-NOV-01 : This is not required as grid values are now exported
C                But don't remove this section as it may be required in the
C                future.
C      IF (KTYP90.EQ.9) THEN  ! Coupled flow and solid mech problem

C        IF(ITYP16(nr,nx).eq.4) THEN   ! Lax-Wendroff
C          IF(ITYP3(nr,nx).eq.1) THEN  ! flow in elastic tube
C            nj=NJ_LOC(NJL_FIEL,2,nr)  ! pressure
C            nj1=NJ_LOC(NJL_FIEL,3,nr) ! velocity
C            nj2=NJ_LOC(NJL_FIEL,1,nr) ! radius
C            DO no_ne=1,NEELEM(0,nr)
C              ne1=NEELEM(no_ne,nr)
C              nb1=NBJ(1,ne1)
C              np1=NPNE(1,nb1,ne1)
C              np2=NPNE(NNT(nb1),nb1,ne1)
C              nq1=NQNE(ne1,1)
C              nq1_pre=CONECT(-1,1,NQ1)
C              nq2=NQNE(ne1,NQET(NQS(ne1)))
C              nq2_pre=CONECT(-1,1,NQ2)
C              XI_DIST1=(DABS(XIP(1,NP1)-YQ(NYNQ(1,nq1,0),10)))
C     '          +(DABS(XIP(2,NP1)-YQ(NYNQ(2,nq1,0),10)))
C     '          +(DABS(XIP(3,NP1)-YQ(NYNQ(3,nq1,0),10)))
C              XI_DIST2=(DABS(XIP(1,NP2)-YQ(NYNQ(1,nq2,0),10)))
C     '          +(DABS(XIP(2,NP2)-YQ(NYNQ(2,nq2,0),10)))
C     '          +(DABS(XIP(3,NP2)-YQ(NYNQ(3,nq2,0),10)))
C              XI_SUM1=DABS(YQ(NYNQ(1,nq1,0),10))
C     '          +DABS(YQ(NYNQ(2,nq1,0),10))
C     '          +DABS(YQ(NYNQ(3,nq1,0),10))
C              XI_SUM2=DABS(YQ(NYNQ(1,nq2,0),10))
C     '          +DABS(YQ(NYNQ(2,nq2,0),10))
C     '          +DABS(YQ(NYNQ(3,nq2,0),10))
C could use radius for this loop as well
C              IF((XI_SUM2+XI_SUM1).GT.LOOSE_TOL) THEN
C                IF((XI_DIST1+XI_DIST2).GT.LOOSE_TOL) THEN
C                  np2=NPNE(1,nb1,ne1)
C                  np1=NPNE(NNT(nb1),nb1,ne1)
C                ENDIF
C              ENDIF
C              IF(NQ1_PRE.EQ.0) THEN
C                ny_p=NYNQ(1,nq1,0)
C                ny_r=NYNQ(2,nq1,0)
C                ny_v=NYNQ(3,nq1,0)
C                IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
C                  XP(1,1,nj,np1)=YQ(ny_p,1)
C                ELSE
C                  XP(1,1,nj,np1)=LOOSE_TOL
C                ENDIF
C                XP(1,1,nj2,np1)=YQ(ny_r,3) !radius left unchanged
C                IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
C                  XP(1,1,nj1,np1)=YQ(ny_v,1)
C                ELSE
C                  XP(1,1,nj1,np1)=LOOSE_TOL
C                ENDIF
C              ELSE
C                IF(CONECT(1,0,NQ1_pre).LT.2) THEN
C not a distal bifuraction
C                  ny_p=NYNQ(1,nq1,0)
C                  ny_r=NYNQ(2,nq1,0)
C                  ny_v=NYNQ(3,nq1,0)
C                  IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
C                    XP(1,1,nj,np1)=YQ(ny_p,1)
C                  ELSE
C                    XP(1,1,nj,np1)=LOOSE_TOL
C                  ENDIF
C                  XP(1,1,nj2,np1)=YQ(ny_r,3) !radius left unchanged
C                  IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
C                    XP(1,1,nj1,np1)=YQ(ny_v,1)
C                  ELSE
C                    XP(1,1,nj1,np1)=LOOSE_TOL
C                  ENDIF
C                ENDIF
C              ENDIF
C              IF(NQ2_PRE.EQ.0) THEN
C                ny_p=NYNQ(1,nq2,0)
C                ny_r=NYNQ(2,nq2,0)
C                ny_v=NYNQ(3,nq2,0)
C                IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
C                  XP(1,1,nj,np2)=YQ(ny_p,1)
C                ELSE
C                  XP(1,1,nj,np2)=LOOSE_TOL
C                ENDIF
C                XP(1,1,nj2,np2)=YQ(ny_r,3) !radius left unchanged
C                IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
C                  XP(1,1,nj1,np2)=YQ(ny_v,1)
C                ELSE
C                  XP(1,1,nj1,np2)=LOOSE_TOL
C                ENDIF
C              ELSE
C                IF(CONECT(1,0,NQ2_pre).LT.2) THEN
C not a distal bifuraction
C                  ny_p=NYNQ(1,nq2,0)
C                  ny_r=NYNQ(2,nq2,0)
C                  ny_v=NYNQ(3,nq2,0)
C                  IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
C                    XP(1,1,nj,np2)=YQ(ny_p,1)
C                  ELSE
C                    XP(1,1,nj,np2)=LOOSE_TOL
C                  ENDIF
C                  XP(1,1,nj2,np2)=YQ(ny_r,3)
C                  IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
C                    XP(1,1,nj1,np2)=YQ(ny_v,1)
C                  ELSE
C                    XP(1,1,nj1,np2)=LOOSE_TOL
C                  ENDIF
C                ENDIF
C              ENDIF
C            ENDDO
C          ENDIF
C        ENDIF

C        CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)

C        IF(DATA_FILE_COUNT.NE.0) THEN
C          DATA_FILE_COUNT=DATA_FILE_COUNT+1
C        ENDIF

C      ENDIF  ! Coupled flow and solid mech problem
      CALL EXITS('MARCH4')
      RETURN
 9999 CALL ERRORS('MARCH4',ERROR)
      CALL EXITS('MARCH4')
      RETURN 1
      END



