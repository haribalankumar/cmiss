      SUBROUTINE MG_COLLOCATE(NBJ,NAQ,NLQ,NWQ,nx,NXQ,
     '  CQ,GCHQ,GUQ,PROPQ,YQ,ERROR,*)

C#### Subroutine: MG_COLLOCATE
C###  Description:
C###    MG_COLLOCATE solves div(kgrad(u))=f by Gauss-Seidel
C###    iteration with multigrid acceleration.
C###    CQ(1,nq) is domain source f
C###    CQ(2,nq) is permeability k11
C###    CQ(3,nq) is permeability k22
C###    CQ(4,nq) is permeability k33
C###
C###    na=1..NMGT are grid levels used in V-cycle (na=1 is finest)
C###    N_GRID1 is # interior grid points in 1D for fine grid (grid 1)
C###    N_GRID2 is # interior grid points in 1D for coarsest grid
C###

C**** NQGE(ng,ne) is global grid pt# nq for local gauss pt ng of ne
C**** NAQ(nq,na)= 0 if nq is on grid level na
C****           = 1  "  "  " betw 2 grid na points in Xi(1) dir.n
C****           = 2  "  "  " betw 2 grid na points in Xi(2) dir.n
C****           = 3  "  "  " betw 2 grid na points in Xi(3) dir.n
C****           = 4  "  "  " betw 4 grid na points in Xi1,2 plane
C****           = 5  "  "  " betw 4 grid na points in Xi2,3 plane
C****           = 6  "  "  " betw 4 grid na points in Xi3,1 plane
C****           = 7  "  "  " betw 8 grid na points in 1,2,3 space
C****           =-1  "  "  " does not belong to grids na or na-1
C**** NWQ defined in IPGRID for each multigrid level na=1..NMGT:
C**** NWQ(1,nq,na)  is 0 if nq is internal
C****    "           mq1 if nq on external bdy (for calc no-flux b.c.)
C**** NWQ(2,nq,na)  is 0 if nq is internal
C****    "           mq2 if nq on external bdy
C**** NWQ(4,nq,na)  is 0 if global point nq is not active
C****    "          1  "    "     "    "  " active
C****    "          2  "    "     "    "  " surrounded by active points
C**** NWQ(5,nq,na)  is 1 for Dirichlet b.c. on boundary point
C****    "             2  "  Neumann    "    "     "      "
C**** nxQ(-3..0..3,i,nq,na) are neighbouring pts for multigrid level na
C**** XQ(nj,nq)     are rect. cart. coords at nq
C**** GUQ(ni1,ni2,nq) are contravariant cpts of metric tensor at nq
C**** GCHQ(nk,nq) are components of GUQ(i,j) * Christoffel_symbol(i,j,k)
C**** dNUdXQ(ni,nj,nq) are derivs of nu coords wrt X coords at nq
C**** dXdXiQ(nj,ni,nq) are derivs of X coords wrt Xi at nq
C**** UQ contain the value of phi(m) and phi(e) at the
C**** grid points of the local quadratic element about nq.
C**** dUdX
C**** d2UdX2

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'solv00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NAQ(NQM,NAM),NLQ(NQM),NWQ(8,0:NQM,NAM),nx,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 CQ(NMM,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),PROPQ(3,3,4,2,NQM),
     '  YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER Iteration,na,nb,NITB,nq
      REAL*8 DT,L2_Norm,L2_Norm_old,L2_Ratio
      LOGICAL ISOTROPIC
      CHARACTER TYPE*7

      CALL ENTERS('MG_COLLOCATE',*9999)
      nb=NBJ(1,1)
      NITB=NIT(nb)

      ISOTROPIC=.false.
      DT=0.0d0
      TYPE='STATIC'

      WRITE(OP_STRING,'('' >Warning, solution type set to static'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C      WRITE(*,'($'' Enter (1)Static,(2)Dynamic,(3)Active: '')')
C      READ(*,*) I
C      IF(I.EQ.1) THEN
C        TYPE='STATIC'
C      ELSE IF(I.EQ.2) THEN
C        TYPE='DYNAMIC'
C      ELSE IF(I.EQ.3) THEN
C        TYPE='ACTIVE'
C      ENDIF
C
C      IF(TYPE(1:7).EQ.'DYNAMIC'.OR.TYPE(1:6).EQ.'ACTIVE') THEN
C        WRITE(*,'($'' Enter time step: '')')
C        READ(*,*) DT
C      ENDIF

      L2_Norm=1.0D0
      L2_Norm_old=0.0D0
      Iteration=0

! Main loop over V-cycles
      DO WHILE(L2_Norm.GT.1.0D-6.AND.Iteration.LT.NTITER)
        Iteration=Iteration+1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_COLLOCATE_1)
          WRITE(OP_STRING,'(/17(''*''),'' V-cycle iteration #'',I3)')
     '      Iteration
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_COLLOCATE_1)
        ENDIF !dop

!     Set RHS=source term at grid level na=1 & RHS=0 for na=2..
        DO nq=1,NQT
          YQ(nq,5,1)=CQ(1,nq)
          DO na=2,NMGT
            YQ(nq,5,na)=0.0D0
          ENDDO !na
        ENDDO !nq

!     Move down grids
        DO na=1,NMGT-1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_COLLOCATE_2)
            WRITE(OP_STRING,'(/17(''*''),'
     '        //''' V-cycle descent: Grid #'',I2,3(''*''))') na
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_COLLOCATE_2)
          ENDIF

!     Relax YQ(nq,1,na) with Relax1 iterations & RHS=YQ(nq,5,na)
          CALL MG_RELAX(Relax1(nx),na,1,NITB,NAQ(1,na),NLQ,
     '      NWQ(1,0,na),NXQ(-NIM,0,0,na),.FALSE.,ISOTROPIC,TYPE,
     '      DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*9999)

!     Soln YQ(nq,1,na)-->Resid YQ(nq,2,na)-->Restr YQ(nq,2,na+1)
          CALL MG_RESIDUAL(na,2,NITB,NAQ(1,na),NLQ,NWQ(1,0,na),
     '      NXQ(-NIM,0,0,na),.FALSE.,ISOTROPIC,TYPE,DT,GCHQ,GUQ,
     '      PROPQ,YQ,ERROR,*9999)
          CALL MG_RESTRICT(na,2,NITB,NAQ,NWQ,NXQ,YQ,ERROR,*9999)

!     Soln YQ(nq,1,na)-->Restr YQ(nq,1,na+1)-->Resid YQ(nq,3,na+1)
          CALL MG_RESTRICT(na  ,1,NITB,NAQ,NWQ,NXQ,YQ,ERROR,*9999)
          CALL MG_RESIDUAL(na+1,3,NITB,NAQ(1,na+1),NLQ,NWQ(1,0,na+1),
     '      NXQ(-NIM,0,0,na+1),.FALSE.,ISOTROPIC,TYPE,DT,GCHQ,GUQ,
     '      PROPQ,YQ,ERROR,*9999)

!     Store restr YQ(nq,1,na+1) in YQ(nq,4,na+1) for upstroke of V
          DO nq=1,NQT
            YQ(nq,4,na+1)=YQ(nq,1,na+1)
          ENDDO

!     Add discretization error to RHS in YQ(nq,5,na+1)
          DO nq=1,NQT !loop over free points of grid na+1
            IF(NAQ(nq,na+1).EQ.0.AND.NWQ(5,nq,na+1).ne.1) THEN
              YQ(nq,5,na+1)=YQ(nq,5,na+1)+YQ(nq,3,na+1)-
     '          YQ(nq,2,na+1)
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_COLLOCATE_3)
                WRITE(OP_STRING,'(''    Error YQ('',I6,'',5,'',I1,'
     '            //''')='',D12.4)') nq,na+1,YQ(nq,5,na+1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_COLLOCATE_3)
              ENDIF
            ENDIF
          ENDDO !nq

        ENDDO !move down grids

!     Put boundary conditions in YQ(nq,2,na)
        DO na=1,NMGT !loop over grids
          DO nq=1,NQT  !loop over boundary points in grid na
            IF(NAQ(nq,na).EQ.0.AND.NWQ(5,nq,na).EQ.1) THEN
              YQ(nq,2,na)=YQ(nq,1,na)
            ENDIF
          ENDDO !nq
        ENDDO !na

!     Solve residual equation on coarsest grid-->Soln YQ(nq,2,na=NMGT)
        na=NMGT
        CALL MG_SOLVE(na,NITB,NAQ(1,na),NWQ(1,0,na),NXQ(-NIM,0,0,na),
     '    ISOTROPIC,TYPE,DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*9999)

!   Move up grids
        DO na=NMGT-1,1,-1
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_COLLOCATE_4)
            WRITE(OP_STRING,'(17(''*''),'
     '        //''' V-cycle  ascent: Grid #'',I2,3(''*''))') na
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_COLLOCATE_4)
          ENDIF

!       Extract coarse error from current coarse sol.n
!                               - projected fine sol.n
          DO nq=1,NQT
            IF(NAQ(nq,na+1).EQ.0) THEN !belongs to grid na+1
              IF(NWQ(1,nq,na+1).EQ.0) THEN !interior point
                YQ(nq,3,na+1)=YQ(nq,2,na+1)-YQ(nq,4,na+1)
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_COLLOCATE_5)
                  WRITE(OP_STRING,'('' YQ('',I6,'',2,'',I1,'')= '','
     '              //'D12.4,'' YQ('',I6,'',4,'',I1,'')= '','
     '              //'D12.4,'' YQ('',I6,'',3,'',I1,'')= '',D12.4)')
     '              nq,na+1,YQ(nq,2,na+1),nq,na+1,YQ(nq,4,na+1),
     '              nq,na+1,YQ(nq,3,na+1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_COLLOCATE_5)
                ENDIF
              ELSE !boundary point
                YQ(nq,3,na+1)=0.0D0 !assume zero error
              ENDIF
            ENDIF
          ENDDO

!       Interpolate error YQ(nq,3,na+1)-->YQ(nq,3,na)
          CALL MG_INTERPOL(na,3,NAQ,NLQ,NXQ,YQ,.FALSE.,ERROR,*9999)

!       Add error YQ(nq,3,na) to previous solution YQ(nq,1,na)
          DO nq=1,NQT !loop over interior points of grid na
            IF(NAQ(nq,na).EQ.0.AND.NWQ(1,nq,na).EQ.0) THEN
              YQ(nq,2,na)=YQ(nq,1,na)+YQ(nq,3,na)
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_COLLOCATE_6)
                WRITE(OP_STRING,'('' Updated soln YQ('',I6,'',2,'','
     '            //'I1,'')= '',D12.4)') nq,na,YQ(nq,2,na)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_COLLOCATE_6)
              ENDIF
            ENDIF
          ENDDO !nq
          DO nq=1,NQT !loop over non-fixed bdry points of grid na
            IF(NAQ(nq,na).EQ.0.AND.NWQ(5,nq,na).EQ.2) THEN
              YQ(nq,2,na)=(4.0d0*YQ(NWQ(1,nq,na),2,na)
     '          -YQ(NWQ(2,nq,na),2,na))/3.0d0
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_COLLOCATE_7)
                WRITE(OP_STRING,'('' Updated soln YQ('',I6,'',2,'','
     '            //'I1,'')= '',D12.4)') nq,na,YQ(nq,2,na)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_COLLOCATE_7)
              ENDIF
            ENDIF
          ENDDO !nq

!       Relax YQ(nq,2,na) with Relax2 iterations & YQ(nq,5,na)
          CALL MG_RELAX(Relax2(nx),na,2,NITB,NAQ(1,na),NLQ,
     '      NWQ(1,0,na),NXQ(-NIM,0,0,na),.FALSE.,ISOTROPIC,TYPE,
     '      DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*9999)

        ENDDO !move up grids

!       Transfer soln to YQ(nq,1,1) & compute fine grid residual norm
        DO nq=1,NQT
          YQ(nq,1,1)=YQ(nq,2,1)
        ENDDO
        CALL MG_RESIDUAL(1,2,NITB,NAQ(1,1),NLQ,NWQ(1,0,1),
     '    NXQ(-NIM,0,0,1),.FALSE.,ISOTROPIC,TYPE,DT,GCHQ,GUQ,
     '    PROPQ,YQ,ERROR,*9999)
        L2_Norm=0.0D0
        DO nq=1,NQT !loop over free points of grid 1
          IF(NWQ(5,nq,1).ne.1) THEN
            L2_Norm=L2_Norm+YQ(nq,2,1)**2
          ENDIF
        ENDDO
        L2_Norm=DSQRT(L2_Norm)
          IF(L2_Norm_old.GT.0.0D0) THEN
            L2_Ratio=L2_Norm/L2_Norm_old
            WRITE(OP_STRING,'( '' V-cycle iteration #'',I3,'
     '        //''' L2_Norm = '',D12.4,'' Ratio = '',F5.3)')
     '        Iteration,L2_Norm,L2_Ratio
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(/'' V-cycle iteration #'',I3,'
     '        //''' L2_Norm = '',D12.4)') Iteration,L2_Norm
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          L2_Norm_old=L2_Norm

      ENDDO !while L2_Norm > 1.d-6

      CALL EXITS('MG_COLLOCATE')
      RETURN
 9999 CALL ERRORS('MG_COLLOCATE',ERROR)
      CALL EXITS('MG_COLLOCATE')
      RETURN 1
      END


