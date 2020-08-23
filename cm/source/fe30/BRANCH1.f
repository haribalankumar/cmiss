      SUBROUTINE BRANCH1(BC_POINTS,BRANCH,CONECT,CQ,ITYP12,nr,
     &  NTIME_POINTS,NTIME_NR,NYNQ,TIME,XQ,YQ,TIME_VALUES,ERROR,*)

C#### Subroutine: BRANCH1
C###  Description:
C###    BRANCH1 calculates the presure, radius and velocity at the grid
C###    points which are at the ends of network segments when
C###    solving the equations for flow in elastic tubes with
C###    finite differences. Two cases are considered, the first
C###    is the inflow or outflow points of the network where the
C###    boundary conditions alone are applied the second is the case for
C###    grid points which are adjacent to a bifurcation. The first uses
C###    equations 3.45 and 3.46 (from the PhD thesis of Nic Smith) to
C###    calculate velocity and radius from a specifed pressure. The
C###    bifurcation equations seek to calculate an initial guess for
C###    the flows in each segment and then using the derivative of
C###    velocity and radius with respect to pressure with the Newton
C###    Method to calculate flows which simultaneously satisfy
C###    conservation of momentum and mass across the bifurcation
C###    (see section 4.1 in the PhD thesis of Nic Smith).

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'time_variable.cmn'

!     Parameter List
      INTEGER  BC_POINTS(3,3),CONECT(-1:1,0:2,NQM),
     '  ITYP12,nr,NYNQ(NHM,NQM,0:NRCM),
     '  NTIME_POINTS(NTIMEVARSM),NTIME_NR(0:NTIMEVARSM,NRM)
C not referenced ,NTIME_INTERP(NTIMEVARSM)
      REAL*8  CQ(NMM,NQM),TIME,XQ(NJM,NQM),YQ(NYQM,NIQM),
     '  TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      CHARACTER ERROR*(*)
C not referenced ,TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
      LOGICAL BRANCH
!     Local Variables
      INTEGER Fsign(3),nj,NJTOT,no_ntv,nq,
     '  nq_a1,nq_a2,nq_a3,nq_b1,nq_b2,nq_b3,nq_c1,nq_c2,nq_c3,
     '  nq_next,ny_p,ny_r,ny_v,
     '  ny_a1p,ny_a1r,ny_a1v,
     '  ny_a2p,ny_a2r,ny_a2v,
     '  ny_a3r,ny_a3v,
     '  ny_b1p,ny_b1r,ny_b1v,
     '  ny_b2p,ny_b2r,ny_b2v,
     '  ny_b3r,ny_b3v,
     '  ny_c1p,ny_c1r,ny_c1v,
     '  ny_c2p,ny_c2r,ny_c2v,
     '  ny_c3r,ny_c3v,
     '  NEW_COUNT,seg,
     '  ntv,ntv1,ntp,NT_CYCLE(NTIMEVARSM)

      REAL*8 Aa,Ab,Ac,B(3),BETA,BETAA1h,BETAB1h,BETAC1h,BETAA1,BETAB1,
     &  BETAC1,DENOM,DENOM1,DELTA(NJT),DELTA_X,DELTA_XA,DELTA_XB,
     &  DELTA_XC,dFA,dFB,dFC,DRA,DRB,DRC,DSA,DSB,DSC,DVA,DVB,DVC,
     &  DG(3,3),DX,DXA,DXB,DXC,ELIP_A,ELIP_B,ELIP_TERM,FA1,FA1_OLD,
     &  FA2,FB1,FB1_OLD,FB2,FC1,FC1_OLD,FC2,fo,G_BC,G_PLEURAL,
     &  G_PLEURAL2,G_TERM,G_TERMA,G_TERMB,G_TERMC,H,HEIGHT(NJT),
     &  HEIGHT2(NJT),I,J,K,L,LAMBDA_1OLD,LAMBDA_1,LAMBDA_2OLD,LAMBDA_2,
     &  LAMBDA_A1,LAMBDA_A2,OLD_LAMBDA_A1,OLD_LAMBDA_A2,LAMBDA_B1,
     &  LAMBDA_B2,OLD_LAMBDA_B1,OLD_LAMBDA_B2,LAMBDA_C1,LAMBDA_C2,
     &  OLD_LAMBDA_C1,OLD_LAMBDA_C2,LAM_TERMA,LAM_TERMB,LAM_TERMC,
     &  LAMBDA_AH,LAMBDA_BH,LAMBDA_CH,LA,LB,LC,LOOP_DENOM,OLDFA1,OLDFA2,
     &  M,OLDFB1,OLDFB2,OLDFC1,OLDFC2,OLDP_INITIAL,PA,PB,PC,P_DIFF,
     &  P_INITIAL,PLEURAL_DENSITY,PRESTERM,ra,rb,rc,RO_ART,RO_VIEN,ROI,
     &  ROIA1,ROIB1,ROIC1,ROIA1h,ROIB1h,ROIC1h,TERM1,TERM2,TERM3,TERM4,
     &  TERM5,TERMREST,TRACE,TRACE_A1,TRACE_B1,TRACE_C1,TRACE_A2,
     &  TRACE_B2,TRACE_C2,OLD_TRACE_A1,OLD_TRACE_B1,OLD_TRACE_C1,
     &  OLD_TRACE_A2,OLD_TRACE_B2,OLD_TRACE_C2,SA,SB,SC,SUM,VA,VB,VC,
     &  Va1_new,Va2_new,Va3_mid,Va1_old,Va2_old,Vb1_new,Vb2_new,Vb3_mid,
     &  Vb1_old,Vb2_old,Vc1_new,Vc2_new,Vc3_mid,Vc1_old,Vc2_old,
     &  vien_ro_a1,vien_ro_a1h,vien_ro_b1,vien_ro_b1h,vien_ro_c1,
     &  vien_ro_c1h,X(3),Ain(NTIMEPOINTSM,NTIMEVARSM),Bin(NTIMEPOINTSM,
     &  NTIMEVARSM),T_frac(NTIMEPOINTSM,NTIMEVARSM),T_period(NTIMEVARSM)
C      CHARACTER TRANS
C      PARAMETER (TRANS='N')
      EXTERNAL DGETRS,DGETRF

      CALL ENTERS('BRANCH1',*9999)
      
      NJTOT=NJ_LOC(NJL_GEOM,0,nr)
C Set Fsign(seg) to be +1 if flow +ve into bifurcation, else -1
      DO seg=1,3
        IF(BC_POINTS(seg,3).EQ.BC_POINTS(seg,2)) THEN
          Fsign(seg)= 1
        ELSE
          Fsign(seg)=-1
        ENDIF
      ENDDO !seg
      nq=BC_POINTS(1,1) !pt a1 (1st segment & adjacent)
      IF(.NOT.BRANCH) THEN !terminal or start point (not bifurcation)
C      IF(NXQ(1,0,nq,1).LE.1.AND.NXQ(-1,0,nq,1).LE.1) THEN
        IF(NQ_START(nr).EQ.nq) THEN !start point
C         Note: terminal pts dealt with in BRANCH2
          nq_next=CONECT(1,1,nq) !finds the grid point adjacent to the
C          nq_next=NXQ(1,1,nq,1)
          ny_p=NYNQ(1,nq,0) !grid point on the boundary and
          ny_r=NYNQ(2,nq,0) !determines the ny values for
          ny_v=NYNQ(3,nq,0) !velocity, radius, and pressure
C calculates the streach ratios for the k and k+1 time step
          LAMBDA_1OLD=(YQ(ny_v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(ny_v,9)*(T1-TIME+TINCR)/(T1-T0))
          LAMBDA_1=(YQ(ny_v,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_v,9)*(T1-TIME)/(T1-T0))
          LAMBDA_2OLD=(YQ(NYNQ(3,nq_next,0),2)
     '      *(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(NYNQ(3,nq_next,0),9)*(T1-TIME+TINCR)/(T1-T0))
          LAMBDA_2=(YQ(NYNQ(3,nq_next,0),2)*(TIME-T0)/(T1-T0))+
     '      (YQ(NYNQ(3,nq_next,0),9)*(T1-TIME)/(T1-T0))
          DELTA_X=0.0d0
          DO nj=1,NJTOT !calculate delta x for half time step
            DELTA_X=DELTA_X+((XQ(nj,nq)-XQ(nj,nq_next))**2.0d0)
          ENDDO !nj
C calculates the space step between grid points
C          DELTA_X=DMAX1((DELTA_X**0.5d0),0.60D0)*
C     '      ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)*0.25D0)
          DELTA_X=(DELTA_X**0.5d0)*
     '      ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)*0.25D0)
C linear interpolation of trace though time
          TRACE=(YQ(ny_p,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_p,9)*(T1-TIME)/(T1-T0))
          
C... KSB Oct-03 : Calculating gravity term
          DX=0.d0
          DO nj=1,NJT
            DX=DX+(XQ(nj,nq_next)-XQ(nj,nq))**2.d0
            DELTA(nj)=XQ(nj,nq_next)-XQ(nj,nq)
            HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq)
          ENDDO
          DX=DSQRT(DX)
          G_BC=0.d0 !gravitational force pressure boundary conditions
          G_TERM=0.d0 !gravity term for pressure equation
          DO nj=1,NJT 
            G_BC=G_BC+CQ(1,nq)*G_VECTOR(nj)*GRAVITY*1000.d0*HEIGHT(nj)
            G_TERM=G_TERM+CQ(1,nq)*G_VECTOR(nj)*GRAVITY*1000.d0
     &        *DELTA(nj)
          ENDDO
C... calculating external pressure force (=TRACE):          
C linear interpolation of trace though time
          IF(ITYP12.EQ.2) THEN !pulmonary flow - calculate pleural pressure
            G_PLEURAL=0.d0 !gravitational force - initialise
            PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
 !This density value effectively represents lung tissue density = 1/5th that of water
            DO nj=1,NJT
              G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &          *GRAVITY*1000.d0*HEIGHT(nj)
            ENDDO
            TRACE=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
            !NB/ currently not time-dependent
          ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
            nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
            TRACE=CQ(nj,nq)/1000.d0 !Pa->kPa                                
          ELSE !coronary or other problems
            TRACE=(YQ(ny_p,2)*(TIME-T0)/(T1-T0))+
     '        (YQ(ny_p,9)*(T1-TIME)/(T1-T0))
          ENDIF
          
C PM 26-JUL-01 : boundary conditions
          IF(KTYP3_INIT(nr).EQ.3) THEN ! from a file
            READ(IOFILE1,*) YQ(ny_p,1)
          ELSEIF(KTYP3_INIT(nr).EQ.4) THEN ! from time variable defined
 ! through iptime
C            DO ntv=1,NTIMEVARST
            DO no_ntv=1,NTIME_NR(0,nr) !KSB region-dep boundary conditions, only loop over variables in region
              ntv=NTIME_NR(no_ntv,nr)                                          
              DO ntp=1,NTIME_POINTS(ntv)-1
                Ain(ntp,ntv)=TIME_VALUES(2,ntp,ntv)-
     '            TIME_VALUES(1,ntp,ntv)*
     '            (TIME_VALUES(2,ntp+1,ntv)-
     '            TIME_VALUES(2,ntp,ntv))/
     '            (TIME_VALUES(1,ntp+1,ntv)-
     '            TIME_VALUES(1,ntp,ntv))
                Bin(ntp,ntv)=(TIME_VALUES(2,ntp+1,ntv)-
     '            TIME_VALUES(2,ntp,ntv))/
     '            (TIME_VALUES(1,ntp+1,ntv)-
     '            TIME_VALUES(1,ntp,ntv))
              ENDDO
            ENDDO
C            DO ntv=1,NTIMEVARST
            DO no_ntv=1,NTIME_NR(0,nr) !NTIMEVARST
              ntv=NTIME_NR(no_ntv,nr)               
              DO ntp=1,NTIME_POINTS(ntv)
                T_frac(ntp,ntv)=TIME_VALUES(1,ntp,ntv)/
     '            (TIME_VALUES(1,NTIME_POINTS(ntv),ntv)-
     '            TIME_VALUES(1,1,ntv))
              ENDDO
            ENDDO
C            DO ntv=1,NTIMEVARST
            DO no_ntv=1,NTIME_NR(0,nr) !NTIMEVARST
              ntv=NTIME_NR(no_ntv,nr)                           
              T_period(ntv)=TIME_VALUES(1,NTIME_POINTS(ntv),ntv)-
     '          TIME_VALUES(1,1,ntv)
              NT_CYCLE(ntv)=INT(TIME/T_period(ntv)+1.0d-7)
            ENDDO
            ntv1=NTIME_NR(1,nr)
            ntv=NTIME_NR(2,nr) !KSB: region dependent boundary conditions            

C... KSB - replaced hard coded ntv values 01/12/05            
            IF((PERIODIC.EQ.'N').OR.(PERIODIC.EQ.'n')) THEN
              IF(TIME.LT.TIME_VALUES(1,1,ntv1)) THEN
                YQ(ny_p,1)=TIME_VALUES(2,0,ntv1)+G_BC
              ELSEIF (TIME.GT.TIME_VALUES(1,NTIME_POINTS(ntv1),ntv1))
     &            THEN
                YQ(ny_p,1)=TIME_VALUES(2,NTIME_POINTS(ntv1)+1,ntv1)+G_BC
              ELSE
                DO ntp=1,NTIME_POINTS(ntv1)-1
                  IF((TIME.GE.TIME_VALUES(1,1,ntv1)+T_period(ntv1)*
     '              T_frac(ntp,ntv1)).AND.(TIME.LT.TIME_VALUES(1,1,ntv1)
     &              +T_period(ntv1)*T_frac(ntp+1,ntv1))) THEN
                    YQ(ny_p,1)=Ain(ntp,ntv1)+Bin(ntp,ntv1)*TIME+G_BC
                  ENDIF
                ENDDO ! ntp
              ENDIF
            ELSE
              DO ntp=1,NTIME_POINTS(ntv1)-1
                IF((TIME.GE.T_period(ntv1)*(NT_CYCLE(ntv1)+
     '            T_frac(ntp,ntv1))).AND.(TIME.LT.T_period(ntv1)*
     '            (NT_CYCLE(ntv1)+T_frac(ntp+1,ntv1)))) THEN
                  YQ(ny_p,1)=Ain(ntp,ntv1)+Bin(ntp,ntv1)*
     '              (TIME-NT_CYCLE(ntv1)*T_period(ntv1))+G_BC
                ENDIF
              ENDDO ! ntp
            ENDIF ! PERIODIC
          ELSE
            
C linear interpolation of the boundary condition pressure through time
C PM 26-JUL-01 : now can be defined through iptime
C              IF(TIME.LT.((T1+T0)/1.0d0)) THEN
C                YQ(ny_p,1)=YQ(ny_p,3)+ENTRY_PRESSURE0+
C     '            ((ENTRY_PRESSURE1-ENTRY_PRESSURE0)/(T1-T0))*
C     '            ((TIME-T0)*1.0d0) !index here as well
C                WRITE(*,*) YQ(ny_p,1)
C              ELSE
C           YQ(ny_p,1)=YQ(ny_p,3)+ENTRY_PRESSURE1
C           ENDIF
          ENDIF ! KTYP_INIT
          
          
C calculates the unstressed arterial radius at time k+1
          ROI=CQ(4,nq)/(LAMBDA_1**0.5d0)
          BETA=DMAX1((CQ(6,nq)*LAMBDA_1+CQ(12,nq)),1.0d0)
C calculates the radius using the pressure boundary condition
          IF((YQ(ny_p,1)+TRACE).GT.0.0D0) THEN
            YQ(ny_r,1)=ROI*((((YQ(ny_p,1)+TRACE)
     '        /CQ(5,nq))+1.0d0)**(1.0d0/BETA))
          ELSE
            fo=BETA*CQ(5,nq)/CQ(7,nq)
            YQ(ny_r,1)=ROI*
     '        (1-((YQ(ny_p,1)+TRACE)/fo))**(-1.0d0/CQ(7,nq))
          ENDIF
C calculates the indvidual terms in the velocity boundary condition
C equation
          TERM1=((2.0d0*CQ(3,nq))-0.0d0)*((YQ(ny_v,5)**2.0d0)
     '      /YQ(ny_r,5))*
     '      (YQ(NYNQ(2,nq_next,0),1)+YQ(NYNQ(2,nq_next,0),8)
     '      -YQ(ny_r,1)-YQ(ny_r,8))
          P_DIFF=YQ(NYNQ(1,nq_next,0),1)+YQ(NYNQ(1,nq_next,0),8)
     '      -YQ(ny_p,1)-YQ(ny_p,8)+2.d0*G_TERM!G_TERM*DX*ANGLE
          
          TERM2=(1.0d0/CQ(1,nq))*P_DIFF          
          TERM3=2.0d0*TINCR*CQ(2,nq)*(CQ(3,nq)/(CQ(3,nq)-1.d0))*
     '      (YQ(ny_v,5)/(YQ(ny_r,5)**2.0d0))
          TERM4=(1.0d0/YQ(ny_r,5))*((DELTA_X/TINCR)+
     '      ((2.0d0*CQ(3,nq)-1.0d0)*YQ(ny_v,5)))*
     '      (((LAMBDA_1**0.5d0)*YQ(ny_r,1))+
     '      ((LAMBDA_2**0.5d0)*YQ(NYNQ(2,nq_next,0),1))-
     '      ((LAMBDA_1OLD**0.5d0)*YQ(ny_r,8))
     '      -((LAMBDA_2OLD**0.5d0)*YQ(NYNQ(2,nq_next,0),8)))
          TERM5=YQ(ny_v,5)/YQ(ny_r,5)
     '      *(YQ(NYNQ(2,nq_next,0),1)
     '      +YQ(NYNQ(2,nq_next,0),8)-YQ(ny_r,1)-YQ(ny_r,8))
C calculated the arterial velocity at the inflow arterial grid point
          YQ(ny_v,1)=YQ(NYNQ(3,nq_next,0),8)+
     '      ((TINCR/(2.0d0*DELTA_X))*(TERM1-TERM2))
     '      -TERM3+TERM4+TERM5
          
          IF (((VENOUS_NETWORK.EQ.'Y').OR.
     '      (VENOUS_NETWORK.EQ.'y')).AND.(N_VENOUS_GEOM.EQ.1)) THEN
            
C at the same grid point the venous quanties are now calculated using
C equivalent boundary condtion equations.
            ny_p=NYNQ(4,nq,0)
            ny_r=NYNQ(5,nq,0)
            ny_v=NYNQ(6,nq,0)
C calculates the unstressed venous radius at time k+1
            ROI=CQ(4,nq)/(LAMBDA_1**0.5d0) !ro at the time k+1 vien
            BETA=CQ(9,nq)*LAMBDA_1+CQ(10,nq) !beta at the time k+1
            BETA=DMAX1(BETA,1.0d0)
C linear interpolation of the venous boundary condition pressure
C through time
C PM 26-JUL-01 boundary conditions
            IF(KTYP3_INIT(nr).EQ.3) THEN ! from a file
              READ(IOFILE1,*) YQ(ny_p,1)
            ELSEIF(KTYP3_INIT(nr).EQ.4) THEN ! defined as time variable
 ! through iptime
C              DO ntv=1,NTIMEVARST
              DO no_ntv=1,NTIME_NR(0,nr) !NTIMEVARST
                ntv=NTIME_NR(no_ntv,nr)      
                DO ntp=1,NTIME_POINTS(ntv)-1
                  Ain(ntp,ntv)=TIME_VALUES(2,ntp,ntv)-
     '              TIME_VALUES(1,ntp,ntv)*
     '              (TIME_VALUES(2,ntp+1,ntv)-
     '              TIME_VALUES(2,ntp,ntv))/
     '              (TIME_VALUES(1,ntp+1,ntv)-
     '              TIME_VALUES(1,ntp,ntv))
                  Bin(ntp,ntv)=(TIME_VALUES(2,ntp+1,ntv)-
     '              TIME_VALUES(2,ntp,ntv))/
     '              (TIME_VALUES(1,ntp+1,ntv)-
     '              TIME_VALUES(1,ntp,ntv))
                ENDDO
              ENDDO
C              DO ntv=1,NTIMEVARST
              DO no_ntv=1,NTIME_NR(0,nr) !NTIMEVARST
                ntv=NTIME_NR(no_ntv,nr)      
                DO ntp=1,NTIME_POINTS(ntv)
                  T_frac(ntp,ntv)=TIME_VALUES(1,ntp,ntv)/
     '              (TIME_VALUES(1,NTIME_POINTS(ntv),ntv)-
     '              TIME_VALUES(1,1,ntv))
                ENDDO
              ENDDO
C              DO ntv=1,NTIMEVARST
              DO no_ntv=1,NTIME_NR(0,nr) !NTIMEVARST
                ntv=NTIME_NR(no_ntv,nr)      
                T_period(ntv)=TIME_VALUES(1,NTIME_POINTS(ntv),ntv)-
     '            TIME_VALUES(1,1,ntv)
                NT_CYCLE(ntv)=INT(TIME/T_period(ntv)+1.0d-7)
              ENDDO
              
              IF((PERIODIC.EQ.'N').OR.(PERIODIC.EQ.'n')) THEN
                ntv1=NTIME_NR(1,nr)
                ntv=NTIME_NR(2,nr) !KSB: region dependent boundary conditions
                IF(TIME.LT.TIME_VALUES(1,1,ntv)) THEN
                  YQ(ny_p,1)=TIME_VALUES(2,0,ntv)+G_BC
                ELSEIF (TIME.GT.TIME_VALUES(1,NTIME_POINTS(ntv1),ntv))
     '              THEN
                  YQ(ny_p,1)=TIME_VALUES(2,NTIME_POINTS(ntv1)+1,ntv)
     &              +G_BC
                ELSE
                  DO ntp=1,NTIME_POINTS(ntv)-1
                    IF((TIME.GE.TIME_VALUES(1,1,ntv)+T_period(ntv)*
     '                T_frac(ntp,ntv)).AND.
     '                (TIME.LT.TIME_VALUES(1,1,ntv)+
     '                T_period(ntv)*T_frac(ntp+1,ntv))) THEN
                      YQ(ny_p,1)=Ain(ntp,ntv)+Bin(ntp,ntv)*TIME+G_BC
                    ENDIF
                  ENDDO ! ntp
                ENDIF
              ELSE
                DO ntp=1,NTIME_POINTS(ntv)-1
                  IF((TIME.GE.T_period(ntv)*(NT_CYCLE(ntv)+
     '              T_frac(ntp,ntv))).AND.(TIME.LT.T_period(ntv)*
     '              (NT_CYCLE(ntv)+T_frac(ntp+1,ntv)))) THEN
                    YQ(ny_p,1)=Ain(ntp,ntv)+Bin(ntp,ntv)*
     '                (TIME-NT_CYCLE(ntv)*T_period(ntv))+G_BC
                  ENDIF
                ENDDO ! ntp
              ENDIF ! PERIODIC
            ELSE
              
C             IF(TIME.LT.((T1+T0)/1.0d0)) THEN
C             YQ(ny_p,1)=YQ(ny_p,3)+EXIT_PRESSURE0+
C     '            ((EXIT_PRESSURE1-EXIT_PRESSURE0)/(T1-T0))*
C     '            ((TIME-T0)*1.0d0)
C                WRITE(*,*) YQ(ny_p,1)
C              ELSE
C                YQ(ny_p,1)=YQ(ny_p,3)+EXIT_PRESSURE1
C              ENDIF
            ENDIF
C calculates the venous radius using the pressure boundary condition
            IF((YQ(ny_p,1)+TRACE).GT.0.0D0) THEN
              YQ(ny_r,1)=(ROI*VIEN_RATIO)*((((YQ(ny_p,1)+TRACE)
     '          /CQ(11,nq))+1.0d0)**(1.0d0/BETA))
            ELSE
              fo=BETA*CQ(11,nq)/CQ(7,nq)
              YQ(ny_r,1)=ROI*VIEN_RATIO*
     '          (1-((YQ(ny_p,1)+TRACE)/fo))
     '          **(-1.0d0/CQ(7,nq))
            ENDIF
C calculates the indvidual terms in the velocity boundary condition
C equation
            TERM1=((2.0d0*CQ(3,nq))-0.0d0)*((YQ(ny_v,5)**2.0d0)
     '        /YQ(ny_r,5))*
     '        (YQ(NYNQ(5,nq_next,0),1)+YQ(NYNQ(5,nq_next,0),8)
     '        -YQ(ny_r,1)-YQ(ny_r,8))            
            P_DIFF=YQ(NYNQ(4,nq_next,0),1)+YQ(NYNQ(4,nq_next,0),8)
     '        -YQ(ny_p,1)-YQ(ny_p,8)+2.d0*G_TERM            
            TERM2=(1.0D0/CQ(1,nq))*P_DIFF
            TERM3=2.0d0*TINCR*CQ(2,nq)*(CQ(3,nq)/(CQ(3,nq)-1))*
     '        (YQ(ny_v,5)/(YQ(ny_r,5)**2.0d0))
            TERM4=(1.0d0/YQ(ny_r,5))*((DELTA_X/TINCR)+
     '        ((2.0d0*CQ(3,nq)-1.0d0)*YQ(ny_v,5)))*
     '        (((LAMBDA_1**0.5d0)*YQ(ny_r,1))+
     '        ((LAMBDA_2**0.5d0)*YQ(NYNQ(5,nq_next,0),1))-
     '        ((LAMBDA_1OLD**0.5d0)*YQ(ny_r,8))
     '        -((LAMBDA_2OLD**0.5d0)*YQ(NYNQ(5,nq_next,0),8)))
            TERM5=YQ(ny_v,5)/YQ(ny_r,5)*(YQ(NYNQ(5,nq_next,0),1)
     '        +YQ(NYNQ(5,nq_next,0),8)-YQ(ny_r,1)-YQ(ny_r,8))
C calculated the venous velocity at the inflow arterial grid point
            YQ(ny_v,1)=YQ(NYNQ(6,nq_next,0),8)+
     '        ((TINCR/(2.0d0*DELTA_X))*(TERM1-TERM2))
     '        -TERM3+TERM4+TERM5
            
          ENDIF !start point
          
        ENDIF ! venous-network & identical venous net_work
        
      ELSE !bifurcation
C -------------- now do arterial bifurcation -----------------
C do arterial bifurcation first
C need to work out the six points that surround a bifurcation
C a1,b1,c1 are pts adjacent to bifurcation
C a2,b2,c2 are pts one step away from bifurcation
C a3,b3,c3 are pts where half-step information is stored
C Calculate initial guess at pressures for bifurcation pts a1.b1,c1
C determine the grid points a1 and a2 and the ny values corresponding
C to velocity, radius and pressure
        nq_a1=BC_POINTS(1,1) !was CONECT(0,1,nq)
        ny_a1p=NYNQ(1,nq_a1,0)
        ny_a1r=NYNQ(2,nq_a1,0)
        ny_a1v=NYNQ(3,nq_a1,0)
        nq_a2=BC_POINTS(1,2) !was CONECT(-1,1,nq_a1)
        ny_a2p=NYNQ(1,nq_a2,0)
        ny_a2r=NYNQ(2,nq_a2,0)
        ny_a2v=NYNQ(3,nq_a2,0)
        nq_a3=BC_POINTS(1,3)
        ny_a3r=NYNQ(2,nq_a3,0)
        ny_a3v=NYNQ(3,nq_a3,0)            
C determining the velocity at the grid a1 and a2 at the bifurcation
        Va1_new=Fsign(1)*YQ(ny_a1v,1) !is new veloc at a1
        Va2_new=Fsign(1)*YQ(ny_a2v,1) !is new veloc at a2
        Va3_mid=Fsign(1)*YQ(ny_a3v,5) !is mid veloc at a3
        Va1_old=Fsign(1)*YQ(ny_a1v,8) !is old veloc at a1
        Va2_old=Fsign(1)*YQ(ny_a2v,8) !is old veloc at a2
C determine the grid points b1 and b2 and the ny values corresponding
C to velocity, radius and pressure
        nq_b1=BC_POINTS(2,1) !was CONECT(1,1,nq)
        ny_b1p=NYNQ(1,nq_b1,0)
        ny_b1r=NYNQ(2,nq_b1,0)
        ny_b1v=NYNQ(3,nq_b1,0)
        nq_b2=BC_POINTS(2,2) !was CONECT(1,1,nq_b1)
        ny_b2p=NYNQ(1,nq_b2,0)
        ny_b2r=NYNQ(2,nq_b2,0)
        ny_b2v=NYNQ(3,nq_b2,0)
        nq_b3=BC_POINTS(2,3)
        ny_b3r=NYNQ(2,nq_b3,0)
        ny_b3v=NYNQ(3,nq_b3,0)
C determining the velocity at the grid b1 and b2 at the bifurcation
        Vb1_new=-Fsign(2)*YQ(ny_b1v,1) !is new veloc at b1
        Vb2_new=-Fsign(2)*YQ(ny_b2v,1) !is new veloc at b2
        Vb3_mid=-Fsign(2)*YQ(ny_b3v,5) !is mid veloc at b3
        Vb1_old=-Fsign(2)*YQ(ny_b1v,8) !is old veloc at b1
        Vb2_old=-Fsign(2)*YQ(ny_b2v,8) !is old veloc at b2
C determine the grid points c1 and c2 and the ny values corresponding
C to velocity, radius and pressure
        nq_c1=BC_POINTS(3,1) !was CONECT(1,2,nq)
        ny_c1p=NYNQ(1,nq_c1,0)
        ny_c1r=NYNQ(2,nq_c1,0)
        ny_c1v=NYNQ(3,nq_c1,0)
        nq_c2=BC_POINTS(3,2) !was CONECT(1,1,nq_c1)
        ny_c2p=NYNQ(1,nq_c2,0)
        ny_c2r=NYNQ(2,nq_c2,0)
        ny_c2v=NYNQ(3,nq_c2,0)
        nq_c3=BC_POINTS(3,3)
        ny_c3r=NYNQ(2,nq_c3,0)
        ny_c3v=NYNQ(3,nq_c3,0)
C determining the velocity at the grid c1 and c2 at the bifurcation
        Vc1_new=-Fsign(3)*YQ(ny_c1v,1) !is new veloc at c1
        Vc2_new=-Fsign(3)*YQ(ny_c2v,1) !is new veloc at c2
        Vc3_mid=-Fsign(3)*YQ(ny_c3v,5) !is mid veloc at c3
        Vc1_old=-Fsign(3)*YQ(ny_c1v,8) !is old veloc at c1
        Vc2_old=-Fsign(3)*YQ(ny_c2v,8) !is old veloc at c2
        
C.. KSB Oct 2003 - calculate gravity terms for each segment at
C.. bifurcation
        IF(ITYP12.EQ.2) PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
        DXA=0.d0
        G_TERMA=0.d0
        DO nj=1,NJT
          DXA=DXA+(XQ(nj,nq_a2)-XQ(nj,nq_a1))**2.d0
          DELTA(nj)=XQ(nj,nq_a2)-XQ(nj,nq_a1)
          HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq_a1)
          HEIGHT2(nj)=XQ_IN(nj)-XQ(nj,nq_a2)
        ENDDO
        DXA=DSQRT(DXA)
        DO nj=1,NJT !sums gravity vector components
          G_TERMA=G_TERMA+CQ(1,nq_a1)*G_VECTOR(nj)*GRAVITY*1000.d0
     &      *DELTA(nj)/DXA
        ENDDO
        IF(ITYP12.EQ.2) THEN !pulmonary flow - calculate pleural pressure
          G_PLEURAL=0.d0 !gravitational force - initialise
          G_PLEURAL2=0.d0
          DO nj=1,NJT
            G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &        *GRAVITY*1000.d0*HEIGHT(nj)
            G_PLEURAL2=G_PLEURAL2+PLEURAL_DENSITY*G_VECTOR(nj)
     &        *GRAVITY*1000.d0*HEIGHT2(nj)
          ENDDO
          TRACE_A1=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
          TRACE_A2=-(PLEURAL_P+G_PLEURAL2)
          OLD_TRACE_A1=TRACE_A1 !currently not time-dependent
          OLD_TRACE_A2=TRACE_A2
        ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
          nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
          TRACE_A1=CQ(nj,nq_a1)/1000.d0 !Pa->kPa 
          TRACE_A2=CQ(nj,nq_a2)/1000.d0 !Pa->kPa
          OLD_TRACE_A1=TRACE_A1 !currently not time-dependent
          OLD_TRACE_A2=TRACE_A2
        ENDIF
        
        DXB=0.d0
        G_TERMB=0.d0
        DO nj=1,NJT
          DXB=DXB+(XQ(nj,nq_b2)-XQ(nj,nq_b1))**2.d0
          DELTA(nj)=XQ(nj,nq_b2)-XQ(nj,nq_b1)
          HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq_b1)
          HEIGHT2(nj)=XQ_IN(nj)-XQ(nj,nq_b2)
        ENDDO
        DXB=DSQRT(DXB)
        DO nj=1,NJT !sums gravity vector components
          G_TERMB=G_TERMB+CQ(1,nq_b1)*G_VECTOR(nj)*GRAVITY*1000.d0
     &      *DELTA(nj)/DXB
        ENDDO
        IF(ITYP12.EQ.2) THEN !pulmonary flow - calculate pleural pressure
          G_PLEURAL=0.d0 !gravitational force - initialise
          G_PLEURAL2=0.d0 
          DO nj=1,NJT
            G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &        *GRAVITY*1000.d0*HEIGHT(nj)
            G_PLEURAL2=G_PLEURAL2+PLEURAL_DENSITY*G_VECTOR(nj)
     &        *GRAVITY*1000.d0*HEIGHT2(nj)
          ENDDO
          TRACE_B1=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
          TRACE_B2=-(PLEURAL_P+G_PLEURAL2) 
          OLD_TRACE_B1=TRACE_B1 !currently not time-dependent
          OLD_TRACE_B2=TRACE_B2
        ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
          nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
          TRACE_B1=CQ(nj,nq_b1)/1000.d0 !Pa->kPa 
          TRACE_B2=CQ(nj,nq_b2)/1000.d0 !Pa->kPa
          OLD_TRACE_B1=TRACE_B1 !currently not time-dependent
          OLD_TRACE_B2=TRACE_B2
        ENDIF

        DXC=0.d0
        G_TERMC=0.d0
        DO nj=1,NJT
          DXC=DXC+(XQ(nj,nq_c2)-XQ(nj,nq_c1))**2.d0
          DELTA(nj)=XQ(nj,nq_c2)-XQ(nj,nq_c1)
          HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq_c1)
          HEIGHT2(nj)=XQ_IN(nj)-XQ(nj,nq_c2)
        ENDDO
        DXC=DSQRT(DXC)
        DO nj=1,NJT !sums gravity vector components
          G_TERMC=G_TERMC+CQ(1,nq_c1)*G_VECTOR(nj)*GRAVITY*1000.d0
     &      *DELTA(nj)/DXC
        ENDDO
        IF(ITYP12.EQ.2) THEN !pulmonary flow - calculate pleural pressure
          G_PLEURAL=0.d0 !gravitational force - initialise
          G_PLEURAL2=0.d0 
          DO nj=1,NJT
            G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &        *GRAVITY*1000.d0*HEIGHT(nj)
            G_PLEURAL2=G_PLEURAL2+PLEURAL_DENSITY*G_VECTOR(nj)
     &        *GRAVITY*1000.d0*HEIGHT2(nj)
          ENDDO
          TRACE_C1=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
          TRACE_C2=-(PLEURAL_P+G_PLEURAL2) 
          OLD_TRACE_C1=TRACE_C1 !currently not time-dependent
          OLD_TRACE_C2=TRACE_C2
        ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
          nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
          TRACE_C1=CQ(nj,nq_c1)/1000.d0 !Pa->kPa 
          TRACE_C2=CQ(nj,nq_c2)/1000.d0 !Pa->kPa
          OLD_TRACE_C1=TRACE_C1 !currently not time-dependent
          OLD_TRACE_C2=TRACE_C2
        ENDIF
        
        IF(ITYP12.EQ.1) THEN !coronary or other problems (not pulmonary)
C calculate the trace at the points surrounding the bifurcation
C by linearly interpolating linearly through time
          TRACE_A1=(YQ(ny_a1p,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_a1p,9)*(T1-TIME)/(T1-T0))
          TRACE_A2=(YQ(ny_a2p,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_a2p,9)*(T1-TIME)/(T1-T0))
          OLD_TRACE_A1=(YQ(ny_a1p,2)*(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(ny_a1p,9)*(T1-TIME+TINCR)/(T1-T0))
          OLD_TRACE_A2=(YQ(ny_a2p,2)*(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(ny_a2p,9)*(T1-TIME+TINCR)/(T1-T0))
          TRACE_B1=(YQ(ny_b1p,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_b1p,9)*(T1-TIME)/(T1-T0))
          TRACE_B2=(YQ(ny_b2p,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_b2p,9)*(T1-TIME)/(T1-T0))
          OLD_TRACE_B1=(YQ(ny_b1p,2)*(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(ny_b1p,9)*(T1-TIME+TINCR)/(T1-T0))
          OLD_TRACE_B2=(YQ(ny_b2p,2)*(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(ny_b2p,9)*(T1-TIME+TINCR)/(T1-T0))
          TRACE_C1=(YQ(ny_c1p,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_c1p,9)*(T1-TIME)/(T1-T0))
          TRACE_C2=(YQ(ny_c2p,2)*(TIME-T0)/(T1-T0))+
     '      (YQ(ny_c2p,9)*(T1-TIME)/(T1-T0))
          OLD_TRACE_C1=(YQ(ny_c1p,2)*(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(ny_c1p,9)*(T1-TIME+TINCR)/(T1-T0))
          OLD_TRACE_C2=(YQ(ny_c2p,2)*(TIME-T0-TINCR)/(T1-T0))+
     '      (YQ(ny_c2p,9)*(T1-TIME+TINCR)/(T1-T0))
        ENDIF
C calculate the extension ratios  at the points surrounding
C the bifurcation by linearly interpolating linearly through time
        LAMBDA_A1=(YQ(ny_a1v,2)*(TIME-T0)/(T1-T0))+
     '    (YQ(ny_a1v,9)*(T1-TIME)/(T1-T0))
        LAMBDA_A2=(YQ(ny_a2v,2)*(TIME-T0)/(T1-T0))+
     '    (YQ(ny_a2v,9)*(T1-TIME)/(T1-T0))
        OLD_LAMBDA_A1=(YQ(ny_a1v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '    (YQ(ny_a1v,9)*(T1-TIME+TINCR)/(T1-T0))
        OLD_LAMBDA_A2=(YQ(ny_a2v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '    (YQ(ny_a2v,9)*(T1-TIME+TINCR)/(T1-T0))
        LAMBDA_AH=(LAMBDA_A1+LAMBDA_A2+
     '    OLD_LAMBDA_A1+OLD_LAMBDA_A2)/4.0d0
        LAMBDA_B1=(YQ(ny_b1v,2)*(TIME-T0)/(T1-T0))+
     '    (YQ(ny_b1v,9)*(T1-TIME)/(T1-T0))
        LAMBDA_B2=(YQ(ny_b2v,2)*(TIME-T0)/(T1-T0))+
     '    (YQ(ny_b2v,9)*(T1-TIME)/(T1-T0))
        OLD_LAMBDA_B1=(YQ(ny_b1v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '    (YQ(ny_b1v,9)*(T1-TIME+TINCR)/(T1-T0))
        OLD_LAMBDA_B2=(YQ(ny_b2v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '    (YQ(ny_b2v,9)*(T1-TIME+TINCR)/(T1-T0))
        LAMBDA_BH=(LAMBDA_B1+LAMBDA_B2+
     '    OLD_LAMBDA_B1+OLD_LAMBDA_B2)/4.0d0
        LAMBDA_C1=(YQ(ny_c1v,2)*(TIME-T0)/(T1-T0))+
     '    (YQ(ny_c1v,9)*(T1-TIME)/(T1-T0))
        LAMBDA_C2=(YQ(ny_c2v,2)*(TIME-T0)/(T1-T0))+
     '    (YQ(ny_c2v,9)*(T1-TIME)/(T1-T0))
        OLD_LAMBDA_C1=(YQ(ny_c1v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '    (YQ(ny_c1v,9)*(T1-TIME+TINCR)/(T1-T0))
        OLD_LAMBDA_C2=(YQ(ny_c2v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '    (YQ(ny_c2v,9)*(T1-TIME+TINCR)/(T1-T0))
        LAMBDA_CH=(LAMBDA_C1+LAMBDA_C2+
     '    OLD_LAMBDA_C1+OLD_LAMBDA_C2)/4.0d0
C calculating the unstressed radii at the points surrounding the
C bifurcation
        ROIA1h=(CQ(4,nq_a1)+CQ(4,nq_a2))*0.5d0/
     '    ((LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)
     '    *0.25d0)**0.5d0
        ROIB1h=(CQ(4,nq_b1)+CQ(4,nq_b2))*0.5d0/
     '    ((LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
     '    *0.25d0)**0.5d0
        ROIC1h=(CQ(4,nq_c1)+CQ(4,nq_c2))*0.5d0/
     '    ((LAMBDA_C1+LAMBDA_C2+OLD_LAMBDA_C1+OLD_LAMBDA_C2)
     '    *0.25d0)**0.5d0
C unstressed radii at half time step i+1/2 and k+1/2
        ROIA1=CQ(4,nq_a1)/(LAMBDA_A1**0.5d0)
        ROIB1=CQ(4,nq_b1)/(LAMBDA_B1**0.5d0)
        ROIC1=CQ(4,nq_c1)/(LAMBDA_C1**0.5d0)
C unstressed radii at full time step i+1 and k+1
        BETAA1h=(CQ(6,nq_a1)+CQ(6,nq_a2))*0.5d0*
     '    ((LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)
     '    *0.25d0)+((CQ(12,nq_a1)+CQ(12,nq_a2))*0.5d0)
        BETAB1h=(CQ(6,nq_b1)+CQ(6,nq_b2))*0.5d0*
     '    ((LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
     '    *0.25d0)+((CQ(12,nq_b1)+CQ(12,nq_b2))*0.5d0)
        BETAC1h=(CQ(6,nq_c1)+CQ(6,nq_c2))*0.5d0*
     '    ((LAMBDA_C1+LAMBDA_C2+OLD_LAMBDA_C1+OLD_LAMBDA_C2)
     '    *0.25d0)+((CQ(12,nq_c1)+CQ(12,nq_c2))*0.5d0)
        BETAA1h=DMAX1(BETAA1h,1.0d0)
        BETAB1h=DMAX1(BETAB1h,1.0d0)
        BETAC1h=DMAX1(BETAC1h,1.0d0)
C calculate the pressure radius exponent  at the half time
C step i+1/2 and k+1/2
        BETAA1=CQ(6,nq_a1)*LAMBDA_A1+CQ(12,nq_a1)
        BETAB1=CQ(6,nq_b1)*LAMBDA_B1+CQ(12,nq_b1)
        BETAC1=CQ(6,nq_c1)*LAMBDA_C1+CQ(12,nq_c1)
        BETAA1=DMAX1(BETAA1,1.0d0)
        BETAB1=DMAX1(BETAB1,1.0d0)
        BETAC1=DMAX1(BETAC1,1.0d0)
C calculate the pressure radius exponent  at the half time
C step i+1 and k+1
        DELTA_XA=0.0d0
        DO nj=1,NJTOT !calculate delta x for half time step
          DELTA_XA=DELTA_XA+((XQ(nj,nq_a1)-XQ(nj,nq_a2))**2.0d0)
        ENDDO
C       DELTA_XA=DMAX1((DELTA_XA**0.5d0),0.6d0)*
C       '      (LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)*0.25d0
        DELTA_XA=(DELTA_XA**0.5d0)*
     '    (LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)*0.25d0
        DELTA_XB=0.0d0
        DO nj=1,NJTOT !calculate delta x for half time step
          DELTA_XB=DELTA_XB+((XQ(nj,nq_b1)-XQ(nj,nq_b2))**2.0d0)
        ENDDO
C       DELTA_XB=DMAX1((DELTA_XB**0.5d0),0.6d0)*
C       '      (LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)*0.25d0
        DELTA_XB=(DELTA_XB**0.5d0)*
     '    (LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)*0.25d0
        DELTA_XC=0.0d0
        DO nj=1,NJTOT !calculate delta x for half time step
          DELTA_XC=DELTA_XC+((XQ(nj,nq_c1)-XQ(nj,nq_c2))**2.0d0)
        ENDDO
C       DELTA_XC=DMAX1((DELTA_XC**0.5d0),0.6d0)*
C       '      (LAMBDA_C1+LAMBDA_C2+OLD_LAMBDA_C1+OLD_LAMBDA_C2)*0.25d0
        DELTA_XC=(DELTA_XC**0.5d0)*
     '    (LAMBDA_C1+LAMBDA_C2+OLD_LAMBDA_C1+OLD_LAMBDA_C2)*0.25d0
C calculate the delta x space steps between the first and second points
C e.g. a1 and a2, surrounding a bifurcation
        OLDP_INITIAL=YQ(NYNQ(1,nq,0),6)
C calculating Po the pressre at centre of the bifurcation
        LA=1.0d-9/((ROIA1/CQ(4,1))**2.0d0)
        LB=1.0d-9/((ROIB1/CQ(4,1))**2.0d0)
        LC=1.0d-9/((ROIC1/CQ(4,1))**2.0d0)
c setting the relative inductance values for the bifurcation junction
        FA2=PI*(YQ(ny_a2r,1)**2.0d0)*Va2_new
        OLDFA1=PI*(YQ(ny_a1r,8)**2.0d0)*Va1_old
        OLDFA2=PI*(YQ(ny_a2r,8)**2.0d0)*Va2_old
        FB2=PI*(YQ(ny_b2r,1)**2.0d0)*Vb2_new
        OLDFB1=PI*(YQ(ny_b1r,8)**2.0d0)*Vb1_old
        OLDFB2=PI*(YQ(ny_b2r,8)**2.0d0)*Vb2_old
        FC2=PI*(YQ(ny_c2r,1)**2.0d0)*Vc2_new
        OLDFC1=PI*(YQ(ny_c1r,8)**2.0d0)*Vc1_old
        OLDFC2=PI*(YQ(ny_c2r,8)**2.0d0)*Vc2_old
C calculating the flows at the bifurcation grid points
        Aa=((CQ(5,nq_a1)+CQ(5,nq_a2))*BETAA1h)
     '    *((YQ(ny_a3r,5)/ROIA1h)
     '    **(BETAA1h-1.0d0))
     '    /(PI*YQ(ny_a3r,5)*ROIA1h)
        Ab=((CQ(5,nq_b1)+CQ(5,nq_b2))*BETAB1h)
     '    *((2.0d0*YQ(ny_b3r,5)/ROIB1h)
     '    **(BETAB1h-1.0d0))
     '    /(PI*YQ(ny_b3r,5)*ROIB1h)
        Ac=((CQ(5,nq_c1)+CQ(5,nq_c2))*BETAC1h)
     '    *((2.0d0*YQ(ny_c3r,5)/ROIC1h)
     '    **(BETAC1h-1.0d0))
     '    /(PI*YQ(ny_c3r,5)*ROIC1h)
        PRESTERM=YQ(ny_a1p,8)+OLD_TRACE_A1
     '    +YQ(ny_a2p,8)+OLD_TRACE_A2-YQ(ny_a2p,1)-TRACE_A2
     '    -TRACE_A1
        DENOM=(1.0d0+(Aa*(TINCR**2.0d0)/(2.0d0*DELTA_XA*LA)))
        TERMREST=(Aa*(TINCR/DELTA_XA))*(FA2+OLDFA2-(2.0d0*OLDFA1)-
     '    ((TINCR/(2.0d0*LA))*(YQ(ny_a1p,8)-OLDP_INITIAL)))
        LAM_TERMA=(PI*Aa*(YQ(ny_a3r,5)**2.0d0)/LAMBDA_AH)*
     '    (LAMBDA_A1+LAMBDA_A2-OLD_LAMBDA_A1-OLD_LAMBDA_A2)
        H=(PRESTERM+TERMREST-LAM_TERMA)/DENOM
        I=(DENOM-1)/DENOM
        PRESTERM=YQ(ny_b1p,8)+YQ(ny_b2p,8)-YQ(ny_b2p,1)
     '    +OLD_TRACE_B1+OLD_TRACE_B2-TRACE_B1-TRACE_B2
        DENOM=(1.0d0+(Ab*(TINCR**2.0d0)/(2.0d0*DELTA_XB*LB)))
        TERMREST=(Ab*(TINCR/DELTA_XB))*(FB2+OLDFB2-(2.0d0*OLDFB1)-
     '    ((TINCR/(2.0d0*LB))*(OLDP_INITIAL-YQ(ny_b1p,8))))
        LAM_TERMB=(PI*Ab*(YQ(ny_b3r,5)**2.0d0)/LAMBDA_BH)*
     '    (LAMBDA_B1+LAMBDA_B2-OLD_LAMBDA_B1-OLD_LAMBDA_B2)
        J=(PRESTERM-TERMREST-LAM_TERMB)/DENOM
        K=(DENOM-1)/DENOM
        PRESTERM=YQ(ny_c1p,8)+YQ(ny_c2p,8)-YQ(ny_c2p,1)
     '    +OLD_TRACE_C1+OLD_TRACE_C2-TRACE_C1-TRACE_C2
        DENOM=(1.0d0+(Ac*(TINCR**2.0d0)/(2.0d0*DELTA_XC*LC)))
        TERMREST=(Ac*(TINCR/DELTA_XC))*(FC2+OLDFC2-(2.0d0*OLDFC1)-
     '    ((TINCR/(2.0d0*LC))*(OLDP_INITIAL-YQ(ny_c1p,8))))
        LAM_TERMC=(PI*Ac*(YQ(ny_c3r,5)**2.0d0)/LAMBDA_CH)*
     '    (LAMBDA_C1+LAMBDA_C2-OLD_LAMBDA_C1-OLD_LAMBDA_C2)
        L=(PRESTERM-TERMREST-LAM_TERMC)/DENOM
        M=(DENOM-1)/DENOM
        DENOM=((1.0d0/LA)+(1.0d0/LB)+(1.0d0/LC))
        DENOM1=1.0d0-(I/(LA*DENOM))-(K/(LB*DENOM))-(M/(LC*DENOM))
        
C!!! initial pressure guess (KSB) Equation 4.7 (Nic's thesis)        
        P_INITIAL=(-OLDP_INITIAL+((((H+YQ(ny_a1p,8)+G_TERMA*LA)
     '    /LA)+((J+YQ(ny_b1p,8)-G_TERMB*LB)/LB)+((L+
     '    YQ(ny_c1p,8)-G_TERMC*LC)/LC))/DENOM))/DENOM1
        YQ(NYNQ(1,nq,0),6)=P_INITIAL
C calculated the initial guess for the pressure at the centre of the
C bifurcation
        PRESTERM=YQ(ny_a1p,8)+YQ(ny_a2p,8)-YQ(ny_a2p,1)
     '    +OLD_TRACE_A1+OLD_TRACE_A2-TRACE_A1-TRACE_A2
        DENOM=(1.0d0+(Aa*(TINCR**2.0d0)/(2.0d0*DELTA_XA*LA)))
        TERMREST=(AA*(TINCR/DELTA_XA))*(FA2+OLDFA2-(2.d0*OLDFA1)-
     '    ((TINCR/(2.0d0*LA))*(YQ(ny_a1p,8)-
     '    P_INITIAL-OLDP_INITIAL)))
        YQ(ny_a1p,1)=(PRESTERM+TERMREST-LAM_TERMA)/DENOM
        PRESTERM=YQ(ny_b1p,8)+YQ(ny_b2p,8)-YQ(ny_b2p,1)
     '    +OLD_TRACE_B1+OLD_TRACE_B2-TRACE_B1-TRACE_B2
        DENOM=(1.0d0+(Ab*(TINCR**2.0d0)/(2.0d0*DELTA_XB*LB)))
        TERMREST=(Ab*(TINCR/DELTA_XB))*(FB2+OLDFB2-(2.0d0*OLDFB1)-
     '    ((TINCR/(2.0d0*LB))*(P_INITIAL+OLDP_INITIAL-YQ(ny_b1p,8))))
        YQ(ny_b1p,1)=(PRESTERM-TERMREST-LAM_TERMB)/DENOM
        PRESTERM=YQ(ny_c1p,8)+YQ(ny_c2p,8)-YQ(ny_c2p,1)
     '    +OLD_TRACE_C1+OLD_TRACE_C2-TRACE_C1-TRACE_C2
        DENOM=(1.0d0+(Ac*(TINCR**2.0d0)/(2.0d0*DELTA_XC*LC)))
        TERMREST=(Ac*(TINCR/DELTA_XC))*(FC2+OLDFC2-(2.0d0*OLDFC1)-
     '    ((TINCR/(2.0d0*LC))*(P_INITIAL+OLDP_INITIAL-YQ(ny_c1p,8))))
        YQ(ny_c1p,1)=(PRESTERM-TERMREST-LAM_TERMC)/DENOM
C calculates the initial guess of the bifurcation pressures
C at grid points a1, b1, and c1
        IF((YQ(ny_a1p,1)+TRACE_A1).GT.0.0d0) THEN
          YQ(ny_a1r,1)=ROIA1*((((YQ(ny_a1p,1)+TRACE_A1)
     '      /CQ(5,nq_a1))+1.0d0)**(1.0d0/BETAA1))
        ELSE
          fo=BETAA1*CQ(5,nq_a1)/CQ(7,nq_a1)
          YQ(ny_a1r,1)=(ROIA1/((1-((YQ(ny_a1p,1)+TRACE_A1)/
     '      fo))**(1.0d0/CQ(7,nq_a1))))
        ENDIF
        IF((YQ(ny_b1p,1)+TRACE_B1).GT.0.0d0) THEN
          YQ(ny_b1r,1)=ROIB1*((((YQ(ny_b1p,1)+TRACE_B1)
     '      /CQ(5,nq_b1))+1.0d0)**(1.0d0/BETAB1))
        ELSE
          fo=BETAB1*CQ(5,nq_b1)/CQ(7,nq_b1)
          YQ(ny_b1r,1)=(ROIB1/((1-((YQ(ny_b1p,1)+TRACE_B1)/
     '      fo))**(1.0d0/CQ(7,nq_b1))))
        ENDIF
        IF((YQ(ny_c1p,1)+TRACE_C1).GT.0.0d0) THEN
          YQ(ny_c1r,1)=ROIC1*((((YQ(ny_c1p,1)+TRACE_C1)
     '      /CQ(5,nq_c1))+1.0d0)**(1.0d0/BETAC1))
        ELSE
          fo=BETAC1*CQ(5,nq_c1)/CQ(7,nq_c1)
          YQ(ny_c1r,1)=(ROIC1/((1-((YQ(ny_c1p,1)+TRACE_C1)/
     '      fo))**(1.0d0/CQ(7,nq_c1))))
        ENDIF
C calculates the inital radius values at the bifurcation grid points
        TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '    ((Va3_mid**2.0d0)/YQ(ny_a3r,5))*
     '    (YQ(ny_a1r,1)+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))        
        P_DIFF=YQ(ny_a1p,1)+YQ(ny_a1p,8)-YQ(ny_a2p,1)-YQ(ny_a2p,8)
     '    -2.d0*G_TERMA*DELTA_XA
        TERM2=(1.0D0/CQ(1,nq_a1))*P_DIFF            
        IF(YQ(ny_a3r,5).GT.
     '    ROIA1h) THEN
          TERM3=2.0d0*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)/(CQ(3,nq_a1)
     &      -1.0d0))*(Va3_mid/(YQ(ny_a3r,5)**2.0d0))
        ELSE
          RO_ART=ROIA1h
          ELIP_A=((RO_ART**2.0d0)+
     '      (RO_ART**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '      **0.5d0)**0.5d0
          ELIP_B=((RO_ART**2.0d0)-
     '      (RO_ART**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '      **0.5d0)**0.5d0
          ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '      (YQ(ny_a3r,5)**2.0d0*(ELIP_A**2.0d0+
     '      ELIP_B**2.0d0)))**0.5D0
          TERM3=2.0d0*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)/(CQ(3,nq_a1)
     &      -1.0d0))*(Va3_mid/(ELIP_TERM**2.0d0))
        ENDIF
        TERM4=(1.0d0/YQ(ny_a3r,5))*((-DELTA_XA/TINCR)+
     '    ((2.0d0*CQ(3,nq_a1)-1.0d0)*Va3_mid))*
     '    ((YQ(ny_a2r,1)*(LAMBDA_A2**0.5D0))+
     '    (YQ(ny_a1r,1)*(LAMBDA_A1**0.5D0))-
     '    (YQ(ny_a2r,8)*(OLD_LAMBDA_A2**0.5D0))-
     '    (YQ(ny_a1r,8)*(OLD_LAMBDA_A1**0.5D0)))
        TERM5=Va3_mid/YQ(ny_a3r,5)
     '    *(YQ(ny_a1r,1)+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
        Va1_new=Va2_old+((TINCR/(2.0d0*DELTA_XA))*
     '    (TERM1-TERM2))-TERM3+TERM4-TERM5
        TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '    ((Vb3_mid**2.0d0)/YQ(ny_b3r,5))*
     '    (YQ(ny_b2r,1)+YQ(ny_b2r,8)-YQ(ny_b1r,1)-YQ(ny_b1r,8))
        P_DIFF=YQ(ny_b2p,1)+YQ(ny_b2p,8)-YQ(ny_b1p,1)-YQ(ny_b1p,8)
     '    +2.d0*G_TERMB*DELTA_XB          
        TERM2=(1.0D0/CQ(1,nq_b1))*P_DIFF         
        IF(YQ(ny_b3r,5).GT.ROIB1h) THEN
          TERM3=2.0d0*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)/(CQ(3,nq_b1)
     &      -1.0d0))*(Vb3_mid/(YQ(ny_b3r,5)**2.0d0))
        ELSE
          RO_ART=ROIB1h
          ELIP_A=((RO_ART**2.0d0)+
     '      (RO_ART**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '      **0.5d0)**0.5d0
          ELIP_B=((RO_ART**2.0d0)-
     '      (RO_ART**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '      **0.5d0)**0.5d0
          ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '      (YQ(ny_b3r,5)**2.0d0*(ELIP_A**2.0d0+
     '      ELIP_B**2.0d0)))**0.5D0
          TERM3=2.0d0*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)/(CQ(3,nq_b1)
     &      -1.0d0))*(Vb3_mid/(ELIP_TERM**2.0d0))
        ENDIF
        TERM4=(1.0d0/YQ(ny_b3r,5))*((DELTA_XB/TINCR)+
     '    ((2.0d0*CQ(3,nq_b1)-1.0d0)*Vb3_mid))*
     '    (((LAMBDA_B1**0.5d0)*YQ(ny_b1r,1))+
     '    ((LAMBDA_B2**0.5d0)*YQ(ny_b2r,1))-
     '    ((OLD_LAMBDA_B2**0.5d0)*YQ(ny_b2r,8))-
     '    ((OLD_LAMBDA_B1**0.5d0)*YQ(ny_b1r,8)))
        TERM5=Vb3_mid/YQ(ny_b3r,5)
     '    *(YQ(ny_b2r,1)+YQ(ny_b2r,8)-YQ(ny_b1r,1)-YQ(ny_b1r,8))
        Vb1_new=Vb2_old+((TINCR/(2.0d0*DELTA_XB))*
     '    (TERM1-TERM2))-TERM3+TERM4+TERM5
        TERM1=((2.0d0*CQ(3,nq_c1))-0.0d0)*
     '    ((Vc3_mid**2.0d0)/YQ(ny_c3r,5))*
     '    (YQ(ny_c2r,1)+YQ(ny_c2r,8)-YQ(ny_c1r,1)-YQ(ny_c1r,8))
        P_DIFF=YQ(ny_c2p,1)+YQ(ny_c2p,8)-YQ(ny_c1p,1)-YQ(ny_c1p,8)
     '    +2.d0*G_TERMC*DELTA_XC
        TERM2=(1.0D0/CQ(1,nq_c1))*P_DIFF
        IF(YQ(ny_c3r,5).GT.ROIC1h) THEN
          TERM3=2.0d0*TINCR*CQ(2,nq_c1)*(CQ(3,nq_c1)/(CQ(3,nq_c1)
     &      -1.0d0))*(Vc3_mid/(YQ(ny_c3r,5)**2.0d0))
        ELSE
          RO_ART=ROIC1h
          ELIP_A=((RO_ART**2.0d0)+
     '      (RO_ART**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '      **0.5d0)**0.5d0
          ELIP_B=((RO_ART**2.0d0)-
     '      (RO_ART**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '      **0.5d0)**0.5d0
          ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '      (YQ(ny_c3r,5)**2.0d0*(ELIP_A**2.0d0+
     '      ELIP_B**2.0d0)))**0.5D0
          TERM3=2.0d0*TINCR*CQ(2,nq_c1)*(CQ(3,nq_c1)/(CQ(3,nq_c1)
     &      -1.0d0))*(Vc3_mid/(ELIP_TERM**2.0d0))
        ENDIF
        TERM4=(1.0d0/YQ(ny_c3r,5))*((DELTA_XC/TINCR)+
     '    ((2.0d0*CQ(3,nq_c1)-1.0d0)*Vc3_mid))*
     '    (((LAMBDA_C1**0.5d0)*YQ(ny_c1r,1))+
     '    ((LAMBDA_C2**0.5d0)*YQ(ny_c2r,1))-
     '    ((OLD_LAMBDA_C2**0.5d0)*YQ(ny_c2r,8))-
     '    ((OLD_LAMBDA_C1**0.5d0)*YQ(ny_c1r,8)))
        TERM5=Vc3_mid/YQ(ny_c3r,5)
     '    *(YQ(ny_c2r,1)+YQ(ny_c2r,8)-YQ(ny_c1r,1)-YQ(ny_c1r,8))
        Vc1_new=Vc2_old+((TINCR/(2.0d0*DELTA_XC))*
     '    (TERM1-TERM2))-TERM3+TERM4+TERM5
C calculates the initial guess for the velocity values at grid points
C surrounding a bifurcation. To account for the increased resistance
C with collapse these terms for each point are calculated using a
C modified term. Calculate flows in 3 segments
        FA1=PI*(YQ(ny_a1r,1)**2.0d0)*Va1_new
        FB1=PI*(YQ(ny_b1r,1)**2.0d0)*Vb1_new
        FC1=PI*(YQ(ny_c1r,1)**2.0d0)*Vc1_new
C Newton-Raphson loop for arterial bifurcation conservation eqns
        LOOP_DENOM=((1.0d0/LA)+(1.0d0/LB)+(1.0d0/LC))
        PA=YQ(ny_a1p,1)
        PB=YQ(ny_B1p,1)
        PC=YQ(ny_C1p,1)
        SUM=1.0D0
        NEW_COUNT=0
        DO WHILE(((SUM.GT.(CONVERG_TOL*1000.0d0)).AND.
     '    ((dabs(FA1-FB1-FC1)).GT.(CONVERG_TOL*1000.0d0))
     '    .OR.(NEW_COUNT.LT.5))
     '    .AND.(NEW_COUNT.LT.50))
          NEW_COUNT=NEW_COUNT+1
          IF((PA+TRACE_A1).GT.0.0d0) THEN
            Sa=PI*(ROIA1**2.0d0)*
     '        ((((PA+TRACE_A1)/CQ(5,nq_a1))+1)
     '        **(2.0D0/BETAA1))
            dSa=(PI*(ROIA1**2.0d0)*2.0d0/
     '        (BETAA1*CQ(5,nq_a1))
     '        )*((((PA+trace_a1)/CQ(5,nq_a1))+1)
     '        **((2.0D0/BETAA1)-1.0d0))
            Ra=ROIA1*((((PA+TRACE_A1)/CQ(5,nq_a1))+1)
     '        **(1.0D0/BETAA1))
            dRa=(ROIA1/
     '        (BETAA1*CQ(5,nq_a1)))
     '        *((((PA+trace_a1)/CQ(5,nq_a1))+1)
     '        **((1.0D0/BETAA1)-1.0d0))
          ELSE
            fo=BETAA1*CQ(5,nq_a1)/CQ(7,nq_a1)
            Sa=PI*((ROIA1/((1-((PA+TRACE_A1)/
     '        fo))**(1.0d0/CQ(7,nq_a1))))**2.0d0)
            dSa=(2.0D0*PI*(ROIA1**2.0d0)/(CQ(7,nq_a1)*fo))
     '        *(1-(PA+TRACE_A1)/Fo)**(-2.0d0/CQ(7,nq_a1)-1.0d0)
            Ra=(ROIA1/((1-((PA+TRACE_A1)/
     '        fo))**(1.0d0/CQ(7,nq_a1))))
            dRa=(ROIA1/(CQ(7,nq_a1)*fo))
     '        *(1-(PA+TRACE_A1)/Fo)**(-1.0d0/CQ(7,nq_a1)-1.0d0)
          ENDIF !(PA+TRACE_A1).GT.0.0d0
C calculates the radius, area, and derivitives with respect to pressure
C at the "a" grid point adjacent to a bifurcation
          IF((PB+TRACE_B1).GT.0.0d0) THEN
            Sb=PI*(ROIB1**2.0d0)*
     '        ((((PB+TRACE_B1)/CQ(5,nq_b1))+1)
     '        **(2.0D0/BETAB1))
            dSb=(PI*(ROIB1**2.0d0)*2.0d0/
     '        (BETAB1*CQ(5,nq_b1))
     '        )*((((Pb+trace_b1)/CQ(5,nq_b1))+1)
     '        **((2.0D0/BETAB1)-1.0d0))
            Rb=ROIB1*((((PB+TRACE_B1)/CQ(5,nq_b1))+1)
     '        **(1.0D0/BETAB1))
            dRb=(ROIB1/
     '        (BETAB1*CQ(5,nq_b1)))
     '        *((((Pb+trace_b1)/CQ(5,nq_b1))+1)
     '        **((1.0D0/BETAB1)-1.0d0))
          ELSE
            fo=BETAB1*CQ(5,nq_b1)/CQ(7,nq_b1)
            Sb=PI*((ROIB1/((1-((PB+TRACE_B1)/
     '        fo))**(1.0d0/CQ(7,nq_b1))))**2.0d0)
            dSb=(2.0D0*PI*(ROIB1**2.0d0)/(CQ(7,nq_b1)*fo))
     '        *(1-(PB+TRACE_B1)/Fo)**(-2.0d0/CQ(7,nq_b1)-1.0d0)
            Rb=(ROIB1/((1-((PB+TRACE_B1)/
     '        fo))**(1.0d0/CQ(7,nq_b1))))
            dRb=(ROIB1/(CQ(7,nq_b1)*fo))
     '        *(1-(PB+TRACE_B1)/Fo)**(-1.0d0/CQ(7,nq_b1)-1.0d0)
          ENDIF !(PB+TRACE_B1).GT.0.0d0
C calculates the radius, area, and derivitives with respect to pressure
C at the "b" grid point adjacent to a bifurcation
          IF((PC+TRACE_C1).GT.0.0d0) THEN
            Sc=PI*(ROIC1**2.0d0)*
     '        ((((PC+TRACE_C1)/CQ(5,nq_c1))+1)
     '        **(2.0D0/BETAC1))
            dSc=(PI*(ROIC1**2.0d0)*2.0d0/
     '        (BETAC1*CQ(5,nq_c1))
     '        )*((((PC+trace_c1)/CQ(5,nq_c1))+1)
     '        **((2.0D0/BETAC1)-1.0d0))
            Rc=ROIC1*((((PC+TRACE_C1)/CQ(5,nq_c1))+1)
     '        **(1.0D0/BETAC1))
            dRc=(ROIC1/
     '        (BETAC1*CQ(5,nq_c1)))
     '        *((((PC+trace_c1)/CQ(5,nq_c1))+1)
     '        **((1.0D0/BETAC1)-1.0d0))
          ELSE
            fo=BETAC1*CQ(5,nq_c1)/CQ(7,nq_c1)
            Sc=PI*((ROIC1/((1-((PC+TRACE_C1)/
     '        fo))**(1.0d0/CQ(7,nq_c1))))**2.0d0)
            dSc=(2.0D0*PI*(ROIC1**2.0d0)/(CQ(7,nq_c1)*fo))
     '        *(1-(PC+TRACE_C1)/Fo)**(-2.0d0/CQ(7,nq_c1)-1.0d0)
            Rc=(ROIC1/((1-((PC+TRACE_C1)/
     '        fo))**(1.0d0/CQ(7,nq_c1))))
            dRc=(ROIC1/(CQ(7,nq_c1)*fo))
     '        *(1-(PC+TRACE_C1)/Fo)**(-1.0d0/CQ(7,nq_c1)-1.0d0)
          ENDIF !(PC+TRACE_C1).GT.0.0d0
C calculates the radius, area, and derivitives with respect to pressure
C at the "c" grid point adjacent to a bifurcation
          TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '      ((Va3_mid**2.0d0)/YQ(ny_a3r,5))*
     '      (RA+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
          P_DIFF=PA+YQ(ny_a1p,8)-YQ(ny_a2p,1)-YQ(ny_a2p,8)-2.d0
     '      *G_TERMA*DELTA_XA
          TERM2=(1.0D0/CQ(1,nq_a1))*P_DIFF
          IF(YQ(ny_a3r,5).GT.ROIA1h) THEN
            TERM3=2.0d0*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)
     '        /(CQ(3,nq_a1)-1.0d0))*
     '        (Va3_mid/(YQ(ny_a3r,5)**2.0d0))
          ELSE
            RO_ART=ROIA1h
            ELIP_A=((RO_ART**2.0d0)+
     '        (RO_ART**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '        **0.5d0)**0.5d0
            ELIP_B=((RO_ART**2.0d0)-
     '        (RO_ART**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '        **0.5d0)**0.5d0
            ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '        (YQ(ny_a3r,5)**2.0d0*(ELIP_A**2.0d0+
     '        ELIP_B**2.0d0)))**0.5D0
            TERM3=2.0d0*TINCR*CQ(2,nq_a1)*
     '        (CQ(3,nq_a1)/(CQ(3,nq_a1)-1.0d0))*
     '        (Va3_mid/(ELIP_TERM**2.0d0))
          ENDIF
          TERM4=(1.0d0/YQ(ny_a3r,5))*((-DELTA_XA/TINCR)+
     '      ((2.0d0*CQ(3,nq_a1)-1.0d0)*Va3_mid))*
     '      ((YQ(ny_a2r,1)*(LAMBDA_A2**0.5D0))+
     '      (RA*(LAMBDA_A1**0.5D0))-
     '      (YQ(ny_a2r,8)*(OLD_LAMBDA_A2**0.5D0))-
     '      (YQ(ny_a1r,8)*(OLD_LAMBDA_A1**0.5D0)))
            TERM5=Va3_mid/YQ(ny_a3r,5)
     '      *(RA+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
            VA=Va2_old+((TINCR/(2.0d0*DELTA_XA))*
     '        (TERM1-TERM2))-TERM3+TERM4-TERM5
C uses the end point bounary condtions to calculate velocity at grid
C point "a" and the derivative with respect to pressure
            TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '        ((Va3_mid**2.0d0)/YQ(ny_a3r,5))*dRa
            TERM2=1.0d0/CQ(1,nq_a1)
            TERM3=((LAMBDA_A1**0.5d0)/YQ(ny_a3r,5))*((-DELTA_XA/TINCR)+
     '        ((2.0d0*CQ(3,nq_a1)-1.0d0)*Va3_mid))*dRa
            TERM4=Va3_mid/YQ(ny_a3r,5)*dRa
            DVA=((TINCR/(2.0d0*DELTA_XA))*
     '        (TERM1-TERM2))+TERM3-TERM4
            TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '        ((Vb3_mid**2.0d0)/YQ(ny_b3r,5))*
     '        (YQ(ny_b2r,1)+YQ(ny_b2r,8)-Rb-YQ(ny_b1r,8))
            P_DIFF=YQ(ny_b2p,1)+YQ(ny_b2p,8)-PB-YQ(ny_b1p,8)+2.d0
     '        *G_TERMB*DELTA_XB
            TERM2=(1.0D0/CQ(1,nq_b1))*P_DIFF
            IF(YQ(ny_b3r,5).GT.ROIB1h) THEN
              TERM3=2.0d0*TINCR*CQ(2,nq_b1)*
     '          (CQ(3,nq_b1)/(CQ(3,nq_b1)-1.0d0))*
     '          (Vb3_mid/(YQ(ny_b3r,5)**2.0d0))
            ELSE
              RO_ART=ROIB1h
              ELIP_A=((RO_ART**2.0d0)+
     '          (RO_ART**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_B=((RO_ART**2.0d0)-
     '          (RO_ART**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '          (YQ(ny_b3r,5)**2.0d0*(ELIP_A**2.0d0+
     '          ELIP_B**2.0d0)))**0.5D0
              TERM3=2.0d0*TINCR*CQ(2,nq_b1)*
     '          (CQ(3,nq_b1)/(CQ(3,nq_b1)-1.0d0))*
     '          (Vb3_mid/(ELIP_TERM**2.0d0))
            ENDIF
            TERM4=(1.0d0/YQ(ny_b3r,5))*((DELTA_XB/TINCR)+
     '        ((2.0d0*CQ(3,nq_b1)-1.0d0)*Vb3_mid))*
     '        (((LAMBDA_B1**0.5d0)*Rb)+
     '        ((LAMBDA_B2**0.5d0)*YQ(ny_b2r,1))-
     '        ((OLD_LAMBDA_B2**0.5d0)*YQ(ny_b2r,8))-
     '        ((OLD_LAMBDA_B1**0.5d0)*YQ(ny_b1r,8)))
            TERM5=Vb3_mid/YQ(ny_b3r,5)
     '        *(YQ(ny_b2r,1)+YQ(ny_b2r,8)-Rb-YQ(ny_b1r,8))
            VB=Vb2_old+((TINCR/(2.0d0*DELTA_XB))*
     '        (TERM1-TERM2))-TERM3+TERM4+TERM5
C uses the end point bounary condtions to calculate velocity at grid
C point "b" and the derivative with respect to pressure
            TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '        ((Vb3_mid**2.0d0)/YQ(ny_b3r,5))*(-dRb)
            TERM2=-1.0d0/CQ(1,nq_b1)
            TERM3=(1.0d0/YQ(ny_b3r,5))*((DELTA_XB/TINCR)+
     '        ((2.0d0*CQ(3,nq_b1)-1.0d0)*Vb3_mid))*
     '        (LAMBDA_B1**0.5d0)*dRb
            TERM4=Vb3_mid/YQ(ny_b3r,5)*dRb
            DVB=((TINCR/(2.0d0*DELTA_XB))*
     '        (TERM1-TERM2))+TERM3-TERM4
            TERM1=((2.0d0*CQ(3,nq_c1))-0.0d0)*
     '        ((Vc3_mid**2.0d0)/YQ(ny_c3r,5))*
     '        (YQ(ny_c2r,1)+YQ(ny_c2r,8)-Rc-YQ(ny_c1r,8))
            P_DIFF=YQ(ny_c2p,1)+YQ(ny_c2p,8)-PC-YQ(ny_c1p,8)+2.d0
     '        *G_TERMC*DELTA_XC
            TERM2=(1.0D0/CQ(1,nq_c1))*P_DIFF
            IF(YQ(ny_c3r,5).GT.
     '        ROIC1h) THEN
              TERM3=2.0d0*TINCR*CQ(2,nq_c1)*(CQ(3,nq_c1)
     '          /(CQ(3,nq_c1)-1.0d0))*
     '          (Vc3_mid/(YQ(ny_c3r,5)**2.0d0))
            ELSE
              RO_ART=ROIC1h
              ELIP_A=((RO_ART**2.0d0)+
     '          (RO_ART**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_B=((RO_ART**2.0d0)-
     '          (RO_ART**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '          (YQ(ny_c3r,5)**2.0d0*(ELIP_A**2.0d0+
     '          ELIP_B**2.0d0)))**0.5D0
              TERM3=2.0d0*TINCR*CQ(2,nq_c1)*
     '          (CQ(3,nq_c1)/(CQ(3,nq_c1)-1.0d0))*
     '          (Vc3_mid/(ELIP_TERM**2.0d0))
            ENDIF
            TERM4=(1.0d0/YQ(ny_c3r,5))*((DELTA_XC/TINCR)+
     '        ((2.0d0*CQ(3,nq_c1)-1.0d0)*Vc3_mid))*
     '        (((LAMBDA_C1**0.5d0)*Rc)+
     '        ((LAMBDA_C2**0.5d0)*YQ(ny_c2r,1))-
     '        ((OLD_LAMBDA_C2**0.5d0)*YQ(ny_c2r,8))-
     '        ((OLD_LAMBDA_C1**0.5d0)*YQ(ny_c1r,8)))
            TERM5=Vc3_mid/YQ(ny_c3r,5)
     '        *(YQ(ny_c2r,1)+YQ(ny_c2r,8)-Rc-YQ(ny_c1r,8))
            VC=Vc2_old+((TINCR/(2.0d0*DELTA_XC))*
     '        (TERM1-TERM2))-TERM3+TERM4+TERM5
C uses the end point bounary condtions to calculate velocity at grid
C point "c" and the derivative with respect to pressure
            TERM1=((2.0d0*CQ(3,nq_c1))-0.0d0)*
     '        ((Vc3_mid**2.0d0)/YQ(ny_c3r,5))*(-dRc)
            TERM2=-1.0d0/CQ(1,nq_c1)
            TERM3=(1.0d0/YQ(ny_c3r,5))*((DELTA_XC/TINCR)+
     '        ((2.0d0*CQ(3,nq_c1)-1.0d0)*Vc3_mid))*
     '        (LAMBDA_C1**0.5d0)*dRc
            TERM4=Vc3_mid/YQ(ny_c3r,5)*dRc
            DVC=((TINCR/(2.0d0*DELTA_XC))*
     '        (TERM1-TERM2))+TERM3-TERM4
            dFA=(Sa*DVA)+(dSa*VA)
            dFb=(Sb*DVb)+(dSb*Vb)
            dFc=(Sc*DVc)+(dSc*Vc)
            FA1=PI*(Ra**2.0d0)*va
            FB1=PI*(Rb**2.0d0)*vb
            FC1=PI*(Rc**2.0d0)*vc
C calculates flow and their derivative with respect to pressure in
C each of the bifurcation segments
            FA1_OLD=PI*(YQ(ny_a1r,8)**2.0d0)*Va1_old
            FB1_OLD=PI*(YQ(ny_b1r,8)**2.0d0)*Vb1_old
            FC1_OLD=PI*(YQ(ny_c1r,8)**2.0d0)*Vc1_old
            DG(1,1)=1-((2.0d0*LA/TINCR)*dFA)
            DG(1,2)=-1-((2.0d0*LB/TINCR)*dFB)
            DG(1,3)=0.0D0
            DG(2,1)=DG(1,1)
            DG(2,2)=0.0D0
            DG(2,3)=-1-((2.0d0*LC/TINCR)*dFC)
            DG(3,1)=dFA
            DG(3,2)=-dFB
            DG(3,3)=-dFC

C... Equation 4.17 - added gravity term (08/04)
            X(1)=((2.0d0*LA/TINCR)*(FA1-FA1_OLD))+
     '        ((2.0d0*Lb/TINCR)*(FB1-FB1_OLD))+
     '        PB+YQ(ny_b1p,8)-PA-YQ(ny_a1p,8)-2.d0*LA*G_TERMA
     '        -2.d0*LB*G_TERMB
            X(2)=((2.0d0*LA/TINCR)*(FA1-FA1_OLD))+
     '        ((2.0d0*Lc/TINCR)*(FC1-FC1_OLD))+
     '        PC+YQ(ny_C1p,8)-PA-YQ(ny_a1p,8)-2.d0*LA*G_TERMA
     '        -2.d0*LC*G_TERMC
            X(3)=FC1+FB1-FA1
            B(1)=X(1)
            B(2)=X(2)
            B(3)=X(3)
            X(1)=(b(3)-(b(1)*DG(3,2)/DG(1,2))-(b(2)*DG(3,3)/DG(2,3)))/
     '        (DG(3,1)-(DG(1,1)*DG(3,3)/DG(2,3))-
     '        (DG(1,1)*DG(3,2)/DG(1,2)))
            X(2)=(b(1)-DG(1,1)*X(1))/DG(1,2)
            X(3)=(b(2)-DG(1,1)*X(1))/DG(2,3)
c            CALL DGETRF(3,3,DG,3,IPIV,INFO)
c            CALL DGETRS(TRANS,3,1,DG,3,IPIV,X,3,INFO)
C solves the Ax=b to determine the Newton Step to satisfy the
C equations of conservation
            PA=PA+X(1)
            PB=PB+X(2)
            PC=PC+X(3)
            SUM=DABS(x(1))+DABS(x(2))+DABS(x(3))
          ENDDO !convergence loop
C update the YQ array with the converged values from the Newton Rapsom
C iterations
          YQ(ny_a1p,1)=PA
          YQ(ny_b1p,1)=PB
          YQ(ny_c1p,1)=PC
          YQ(ny_a1r,1)=RA
          YQ(ny_b1r,1)=RB
          YQ(ny_c1r,1)=RC
          YQ(ny_a1v,1)= Fsign(1)*VA !is new arterial veloc at a1
          YQ(ny_b1v,1)=-Fsign(2)*VB !is new arterial veloc at b1
          YQ(ny_c1v,1)=-Fsign(3)*VC !is new arterial veloc at c1
          YQ(NYNQ(1,nq,0),6)=-OLDP_INITIAL+
     '      ((((YQ(ny_a1p,1)+YQ(ny_a1p,8)+G_TERMA*LA)/LA)
     '      +((YQ(ny_b1p,1)+YQ(ny_b1p,8)-G_TERMB*LB)/LB)
     '      +((YQ(ny_c1p,1)+YQ(ny_c1p,8)-G_TERMC*LC)/LC))
     '      /LOOP_DENOM)
C ------------------- end of arterial ----------------------

          IF (((VENOUS_NETWORK.EQ.'Y').OR.
     '      (VENOUS_NETWORK.EQ.'y')).AND.(N_VENOUS_GEOM.EQ.1)) THEN
            
C -------------- now do venous bifurcation -----------------
C The code below represents the process performed above for the arterial
C junction with the conservation equations reversed to account for
C the fact that velocities are now negative.
C Determine the grid points a1 and a2 and the ny values corresponding
C to velocity, radius and pressure
            nq_a1=BC_POINTS(1,1) !was CONECT(0,1,nq)
            ny_a1p=NYNQ(4,nq_a1,0)
            ny_a1r=NYNQ(5,nq_a1,0)
            ny_a1v=NYNQ(6,nq_a1,0)
            nq_a2=BC_POINTS(1,2) !was CONECT(-1,1,nq_a1)
            ny_a2p=NYNQ(4,nq_a2,0)
            ny_a2r=NYNQ(5,nq_a2,0)
            ny_a2v=NYNQ(6,nq_a2,0)
            nq_a3=BC_POINTS(1,3)
            ny_a3r=NYNQ(5,nq_a3,0)
            ny_a3v=NYNQ(6,nq_a3,0)
C determining the velocity at the grid a1 and a2 at the bifurcation
            Va1_new=Fsign(1)*YQ(ny_a1v,1) !is new veloc at a1
            Va2_new=Fsign(1)*YQ(ny_a2v,1) !is new veloc at a2
            Va3_mid=Fsign(1)*YQ(ny_a3v,5) !is mid veloc at a3
            Va1_old=Fsign(1)*YQ(ny_a1v,8) !is old veloc at a1
            Va2_old=Fsign(1)*YQ(ny_a2v,8) !is old veloc at a2
C determine the grid points b1 and b2 and the ny values corresponding
C to velocity, radius and pressure
            nq_b1=BC_POINTS(2,1) !was CONECT(1,1,nq)
            ny_b1p=NYNQ(4,nq_b1,0)
            ny_b1r=NYNQ(5,nq_b1,0)
            ny_b1v=NYNQ(6,nq_b1,0)
            nq_b2=BC_POINTS(2,2) !was CONECT(1,1,nq_b1)
            ny_b2p=NYNQ(4,nq_b2,0)
            ny_b2r=NYNQ(5,nq_b2,0)
            ny_b2v=NYNQ(6,nq_b2,0)
            nq_b3=BC_POINTS(2,3)
            ny_b3r=NYNQ(5,nq_b3,0)
            ny_b3v=NYNQ(6,nq_b3,0)
C determining the velocity at the grid b1 and b2 at the bifurcation
            Vb1_new=-Fsign(2)*YQ(ny_b1v,1) !is new veloc at b1
            Vb2_new=-Fsign(2)*YQ(ny_b2v,1) !is new veloc at b2
            Vb3_mid=-Fsign(2)*YQ(ny_b3v,5) !is mid veloc at b3
            Vb1_old=-Fsign(2)*YQ(ny_b1v,8) !is old veloc at b1
            Vb2_old=-Fsign(2)*YQ(ny_b2v,8) !is old veloc at b2
C determine the grid points c1 and c2 and the ny values corresponding
C to velocity, radius and pressure
            nq_c1=BC_POINTS(3,1) !was CONECT(1,2,nq)
            ny_c1p=NYNQ(4,nq_c1,0)
            ny_c1r=NYNQ(5,nq_c1,0)
            ny_c1v=NYNQ(6,nq_c1,0)
            nq_c2=BC_POINTS(3,2) !was CONECT(1,1,nq_c1)
            ny_c2p=NYNQ(4,nq_c2,0)
            ny_c2r=NYNQ(5,nq_c2,0)
            ny_c2v=NYNQ(6,nq_c2,0)
            nq_c3=BC_POINTS(3,3)
            ny_c3r=NYNQ(5,nq_c3,0)
            ny_c3v=NYNQ(6,nq_c3,0)
C determining the velocity at the grid c1 and c2 at the bifurcation
            Vc1_new=-Fsign(3)*YQ(ny_c1v,1) !is new veloc at c1
            Vc2_new=-Fsign(3)*YQ(ny_c2v,1) !is new veloc at c2
            Vc3_mid=-Fsign(3)*YQ(ny_c3v,5) !is mid veloc at c3
            Vc1_old=-Fsign(3)*YQ(ny_c1v,8) !is old veloc at c1
            Vc2_old=-Fsign(3)*YQ(ny_c2v,8) !is old veloc at c2
C calculating the unstressed radii at the points surrounding the
C bifurcation
            vien_ro_a1h=(CQ(4,nq_a1)+CQ(4,nq_a2))*VIEN_RATIO*0.5d0/
     '        ((LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)
     '        *0.25d0)**0.5d0
            vien_ro_b1h=(CQ(4,nq_b1)+CQ(4,nq_b2))*VIEN_RATIO*0.5d0/
     '        ((LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
     '        *0.25d0)**0.5d0

            vien_ro_c1h=(CQ(4,nq_c1)+CQ(4,nq_c2))*VIEN_RATIO*0.5d0/
     '        ((LAMBDA_C1+LAMBDA_C2+OLD_LAMBDA_C1+OLD_LAMBDA_C2)
     '        *0.25d0)**0.5d0
C unstressed radii at half time step i+1/2 and k+1/2
            vien_ro_a1=CQ(4,nq_a2)*VIEN_RATIO/(LAMBDA_A1**0.5d0)
            vien_ro_b1=CQ(4,nq_b2)*VIEN_RATIO/(LAMBDA_B1**0.5d0)
            vien_ro_c1=CQ(4,nq_c2)*VIEN_RATIO/(LAMBDA_C1**0.5d0)
C unstress radii at full time step i+1 and k+1
            BETAA1h=(CQ(9,nq_a1)+CQ(9,nq_a2))*0.5d0*
     '        ((LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)
     '        *0.25d0)+((CQ(10,nq_a1)+CQ(10,nq_a2))*0.5d0)
            BETAB1h=(CQ(9,nq_b1)+CQ(9,nq_b2))*0.5d0*
     '        ((LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
     '        *0.25d0)+((CQ(10,nq_b1)+CQ(10,nq_b2))*0.5d0)
            BETAC1h=(CQ(9,nq_c1)+CQ(9,nq_c2))*0.5d0*
     '        ((LAMBDA_C1+LAMBDA_C2+OLD_LAMBDA_C1+OLD_LAMBDA_C2)
     '        *0.25d0)+((CQ(10,nq_c1)+CQ(10,nq_c2))*0.5d0)
            BETAA1h=DMAX1(BETAA1h,1.0d0)
            BETAB1h=DMAX1(BETAB1h,1.0d0)
            BETAC1h=DMAX1(BETAC1h,1.0d0)
C calculate the pressure radius exponent  at the half time
C  step i+1/2 and k+1/2
            BETAA1=CQ(9,nq_a1)*LAMBDA_A1+CQ(10,nq_a1)
            BETAB1=CQ(9,nq_b1)*LAMBDA_B1+CQ(10,nq_b1)
            BETAC1=CQ(9,nq_c1)*LAMBDA_C1+CQ(10,nq_c1)
            BETAA1=DMAX1(BETAA1,1.0d0)
            BETAB1=DMAX1(BETAB1,1.0d0)
            BETAC1=DMAX1(BETAC1,1.0d0)
C calculate the pressure radius exponent at the half time
C step i+1 and k+1
            OLDP_INITIAL=YQ(NYNQ(4,nq,0),6)
C calculating Po
            LA=1.0d-9/((vien_ro_a1/CQ(4,1))**2.0d0)
            LB=1.0d-9/((vien_ro_b1/CQ(4,1))**2.0d0)
            LC=1.0d-9/((vien_ro_c1/CQ(4,1))**2.0d0)
C setting the realative inductance values k=1.0E-4
            FA2=PI*(YQ(ny_a2r,1)**2.0d0)*Va2_new
            OLDFA1=PI*(YQ(ny_a1r,8)**2.0d0)*Va1_old
            OLDFA2=PI*(YQ(ny_a2r,8)**2.0d0)*Va2_old
            FB2=PI*(YQ(ny_b2r,1)**2.0d0)*Vb2_new
            OLDFB1=PI*(YQ(ny_b1r,8)**2.0d0)*Vb1_old
            OLDFB2=PI*(YQ(ny_b2r,8)**2.0d0)*Vb2_old
            FC2=PI*(YQ(ny_c2r,1)**2.0d0)*Vc2_new
            OLDFC1=PI*(YQ(ny_c1r,8)**2.0d0)*Vc1_old
            OLDFC2=PI*(YQ(ny_c2r,8)**2.0d0)*Vc2_old
C calculating the flows at the bifurcation grid points
            Aa=((CQ(11,nq_a1)+CQ(11,nq_a2))*BETAA1h)
     '        *((YQ(ny_a3r,5)/vien_ro_a1h)
     '        **(BETAA1h-1.0d0))
     '        /(PI*YQ(ny_a3r,5)*vien_ro_a1h)
            Ab=((CQ(11,nq_b1)+CQ(11,nq_b2))*BETAB1h)
     '        *((YQ(ny_b3r,5)/vien_ro_b1h)
     '        **(BETAB1h-1.0d0))
     '        /(PI*YQ(ny_b3r,5)*vien_ro_b1h)
            Ac=((CQ(11,nq_c1)+CQ(11,nq_c2))*BETAC1h)
     '        *((YQ(ny_c3r,5)/vien_ro_c1h)
     '        **(BETAC1h-1.0d0))
     '        /(PI*YQ(ny_c3r,5)*vien_ro_c1h)
            PRESTERM=YQ(ny_a1p,8)+YQ(ny_a2p,8)-YQ(ny_a2p,1)
     '        +OLD_TRACE_A1+OLD_TRACE_A2-TRACE_A1-TRACE_A2
            DENOM=(1.0d0+(Aa*TINCR**2/(2.0d0*DELTA_XA*LA)))
            TERMREST=(Aa*(TINCR/DELTA_XA))*(FA2+OLDFA2-(2.0d0*OLDFA1)-
     '        ((TINCR/(2.0d0*LA))*(YQ(ny_a1p,8)-OLDP_INITIAL)))
            LAM_TERMA=(PI*Aa*(YQ(ny_a3r,5)**2.0d0)/LAMBDA_AH)*
     '        (LAMBDA_A1+LAMBDA_A2-OLD_LAMBDA_A1-OLD_LAMBDA_A2)
            H=(PRESTERM+TERMREST-LAM_TERMA)/DENOM
            I=(DENOM-1)/DENOM
            PRESTERM=YQ(ny_b1p,8)+YQ(ny_b2p,8)-YQ(ny_b2p,1)
     '        +OLD_TRACE_B1+OLD_TRACE_B2-TRACE_B1-TRACE_B2
            DENOM=(1.0d0+(Ab*(TINCR**2.0d0)/(2.0d0*DELTA_XB*LB)))
            TERMREST=(Ab*(TINCR/DELTA_XB))*(FB2+OLDFB2-(2.0d0*OLDFB1)-
     '        ((TINCR/(2.0d0*LB))*(OLDP_INITIAL-YQ(ny_b1p,8))))
            LAM_TERMB=(PI*Ab*(YQ(ny_b3r,5)**2.0d0)/LAMBDA_BH)*
     '        (LAMBDA_B1+LAMBDA_B2-OLD_LAMBDA_B1-OLD_LAMBDA_B2)
            J=(PRESTERM-TERMREST-LAM_TERMB)/DENOM
            K=(DENOM-1)/DENOM
            PRESTERM=YQ(ny_c1p,8)+YQ(ny_c2p,8)-YQ(ny_c2p,1)
     '        +OLD_TRACE_C1+OLD_TRACE_C2-TRACE_C1-TRACE_C2
            DENOM=(1.0d0+(Ac*(TINCR**2.0d0)/(2.0d0*DELTA_XC*LC)))
            TERMREST=(Ac*(TINCR/DELTA_XC))*(FC2+OLDFC2-(2.0d0*OLDFC1)-
     '        ((TINCR/(2.0d0*LC))*(OLDP_INITIAL-YQ(ny_c1p,8))))
            LAM_TERMC=(PI*Ac*(YQ(ny_c3r,5)**2.0d0)/LAMBDA_CH)*
     '        (LAMBDA_C1+LAMBDA_C2-OLD_LAMBDA_C1-OLD_LAMBDA_C2)
            L=(PRESTERM-TERMREST-LAM_TERMC)/DENOM
            M=(DENOM-1)/DENOM
            DENOM=((1.0d0/LA)+(1.0d0/LB)+(1.0d0/LC))
            DENOM1=1.0d0-(I/(LA*DENOM))-(K/(LB*DENOM))-(M/(LC*DENOM))
            P_INITIAL=(-OLDP_INITIAL+((((H+YQ(ny_a1p,8))/LA)+
     '        ((J+YQ(ny_b1p,8))/LB)+ ((L+YQ(ny_c1p,8))/LC))
     '        /DENOM))/DENOM1
            YQ(NYNQ(4,nq,0),6)=P_INITIAL
C calculated the initial guess for the pressure at the centre of the
C bifurcation
            PRESTERM=YQ(ny_a1p,8)+YQ(ny_a2p,8)-YQ(ny_a2p,1)
     '        +OLD_TRACE_A1+OLD_TRACE_A2-TRACE_A1-TRACE_A2
            DENOM=(1.0d0+(Aa*(TINCR**2.0d0)/(2.0d0*DELTA_XA*LA)))
            TERMREST=(AA*(TINCR/DELTA_XA))*(FA2+OLDFA2-(2.0d0*OLDFA1)-
     '        ((TINCR/(2.0d0*LA))*(YQ(ny_a1p,8)-
     '        P_INITIAL-OLDP_INITIAL)))
            YQ(ny_a1p,1)=(PRESTERM+TERMREST-LAM_TERMA)/DENOM
            PRESTERM=YQ(ny_b1p,8)+YQ(ny_b2p,8)-YQ(ny_b2p,1)
     '        +OLD_TRACE_B1+OLD_TRACE_B2-TRACE_B1-TRACE_B2
            DENOM=(1.0d0+(Ab*(TINCR**2.0d0)/(2.0d0*DELTA_XB*LB)))
            TERMREST=(Ab*(TINCR/DELTA_XB))*(FB2+OLDFB2-(2.0d0*OLDFB1)-
     '        ((TINCR/(2.0d0*LB))*
     '        (P_INITIAL+OLDP_INITIAL-YQ(ny_b1p,8))))
            YQ(ny_b1p,1)=(PRESTERM-TERMREST-LAM_TERMB)/DENOM
            PRESTERM=YQ(ny_c1p,8)+YQ(ny_c2p,8)-YQ(ny_c2p,1)
     '        +OLD_TRACE_C1+OLD_TRACE_C2-TRACE_C1-TRACE_C2
            DENOM=(1.0d0+(Ac*(TINCR**2.0d0)/(2.0d0*DELTA_XC*LC)))
            TERMREST=(Ac*(TINCR/DELTA_XC))*
     '        (FC2+OLDFC2-(2.0d0*OLDFC1)-((TINCR/(2.0d0*LC))*
     '        (P_INITIAL+OLDP_INITIAL-YQ(ny_c1p,8))))
            YQ(ny_c1p,1)=(PRESTERM-TERMREST-LAM_TERMC)/DENOM
C calculates the initial guess of the bifurcation pressures
C at grid points a1, b1, and c1
            IF((YQ(ny_a1p,1)+TRACE_A1).GT.0.0d0) THEN
              YQ(ny_a1r,1)=vien_ro_a1*((((YQ(ny_a1p,1)+TRACE_A1)
     '          /CQ(11,nq_a1))+1.0d0)**(1.0d0/BETAA1))
            ELSE
              fo=BETAA1*CQ(11,nq_a1)/CQ(7,nq_a1)
              YQ(ny_a1r,1)=(vien_ro_a1/((1-((YQ(ny_a1p,1)+TRACE_A1)/
     '          fo))**(1.0d0/CQ(7,nq_a1))))
            ENDIF
            IF((YQ(ny_b1p,1)+TRACE_B1).GT.0.0d0) THEN
              YQ(ny_b1r,1)=vien_ro_b1*((((YQ(ny_b1p,1)+TRACE_B1)
     '          /CQ(11,nq_b1))+1.0d0)**(1.0d0/BETAB1))
            ELSE
              fo=BETAB1*CQ(11,nq_b1)/CQ(7,nq_b1)
              YQ(ny_b1r,1)=(vien_ro_b1/((1-((YQ(ny_b1p,1)+TRACE_B1)/
     '        fo))**(1.0d0/CQ(7,nq_b1))))
            ENDIF
            IF((YQ(ny_c1p,1)+TRACE_C1).GT.0.0d0) THEN
              YQ(ny_c1r,1)=vien_ro_c1*((((YQ(ny_c1p,1)+TRACE_C1)
     '          /CQ(11,nq_c1))+1.0d0)**(1.0d0/BETAC1))
            ELSE
              fo=BETAC1*CQ(11,nq_c1)/CQ(7,nq_c1)
              YQ(ny_c1r,1)=(vien_ro_c1/((1-((YQ(ny_c1p,1)+TRACE_C1)/
     '          fo))**(1.0d0/CQ(7,nq_c1))))
            ENDIF
C calculates the inital radius values at the bifurcation grid points
            TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '        ((Va3_mid**2.0d0)/YQ(ny_a3r,5))*
     '        (YQ(ny_a1r,1)+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
            P_DIFF=YQ(ny_a1p,1)+YQ(ny_a1p,8)-YQ(ny_a2p,1)-YQ(ny_a2p,8)
     '        -2.d0*G_TERMA*DELTA_XA
            TERM2=(1.0D0/CQ(1,nq_a1))*P_DIFF
            IF(YQ(ny_a3r,5).GT.vien_ro_a1h) THEN
              TERM3=2.0d0*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)
     '          /(CQ(3,nq_a1)-1.0d0))*(Va3_mid/(YQ(ny_a3r,5)**2.0d0))
            ELSE
              RO_VIEN=vien_ro_a1h
              ELIP_A=((RO_VIEN**2.0d0)+
     '          (RO_VIEN**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_B=((RO_VIEN**2.0d0)-
     '          (RO_VIEN**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '          (YQ(ny_a3r,5)**2.0d0*(ELIP_A**2.0d0+
     '          ELIP_B**2.0d0)))**0.5D0
              TERM3=2.0d0*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)
     '          /(CQ(3,nq_a1)-1.0d0))*
     '          (Va3_mid/(ELIP_TERM**2.0d0))
            ENDIF
            TERM4=(1.0d0/YQ(ny_a3r,5))*((-DELTA_XA/TINCR)+
     '        ((2.0d0*CQ(3,nq_a1)-1.0d0)*Va3_mid))*
     '        ((YQ(ny_a2r,1)*(LAMBDA_A2**0.5D0))+
     '        (YQ(ny_a1r,1)*(LAMBDA_A1**0.5D0))-
     '        (YQ(ny_a2r,8)*(OLD_LAMBDA_A2**0.5D0))-
     '        (YQ(ny_a1r,8)*(OLD_LAMBDA_A1**0.5D0)))
            TERM5=Va3_mid/YQ(ny_a3r,5)
     '        *(YQ(ny_a1r,1)+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
            Va1_new=Va2_old+((TINCR/(2.0d0*DELTA_XA))*
     '        (TERM1-TERM2))-TERM3+TERM4-TERM5
            TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '        ((Vb3_mid**2.0d0)/YQ(ny_b3r,5))*
     '        (YQ(ny_b2r,1)+YQ(ny_b2r,8)-YQ(ny_b1r,1)-YQ(ny_b1r,8))
            P_DIFF=YQ(ny_b2p,1)+YQ(ny_b2p,8)-YQ(ny_b1p,1)-YQ(ny_b1p,8)
     '        +2.d0*G_TERMB*DELTA_XB
            TERM2=(1.0D0/CQ(1,nq_b1))*P_DIFF
            IF(YQ(ny_b3r,5).GT.vien_ro_b1h) THEN
              TERM3=2.0d0*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)/
     '          (CQ(3,nq_b1)-1.0d0))*(Vb3_mid/(YQ(ny_b3r,5)**2.0d0))
            ELSE
              RO_VIEN=vien_ro_b1h
              ELIP_A=((RO_VIEN**2.0d0)+
     '          (RO_VIEN**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_B=((RO_VIEN**2.0d0)-
     '          (RO_VIEN**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '          (YQ(ny_b3r,5)**2.0d0*(ELIP_A**2.0d0+
     '          ELIP_B**2.0d0)))**0.5D0
              TERM3=2.0d0*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)/
     '          (CQ(3,nq_b1)-1.0d0))*(Vb3_mid/(ELIP_TERM**2.0d0))
            ENDIF
            TERM4=(1.0d0/YQ(ny_b3r,5))*((DELTA_XB/TINCR)+
     '        ((2.0d0*CQ(3,nq_b1)-1.0d0)*Vb3_mid))*
     '        (((LAMBDA_B1**0.5d0)*YQ(ny_b1r,1))+
     '        ((LAMBDA_B2**0.5d0)*YQ(ny_b2r,1))-
     '        ((OLD_LAMBDA_B2**0.5d0)*YQ(ny_b2r,8))-
     '        ((OLD_LAMBDA_B1**0.5d0)*YQ(ny_b1r,8)))
            TERM5=Vb3_mid/YQ(ny_b3r,5)
     '        *(YQ(ny_b2r,1)+YQ(ny_b2r,8)-YQ(ny_b1r,1)-YQ(ny_b1r,8))
            Vb1_new=Vb2_old+((TINCR/(2.0d0*DELTA_XB))*
     '        (TERM1-TERM2))-TERM3+TERM4+TERM5
            TERM1=((2.0d0*CQ(3,nq_c1))-0.0d0)*
     '        ((Vc3_mid**2.0d0)/YQ(ny_c3r,5))*
     '        (YQ(ny_c2r,1)+YQ(ny_c2r,8)-YQ(ny_c1r,1)-YQ(ny_c1r,8))
            P_DIFF=YQ(ny_c2p,1)+YQ(ny_c2p,8)-YQ(ny_c1p,1)-YQ(ny_c1p,8)
     '        +2.d0*G_TERMC*DELTA_XC
            TERM2=(1.0D0/CQ(1,nq_c1))*P_DIFF
            IF(YQ(ny_c3r,5).GT.vien_ro_c1h) THEN
              TERM3=2.0d0*TINCR*CQ(2,nq_c1)*(CQ(3,nq_c1)
     '          /(CQ(3,nq_c1)-1.0d0))*
     '          (Vc3_mid/(YQ(ny_c3r,5)**2.0d0))
            ELSE
              RO_VIEN=vien_ro_c1h
              ELIP_A=((RO_VIEN**2.0d0)+
     '          (RO_VIEN**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_B=((RO_VIEN**2.0d0)-
     '          (RO_VIEN**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '          (YQ(ny_c3r,5)**2.0d0*(ELIP_A**2.0d0+
     '          ELIP_B**2.0d0)))**0.5D0
              TERM3=2.0d0*TINCR*CQ(2,nq_c1)*
     '          (CQ(3,nq_c1)/(CQ(3,nq_c1)-1.0d0))*
     '          (Vc3_mid/(ELIP_TERM**2.0d0))
            ENDIF
            TERM4=(1.0d0/YQ(ny_c3r,5))*((DELTA_XC/TINCR)+
     '        ((2.0d0*CQ(3,nq_c1)-1.0d0)*Vc3_mid))*
     '        (((LAMBDA_C1**0.5d0)*YQ(ny_c1r,1))+
     '        ((LAMBDA_C2**0.5d0)*YQ(ny_c2r,1))-
     '        ((OLD_LAMBDA_C2**0.5d0)*YQ(ny_c2r,8))-
     '        ((OLD_LAMBDA_C1**0.5d0)*YQ(ny_c1r,8)))
            TERM5=Vc3_mid/YQ(ny_c3r,5)
     '        *(YQ(ny_c2r,1)+YQ(ny_c2r,8)-YQ(ny_c1r,1)-YQ(ny_c1r,8))
            Vc1_new=Vc2_old+((TINCR/(2.0d0*DELTA_XC))*
     '        (TERM1-TERM2))-TERM3+TERM4+TERM5
C calculates the initial guess for the velocity values at grid points
C surrouding a bifurcation. To account for the increased resistance with
C collapse these terms for each point are calculated using a modifed
C term.
C Calculate flows in 3 segments
            FA1=PI*(YQ(ny_a1r,1)**2.0d0)*Va1_new
            FB1=PI*(YQ(ny_b1r,1)**2.0d0)*Vb1_new
            FC1=PI*(YQ(ny_c1r,1)**2.0d0)*Vc1_new
            LOOP_DENOM=((1.0d0/LA)+(1.0d0/LB)+(1.0d0/LC))
C new for veins NPS 12/11/98
C Newton-Raphson loop for venous bifurcation conservation eqns
            PA=YQ(ny_a1p,1)
            PB=YQ(ny_B1p,1)
            PC=YQ(ny_C1p,1)
            SUM=1.0D0
            NEW_COUNT=0
            DO WHILE(((SUM.GT.(CONVERG_TOL*1000.0d0)).AND.
     '        ((dabs(FA1-FB1-FC1)).GT.(CONVERG_TOL*1000.0d0))
     '        .OR.(NEW_COUNT.LT.5)).AND.(NEW_COUNT.LT.50))
              NEW_COUNT=NEW_COUNT+1
              IF((PA+TRACE_A1).GT.0.0d0) THEN
                Sa=PI*(vien_ro_a1**2.0d0)*
     '            ((((PA+trace_a1)/CQ(11,nq_a1))+1)
     '            **(2.0D0/BETAA1))
                dSa=(PI*(vien_ro_a1**2.0d0)*2.0d0/
     '            (BETAA1*CQ(11,nq_a1))
     '            )*((((PA+trace_a1)/CQ(11,nq_a1))+1)
     '            **((2.0D0/BETAA1)-1.0d0))
                Ra=vien_ro_a1*((((PA+trace_a1)/CQ(11,nq_a1))+1)
     '            **(1.0D0/BETAA1))
                dRa=(vien_ro_a1/
     '            (BETAA1*CQ(11,nq_a1)))
     '            *((((PA+trace_a1)/CQ(11,nq_a1))+1)
     '            **((1.0D0/BETAA1)-1.0d0))
              ELSE
                fo=BETAA1*CQ(11,nq_a1)/CQ(7,nq_a1)
                Sa=PI*((vien_ro_a1/((1-((PA+TRACE_A1)/
     '            fo))**(1.0d0/CQ(7,nq_a1))))**2.0d0)
                dSa=(2.0D0*PI*(vien_ro_a1**2.0d0)/(CQ(7,nq_a1)*fo))
     '            *(1-(PA+TRACE_A1)/Fo)**(-2.0d0/CQ(7,nq_a1)-1.0d0)
                Ra=(vien_ro_a1/((1-((PA+TRACE_A1)/
     '            fo))**(1.0d0/CQ(7,nq_a1))))
                dRa=(vien_ro_a1/(CQ(7,nq_a1)*fo))
     '            *(1-(PA+TRACE_A1)/Fo)**(-1.0d0/CQ(7,nq_a1)-1.0d0)
              ENDIF !((PA+TRACE_A1).GT.0.0d0)
C calculates the radius, area, and derivitives with respect to pressure
C at the  "a" grid point adjacent to a bifurcation
              IF((PB+TRACE_B1).GT.0.0d0) THEN
                Sb=PI*(vien_ro_b1**2.0d0)*
     '            ((((PB+trace_b1)/CQ(11,nq_b1))+1)
     '            **(2.0D0/BETAB1))
                dSb=(PI*(vien_ro_b1**2.0d0)*2.0d0/
     '            (BETAB1*CQ(11,nq_b1))
     '            )*((((Pb+trace_b1)/CQ(11,nq_b1))+1)
     '            **((2.0D0/BETAB1)-1.0d0))
                Rb=vien_ro_b1*((((PB+trace_b1)/CQ(11,nq_b1))+1)
     '            **(1.0D0/BETAB1))
                dRb=(vien_ro_b1/
     '            (BETAB1*CQ(11,nq_b1)))
     '            *((((Pb+trace_b1)/CQ(11,nq_b1))+1)
     '            **((1.0D0/BETAB1)-1.0d0))
              ELSE
                fo=BETAB1*CQ(11,nq_b1)/CQ(7,nq_b1)
                Sb=PI*((vien_ro_b1/((1-((PB+TRACE_B1)/
     '            fo))**(1.0d0/CQ(7,nq_b1))))**2.0d0)
                dSb=(2.0D0*PI*(vien_ro_b1**2.0d0)/(CQ(7,nq_b1)*fo))
     '            *(1-(PB+TRACE_B1)/Fo)**(-2.0d0/CQ(7,nq_b1)-1.0d0)
                Rb=(vien_ro_b1/((1-((PB+TRACE_B1)/
     '            fo))**(1.0d0/CQ(7,nq_b1))))
                dRb=(vien_ro_b1/(CQ(7,nq_b1)*fo))
     '            *(1-(PB+TRACE_B1)/Fo)**(-1.0d0/CQ(7,nq_b1)-1.0d0)
              ENDIF !((PB+TRACE_B1).GT.0.0d0)
C calculates the radius, area, and derivitives with respect to pressue
C at the  "b" grid point adjacent to a bifurcation
              IF((PC+TRACE_C1).GT.0.0d0) THEN
                Sc=PI*(vien_ro_c1**2.0d0)*
     '            ((((PC+trace_c1)/CQ(11,nq_c1))+1)
     '            **(2.0D0/BETAC1))
                dSc=(PI*(vien_ro_c1**2.0d0)*2.0d0/
     '            (BETAC1*CQ(11,nq_c1))
     '            )*((((PC+trace_c1)/CQ(11,nq_c1))+1)
     '            **((2.0D0/BETAC1)-1.0d0))
                Rc=vien_ro_c1*((((PC+trace_c1)/CQ(11,nq_c1))+1)
     '            **(1.0D0/BETAC1))
                dRc=(vien_ro_c1/
     '            (BETAC1*CQ(11,nq_c1)))
     '            *((((PC+trace_c1)/CQ(11,nq_c1))+1)
     '            **((1.0D0/BETAC1)-1.0d0))
              ELSE
                fo=BETAC1*CQ(11,nq_c1)/CQ(7,nq_c1)
                Sc=PI*((vien_ro_c1/((1-((PC+TRACE_C1)/
     '            fo))**(1.0d0/CQ(7,nq_c1))))**2.0d0)
                dSc=(2.0D0*PI*(vien_ro_c1**2.0d0)/(CQ(7,nq_c1)*fo))
     '            *(1-(PC+TRACE_C1)/Fo)**(-2.0d0/CQ(7,nq_c1)-1.0d0)
                Rc=(vien_ro_c1/((1-((PC+TRACE_C1)/
     '            fo))**(1.0d0/CQ(7,nq_c1))))
                dRc=(vien_ro_c1/(CQ(7,nq_c1)*fo))
     '            *(1-(PC+TRACE_C1)/Fo)**(-1.0d0/CQ(7,nq_c1)-1.0d0)
              ENDIF !((PC+TRACE_C1).GT.0.0d0)
C calculates the radius, area, and derivitives with respect to pressue
C at the  "c" grid point adjacent to a bifurcation
              TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '          ((Va3_mid**2.0d0)/YQ(ny_a3r,5))*
     '          (RA+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
              TERM2=(1.0D0/CQ(1,nq_a1))*(PA+YQ(ny_a1p,8)-YQ(ny_a2p,1)
     &          -YQ(ny_a2p,8)-2.d0*G_TERMA*DELTA_XA)
              IF(YQ(ny_a3r,5).GT.
     '          vien_ro_a1h) THEN
                TERM3=2.0d0*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)
     '            /(CQ(3,nq_a1)-1.0d0))*
     '            (Va3_mid/(YQ(ny_a3r,5)**2.0d0))
              ELSE
                RO_VIEN=vien_ro_a1h
                ELIP_A=((RO_VIEN**2.0d0)+
     '            (RO_VIEN**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_B=((RO_VIEN**2.0d0)-
     '            (RO_VIEN**4.0d0-YQ(ny_a3r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '            (YQ(ny_a3r,5)**2.0d0*(ELIP_A**2.0d0+
     '            ELIP_B**2.0d0)))**0.5D0
                TERM3=2.0d0*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)
     '            /(CQ(3,nq_a1)-1.0d0))*
     '            (Va3_mid/(ELIP_TERM**2.0d0))
              ENDIF
              TERM4=(1.0d0/YQ(ny_a3r,5))*((-DELTA_XA/TINCR)+
     '          ((2.0d0*CQ(3,nq_a1)-1.0d0)*Va3_mid))*
     '          ((YQ(ny_a2r,1)*(LAMBDA_A2**0.5D0))+
     '          (RA*(LAMBDA_A1**0.5D0))-
     '          (YQ(ny_a2r,8)*(OLD_LAMBDA_A2**0.5D0))-
     '          (YQ(ny_a1r,8)*(OLD_LAMBDA_A1**0.5D0)))
              TERM5=Va3_mid/YQ(ny_a3r,5)
     '          *(RA+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
              VA=Va2_old+((TINCR/(2.0d0*DELTA_XA))*
     '          (TERM1-TERM2))-TERM3+TERM4-TERM5
C uses the end point boundary condtions to calculate velocity at grid
C point "a" and the derivative with respect to pressure
              TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '          ((Va3_mid**2.0d0)/YQ(ny_a3r,5))*dRa
              TERM2=1.0d0/CQ(1,nq_a1)
              TERM3=((LAMBDA_A1**0.5d0)/YQ(ny_a3r,5))*
     '          ((-DELTA_XA/TINCR)+
     '          ((2.0d0*CQ(3,nq_a1)-1.0d0)*Va3_mid))*dRa
              TERM4=Va3_mid/YQ(ny_a3r,5)*dRa
              DVA=((TINCR/(2.0d0*DELTA_XA))*
     '          (TERM1-TERM2))+TERM3-TERM4
              TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '          ((Vb3_mid**2.0d0)/YQ(ny_b3r,5))*
     '          (YQ(ny_b2r,1)+YQ(ny_b2r,8)-Rb-YQ(ny_b1r,8))
              TERM2=(1.0D0/CQ(1,nq_b1))*(YQ(ny_b2p,1)
     '          +YQ(ny_b2p,8)-PB-YQ(ny_b1p,8)+2.d0*G_TERMB*DELTA_XB)
              IF(YQ(ny_b3r,5).GT.
     '          vien_ro_b1h) THEN
                TERM3=2.0d0*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)
     '            /(CQ(3,nq_b1)-1.0d0))*
     '            (Vb3_mid/(YQ(ny_b3r,5)**2.0d0))
              ELSE
                RO_VIEN=vien_ro_b1h
                ELIP_A=((RO_VIEN**2.0d0)+
     '            (RO_VIEN**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_B=((RO_VIEN**2.0d0)-
     '            (RO_VIEN**4.0d0-YQ(ny_b3r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '            (YQ(ny_b3r,5)**2.0d0*(ELIP_A**2.0d0+
     '            ELIP_B**2.0d0)))**0.5D0
                TERM3=2.0d0*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)
     '            /(CQ(3,nq_b1)-1.0d0))*
     '            (Vb3_mid/(ELIP_TERM**2.0d0))
              ENDIF
              TERM4=(1.0d0/YQ(ny_b3r,5))*((DELTA_XB/TINCR)+
     '          ((2.0d0*CQ(3,nq_b1)-1.0d0)*Vb3_mid))*
     '          (((LAMBDA_B1**0.5d0)*Rb)+
     '          ((LAMBDA_B2**0.5d0)*YQ(ny_b2r,1))-
     '          ((OLD_LAMBDA_B2**0.5d0)*YQ(ny_b2r,8))-
     '          ((OLD_LAMBDA_B1**0.5d0)*YQ(ny_b1r,8)))
              TERM5=Vb3_mid/YQ(ny_b3r,5)
     '          *(YQ(ny_b2r,1)+YQ(ny_b2r,8)-Rb-YQ(ny_b1r,8))
              VB=Vb2_old+((TINCR/(2.0d0*DELTA_XB))*
     '          (TERM1-TERM2))-TERM3+TERM4+TERM5
C uses the end point bounary condtions to calculate velocity at grid
C point "b" and the derivative with respect to pressure
              TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '          ((Vb3_mid**2.0d0)/YQ(ny_b3r,5))*(-dRb)
              TERM2=-1.0d0/CQ(1,nq_b1)
              TERM3=(1.0d0/YQ(ny_b3r,5))*((DELTA_XB/TINCR)+
     '          ((2.0d0*CQ(3,nq_b1)-1.0d0)*Vb3_mid))*
     '          (LAMBDA_B1**0.5d0)*dRb
              TERM4=Vb3_mid/YQ(ny_b3r,5)*dRb
              DVB=((TINCR/(2.0d0*DELTA_XB))*
     '          (TERM1-TERM2))+TERM3-TERM4
              TERM1=((2.0d0*CQ(3,nq_c1))-0.0d0)*
     '          ((Vc3_mid**2.0d0)/YQ(ny_c3r,5))*
     '          (YQ(ny_c2r,1)+YQ(ny_c2r,8)-Rc-YQ(ny_c1r,8))
              TERM2=(1.0D0/CQ(1,nq_c1))*(YQ(ny_c2p,1)
     '          +YQ(ny_c2p,8)-PC-YQ(ny_c1p,8)+2.d0*G_TERMB*DELTA_XB)
              IF(YQ(ny_c3r,5).GT.vien_ro_c1h) THEN
                TERM3=2.0d0*TINCR*CQ(2,nq_c1)*(CQ(3,nq_c1)
     '            /(CQ(3,nq_c1)-1.0d0))*
     '            (Vc3_mid/(YQ(ny_c3r,5)**2.0d0))
              ELSE
                RO_VIEN=vien_ro_c1h
                ELIP_A=((RO_VIEN**2.0d0)+
     '            (RO_VIEN**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_B=((RO_VIEN**2.0d0)-
     '            (RO_VIEN**4.0d0-YQ(ny_c3r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '            (YQ(ny_c3r,5)**2.0d0*(ELIP_A**2.0d0+
     '            ELIP_B**2.0d0)))**0.5D0
                TERM3=2.0d0*TINCR*CQ(2,nq_c1)*
     '            (CQ(3,nq_c1)/(CQ(3,nq_c1)-1.0d0))*
     '            (Vc3_mid/(ELIP_TERM**2.0d0))
              ENDIF
              TERM4=(1.0d0/YQ(ny_c3r,5))*((DELTA_XC/TINCR)+
     '          ((2.0d0*CQ(3,nq_c1)-1.0d0)*Vc3_mid))*
     '          (((LAMBDA_C1**0.5d0)*Rc)+
     '          ((LAMBDA_C2**0.5d0)*YQ(ny_c2r,1))-
     '          ((OLD_LAMBDA_C2**0.5d0)*YQ(ny_c2r,8))-
     '          ((OLD_LAMBDA_C1**0.5d0)*YQ(ny_c1r,8)))
              TERM5=Vc3_mid/YQ(ny_c3r,5)
     '          *(YQ(ny_c2r,1)+YQ(ny_c2r,8)-Rc-YQ(ny_c1r,8))
              VC=Vc2_old+((TINCR/(2.0d0*DELTA_XC))*
     '          (TERM1-TERM2))-TERM3+TERM4+TERM5
C uses the end point bounary condtions to calculate velocity at grid
C point "c" and the derivative with respect to pressure
              TERM1=((2.0d0*CQ(3,nq_c1))-0.0d0)*
     '          ((Vc3_mid**2.0d0)/YQ(ny_c3r,5))*(-dRc)
              TERM2=-1.0d0/CQ(1,nq_c1)
              TERM3=(1.0d0/YQ(ny_c3r,5))*((DELTA_XC/TINCR)+
     '          ((2.0d0*CQ(3,nq_c1)-1.0d0)*Vc3_mid))*
     '          (LAMBDA_C1**0.5d0)*dRc
              TERM4=Vc3_mid/YQ(ny_c3r,5)*dRc
              DVC=((TINCR/(2.0d0*DELTA_XC))*
     '          (TERM1-TERM2))+TERM3-TERM4
              dFA=(Sa*DVA)+(dSa*VA)
              dFb=(Sb*DVb)+(dSb*Vb)
              dFc=(Sc*DVc)+(dSc*Vc)
              FA1=PI*(Ra**2.0d0)*va
              FB1=PI*(Rb**2.0d0)*vb
              FC1=PI*(Rc**2.0d0)*vc
C calculates flow and their derivative with respect to pressure in
C in each of the bifurcation segments
              FA1_OLD=PI*(YQ(ny_a1r,8)**2.0d0)*Va1_old
              FB1_OLD=PI*(YQ(ny_b1r,8)**2.0d0)*Vb1_old
              FC1_OLD=PI*(YQ(ny_c1r,8)**2.0d0)*Vc1_old
              DG(1,1)=1-((2.0d0*LA/TINCR)*dFA)
              DG(1,2)=-1-((2.0d0*LB/TINCR)*dFB)
              DG(1,3)=0.0D0
              DG(2,1)=DG(1,1)
              DG(2,2)=0.0D0
              DG(2,3)=-1-((2.0d0*LC/TINCR)*dFC)
              DG(3,1)=dFA
              DG(3,2)=-dFB
              DG(3,3)=-dFC
              X(1)=((2.0d0*LA/TINCR)*(FA1-FA1_OLD))+
     '          ((2.0d0*Lb/TINCR)*(FB1-FB1_OLD))+
     '          PB+YQ(ny_b1p,8)-PA-YQ(ny_a1p,8)
              X(2)=((2.0d0*LA/TINCR)*(FA1-FA1_OLD))+
     '          ((2.0d0*Lc/TINCR)*(FC1-FC1_OLD))+
     '          PC+YQ(ny_C1p,8)-PA-YQ(ny_a1p,8)
              X(3)=FC1+FB1-FA1
              B(1)=X(1)
              B(2)=X(2)
              B(3)=X(3)
              X(1)=(B(3)-(B(1)*DG(3,2)/DG(1,2))-(B(2)*DG(3,3)/DG(2,3)))/
     '          (DG(3,1)-(DG(1,1)*DG(3,3)/DG(2,3))-
     '          (DG(1,1)*DG(3,2)/DG(1,2)))
              X(2)=(B(1)-DG(1,1)*X(1))/DG(1,2)
              X(3)=(B(2)-DG(1,1)*X(1))/DG(2,3)
c              CALL DGETRF(3,3,DG,3,IPIV,INFO)
c              CALL DGETRS(TRANS,3,1,DG,3,IPIV,X,3,INFO)
C solves the Ax=b to determine the Newton Step to satisfy the
C equations of conservation
              PA=PA+X(1)
              PB=PB+X(2)
              PC=PC+X(3)
              SUM=DABS(x(1))+DABS(x(2))+DABS(x(3))
            ENDDO !convergence loop
C update the YQ array with the converged values from the Newton Rapsom
C iterations
            YQ(ny_a1p,1)=PA
            YQ(ny_b1p,1)=PB
            YQ(ny_c1p,1)=PC
            YQ(ny_a1r,1)=RA
            YQ(ny_b1r,1)=RB
            YQ(ny_c1r,1)=RC
            YQ(ny_a1v,1)= Fsign(1)*VA !is new venous veloc at a1
            YQ(ny_b1v,1)=-Fsign(2)*VB !is new venous veloc at b1
            YQ(ny_c1v,1)=-Fsign(3)*VC !is new venous veloc at c1
            YQ(NYNQ(4,nq,0),6)=-OLDP_INITIAL+
     '        ((((YQ(ny_a1p,1)+YQ(ny_a1p,8))/LA)
     '        +((YQ(ny_b1p,1)+YQ(ny_b1p,8))/LB)
     '        +((YQ(ny_c1p,1)+YQ(ny_c1p,8))/LC))
     '        /LOOP_DENOM)

C ----------------- end of venous section  ------------------

          ENDIF  ! venous_network & identical venous net_work

        ENDIF !nonterminal pt/bifurcation


      CALL EXITS('BRANCH1')
      RETURN
 9999 CALL ERRORS('BRANCH1',ERROR)
      CALL EXITS('BRANCH1')
      RETURN 1
      END


