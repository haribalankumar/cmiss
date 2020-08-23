      SUBROUTINE BRANCH2(BC_POINTS,BRANCH,CONECT,ITYP12,CQ,no_nq,nr,
     '  NTIME_POINTS,NTIME_NR,NYNQ,TIME,TIME_VALUES,XQ,YQ,ERROR,*)


C#### Subroutine: BRANCH2
C###  Description:
C###    BRANCH2 calculates the flow through the inflow and outflow grid
C###    points of a mesh which are connected by the lumped parameter
C###    microcirculation model. Each of the lumped parameter model
C###    components is represented using a rational polynomial where the
C###    individual components are entered in the .ipmate file. The
C###    mathematical technique is similar to that used in BRANCH1 where
C###    the Newton Raphson method is used to find a simultaneous
C###    solution to the boundary conditions of the discrete model and
C###    the lumped parameter model equations(see section 4.2 in the PhD
C###    thesis of Nic Smith). This routine is parallelised using
C###    dynamic scheduling over the number of vessel endings. Due to
C###    the even work associated with the iterative solution procedure
C###    scalability is reduced for a large number of processors

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
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
      INTEGER  BC_POINTS,CONECT(-1:1,0:2,NQM),ITYP12,nr,
     &  NTIME_POINTS(NTIMEVARSM),NTIME_NR(0:NTIMEVARSM,NRM),
     &  NYNQ(NHM,NQM,0:NRCM)
      REAL*8  CQ(NMM,NQM),TIME,
     '  TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM),
     '  XQ(NJM,NQM),YQ(NYQM,NIQM)
      CHARACTER ERROR*(*)
      LOGICAL BRANCH(NQM)
!     Local Variables
      INTEGER iterate_count,nj,
     '  NJTOT,no_nq,nq,nq_a1,nq_a2,nq_b1,nq_b2,
     '  nq_previous,NT_CYCLE(NTIMEVARSM),ntp,no_ntv,ntv,ntv1,ny_a1p,
     &  ny_a1r,ny_a1v,ny_a2p,ny_a2r,ny_a2v,ny_b1p,ny_b1r,ny_b1v,
     '  ny_b2p,ny_b2r,ny_b2v,ny_p,ny_r,ny_v,ny_p2,ny_r2,
     '  ny_v2
      REAL*8 A(2,2),Ain(NTIMEPOINTSM,NTIMEVARSM),B(2),
     '  Bin(NTIMEPOINTSM,NTIMEVARSM),BSUM,BSUM_OLD,BCONVERGE,
     '  BETAA1h,BETAB1h,BETAA1,BETAB1,BETA,BETAh,COM1,COM1_OLD,COM2,
     '  COM2_OLD,DCOM1,DCOM2,DELTA(NJT),DELTA_XA,DELTA_XB,DELTA_X,DFA,
     &  DFB,DFCPA,DFCPB,DP1DPA,DP2DPV,DRA,DRART,DRB,DRCAP_PA,DRCAP_PB,
     '  DRVIE,dSa,DSB,DVA,DVB,DX,ELIP_A,ELIP_B,ELIP_TERM,
     '  FA1,FA1_OLD,FB1,FB1_OLD,FC,FC_OLD,fo,G_BC,G_PLEURAL,G_TERM,
     &  HEIGHT(NJT),HEIGHT_A1(NJT),HEIGHT_B1(NJT),OLD_TRACE_A1,
     &  LAMBDA_1OLD,LAMBDA_1,LAMBDA_2OLD,LAMBDA_2,LAMBDA_A1,LAMBDA_A2,
     &  OLD_LAMBDA_A1,OLD_LAMBDA_A2,LAMBDA_B1,LAMBDA_B2,OLD_LAMBDA_B1,
     &  OLD_LAMBDA_B2,P1,P2,PA,PA_OLD,PB,PB_OLD,PLEURAL_DENSITY,Ra,RART,
     &  RART_OLD,RB,RCAP,RO_ART,RO_VIEN,RVIE,ROIA1h,ROIB1h,ROIA1,ROI,
     &  ROIh,RVIE_OLD,SA,SB,T_frac(NTIMEPOINTSM,NTIMEVARSM), !SMAR009 22/12/98 TERM6
     &  T_period(NTIMEVARSM),TRACE_A1,TERM1,TERM2,TERM3,TERM4,TERM5,
     &  TRACE,trace_b1,VA,VB,vien_ro_b1,X(2),XSUM,XSUM_OLD,XCONVERGE
C      CHARACTER TRANS
C      PARAMETER (TRANS='N')
      EXTERNAL DGETRS,DGETRF

      CALL ENTERS('BRANCH2',*9999)
      
      NJTOT=NJ_LOC(NJL_GEOM,0,nr)
      nq=BC_POINTS
      IF (.NOT.BRANCH(no_nq)) THEN ! not a bifurcation
C      IF(NXQ(1,0,nq,1).LE.1.AND.NXQ(-1,0,nq,1).LE.1) THEN
        IF (NQ_START(nr).ne.nq) THEN !arterial end point
C PM 26-JUL-01
C         IF(FLOW_CONST.GT.ZERO_TOL) THEN

          IF(((VENOUS_NETWORK.EQ.'Y').OR.(VENOUS_NETWORK.EQ.'y'))
     '     .AND.(N_VENOUS_GEOM.EQ.1)) THEN
 !venous network with identical geometry
            iterate_count=0
            nq_a1=CONECT(0,1,nq)
C            nq_a1=NXQ(0,1,nq,1)
            ny_a1p=NYNQ(1,nq_a1,0)
            ny_a1r=NYNQ(2,nq_a1,0)
            ny_a1v=NYNQ(3,nq_a1,0)
            nq_a2=CONECT(-1,1,nq_a1)
C            nq_a2=NXQ(-1,1,nq_a1,1)
            
            ny_a2p=NYNQ(1,nq_a2,0)
            ny_a2r=NYNQ(2,nq_a2,0)
            ny_a2v=NYNQ(3,nq_a2,0)
            nq_b1=CONECT(0,1,nq)
C            nq_b1=NXQ(0,1,nq,1)
            
            ny_b1p=NYNQ(4,nq_b1,0)
            ny_b1r=NYNQ(5,nq_b1,0)
            ny_b1v=NYNQ(6,nq_b1,0)
            nq_b2=CONECT(-1,1,nq_b1)
C            nq_b2=NXQ(-1,1,nq_b1,1)            
            
            ny_b2p=NYNQ(4,nq_b2,0)
            ny_b2r=NYNQ(5,nq_b2,0)
            ny_b2v=NYNQ(6,nq_b2,0)
C determing the YQ indeices for the flow quanties at the grid points
C adjacent to a bifurcation

            LAMBDA_A1=YQ(ny_a1v,2)*(TIME-T0)/(T1-T0)+
     '        YQ(ny_a1v,9)*(T1-TIME)/(T1-T0)
            LAMBDA_A2=YQ(ny_a2v,2)*(TIME-T0)/(T1-T0)+
     '        YQ(ny_a2v,9)*(T1-TIME)/(T1-T0)
            OLD_LAMBDA_A1=(YQ(ny_a1v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '        (YQ(ny_a1v,9)*(T1-TIME+TINCR)/(T1-T0))
            OLD_LAMBDA_A2=(YQ(ny_a2v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '        (YQ(ny_a2v,9)*(T1-TIME+TINCR)/(T1-T0))

            LAMBDA_B1=(YQ(ny_b1v,2)*(TIME-T0)/(T1-T0))+
     '        (YQ(ny_b1v,9)*(T1-TIME)/(T1-T0))
            LAMBDA_B2=(YQ(ny_b2v,2)*(TIME-T0)/(T1-T0))+
     '        (YQ(ny_b2v,9)*(T1-TIME)/(T1-T0))
            OLD_LAMBDA_B1=(YQ(ny_b1v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '        (YQ(ny_b1v,9)*(T1-TIME+TINCR)/(T1-T0))
            OLD_LAMBDA_B2=(YQ(ny_b2v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '        (YQ(ny_b2v,9)*(T1-TIME+TINCR)/(T1-T0))

C calculating the vessel stretch ratios at the grid points
C adjacent to a bifurcation by interpolating between the values
C at the begining and end of the mechanics step


            ROIA1h=(CQ(4,nq_a1)+CQ(4,nq_a2))*0.5d0/
     '        ((LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)
     '        *0.25d0)**0.5d0
            ROIB1h=(CQ(4,nq_b1)+CQ(4,nq_b2))*0.5d0*VIEN_RATIO/
     '        ((LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
     '        *0.25d0)**0.5d0

            ROIA1=CQ(4,nq_a1)/(LAMBDA_A1**0.5d0)
            vien_ro_b1=CQ(4,nq_b1)*VIEN_RATIO
     '        /(LAMBDA_B1**0.5d0)

            BETAA1h=(CQ(6,nq_a1)+CQ(6,nq_a2))*0.5d0*
     '        ((LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+OLD_LAMBDA_A2)
     '        *0.25d0)+((CQ(12,nq_a1)+CQ(12,nq_a2))*0.5d0)

            BETAB1h=(CQ(9,nq_b1)+CQ(9,nq_b2))*0.5d0*
     '        ((LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
     '        *0.25d0)+((CQ(10,nq_b1)+CQ(10,nq_b2))*0.5d0)

            BETAA1=CQ(6,nq_a1)*LAMBDA_A1+CQ(12,nq_a1)
            BETAB1=CQ(9,nq_b1)*LAMBDA_B1+CQ(10,nq_b1)
            BETAA1h=DMAX1(BETAA1h,1.0d0)
            BETAB1h=DMAX1(BETAB1h,1.0d0)
            BETAA1=DMAX1(BETAA1,1.0d0)
            BETAB1=DMAX1(BETAB1,1.0d0)

C determining the unstressed radius and vessel wall exponents

            FA1_OLD=PI*(YQ(ny_a1r,8)**2.0d0)*YQ(ny_a1v,8)
            FB1_OLD=-PI*(YQ(ny_b1r,8)**2.0d0)*YQ(ny_b1v,8)

C.. TRACE = external pressure force
            IF(ITYP12.EQ.2) THEN !pulmonary flow - calculate pleural pressure
              PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
              DO nj=1,NJT
                HEIGHT_A1(nj)=XQ_IN(nj)-XQ(nj,nq_a1)
                HEIGHT_B1(nj)=XQ_IN(nj)-XQ(nj,nq_b1)
              ENDDO
              G_PLEURAL=0.d0
              DO nj=1,NJT
                G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &            *GRAVITY*1000.d0*HEIGHT_A1(nj)
              ENDDO
              TRACE_A1=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
              OLD_TRACE_A1=TRACE_A1 !currently not time-dependent
            ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
              nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
              TRACE_A1=CQ(nj,nq)/1000.d0 !Pa->kPa
              OLD_TRACE_A1=TRACE_A1 !currently not time-dependent
            ELSE
              TRACE_A1=(YQ(ny_a1p,2)*(TIME-T0)/(T1-T0))+
     '          (YQ(ny_a1p,9)*(T1-TIME)/(T1-T0))
              OLD_TRACE_A1=(YQ(ny_a1p,2)*(TIME-T0-TINCR)/(T1-T0))+
     '          (YQ(ny_a1p,9)*(T1-TIME+TINCR)/(T1-T0))
            ENDIF

            DELTA_XA=0.0d0
            DO nj=1,NJTOT !calculate delta x for half time step
              DELTA_XA=DELTA_XA+((XQ(nj,nq_a1)-XQ(nj,nq_a2))**2.0d0)
            ENDDO
C            DELTA_XA=DMAX1((DELTA_XA**0.5d0),0.6d0)*
C     '        (LAMBDA_A1+LAMBDA_A2+OLD_LAMBDA_A1+
C     '        OLD_LAMBDA_A2)*0.25d0
            DELTA_XA=(DELTA_XA**0.5d0)*(LAMBDA_A1+LAMBDA_A2
     &        +OLD_LAMBDA_A1+OLD_LAMBDA_A2)*0.25d0

            IF(ITYP12.EQ.2) THEN !pulmonary flow - calculate pleural pressure
              G_PLEURAL=0.d0
              DO nj=1,NJT
                G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &            *GRAVITY*1000.d0*HEIGHT_B1(nj)
              ENDDO
              TRACE_B1=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
            ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
              nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
              TRACE_B1=CQ(nj,nq)/1000.d0 !Pa->kPa
            ELSE
              TRACE_B1=(YQ(ny_b1p,2)*(TIME-T0)/(T1-T0))+
     '          (YQ(ny_b1p,9)*(T1-TIME)/(T1-T0))
            ENDIF

            DELTA_XB=0.0d0
            DO nj=1,NJTOT !calculate delta x for half time step
              DELTA_XB=DELTA_XB+((XQ(nj,nq_b1)-XQ(nj,nq_b2))**2.0d0)
            ENDDO
C           DELTA_XB=DMAX1((DELTA_XB**0.5d0),0.6d0)*
C    '        (LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
C    '        *0.25d0
            DELTA_XB=(DELTA_XB**0.5d0)*
     '        (LAMBDA_B1+LAMBDA_B2+OLD_LAMBDA_B1+OLD_LAMBDA_B2)
     '        *0.25d0

            PA=YQ(ny_a1p,8)-TRACE_A1+OLD_TRACE_A1
! inital starting guess for pa
            PA_OLD=YQ(ny_a1p,8)

            PB=YQ(ny_b1p,8)-TRACE_A1+OLD_TRACE_A1
! inital starting guess for pb, trace values are the same
            PB_OLD=YQ(ny_b1p,8)
            XSUM=0.0d0
            BSUM=0.0d0
            BCONVERGE=1.0d0
            XCONVERGE=1.0d0

            DO WHILE(((XCONVERGE.GT.(CONVERG_TOL*100.0d0)).AND.
     '        (BCONVERGE.GT.(CONVERG_TOL*100.0d0)).AND.
     '        (iterate_count.LT.100))
     '        .OR.(iterate_count.LT.5))

              iterate_count=iterate_count+1

              IF((PA+TRACE_A1).GT.0.0d0) THEN
                Sa=PI*(ROIA1**2.0d0)*
     '            ((((PA+TRACE_A1)/CQ(5,nq_a1))+1)
     '            **(2.0D0/BETAA1))
                dSa=(PI*(ROIA1**2.0d0)*2.0d0/
     '            (BETAA1*CQ(5,nq_a1))
     '            )*((((PA+trace_a1)/CQ(5,nq_a1))+1)
     '            **((2.0D0/BETAA1)-1.0d0))
                Ra=ROIA1*((((PA+TRACE_A1)/CQ(5,nq_a1))+1)
     '            **(1.0D0/BETAA1))
                dRa=(ROIA1/
     '            (BETAA1*CQ(5,nq_a1)))
     '            *((((PA+trace_a1)/CQ(5,nq_a1))+1)
     '            **((1.0D0/BETAA1)-1.0d0))
              ELSE
                fo=BETAA1*CQ(5,nq_a1)/CQ(7,nq_a1)
                Sa=PI*((ROIA1/((1-((PA+TRACE_A1)/
     '            fo))**(1.0d0/CQ(7,nq_a1))))**2.0d0)
                dSa=(2.0D0*PI*(ROIA1**2.0d0)/(CQ(7,nq_a1)*fo))
     '            *(1-(PA+TRACE_A1)/Fo)**(-2.0d0/CQ(7,nq_a1)-1.0d0)
                Ra=(ROIA1/((1-((PA+TRACE_A1)/
     '            fo))**(1.0d0/CQ(7,nq_a1))))
                dRa=(ROIA1/(CQ(7,nq_a1)*fo))
     '            *(1-(PA+TRACE_A1)/Fo)**(-1.0d0/CQ(7,nq_a1)-1.0d0)

              ENDIF

              TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '          ((YQ(ny_a2v,5)**2)/YQ(ny_a2r,5))*
     '          (RA+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))
              TERM2=(1.0D0/CQ(1,nq_a1))*(PA
     '          +YQ(ny_a1p,8)-YQ(ny_a2p,1)-YQ(ny_a2p,8))!-2.d0*G_TERM*DX
!     '          *DSIN(ANGLE))

              IF(YQ(ny_a2r,5).GT.ROIA1h) THEN
                TERM3=2*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)/
     '            (CQ(3,nq_a1)-1.0d0))*
     '            (YQ(ny_a2v,5)/(YQ(ny_a2r,5)**2.0d0))
              ELSE
                RO_ART=ROIA1h
                ELIP_A=((RO_ART**2.0d0)+
     '            (RO_ART**4.0d0-YQ(ny_a2r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_B=((RO_ART**2.0d0)-
     '            (RO_ART**4.0d0-YQ(ny_a2r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '            (YQ(ny_a2r,5)**2.0d0*(ELIP_A**2.0d0+
     '            ELIP_B**2.0d0)))**0.5D0
                TERM3=2*TINCR*CQ(2,nq_a1)*(CQ(3,nq_a1)
     '            /(CQ(3,nq_a1)-1.0d0))*
     '            (YQ(ny_a2v,5)/(ELIP_TERM**2.0d0))
              ENDIF


              TERM4=(1.0d0/YQ(ny_a2r,5))*((-DELTA_XA/TINCR)+
     '          ((2.0d0*CQ(3,nq_a1)-1.0d0)*YQ(ny_a2v,5)))*
     '          (((LAMBDA_A2**0.5d0)*YQ(ny_a2r,1))+
     '          ((LAMBDA_A1**0.5d0)*RA)-
     '          ((OLD_LAMBDA_A2**0.5d0)*YQ(ny_a2r,8))
     '          -((OLD_LAMBDA_A1**0.5d0)*YQ(ny_a1r,8)))

              TERM5=YQ(ny_a2v,5)/YQ(ny_a2r,5)
     '          *(RA+YQ(ny_a1r,8)-YQ(ny_a2r,1)-YQ(ny_a2r,8))

              VA=YQ(ny_a2v,8)+((TINCR/(2.0d0*DELTA_XA))*
     '          (TERM1-TERM2))-TERM3+TERM4-TERM5


              TERM1=((2.0d0*CQ(3,nq_a1))-0.0d0)*
     '          ((YQ(ny_a2v,5)**2)/YQ(ny_a2r,5))*dRa
              TERM2=1.0d0/CQ(1,nq_a1)
              TERM3=((LAMBDA_A1**0.5d0)
     '          /YQ(ny_a2r,5))*((-DELTA_XA/TINCR)+
     '          ((2.0d0*CQ(3,nq_a1)-1.0d0)*YQ(ny_a2v,5)))*dRa
              TERM4=YQ(ny_a2v,5)/YQ(ny_a2r,5)*(-dRa)
              DVA=((TINCR/(2.0d0*DELTA_XA))*
     '          (TERM1-TERM2))+TERM3+TERM4

              dFA=(Sa*DVA)+(dSa*VA)
              FA1=PI*(Ra**2.0d0)*va

              IF((PB+TRACE_B1).GT.0.0d0) THEN
                Sb=PI*(vien_ro_b1**2.0d0)*
     '            ((((PB+TRACE_B1)/CQ(11,nq_b1))+1)
     '            **(2.0D0/BETAB1))
                dSb=(PI*(vien_ro_b1**2.0d0)*2.0d0/
     '            (BETAB1*CQ(11,nq_b1))
     '            )*((((PB+trace_b1)/CQ(11,nq_b1))+1)
     '            **((2.0D0/BETAB1)-1.0d0))
                Rb=vien_ro_b1*((((PB+TRACE_B1)/CQ(11,nq_b1))+1)
     '            **(1.0D0/BETAB1))
                dRb=(vien_ro_b1/
     '            (BETAB1*CQ(11,nq_b1)))
     '            *((((PB+trace_b1)/CQ(11,nq_b1))+1)
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
              ENDIF

              TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '          ((YQ(ny_b2v,5)**2)/YQ(ny_b2r,5))*
     '          (RB+YQ(ny_b1r,8)-YQ(ny_b2r,1)-YQ(ny_b2r,8))
              TERM2=(1.0D0/CQ(1,nq_b1))*(PB
     '          +YQ(ny_b1p,8)-YQ(ny_b2p,1)-YQ(ny_b2p,8))!-2.d0*G_TERM*DX
!             '          *DSIN(ANGLE))

              IF(YQ(ny_b2r,5).GT.ROIB1h) THEN
                TERM3=2*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)/
     '            (CQ(3,nq_b1)-1.0d0))*
     '            (YQ(ny_b2v,5)/(YQ(ny_b2r,5)**2.0d0))
              ELSE
                RO_VIEN=ROIB1h
                ELIP_A=((RO_VIEN**2.0d0)+
     '            (RO_VIEN**4.0d0-YQ(ny_b2r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_B=((RO_VIEN**2.0d0)-
     '            (RO_VIEN**4.0d0-YQ(ny_b2r,5)**4.0d0)
     '            **0.5d0)**0.5d0
                ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '            (YQ(ny_b2r,5)**2.0d0*(ELIP_A**2.0d0+
     '            ELIP_B**2.0d0)))**0.5D0
                TERM3=2*TINCR*CQ(2,nq_b1)*(CQ(3,nq_b1)
     '            /(CQ(3,nq_b1)-1.0d0))*
     '            (YQ(ny_b2v,5)/(ELIP_TERM**2.0d0))
              ENDIF

              TERM4=(1.0d0/YQ(ny_b2r,5))*((-DELTA_XB/TINCR)+
     '          ((2.0d0*CQ(3,nq_b1)-1.0d0)*YQ(ny_b2v,5)))*
     '          (((LAMBDA_B2**0.5d0)*YQ(ny_b2r,1))+
     '          ((LAMBDA_B1**0.5d0)*RB)-
     '          ((OLD_LAMBDA_B2**0.5d0)*YQ(ny_b2r,8))
     '          -((OLD_LAMBDA_B1**0.5d0)*YQ(ny_b1r,8)))

              TERM5=YQ(ny_b2v,5)/YQ(ny_b2r,5)
     '          *(RB+YQ(ny_b1r,8)-YQ(ny_b2r,1)-YQ(ny_b2r,8))

              VB=YQ(ny_b2v,8)+((TINCR/(2.0d0*DELTA_XB))*
     '          (TERM1-TERM2))-TERM3+TERM4-TERM5


              TERM1=((2.0d0*CQ(3,nq_b1))-0.0d0)*
     '          ((YQ(ny_b2v,5)**2)/YQ(ny_b2r,5))*dRb
              TERM2=1.0d0/CQ(1,nq_b1)
              TERM3=((LAMBDA_B1**0.5d0)
     '          /YQ(ny_b2r,5))*((-DELTA_XB/TINCR)+
     '          ((2.0d0*CQ(3,nq_b1)-1.0d0)*YQ(ny_b2v,5)))*dRb
              TERM4=YQ(ny_b2v,5)/YQ(ny_b2r,5)*(-dRb)
              DVB=((TINCR/(2.0d0*DELTA_XB))*
     '          (TERM1-TERM2))+TERM3+TERM4

              dFB=(Sb*DVB)+(dSb*VB)
              FB1=PI*(Rb**2.0d0)*vb

              dFB=-dFB
              FB1=-FB1

              FC_OLD=YQ(ny_a1v,4)

              CALL USER3_CORONARY1(DENOM_BOX(1,1),DRART,
     '          (dFA/1000.0d0),NUM_BOX(1,1),RART,(PA+TRACE_A1),
     '          (FA1/1000.0d0),ERROR,*9999)

              RART=(RART/1000.0d0)
     '          *(LAMBDA_A1**2.0D0) ! convert from mls to mm^3
              DRART=(DRART/1000.0d0)
     '          *(LAMBDA_A1**2.0D0) !and account for length change

              RART_OLD=YQ(ny_a1r,4)

              CALL USER3_CORONARY1(DENOM_BOX(1,3),DRVIE,
     '          (-dFB/1000.0d0),
     '          NUM_BOX(1,3),RVIE,(PB+TRACE_B1),(-FB1/1000.0d0),
     '          ERROR,*9999)

              RVIE= (RVIE/1000.0d0)
     '          *(LAMBDA_B1**2.0D0) ! convert from mls to mm^3
              DRVIE=(DRVIE/1000.0d0)
     '          *(LAMBDA_B1**2.0D0) ! and account for length change

              RVIE_OLD=YQ(ny_b1r,4)

              P1=PA-FA1*RART
              P2=PB+FB1*RVIE
              DP1DPA=1-((DFA*RART)+(FA1*DRART))
              DP2DPV=1+((DFB*RVIE)+(FB1*DRVIE))

              CALL USER3_CORONARY2(DENOM_BOX(1,2),
     '          NUM_BOX(1,2),DRCAP_PA,DRCAP_PB,
     '          RCAP,(P1+TRACE_A1),(P2+TRACE_B1),ERROR,*9999)

              RCAP=(RCAP/1000.0D0)*(LAMBDA_B1**2.0D0)
C convert from mls to mm^3
              DRCAP_PA=((DRCAP_PA*DP1DPA)/1000.0D0)
     '          *(LAMBDA_B1**2.0D0)
C convert from mls to mm^3 and change derivative from drcdp1 to drcdpa
C also now take imto account length change
              DRCAP_PB=((DRCAP_PB*DP2DPV)/1000.0D0)
     '          *(LAMBDA_B1**2.0D0)
C as above

              CALL USER3_CORONARY1(DENOM_BOX(1,4),DCOM1,
     '          DP1DPA,NUM_BOX(1,4),COM1,(PA+TRACE_A1)
     '          ,(P1+TRACE_A1),ERROR,*9999)

              COM1=COM1*1000.0D0 !convert from mls to mm^3
              DCOM1=DCOM1*1000.0D0 !convert from mls to mm^3

              COM1_OLD=YQ(ny_a1p,4)

              IF(COM1.lt.0.0001D0) THEN
                COM1=0.0001D0
                DCOM1=0.0D0
              ENDIF

              CALL USER3_CORONARY1(DENOM_BOX(1,5),DCOM2,
     '          DP2DPV,NUM_BOX(1,5),COM2,
     '          (PB+TRACE_B1),(P2+TRACE_B1),ERROR,*9999)

              COM2=COM2*1000.0D0 !convert from mls to mm^3
              DCOM2=DCOM2*1000.0D0 !convert from mls to mm^3

              COM2_OLD=YQ(ny_B1p,4)

              IF(COM2.lt.0.0001D0) THEN
                COM2=0.0001D0
                DCOM2=0.0D0
              ENDIF

              DFCPA=((1.0d0-FA1*DRART-RART*DFA)*RCAP-
     '          (PA-FA1*RART-FB1*RVIE-PB)*DRCAP_PA)/(RCAP**2.0d0)

              DFCPB=((-1.0d0-FB1*DRVIE-RVIE*DFB)*RCAP-
     '          (PA-FA1*RART-FB1*RVIE-PB)*DRCAP_PB)/(RCAP**2.0d0)

              FC=(PA-FA1*RART-FB1*RVIE-PB)/RCAP

              BSUM_OLD=BSUM

              B(1)=FA1*RART-PA+(FA1+FA1_OLD-FC-FC_OLD)*
     '          (1.0D0*TINCR/(COM1+COM1_OLD))+
     '          (PA_OLD-FA1_OLD*RART_OLD)-TRACE_A1+OLD_TRACE_A1

              B(2)=FB1*RVIE+PB-(FC+FC_OLD-FB1-FB1_OLD)*
     '          (1.0D0*TINCR/(COM2+COM2_OLD))-
     '          (PB_OLD+FB1_OLD*RVIE_OLD)+TRACE_A1-OLD_TRACE_A1

              BSUM=DABS(B(1))+DABS(B(2))

              X(1)=B(1)
              X(2)=B(2)

              A(1,1)=-1.0d0+DFA*RART+DRART*FA1+(DFA-DFCPA)*
     '          (1.0d0*TINCR/(COM1+COM1_OLD))+
     '          (FA1+FA1_OLD-FC-FC_OLD)*(-1.0D0*TINCR)/
     '          ((COM1+COM1_OLD)**2.0d0)*DCOM1

              A(1,2)=(-1.0d0*TINCR)/(COM1+COM1_OLD)*DFCPB

              A(2,1)=(-1.0d0*TINCR)/(COM2+COM2_OLD)*DFCPA

              A(2,2)=1.0d0+DFB*RVIE+DRVIE*FB1-(DFCPB-DFB)*
     '          (1.0d0*TINCR/(COM2+COM2_OLD))-
     '          (FC+FC_OLD-FA1-FA1_OLD)*(-1.0D0*TINCR)/
     '          ((COM2+COM2_OLD)**2.0d0)*DCOM2

              X(1)=(B(2)-(A(2,2)*B(1)/A(1,2)))/
     '          (A(2,1)-(A(1,1)*A(2,2)/A(1,2)))
              X(2)=(B(1)-(A(1,1)*X(1)))/A(1,2)

c              CALL DGETRF(2,2,A,2,IPIV,INFO)
c              CALL DGETRS(TRANS,2,1,A,2,IPIV,X,2,INFO)

              PA=PA-X(1)
              PB=PB-X(2)
              XSUM_OLD=XSUM
              XSUM=DABS(x(1))+DABS(x(2))
              XCONVERGE=DABS(XSUM_OLD-XSUM)/(DABS(XSUM)+1.0D0)
              BCONVERGE=DABS(BSUM_OLD-BSUM)/(DABS(BSUM)+1.0D0)

            ENDDO

            YQ(ny_a1p,1)=PA
            YQ(ny_B1p,1)=PB
            YQ(ny_a1R,1)=RA
            YQ(ny_B1R,1)=RB
            YQ(ny_a1V,1)=VA
            YQ(ny_B1V,1)=VB

C storing lumped parameter values for the next time step

            YQ(ny_a1p,4)=COM1
            YQ(ny_a1r,4)=RART
            YQ(ny_a1v,4)=FC
            YQ(ny_B1p,4)=COM2
            YQ(ny_b1r,4)=RVIE

          ELSE !no identical venous network OR no venous network

            nq_previous=CONECT(-1,1,nq)
C            nq_previous=NXQ(-1,1,nq,1)

            ny_p=NYNQ(1,nq,0)
            ny_r=NYNQ(2,nq,0)
            ny_v=NYNQ(3,nq,0)
            ny_p2=NYNQ(1,nq_previous,0)
            ny_r2=NYNQ(2,nq_previous,0)
            ny_v2=NYNQ(3,nq_previous,0)

            LAMBDA_1OLD=(YQ(ny_v,2)*(TIME-T0-TINCR)/(T1-T0))+
     '        (YQ(ny_v,9)*(T1-TIME+TINCR)/(T1-T0))
            LAMBDA_1=(YQ(ny_v,2)*(TIME-T0)/(T1-T0))+
     '        (YQ(ny_v,9)*(T1-TIME)/(T1-T0))

            LAMBDA_2OLD=(YQ(NYNQ(3,nq_previous,0),2)
     '        *(TIME-T0-TINCR)/(T1-T0))+
     '        (YQ(NYNQ(3,nq_previous,0),9)*(T1-TIME+TINCR)/(T1-T0))
            LAMBDA_2=(YQ(NYNQ(3,nq_previous,0),2)*
     '        (TIME-T0)/(T1-T0))+
     '        (YQ(NYNQ(3,nq_previous,0),9)*(T1-TIME)/(T1-T0))

            DELTA_X=0.0d0
            DO nj=1,NJTOT !calculate delta x for half time step
              DELTA_X=DELTA_X+
     '          ((XQ(nj,nq)-XQ(nj,nq_previous))**2.0d0)
            ENDDO
C            DELTA_X=DMAX1((DELTA_X**0.5d0),0.60d0)*
C     '        ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)*0.25D0)
            DELTA_X=(DELTA_X**0.5d0)*
     '        ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)*0.25D0)

            ROI=CQ(4,nq)/(LAMBDA_1**0.5d0)
            ROIh=(CQ(4,nq)+CQ(4,nq))*0.5d0/
     '        ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)*0.25d0)

            BETA=CQ(6,nq)*LAMBDA_1+CQ(12,nq)
            BETAh=CQ(6,nq)*
     '        ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)
     '        *0.25d0)+CQ(12,nq)

            BETA=DMAX1(BETA,1.0d0)
            BETAh=DMAX1(BETAh,1.0d0)

C... KSB Oct-03 : Calculating gravity term
            DX=0.d0
            DO nj=1,NJT
              DX=DX+(XQ(nj,nq_previous)-XQ(nj,nq))**2.d0
              HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq)
              DELTA(nj)=XQ(nj,nq_previous)-XQ(nj,nq)
            ENDDO
            DX=DSQRT(DX)
            G_BC=0.d0
            G_TERM=0.d0
            DO nj=1,NJT
              G_BC=G_BC+CQ(1,nq)*G_VECTOR(nj)*GRAVITY*1000.d0*HEIGHT(nj)
              G_TERM=G_TERM+CQ(1,nq)*G_VECTOR(nj)*GRAVITY*1000.d0
     &          *DELTA(nj)
            ENDDO

C... Calculating external pressure force            
            IF(ITYP12.EQ.2) THEN !pulmonary flow
              G_PLEURAL=0.d0 !gravitational force - initialise
              PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
 !This density value effectively represents lung tissue density = 1/5th that of water
              DO nj=1,NJT
                G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &            *GRAVITY*1000.d0*HEIGHT(nj)
              ENDDO
              TRACE=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
              !NB/ currently not time-dependent
            ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
              nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
              TRACE=CQ(nj,nq)/1000.d0 !Pa->kPa
            ELSE              
              TRACE=(YQ(ny_p,2)*(TIME-T0)/(T1-T0))+
     '          (YQ(ny_p,9)*(T1-TIME)/(T1-T0))
            ENDIF
            
C PM 29-NOV-01: Boundary conditions from iptime or file

            IF(KTYP3_INIT(nr).EQ.3) THEN         ! from a file
              READ(IOFILE1,*) YQ(ny_p,1)
            ELSEIF(KTYP3_INIT(nr).EQ.4) THEN     ! defined as time variable
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
              ntv1=NTIME_NR(1,nr)
              ntv=NTIME_NR(2,nr) !KSB: region dependent boundary conditions
C... KSB replacing hard coded ntv values with actual region dep values 1/12/05
              IF((PERIODIC.EQ.'N').OR.(PERIODIC.EQ.'n')) THEN
                IF(TIME.LT.TIME_VALUES(1,1,ntv)) THEN
                  YQ(ny_p,1)=TIME_VALUES(2,0,ntv)+G_BC
                ELSEIF (TIME.GT.TIME_VALUES(1,NTIME_POINTS(ntv1),ntv))
     &              THEN
                  YQ(ny_p,1)=TIME_VALUES(2,NTIME_POINTS(ntv1)+1,ntv)
     &              +G_BC
                ELSE
                  DO ntp=1,NTIME_POINTS(ntv)-1
                    IF((TIME.GE.TIME_VALUES(1,1,ntv)+T_period(ntv)*
     '                T_frac(ntp,ntv)).AND.(TIME.LT.TIME_VALUES(1,1,ntv)
     &                +T_period(ntv)*T_frac(ntp+1,ntv))) THEN
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
     '              (TIME-NT_CYCLE(ntv)*T_period(ntv))+G_BC
                  ENDIF
                ENDDO ! ntp
              ENDIF  ! PERIODIC
            ELSE
            ENDIF  ! KTYP_INIT
            
C            IF(TIME.LT.((T1+T0)/1.0d0)) THEN
C              YQ(ny_p,1)=YQ(ny_p,3)+EXIT_PRESSURE0+
C     '          ((EXIT_PRESSURE1-EXIT_PRESSURE0)/(T1-T0))*
C     '          ((TIME-T0)*1.0d0)
C            ELSE
C              YQ(ny_p,1)=YQ(ny_p,3)+EXIT_PRESSURE1
C            ENDIF

            YQ(NYNQ(2,nq,0),4)=0
            YQ(NYNQ(3,nq,0),4)=0

            IF((YQ(ny_p,1)+TRACE).GT.0) THEN
              YQ(ny_r,1)=ROI*((((YQ(ny_p,1)+TRACE)
     '          /CQ(5,nq))+1.0d0)**(1.0d0/BETA))
            ELSE
              fo=BETA*CQ(5,nq)/CQ(7,nq)
              YQ(ny_r,1)=ROI*(1-((YQ(ny_p,1)+TRACE)/fo))
     '          **(-1.0d0/CQ(7,nq))
            ENDIF
            TERM1=((2.0d0*CQ(3,nq))-0.0d0)*
     '        ((YQ(ny_v2,5)**2)/YQ(ny_r2,5))*
     '        (YQ(ny_r,1)+YQ(ny_r,8)-YQ(ny_r2,1)-YQ(ny_r2,8))
            TERM2=(1.0D0/CQ(1,nq))*(YQ(ny_p,1)
     '        +YQ(ny_p,8)-YQ(ny_p2,1)-YQ(ny_p2,8)-2.d0*G_TERM)

            IF(YQ(ny_r2,5).GT.ROIh) THEN

              TERM3=2*TINCR*CQ(2,nq)*(CQ(3,nq)/(CQ(3,nq)-1.0d0))*
     '          (YQ(ny_v2,5)/(YQ(ny_r2,5)**2.0d0))
            ELSE
              RO_ART=ROIh
              ELIP_A=((RO_ART**2.0d0)+
     '          (RO_ART**4.0d0-YQ(ny_r2,5)**4.0d0)
     '          **0.5d0)**0.5d0
              ELIP_B=((RO_ART**2.0d0)-
     '          (RO_ART**4.0d0-YQ(ny_r2,5)**4.0d0)
     '          **0.5d0)**0.5d0


              ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '          (YQ(ny_r2,5)**2.0d0*(ELIP_A**2.0d0+
     '          ELIP_B**2.0d0)))**0.5D0
              TERM3=2*TINCR*CQ(2,nq)*(CQ(3,nq)
     '          /(CQ(3,nq)-1.0d0))*
     '          (YQ(ny_v2,5)/(ELIP_TERM**2.0d0))
            ENDIF

            TERM4=(1.0d0/YQ(ny_r2,5))*((-DELTA_X/TINCR)+
     '        ((2.0d0*CQ(3,nq)-1.0d0)*YQ(ny_v2,5)))*
     '        (((LAMBDA_2**0.5d0)*YQ(ny_r2,1))+
     '        ((LAMBDA_1**0.5d0)*YQ(ny_r,1))-
     '        ((LAMBDA_2OLD**0.5d0)*YQ(ny_r2,8))
     '        -((LAMBDA_1OLD**0.5d0)*YQ(ny_r,8)))

            TERM5=YQ(ny_v2,5)/YQ(ny_r2,5)
     '        *(YQ(ny_r,1)+YQ(ny_r,8)-YQ(ny_r2,1)-YQ(ny_r2,8))

            YQ(ny_v,1)=YQ(ny_v2,8)+((TINCR/(2.0d0*DELTA_X))*
     '        (TERM1-TERM2))-TERM3+TERM4-TERM5

            IF(((VENOUS_NETWORK.EQ.'Y').OR.(VENOUS_NETWORK.EQ.'y'))
     '        .AND.(N_VENOUS_GEOM.EQ.1)) THEN
              !venous network with identical geometry

              ny_p=NYNQ(4,nq,0)
              ny_r=NYNQ(5,nq,0)
              ny_v=NYNQ(6,nq,0)
              ny_p2=NYNQ(4,nq_previous,0)
              ny_r2=NYNQ(5,nq_previous,0)
              ny_v2=NYNQ(6,nq_previous,0)

              ROI=CQ(4,nq)/(LAMBDA_1**0.5d0)
              ROIh=(CQ(4,nq)+CQ(4,nq))*0.5d0/
     '          ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)*0.25d0)

              BETA=CQ(9,nq)*LAMBDA_1+CQ(10,nq)
              BETAh=CQ(9,nq)*
     '          ((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)
     '          *0.25d0)+CQ(10,nq)

              BETA=DMAX1(BETA,1.0d0)
              BETAh=DMAX1(BETAh,1.0d0)

              IF(ITYP12.EQ.2) THEN !pulmonary flow
                DO nj=1,NJT
                  HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq)
                ENDDO
                G_PLEURAL=0.d0 !gravitational force - initialise
                PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
 !This density value effectively represents lung tissue density = 1/5th that of water
                DO nj=1,NJT
                  G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &              *GRAVITY*1000.d0*HEIGHT(nj)
                ENDDO
                TRACE=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
                !NB/ currently not time-dependent                
              ELSEIF(ITYP12.EQ.3) THEN !external pressure field read in
                nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
                TRACE=CQ(nj,nq)/1000.d0 !Pa->kPa
              ELSE !all other flow types
                TRACE=(YQ(ny_p,2)*(TIME-T0)/(T1-T0))+
     '            (YQ(ny_p,9)*(T1-TIME)/(T1-T0))
              ENDIF

              IF(TIME.LT.((T1+T0)/1.0d0)) THEN
                YQ(ny_p,1)=YQ(ny_p,3)+EXIT_PRESSURE0+
     '            ((EXIT_PRESSURE1-EXIT_PRESSURE0)/(T1-T0))*
     '            ((TIME-T0)*1.0d0)+G_BC
              ELSE
                YQ(ny_p,1)=YQ(ny_p,3)+EXIT_PRESSURE1+G_BC
              ENDIF

              YQ(NYNQ(2,nq,0),4)=0
              YQ(NYNQ(3,nq,0),4)=0

              IF((YQ(ny_p,1)+TRACE).GT.0) THEN
                YQ(ny_r,1)=ROI*VIEN_RATIO*((((YQ(ny_p,1)+TRACE)
     '            /CQ(11,nq))+1.0d0)**(1.0d0/BETA))
              ELSE
                fo=BETA*CQ(11,nq)/CQ(7,nq)
                YQ(ny_r,1)=ROI*VIEN_RATIO*
     '            (1-((YQ(ny_p,1)+TRACE)/fo))
     '            **(-1.0d0/CQ(7,nq))
              ENDIF
              TERM1=((2.0d0*CQ(3,nq))-0.0d0)*
     '          ((YQ(ny_v2,5)**2)/YQ(ny_r2,5))*
     '          (YQ(ny_r,1)+YQ(ny_r,8)-YQ(ny_r2,1)-YQ(ny_r2,8))
              TERM2=(1.0D0/CQ(1,nq))*(YQ(ny_p,1)
     '          +YQ(ny_p,8)-YQ(ny_p2,1)-YQ(ny_p2,8)-2.d0*G_TERM)

              IF(YQ(ny_r2,5).GT.
     '          (ROIh*VIEN_RATIO)) THEN

                TERM3=2*TINCR*CQ(2,nq)*(CQ(3,nq)/(CQ(3,nq)-1.0d0))*
     '            (YQ(ny_v2,5)/(YQ(ny_r2,5)**2.0d0))
              ELSE
                RO_VIEN=ROIh*VIEN_RATIO
                ELIP_A=(RO_VIEN**2.0d0+
     '            (RO_VIEN**4.0d0-YQ(ny_r2,5)**4.0d0)
     '            **0.5d0)**0.5d0

                ELIP_B=(RO_VIEN**2.0d0-
     '            (RO_VIEN**4.0d0-YQ(ny_r2,5)**4.0d0)
     '            **0.5d0)**0.5d0

                ELIP_TERM=((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '            (YQ(ny_r2,5)**2.0d0*(ELIP_A**2.0d0+
     '            ELIP_B**2.0d0)))**0.5D0
                TERM3=2*TINCR*CQ(2,nq)*(CQ(3,nq)/(CQ(3,nq)-1.0d0))*
     '            (YQ(ny_v2,5)/(ELIP_TERM**2.0d0))
              ENDIF

              TERM4=(1.0d0/YQ(ny_r2,5))*((-DELTA_X/TINCR)+
     '          ((2.0d0*CQ(3,nq)-1.0d0)*YQ(ny_v2,5)))*
     '          (((LAMBDA_2**0.5d0)*YQ(ny_r2,1))+
     '          ((LAMBDA_1**0.5d0)*YQ(ny_r,1))-
     '          ((LAMBDA_2OLD**0.5d0)*YQ(ny_r2,8))
     '          -((LAMBDA_1OLD**0.5d0)*YQ(ny_r,8)))

              TERM5=YQ(ny_v2,5)/YQ(ny_r2,5)
     '          *(YQ(ny_r,1)+YQ(ny_r,8)-YQ(ny_r2,1)-YQ(ny_r2,8))
              YQ(ny_v,1)=YQ(ny_v2,8)+((TINCR/(2.0d0*DELTA_X))*
     '          (TERM1-TERM2))-TERM3+TERM4-TERM5

            ENDIF  ! venous network with identical geometry

          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('BRANCH2')
      RETURN
 9999 CALL ERRORS('BRANCH2',ERROR)
      CALL EXITS('BRANCH2')
      RETURN 1
      END


