      SUBROUTINE XPFD30(NXQ,N_TERM_Q,CONECT,CQ,ES,HALF_TIME_STEP,NHQ,nq,
     '  nr,nx,NYNQ,UPDATE_MATRIX,UPDATE_VECTOR,TIME,XQ,YQ,ERROR,*)


C#### Subroutine: XPFD30
C###  Description:
C###    XPFD30 calculates difference grid point matrices
C###    (ES(nhs1,nhs2))
C****  NPS 2/11/96 Currently only implemented for Lax Wendroff
C****  TDS 13.07.06 Changed all x**0.5 terms to DSQRT(x) and y=x**2 terms
C****  to tmp=x; y=tmp*tmp or y=x*x depending on the original usage.
C****  Also formatted the original code slightly and added some new
C****  comments. Introduced the new constants - density, viscosity and
C****  alpha so that we don't have to use CQ(xx,nq) where required.
      
      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER  CONECT(-1:1,0:2,NQM),NHQ(NRM),N_TERM_Q(-1:1,0:NEM),
     '  nq,nr,nx,NYNQ(NHM,NQM,0:NRCM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 CQ(NMM,NQM),ES(NHM*NSM,NHM*NSM),TIME,
     '  XQ(NJM,NQM),YQ(NYQM,NIQM)
      CHARACTER ERROR*(*)
      LOGICAL HALF_TIME_STEP,UPDATE_MATRIX,UPDATE_VECTOR

!     Local Variables
      INTEGER i,nhs1,nhs2,nj,NJTOT,nonq,nq_next,nq_term,nq_pre,
     '  ny_p1,ny_p2,ny_r1,ny_r2,ny_v1,ny_v2,VIEN_OFFSET,
     '  yp_index,yp_value,NET_WORK_NO
      REAL*8 BETA_A,BETA_NEXT_A,BETA_V,BETA_NEXT_V,DELTA(NJT),
     &  DELTA_X,DELTA_X1,DELTA_X2,DX,ELIP_A,ELIP_B,ELIP_TERM,fo1,fo2,
     &  G_PLEURAL,G_TERM,HEIGHT(NJT),LAMBDA_1OLD,LAMBDA_1,LAMBDA_2OLD,
     &  LAMBDA_2,LAMBDA_3OLD,LAMBDA_3,LAMBDA_1N,LAMBDA_2N,LAMBDA_1H,
     &  LAMBDA_2H,P1,P2,P_DIFF,PLEURAL_DENSITY,R1,R2,RO,ROI,ROI_NEXT,
     &  ROI_PRE,TERM1,TERM2,TERM3,TERM4,TERM5,test1,test2,V1,V2,ALPHA,
     &  DENSITY,VISCOSITY,TMP
      LOGICAL TERMINAL

      CALL ENTERS('XPFD30',*9999)

      IF(HALF_TIME_STEP) THEN
        yp_index=1 !previous time step
        yp_value=5 !previous time step
      ELSE
        yp_index=5 !half time step
        yp_value=1 !half time step
      ENDIF

      NJTOT=NJ_LOC(NJL_GEOM,0,nr)

      DENSITY=CQ(1,nq)
      VISCOSITY=CQ(2,nq)
      ALPHA=CQ(3,nq)
      
      IF (UPDATE_MATRIX) THEN
        IF(ITYP16(nr,nx).eq.1) THEN
C ES is the idenity for explicit schemes
          DO nhs1=1,NHQ(nr)
            DO nhs2=1,NHQ(nr)
              IF(nhs1.eq.nhs2) THEN
                ES(nhs1,nhs2)=1.0d0
              ELSE
                ES(nhs1,nhs2)=0.0d0
              ENDIF
            ENDDO
          ENDDO
        ELSE
C implicit  schemes
        ENDIF
      ENDIF
      
C PM 13-12-01: Transition and neighbouring point
      TERMINAL=.FALSE.
      
      IF(((VENOUS_NETWORK.EQ.'Y').OR.(VENOUS_NETWORK.EQ.'y'))
     '  .AND.(N_VENOUS_GEOM.EQ.2)) THEN
        DO nonq=1,N_TERM_Q(0,0)
          nq_term=N_TERM_Q(0,nonq)
          IF(nq_term.EQ.nq.OR.
     &      NXQ(-1,1,nq_term,1).EQ.nq.OR.
     &      NXQ(1,1,nq_term,1).EQ.nq.OR.
     &      NXQ(-1,1,NXQ(-1,1,nq_term,1),1).EQ.nq) THEN
            TERMINAL=.TRUE.
          ENDIF
        ENDDO
      ENDIF

      IF (UPDATE_VECTOR) THEN
        IF (ITYP16(nr,nx).eq.4) THEN !Lax-Wendroff
          IF (ITYP3(nr,nx).eq.1) THEN !flow in elastic tube
            VIEN_OFFSET=3
            IF (((CONECT(-1,0,nq).EQ.1).OR.
     &        HALF_TIME_STEP).AND.(CONECT(1,0,nq).EQ.1)) THEN
              
C The grid point nq is not an end point or a bifrucation
              IF (HALF_TIME_STEP) THEN !Lax-Wendroff 1/2 step
c not endpoint-need to also exclude
C points around bifurcation

C PM 13-12-01 : CONECT(1,1,nq) is supposed to give the grid point number one
C point downstream of nq. It seems that it does not do it correctly for some
C elements. Thus NXQ is used to do this for the grid points in the vicinity of
C the transition points. For details how the transition points are handled see
C BRANCH3. However, NXQ cannot be used for all grid points as two of the three
C grid points in a bifurcation is defined by CONECT. Needs to do a comprehensive
C check on CONECT.
                IF(TERMINAL) THEN
                  nq_next=NXQ(1,1,nq,1)
                ELSE
                  nq_next=CONECT(1,1,nq)
                ENDIF

                CALL ASSERT(NIQM.GE.9,'>>Increase NIQM, must be >= 9',
     '            ERROR,*9999)
                
                LAMBDA_1=(YQ(NYNQ(3,nq,0),2)*
     '            (TIME-T0-TINCR/2.0d0)/(T1-T0))+
     '            (YQ(NYNQ(3,nq,0),9)*(T1-TIME+TINCR/2.0d0)/(T1-T0))
                LAMBDA_1N=(YQ(NYNQ(3,nq,0),2)*(TIME-T0+TINCR/2.0d0)
     '            /(T1-T0))+
     '            (YQ(NYNQ(3,nq,0),9)*(T1-TIME-TINCR/2.0d0)/(T1-T0))
                
                LAMBDA_2=(YQ(NYNQ(3,nq_next,0),2)
     '            *(TIME-T0-TINCR/2.0d0)/(T1-T0))+
     '            (YQ(NYNQ(3,nq_next,0),9)*(T1-TIME+TINCR/2.0d0)
     '            /(T1-T0))
                LAMBDA_2N=(YQ(NYNQ(3,nq_next,0),2)
     '            *(TIME-T0+TINCR/2.0d0)/(T1-T0))+
     '            (YQ(NYNQ(3,nq_next,0),9)*(T1-TIME-TINCR/2.0d0)
     '            /(T1-T0))
                
                LAMBDA_1H=(YQ(NYNQ(3,nq,0),2)*(TIME-T0)
     '            /(T1-T0))+(YQ(NYNQ(3,nq,0),9)*
     '            (T1-TIME)/(T1-T0))
                LAMBDA_2H=(YQ(NYNQ(3,nq_next,0),2)
     '            *(TIME-T0)/(T1-T0))+(YQ(NYNQ(3,nq_next,0),9)*(T1-TIME)
     '            /(T1-T0))

                DELTA_X=0.0d0
                DO nj=1,NJTOT !calculate delta x for half time step
                  TMP=XQ(nj,nq)-XQ(nj,nq_next)
                  DELTA_X=DELTA_X+(TMP*TMP)
                ENDDO
C calculates spatial step
                DELTA_X=DSQRT(DELTA_X)*((LAMBDA_1+LAMBDA_2)/2.0d0)

C PM 29-NOV-01 : Find whether computaion on venous network is required.
                IF(((VENOUS_NETWORK.EQ.'Y').OR.
     '            (VENOUS_NETWORK.EQ.'y')).AND.(N_VENOUS_GEOM.EQ.1))
     '            THEN
                  NET_WORK_NO=1
                ELSE
                  NET_WORK_NO=0
                ENDIF

                DO I=0,NET_WORK_NO ! calculates viens and arteries
                  ROI=CQ(4,nq)/DSQRT(LAMBDA_1)
! calculated at pre full time step
                  ROI_NEXT=CQ(4,nq_next)/DSQRT(LAMBDA_2)
                  
! for the half time step                  
                  ny_p1=NYNQ(1+(I*VIEN_OFFSET),nq,0)
                  ny_r1=NYNQ(2+(I*VIEN_OFFSET),nq,0)
                  ny_v1=NYNQ(3+(I*VIEN_OFFSET),nq,0)
                  ny_p2=NYNQ(1+(I*VIEN_OFFSET),nq_next,0)
                  ny_r2=NYNQ(2+(I*VIEN_OFFSET),nq_next,0)
                  ny_v2=NYNQ(3+(I*VIEN_OFFSET),nq_next,0)
                
C finds ny's for two grid points to calculate 1/2 time step
                  V1=YQ(ny_v1,yp_index)
                  V2=YQ(ny_v2,yp_index)
                  R1=YQ(ny_r1,yp_index)
                  R2=YQ(ny_r2,yp_index)
                  P1=YQ(ny_p1,yp_index)
                  P2=YQ(ny_p2,yp_index)
                  
C radius at 1/2 time step
                  TERM1=(R2+R1)/2.0d0
                  TERM2=(R2+R1)*(V2-V1)/4.0d0
                  TERM3=(V2+V1)*(R2-R1)/2.0d0
                  TERM4=1.0d0/
     '              (DSQRT((LAMBDA_1N+LAMBDA_2N+LAMBDA_1+LAMBDA_2)
     &              /4.0d0))
                  TERM5=DSQRT((LAMBDA_1+LAMBDA_2)/2.0D0)
                  YQ(ny_r1,yp_value)=TERM4*(TERM5*TERM1-
     '              ((TINCR/(2.d0*DELTA_X))*(TERM2+TERM3)))

C velocity at 1/2 time step
                  TERM1=(V1+V2)/2.0d0
                  TERM2=(((2.d0*ALPHA)-1.0d0)/2.0d0)*(V2+V1)*(V2-V1)
                  TERM3=(ALPHA-1.0d0)*
     '              (((V2+V1)*(V2+V1))/(R2+R1)) *(R2-R1)
                  
C KSB adding gravity term to TERM4 (velocity)
                  DX=0.d0
                  DO nj=1,NJT
                    TMP=XQ(nj,nq_next)-XQ(nj,nq)
                    DX=DX+TMP*TMP
                    DELTA(nj)=TMP
                  ENDDO
                  DX=DSQRT(DX) !distance between grid points
                  G_TERM=0.d0
                  DO nj=1,NJT
                    G_TERM=G_TERM+DX*DENSITY*G_VECTOR(nj)*GRAVITY
     &                *1000.d0*(DELTA(nj)/DX)
                  ENDDO
                  P_DIFF=P2-P1+G_TERM                  
                  TERM4=(1.0d0/DENSITY)*P_DIFF
                  IF((R2+R1).GE.(ROI+ROI_NEXT)) THEN
                    TERM5=2.0d0*TINCR*VISCOSITY*(ALPHA/(ALPHA-1))*
     '                ((V2+V1)/((R2+R1)*(R2+R1)))
                  ELSE
                    IF(I.EQ.0) THEN !artery
                      RO=(ROI+ROI_NEXT)/2.0D0
                    ELSE !vien thus Ro=CQ(4,nq)*sqrt(1.5)
                      RO=((ROI+ROI_NEXT)/2.0d0)*VIEN_RATIO
                    ENDIF
                    ELIP_A=DSQRT((RO*RO)+
     '                DSQRT(RO**4.0d0-((R2+R1)/2.0d0)**4.0d0))
                    ELIP_B=DSQRT((RO*RO)-
     '                DSQRT(RO**4.0d0-((R2+R1)/2.0d0)**4.0d0))
                    ELIP_TERM=DSQRT((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '                (((R2+R1)/2.0d0)*((R2+R1)/2.0d0)*(ELIP_A*ELIP_A+
     '                ELIP_B*ELIP_B)))
                    TERM5=2.0d0*TINCR*VISCOSITY*
     '                ((ALPHA/(ALPHA-1))*
     '                ((V2+V1)/((2.0d0*ELIP_TERM)*(2.0d0*ELIP_TERM))))
                  ENDIF
                  YQ(ny_v1,yp_value)=TERM1-((TINCR/(2.d0*DELTA_X))*
     '              (TERM2+TERM3+TERM4))-TERM5                  

C pressure at 1/2 time step
                  IF(ITYP12(nr,nx).EQ.2) THEN !pulmonary flow
C.. Also calculating pleural pressure (including gravity effect) for
C.. TERM4 in pressure calculation below                  
                    DO nj=1,NJT
                      HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq)
                    ENDDO
                    G_PLEURAL=0.d0 !gravitational force
                    PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
                    !This density value effectively represents lung tissue density = 1/4th that of water
                    DO nj=1,NJT
                      G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &                  *GRAVITY*1000.d0*HEIGHT(nj)
                    ENDDO
                    TERM4=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
                  ELSEIF(ITYP12(nr,nx).EQ.3) THEN !external pressure field read in
                    !linearly interpolate nodal pressure field to grid points
                    nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
                    TERM4=CQ(nj,nq)/1000.d0 !Pa->kPa
                  ELSE !coronary or other flow
                    TERM4=(0.5d0*(YQ(ny_p1,2)+YQ(ny_p2,2)))
     '                *((TIME-T0)/(T1-T0))+
     '                (0.5d0*(YQ(ny_p1,9)+YQ(ny_p2,9)))
     '                *((T1-TIME)/(T1-T0)) !t_wall (external pressure)
                  ENDIF

                  BETA_A=CQ(6,nq)*LAMBDA_1H+CQ(12,nq)
!calc beta as half time step
                  BETA_NEXT_A=CQ(6,nq_next)*LAMBDA_2H+CQ(12,nq_next)
                  BETA_V=CQ(9,nq)*LAMBDA_1H+CQ(10,nq)
!calc beta as half time step
                  
                  BETA_NEXT_V=CQ(9,nq_next)*LAMBDA_2H+CQ(10,nq_next)                  
                  BETA_A=DMAX1(BETA_A,1.0d0)
                  BETA_NEXT_A=DMAX1(BETA_NEXT_A,1.0d0)
                  BETA_V=DMAX1(BETA_V,1.0d0)
                  BETA_NEXT_V=DMAX1(BETA_NEXT_V,1.0d0)

                  ROI=CQ(4,nq)/DSQRT(LAMBDA_1H) ! calculated at half time step
                  ROI_NEXT=CQ(4,nq_next)/DSQRT(LAMBDA_2H)
! for the half time step
                  
                  IF(I.EQ.0) THEN !artery
                    test1=2.0D0*YQ(ny_r1,yp_value)
                    test2=ROI+ROI_NEXT
                    IF(test1.GT.test2)THEN
                      TERM1=0.5d0*(CQ(5,nq)+CQ(5,nq_next))
                      TERM2=(((YQ(ny_r1,yp_value)*2.0d0)
     '                  /(ROI+ROI_NEXT))
     '                  **((BETA_A+BETA_NEXT_A)/2.0D0))-1.0d0
                    ELSE  !negative pressure
                      fo1=BETA_A*CQ(5,nq)/CQ(7,nq)
                      fo2=BETA_NEXT_A*CQ(5,nq_next)/CQ(7,nq_next)
                      TERM1=0.5D0*(FO1+FO2)
                      TERM2=1.0D0-(((ROI+ROI_NEXT)/
     '                  (YQ(ny_r1,yp_value)*2.0d0))**
     '                  (0.5D0*(CQ(7,NQ)+CQ(7,NQ_NEXT))))
                    ENDIF
C term4 is the trace term added on to veesel pressure to include
C transmural pressure as a function of time
                    YQ(ny_p1,yp_value)=(TERM1*TERM2)-TERM4

                  ELSE !vien thus Ro=CQ(4,nq)*sqrt(1.5)                    
                    IF ((2.0D0*YQ(ny_r1,yp_value)).GT.
     '                ((ROI+ROI_NEXT)*VIEN_RATIO)) THEN
C Positive pressure
                      TERM1=(CQ(11,nq)+CQ(11,nq_next))/2.0d0
                      TERM2=(((YQ(ny_r1,yp_value)*2.0d0)
     '                  /((ROI+ROI_NEXT)*
     '                  VIEN_RATIO))**
     '                  ((BETA_V+BETA_NEXT_V)/2.0d0))-1.0d0                      
                    ELSE !negitive pressure
                      fo1=BETA_V*CQ(11,nq)/CQ(7,nq)
                      fo2=BETA_NEXT_V*CQ(11,nq_next)/CQ(7,nq_next)
                      TERM1=0.5D0*(FO1+FO2)
                      TERM2=1.0D0-((((ROI+ROI_NEXT)*
     '                  VIEN_RATIO)/
     '                  (YQ(ny_r1,yp_value)*2.0d0))**
     '                  (0.5D0*(CQ(7,NQ)+CQ(7,NQ_NEXT))))
                    ENDIF
                  ENDIF
                  YQ(ny_p1,yp_value)=(TERM1*TERM2)-TERM4
                ENDDO !calculates viens and arteries (for half timestep)
                
              ELSE !Lax-Wendroff full time step

C PM 13-DEC-01 : See the note added on the same day for half time step.
                IF(TERMINAL) THEN
                  nq_pre=NXQ(-1,1,nq,1)
                  nq_next=NXQ(1,1,nq,1)
                ELSE
                  nq_pre=CONECT(-1,1,nq)
                  nq_next=CONECT(1,1,nq)
                ENDIF

                DELTA_X1=0.0d0
                DELTA_X2=0.0d0

                LAMBDA_1OLD=(YQ(NYNQ(3,nq_pre,0),2)
     '            *(TIME-T0-TINCR)/(T1-T0))+
     '            (YQ(NYNQ(3,nq_pre,0),9)*(T1-TIME+TINCR)/(T1-T0))
                LAMBDA_1=(YQ(NYNQ(3,nq_pre,0),2)*(TIME-T0)/(T1-T0))+
     '            (YQ(NYNQ(3,nq_pre,0),9)*(T1-TIME)/(T1-T0))

                LAMBDA_2OLD=(YQ(NYNQ(3,nq,0),2)
     '            *(TIME-T0-TINCR)/(T1-T0))+
     '            (YQ(NYNQ(3,nq,0),9)*(T1-TIME+TINCR)/(T1-T0))
                LAMBDA_2=(YQ(NYNQ(3,nq,0),2)*(TIME-T0)/(T1-T0))+
     '            (YQ(NYNQ(3,nq,0),9)*(T1-TIME)/(T1-T0))

                LAMBDA_3OLD=(YQ(NYNQ(3,nq_next,0),2)
     '            *(TIME-T0-TINCR)/(T1-T0))+
     '            (YQ(NYNQ(3,nq_next,0),9)*(T1-TIME+TINCR)/(T1-T0))
                LAMBDA_3=(YQ(NYNQ(3,nq_next,0),2)*(TIME-T0)/(T1-T0))+
     '            (YQ(NYNQ(3,nq_next,0),9)*(T1-TIME)/(T1-T0))

                DO nj=1,NJTOT!calculate delta x for half time step
                  TMP=XQ(nj,nq_pre)-XQ(nj,nq)
                  DELTA_X1=DELTA_X1+(TMP*TMP)

                  TMP=XQ(nj,nq)-XQ(nj,nq_next)
                  DELTA_X2=DELTA_X2+(TMP*TMP)
                ENDDO
                
C calculates spatial step
                DELTA_X1=DSQRT(DELTA_X1)
                DELTA_X2=DSQRT(DELTA_X2)
                DELTA_X=((DELTA_X1+DELTA_X2)/2.0d0)*
     '            ((LAMBDA_1OLD+LAMBDA_1+(LAMBDA_2OLD*2.0d0)
     '            +(LAMBDA_2*2.0d0)+LAMBDA_3OLD+LAMBDA_3)*0.125D0)

C PM 29-NOV-01 : Find whether computaion on venous network is required.
                IF(((VENOUS_NETWORK.EQ.'Y').OR.
     '            (VENOUS_NETWORK.EQ.'y')).AND.(N_VENOUS_GEOM.EQ.1))
     '            THEN
                  NET_WORK_NO=1
                ELSE
                  NET_WORK_NO=0
                ENDIF
                
                DO I=0,NET_WORK_NO ! calculates viens and arteries
                  ROI=CQ(4,nq)*
     '              DSQRT((LAMBDA_2OLD+LAMBDA_2+LAMBDA_3OLD+LAMBDA_3)
     &              /4.0d0)
                  
                  ! calculated at half time step
                  ROI_PRE=CQ(4,nq_pre)*
     '              DSQRT((LAMBDA_1OLD+LAMBDA_1+LAMBDA_2OLD+LAMBDA_2)
     '              /4.0d0)
                  
! for the full time step
                  ny_p1=NYNQ(1+(I*VIEN_OFFSET),nq_pre,0)
                  ny_r1=NYNQ(2+(I*VIEN_OFFSET),nq_pre,0)
                  ny_v1=NYNQ(3+(I*VIEN_OFFSET),nq_pre,0)
                  ny_p2=NYNQ(1+(I*VIEN_OFFSET),nq,0)
                  ny_r2=NYNQ(2+(I*VIEN_OFFSET),nq,0)
                  ny_v2=NYNQ(3+(I*VIEN_OFFSET),nq,0)
                  V1=YQ(ny_v1,yp_index)
                  V2=YQ(ny_v2,yp_index)
                  R1=YQ(ny_r1,yp_index)
                  R2=YQ(ny_r2,yp_index)
                  P1=YQ(ny_p1,yp_index)
                  P2=YQ(ny_p2,yp_index)
C finds ny's and then YP values associated with those ny's for the two                
C 1/2 step around grid point nq

C radius at full time step                  
                  TERM1=YQ(ny_r2,1)
                  TERM2=(R2+R1)*(V2-V1)/4.0d0
                  TERM3=(V2+V1)*(R2-R1)/2.0d0
                  YQ(ny_r2,yp_value)=
     '              ((DSQRT(LAMBDA_2OLD)
     '              *TERM1)-((TINCR/DELTA_X)*(TERM2+TERM3)))/
     '              DSQRT(LAMBDA_2)

C velocity at full time step
                  TERM1=YQ(ny_v2,1)
                  TERM2=(((2.0d0*ALPHA)-1.0d0)/2.0d0)*(V2+V1)*(V2-V1)
                  TERM3=(ALPHA-1.0d0)*
     '              (((V2+V1)*(V2+V1))/(R2+R1))*(R2-R1)
                  
C KSB adding gravity term to TERM4 (velocity)
                  DX=0.d0
                  DO nj=1,NJT
                    TMP=XQ(nj,nq)-XQ(nj,nq_pre)
                    DX=DX+TMP*TMP
                    DELTA(nj)=TMP
                  ENDDO
                  DX=DSQRT(DX)
                  G_TERM=0.d0
                  DO nj=1,NJT
                    IF(DX.NE.0.d0) G_TERM=G_TERM+DX*DENSITY
     &                *G_VECTOR(nj)*GRAVITY*1000.d0*(DELTA(nj)/DX)
                  ENDDO
                  P_DIFF=P2-P1+G_TERM
                  TERM4=(1.0d0/DENSITY)*P_DIFF
                  IF((R1+R2).GE.(ROI+ROI_PRE)) THEN
                    TERM5=4.0d0*TINCR*VISCOSITY*((ALPHA/(ALPHA-1))*
     '                ((V2+V1)/((R2+R1)*(R2+R1))))
                  ELSE
                    IF(I.EQ.0) THEN !artery
                      RO=(ROI+ROI_PRE)/2.0D0
                    ELSE !vien thus Ro=CQ(4,nq)*sqrt(1.5)
                      RO=((ROI+ROI_PRE)/2.0D0)*VIEN_RATIO
                    ENDIF
                    ELIP_A=DSQRT((RO*RO)+
     '                DSQRT(RO**4.0d0-((R2+R1)/2.0d0)**4.0d0))
                    ELIP_B=DSQRT((RO*RO)-
     '                DSQRT(RO**4.0d0-((R2+R1)/2.0d0)**4.0d0))
                    ELIP_TERM=DSQRT((2.0d0*ELIP_A**3.0d0*ELIP_B**3.0d0)/
     '                (((R2+R1)/2.0d0)*((R2+R1)/2.0d0)*(ELIP_A*ELIP_A+
     '                ELIP_B*ELIP_B)))
                    TERM5=4.0d0*TINCR*VISCOSITY*((ALPHA
     '              /(ALPHA-1))*
     '                ((V2+V1)/((2.0d0*ELIP_TERM)*(2.0d0*ELIP_TERM))))
                  ENDIF
                  YQ(ny_v2,yp_value)=TERM1-((TINCR/DELTA_X)*
     '              (TERM2+TERM3+TERM4))-TERM5

C pressure at full time step
                  BETA_A=CQ(6,nq)*LAMBDA_2+CQ(12,nq)
!calc beta at full time step arteries
                  BETA_V=CQ(9,nq)*LAMBDA_2+CQ(10,nq)
                  BETA_A=DMAX1(BETA_A,1.0d0)
                  BETA_V=DMAX1(BETA_V,1.0d0)
                  
!calc beta at full time step viens
                  ROI=CQ(4,nq)/DSQRT(LAMBDA_2)
! calculated at FULL time step
                  IF(ITYP12(nr,nx).EQ.2) THEN !pulmonary flow - pleural pressure term
                    DO nj=1,NJT
                      HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq)
                    ENDDO
                    G_PLEURAL=0.d0 !gravitational force
                    PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
                    DO nj=1,NJT
                      G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*G_VECTOR(nj)
     &                  *GRAVITY*1000.d0*HEIGHT(nj)
                    ENDDO
                    TERM4=-(PLEURAL_P+G_PLEURAL) !pleural pressure term
                  ELSEIF(ITYP12(nr,nx).EQ.3) THEN !external pressure field read in
                    nj=13 !Now uses material parameter array CQ nj=13 see example 9823b
                    TERM4=CQ(nj,nq)/1000.d0 !Pa->kPa
                  ELSE !coronary or other flow
                    TERM4=(YQ(ny_p2,2)*((TIME-T0)/(T1-T0)))
     '                +(YQ(ny_p2,9)*((T1-TIME)/(T1-T0)))
                  ENDIF

                  IF(I.EQ.0) THEN !artery
                    IF (YQ(ny_r2,yp_value).GT.ROI) THEN !Positive pressure
                      TERM1=CQ(5,nq)
                      TERM2=(((YQ(ny_r2,yp_value)*2.0d0)
     '                  /((2.d0*ROI
     '                  )))**BETA_A)-1.0d0
                    ELSE !negitive pressure
                      TERM1=BETA_A*CQ(5,nq)/CQ(7,nq)
                      TERM2=1.0D0-((ROI/YQ(ny_r2,yp_value))
     '                  **CQ(7,nq))
                    ENDIF
                  ELSE !vien
                    IF (YQ(ny_r2,yp_value).GT.
     '                (ROI*VIEN_RATIO)) THEN !Positive pressure                      
                      TERM1=CQ(11,nq)
                      TERM2=(((YQ(ny_r2,yp_value)*2.0d0)
     '                  /((2.d0*ROI)*VIEN_RATIO))
     '                  **BETA_V)-1.0d0
                    ELSE !negitive pressure
                      TERM1=BETA_V*CQ(11,nq)/CQ(7,nq)
                      TERM2=1.0D0-(((ROI*VIEN_RATIO)
     '                  /YQ(ny_r2,yp_value))**CQ(7,nq))                      
                    ENDIF
                  ENDIF                  
C term4 is the trace term added on to veesel pressure to include
C transmural pressure. It is interploated between the two mecahnics
C steps
                  YQ(ny_p2,yp_value)=(TERM1*TERM2)-TERM4
                ENDDO ! calculates viens and arteries

              ENDIF ! full/half timestep
            ENDIF
          ENDIF !flow in elastic tube
        ENDIF !Lax-Wendroff
      ENDIF! (UPDATE_VECTOR)

      CALL EXITS('XPFD30')
      RETURN
 9999 CALL ERRORS('XPFD30',ERROR)
      CALL EXITS('XPFD30')
      RETURN 1
      END


