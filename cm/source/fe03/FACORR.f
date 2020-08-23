      SUBROUTINE FACORR(XD,ERROR,*)

C#### Subroutine: FACORR
C###  Description:
C###    FACORR Corrects fibre angle from co-ordinate rig
C###    by adjusting for sloping surface.
C**** Created by CS 16 May 2001
C****    This is taken from Ian LeGrice's geom program
C****    with some modification to work in cmiss

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b00.cmn'

!     Parameter List
      REAL*8 XD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
C      IMPLICIT REAL (A-H,O-Z)
      INTEGER  ID,L,L2,LG,LISTS,M,MB,MG,N,ND,NDD,NLEVEL,NX,LINES(75)
      REAL*8 ALPHA,ALPHAR,AXIS,AXISNG,adj_an_dist,BANGLE,BANGLR,
     '  BETA,BETAR,BOT,D(75,100,5),DUM,DUM1,DUM2,DTG,
     '  DTGR,DTHETA,DTHETR,GAMMAR,LFA,OFA,
     '  OFAR,OMEGA,OMEGAR,PRC,PRNG,RC,RNB,RNG,
     '  RNG1,RNG2,RR,RRG,THETA,TOB,TOP
C      CHARACTER ZLABEL(NDMAX)*12,DLABEL*9

      CALL ENTERS('FACORR',*9999)

      ND=0
      NLEVEL=0

C     ** Fill array for new level **

 1    L=1
      M=1
      ID=0
C 5    ND=ND+1
      ND=ND+1

C**** CS This isn't necessary as we will assume we only
C****    loaded points to correct
C8    IF(((INDEX(ZLABEL(ND),'EPIC').GT.0)
C     ,   .AND.(INDEX(ZLABEL(ND),'BASE').EQ.0))
C     ,   .OR.(INDEX(ZLABEL(ND),'MY').GT.0)
C     ,   .OR.((INDEX(ZLABEL(ND),'EN').GT.0)
C     ,   .AND.(INDEX(ZLABEL(ND),'VA').EQ.0))) THEN
C        DLABEL=ZLABEL(ND)(1:9)
C        NLEVEL=NLEVEL+1
        NDD=ND-1
C        WRITE(3,'(2X,I5,3X,A9)')ND,DLABEL
        GO TO 10
C      ELSE
C        IF(ND.GE.NDT) GO TO 901
C        GO TO 5
C      ENDIF

 9    ND=ND+1

 10   M=M+1
      DO 20 N=1,4
        D(L,M,N)=XD(N,ND)
 20   CONTINUE

C**** CS We will only work with one level at a time since we don't have labels
CC *** Test for end of level
C      IF(ZLABEL(ND+1)(1:9).NE.DLABEL) THEN
C        ID=1
C        LISTS=L
C        ND=NDD
C        GO TO 25
C      ENDIF
      IF(ND.EQ.NDT) THEN
        ID=1
        LISTS=L
        ND=NDD
        GO TO 25
      ENDIF

C *** Test for end of annulus

      IF(DABS(XD(3,ND)-XD(3,ND+1)).LT.1.0d0) THEN
        GO TO 9
      ENDIF

C *** Dummy first and last points

 25   DO 30  N=1,4
        D(L,1,N)=D(L,M,N)
        D(L,M+1,N)=D(L,2,N)
 30   CONTINUE
      LINES(L)=M
      IF(ID.EQ.1) GO TO 100
      L=L+1
      M=1
      GO TO 9

C *** CORRECTION
C *** ----------

 100  CONTINUE
      DO 900 L=1,LISTS
        DO 800 M=2,LINES(L)
          ND=ND+1
          RC=D(L,M,1)
          THETA=D(L,M,2)
          AXIS=D(L,M,3)
          OFA=D(L,M,4)
          LFA=OFA

C ***     BETA  CORRECTION  (B)

C**** CS We assume all data coming though here has a fibre angle
C          IF((LFA.EQ.0).OR.(INDEX(ZLABEL(ND),'NFA').GT.0))GO TO 800
          IF(LFA.EQ.0)GO TO 800
          IF(DABS(OFA).LE.90.0d0) THEN
            MB=M-1
          ELSE
            MB=M+1
          ENDIF
          RNB=D(L,MB,1)
          RR =RC/RNB
C         WRITE(*,*) 'ND= ',ND,(D(L,M,NX),NX=1,4)
          DTHETA=DABS(D(L,M,2)-D(L,MB,2))

C ***     ? Past 360 ?  (DTHETA > 200  say)

          IF(DTHETA.GT.200.0d0) DTHETA=360.0d0-DTHETA

C ***     Step  too big ? -No  correction

          IF(DTHETA .GT. 60.0d0) THEN
            BETA=0.0d0
            GO TO 550
          ENDIF

C ***     Calculate  BETA

          DTHETR=DTHETA*PI/180.0d0
          TOP=SIN(DTHETR)
          BOT=SQRT(1.0d0+RR*(RR-2.0d0*COS(DTHETR)))
          IF(BOT.EQ.0.0d0) THEN
            WRITE(*,*) ' Divide by zero. Check data point for'
            WRITE(*,*) ' duplication or single point in annulus.'
            WRITE(*,*) ' at point.     ',ND,(D(L,M,NX),NX=1,4)
            GOTO 800
          ENDIF
          TOB=TOP/BOT
          IF(DABS(TOB).GT.1.0d0) THEN
            TOB=1.0d0*DABS(TOB)/TOB
          ENDIF
C  150     BANGLR=ASIN(TOB)
          BANGLR=ASIN(TOB)
          BANGLE=BANGLR*180.0d0/PI
          IF(DABS(OFA).LE.90.0d0) THEN
            IF(RC/COS(DTHETR).GE.RNB) THEN
              BETA=90.0d0-BANGLE
            ELSE
              BETA=-(90.0d0-BANGLE)
            ENDIF
          ELSE
            IF(RC/COS(DTHETR).GE.RNB) THEN
              BETA=-(90.0d0-BANGLE)
            ELSE
              BETA=90.0d0-BANGLE
            ENDIF
          ENDIF

C ***      GAMMA  CORRECTION  (G)

C ***     Find Radius for GAMMA  Calculation
C ***     1) Special Cases - First and Last Annulae


C**** CS This seems to assume the data is collected in annulae
C**** from base to apex. That is not always the case, so
C**** I have generalised this.
C 550      IF(L.EQ.1) THEN
C            LG=2
C            GO TO 580
C          ENDIF
C          IF(L.EQ.LISTS) THEN
C            LG=LISTS-1
C            GO TO 580
C          ENDIF
CC ***     2) Middle  Annulae
C
C          OFADUM=(DABS(OFA))/OFA
C          IF(OFADUM.LT.0.0d0) THEN
C            LG=L-1
C            ELSE
C            LG=L+1
C          ENDIF

 550      adj_an_dist=RMAX
          DO L2=1,LISTS
            IF(D(L2,1,3).NE.D(L,M,3)) THEN
              IF(DABS(D(L2,1,3)-D(L,1,3)).LT.adj_an_dist) THEN
                adj_an_dist=DABS(D(L2,1,3)-D(L,1,3))
                LG=L2
              ENDIF
            ENDIF
          ENDDO


C ***     Search for THETA

C 580      MG=1
          MG=1
 600      MG=MG+1
          DUM1=D(LG,MG,2)
          DUM2=D(LG,MG-1,2)
          IF(DUM1.LT.DUM2) THEN
            IF(DUM1.GE.THETA.AND.THETA.GE.0.0d0) GO TO 650
            IF(DUM2.LT.THETA.AND.THETA.LE.360.0d0) GO TO 650
          ENDIF
          IF(DUM1.GE.THETA.AND.DUM2.LT.THETA) GO TO 650
          GO TO 600
 650      IF(D(LG,MG,2).EQ.THETA) THEN
            RNG=D(LG,MG,1)
            GO TO 700
          ELSE

C ***       Calculate  New  Radius

            DTG=(D(LG,MG,2)-D(LG,MG-1,2))
            IF(DTG/DABS(DTG).LT.0.0d0) DTG=360.0d0+DTG

C ***       Step  too  big?  No  Correction

            IF(DTG.GT.60.0d0) THEN
              GAMMAR=0.0d0
              BETAR=BETA*PI/180.0d0
              GO TO 720
            ENDIF
            RNG1=D(LG,MG,1)
            RNG2=D(LG,MG-1,1)
            RRG=RNG1/RNG2
            DTGR=DTG*PI/180.0d0
            TOP=SIN(DTGR)
            BOT=SQRT(1.0d0+RRG*(RRG-2.0d0*COS(DTGR)))
            IF(BOT.EQ.0.0d0) THEN
              WRITE(*,*) ' Divide by zero. Gamma calculation.'
              WRITE(*,*) ' Check data on adjacent annulus.'
              WRITE(*,*) ' point = ',ND
              GOTO 800
            ENDIF
            TOB=TOP/BOT
            IF(DABS(TOB).GT.1.0d0) THEN
              TOB=1.0d0*DABS(TOB)/TOB
            ENDIF
C 660        ALPHAR=ASIN(TOB)
            ALPHAR=ASIN(TOB)
            ALPHA=ALPHAR*180.0d0/PI
            IF(RNG2.GE.RNG1/COS(DTGR)) ALPHA=180.0d0-ALPHA
            OMEGA=180.0d0-(ALPHA-THETA+D(LG,MG,2))
            OMEGAR=OMEGA*PI/180.0d0
            RNG=(RNG1*SIN(ALPHAR))/SIN(OMEGAR)
          ENDIF

C ***     Calculate  GAMMA

 700      BETAR=BETA*PI/180.0d0
          PRC=RC*COS(BETAR)
          PRNG=RNG*COS(BETAR)
          AXISNG=D(LG,MG,3)
          TOP=PRNG-PRC
          BOT=AXISNG-AXIS
          GAMMAR=ATAN(TOP/BOT)

C ***     FIBRE  ANGLE  CORRECTION

 720      OFAR=OFA*PI/180.0d0
          TOP=(COS(BETAR))*(SIN(OFAR))
          BOT=COS(GAMMAR)*COS(OFAR)+SIN(GAMMAR)*SIN(BETAR)*SIN(OFAR)
          IF(BOT.EQ.0.0d0) THEN
            D(L,M,5)=90.0d0
          ELSE
            D(L,M,5)=(ATAN (TOP/BOT))*180.0d0/PI
          ENDIF
          IF(D(L,M,5).EQ.0.0d0) THEN
            D(L,M,5)=1.0d-1
            WRITE (*,*) 'NFA=0.0 at point',ND
          ENDIF
          DUM=D(L,M,4)/D(L,M,5)
          IF(DUM.LT.0.0d0) THEN
            IF(D(L,M,5).LT.0.0d0) THEN
              D(L,M,5)=D(L,M,5)+180.0d0
            ELSE
              D(L,M,5)=-(180.0d0-D(L,M,5))
            ENDIF
          ENDIF
          XD(4,ND)=D(L,M,5)
 800    CONTINUE
 900  CONTINUE
      IF(ND.LT.NDT) GO TO 1
C 901  CONTINUE

      WRITE(*,*) ' Fibre angles corrected. '
      WRITE(*,*) ' ND = ',ND , '    Should = NDT = ',NDT
      WRITE(*,*) ' Number of levels used   = ',NLEVEL

      CALL EXITS('FACORR')
      RETURN
 9999 CALL ERRORS('FACORR',ERROR)
      CALL EXITS('FACORR')
      RETURN 1
      END


