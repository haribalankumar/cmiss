      SUBROUTINE NQDS(APPROX,na,NENQ,ni,niq,nq,NQS,NQXI,NXQ,AQ,CQ,
     &  DNUDXQ,DS,DXDXIQ,DXDXIQ2,XQ,YQ,FLUX,ERROR,*)

C#### Subroutine: NQDS
C###  Description:
C###    NQDS finds the arc length derivative about a grid point
C###    for a given direction. If FLUX is true then the arc
C###    length derivative of flux is found. If FLUX is false
C###    the arc length derivative of potential is found.
C###    If APPROX is 1 then a two point difference is calculated.
C###    If APPROX is 2 then a four point difference is calculated.
C**** Created by Martin Buist 11-Dec-1998

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'maqloc00.inc'

!     Parameter List
      INTEGER APPROX,na,NENQ(0:8,NQM),ni,niq,nq,NQS(NEQM),
     &  NQXI(0:NIM,NQSCM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM),DNUDXQ(3,3,NQM),DS,DXDXIQ(3,3,NQM),
     &  DXDXIQ2(3,3,NQM),XQ(NJM,NQM),YQ(NYQM,NIQM,NAM)
      LOGICAL FLUX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,nq0,nq1,nq2,nq3,nq4
      REAL*8 DQ1,DQ2,DQ3,DQ4,DYQ1,LOCFLUX(5),LOCPOTE(5)

      CALL ENTERS('NQDS',*9999)

      IF(APPROX.EQ.1) THEN
        nq1=NXQ(-ni,1,nq,na)
        nq2=NXQ(ni,1,nq,na)

        IF(nq1.EQ.0) THEN
          DQ1=0.0d0
          DO nj=1,NJT
            DQ1=DQ1+(XQ(nj,nq)-XQ(nj,nq2))**2.0d0
          ENDDO
          DQ1=DSQRT(DQ1)

          IF(FLUX) THEN
            CALL GGRADPHIQDN(NENQ,niq,nq,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(1),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq2,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(3),YQ(1,1,na),ERROR,*9999)

            DYQ1=LOCFLUX(3)-LOCFLUX(1)
          ELSE
            DYQ1=YQ(nq2,niq,na)-YQ(nq,niq,na)
          ENDIF

          CALL ASSERT(DQ1.GT.0,' >>Zero arc length',ERROR,*9999)
          DS=DYQ1/DQ1
        ELSE IF(nq2.EQ.0) THEN
          DQ1=0.0d0
          DO nj=1,NJT
            DQ1=DQ1+(XQ(nj,nq)-XQ(nj,nq1))**2.0d0
          ENDDO
          DQ1=DSQRT(DQ1)

          IF(FLUX) THEN
            CALL GGRADPHIQDN(NENQ,niq,nq1,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(1),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(3),YQ(1,1,na),ERROR,*9999)

            DYQ1=LOCFLUX(3)-LOCFLUX(1)
          ELSE
            DYQ1=YQ(nq,niq,na)-YQ(nq1,niq,na)
          ENDIF

          CALL ASSERT(DQ1.GT.0,' >>Zero arc length',ERROR,*9999)
          DS=DYQ1/DQ1
        ELSE
          DQ1=0.0d0
          DQ2=0.0d0
          DO nj=1,NJT
            DQ1=DQ1+(XQ(nj,nq)-XQ(nj,nq1))**2.0d0
            DQ2=DQ2+(XQ(nj,nq)-XQ(nj,nq2))**2.0d0
          ENDDO
          DQ1=DSQRT(DQ1)
          DQ2=DSQRT(DQ2)

          IF(FLUX) THEN
            CALL GGRADPHIQDN(NENQ,niq,nq1,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(1),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq2,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(3),YQ(1,1,na),ERROR,*9999)

            DYQ1=LOCFLUX(3)-LOCFLUX(1)
          ELSE
            DYQ1=YQ(nq2,niq,na)-YQ(nq1,niq,na)
          ENDIF

          CALL ASSERT((DQ1+DQ2).GT.0,' >>Zero arc length',ERROR,*9999)

          DS=DYQ1/(DQ1+DQ2)
        ENDIF

      ELSEIF(APPROX.EQ.2) THEN
        nq2=nq
        nq1=NXQ(-ni,1,nq2,na)
        nq0=NXQ(-ni,1,nq1,na)
        nq3=NXQ(ni,1,nq2,na)
        nq4=NXQ(ni,1,nq3,na)
C        nq1=NXQ(-ni,1,nq2,na)
C        IF(NXQ(-ni,0,nq2,na).GT.1) nq1=NXQ(-ni,2,nq2,na)
C        nq0=NXQ(-ni,1,nq1,na)
C        IF(NXQ(-ni,0,nq1,na).GT.1) nq0=NXQ(-ni,2,nq1,na)
C        nq3=NXQ(ni,1,nq2,na)
C        IF(NXQ(ni,0,nq2,na).GT.1) nq3=NXQ(ni,2,nq2,na)
C        nq4=NXQ(ni,1,nq3,na)
C        IF(NXQ(ni,0,nq3,na).GT.1) nq4=NXQ(ni,2,nq3,na)

        IF(nq1.EQ.0.OR.nq0.EQ.0) THEN
          DQ1=0.0d0
          DQ2=0.0d0
          DO nj=1,NJT
            DQ1=DQ1+(XQ(nj,nq2)-XQ(nj,nq3))**2.0d0
            DQ2=DQ2+(XQ(nj,nq3)-XQ(nj,nq4))**2.0d0
          ENDDO
          DQ1=DSQRT(DQ1)
          DQ2=DSQRT(DQ2)

          IF(FLUX) THEN
            CALL GGRADPHIQDN(NENQ,niq,nq2,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(3),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq3,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(4),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq4,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(5),YQ(1,1,na),ERROR,*9999)

            DYQ1=-((3.0d0*LOCFLUX(3))-(4.0d0*LOCFLUX(4))+
     &        (1.0d0*LOCFLUX(5)))
          ELSE
            LOCPOTE(3)=YQ(nq2,niq,na)
            LOCPOTE(4)=YQ(nq3,niq,na)
            LOCPOTE(5)=YQ(nq4,niq,na)

            DYQ1=-((3.0d0*LOCPOTE(3))-(4.0d0*LOCPOTE(4))+
     &        (1.0d0*LOCPOTE(5)))
          ENDIF

          DS=DYQ1/(DQ1+DQ2)
        ELSE IF(nq3.EQ.0.OR.nq4.EQ.0) THEN
          DQ1=0.0d0
          DQ2=0.0d0
          DO nj=1,NJT
            DQ1=DQ1+(XQ(nj,nq0)-XQ(nj,nq1))**2.0d0
            DQ2=DQ2+(XQ(nj,nq1)-XQ(nj,nq2))**2.0d0
          ENDDO
          DQ1=DSQRT(DQ1)
          DQ2=DSQRT(DQ2)

          IF(FLUX) THEN
            CALL GGRADPHIQDN(NENQ,niq,nq0,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(1),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq1,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(2),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq2,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(3),YQ(1,1,na),ERROR,*9999)

            DYQ1=((3.0d0*LOCFLUX(3))-(4.0d0*LOCFLUX(2))+
     &        (1.0d0*LOCFLUX(1)))
          ELSE
            LOCPOTE(1)=YQ(nq0,niq,na)
            LOCPOTE(2)=YQ(nq1,niq,na)
            LOCPOTE(3)=YQ(nq2,niq,na)

            DYQ1=((3.0d0*LOCPOTE(3))-(4.0d0*LOCPOTE(2))+
     &        (1.0d0*LOCPOTE(1)))
          ENDIF

          DS=DYQ1/(DQ1+DQ2)
        ELSE
          DQ1=0.0d0
          DQ2=0.0d0
          DQ3=0.0d0
          DQ4=0.0d0
          DO nj=1,NJT
            DQ1=DQ1+(XQ(nj,nq0)-XQ(nj,nq1))**2.0d0
            DQ2=DQ2+(XQ(nj,nq1)-XQ(nj,nq2))**2.0d0
            DQ3=DQ3+(XQ(nj,nq2)-XQ(nj,nq3))**2.0d0
            DQ4=DQ4+(XQ(nj,nq3)-XQ(nj,nq4))**2.0d0
          ENDDO
          DQ1=DSQRT(DQ1)
          DQ2=DSQRT(DQ2)
          DQ3=DSQRT(DQ3)
          DQ4=DSQRT(DQ4)

          IF(FLUX) THEN
            CALL GGRADPHIQDN(NENQ,niq,nq0,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(1),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq1,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(2),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq3,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(4),YQ(1,1,na),ERROR,*9999)
            CALL GGRADPHIQDN(NENQ,niq,nq4,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &        DXDXIQ,DXDXIQ2,LOCFLUX(5),YQ(1,1,na),ERROR,*9999)

            DYQ1=LOCFLUX(1)-(8.0d0*LOCFLUX(2))+(8.0d0*LOCFLUX(4))-
     &        LOCFLUX(5)
C            DYQ1=-DYQ1
          ELSE
            LOCPOTE(1)=YQ(nq0,niq,na)
            LOCPOTE(2)=YQ(nq1,niq,na)
            LOCPOTE(3)=YQ(nq2,niq,na)
            LOCPOTE(4)=YQ(nq3,niq,na)
            LOCPOTE(5)=YQ(nq4,niq,na)

            DYQ1=LOCPOTE(1)-(8.0d0*LOCPOTE(2))+(8.0d0*LOCPOTE(4))-
     &        LOCPOTE(5)
C            DYQ1=-DYQ1
          ENDIF

          DS=DYQ1/(3.0d0*(DQ1+DQ2+DQ3+DQ4)) !12h
        ENDIF
      ELSE
        ERROR=' >>Invalid approximation type'
        GOTO 9999
      ENDIF

      CALL EXITS('NQDS')
      RETURN
 9999 CALL ERRORS('NQDS',ERROR)
      CALL EXITS('NQDS')
      RETURN 1
      END
