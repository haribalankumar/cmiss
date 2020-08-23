      SUBROUTINE GRADPHIQRC(NENQ,niqV,nq,NQS,NQXI,NXQ,DPHIDX,DXDXIQ,YQ,
     &  ERROR,*)

C#### Subroutine: GRADPHIQRC
C###  Description:
C###    GRADPHIQRC calculates the gradient of Phi at grid point
C###    nq in rectangular cartesian coordinates.
C***    APPROX=1 gets 2 point 1 sided and 2 point 2 sided approximations
C***    APPROX=2 gets 3 point 1 sided and 2 point 2 sided approximations
C***    APPROX=3 gets 4 point 1 sided and 4 point 2 sided approximations
C***  Created by Martin Buist 18-Dec-1998

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NENQ(0:8,NQM),niqV,nq,NQS(NEQM),NQXI(0:NIM,NQSCM),
     &  NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 DPHIDX(3),DXDXIQ(3,3,NQM),YQ(NYQM,NIQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER APPROX,ne,ni,nj,nq1,nq2,nq3,nq4,SCHEME
      REAL*8 DET,DPHIDXI(3),DXDXI(3,3),DXI,DXIDX(3,3)

      CALL ENTERS('GRADPHIQRC',*9999)

      ne=NENQ(1,nq)
      SCHEME=NQS(ne)
      APPROX=2
      DXI=0.5d0

      IF(USE_LAT.EQ.0) THEN
        !calculate dphi/dxi
        DO ni=1,NQXI(0,SCHEME)
          IF((NXQ(-ni,1,nq).EQ.0).OR.(NXQ(-ni,0,nq).GT.1)) THEN
            !no grid point in -ni
            !use onde sided difference
            IF(APPROX.EQ.1) THEN
              nq1=NXQ(ni,1,nq)
              DPHIDXI(ni)=-(YQ(nq,niqV)-YQ(nq1,niqV))/DXI
            ELSEIF(APPROX.EQ.2) THEN
              nq1=NXQ(ni,1,nq)
              nq2=NXQ(ni,1,nq1)
              IF(nq2.EQ.0) THEN
                DPHIDXI(ni)=-(YQ(nq,niqV)-YQ(nq1,niqV))/DXI
              ELSE
                DPHIDXI(ni)=-((3.0d0*YQ(nq,niqV))-(4.0d0*YQ(nq1,niqV))+
     &            (1.0d0*YQ(nq2,niqV)))/(2.0d0*DXI)
              ENDIF
            ELSEIF(APPROX.EQ.3) THEN
              nq1=NXQ(ni,1,nq)
              nq2=NXQ(ni,1,nq1)
              nq3=NXQ(ni,1,nq2)
              DPHIDXI(ni)=-((11.0d0*YQ(nq,niqV))-(18.0d0*YQ(nq1,niqV))+
     &          (9.0d0*YQ(nq2,niqV))-(2.0d0*YQ(nq3,niqV)))/(6.0d0*DXI)
            ENDIF
          ELSEIF((NXQ(ni,1,nq).EQ.0).OR.(NXQ(ni,0,nq).GT.1)) THEN
            !no grid point in +ni
            !use one sided difference
            IF(APPROX.EQ.1) THEN
              nq1=NXQ(-ni,1,nq)
              DPHIDXI(ni)=(YQ(nq,niqV)-YQ(nq1,niqV))/DXI
            ELSEIF(APPROX.EQ.2) THEN
              nq1=NXQ(-ni,1,nq)
              nq2=NXQ(-ni,1,nq1)
              IF(nq2.EQ.0) THEN
                DPHIDXI(ni)=(YQ(nq,niqV)-YQ(nq1,niqV))/DXI
              ELSE
                DPHIDXI(ni)=((3.0d0*YQ(nq,niqV))-(4.0d0*YQ(nq1,niqV))+
     &            (1.0d0*YQ(nq2,niqV)))/(2.0d0*DXI)
              ENDIF
            ELSEIF(APPROX.EQ.3) THEN
              nq1=NXQ(-ni,1,nq)
              nq2=NXQ(-ni,1,nq1)
              nq3=NXQ(-ni,1,nq2)
              DPHIDXI(ni)=((11.0d0*YQ(nq,niqV))-(18.0d0*YQ(nq1,niqV))+
     &          (9.0d0*YQ(nq2,niqV))-(2.0d0*YQ(nq3,niqV)))/(6.0d0*DXI)
            ENDIF
          ELSE
            !points in both +ni and -ni
            !we can use a two sided difference
            IF((APPROX.EQ.1).OR.(APPROX.EQ.2)) THEN
              nq1=NXQ(-ni,1,nq)
              nq2=NXQ(ni,1,nq)
              DPHIDXI(ni)=(YQ(nq2,niqV)-YQ(nq1,niqV))/(2.0d0*DXI)
            ELSEIF(APPROX.EQ.3) THEN
              nq1=NXQ(-ni,1,nq)
              nq2=NXQ(ni,1,nq)
              nq3=NXQ(-ni,1,nq1)
              nq4=NXQ(ni,1,nq2)
              DPHIDXI(ni)=(YQ(nq3,niqV)-(8.0d0*YQ(nq1,niqV))+(8.0d0*
     &          YQ(nq2,niqV))-YQ(nq4,niqV))/(12.0d0*DXI)
            ENDIF
          ENDIF
        ENDDO

        !calculate dxi/dx
        IF(NQXI(0,SCHEME).EQ.1) THEN
          IF(DABS(DXDXIQ(1,1,nq)).GT.ZERO_TOL) THEN
            DXIDX(1,1)=1.0d0/DXDXIQ(1,1,nq)
            DXIDX(1,2)=0.0d0
            DXIDX(1,3)=0.0d0
          ELSE
            ERROR='>>Divide by zero'
            GOTO 9999
          ENDIF
        ELSE
          DO ni=1,NQXI(0,SCHEME)
            DO nj=1,NJT
              DXDXI(nj,ni)=DXDXIQ(nj,ni,nq)
            ENDDO
          ENDDO
          CALL ASSERT(NJT.EQ.NQXI(0,SCHEME),
     &      ' >>NJT.NE.NIT - cannot invert',ERROR,*9999)
          CALL INVERT(NJT,DXDXI,DXIDX,DET)
        ENDIF

        !calculate dphi/dx
        DO nj=1,NJT
          DPHIDX(nj)=0.0d0
          DO ni=1,NQXI(0,SCHEME)
            DPHIDX(nj)=DPHIDX(nj)+DPHIDXI(ni)*DXIDX(ni,nj)
          ENDDO
        ENDDO
      ELSE !use MP inverse to find dphi/dx
        ERROR='Need MP to be passed to GRADPHIQRC'
        GOTO 9999
      ENDIF

      CALL EXITS('GRADPHIQRC')
      RETURN
 9999 CALL ERRORS('GRADPHIQRC',ERROR)
      CALL EXITS('GRADPHIQRC')
      RETURN 1
      END
