      SUBROUTINE CALC_TRANSMEMBRANE_RHS(NENQ,niqV,NQS,NQXI,nr,NWQ,
     '  nx,nx_ext,NXQ,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,RHS,YQ,ERROR,*)

C#### Subroutine: CALC_TRANSMEMBRANE_RHS
C###  Description:
C###    CALC_TRANSMEMBRANE_RHS
C**** Created by Martin Buist, February 2000.

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
C      INCLUDE 'cmiss$reference:cell02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ktyp30.cmn'

!     Parameter list
      INTEGER NENQ(0:8,NQM),niqV,NQS(NEQM),NQXI(0:NIM,NQSCM),nr,
     '  NWQ(8,0:NQM,NAM),nx,nx_ext,NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     '  DXDXIQ2(3,3,NQM),RHS(NQM),
     '  YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nq
      LOGICAL ERROR_FLAG

      CALL ENTERS('CALC_TRANSMEMBRANE_RHS',*9999)

      IF(KTYP38.EQ.1) THEN !no Vm flux
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(niqV,nr,NWQ,nx,RHS,YQ)
        DO nq=NQR(1,nr),NQR(2,nr)
          IF(NWQ(1,nq,1).GT.0) THEN !external
            RHS(nq)=0.0d0
          ELSE !internal
            RHS(nq)=YQ(nq,niqV,1,nx)
          ENDIF !internal/external
        ENDDO !nq
C$OMP   END PARALLEL DO
      ELSE IF(KTYP38.EQ.2) THEN !no PHIi flux
        ERROR_FLAG=.FALSE.
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(nq),
C$OMP&  SHARED(CQ,DNUDXQ,DXDXIQ,DXDXIQ2,RHS,NENQ,niqV,nr,NQS,
C$OMP&    NQXI,nx,NXQ,NWQ,YQ)
        DO nq=NQR(1,nr),NQR(2,nr)
          IF(.NOT.ERROR_FLAG) THEN
            IF(NWQ(1,nq,1).GT.0) THEN !external
              CALL GGRADPHIQDN(NENQ,niqV,nq,NQS,NQXI,NXQ(-NIM,0,0,1),AQ,
     &          CQ(3,nq),DNUDXQ,DXDXIQ,DXDXIQ2,RHS(nq),YQ(1,1,1,nx_ext),
     &          ERROR,*100)
              RHS(nq)=-RHS(nq)
              IF(DABS(RHS(nq)).LT.1.0d-6) RHS(nq)=0.0d0
C              RHS(nq)=-RHS(nq)/1000.0d0 !mS - S

C MLB please leave
C              RHS(nq)=RHS(nq)/1000.0d0 !mS - S
C              IF(DABS(RHS(nq).LT.1.0d-6)) RHS(nq)=0.0d0

            ELSE !internal
              RHS(nq)=YQ(nq,niqV,1,nx)
            ENDIF !internal/external

            GOTO 102
 100        CONTINUE
C$OMP       CRITICAL(CALC_TRANSMEMBRANE_RHS_1)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*101)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '        //'results may be unreliable'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP       END CRITICAL(CALC_TRANSMEMBRANE_RHS_1)
 102        CONTINUE
          ENDIF !.NOT.ERROR_FLAG
        ENDDO !nq
C$OMP   END PARALLEL DO
      ENDIF

      CALL EXITS('CALC_TRANSMEMBRANE_RHS')
      RETURN
 9999 CALL ERRORS('CALC_TRANSMEMBRANE_RHS',ERROR)
      CALL EXITS('CALC_TRANSMEMBRANE_RHS')
      RETURN 1
      END
