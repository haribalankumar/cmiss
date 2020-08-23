      SUBROUTINE GGRADPHIQDN(NENQ,niqV,nq,NQS,NQXI,NXQ,AQ,CQ,DNUDXQ,
     &  DXDXIQ,DXDXIQ2,GDPHIDN,YQ,ERROR,*)

C#### Subroutine: GGRADPHIQDN
C###  Description:
C###    GGRADPHIQDN calculates g*(grad phi).n where g is a conductivity
C###    tensor and n is a normal. This represents the current flow
C###    into or out of the grid domain.
C***  Created by Martin Buist 21-Dec-1998

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'maqloc00.inc'

!     Parameter List
      INTEGER NENQ(0:8,NQM),niqV,nq,NQS(NEQM),NQXI(0:NIM,NQSCM),
     &  NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     &  DXDXIQ2(3,3,NQM),GDPHIDN,YQ(NYQM,NIQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER maq,nj1,nj2,nj3
      REAL*8 DET,DNUDX(3,3),DPHIDX(3),DXDNU(3,3),G(3,3),GGRADPHI(3),
     &  XNLOCAL(3)

      CALL ENTERS('GGRADPHIQDN',*9999)

      !calculate dphi/dx
      CALL GRADPHIQRC(NENQ,niqV,nq,NQS,NQXI,NXQ,DPHIDX,DXDXIQ,YQ,
     &  ERROR,*9999)

      !calculate g(nj) conductivities
      DO nj2=1,NJT
        DO nj1=1,NJT
          DNUDX(nj1,nj2)=DNUDXQ(nj1,nj2,nq)
        ENDDO
      ENDDO
      CALL INVERT(NJT,DNUDX,DXDNU,DET)
      DO nj1=1,NJT
        DO nj2=1,NJT
          G(nj1,nj2)=0.0d0
          DO nj3=1,NJT
            G(nj1,nj2)=G(nj1,nj2)+
     &        (DXDNU(nj1,nj3)*DNUDX(nj3,nj2)*CQ(nj3))
          ENDDO
        ENDDO
      ENDDO
      
      !calculate a normal vector at nq
      IF(USE_LAT.EQ.0) THEN
        CALL NORM31(NJT,nq,NXQ,DXDXIQ,DXDXIQ2,XNLOCAL,ERROR,*9999)
      ELSE
C        CALL NORM_LATTICE(NJT,nq,NWQ,DXDXIQ,DXDXIQ2,XNLOCAL,ERROR,
C     &    *9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_X,ERROR,*9999)
        XNLOCAL(1)=AQ(maq,nq)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Y,ERROR,*9999)
        XNLOCAL(2)=AQ(maq,nq)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Z,ERROR,*9999)
        XNLOCAL(3)=AQ(maq,nq)
      ENDIF
      !calculate g*(grad phi)
      DO nj1=1,NJT
        GGRADPHI(nj1)=0.0d0
        DO nj2=1,NJT
          GGRADPHI(nj1)=GGRADPHI(nj1)+(G(nj1,nj2)*DPHIDX(nj2))
        ENDDO
      ENDDO

      !calculate g*(grad phi).n
      GDPHIDN=0.0d0
      DO nj1=1,NJT
        GDPHIDN=GDPHIDN+(GGRADPHI(nj1)*XNLOCAL(nj1))
      ENDDO

      CALL EXITS('GGRADPHIQDN')
      RETURN
 9999 CALL ERRORS('GGRADPHIQDN',ERROR)
      CALL EXITS('GGRADPHIQDN')
      RETURN 1
      END
