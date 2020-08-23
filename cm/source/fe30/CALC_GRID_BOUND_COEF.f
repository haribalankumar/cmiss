      SUBROUTINE CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,NXQ,COEFF,CQ,
     '  DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,*)

C#### Subroutine: CALC_GRID_BOUND_COEF
C###  Description:
C###    CALC_GRID_BOUND_COEF calculates the flux coefficients for
C###    g*(grad phi).(n)
C***  14/7/98 MLB replacing all NJT with NITB as nit must = njt
C***  Created by Martin Buist, January 1999

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER NENQ(0:8,NQM),nq,NQS(NEQM),NQXI(0:NIM,NQSCM),
     '  NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 COEFF(NQGM),CQ(NMM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     '  DXDXIQ2(3,3,NQM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,n,ne,ni,NITB,nj,nj1,nj2,nj3,nk,SCHEME
      REAL*8 DET,DNUDX(3,3),DXDNU(3,3),DXDXI(3,3),DXIDX(3,3),
     '  G(3,3),TEMPCOEF,XNLOCAL(3)
C MHT 24-08-99 unreferenced
C      INTEGER nq1,nq2,nq3,nq4
C      REAL*8 DXI,
C      LOGICAL FOUND

      CALL ENTERS('CALC_GRID_BOUND_COEF',*9999)

      ne=NENQ(1,nq)
      SCHEME=NQS(ne)
      NITB=NQXI(0,SCHEME)
      DO i=1,NQGM
        COEFF(i)=0.0d0
      ENDDO !i
      DO ni=1,3
        DO nj=1,3
          DXIDX(ni,nj)=0.0d0
        ENDDO
      ENDDO

      !calculate dxi/dx
      DO nj=1,NJT
        DO ni=1,NITB
          DXDXI(nj,ni)=DXDXIQ(nj,ni,nq)
        ENDDO !ni
      ENDDO !nj
      IF(NITB.EQ.1) THEN
        IF(DABS(DXDXI(1,1)).GT.ZERO_TOL) THEN
          DXIDX(1,1)=1.0d0/DXDXI(1,1)
        ELSE
          DXIDX(1,1)=1.0d0
        ENDIF
      ELSE
        CALL ASSERT(NJT.EQ.NITB,
     '    ' >>NJT.NE.NIT - cannot invert',ERROR,*9999)
        CALL INVERT(NITB,DXDXI,DXIDX,DET)
      ENDIF

      !calculate g(nj) conductivities
      DO nj2=1,NITB
        DO nj1=1,NITB
          DNUDX(nj1,nj2)=DNUDXQ(nj1,nj2,nq)
        ENDDO
      ENDDO
      CALL INVERT(NITB,DNUDX,DXDNU,DET)
      DO nj1=1,NITB
        DO nj2=1,NITB
          G(nj1,nj2)=0.0d0
          DO nj3=1,NITB
            G(nj1,nj2)=G(nj1,nj2)+
     '        (DXDNU(nj1,nj3)*DNUDX(nj3,nj2)*CQ(nj3))
          ENDDO
        ENDDO
      ENDDO

      !calculate a normal vector at nq
      CALL NORM31(NITB,nq,NXQ,DXDXIQ,DXDXIQ2,XNLOCAL,ERROR,*9999)

      !calculate coefficients and grid points
      n=1
      DO ni=1,NITB
        TEMPCOEF=0.0d0
        DO nj=1,NITB
          DO nk=1,NITB
            TEMPCOEF=TEMPCOEF+(XNLOCAL(nk)*G(nk,nj)*DXIDX(ni,nj))
          ENDDO
        ENDDO
        IF(NXQ(-ni,1,nq).EQ.0) THEN
          !no grid point in -ni -> use one sided difference
          COEFF(1)=COEFF(1)-(3.0d0*TEMPCOEF)
          n=n+1
          COEFF(n)=COEFF(n)+(4.0d0*TEMPCOEF)
          n=n+1
          COEFF(n)=COEFF(n)-(1.0d0*TEMPCOEF)
        ELSEIF(NXQ(ni,1,nq).EQ.0) THEN
          !no grid point in +ni -> use one sided difference
          COEFF(1)=COEFF(1)+(3.0d0*TEMPCOEF)
          n=n+1
          COEFF(n)=COEFF(n)-(4.0d0*TEMPCOEF)
          n=n+1
          COEFF(n)=COEFF(n)+(1.0d0*TEMPCOEF)
        ELSE
          !points in both +ni and -ni -> use a two sided difference
          n=n+1
          COEFF(n)=COEFF(n)-TEMPCOEF
          n=n+1
          COEFF(n)=COEFF(n)+TEMPCOEF
        ENDIF !NXQ non-zero
      ENDDO !ni

C      DXI=0.0d0
C      DO ni=2,NQGP(0,nq)
C        DXI=DXI+COEFF(ni)
C      ENDDO
C      COEFF(1)=-DXI

      CALL EXITS('CALC_GRID_BOUND_COEF')
      RETURN
 9999 CALL ERRORS('CALC_GRID_BOUND_COEF',ERROR)
      CALL EXITS('CALC_GRID_BOUND_COEF')
      RETURN 1
      END
