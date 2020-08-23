      SUBROUTINE CALC_LATT_COEF_TRANS(NITB,nq,NWQ,SUPPORT,AQ,M,NQGW,
     &  PROPQ,FIXQ,ERROR,*)

C#### Subroutine: CALC_LATT_COEF_TRANS
C###  Description:
C###    CALC_LATT_COEF_TRANS calculates the coefficients for matrix
C###    rows when using the lattice based grid points. These
C###    coefficients are based upon taking appropriate weights
C###    (found using a Moore-Penrose Pseudo Inverse) for the
C###    grid points supporting each nq.

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'solv00.cmn'

!     Parameter List
      INTEGER NITB,nq,NWQ(8),SUPPORT(0:NQGM)
      REAL*8 AQ(NMAQM),M(9,NQGM),NQGW(NQGM),PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM)

!     Local Variables
      INTEGER CROSS_LOOK(3,3),i,j,k
      REAL*8 GDOTN(3),DSIGMADX(3)

      DATA CROSS_LOOK /4,7,8,
     &                 7,5,9,
     &                 8,9,6/

      CALL ENTERS('CALC_LATT_COEF_TRANS',*9999)

      IF(NWQ(1).EQ.0) THEN !internal
        DO j=1,NITB
          DSIGMADX(j)=0.0d0
          DO i=1,SUPPORT(0)
            DO k=1,NITB
              DSIGMADX(j)=DSIGMADX(j)+M(k,i)*PROPQ(k,j,1,1,SUPPORT(i))
            ENDDO
          ENDDO
        ENDDO
        DO i=1,SUPPORT(0)
          NQGW(i)=0.0d0
          DO j=1,NITB
            NQGW(i)=NQGW(i)+DSIGMADX(j)*M(j,i)
            DO k=1,NITB
              NQGW(i)=NQGW(i)+PROPQ(k,j,1,1,nq)*M(CROSS_LOOK(k,j),i)
            ENDDO !k
          ENDDO !j
        ENDDO
      ELSE !external
        IF(FIXQ(nq,1)) THEN !potential boundary condition
          DO i=1,SUPPORT(0)
            IF (SUPPORT(i).EQ.nq) THEN
              NQGW(i)=1.0d0
            ELSE
              NQGW(i)=0.0d0
            ENDIF
          ENDDO
        ELSEIF(FIXQ(nq,2)) THEN !flux boundary condition
          CALL CALC_LATT_GDOTN(NITB,AQ,GDOTN,PROPQ(1,1,1,1,nq),
     &      ERROR,*9999)
          DO i=1,SUPPORT(0)
            NQGW(i)=0.0d0
            DO j=1,NITB
              NQGW(i)=NQGW(i)+M(j,i)*GDOTN(j)
            ENDDO
          ENDDO !i
        ENDIF !FIXQ
      ENDIF !NWQ

      CALL EXITS('CALC_LATT_COEF_TRANS')
      RETURN
 9999 CALL ERRORS('CALC_LATT_COEF_TRANS',ERROR)
      CALL EXITS('CALC_LATT_COEF_TRANS')
      RETURN 1
      END


