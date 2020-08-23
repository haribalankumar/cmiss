      SUBROUTINE CALC_LATT_COEF_EXT(NITB,nq,NWQ,SUPPORT,AQ,COEFFSEXT,M,
     &  PROPQ,FIXQ,SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: CALC_LATT_COEF_EXT
C###  Description:
C###    CALC_LATT_COEF_EXT calculates the coefficients for matrix
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
      REAL*8 AQ(NMAQM),COEFFSEXT(NQGM),M(9,NQGM),PROPQ(3,3,4,2,NQM)
      CHARACTER ERROR*(*)
      LOGICAL FIXQ(NYQM,NIYFIXM),SOLVEEIGHTPROBLEM

!     Local Variables
      INTEGER i,j
      REAL*8 GDOTN(3)
      LOGICAL FLUXBC

      CALL ENTERS('CALC_LATT_COEF_EXT',*9999)

      IF(NWQ(1).GT.0) THEN !external grid point
        IF(SOLVEEIGHTPROBLEM) THEN
          IF(SOLVE8_FLUXBC) THEN
            IF(nq.EQ.1) THEN !Need one potential BC
              FLUXBC=.FALSE.
            ELSE
              FLUXBC=.TRUE.
            ENDIF
          ELSEIF(FIXQ(nq,1)) THEN !potential boundary condition
              FLUXBC=.FALSE.
          ELSEIF(FIXQ(nq,2)) THEN !flux boundary condition
              FLUXBC=.TRUE.
          ENDIF
        ELSE
          IF(FIXQ(nq,1)) THEN !potential boundary condition
            FLUXBC=.FALSE.
          ELSEIF(FIXQ(nq,2)) THEN !flux boundary condition
            FLUXBC=.TRUE.
          ELSEIF(FIXQ(nq,3)) THEN
            FLUXBC=.FALSE.
          ELSE
            WRITE(ERROR,'(''>>No boundary condition set at point'
     &        //' (ext.) '',I8)') nq
            GOTO 9999
          ENDIF
        ENDIF
        IF(FLUXBC) THEN
          CALL CALC_LATT_GDOTN(NITB,AQ,GDOTN,PROPQ(1,1,1,2,nq),
     &      ERROR,*9999)
          DO i=1,SUPPORT(0)
            COEFFSEXT(i)=0.0d0
            DO j=1,NITB
              COEFFSEXT(i)=COEFFSEXT(i) + M(j,i)*GDOTN(j)
            ENDDO
          ENDDO
        ELSE
          DO i=1,SUPPORT(0)
            IF (SUPPORT(i).EQ.nq) THEN
              COEFFSEXT(i)=1.0d0
            ELSE
              COEFFSEXT(i)=0.0d0
            ENDIF
          ENDDO
        ENDIF
      ELSE !internal grid point
        IF(NITB.EQ.1) THEN
          DO i=1,SUPPORT(0)
            COEFFSEXT(i)=M(4,i)
          ENDDO
        ELSEIF(NITB.EQ.2) THEN
          DO i=1,SUPPORT(0)
            COEFFSEXT(i)=M(4,i)+M(5,i)
          ENDDO
        ELSE !NITB==3
          DO i=1,SUPPORT(0)
            COEFFSEXT(i)=M(4,i)+M(5,i)+M(6,i)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('CALC_LATT_COEF_EXT')
      RETURN
 9999 CALL ERRORS('CALC_LATT_COEF_EXT',ERROR)
      CALL EXITS('CALC_LATT_COEF_EXT')
      RETURN 1
      END


