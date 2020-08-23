      SUBROUTINE CALC_CONTRIB_COEFF(COEFF,nb_e,nb_f,nt,NPL,PG,WG,
     '  ERROR,*)

C#### Subroutine: CALC_CONTRIB_COEFF
C###  Description:
C###    CALC_CONTRIB_COEFF calculates the nodal contribution
C###    coefficients for integrated natural boundary conditions.
C**** Created by Carey Stevens 22 Oct 1997

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER nb_e,nb_f,NPL(5,0:3),nt
      REAL*8 COEFF(4),PG(NSM,NUM,NGM,NBM),WG(NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IG(4),N,ng,NGA
      REAL*8 PG_LOCAL,PL1,PL2,PL3,W,WG_LOCAL(10),XI,
     '  XIGG(10)

      DATA XIGG/0.5000000000000D0,0.2113248654051D0,0.7886751345948D0,
     '          0.1127016653792D0,0.5000000000000D0,0.8872983346207D0,
     '          0.0694318442029D0,0.3300094782075D0,0.6699905217924D0,
     '          0.9305681557970D0/
      DATA WG_LOCAL/1.0000000000000D0,0.5000000000000D0,0.50000000000D0,
     '          0.2777777777778D0,0.4444444444444D0,0.2777777777778D0,
     '          0.1739274225687D0,0.3260725774313D0,0.3260725774313D0,
     '          0.1739274225687D0/

      DATA   IG/0,1,3,6/,NGA/4/

      CALL ENTERS('CALC_CONTRIB_COEFF',*9999)


      IF(NIT(nb_e).EQ.2) THEN
        IF(JTYP2B.EQ.1.AND.NPL(4,0).LT.0) GOTO 9998
        IF(NPL(2,1).NE.NPL(3,1)) THEN
C       Calculate equivalent nodal load contribution by
C       Gaussian quadrature
          IF(NPL(1,1).EQ.2) nt=3
          IF(NPL(1,1).EQ.3) nt=4
          DO N=1,nt
            COEFF(N)=0.0d0
            DO ng=1,NGA
              XI=XIGG(IG(NGA)+ng)
              W=WG_LOCAL(IG(NGA)+ng)

              IF(NPL(1,1).EQ.1) THEN
                PG_LOCAL=PL1(N,1,XI)
              ELSE IF(NPL(1,1).EQ.2) THEN
                PG_LOCAL=PL2(N,1,XI)
              ELSE IF(NPL(1,1).EQ.3) THEN
                PG_LOCAL=PL3(N,1,XI)
              ELSE
                CALL ASSERT(.FALSE.,'>>Not implemented for '
     '            //'interpolation type',ERROR,*9999)
              ENDIF
              COEFF(N)=COEFF(N)+PG_LOCAL*W
            ENDDO
          ENDDO
        ENDIF
      ELSE
        DO N=1,nt
            COEFF(N)=0.0d0
            DO ng=1,NGT(nb_f)
              COEFF(N)=COEFF(N)+PG(N,1,ng,nb_f)*WG(ng,nb_f)
            ENDDO
          ENDDO
      ENDIF

 9998 CALL EXITS('CALC_CONTRIB_COEFF')
      RETURN
 9999 CALL ERRORS('CALC_CONTRIB_COEFF',ERROR)
      CALL EXITS('CALC_CONTRIB_COEFF')
      RETURN 1
      END


