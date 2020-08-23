      SUBROUTINE CALC_STIMULUS2(maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,
     '  maqp2i,niqV,NQGP,NQGP_PIVOT,nx_ext,AQ,NQGW,
     '  PSTIMULUS,STIMULUS,T,YQ,CALCULATE_STIMULUS,ERROR,*)

C#### Subroutine: CALC_STIMULUS2
C###  Description:
C###    CALC_STIMULUS2

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,niqV,
     '  NQGP(0:NQGM),NQGP_PIVOT(NQGM),nx_ext
      REAL*8 AQ(NMAQM),NQGW(NQGM),PSTIMULUS,STIMULUS,T,
     '  YQ(NYQM,NIQM,NAM,NXM)
      LOGICAL CALCULATE_STIMULUS
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i

      CALL ENTERS('CALC_STIMULUS2',*9999)

      !Calculate Is from current injection conditions
      STIMULUS=0.0d0
      IF(CALCULATE_STIMULUS) THEN
        IF(DABS(AQ(maqp1i)).GT.ZERO_TOL.AND.
     '    T.GE.AQ(maqp1t0).AND.T.LT.AQ(maqp1t1)) THEN
          STIMULUS=AQ(maqp1i)
        ELSE IF(DABS(AQ(maqp2i)).GT.ZERO_TOL.AND.
     '    T.GE.AQ(maqp2t0).AND.T.LT.AQ(maqp2t1)) THEN
          STIMULUS=AQ(maqp2i)
        ENDIF
      ENDIF

      !Calculate the diffusive parts of the pseudo stimulus
      PSTIMULUS=0.0d0
      IF(NQGP(0).EQ.3) THEN !1d internal
        PSTIMULUS=PSTIMULUS+
     '    NQGW(NQGP_PIVOT(1))*YQ(NQGP(1),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(2))*YQ(NQGP(2),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(3))*YQ(NQGP(3),niqV,1,nx_ext)
      ELSE IF(NQGP(0).EQ.9) THEN !2d internal
        PSTIMULUS=PSTIMULUS+
     '    NQGW(NQGP_PIVOT(1))*YQ(NQGP(1),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(2))*YQ(NQGP(2),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(3))*YQ(NQGP(3),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(4))*YQ(NQGP(4),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(5))*YQ(NQGP(5),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(6))*YQ(NQGP(6),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(7))*YQ(NQGP(7),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(8))*YQ(NQGP(8),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(9))*YQ(NQGP(9),niqV,1,nx_ext)
      ELSE IF(NQGP(0).EQ.19) THEN !3d internal
        PSTIMULUS=PSTIMULUS+
     '    NQGW(NQGP_PIVOT( 1))*YQ(NQGP( 1),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 2))*YQ(NQGP( 2),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 3))*YQ(NQGP( 3),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 4))*YQ(NQGP( 4),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 5))*YQ(NQGP( 5),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 6))*YQ(NQGP( 6),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 7))*YQ(NQGP( 7),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 8))*YQ(NQGP( 8),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT( 9))*YQ(NQGP( 9),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(10))*YQ(NQGP(10),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(11))*YQ(NQGP(11),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(12))*YQ(NQGP(12),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(13))*YQ(NQGP(13),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(14))*YQ(NQGP(14),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(15))*YQ(NQGP(15),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(16))*YQ(NQGP(16),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(17))*YQ(NQGP(17),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(18))*YQ(NQGP(18),niqV,1,nx_ext)+
     '    NQGW(NQGP_PIVOT(19))*YQ(NQGP(19),niqV,1,nx_ext)
      ELSE                       !All other grid schemes
        DO i=1,NQGP(0)
          PSTIMULUS=PSTIMULUS+(NQGW(NQGP_PIVOT(i))*
     '      YQ(NQGP(i),niqV,1,nx_ext))
        ENDDO !i
      ENDIF

      CALL EXITS('CALC_STIMULUS2')
      RETURN
 9999 CALL ERRORS('CALC_STIMULUS2',ERROR)
      CALL EXITS('CALC_STIMULUS2')
      RETURN 1
      END



