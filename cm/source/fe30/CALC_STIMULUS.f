      SUBROUTINE CALC_STIMULUS(maqp1t0,maqp1t1,maqp1i,maqp2t0,
     '  maqp2t1,maqp2i,niqOLD,niqV,nq,NQGP,NQGP_PIVOT,nr,nx_ext,
     '  nx_trans,AQ,CQ,NQGW,PSTIMULUS,STIMULUS,T,THETA,YQ,BIDOMAIN,
     '  CALCULATE_STIMULUS,ERROR,*)

C#### Subroutine: CALC_STIMULUS
C###  Description:
C###    CALC_STIMULUS calculates the applied and pseudo stimulus
C###    currents to be applied to a cell. The pseudo stimulus
C###    current is made up of the explicit component of the
C###    intracellular diffusive current, the contribution
C###    from the bidomain extracellular diffusive current if they
C###    exist and any stretch activated channels. NOTE: The
C###    psuedo stimulus has units of uA/mm^3.

C*** NOTE: This routine replaces CALC_DIFFUSION and GET_EXT_CONTRIB
C***       from previous solution techniques - MLB 8-April-1999.

      IMPLICIT NONE

c     INCLUDE 'cmiss$reference:b01.cmn'
c     INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,niqOLD,
     '  niqV,nq,NQGP(0:NQGM),NQGP_PIVOT(NQGM),nr,nx_ext,nx_trans
      REAL*8 AQ(NMAQM),CQ(NMM),NQGW(NQGM),PSTIMULUS,STIMULUS,T,THETA,
     '  YQ(NYQM,NIQM,NAM,NXM)
      LOGICAL BIDOMAIN,CALCULATE_STIMULUS
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,niqSAC
      LOGICAL SAC
      CHARACTER ERR_DUMMY*32

      CALL ENTERS('CALC_STIMULUS',*9999)

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

C SGM 29 Nov 2000 The explicit component of the intracellular diffusive
C     current for grid-based FE is handled in MARCH8
C MLT 2Dec02 Added grid finite volumes
      IF(ITYP4(nr,nx_trans).NE.6.AND.
     '   ITYP4(nr,nx_trans).NE.7) THEN !Not grid-based FE
        IF(NQGP(0).EQ.3) THEN       !1d internal
          IF(THETA.LT.1.0d0) THEN   !Some explicit component
C SGM 18Dec2000 replace all occurances of niqOLD with niqV
C           as YQ(niqV,nx_trans) hasn't yet been updated to new timestep
C           in MARCH8 and can be used as YQ(niqOLD,nx_trans).
C           YQ(niqOLD,nx_trans) is now being used for temporary storage
C           in MARCH8.
            IF(USE_LAT.EQ.0) THEN
              PSTIMULUS=PSTIMULUS+
     '          (NQGW(NQGP_PIVOT(1))*YQ(NQGP(1),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(2))*YQ(NQGP(2),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(3))*YQ(NQGP(3),niqV,1,nx_trans))*
     '          (1.0d0-THETA)
            ELSE !lattice method
                PSTIMULUS=PSTIMULUS+
     '          (NQGW(1)*YQ(NQGP(1),niqV,1,nx_trans)+
     '          NQGW(2)*YQ(NQGP(2),niqV,1,nx_trans)+
     '          NQGW(3)*YQ(NQGP(3),niqV,1,nx_trans))*
     '          (1.0d0-THETA)                
            ENDIF
          ENDIF
          IF(BIDOMAIN) THEN         !Extracellular contribution
            IF(USE_LAT.EQ.0) THEN
              PSTIMULUS=PSTIMULUS+
     '          NQGW(NQGP_PIVOT(1))*YQ(NQGP(1),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(2))*YQ(NQGP(2),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(3))*YQ(NQGP(3),niqV,1,nx_ext)
            ELSE !lattice method
              PSTIMULUS=PSTIMULUS+
     '          NQGW(1)*YQ(NQGP(1),niqV,1,nx_ext)+
     '          NQGW(2)*YQ(NQGP(2),niqV,1,nx_ext)+
     '          NQGW(3)*YQ(NQGP(3),niqV,1,nx_ext)              
            ENDIF
          ENDIF
        ELSE IF(NQGP(0).EQ.9) THEN  !2d internal
          IF(THETA.LT.1.0d0) THEN   !Some explicit component
            IF(USE_LAT.EQ.0) THEN
              PSTIMULUS=PSTIMULUS+
     '          (NQGW(NQGP_PIVOT(1))*YQ(NQGP(1),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(2))*YQ(NQGP(2),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(3))*YQ(NQGP(3),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(4))*YQ(NQGP(4),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(5))*YQ(NQGP(5),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(6))*YQ(NQGP(6),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(7))*YQ(NQGP(7),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(8))*YQ(NQGP(8),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(9))*YQ(NQGP(9),niqV,1,nx_trans))*
     '          (1.0d0-THETA)
            ELSE
              PSTIMULUS=PSTIMULUS+
     '          (NQGW(1)*YQ(NQGP(1),niqV,1,nx_trans)+
     '          NQGW(2)*YQ(NQGP(2),niqV,1,nx_trans)+
     '          NQGW(3)*YQ(NQGP(3),niqV,1,nx_trans)+
     '          NQGW(4)*YQ(NQGP(4),niqV,1,nx_trans)+
     '          NQGW(5)*YQ(NQGP(5),niqV,1,nx_trans)+
     '          NQGW(6)*YQ(NQGP(6),niqV,1,nx_trans)+
     '          NQGW(7)*YQ(NQGP(7),niqV,1,nx_trans)+
     '          NQGW(8)*YQ(NQGP(8),niqV,1,nx_trans)+
     '          NQGW(9)*YQ(NQGP(9),niqV,1,nx_trans))*
     '          (1.0d0-THETA)              
            ENDIF
          ENDIF
          IF(BIDOMAIN) THEN         !Extracellular contribution
            IF(USE_LAT.EQ.0) THEN
              PSTIMULUS=PSTIMULUS+
     '          NQGW(NQGP_PIVOT(1))*YQ(NQGP(1),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(2))*YQ(NQGP(2),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(3))*YQ(NQGP(3),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(4))*YQ(NQGP(4),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(5))*YQ(NQGP(5),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(6))*YQ(NQGP(6),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(7))*YQ(NQGP(7),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(8))*YQ(NQGP(8),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(9))*YQ(NQGP(9),niqV,1,nx_ext)
            ELSE
              PSTIMULUS=PSTIMULUS+
     '          NQGW(1)*YQ(NQGP(1),niqV,1,nx_ext)+
     '          NQGW(2)*YQ(NQGP(2),niqV,1,nx_ext)+
     '          NQGW(3)*YQ(NQGP(3),niqV,1,nx_ext)+
     '          NQGW(4)*YQ(NQGP(4),niqV,1,nx_ext)+
     '          NQGW(5)*YQ(NQGP(5),niqV,1,nx_ext)+
     '          NQGW(6)*YQ(NQGP(6),niqV,1,nx_ext)+
     '          NQGW(7)*YQ(NQGP(7),niqV,1,nx_ext)+
     '          NQGW(8)*YQ(NQGP(8),niqV,1,nx_ext)+
     '          NQGW(9)*YQ(NQGP(9),niqV,1,nx_ext)              
            ENDIF
          ENDIF
        ELSE IF(NQGP(0).EQ.19) THEN !3d internal
          IF(THETA.LT.1.0d0) THEN   !Some explicit component
            IF(USE_LAT.EQ.0) THEN
              PSTIMULUS=PSTIMULUS+
     '          (NQGW(NQGP_PIVOT( 1))*YQ(NQGP( 1),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 2))*YQ(NQGP( 2),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 3))*YQ(NQGP( 3),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 4))*YQ(NQGP( 4),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 5))*YQ(NQGP( 5),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 6))*YQ(NQGP( 6),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 7))*YQ(NQGP( 7),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 8))*YQ(NQGP( 8),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT( 9))*YQ(NQGP( 9),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(10))*YQ(NQGP(10),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(11))*YQ(NQGP(11),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(12))*YQ(NQGP(12),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(13))*YQ(NQGP(13),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(14))*YQ(NQGP(14),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(15))*YQ(NQGP(15),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(16))*YQ(NQGP(16),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(17))*YQ(NQGP(17),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(18))*YQ(NQGP(18),niqV,1,nx_trans)+
     '          NQGW(NQGP_PIVOT(19))*YQ(NQGP(19),niqV,1,nx_trans))*
     '          (1.0d0-THETA)
            ELSE
              PSTIMULUS=PSTIMULUS+
     '          (NQGW( 1)*YQ(NQGP( 1),niqV,1,nx_trans)+
     '          NQGW( 2)*YQ(NQGP( 2),niqV,1,nx_trans)+
     '          NQGW( 3)*YQ(NQGP( 3),niqV,1,nx_trans)+
     '          NQGW( 4)*YQ(NQGP( 4),niqV,1,nx_trans)+
     '          NQGW( 5)*YQ(NQGP( 5),niqV,1,nx_trans)+
     '          NQGW( 6)*YQ(NQGP( 6),niqV,1,nx_trans)+
     '          NQGW( 7)*YQ(NQGP( 7),niqV,1,nx_trans)+
     '          NQGW( 8)*YQ(NQGP( 8),niqV,1,nx_trans)+
     '          NQGW( 9)*YQ(NQGP( 9),niqV,1,nx_trans)+
     '          NQGW(10)*YQ(NQGP(10),niqV,1,nx_trans)+
     '          NQGW(11)*YQ(NQGP(11),niqV,1,nx_trans)+
     '          NQGW(12)*YQ(NQGP(12),niqV,1,nx_trans)+
     '          NQGW(13)*YQ(NQGP(13),niqV,1,nx_trans)+
     '          NQGW(14)*YQ(NQGP(14),niqV,1,nx_trans)+
     '          NQGW(15)*YQ(NQGP(15),niqV,1,nx_trans)+
     '          NQGW(16)*YQ(NQGP(16),niqV,1,nx_trans)+
     '          NQGW(17)*YQ(NQGP(17),niqV,1,nx_trans)+
     '          NQGW(18)*YQ(NQGP(18),niqV,1,nx_trans)+
     '          NQGW(19)*YQ(NQGP(19),niqV,1,nx_trans))*
     '          (1.0d0-THETA)              
            ENDIF
          ENDIF
          IF(BIDOMAIN) THEN         !Extracellular contribution
            IF(USE_LAT.EQ.0) THEN
              PSTIMULUS=PSTIMULUS+
     '          NQGW(NQGP_PIVOT( 1))*YQ(NQGP( 1),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 2))*YQ(NQGP( 2),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 3))*YQ(NQGP( 3),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 4))*YQ(NQGP( 4),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 5))*YQ(NQGP( 5),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 6))*YQ(NQGP( 6),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 7))*YQ(NQGP( 7),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 8))*YQ(NQGP( 8),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT( 9))*YQ(NQGP( 9),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(10))*YQ(NQGP(10),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(11))*YQ(NQGP(11),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(12))*YQ(NQGP(12),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(13))*YQ(NQGP(13),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(14))*YQ(NQGP(14),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(15))*YQ(NQGP(15),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(16))*YQ(NQGP(16),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(17))*YQ(NQGP(17),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(18))*YQ(NQGP(18),niqV,1,nx_ext)+
     '          NQGW(NQGP_PIVOT(19))*YQ(NQGP(19),niqV,1,nx_ext)
            ELSE
              PSTIMULUS=PSTIMULUS+
     '          NQGW( 1)*YQ(NQGP( 1),niqV,1,nx_ext)+
     '          NQGW( 2)*YQ(NQGP( 2),niqV,1,nx_ext)+
     '          NQGW( 3)*YQ(NQGP( 3),niqV,1,nx_ext)+
     '          NQGW( 4)*YQ(NQGP( 4),niqV,1,nx_ext)+
     '          NQGW( 5)*YQ(NQGP( 5),niqV,1,nx_ext)+
     '          NQGW( 6)*YQ(NQGP( 6),niqV,1,nx_ext)+
     '          NQGW( 7)*YQ(NQGP( 7),niqV,1,nx_ext)+
     '          NQGW( 8)*YQ(NQGP( 8),niqV,1,nx_ext)+
     '          NQGW( 9)*YQ(NQGP( 9),niqV,1,nx_ext)+
     '          NQGW(10)*YQ(NQGP(10),niqV,1,nx_ext)+
     '          NQGW(11)*YQ(NQGP(11),niqV,1,nx_ext)+
     '          NQGW(12)*YQ(NQGP(12),niqV,1,nx_ext)+
     '          NQGW(13)*YQ(NQGP(13),niqV,1,nx_ext)+
     '          NQGW(14)*YQ(NQGP(14),niqV,1,nx_ext)+
     '          NQGW(15)*YQ(NQGP(15),niqV,1,nx_ext)+
     '          NQGW(16)*YQ(NQGP(16),niqV,1,nx_ext)+
     '          NQGW(17)*YQ(NQGP(17),niqV,1,nx_ext)+
     '          NQGW(18)*YQ(NQGP(18),niqV,1,nx_ext)+
     '          NQGW(19)*YQ(NQGP(19),niqV,1,nx_ext)              
            ENDIF
          ENDIF
        ELSE                       !All other grid schemes
          IF(THETA.LT.1.0d0) THEN  !Some explicit component
            IF(USE_LAT.EQ.0) THEN
              DO i=1,NQGP(0)
                PSTIMULUS=PSTIMULUS+(NQGW(NQGP_PIVOT(i))*
     '            YQ(NQGP(i),niqV,1,nx_trans))
              ENDDO !i
              PSTIMULUS=PSTIMULUS*(1.0d0-THETA)
            ELSE
              DO i=1,NQGP(0)
                PSTIMULUS=PSTIMULUS+(NQGW(i)*
     '            YQ(NQGP(i),niqV,1,nx_trans))
              ENDDO !i
              PSTIMULUS=PSTIMULUS*(1.0d0-THETA)              
            ENDIF
          ENDIF
          IF(BIDOMAIN) THEN        !Extracellular Contribution
            IF(USE_LAT.EQ.0) THEN
              DO i=1,NQGP(0)
                PSTIMULUS=PSTIMULUS+(NQGW(NQGP_PIVOT(i))*
     '            YQ(NQGP(i),niqV,1,nx_ext))
              ENDDO !i
            ELSE
              DO i=1,NQGP(0)
                PSTIMULUS=PSTIMULUS+(NQGW(i)*
     '            YQ(NQGP(i),niqV,1,nx_ext))
              ENDDO !i              
            ENDIF
          ENDIF
        ENDIF
      ENDIF !Not grid-based FE

      !Adjust pstimulus for stretch activated channels
      SAC=.FALSE.
      IF(ITYP19(nr,nx_trans).EQ.1.AND.ITYP3(nr,nx_trans).EQ.2) THEN !FHN
        IF(DABS(CQ(18)).GT.ZERO_TOL) SAC=.TRUE.
      ELSEIF(ITYP19(nr,nx_trans).EQ.1.AND.ITYP3(nr,nx_trans).EQ.3)THEN !VCD
        IF(DABS(CQ(13)).GT.ZERO_TOL) SAC=.TRUE.
      ENDIF
      IF(SAC) THEN
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,ERR_DUMMY,*9999)
        PSTIMULUS=PSTIMULUS-YQ(nq,niqSAC,1,nx_trans)
      ENDIF

C *** DPN 16 July 1999 - if solving a single cell (grid point)
      IF (NQT.EQ.1) PSTIMULUS=0.0d0

      CALL EXITS('CALC_STIMULUS')
      RETURN
 9999 CALL ERRORS('CALC_STIMULUS',ERROR)
      CALL EXITS('CALC_STIMULUS')
      RETURN 1
      END



