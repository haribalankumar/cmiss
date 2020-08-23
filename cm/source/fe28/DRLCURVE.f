      SUBROUTINE DRLCURVE(ISEG,ISIZE_TBH,ISPLOTXY,PHI,REG_PARAMETER,
     '  SIGMA_T_BH,U_PHI,U_T_BH,WK1_INV,CSEG,STRING,ERROR,*)

C#### Subroutine: DRLCUR
C###  Description:
C###    DRLCUR draws the L-curve.
CC JMB 13-MAR-2000

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'plot02.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER ISEG(*),ISIZE_TBH(*),ISPLOTXY(2)
      REAL*8 PHI(NY_TRANSFER_M,NY_TRANSFER_M),
     '  REG_PARAMETER(0:NTSM),SIGMA_T_BH(NY_TRANSFER_M),
     '  U_PHI(NY_TRANSFER_M,NY_TRANSFER_M),
     '  U_T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M)
      CHARACTER CSEG(*)*(*),STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,INDEX,ISIZE_SM_TBH(2),iw,j,N3CO,nts
      REAL*8 AX,BETA2,ETA(NPOINTS),ETA_C,D_PTS(3,NPOINTS),
     '  F,LAMBDA,ONE,RATIO,
     '  RHO(NPOINTS),RHO_C,X(3),Y(3),ZERO
      CHARACTER STR*80
      LOGICAL CBBREV
      EXTERNAL POLYLINE_DYNAM,POLYMARKER_DYNAM
      PARAMETER (iw=10,ONE=1.0d0,ZERO=0.0d0)
!     Functions
      INTEGER IFROMC,INDEX_POLYLINE,INDEX_POLYMARKER
      REAL*8 DDOT,DNRM2

      CALL ENTERS('DRLCURVE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C-----------------------------------------------------------------------

C#### Command: FEM draw lcurve
C###  Parameter:    <index #[1]>
C###    Specify the time/equation index of the L-curve.
C###  Parameter:    <rgb=RGB[blue]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draw the L-curve plot on workstation 10.

        OP_STRING(1)=STRING(1:IEND)//' <index #[1]>'
        OP_STRING(2)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C-----------------------------------------------------------------------

      ELSEIF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRLCURVE',ERROR,*9999)
      ELSE
        CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
     '    ERROR,*9999)
        CALL ASSERT(USE_GRAPHICS.EQ.1,'>>Set USE_GRAPHICS=1',ERROR,
     '    *9999)
        CALL ASSERT(IREGULARISE.EQ.2,'>>Define L-curve criterion',ERROR,
     '    *9999)

        IF(CBBREV(CO,'INDEX',1,noco+1,NTCO,N3CO)) THEN
          nts=IFROMC(CO(N3CO+1))
          CALL ASSERT((nts.GE.1).AND.(nts.LE.
     '      INT(REG_PARAMETER(0))),'>>Non valid index',ERROR,*9999)
        ELSE
          nts=1
        ENDIF

        IF(ISTABILISE.EQ.2) THEN
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1',CO(N3CO+1))
          ELSE
            INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1','BLUE')
          ENDIF
        ELSEIF(ISTABILISE.EQ.3) THEN
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
          ELSE
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
          ENDIF
        ENDIF

        ! Initialisation
        ISIZE_SM_TBH(1)=MIN(ISIZE_TBH(1),ISIZE_TBH(2))
        ISHOLD=.TRUE.
        NEXTPLOT='REPLACE'
        XSCALE='LOG'
        YSCALE='LOG'
        CALL ACWK(iw,0,ERROR,*9999)

        IF(ICOUPLING.EQ.2) THEN
          CALL DGEMV('T',ISIZE_TBH(1),ISIZE_SM_TBH(1),ONE,U_T_BH,
     '      NY_TRANSFER_M,U_PHI(1,nts),1,ZERO,WK1_INV(1,1),1)
          BETA2=DDOT(ISIZE_TBH(1),U_PHI(1,nts),1,U_PHI(1,nts),1)
     '      -DDOT(ISIZE_SM_TBH(1),WK1_INV(1,1),1,WK1_INV(1,1),1)
        ELSE
          CALL DGEMV('T',ISIZE_TBH(1),ISIZE_SM_TBH(1),ONE,U_T_BH,
     '      NY_TRANSFER_M,PHI(1,nts),1,ZERO,WK1_INV(1,1),1)
          BETA2=DDOT(ISIZE_TBH(1),PHI(1,nts),1,PHI(1,nts),1)
     '      -DDOT(ISIZE_SM_TBH(1),WK1_INV(1,1),1,WK1_INV(1,1),1)
        ENDIF

        IF(ISTABILISE.EQ.2) THEN
          ! TGSVD
          ETA(1)=(WK1_INV(1,1)/SIGMA_T_BH(1))**2
          DO i=2,ISIZE_SM_TBH(1)
            ETA(i)=ETA(i-1)+(WK1_INV(i,1)/SIGMA_T_BH(i))**2
          ENDDO
          IF(ISIZE_TBH(1).GT.ISIZE_TBH(2)) THEN
            IF(BETA2.GT.ZERO) THEN
              RHO(ISIZE_SM_TBH(1))=BETA2
            ELSE
              RHO(ISIZE_SM_TBH(1))=ZERO_TOL
            ENDIF
          ELSE
            RHO(ISIZE_SM_TBH(1))=ZERO_TOL
          ENDIF
          DO i=(ISIZE_SM_TBH(1)-1),1,-1
            RHO(i)=RHO(i+1)+WK1_INV(i+1,1)**2
          ENDDO
          DO i=1,ISIZE_SM_TBH(1)
            ETA(i)=DSQRT(ETA(i))
            RHO(i)=DSQRT(RHO(i))
          ENDDO
          RHO_C=RHO(INT(REG_PARAMETER(nts)))
          ETA_C=ETA(INT(REG_PARAMETER(nts)))

          ! Plot the L-curve
          XMAX=RHO(1)
          XMIN=RHO(ISIZE_SM_TBH(1))
          YMAX=ETA(ISIZE_SM_TBH(1))
          YMIN=ETA(1)

          CALL SGPLOTXY(POLYMARKER_DYNAM,INDEX,ISEG,1,ISPLOTXY(1),iw,
     '      ISIZE_SM_TBH(1),D_PTS,RHO,ETA,CSEG,ERROR,*9999)

          IF(ISTABILISE.EQ.2) THEN
            WRITE(STR,'(I5)') INT(REG_PARAMETER(nts))
            CALL STRING_TRIM(STR,IBEG,IEND)
            WRITE(TITLE,*) 'L-curve TSVD corner at ',STR(IBEG:IEND)
          ELSE
            WRITE(STR,'(I5)') INT(REG_PARAMETER(nts))
            CALL STRING_TRIM(STR,IBEG,IEND)
            WRITE(TITLE,*) 'L-curve Greensite corner at ',STR(IBEG:IEND)
          ENDIF

        ELSEIF(ISTABILISE.EQ.3) THEN
          ! Tikhonov
          AX=MAX(SIGMA_T_BH(ISIZE_SM_TBH(1)),DSQRT(LOOSE_TOL))
          RATIO=(SIGMA_T_BH(1)/AX)**(ONE/(DBLE(NPOINTS)-ONE))
          LAMBDA=AX
          DO i=NPOINTS,1,-1
            DO j=1,ISIZE_SM_TBH(1)
              F=SIGMA_T_BH(j)**2/(SIGMA_T_BH(j)**2+LAMBDA**2)
              WK1_INV(j,2)=F*WK1_INV(j,1)/SIGMA_T_BH(j)
              WK1_INV(j,3)=(1-F)*WK1_INV(j,1)
            ENDDO
            ETA(i)=DNRM2(ISIZE_SM_TBH(1),WK1_INV(1,2),1)
            RHO(i)=DNRM2(ISIZE_SM_TBH(1),WK1_INV(1,3),1)
            LAMBDA=RATIO*LAMBDA
          ENDDO

          IF((ISIZE_TBH(1).GT.ISIZE_TBH(2)).AND.BETA2.GT.ZERO) THEN
            DO i=1,NPOINTS
              RHO(i)=DSQRT(RHO(i)**2+BETA2)
            ENDDO
          ENDIF
          DO j=1,ISIZE_SM_TBH(1)
            F=SIGMA_T_BH(j)**2/(SIGMA_T_BH(j)**2+
     '        REG_PARAMETER(nts)**2)
            WK1_INV(j,2)=F*WK1_INV(j,1)/SIGMA_T_BH(j)
            WK1_INV(j,3)=(1-F)*WK1_INV(j,1)
          ENDDO
          ETA_C=DNRM2(ISIZE_SM_TBH(1),WK1_INV(1,2),1)
          RHO_C=DNRM2(ISIZE_SM_TBH(1),WK1_INV(1,3),1)
          IF((ISIZE_TBH(1).GT.ISIZE_TBH(2)).AND.BETA2.GT.ZERO) THEN
            RHO_C=DSQRT(RHO_C**2+BETA2)
          ENDIF

          ! Plot the L-curve
          XMAX=RHO(1)
          XMIN=RHO(NPOINTS)
          YMAX=ETA(NPOINTS)
          YMIN=ETA(1)

          CALL SGPLOTXY(POLYLINE_DYNAM,INDEX,ISEG,1,ISPLOTXY(1),iw,
     '      NPOINTS,D_PTS,RHO,ETA,CSEG,ERROR,*9999)
          WRITE(STR,'(E10.4)') REG_PARAMETER(nts)
          CALL STRING_TRIM(STR,IBEG,IEND)
          WRITE(TITLE,*) 'L-curve Tikhonov corner at ',STR(IBEG:IEND)
        ENDIF

        ! Plot the L-corner
        X(1)=XMIN
        X(2)=RHO_C
        X(3)=RHO_C
        Y(1)=ETA_C
        Y(2)=ETA_C
        Y(3)=YMIN
        ISHOLD=.FALSE.
        NEXTPLOT='ADD'
        INDEX=INDEX_POLYLINE(0,'DASHED','WIDTH1','RED')
        CALL SGPLOTXY(POLYLINE_DYNAM,INDEX,ISEG,1,ISPLOTXY(1),iw,3,
     '    D_PTS,X,Y,CSEG,ERROR,*9999)
        XGRID=.TRUE.
        YGRID=.TRUE.
C LKC 6-NOV-2000 Zero length strings
C        XLABEL=''
C        YLABEL=''
        XLABEL='-'
        YLABEL='-'
        CALL SGPLOTXY(%VAL(0),INDEX,ISEG,2,ISPLOTXY(2),iw,%VAL(0),
     '    D_PTS,%VAL(0),%VAL(0),CSEG,ERROR,*9999)
        CALL DAWK(iw,0,ERROR,*9999)
        CALL REFRESH_GRAPHICS(0.0,ERROR,*9999)
      ENDIF

      CALL EXITS('DRLCURVE')
      RETURN
 9999 CALL ERRORS('DRLCURVE',ERROR)
      CALL EXITS('DRLCURVE')
      RETURN 1
      END

