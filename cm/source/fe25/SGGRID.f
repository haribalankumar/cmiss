      SUBROUTINE SGGRID(INDEX,ISEG,ISGRID,iw,NLQ,NQLIST,
     '   CSEG,ADAPTIVE,TYPE,XQ,YQ,ERROR,*)

C#### Subroutine: SGGRID
C###  Description:
C###    SGGRID creates grid point or field segment ISGRID.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISGRID,iw,NLQ(NQM),
     '  NQLIST(0:NQM)
      REAL*8 XQ(NJM,NQM),YQ(NYQM,NIQM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TYPE*(*)
      LOGICAL ADAPTIVE
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD,INDEX_POLYMARKER,nj,NLQQ,nq,nqa
      REAL*8 POINT(3),ZVAL
      CHARACTER STR*6,TYPE_MARKER*5,COLOUR*6

      CALL ENTERS('SGGRID',*9999)

      CALL OPEN_SEGMENT(ISGRID,ISEG,iw,'GRID',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)
      DO nqa=1,NQLIST(0)
        nq=NQLIST(nqa)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' nqa='',I5,'' nq='',I5)') nqa,nq
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !dop

        DO nj=1,NJT
          POINT(nj)=XQ(nj,nq)
        ENDDO
        DO nj=NJT+1,3
          POINT(nj)=0.0d0
        ENDDO
        IF(TYPE(1:6).EQ.'POINTS') THEN
C KAT 26Jan00
          IF(ADAPTIVE) THEN
            NLQQ=NLQ(nq)
            IF(NLQQ.LT.10) THEN      !nq is residual point
              TYPE_MARKER='PLUS'
            ELSE                     !nq is interpolated point
              TYPE_MARKER='POINT'
              NLQQ=NLQQ/10
            ENDIF
            IF(NLQQ.EQ.3) THEN      !operation done on grid level 3
              COLOUR='GREEN'
            ELSE IF(NLQQ.EQ.2) THEN !operation done on grid level 2
              COLOUR='BLUE'
            ELSE IF(NLQQ.EQ.1) THEN !operation done on grid level 1
              COLOUR='RED'
            ELSE                         !operation done on grid level >3
              COLOUR='BLACK'
            ENDIF !NLQQ
            INDEX=INDEX_POLYMARKER(0,TYPE_MARKER,'SIZE1',COLOUR)
C           ELSE supplied value for INDEX is used
          ENDIF !adaptive
C            IF(NAQ(nq,3).EQ.0) THEN      !nq belongs to grid level 3
C              COLOUR='GREEN'
C            ELSE IF(NAQ(nq,2).EQ.0) THEN !nq belongs to grid level 2
C              COLOUR='BLUE'
C            ELSE IF(NAQ(nq,1).EQ.0) THEN !nq belongs to grid level 1
C              COLOUR='RED'
C            ELSE                         !nq belongs to grid level >3
C              COLOUR='BLACK'
C            ENDIF !NAQ
C            IF(NLQ(nq).LT.10) THEN       !nq is residual point
C              INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1',COLOUR)
C            ELSE IF(NLQ(nq).GT.10) THEN  !nq is interpolated point
C              INDEX=INDEX_POLYMARKER(0,'POINT','SIZE1',COLOUR)
C            ENDIF !NLQ
C          ENDIF !adaptive
          CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)

        ELSE IF(TYPE(1:7).EQ.'NUMBERS') THEN
          WRITE(STR,'(I6)') nq
          CALL STRING_TRIM(STR,IBEG,IEND)
          CALL TEXT(1,iw,STR(IBEG:IEND),POINT,ERROR,*9999)

        ELSE IF(TYPE(1:6).EQ.'VALUES') THEN
          ZVAL=YQ(nq,1)
          IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
          IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
          IF(COLOUR_WS) THEN
! PJH 9Jan96 INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0d0)
C         Keep index between 45 and 200
            IF(DABS(ZDIFF).GT.1.D-6) THEN
              INDEX=200-NINT((ZVAL-ZMINI)/ZDIFF*155.0d0)
            ELSE
              INDEX=15
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' ZVAL='',E12.3,'' Index='',I3)')
     '          ZVAL,INDEX
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !dop
          ENDIF
          CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)
        ENDIF
      ENDDO !nq
      CALL CLOSE_SEGMENT(ISGRID,iw,ERROR,*9999)

      CALL EXITS('SGGRID')
      RETURN
 9999 CALL ERRORS('SGGRID',ERROR)
      CALL EXITS('SGGRID')
      RETURN 1
      END


