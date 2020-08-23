      SUBROUTINE CalculateNLQ(NLQ,NWQ,NXQ,AQ,YQ,ERROR,*)

C#### Subroutine: CalculateNLQ
C###  Description:
C###    CalculateNLQ calculates grid point gradients to determine which
C###    grid points should be active at the fine grid level (NLQ(nq)=1).
C###    Grid pts where a stimulus current is applied are also made active.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NLQ(NQM),NWQ(8,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IJ,IK,maqp1i,mq,na,nii,nij,nik,nq,NITB

      CALL ENTERS('CalculateNLQ',*9999)
      na=1   !fine grid level
      NITB=2 !temporary
!Find position maqp1i of stimulus current in AQ array
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,MAQ_CURRENT,
     '  ERROR,*9999)
      write(*,'('' maqp1i='',I3)') maqp1i
      IF(maqp1i.EQ.0) maqp1i=3 !temporary

      DO nq=1,NQT !loop over fine grid points (level 1)

        IF(NWQ(1,nq,na).EQ.0) THEN !interior g.p.
          IK=MAX(0,NITB-2) !zero for 1,2D, one for 3D
          IJ=MIN(NITB-1,1) !zero for 1D, one for 2,3D

          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' grid pts adjacent to nq='',I8)') nq
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !DOP

          IF(DABS(AQ(maqp1i,nq)).GT.ZERO_TOL) THEN !stimulus current exists
            NLQ(nq)=1
            WRITE(OP_STRING,'('' stimulus current applied at nq='',I8)')
     '        nq
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          DO nik=-IK,IK
            DO nij=-IJ,IJ
              DO nii=-1,1
                mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,na),na),na) !adj g.p.
C                write(*,'('' nik='',I2,'' nij='',I2,'
C     '            //''' nii='',I2,'' mq='',I8)') nik,nij,nii,mq
                IF(mq.EQ.0) THEN
                  WRITE(OP_STRING,'('' Error: nq='',I6,'
     '              //''' nii='',I2,'' nij='',I2,'' nik='',I2,'
     '              //''' outside grid in CalculateNLQ'')')
     '              nq,nii,nij,nik
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF !mq=0

                IF(mq.NE.nq
     '            .AND.DABS(YQ(mq,1,1)-YQ(nq,1,1)).GT.1.d-1) THEN !grad exceeded
                  NLQ(nq)=1
                  WRITE(OP_STRING,'('' grad > threshold at mq='',I8)')
     '              mq
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  GOTO 100
                ENDIF !difference
              ENDDO !nii
            ENDDO !nij
          ENDDO !nik
        ENDIF !interior residual g.p.

 100    CONTINUE
      ENDDO !nq

      CALL EXITS('CalculateNLQ')
      RETURN
 9999 CALL ERRORS('CalculateNLQ',ERROR)
      CALL EXITS('CalculateNLQ')
      RETURN 1
      END


