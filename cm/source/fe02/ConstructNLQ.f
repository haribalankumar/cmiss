      SUBROUTINE ConstructNLQ(NAQ,NLQ,NWQ,NXQ,ERROR,*)

C#### Subroutine: ConstructNLQ
C###  Description:
C###    Construct NLQ interpolating connections.

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NAQ(NQM,NAM),NLQ(NQM),
     '  NWQ(8,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IJ,IK,mq,na,nq,nii,nij,nik,NITB

      CALL ENTERS('ConstructNLQ',*9999)

      NITB=2 !temporary

      DO na=1,NMGT !loops over all grid levels
        IF(DOP) THEN
          WRITE(OP_STRING,'(//'' Grid level na='',I2)') na
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !DOP
        DO nq=1,NQT
          IF(NLQ(nq).EQ.na.AND.NWQ(1,nq,na).EQ.0) THEN !interior residual g.p.
            IK=MAX(0,NITB-2) !zero for 1,2D, one for 3D
            IJ=MIN(NITB-1,1) !zero for 1D, one for 2,3D

C         Do adjacent grid pts to left & right and above & below
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' grid pts adjacent to nq='',I8)') nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !DOP
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,na),na),na)!adj g.p.
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' nik='',I2,'' nij='',I2,'
     '                //''' nii='',I2,'' mq='',I8)') nik,nij,nii,mq
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF !DOP
                  IF(mq.EQ.0) THEN
                    WRITE(OP_STRING,'('' Error: nq='',I6,'
     '                //''' nii='',I2,'' nij='',I2,'' nik='',I2,'
     '                //''' outside grid in ConstructNLQ'')')
     '                nq,nii,nij,nik
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF !mq=0
                  IF(NLQ(mq).EQ.0) THEN !not residual or interpolated g.p.
                    IF(NAQ(mq,na+1).EQ.0) THEN  !belongs to level na+1
                      NLQ(mq)=na+1 !..so set as residual in level na+1
                    ELSE !interpolated from level na+1 using NAQ
                      NLQ(mq)=11*na+1      !.. so set connectivity
                    ENDIF
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' mq='',I8,'' new NLQ='',I3)')
     '                  mq,NLQ(mq)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF !DOP
                  ENDIF !NLQ
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik

          ENDIF !interior residual g.p.
        ENDDO !nq
      ENDDO !na

      CALL EXITS('ConstructNLQ')
      RETURN
 9999 CALL ERRORS('ConstructNLQ',ERROR)
      CALL EXITS('ConstructNLQ')
      RETURN 1
      END


