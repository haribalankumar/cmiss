      SUBROUTINE ConstructNXQ(NAQ,nr,NWQ,NXQ,ERROR,*)

C#### Subroutine: ConstructNXQ
C###  Description:
C###    ConstructNXQ creates the grid connectivity array NXQ.

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

! Parameter list
      INTEGER NAQ(NQM,NAM),nr,NWQ(8,0:NQM,NAM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      CHARACTER ERROR*(*)
! Local variables
      INTEGER mq1,mq2,mq3,na,ni,nq
      CHARACTER ROW*6

      CALL ENTERS('ConstructNXQ',*9999)
      CALL ASSERT(CALL_GRID,'>>Must have a grid defined',ERROR,*9999)

C MHT 30-06-99 set but not used
C      NITB=NQXI(0,NQS(NEELEM(1,nr)))

! Calculate NXQ(-3:3,1,nq,na+1) array
      DO na=1,NMGT-1
        DO nq=NQR(1,nr),NQR(2,nr)
          IF(NAQ(nq,na+1).GE.0) THEN !nq is in grid na+1
            NXQ(0,1,nq,na+1)=nq
            DO ni=-1,1
              mq3=NXQ(3*ni,1,nq,na) !Xi3 neighbour at grid level na
              IF(mq3.ne.0) NXQ(3*ni,1,nq,na+1)=NXQ(3*ni,1,mq3,na)
              mq2=NXQ(2*ni,1,nq,na) !Xi2 neighbour at grid level na
              IF(mq2.ne.0) NXQ(2*ni,1,nq,na+1)=NXQ(2*ni,1,mq2,na)
              mq1=NXQ(  ni,1,nq,na) !Xi1 neighbour at grid level na
              IF(mq1.ne.0) NXQ(  ni,1,nq,na+1)=NXQ(  ni,1,mq1,na)
            ENDDO !ni
          ENDIF !nq is in grid na+1
        ENDDO !nq
      ENDDO !na

! Calculate NWQ(ni,nq,na=2..) array
      DO na=2,NMGT
        IF(DOP) THEN
          WRITE(OP_STRING,'('' Calculate NWQ for grid na='',I2)') na
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nq=NQR(1,nr),NQR(2,nr)
          DO ni=1,6
            NWQ(ni,nq,na)=0
          ENDDO
          IF(NAQ(nq,na).GE.0) THEN !nq is in grid na
            IF(NXQ(-2,1,nq,na).eq.0) THEN !bottom row
              ROW='Lower'
              mq1=NXQ(2,1,nq,na) !is one row  above
              mq2=NXQ(2,1,mq1,na) !is two rows above
            ELSE IF( NXQ(2,1,nq,na).eq.0) THEN !top row
              ROW='Upper'
              mq1=NXQ(-2,1,nq,na) !is one row  below
              mq2=NXQ(-2,1,mq1,na) !is two rows below
            ELSE
              ROW='Centre'
            ENDIF
            IF(ROW(1:5).EQ.'Lower'.OR.ROW(1:5).EQ.'Upper') THEN
              IF(NXQ(-1,1,nq,na).EQ.0) THEN !LH corner
                NWQ(1,nq,na)=NXQ(1,1,mq1,na)
                NWQ(2,nq,na)=NXQ(1,1,NXQ(1,1,mq2,na),na)
              ELSE IF(NXQ(1,1,nq,na).EQ.0) THEN !RH corner
                NWQ(1,nq,na)=NXQ(-1,1,mq1,na)
                NWQ(2,nq,na)=NXQ(-1,1,NXQ(-1,1,mq2,na),na)
              ELSE !midside
                NWQ(1,nq,na)=mq1
                NWQ(2,nq,na)=mq2
              ENDIF
            ELSE IF(ROW(1:6).EQ.'Centre') THEN
              IF(NXQ(-1,1,nq,na).eq.0) THEN !LH edge
                NWQ(1,nq,na)=NXQ(1,1,nq,na)
                NWQ(2,nq,na)=NXQ(1,1,NXQ(1,1,nq,na),na)
              ELSE IF(NXQ(1,1,nq,na).eq.0) THEN !RH edge
                NWQ(1,nq,na)=NXQ(-1,1,nq,na)
                NWQ(2,nq,na)=NXQ(-1,1,NXQ(-1,1,nq,na),na)
              ENDIF
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'(''***** nq='',I6,'' ROW='',A,'
     '          //''' NWQ='',2I6)') nq,ROW,NWQ(1,nq,na),NWQ(2,nq,na)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF !nq is in grid na
        ENDDO !nq
      ENDDO !na

      CALL EXITS('ConstructNXQ')
      RETURN
 9999 CALL ERRORS('ConstructNXQ',ERROR)
      CALL EXITS('ConstructNXQ')
      RETURN 1
      END


