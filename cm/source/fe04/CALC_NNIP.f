      SUBROUTINE CALC_NNIP(COLLAPSED_XI,IBT,IDRN,IDR2,IDR3,IMAP,INP,
     '  nb,nb2,NNIP,NNIP2,NUM_COLLAPSED,NUMI1,NUMI2,NUMI3,PERP_XI,
     '  ATCOLLAPSE,DIFFELEMS,SECTOR,ERROR,*)

C#### Subroutine: CALC_NNIP
C###  Description:
C###    CALC_NNIP calculates the local node numbers (nn) of the element
C###    vertices for the element (NNIP) the number of vertices in each
C###    and other information required for refining. If DIFFELEMS is
C###    .TRUE. then the new element (with basis number nb2) created
C###    with refinement is different than the current element (basis
C###    type nb) (as is the case with sectors sometimes) then NNIP2
C###    (for this new element) will also be returned.
C###    SECTOR indicates whether sector calculations will be used.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER COLLAPSED_XI,IBT(3,NIM,NBFM),IDRN,IDR2,IDR3,IMAP(3,4,4,4),
     '  INP(NNM,NIM,NBFM),nb,nb2,NNIP(4,4,4),NNIP2(4,4,4),
     '  NUM_COLLAPSED,NUMI1,NUMI2,NUMI3,PERP_XI
      CHARACTER ERROR*(*)
      LOGICAL ATCOLLAPSE(4,4,4),DIFFELEMS,SECTOR
!     Local Variables
      INTEGER i1,i2,i3,ii1,ii2,ii3,nn
      LOGICAL COLLAPSE

      CALL ENTERS('CALC_NNIP',*9999)

      IF(NIT(nb).EQ.1) THEN
        NUMI2=1
        NUMI3=1
      ELSE IF(NIT(nb).EQ.2) THEN
        IF(IBT(1,IDR2,nb).EQ.5.OR.IBT(1,IDR2,nb).EQ.6) THEN
          IF(IBT(2,IDR2,nb).EQ.4) THEN !Hermite
            NUMI2=2
          ELSE
            NUMI2=IBT(2,IDR2,nb)+1
          ENDIF
        ELSE IF(IBT(1,IDR2,nb).EQ.1) THEN !Lagrange
          NUMI2=IBT(2,IDR2,nb)+1
        ELSE IF(IBT(1,IDR2,nb).EQ.2) THEN !Hermite
          NUMI2=2
        ELSE
          ERROR='>>Invalid basis type'
          GOTO 9999
        ENDIF
        NUMI3=1
      ELSE
        IF(IBT(1,IDR2,nb).EQ.5.OR.IBT(1,IDR2,nb).EQ.6) THEN
          IF(IBT(2,IDR2,nb).EQ.4) THEN !Hermite
            NUMI2=2
          ELSE
            NUMI2=IBT(2,IDR2,nb)+1
          ENDIF
        ELSE IF(IBT(1,IDR2,nb).EQ.1) THEN !Lagrange
          NUMI2=IBT(2,IDR2,nb)+1
        ELSE IF(IBT(1,IDR2,nb).EQ.2) THEN !Hermite
          NUMI2=2
        ELSE
          ERROR='>>Invalid basis type'
          GOTO 9999
        ENDIF
        IF(IBT(1,IDR3,nb).EQ.5.OR.IBT(1,IDR3,nb).EQ.6) THEN
          IF(IBT(2,IDR3,nb).EQ.4) THEN !Hermite
            NUMI3=2
          ELSE
            NUMI3=IBT(2,IDR3,nb)+1
          ENDIF
        ELSE IF(IBT(1,IDR3,nb).EQ.1) THEN !Lagrange
          NUMI3=IBT(2,IDR3,nb)+1
        ELSE IF(IBT(1,IDR3,nb).EQ.2) THEN !Hermite
          NUMI3=2
        ELSE
          ERROR='>>Invalid basis type'
          GOTO 9999
        ENDIF
      ENDIF
      IF(IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6) THEN
        IF(IBT(2,IDRN,nb).EQ.4) THEN !Hermite
          NUMI1=2
        ELSE
          NUMI1=IBT(2,IDRN,nb)+1
        ENDIF
      ELSE IF(IBT(1,IDRN,nb).EQ.1) THEN !Lagrange
        NUMI1=IBT(2,IDRN,nb)+1
      ELSE IF(IBT(1,IDRN,nb).EQ.2) THEN !Hermite
        NUMI1=2
      ELSE
        ERROR='>>Invalid basis type'
        GOTO 9999
      ENDIF
      IF(SECTOR) THEN
C*** Find IMAP(i1,i2,i3,j) which stores the ii values for i_j at the
C*** vertex point i1,i2,i3.
        DO i3=1,NUMI3
          DO i2=1,NUMI2
            DO i1=1,NUMI1
              COLLAPSE=.FALSE.
              IF(NUM_COLLAPSED.EQ.2) THEN
                IF(PERP_XI.EQ.IDR2) THEN
                  IF((IBT(1,IDRN,nb).EQ.5.AND.i2.EQ.1).OR.
     '              (IBT(1,IDRN,nb).EQ.6.AND.i2.EQ.NUMI2)) THEN
                    ii1=1
                    ii3=1
                    COLLAPSE=.TRUE.
                  ELSE
                    ii1=i1
                    ii3=i3
                  ENDIF
                ELSE
                  IF((IBT(1,IDRN,nb).EQ.5.AND.i3.EQ.1).OR.
     '              (IBT(1,IDRN,nb).EQ.6.AND.i3.EQ.NUMI3)) THEN
                    ii1=1
                    ii2=1
                    COLLAPSE=.TRUE.
                  ELSE
                    ii1=i1
                    ii2=i2
                  ENDIF
                ENDIF
              ELSE
                IF(COLLAPSED_XI.EQ.IDRN) THEN
                  IF(PERP_XI.EQ.IDR2) THEN
                    IF((IBT(1,IDRN,nb).EQ.5.AND.i2.EQ.1).OR.
     '                (IBT(1,IDRN,nb).EQ.6.AND.i2.EQ.NUMI2)) THEN
                      ii1=1
                      COLLAPSE=.TRUE.
                    ELSE
                      ii1=i1
                    ENDIF
                  ELSE
                    IF((IBT(1,IDRN,nb).EQ.5.AND.i3.EQ.1).OR.
     '                (IBT(1,IDRN,nb).EQ.6.AND.i3.EQ.NUMI3)) THEN
                      ii1=1
                      COLLAPSE=.TRUE.
                    ELSE
                      ii1=i1
                    ENDIF
                  ENDIF
                  ii2=i2
                  ii3=i3
                ELSE IF(COLLAPSED_XI.EQ.IDR2) THEN
                  IF(PERP_XI.EQ.IDRN) THEN
                    IF((IBT(1,IDR2,nb).EQ.5.AND.i1.EQ.1).OR.
     '                (IBT(1,IDR2,nb).EQ.6.AND.i1.EQ.NUMI1)) THEN
                      ii2=1
                      COLLAPSE=.TRUE.
                    ELSE
                      ii2=i2
                    ENDIF
                  ELSE
                    IF((IBT(1,IDR2,nb).EQ.5.AND.i3.EQ.1).OR.
     '                (IBT(1,IDR2,nb).EQ.6.AND.i3.EQ.NUMI3)) THEN
                      ii2=1
                      COLLAPSE=.TRUE.
                    ELSE
                      ii2=i2
                    ENDIF
                  ENDIF
                  ii1=i1
                  ii3=i3
                ELSE
                  IF(PERP_XI.EQ.IDRN) THEN
                    IF((IBT(1,IDR3,nb).EQ.5.AND.i1.EQ.1).OR.
     '                (IBT(1,IDR3,nb).EQ.6.AND.i1.EQ.NUMI1)) THEN
                      ii3=1
                      COLLAPSE=.TRUE.
                    ELSE
                      ii3=i3
                    ENDIF
                  ELSE
                    IF((IBT(1,IDR3,nb).EQ.5.AND.i2.EQ.1).OR.
     '                (IBT(1,IDR3,nb).EQ.6.AND.i2.EQ.NUMI2)) THEN
                      ii3=1
                      COLLAPSE=.TRUE.
                    ELSE
                      ii3=i3
                    ENDIF
                  ENDIF
                  ii1=i1
                  ii2=i2
                ENDIF
              ENDIF
              ATCOLLAPSE(i1,i2,i3)=COLLAPSE
              IMAP(1,i1,i2,i3)=ii1
              IMAP(2,i1,i2,i3)=ii2
              IMAP(3,i1,i2,i3)=ii3
            ENDDO !i1
          ENDDO !i2
        ENDDO !i3
      ELSE
        DO i3=1,NUMI3
          DO i2=1,NUMI2
            DO i1=1,NUMI1
              IMAP(1,i1,i2,i3)=i1
              IMAP(2,i1,i2,i3)=i2
              IMAP(3,i1,i2,i3)=i3
            ENDDO !i1
          ENDDO !i2
        ENDDO !i3
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' NUMI1='',I1,'', NUMI2='',I1,'
     '    //''', NUMI3='',I1)') NUMI1,NUMI2,NUMI3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO i3=1,NUMI3
        DO i2=1,NUMI2
          DO i1=1,NUMI1
            ii1=IMAP(1,i1,i2,i3)
            ii2=IMAP(2,i1,i2,i3)
            ii3=IMAP(3,i1,i2,i3)
            NNIP(ii1,ii2,ii3)=0
            IF(DIFFELEMS) NNIP2(i1,i2,i3)=0
            IF(NUMI2.EQ.1.AND.NUMI3.EQ.1) THEN !1 xi directions only
              DO nn=1,NNT(nb)
                IF(INP(nn,IDRN,nb).EQ.ii1) NNIP(ii1,ii2,ii3)=nn
              ENDDO !nn
              IF(DIFFELEMS) THEN
                DO nn=1,NNT(nb2)
                  IF(INP(nn,IDRN,nb2).EQ.i1) NNIP2(i1,i2,i3)=nn
                ENDDO !nn
              ENDIF
            ELSE IF(NUMI3.EQ.1) THEN !2 xi directions only
              DO nn=1,NNT(nb)
                IF(INP(nn,IDRN,nb).EQ.ii1.AND.
     '            INP(nn,IDR2,nb).EQ.ii2) NNIP(ii1,ii2,ii3)=nn
              ENDDO !nn
              IF(DIFFELEMS) THEN
                DO nn=1,NNT(nb2)
                  IF(INP(nn,IDRN,nb2).EQ.i1.AND.
     '              INP(nn,IDR2,nb2).EQ.i2) NNIP2(i1,i2,i3)=nn
                ENDDO !nn
              ENDIF
            ELSE !3 xi directions
              DO nn=1,NNT(nb)
                IF(INP(nn,IDRN,nb).EQ.ii1.AND.
     '            INP(nn,IDR2,nb).EQ.ii2.AND.
     '            INP(nn,IDR3,nb).EQ.ii3) NNIP(ii1,ii2,ii3)=nn
              ENDDO !nn
              IF(DIFFELEMS) THEN
                DO nn=1,NNT(nb2)
                  IF(INP(nn,IDRN,nb2).EQ.i1.AND.
     '              INP(nn,IDR2,nb2).EQ.i2.AND.
     '              INP(nn,IDR3,nb2).EQ.i3) NNIP2(i1,i2,i3)=nn
                ENDDO !nn
              ENDIF
            ENDIF
            CALL ASSERT(NNIP(ii1,ii2,ii3).NE.0,
     '        '>>Could not find local node',ERROR,*9999)
            NNIP(i1,i2,i3)=NNIP(ii1,ii2,ii3)
            IF(DIFFELEMS) CALL ASSERT(NNIP2(i1,i2,i3).NE.0,
     '        '>>Could not find local node',ERROR,*9999)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' NNIP1('',I1,'','',I1,'','
     '          //''',I1,'')='',I2)') ii1,ii2,ii3,NNIP(ii1,ii2,ii3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              IF(DIFFELEMS) THEN
                WRITE(OP_STRING,'('' NNIP2('',I1,'','',I1,'','
     '            //''',I1,'')='',I2)') i1,i2,i3,NNIP2(i1,i2,i3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
CC$            call mp_unsetlock()
            ENDIF
          ENDDO !i1
        ENDDO !i2
      ENDDO !i3

      CALL EXITS('CALC_NNIP')
      RETURN
 9999 CALL ERRORS('CALC_NNIP',ERROR)
      CALL EXITS('CALC_NNIP')
      RETURN 1
      END


