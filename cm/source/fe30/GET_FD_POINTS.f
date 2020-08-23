      SUBROUTINE GET_FD_POINTS(NEELEM,NENQ,NLQ,NQGP,NQGP_PIVOT,
     '  NQS,NQXI,NRLIST,NWQ,NXQ,RET_ERROR,*)

C#### Subroutine: GET_FD_POINTS
C###  Description:
C###    This routine creates a list of the local quadratic element
C###    grid points for each internal grid point.
C***  Created by Martin Buist 20 August 1997

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter list
      INTEGER NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLQ(NQM),
     '  NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NWQ(8,0:NQM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      CHARACTER RET_ERROR*(*)
!     Local variables
      INTEGER CASE,i,IEND,j,n,ne1,ne2,na,ni,NITB,nq,nq1,nq2,NQGT,
     '  nqHOLD1,nr,nrr
      LOGICAL ADAPTIVE,BRANCH,ERROR_FLAG,FOUND
      CHARACTER ERROR*(ERRSTRLEN)
!     Functions
      INTEGER LEN_TRIM

      CALL ENTERS('GET_FD_POINTS',*9999)

      ERROR_FLAG=.FALSE.

      IF(NMGT.GT.1) THEN !adaptive grid
        ADAPTIVE=.TRUE.
      ELSE               !not adaptive
        ADAPTIVE=.FALSE.
      ENDIF
      nqHOLD1=0

      IF(ADAPTIVE) THEN
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(i,nq)
C$OMP&  SHARED(NQGP)
        DO nq=1,NQM
          DO i=0,19
            NQGP(i,nq)=0
          ENDDO !i
        ENDDO !nq
C$OMP   END PARALLEL DO
      ELSE
        na=1
      ENDIF

      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
        NITB=NQXI(0,NQS(NEELEM(1,nr)))
C$OMP PARALLEL DO
C$OMP&PRIVATE(BRANCH,CASE,FOUND,i,IEND,j,n,na,ne1,ne2,ni,
C$OMP&        nq,nq1,nq2,NQGT)
C$OMP&SHARED(NITB,NENQ,NQGP,NQGP_PIVOT,nqHOLD1,NWQ,NXQ)
        DO nq=NQR(1,nr),NQR(2,nr)

C *** DPN 06 September 1999 - adding alternate return for MP
          IF (.NOT.ERROR_FLAG) THEN

            IF(ADAPTIVE) THEN
              na=NLQ(nq) !is grid level for nq
C DPN 02 July 1999 - needs to be fixed !!!
c              IF(DOP) THEN
c                WRITE(OP_STRING,'('' nq='',I6,'' na='',I2)') nq,na
c                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              ENDIF
            ELSE
              na=1 !is fine grid level
            ENDIF

            IF(NWQ(1,nq).EQ.0 !internal grid pt
     '        .AND.(.NOT.ADAPTIVE
     '        .OR.ADAPTIVE.AND.NLQ(nq).GT.0
     '        .AND.NLQ(nq).LE.NMGT)) THEN !..& need residual

              IF(ADAPTIVE) THEN
C KAT 23Nov00: Check that NQGP is big enough
C               Get standard template sizes
                IF(NITB.EQ.1) THEN
                  NQGT=1+NXQ(-1,i,nq,na)+NXQ(-1,i,nq,na)
                ELSEIF(NITB.EQ.2) THEN
                  NQGT=9
                ELSE !IF(NITB.EQ.3) THEN
                  NQGT=27-8 ! no corners
                ENDIF
                IF(NQGT.GT.NQGM) THEN
                  IEND=0
                  CALL APPENDC(IEND,'NQGM needs to be at least ',ERROR)
                  CALL APPENDI(IEND,NQGT,ERROR)
                  GOTO 100
                ENDIF ! nqquadnum > NQGM
                IF(NITB.EQ.1) THEN
!one xi direction - may be branching
                  n=0
                  DO i=1,NXQ(-1,0,nq,na)
                    n=n+1
                    NQGP(n,nq)=NXQ(-1,i,nq,na)
                  ENDDO
                  n=n+1
                  NQGP(n,nq)=nq
                  DO i=1,NXQ(1,0,nq,na)
                    n=n+1
                    NQGP(n,nq)=NXQ( 1,i,nq,na)
                  ENDDO
                  NQGP(0,nq)=n

                ELSE IF(NITB.EQ.2) THEN
!two xi directions - branching?
                  BRANCH=.FALSE.
                  IF(NXQ(-1,0,NXQ(-2,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=1
                  ENDIF
                  IF(NXQ(-1,0,NXQ(2,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=2
                  ENDIF
                  IF(NXQ(1,0,NXQ(-2,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=3
                  ENDIF
                  IF(NXQ(1,0,NXQ(2,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=4
                  ENDIF
                  IF(NXQ(-2,0,NXQ(-1,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=5
                  ENDIF
                  IF(NXQ(-2,0,NXQ(1,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=6
                  ENDIF
                  IF(NXQ(2,0,NXQ(-1,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=7
                  ENDIF
                  IF(NXQ(2,0,NXQ(1,1,nq,na),na).GT.1) THEN
                    BRANCH=.TRUE.
                    CASE=8
                  ENDIF

                  IF(BRANCH) THEN
!more checks needed
                    IF(CASE.EQ.1) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(-1,0,NXQ(-2,1,nq,na),na)
                        nq2=NXQ(-1,n,NXQ(-2,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(1,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE IF(CASE.EQ.5) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(-2,0,NXQ(-1,1,nq,na),na)
                        nq2=NXQ(-2,n,NXQ(-1,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(1,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE
                      NQGP(1,nq)=NXQ(-1,1,NXQ(-2,1,nq,na),na)
                    ENDIF

                    IF(CASE.EQ.3) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(1,0,NXQ(-2,1,nq,na),na)
                        nq2=NXQ(1,n,NXQ(-2,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(3,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE IF(CASE.EQ.6) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(-2,0,NXQ(1,1,nq,na),na)
                        nq2=NXQ(-2,n,NXQ(1,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(3,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE
                      NQGP(3,nq)=NXQ( 1,1,NXQ(-2,1,nq,na),na)
                    ENDIF

                    IF(CASE.EQ.2) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(-1,0,NXQ(2,1,nq,na),na)
                        nq2=NXQ(-1,n,NXQ(2,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(7,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE IF(CASE.EQ.7) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(2,0,NXQ(-1,1,nq,na),na)
                        nq2=NXQ(2,n,NXQ(-1,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(7,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE
                      NQGP(7,nq)=NXQ(-1,1,NXQ(2,1,nq,na),na)
                    ENDIF

                    IF(CASE.EQ.4) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(1,0,NXQ(2,1,nq,na),na)
                        nq2=NXQ(1,n,NXQ(2,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(9,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE IF(CASE.EQ.8) THEN
                      FOUND=.FALSE.
                      DO n=1,NXQ(2,0,NXQ(1,1,nq,na),na)
                        nq2=NXQ(2,n,NXQ(1,1,nq,na),na)
                        DO i=1,NENQ(0,nq)
                          ne1=NENQ(i,nq)
                          DO j=1,NENQ(0,nq2)
                            ne2=NENQ(j,nq2)
                            IF(ne1.EQ.ne2) THEN
                              FOUND=.TRUE.
                              NQGP(9,nq)=nq2
                            ENDIF
                          ENDDO !j
                        ENDDO !i
                      ENDDO !n
                      IF(.NOT.FOUND) nqHOLD1=nq
                    ELSE
                      NQGP(9,nq)=NXQ( 1,1,NXQ(2,1,nq,na),na)
                    ENDIF

!point should be the same
                    NQGP(2,nq)=NXQ(-2,1,nq,na)
                    NQGP(4,nq)=NXQ(-1,1,nq,na)
                    NQGP(5,nq)=nq
                    NQGP(6,nq)=NXQ( 1,1,nq,na)
                    NQGP(8,nq)=NXQ( 2,1,nq,na)
                  ELSE !standard case
                    NQGP(1,nq)=NXQ(-1,1,NXQ(-2,1,nq,na),na)
                    NQGP(2,nq)=NXQ(-2,1,nq,na)
                    NQGP(3,nq)=NXQ( 1,1,NXQ(-2,1,nq,na),na)
                    NQGP(4,nq)=NXQ(-1,1,nq,na)
                    NQGP(5,nq)=nq
                    NQGP(6,nq)=NXQ( 1,1,nq,na)
                    NQGP(7,nq)=NXQ(-1,1,NXQ(2,1,nq,na),na)
                    NQGP(8,nq)=NXQ( 2,1,nq,na)
                    NQGP(9,nq)=NXQ( 1,1,NXQ(2,1,nq,na),na)
                  ENDIF !branches
                  NQGP(0,nq)=NQGT

                ELSE !IF(NITB.EQ.3) THEN
!three xi directions - no branching
                  NQGP(1,nq)=NXQ(-2,1,NXQ(-3,1,nq,na),na)
                  NQGP(2,nq)=NXQ(-1,1,NXQ(-3,1,nq,na),na)
                  NQGP(3,nq)=NXQ(-3,1,nq,na)
                  NQGP(4,nq)=NXQ( 1,1,NXQ(-3,1,nq,na),na)
                  NQGP(5,nq)=NXQ( 2,1,NXQ(-3,1,nq,na),na)
                  NQGP(6,nq)=NXQ(-1,1,NXQ(-2,1,nq,na),na)
                  NQGP(7,nq)=NXQ(-2,1,nq,na)
                  NQGP(8,nq)=NXQ( 1,1,NXQ(-2,1,nq,na),na)
                  NQGP(9,nq)=NXQ(-1,1,nq,na)
                  NQGP(10,nq)=nq
                  NQGP(11,nq)=NXQ( 1,1,nq,na)
                  NQGP(12,nq)=NXQ(-1,1,NXQ(2,1,nq,na),na)
                  NQGP(13,nq)=NXQ( 2,1,nq,na)
                  NQGP(14,nq)=NXQ( 1,1,NXQ(2,1,nq,na),na)
                  NQGP(15,nq)=NXQ(-2,1,NXQ(3,1,nq,na),na)
                  NQGP(16,nq)=NXQ(-1,1,NXQ(3,1,nq,na),na)
                  NQGP(17,nq)=NXQ( 3,1,nq,na)
                  NQGP(18,nq)=NXQ( 1,1,NXQ(3,1,nq,na),na)
                  NQGP(19,nq)=NXQ( 2,1,NXQ(3,1,nq,na),na)
                  NQGP(0,nq)=NQGT
                ENDIF !NITB
              ENDIF !ADAPTIVE
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP           CRITICAL(GET_FD_POINTS_2)
                WRITE(OP_STRING,'('' nq='',I6,'' NQGP:'',19I7)')
     '            nq,(NQGP(i,nq),i=1,19)
                CALL WRITES(IODI,OP_STRING,ERROR,*100)
CC$OMP           END CRITICAL(GET_FD_POINTS_2)
              ENDIF

            ELSE !boundary grid point

              NQGP(0,nq)=1
              NQGP(1,nq)=nq
              DO ni=1,NITB
                IF(NXQ(-ni,1,nq,na).EQ.0) THEN
!no grid point in -ni
!use one sided difference
                  nq1=NXQ(ni,1,nq,na)
                  nq2=NXQ(ni,1,nq1,na)

                  NQGP(0,nq)=NQGP(0,nq)+1
                  NQGP(NQGP(0,nq),nq)=nq1
                  NQGP(0,nq)=NQGP(0,nq)+1
                  NQGP(NQGP(0,nq),nq)=nq2
                ELSE IF(NXQ(ni,1,nq,na).EQ.0) THEN
!no grid point in +ni
!use one sided difference
                  nq1=NXQ(-ni,1,nq,na)
                  nq2=NXQ(-ni,1,nq1,na)

                  NQGP(0,nq)=NQGP(0,nq)+1
                  NQGP(NQGP(0,nq),nq)=nq1
                  NQGP(0,nq)=NQGP(0,nq)+1
                  NQGP(NQGP(0,nq),nq)=nq2
                ELSE
!points in both +ni and -ni
!we can use a two sided difference
                  nq1=NXQ(-ni,1,nq,na)
                  nq2=NXQ(ni,1,nq,na)

                  NQGP(0,nq)=NQGP(0,nq)+1
                  NQGP(NQGP(0,nq),nq)=nq1
                  NQGP(0,nq)=NQGP(0,nq)+1
                  NQGP(NQGP(0,nq),nq)=nq2
                ENDIF !NXQ non-zero
              ENDDO !ni

CMLB - for reference
C            NQGP(1,nq)=nq
C            NQGP(2,nq)=NWQ(1,nq)
C            NQGP(3,nq)=NWQ(2,nq)
C            NQGP(0,nq)=3

            ENDIF !internal/bdry point

            CALL ISORTP(NQGP(0,nq),NQGP(1,nq),NQGP_PIVOT(1,nq))

            DO ni=2,NQGP(0,nq)
              IF(NQGP(ni-1,nq).EQ.NQGP(ni,nq)) THEN
                ERROR='>>Duplicate matrix entries in GET_FD_POINTS'
                GOTO 100
              ENDIF
            ENDDO

C *** DPN 06 September 1999 - adding alternate return for MP
            GOTO 102
 100          CONTINUE
              ERROR_FLAG=.TRUE.
              IF(ERROR.NE.' ') THEN
                CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
              ENDIF
 102        CONTINUE

          ENDIF !.NOT.ERROR_FLAG

        ENDDO !nq
C$OMP END PARALLEL DO
        IF(ERROR_FLAG) THEN
          ERROR=' '
          GOTO 9999
        ENDIF
      ENDDO !nr

      IF(nqHOLD1.NE.0) THEN
        CALL ASSERT(.FALSE.,' >>A cusp point connection was not found',
     '    ERROR,*9999)
      ENDIF

      CALL EXITS('GET_FD_POINTS')
      RETURN
 9999 CALL ERRORS('GET_FD_POINTS',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('GET_FD_POINTS')
      RETURN 1
      END
