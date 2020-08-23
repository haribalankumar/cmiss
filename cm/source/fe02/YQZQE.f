      SUBROUTINE YQZQE(NENQ,NITB,nq,NXQ,ZQEDIM,YQ,ZQE,ERROR,*)

C#### Subroutine: YQZQE
C###  Description:
C###    YQZQE is the grid point equivalent of YPZP+ZPZE and transfers
C###    dependent variables from the YQ array to a ZQE array which
C###    is a local quadratic patch. 
C**** Written by Martin Buist, 5 March 2003

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'maqloc00.inc'

!     Parameter List
      INTEGER NENQ(0:8,NQM),NITB,nq,NXQ(-NIM:NIM,0:4,0:NQM),ZQEDIM
      REAL*8 YQ(NYQM),ZQE(ZQEDIM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CASE,CASES(4),i,icase,j,LNQLIST(0:99),n,NCASE,ne1,ne2,
     '  nonq,nq1,nq2
      LOGICAL ACASE

      CALL ENTERS('YQZQE',*9999)

      IF(NITB.EQ.1) THEN
        CALL ASSERT(ZQEDIM.GE.3,'>>Increase ZQE array >= 3',ERROR,*9999)
      ELSE IF(NITB.EQ.2) THEN
        CALL ASSERT(ZQEDIM.GE.9,'>>Increase ZQE array >= 9',ERROR,*9999)
      ELSE IF(NITB.EQ.3) THEN
        CALL ASSERT(ZQEDIM.GE.27,'>>Increase ZQE array >= 27',
     '    ERROR,*9999)
      ENDIF

      ! Computation of ZQE
      IF(NITB.EQ.1) THEN !1 xi direction

        !Average point in -xi 1
        ZQE(1)=0.0d0
        DO n=1,NXQ(-1,0,nq)
          ZQE(1)=ZQE(1)+YQ(NXQ(-1,n,nq))
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.1) THEN
          ZQE(1)=ZQE(1)/DBLE(NXQ(-1,0,nq))
        ENDIF !denom

        !Centre point
        ZQE(2)=YQ(nq)

        !Average point in +xi 1
        ZQE(3)=0.0d0
        DO n=1,NXQ(1,0,nq)
          ZQE(3)=ZQE(3)+YQ(NXQ(1,n,nq))
        ENDDO !n
        IF(NXQ(1,0,nq).GT.1) THEN
          ZQE(3)=ZQE(3)/DBLE(NXQ(1,0,nq))
        ENDIF !denom

      ELSE IF(NITB.EQ.2) THEN !2 xi directions

        !Check to see if corners are next to cusps
        CASE=0
        IF(NXQ(-1,0,NXQ(-2,1,nq)).GT.1) THEN
          CASE=1
        ENDIF
        IF(NXQ(-1,0,NXQ(2,1,nq)).GT.1) THEN
          CASE=2
        ENDIF
        IF(NXQ(1,0,NXQ(-2,1,nq)).GT.1) THEN
          CASE=3
        ENDIF
        IF(NXQ(1,0,NXQ(2,1,nq)).GT.1) THEN
          CASE=4
        ENDIF
        IF(NXQ(-2,0,NXQ(-1,1,nq)).GT.1) THEN
          CASE=5
        ENDIF
        IF(NXQ(-2,0,NXQ(1,1,nq)).GT.1) THEN
          CASE=6
        ENDIF
        IF(NXQ(2,0,NXQ(-1,1,nq)).GT.1) THEN
          CASE=7
        ENDIF
        IF(NXQ(2,0,NXQ(1,1,nq)).GT.1) THEN
          CASE=8
        ENDIF

        !Point in -xi 1 & -xi 2
        IF(CASE.EQ.1) THEN
          DO n=1,NXQ(-1,0,NXQ(-2,1,nq))
            nq2=NXQ(-1,n,NXQ(-2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE IF(CASE.EQ.5) THEN
          DO n=1,NXQ(-2,0,NXQ(-1,1,nq))
            nq2=NXQ(-2,n,NXQ(-1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          nq1=NXQ(-1,1,NXQ(-2,1,nq))
        ENDIF
        ZQE(1)=YQ(nq1)

        !Average point in -xi 2
        ZQE(2)=0.0d0
        DO n=1,NXQ(-2,0,nq)
          ZQE(2)=ZQE(2)+YQ(NXQ(-2,n,nq))
        ENDDO !n
        IF(NXQ(-2,0,nq).GT.1) THEN
          ZQE(2)=ZQE(2)/DBLE(NXQ(-2,0,nq))
        ENDIF !denom

        !Point in +xi 1 & -xi 2
        IF(CASE.EQ.3) THEN
          DO n=1,NXQ(1,0,NXQ(-2,1,nq))
            nq2=NXQ(1,n,NXQ(-2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE IF(CASE.EQ.6) THEN
          DO n=1,NXQ(-2,0,NXQ(1,1,nq))
            nq2=NXQ(-2,n,NXQ(1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          nq1=NXQ(1,1,NXQ(-2,1,nq))
        ENDIF
        ZQE(3)=YQ(nq1)

        !Average point in -xi 1
        ZQE(4)=0.0d0
        DO n=1,NXQ(-1,0,nq)
          ZQE(4)=ZQE(4)+YQ(NXQ(-1,n,nq))
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.1) THEN
          ZQE(4)=ZQE(4)/DBLE(NXQ(-1,0,nq))
        ENDIF !denom

        !Centre point
        ZQE(5)=YQ(nq)

        !Average point in +xi 1
        ZQE(6)=0.0d0
        DO n=1,NXQ(1,0,nq)
          ZQE(6)=ZQE(6)+YQ(NXQ(1,n,nq))
        ENDDO !n
        IF(NXQ(1,0,nq).GT.1) THEN
          ZQE(6)=ZQE(6)/DBLE(NXQ(1,0,nq))
        ENDIF !denom

        !Point in -xi1 & +xi 2
        IF(CASE.EQ.2) THEN
          DO n=1,NXQ(-1,0,NXQ(2,1,nq))
            nq2=NXQ(-1,n,NXQ(2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE IF(CASE.EQ.7) THEN
          DO n=1,NXQ(2,0,NXQ(-1,1,nq))
            nq2=NXQ(2,n,NXQ(-1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          nq1=NXQ(-1,1,NXQ(2,1,nq))
        ENDIF
        ZQE(7)=YQ(nq1)

        !Average point in +xi 2
        ZQE(8)=0.0d0
        DO n=1,NXQ(2,0,nq)
          ZQE(8)=ZQE(8)+YQ(NXQ(2,n,nq))
        ENDDO !n
        IF(NXQ(2,0,nq).GT.1) THEN
          ZQE(8)=ZQE(8)/DBLE(NXQ(2,0,nq))
        ENDIF !denom

        !Point in +xi 1 & +xi 2
        IF(CASE.EQ.4) THEN
          DO n=1,NXQ(1,0,NXQ(2,1,nq))
            nq2=NXQ(1,n,NXQ(2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE IF(CASE.EQ.8) THEN                  
          DO n=1,NXQ(2,0,NXQ(1,1,nq))
            nq2=NXQ(2,n,NXQ(1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  nq1=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          nq1=NXQ(1,1,NXQ(2,1,nq))  
        ENDIF
        ZQE(9)=YQ(nq1)

      ELSE IF(NITB.EQ.3) THEN

        NCASE=0

        ! Node 2
        IF(NXQ(-2,0,NXQ(-3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=1
        ENDIF
        IF(NXQ(-3,0,NXQ(-2,1,nq)).GT.1) THEN
          NCASE=NCASE+1          
          CASES(NCASE)=2
        ENDIF

        ! Node 4
        IF(NXQ(-3,0,NXQ(-1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=3
        ENDIF
        IF(NXQ(-1,0,NXQ(-3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=4
        ENDIF

        ! Node 6
        IF(NXQ(-3,0,NXQ(1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=5
        ENDIF
        IF(NXQ(1,0,NXQ(-3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=6
        ENDIF

        ! Node 8
        IF(NXQ(-3,0,NXQ(2,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=7
        ENDIF
        IF(NXQ(2,0,NXQ(-3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=8
        ENDIF

        ! Node 10
        IF(NXQ(-2,0,NXQ(-1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=9
        ENDIF
        IF(NXQ(-1,0,NXQ(-2,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=10
        ENDIF

        ! Node 12
        IF(NXQ(1,0,NXQ(-2,1,nq)).GT.1) THEN
          NCASE=NCASE+1
         CASES(NCASE)=11
        ENDIF
        IF(NXQ(-2,0,NXQ(1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=12
        ENDIF

        ! Node 16
        IF(NXQ(-1,0,NXQ(2,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=13
        ENDIF
        IF(NXQ(2,0,NXQ(-1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=14
        ENDIF

        ! Node 18
        IF(NXQ(1,0,NXQ(2,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=15
        ENDIF
        IF(NXQ(2,0,NXQ(1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=16
        ENDIF

        ! Node 20
        IF(NXQ(3,0,NXQ(-2,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=17
        ENDIF
        IF(NXQ(-2,0,NXQ(3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=18
        ENDIF

        ! Node 22
        IF(NXQ(3,0,NXQ(-1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=19
        ENDIF
        IF(NXQ(-1,0,NXQ(3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=20
        ENDIF

        ! Node 24
        IF(NXQ(3,0,NXQ(1,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=21
        ENDIF
        IF(NXQ(1,0,NXQ(3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=22
        ENDIF

        ! Node 26
        IF(NXQ(3,0,NXQ(2,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=23
        ENDIF
        IF(NXQ(2,0,NXQ(3,1,nq)).GT.1) THEN
          NCASE=NCASE+1
          CASES(NCASE)=24
        ENDIF

        ! Grid point 1
        IF(NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq))).EQ.0) THEN
          ZQE(1)=0.0d0
        ELSE
          ZQE(1)=YQ(NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq))))
        ENDIF
        
        ! Grid point 2
        ZQE(2)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.1) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-2,0,NXQ(-3,1,nq))
              nq2=NXQ(-2,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.2) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-3,0,NXQ(-2,1,nq))
              nq2=NXQ(-3,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-2,1,NXQ(-3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-2,1,NXQ(-3,1,nq))
          ELSEIF(NXQ(-3,1,NXQ(-2,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-3,1,NXQ(-2,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(2)=ZQE(2)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(2)=ZQE(2)/DBLE(LNQLIST(0))
        ENDIF
        
        ! Grid point 3
        IF(NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq))).EQ.0) THEN
          ZQE(3)=0.0d0
        ELSE
          ZQE(3)=YQ(NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq))))
        ENDIF
        
        ! Grid point 4
        ZQE(4)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.3) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-3,0,NXQ(-1,1,nq))
              nq2=NXQ(-3,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.4) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-1,0,NXQ(-3,1,nq))
              nq2=NXQ(-1,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF 
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(-1,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-3,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(-3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-1,1,NXQ(-3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(4)=ZQE(4)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(4)=ZQE(4)/DBLE(LNQLIST(0))
        ENDIF

        ! Grid point 5 - average point in -xi3
        ZQE(5)=0.0d0
        DO n=1,NXQ(-3,0,nq)
          ZQE(5)=ZQE(5)+YQ(NXQ(-3,n,nq))
        ENDDO !n
        IF(NXQ(-3,0,nq).GT.0) THEN
          ZQE(5)=ZQE(5)/DBLE(NXQ(-3,0,nq))
        ENDIF !denom
        
        ! Grid point 6
        ZQE(6)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.5) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-3,0,NXQ(1,1,nq))
              nq2=NXQ(-3,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.6) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(1,0,NXQ(-3,1,nq))
              nq2=NXQ(1,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(1,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-3,1,NXQ(1,1,nq))
          ELSEIF(NXQ(1,1,NXQ(-3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(1,1,NXQ(-3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(6)=ZQE(6)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(6)=ZQE(6)/DBLE(LNQLIST(0))
        ENDIF

        ! Grid point 7
        IF(NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq))).EQ.0) THEN
          ZQE(7)=0.0d0
        ELSE
          ZQE(7)=YQ(NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq))))
        ENDIF
        
        ! Grid point 8
        ZQE(8)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.7) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-3,0,NXQ(2,1,nq))
              nq2=NXQ(-3,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.8) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(2,0,NXQ(-3,1,nq))
              nq2=NXQ(2,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(2,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-3,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(-3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(2,1,NXQ(-3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(8)=ZQE(8)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(8)=ZQE(8)/DBLE(LNQLIST(0))
        ENDIF

        ! Grid point 9
        IF(NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq))).EQ.0) THEN
          ZQE(9)=0.0d0
        ELSE
          ZQE(9)=YQ(NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq))))
        ENDIF
        
        ! Grid point 10
        ZQE(10)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.9) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-2,0,NXQ(-1,1,nq))
              nq2=NXQ(-2,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.10) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-1,0,NXQ(-2,1,nq))
              nq2=NXQ(-1,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-2,1,NXQ(-1,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-2,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(-2,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-1,1,NXQ(-2,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(10)=ZQE(10)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(10)=ZQE(10)/DBLE(LNQLIST(0))
        ENDIF
          
        !Grid point 11 take average in -xi2
        ZQE(11)=0.0d0
        DO n=1,NXQ(-2,0,nq)
          ZQE(11)=ZQE(11)+YQ(NXQ(-2,n,nq))
        ENDDO !n
        IF(NXQ(-2,0,nq).GT.0) THEN
          ZQE(11)=ZQE(11)/DBLE(NXQ(-2,0,nq))
        ENDIF !denom
        
        ! Grid point 12
        ZQE(12)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.11) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(1,0,NXQ(-2,1,nq))
              nq2=NXQ(1,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.12) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-2,0,NXQ(1,1,nq))
              nq2=NXQ(-2,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(1,1,NXQ(-2,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(1,1,NXQ(-2,1,nq))
          ELSEIF(NXQ(-2,1,NXQ(1,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-2,1,NXQ(1,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(12)=ZQE(12)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(12)=ZQE(12)/DBLE(LNQLIST(0))
        ENDIF

        ! Grid point 13 - average -xi1
        ZQE(13)=0.0d0
        DO n=1,NXQ(-1,0,nq)
          ZQE(13)=ZQE(13)+YQ(NXQ(-1,n,nq))
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.0) THEN
          ZQE(13)=ZQE(13)/DBLE(NXQ(-1,0,nq))
        ENDIF !denom
        
        !Grid point 14 - centre point
        ZQE(14)=YQ(nq)

        !Grid point 15 average xi1
        ZQE(15)=0.0d0
        DO n=1,NXQ(1,0,nq)
          ZQE(15)=ZQE(15)+YQ(NXQ(1,n,nq))
        ENDDO !n
        IF(NXQ(1,0,nq).GT.0) THEN
          ZQE(15)=ZQE(15)/DBLE(NXQ(1,0,nq))
        ENDIF !denom
        
        !Grid point 16
        ZQE(16)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.13) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-1,0,NXQ(2,1,nq))
              nq2=NXQ(-1,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2))
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.14) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(2,0,NXQ(-1,1,nq))
              nq2=NXQ(2,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-1,1,NXQ(2,1,nq)).GT.0)THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-1,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(-1,1,nq)).GT.0)THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(2,1,NXQ(-1,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(16)=ZQE(16)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(16)=ZQE(16)/DBLE(LNQLIST(0))
        ENDIF

        ! Grid point 17 average xi2
        ZQE(17)=0.0d0
        DO n=1,NXQ(2,0,nq)
          ZQE(17)=ZQE(17)+YQ(NXQ(2,n,nq))
        ENDDO !n
        IF(NXQ(2,0,nq).GT.0) THEN
          ZQE(17)=ZQE(17)/DBLE(NXQ(2,0,nq))
        ENDIF !denom
        
        !Grid point 18
        ZQE(18)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.15) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(1,0,NXQ(2,1,nq))
              nq2=NXQ(1,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.16) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(2,0,NXQ(1,1,nq))
              nq2=NXQ(2,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(1,1,NXQ(2,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(1,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(1,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(2,1,NXQ(1,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(18)=ZQE(18)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(18)=ZQE(18)/DBLE(LNQLIST(0))
        ENDIF

        !Grid point 19
        IF(NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq))).EQ.0) THEN
          ZQE(19)=0.0d0
        ELSE
          ZQE(19)=YQ(NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq))))
        ENDIF
        
        !Grid point 20
        ZQE(20)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.17) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(3,0,NXQ(-2,1,nq))
              nq2=NXQ(3,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.18) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(3,0,NXQ(-2,1,nq))
              nq2=NXQ(3,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(-2,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(3,1,NXQ(-2,1,nq))
          ELSEIF(NXQ(-2,1,NXQ(3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-2,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(20)=ZQE(20)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(20)=ZQE(20)/DBLE(LNQLIST(0))
        ENDIF

        !Grid point 21
        IF(NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq))).EQ.0) THEN
          ZQE(21)=0.0d0
        ELSE
          ZQE(21)=YQ(NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq))))
        ENDIF
        
        !Grid point 22
        ZQE(22)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.19) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(3,0,NXQ(-1,1,nq))
              nq2=NXQ(3,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.20) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(-1,0,NXQ(3,1,nq))
              nq2=NXQ(-1,n,NXQ(3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(-1,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(3,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(-1,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(22)=ZQE(22)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(22)=ZQE(22)/DBLE(LNQLIST(0))
        ENDIF

        ! Grid point 23 average xi3
        ZQE(23)=0.0d0
        DO n=1,NXQ(3,0,nq)
          ZQE(23)=ZQE(23)+YQ(NXQ(3,n,nq))
        ENDDO !n
        IF(NXQ(3,0,nq).GT.0) THEN
          ZQE(23)=ZQE(23)/DBLE(NXQ(3,0,nq))
        ENDIF !denom

        !Grid point 24
        ZQE(24)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.21) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(3,0,NXQ(1,1,nq))
              nq2=NXQ(3,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.22) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(1,0,NXQ(3,1,nq))
              nq2=NXQ(1,n,NXQ(3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(1,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(3,1,NXQ(1,1,nq))
          ELSEIF(NXQ(1,1,NXQ(3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(1,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(24)=ZQE(24)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(24)=ZQE(24)/DBLE(LNQLIST(0))
        ENDIF

        !Grid point 25
        IF(NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq))).EQ.0) THEN
          ZQE(25)=0.0d0
        ELSE
          ZQE(25)=YQ(NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq))))
        ENDIF
        
        !Grid point 26
        ZQE(26)=0.0d0
        ACASE=.FALSE.
        LNQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.23) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(3,0,NXQ(2,1,nq))
              nq2=NXQ(3,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ELSEIF(CASE.EQ.24) THEN
            ACASE=.TRUE.
            DO n=1,NXQ(2,0,NXQ(3,1,nq))
              nq2=NXQ(2,n,NXQ(3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    LNQLIST(0)=LNQLIST(0)+1
                    LNQLIST(LNQLIST(0))=nq2
                    IF(LNQLIST(0).EQ.2) THEN
                      IF(LNQLIST(1).EQ.LNQLIST(2)) 
     '                  LNQLIST(0)=LNQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(2,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(3,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(3,1,nq)).GT.0) THEN
            LNQLIST(0)=LNQLIST(0)+1
            LNQLIST(LNQLIST(0))=NXQ(2,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,LNQLIST(0)
          ZQE(26)=ZQE(26)+YQ(LNQLIST(nonq))
        ENDDO !nonq
        IF(LNQLIST(0).GT.1) THEN
          ZQE(26)=ZQE(26)/DBLE(LNQLIST(0))
        ENDIF

        !Grid point 27
        IF(NXQ(3,1,NXQ(2,1,NXQ(1,1,nq))).EQ.0) THEN
          ZQE(27)=0.0d0
        ELSE
          ZQE(27)=YQ(NXQ(3,1,NXQ(2,1,NXQ(1,1,nq))))
        ENDIF
      ENDIF !nitb

      CALL EXITS('YQZQE')
      RETURN
 9999 CALL ERRORS('YQZQE',ERROR)
      CALL EXITS('YQZQE')
      RETURN 1
      END


