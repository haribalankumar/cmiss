      SUBROUTINE XQXE_REF(NENQ,NITB,nq,NQLIST,nr,NXQ,NXEDIM,AQ,XZE,
     '  ERROR,*)

C#### Subroutine: XQXE_REF
C###  Description:
C###    XQXE_REF transfers global parameters XQ to element parameters
C###    XE (geom,fibre,field). This routine transfers the stored 
C###    undeformed reference coordinates for mechanics problems. Note
C###    that lots of the checks have been removed as XQXE will always
C###    be called at least once before calling this routine and the 
C###    checks in there would be duplicated here. Also the returned
C###    variable is XZE(NSM,NXEDIM) so either XE or ZE can be passed
C###    to the routine dimensioned to NJM or NHM. 
C**** Written by Martin Buist, 23 Jan 2003

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'maqloc00.inc'

!     Parameter List
      INTEGER NENQ(0:8,NQM),NITB,nq,NQLIST(0:NQM),nr,NXEDIM,
     '  NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 AQ(NMAQM,NQM),XZE(NSM,NXEDIM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER AQNJ(3),CASE,CASES(4),i,icase,j,n,NCASE,ne1,ne2,nj,njj,
     '  nonq,nq2,ns,nsmax
      REAL*8 DS1(9),DS2(9)
      LOGICAL ACASE

      CALL ENTERS('XQXE_REF',*9999)

      ! Get the indices from AQ where the undeformed coordinates are
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,AQNJ(1),MAQ_X_UNDEF,
     '  ERROR,*9999)
      IF(AQNJ(1).EQ.0) THEN
        CALL ASSERT(.FALSE.,
     '    '>>Must update grid geometry deformed first',ERROR,*9999)
      ENDIF
      IF(NJT.GE.2) THEN
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,AQNJ(2),MAQ_Y_UNDEF,
     '    ERROR,*9999)
        IF(AQNJ(2).EQ.0) THEN
          CALL ASSERT(.FALSE.,
     '      '>>Must update grid geometry deformed first',ERROR,*9999)
        ENDIF
      ENDIF
      IF(NJT.GE.3) THEN
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,AQNJ(3),MAQ_Z_UNDEF,
     '    ERROR,*9999)
        IF(AQNJ(3).EQ.0) THEN
          CALL ASSERT(.FALSE.,
     '      '>>Must update grid geometry deformed first',ERROR,*9999)
        ENDIF
      ENDIF

      ! Initialisation
      IF(NITB.EQ.1) THEN
        CALL ASSERT(NSM.GE.3,'>>Increase NSM >= 3',ERROR,*9999)
        nsmax=3
      ELSE IF(NITB.EQ.2) THEN
        CALL ASSERT(NSM.GE.9,'>>Increase NSM >= 9',ERROR,*9999)
        nsmax=9
      ELSE IF(NITB.EQ.3) THEN
        CALL ASSERT(NSM.GE.27,'>>Increase NSM >= 27',ERROR,*9999)
        nsmax=27
      ENDIF

      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        nj=NJ_LOC(NJL_GEOM,njj,nr)
        DO ns=1,nsmax
          XZE(ns,nj)=0.0d0
        ENDDO !ns
      ENDDO !njj

      ! Computation of XE
      IF(NITB.EQ.1) THEN !1 xi direction

        !Average point in -xi 1
        DO n=1,NXQ(-1,0,nq)
          DS1(n)=0.0d0
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DS1(n)=DS1(n)+(AQ(AQNJ(nj),NXQ(-1,n,nq))-
     '        AQ(AQNJ(nj),nq))**2.0d0
          ENDDO !njj
          DS1(n)=DSQRT(DS1(n))
        ENDDO !n
        DO n=2,NXQ(-1,0,nq)
          DS1(1)=DS1(1)+DS1(n)
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.0) THEN
          DS1(1)=DS1(1)/DBLE(NXQ(-1,0,nq))
          XZE(1,1)=AQ(AQNJ(1),nq)-DS1(1)
        ENDIF !denom

        !Centre point
        XZE(2,1)=AQ(AQNJ(1),nq)

        !Average point in +xi 1
        DO n=1,NXQ(1,0,nq)
          DS2(n)=0.0d0
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DS2(n)=DS2(n)+(AQ(AQNJ(nj),NXQ(1,n,nq))
     '        -AQ(AQNJ(nj),nq))**2.0d0
          ENDDO !njj
          DS2(n)=DSQRT(DS2(n))
        ENDDO !n
        DO n=2,NXQ(1,0,nq)
          DS2(1)=DS2(1)+DS2(n)
        ENDDO !n
        IF(NXQ(1,0,nq).GT.0) THEN
          DS2(1)=DS2(1)/DBLE(NXQ(1,0,nq))
          XZE(3,1)=AQ(AQNJ(1),nq)+DS2(1)
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
                  NQLIST(1)=nq2
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
                  NQLIST(1)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          NQLIST(1)=NXQ(-1,1,NXQ(-2,1,nq))
        ENDIF
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(1,nj)=AQ(AQNJ(nj),NQLIST(1))
        ENDDO !njj

        !Average point in -xi 2
        DO n=1,NXQ(-2,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(2,nj)=XZE(2,nj)+AQ(AQNJ(nj),NXQ(-2,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(-2,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(2,nj)=XZE(2,nj)/DBLE(NXQ(-2,0,nq))
          ENDDO !njj
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
                  NQLIST(3)=nq2
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
                  NQLIST(3)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          NQLIST(3)=NXQ(1,1,NXQ(-2,1,nq))
        ENDIF
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(3,nj)=AQ(AQNJ(nj),NQLIST(3))
        ENDDO !njj

        !Average point in -xi 1
        DO n=1,NXQ(-1,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(4,nj)=XZE(4,nj)+AQ(AQNJ(nj),NXQ(-1,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(4,nj)=XZE(4,nj)/DBLE(NXQ(-1,0,nq))
          ENDDO !njj
        ENDIF !denom

        !Centre point
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(5,nj)=AQ(AQNJ(nj),nq)
        ENDDO !njj

        !Average point in +xi 1
        DO n=1,NXQ(1,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(6,nj)=XZE(6,nj)+AQ(AQNJ(nj),NXQ(1,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(1,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(6,nj)=XZE(6,nj)/DBLE(NXQ(1,0,nq))
          ENDDO !njj
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
                  NQLIST(7)=nq2
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
                  NQLIST(7)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          NQLIST(7)=NXQ(-1,1,NXQ(2,1,nq))
        ENDIF
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(7,nj)=AQ(AQNJ(nj),NQLIST(7))
        ENDDO !njj

        !Average point in +xi 2
        DO n=1,NXQ(2,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(8,nj)=XZE(8,nj)+AQ(AQNJ(nj),NXQ(2,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(2,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(8,nj)=XZE(8,nj)/DBLE(NXQ(2,0,nq))
          ENDDO !njj
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
                  NQLIST(9)=nq2
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
                  NQLIST(9)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
        ELSE
          NQLIST(9)=NXQ(1,1,NXQ(2,1,nq))  
        ENDIF
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(9,nj)=AQ(AQNJ(nj),NQLIST(9))
        ENDDO !njj

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
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq))).GT.0) THEN
            XZE(1,nj)=AQ(AQNJ(nj),NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq))))
          ENDIF
        ENDDO !njj
        
        ! Grid point 2
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-2,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(-3,1,nq))
          ELSEIF(NXQ(-3,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(-2,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(2,nj)=XZE(2,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(2,nj)=XZE(2,nj)/DBLE(NQLIST(0))
        ENDDO !njj
        
        ! Grid point 3
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq))).GT.0) THEN
            XZE(3,nj)=AQ(AQNJ(nj),NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq))))
          ENDIF
        ENDDO !njj
        
        ! Grid point 4
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF 
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(-1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(-3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(4,nj)=XZE(4,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(4,nj)=XZE(4,nj)/DBLE(NQLIST(0))
        ENDDO !njj

        ! Grid point 5 - average point in -xi3
        DO n=1,NXQ(-3,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(5,nj)=XZE(5,nj)+AQ(AQNJ(nj),NXQ(-3,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(-3,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(5,nj)=XZE(5,nj)/DBLE(NXQ(-3,0,nq))
          ENDDO !njj
        ENDIF !denom
        
        ! Grid point 6
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(1,1,nq))
          ELSEIF(NXQ(1,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(-3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(6,nj)=XZE(6,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(6,nj)=XZE(6,nj)/DBLE(NQLIST(0))
        ENDDO !njj

        ! Grid point 7
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq))).GT.0) THEN
            XZE(7,nj)=AQ(AQNJ(nj),NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq))))
          ENDIF
        ENDDO !njj
        
        ! Grid point 8
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(-3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(8,nj)=XZE(8,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(8,nj)=XZE(8,nj)/DBLE(NQLIST(0))
        ENDDO !njj
        
        ! Grid point 9
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq))).GT.0) THEN
            XZE(9,nj)=AQ(AQNJ(nj),NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq))))
          ENDIF
        ENDDO !njj
        
        ! Grid point 10
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-2,1,NXQ(-1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(-2,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(10,nj)=XZE(10,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(10,nj)=XZE(10,nj)/DBLE(NQLIST(0))
        ENDDO !njj
          
        !Grid point 11 take average in -xi2
        DO n=1,NXQ(-2,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(11,nj)=XZE(11,nj)+AQ(AQNJ(nj),NXQ(-2,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(-2,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(11,nj)=XZE(11,nj)/DBLE(NXQ(-2,0,nq))
          ENDDO !njj
        ENDIF !denom
        
        ! Grid point 12
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(1,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(-2,1,nq))
          ELSEIF(NXQ(-2,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(1,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(12,nj)=XZE(12,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(12,nj)=XZE(12,nj)/DBLE(NQLIST(0))
        ENDDO !njj
          
        ! Grid point 13 - average -xi1
        DO n=1,NXQ(-1,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(13,nj)=XZE(13,nj)+AQ(AQNJ(nj),NXQ(-1,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(13,nj)=XZE(13,nj)/DBLE(NXQ(-1,0,nq))
          ENDDO !njj
        ENDIF !denom
        
        !Grid point 14 - centre point
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(14,nj)=AQ(AQNJ(nj),nq)
        ENDDO !njj

        !Grid point 15 average xi1
        DO n=1,NXQ(1,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(15,nj)=XZE(15,nj)+AQ(AQNJ(nj),NXQ(1,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(1,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(15,nj)=XZE(15,nj)/DBLE(NXQ(1,0,nq))
          ENDDO !njj
        ENDIF !denom
        
        !Grid point 16
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-1,1,NXQ(2,1,nq)).GT.0)THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(-1,1,nq)).GT.0)THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(-1,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(16,nj)=XZE(16,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(16,nj)=XZE(16,nj)/DBLE(NQLIST(0))
        ENDDO !njj
          
        ! Grid point 17 average xi2
        DO n=1,NXQ(2,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(17,nj)=XZE(17,nj)+AQ(AQNJ(nj),NXQ(2,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(2,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(17,nj)=XZE(17,nj)/DBLE(NXQ(2,0,nq))
          ENDDO !njj
        ENDIF !denom
        
        !Grid point 18
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(1,1,NXQ(2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(1,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(18,nj)=XZE(18,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(18,nj)=XZE(18,nj)/DBLE(NQLIST(0))
        ENDDO !njj
       
        !Grid point 19
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq))).GT.0) THEN
            XZE(19,nj)=AQ(AQNJ(nj),NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq))))
          ENDIF
        ENDDO !njj
        
        !Grid point 20
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(-2,1,nq))
          ELSEIF(NXQ(-2,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(20,nj)=XZE(20,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(20,nj)=XZE(20,nj)/DBLE(NQLIST(0))
        ENDDO !njj
        
        !Grid point 21
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq))).GT.0) THEN
            XZE(21,nj)=AQ(AQNJ(nj),NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq))))
          ENDIF
        ENDDO !njj
        
        !Grid point 22
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(-1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(22,nj)=XZE(22,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(22,nj)=XZE(22,nj)/DBLE(NQLIST(0))
        ENDDO !njj

        ! Grid point 23 average xi3
        DO n=1,NXQ(3,0,nq)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(23,nj)=XZE(23,nj)+AQ(AQNJ(nj),NXQ(3,n,nq))
          ENDDO !njj
        ENDDO !n
        IF(NXQ(3,0,nq).GT.0) THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(23,nj)=XZE(23,nj)/DBLE(NXQ(3,0,nq))
          ENDDO !njj
        ENDIF !denom

        !Grid point 24
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(1,1,nq))
          ELSEIF(NXQ(1,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(24,nj)=XZE(24,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(24,nj)=XZE(24,nj)/DBLE(NQLIST(0))
        ENDDO !njj
        
        !Grid point 25
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq))).GT.0) THEN
            XZE(25,nj)=AQ(AQNJ(nj),NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq))))
          ENDIF
        ENDDO !njj
        
        !Grid point 26
        ACASE=.FALSE.
        NQLIST(0)=0
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
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
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(3,1,nq))
          ENDIF
        ENDIF
        DO nonq=1,NQLIST(0)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XZE(26,nj)=XZE(26,nj)+AQ(AQNJ(nj),NQLIST(nonq))
          ENDDO !njj
        ENDDO !nonq
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          XZE(26,nj)=XZE(26,nj)/DBLE(NQLIST(0))
        ENDDO !njj
      
        !Grid point 27
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          IF(NXQ(3,1,NXQ(2,1,NXQ(1,1,nq))).GT.0) THEN
            XZE(27,nj)=AQ(AQNJ(nj),NXQ(3,1,NXQ(2,1,NXQ(1,1,nq))))
          ENDIF
        ENDDO !njj

      ENDIF !nit(nb)

      CALL EXITS('XQXE_REF')
      RETURN
 9999 CALL ERRORS('XQXE_REF',ERROR)
      CALL EXITS('XQXE_REF')
      RETURN 1
      END


