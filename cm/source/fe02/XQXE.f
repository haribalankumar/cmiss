      SUBROUTINE XQXE(NBJ,NENQ,nq,NQGP,NQLIST,nr,NXQ,XE,XQ,DEFORMED,
     '  ERROR,*)

C#### Subroutine: XQXE
C###  Description:
C###    XQXE transfers global parameters XQ to element parameters
C###    XE (geom,fibre,field). The XE values are generated for a
C###    local quadratic element in XYZ space and should not be used
C###    for global finite element interpolation.
C**** Written by Martin Buist, 21 June 1999

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'

!     Parameter List
      INTEGER NBJ(NJM),NENQ(0:8,NQM),nq,NQGP(0:NQGM,NQM),NQLIST(0:NQM),
     '  nr,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 XE(NSM,NJM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL DEFORMED
!     Local Variables
      INTEGER CASE,CASES(4),i,icase,j,n,NCASE,ne1,ne2,ni,NITB,
     '  nj,njj1,njj2,njjstep,nonq,nq2,nqgt,nqHOLD,nqquadnum
      REAL*8 DS1(5),DS2(5)
      LOGICAL ACASE,FOUND

      CALL ENTERS('XQXE',*9999)

      nqHOLD=0
      NJJSTEP=1 !Loop over geometry,fibres,field
      NITB=NIT(NBJ(1))
      IF(.NOT.DEFORMED) THEN
C KAT 23Nov00: Check that NQGP is big enough
C       Get standard template sizes
        IF(NITB.EQ.1) THEN
          nqgt=3
        ELSEIF(NITB.EQ.2) THEN
          nqgt=9
        ELSE !IF(NITB.EQ.3) THEN
          nqgt=27-8 ! no corners
        ENDIF
C       Check grid for branches at adjacent grid points.
C       (Branches at diagonally adjacent points are checked later.)
        DO ni=1,NITB
          nqgt=nqgt+NXQ(-ni,0,nq)+NXQ(ni,0,nq)-2 !branches
        ENDDO !ni
        IF(nqgt.GT.NQGM) THEN
          WRITE(ERROR,'(''NQGM needs to be at least1 '',I12)') nqgt
C          IEND=0
C          CALL APPENDC(IEND,'NQGM needs to be at least ',ERROR)
C          CALL APPENDI(IEND,nqgt,ERROR)
          GOTO 9999
        ENDIF ! nqquadnum > NQGM
C KAT 23Nov00: I don't think this needs to be initialized.
C        DO i=0,22
C          NQGP(i,nq)=0
C        ENDDO !i
      ENDIF !deformed

      IF(NITB.EQ.1) THEN !1 xi direction
        nqquadnum=0
        !Average point in -xi 1
        DO nj=1,NJM
          XE(1,nj)=0.0d0
        ENDDO !nj
        DO n=1,NXQ(-1,0,nq)
          DS1(n)=0.0d0
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              DS1(n)=DS1(n)+(XQ(nj,NXQ(-1,n,nq))-XQ(nj,nq))**2.0d0
            ENDDO !njj2
          ENDDO !njj1
          DS1(n)=DSQRT(DS1(n))
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(-1,n,nq)
          ENDIF !deformed
        ENDDO !n
        DO n=2,NXQ(-1,0,nq)
          DS1(1)=DS1(1)+DS1(n)
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.0) THEN
          DS1(1)=DS1(1)/DBLE(NXQ(-1,0,nq))
          XE(1,1)=XQ(1,nq)-DS1(1)
        ENDIF !denom

        !Centre point
        DO nj=1,NJM
          XE(2,nj)=0.0d0
        ENDDO !nj
        XE(2,1)=XQ(1,nq)
        IF(.NOT.DEFORMED) THEN
          nqquadnum=nqquadnum+1
          NQGP(nqquadnum,nq)=nq
        ENDIF !deformed

        !Average point in +xi 1
        DO nj=1,NJM
          XE(3,nj)=0.0d0
        ENDDO !nj
        DO n=1,NXQ(1,0,nq)
          DS2(n)=0.0d0
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              DS2(n)=DS2(n)+(XQ(nj,NXQ(1,n,nq))-XQ(nj,nq))**2.0d0
            ENDDO !njj2
          ENDDO !njj1
          DS2(n)=DSQRT(DS2(n))
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(1,n,nq)
          ENDIF !deformed
        ENDDO !n
        DO n=2,NXQ(1,0,nq)
          DS2(1)=DS2(1)+DS2(n)
        ENDDO !n
        IF(NXQ(1,0,nq).GT.0) THEN
          DS2(1)=DS2(1)/DBLE(NXQ(1,0,nq))
          XE(3,1)=XQ(1,nq)+DS2(1)
        ENDIF !denom

        IF(.NOT.DEFORMED) NQGP(0,nq)=nqquadnum

      ELSE IF(NITB.EQ.2) THEN !2 xi directions

C KAT 11Oct00: Initializing
        CASE=0
        !Check to see if corners are next to cusps
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

        nqquadnum=0

        !Point in -xi 1 & -xi 2
        IF(CASE.EQ.1) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(-1,0,NXQ(-2,1,nq))
            nq2=NXQ(-1,n,NXQ(-2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(1)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE IF(CASE.EQ.5) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(-2,0,NXQ(-1,1,nq))
            nq2=NXQ(-2,n,NXQ(-1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(1)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE
          NQLIST(1)=NXQ(-1,1,NXQ(-2,1,nq))
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqquadnum=nqquadnum+1
          NQGP(nqquadnum,nq)=NQLIST(1)
        ENDIF !deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(1,nj)=XQ(nj,NQLIST(1))
          ENDDO !njj2
        ENDDO !njj1

        !Average point in -xi 2
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(2,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(-2,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(2,nj)=XE(2,nj)+XQ(nj,NXQ(-2,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(-2,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(-2,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(2,nj)=XE(2,nj)/DBLE(NXQ(-2,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Point in +xi 1 & -xi 2
        IF(CASE.EQ.3) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(1,0,NXQ(-2,1,nq))
            nq2=NXQ(1,n,NXQ(-2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(3)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE IF(CASE.EQ.6) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(-2,0,NXQ(1,1,nq))
            nq2=NXQ(-2,n,NXQ(1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(3)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE
          NQLIST(3)=NXQ(1,1,NXQ(-2,1,nq))
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqquadnum=nqquadnum+1
          NQGP(nqquadnum,nq)=NQLIST(3)
        ENDIF !deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(3,nj)=XQ(nj,NQLIST(3))
          ENDDO !njj2
        ENDDO !njj1

        !Average point in -xi 1
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(4,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(-1,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(4,nj)=XE(4,nj)+XQ(nj,NXQ(-1,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(-1,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(4,nj)=XE(4,nj)/DBLE(NXQ(-1,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Centre point
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(5,nj)=XQ(nj,nq)
          ENDDO !njj2
        ENDDO !njj1
        IF(.NOT.DEFORMED) THEN
          nqquadnum=nqquadnum+1
          NQGP(nqquadnum,nq)=nq
        ENDIF !deformed

        !Average point in +xi 1
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(6,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(1,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(6,nj)=XE(6,nj)+XQ(nj,NXQ(1,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(1,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(1,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(6,nj)=XE(6,nj)/DBLE(NXQ(1,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Point in -xi1 & +xi 2
        IF(CASE.EQ.2) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(-1,0,NXQ(2,1,nq))
            nq2=NXQ(-1,n,NXQ(2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(7)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE IF(CASE.EQ.7) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(2,0,NXQ(-1,1,nq))
            nq2=NXQ(2,n,NXQ(-1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(7)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE
          NQLIST(7)=NXQ(-1,1,NXQ(2,1,nq))
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqquadnum=nqquadnum+1
          NQGP(nqquadnum,nq)=NQLIST(7)
        ENDIF !deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(7,nj)=XQ(nj,NQLIST(7))
          ENDDO !njj2
        ENDDO !njj1

        !Average point in +xi 2
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(8,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(2,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(8,nj)=XE(8,nj)+XQ(nj,NXQ(2,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(2,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(2,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(8,nj)=XE(8,nj)/DBLE(NXQ(2,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Point in +xi 1 & +xi 2
        IF(CASE.EQ.4) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(1,0,NXQ(2,1,nq))
            nq2=NXQ(1,n,NXQ(2,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(9)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE IF(CASE.EQ.8) THEN
          FOUND=.FALSE.
          DO n=1,NXQ(2,0,NXQ(1,1,nq))
            nq2=NXQ(2,n,NXQ(1,1,nq))
            DO i=1,NENQ(0,nq)
              ne1=NENQ(i,nq)
              DO j=1,NENQ(0,nq2)
                ne2=NENQ(j,nq2)
                IF(ne1.EQ.ne2) THEN
                  FOUND=.TRUE.
                  NQLIST(9)=nq2
                ENDIF
              ENDDO !j
            ENDDO !i
          ENDDO !n
          IF(.NOT.FOUND) nqHOLD=nq
        ELSE
          NQLIST(9)=NXQ(1,1,NXQ(2,1,nq))
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqquadnum=nqquadnum+1
          NQGP(nqquadnum,nq)=NQLIST(9)
        ENDIF !deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(9,nj)=XQ(nj,NQLIST(9))
          ENDDO !njj2
        ENDDO !njj1

        IF(.NOT.DEFORMED) NQGP(0,nq)=nqquadnum

      ELSE IF(NITB.EQ.3) THEN

C rgb 14/09/99
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

        nqquadnum=0

        ! Grid point 1
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq))).GT.0) THEN
              XE(1,nj)=XQ(nj,NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq))))
            ELSE
              XE(1,nj)=0.d0
            ENDIF
          ENDDO
        ENDDO

        ! Grid point 2
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.1) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-2,0,NXQ(-3,1,nq))
              nq2=NXQ(-2,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.2) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-3,0,NXQ(-2,1,nq))
              nq2=NXQ(-3,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-2,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(-3,1,nq))
          ELSEIF(NXQ(-3,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(-2,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(2,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(2,nj)=XE(2,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(2,nj)=XE(2,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2


        ! Grid point 3
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq))).GT.0) THEN
              XE(3,nj)=XQ(nj,NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq))))
            ELSE
              XE(3,nj)=0.d0
            ENDIF
          ENDDO !njj1
        ENDDO !njj2

        ! Grid point 4
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.3) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-3,0,NXQ(-1,1,nq))
              nq2=NXQ(-3,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.4) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-1,0,NXQ(-3,1,nq))
              nq2=NXQ(-1,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(-1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(-3,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(4,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(4,nj)=XE(4,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(4,nj)=XE(4,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2


        ! Grid point 5 - average point in -xi3
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(5,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(-3,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(5,nj)=XE(5,nj)+XQ(nj,NXQ(-3,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(-3,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(-3,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(5,nj)=XE(5,nj)/DBLE(NXQ(-3,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        ! Grid point 6
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.5) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-3,0,NXQ(1,1,nq))
              nq2=NXQ(-3,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.6) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(1,0,NXQ(-3,1,nq))
              nq2=NXQ(1,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(1,1,nq))
          ELSEIF(NXQ(1,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(-3,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(6,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(6,nj)=XE(6,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(6,nj)=XE(6,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        ! Grid point 7
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq))).GT.0) THEN
              XE(7,nj)=XQ(nj,NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq))))
            ELSE
              XE(7,nj)=0.d0
            ENDIF
          ENDDO !njj1
        ENDDO !njj2

        ! Grid point 8
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.7) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-3,0,NXQ(2,1,nq))
              nq2=NXQ(-3,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.8) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(2,0,NXQ(-3,1,nq))
              nq2=NXQ(2,n,NXQ(-3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-3,1,NXQ(2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-3,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(-3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(-3,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(8,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(8,nj)=XE(8,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(8,nj)=XE(8,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        ! Grid point 9
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq))).GT.0) THEN
              XE(9,nj)=XQ(nj,NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq))))
            ELSE
              XE(9,nj)=0.d0
            ENDIF
          ENDDO !njj1
        ENDDO !njj2

        ! Grid point 10
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.9) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-2,0,NXQ(-1,1,nq))
              nq2=NXQ(-2,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.10) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-1,0,NXQ(-2,1,nq))
              nq2=NXQ(-1,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-2,1,NXQ(-1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(-2,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(10,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(10,nj)=XE(10,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(10,nj)=XE(10,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 11 take average in -xi2
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(11,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(-2,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(11,nj)=XE(11,nj)+XQ(nj,NXQ(-2,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(-2,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(-2,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(11,nj)=XE(11,nj)/DBLE(NXQ(-2,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Grid point 12
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.11) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(1,0,NXQ(-2,1,nq))
              nq2=NXQ(1,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.12) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-2,0,NXQ(1,1,nq))
              nq2=NXQ(-2,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(1,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(-2,1,nq))
          ELSEIF(NXQ(-2,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(1,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(12,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(12,nj)=XE(12,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(12,nj)=XE(12,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        ! grid point 13 - average -xi1
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(13,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(-1,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(13,nj)=XE(13,nj)+XQ(nj,NXQ(-1,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(-1,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(-1,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(13,nj)=XE(13,nj)/DBLE(NXQ(-1,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Grid point 14 - centre point
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(14,nj)=XQ(nj,nq)
          ENDDO !njj2
        ENDDO !njj1
        IF(.NOT.DEFORMED) THEN
          nqquadnum=nqquadnum+1
          NQGP(nqquadnum,nq)=nq
        ENDIF !deformed

        !Grid point 15 average xi1
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(15,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(1,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(15,nj)=XE(15,nj)+XQ(nj,NXQ(1,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(1,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(1,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(15,nj)=XE(15,nj)/DBLE(NXQ(1,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Grid point 16
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.13) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-1,0,NXQ(2,1,nq))
              nq2=NXQ(-1,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.14) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(2,0,NXQ(-1,1,nq))
              nq2=NXQ(2,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(-1,1,NXQ(2,1,nq)).GT.0)THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(-1,1,nq)).GT.0)THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(-1,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(16,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(16,nj)=XE(16,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(16,nj)=XE(16,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        ! Grid point 17 average xi2
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(17,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(2,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(17,nj)=XE(17,nj)+XQ(nj,NXQ(2,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(2,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(2,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(17,nj)=XE(17,nj)/DBLE(NXQ(2,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Grid point 18
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.15) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(1,0,NXQ(2,1,nq))
              nq2=NXQ(1,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.16) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(2,0,NXQ(1,1,nq))
              nq2=NXQ(2,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(1,1,NXQ(2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(1,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(18,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(18,nj)=XE(18,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(18,nj)=XE(18,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 19
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq))).GT.0) THEN
              XE(19,nj)=XQ(nj,NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq))))
            ELSE
              XE(19,nj)=0.d0
            ENDIF
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 20
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.17) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(3,0,NXQ(-2,1,nq))
              nq2=NXQ(3,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.18) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(3,0,NXQ(-2,1,nq))
              nq2=NXQ(3,n,NXQ(-2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(-2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(-2,1,nq))
          ELSEIF(NXQ(-2,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-2,1,NXQ(3,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(20,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(20,nj)=XE(20,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(20,nj)=XE(20,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 21
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq))).GT.0) THEN
              XE(21,nj)=XQ(nj,NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq))))
            ELSE
              XE(21,nj)=0.d0
            ENDIF
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 22
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.19) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(3,0,NXQ(-1,1,nq))
              nq2=NXQ(3,n,NXQ(-1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.20) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(-1,0,NXQ(3,1,nq))
              nq2=NXQ(-1,n,NXQ(3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(-1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(-1,1,nq))
          ELSEIF(NXQ(-1,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(-1,1,NXQ(3,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(22,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(22,nj)=XE(22,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(22,nj)=XE(22,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        ! Grid point 23 average xi3
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(23,nj)=0.0d0
          ENDDO !njj2
        ENDDO !njj1
        DO n=1,NXQ(3,0,nq)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(23,nj)=XE(23,nj)+XQ(nj,NXQ(3,n,nq))
            ENDDO !njj2
          ENDDO !njj1
          IF(.NOT.DEFORMED) THEN
            nqquadnum=nqquadnum+1
            NQGP(nqquadnum,nq)=NXQ(3,n,nq)
          ENDIF !deformed
        ENDDO !n
        IF(NXQ(3,0,nq).GT.0) THEN
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(23,nj)=XE(23,nj)/DBLE(NXQ(3,0,nq))
            ENDDO !njj2
          ENDDO !njj1
        ENDIF !denom

        !Grid point 24
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.21) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(3,0,NXQ(1,1,nq))
              nq2=NXQ(3,n,NXQ(1,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.22) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(1,0,NXQ(3,1,nq))
              nq2=NXQ(1,n,NXQ(3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(1,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(1,1,nq))
          ELSEIF(NXQ(1,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(1,1,NXQ(3,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(24,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(24,nj)=XE(24,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(24,nj)=XE(24,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 25
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq))).GT.0) THEN
              XE(25,nj)=XQ(nj,NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq))))
            ELSE
              XE(25,nj)=0.d0
            ENDIF
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 26
        ACASE=.FALSE.
        NQLIST(0)=0
        DO icase=1,NCASE
          CASE=CASES(icase)
          IF(CASE.EQ.23) THEN
            ACASE=.TRUE.
            FOUND=.FALSE.
            DO n=1,NXQ(3,0,NXQ(2,1,nq))
              nq2=NXQ(3,n,NXQ(2,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ELSEIF(CASE.EQ.24) THEN
            FOUND=.FALSE.
            ACASE=.TRUE.
            DO n=1,NXQ(2,0,NXQ(3,1,nq))
              nq2=NXQ(2,n,NXQ(3,1,nq))
              DO i=1,NENQ(0,nq)
                ne1=NENQ(i,nq)
                DO j=1,NENQ(0,nq2)
                  ne2=NENQ(j,nq2)
                  IF(ne1.EQ.ne2) THEN
                    FOUND=.TRUE.
                    NQLIST(0)=NQLIST(0)+1
                    NQLIST(NQLIST(0))=nq2
                    IF(NQLIST(0).EQ.2) THEN
                      IF(NQLIST(1).EQ.NQLIST(2)) NQLIST(0)=NQLIST(0)-1
                    ENDIF
                  ENDIF
                ENDDO !j
              ENDDO !i
            ENDDO !n
            IF(.NOT.FOUND) nqHOLD=nq
          ENDIF
        ENDDO
        IF(.NOT.ACASE) THEN
          IF(NXQ(3,1,NXQ(2,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(3,1,NXQ(2,1,nq))
          ELSEIF(NXQ(2,1,NXQ(3,1,nq)).GT.0) THEN
            NQLIST(0)=NQLIST(0)+1
            NQLIST(NQLIST(0))=NXQ(2,1,NXQ(3,1,nq))
          ELSE
            CALL ASSERT(.FALSE.,'>>Unable to find adjoining grid'//
     '        ' point',ERROR,*9999)
          ENDIF
        ENDIF
        IF(.NOT.DEFORMED) THEN
          nqgt=nqgt+NQLIST(0)-1
          IF(nqgt.LE.NQGM) THEN
            DO nonq=1,NQLIST(0)
              nqquadnum=nqquadnum+1
              NQGP(nqquadnum,nq)=NQLIST(nonq)
            ENDDO
          ENDIF
        ENDIF !.not.deformed
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(26,nj)=0.d0
          ENDDO !njj1
        ENDDO !njj2
        DO nonq=1,NQLIST(0)
          DO njj1=1,3,NJJSTEP
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              XE(26,nj)=XE(26,nj)+XQ(nj,NQLIST(nonq))
            ENDDO !njj1
          ENDDO !njj2
        ENDDO !nonq
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            XE(26,nj)=XE(26,nj)/DBLE(NQLIST(0))
          ENDDO !njj1
        ENDDO !njj2

        !Grid point 27
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NXQ(3,1,NXQ(2,1,NXQ(1,1,nq))).GT.0) THEN
              XE(27,nj)=XQ(nj,NXQ(3,1,NXQ(2,1,NXQ(1,1,nq))))
            ELSE
              XE(27,nj)=0.d0
            ENDIF
          ENDDO !njj1
        ENDDO !njj2

        IF(.NOT.DEFORMED) NQGP(0,nq)=nqquadnum

C        NQLIST(0)=27
C        NQLIST(1)=NXQ(-3,1,NXQ(-2,1,NXQ(-1,1,nq)))
C        NQLIST(2)=NXQ(-3,1,NXQ(-2,1,nq))
C        NQLIST(3)=NXQ(-3,1,NXQ(-2,1,NXQ(1,1,nq)))
C        NQLIST(4)=NXQ(-3,1,NXQ(-1,1,nq))
C        NQLIST(5)=NXQ(-3,1,nq)
C        NQLIST(6)=NXQ(-3,1,NXQ(1,1,nq))
C        NQLIST(7)=NXQ(-3,1,NXQ(2,1,NXQ(-1,1,nq)))
C        NQLIST(8)=NXQ(-3,1,NXQ(2,1,nq))
C        NQLIST(9)=NXQ(-3,1,NXQ(2,1,NXQ(1,1,nq)))
C        NQLIST(10)=NXQ(-2,1,NXQ(-1,1,nq))
C        NQLIST(11)=NXQ(-2,1,nq)
C        NQLIST(12)=NXQ(-2,1,NXQ(1,1,nq))
C        NQLIST(13)=NXQ(-1,1,nq)
C        NQLIST(14)=nq
C        NQLIST(15)=NXQ(1,1,nq)
C        NQLIST(16)=NXQ(2,1,NXQ(-1,1,nq))
C        NQLIST(17)=NXQ(2,1,nq)
C        NQLIST(18)=NXQ(2,1,NXQ(1,1,nq))
C        NQLIST(19)=NXQ(3,1,NXQ(-2,1,NXQ(-1,1,nq)))
C        NQLIST(20)=NXQ(3,1,NXQ(-2,1,nq))
C        NQLIST(21)=NXQ(3,1,NXQ(-2,1,NXQ(1,1,nq)))
C        NQLIST(22)=NXQ(3,1,NXQ(-1,1,nq))
C        NQLIST(23)=NXQ(3,1,nq)
C        NQLIST(24)=NXQ(3,1,NXQ(1,1,nq))
C        NQLIST(25)=NXQ(3,1,NXQ(2,1,NXQ(-1,1,nq)))
C        NQLIST(26)=NXQ(3,1,NXQ(2,1,nq))
C        NQLIST(27)=NXQ(3,1,NXQ(2,1,NXQ(1,1,nq)))
C
C        IF(.NOT.DEFORMED) THEN
C          NQGP(0,nq)=19
C          NQGP(1,nq)=NQLIST(2)
C          NQGP(2,nq)=NQLIST(4)
C          NQGP(3,nq)=NQLIST(5)
C          NQGP(4,nq)=NQLIST(6)
C          NQGP(5,nq)=NQLIST(8)
C          NQGP(6,nq)=NQLIST(10)
C          NQGP(7,nq)=NQLIST(11)
C          NQGP(8,nq)=NQLIST(12)
C          NQGP(9,nq)=NQLIST(13)
C          NQGP(10,nq)=NQLIST(14)
C          NQGP(11,nq)=NQLIST(15)
C          NQGP(12,nq)=NQLIST(16)
C          NQGP(13,nq)=NQLIST(17)
C          NQGP(14,nq)=NQLIST(18)
C          NQGP(15,nq)=NQLIST(20)
C          NQGP(16,nq)=NQLIST(22)
C          NQGP(17,nq)=NQLIST(23)
C          NQGP(18,nq)=NQLIST(24)
C          NQGP(19,nq)=NQLIST(26)
C        ENDIF !deformed
C
C        DO njj1=1,3,NJJSTEP
C          DO njj2=1,NJ_LOC(njj1,0,nr)
C            nj=NJ_LOC(njj1,njj2,nr)
C            DO nqq=1,NQLIST(0)
C              nq1=NQLIST(nqq)
C              IF(nq1.GT.0) THEN
C                XE(nqq,nj)=XQ(nj,nq1)
C              ELSE
C                XE(nqq,nj)=0.0d0
C              ENDIF
C            ENDDO !nqq
C          ENDDO !njj2
C        ENDDO !njj1

        IF(.NOT.DEFORMED) THEN
          IF(nqgt.GT.NQGM) THEN
            WRITE(ERROR,'(''NQGM needs to be at least2 '',I12)') nqgt
C            WRITE(*,*) '2 nqgt = ',nqgt
C            IEND=0
C            CALL APPENDC(IEND,'NQGM needs to be at least ',ERROR)
C            CALL APPENDI(IEND,nqgt,ERROR)
            GOTO 9999
          ENDIF ! nqquadnum > NQGM
        ENDIF

      ENDIF !nit(nb)


      IF(nqHOLD.NE.0) THEN
        CALL ASSERT(.FALSE.,' >>A cusp point connection was not found',
     '    ERROR,*9999)
      ENDIF

CMLB - for reference
C      IF(NIT(NBJ(nj)).EQ.1) THEN
C        NQLIST(0)=3
C        NQLIST(1)=NXQ(-1,1,nq)
C        NQLIST(2)=nq
C        NQLIST(3)=NXQ(1,1,nq)
C      ELSE IF(NIT(NBJ(nj)).EQ.2) THEN
C        NQLIST(0)=9
C        NQLIST(1)=NXQ(-2,1,NXQ(-1,1,nq))
C        NQLIST(2)=NXQ(-2,1,nq)
C        NQLIST(3)=NXQ(-2,1,NXQ(1,1,nq))
C        NQLIST(4)=NXQ(-1,1,nq)
C        NQLIST(5)=nq
C        NQLIST(6)=NXQ(1,1,nq)
C        NQLIST(7)=NXQ(2,1,NXQ(-1,1,nq))
C        NQLIST(8)=NXQ(2,1,nq)
C        NQLIST(9)=NXQ(2,1,NXQ(1,1,nq))
C      ELSE IF(NIT(NBJ(nj)).EQ.3) THEN
C        !unchanged
C      ENDIF !nit(nb)

      CALL EXITS('XQXE')
      RETURN
 9999 CALL ERRORS('XQXE',ERROR)
      CALL EXITS('XQXE')
      RETURN 1
      END


