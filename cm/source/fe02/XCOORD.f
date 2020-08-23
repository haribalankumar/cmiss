      SUBROUTINE XCOORD(IBT,IDO,INP,NBJ,ne,NEELEM,
     '  NKJE,NPF,NPNE,NPNODE,nr,NVJE,
     '  SE,TOL_XCOORD,XA,XE,XI,XP,XP1,ERROR,*)

C#### Subroutine: XCOORD
C###  Description:
C###    XCOORD evaluates the Xi-coords (XI(ni)) and element number
C###    (ne)  of a specified Xj-coord point (XP1(nj)).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),ne,NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),TOL_XCOORD,XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XI(3),XP(NKM,NVM,NJM,NPM),XP1(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ICHECK,ie,IFLAG(3),nb,nbb,nee,NEP(8),ni,nj,NMIN,nn,
     '  noelem,nonode,np
      REAL*8 D,D2,DMIN,D2MIN,DELTA(3),DIST,PXI,XN_LOCAL(3),XS(3),YN(3),
     '  YP1(3),YS(3)

      CALL ENTERS('XCOORD',*9999)

C GBS 20/10/94
C!!!  Fix this (maybe pass NRE to this subroutine)
C   nr=1
C   Fixed by RGB 7/4/98 by passing in nr from OPC30,OPC40
      DMIN=1.0D6
      D2MIN=DMIN
      CALL COORD(ITYP10(1),1,XP1,YP1,ERROR,*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' Point XP1 has curvilinear coords:'','
     '    //'3E11.4,'' & r.c. coords:'',3E11.4)')
     '    (XP1(I),I=1,3),(YP1(I),I=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DO nonode=1,NPNODE(0,1)
        np=NPNODE(nonode,1)
        DO nj=1,NJT
          XN_LOCAL(nj)=XP(1,1,nj,np)
        ENDDO
        CALL COORD(ITYP10(1),1,XN_LOCAL,YN,ERROR,*9999)
        D=0.0D0
        DO nj=1,NJT
          D=D+(YP1(nj)-YN(nj))**2
        ENDDO
        IF(ITYP10(1).EQ.1) THEN      !rectangular cartesian coords
          IF(D.LT.DMIN) THEN
            DMIN=D
            NMIN=np
          ENDIF
        ELSE IF(ITYP10(1).GT.1) THEN !curvilinear coords
!         Checking for coincident nodes,i.e choose the correct sector.
          D2= DABS(XP1(2)-XN_LOCAL(2))
          IF( DABS(XN_LOCAL(3)-0.50D0*PI).LT.TOL_XCOORD) THEN
            IF(D.LT.DMIN.OR. DABS(D-DMIN).LT.TOL_XCOORD) THEN
              IF(D2.LT.D2MIN) THEN
                DMIN=D
                D2MIN=D2
                NMIN=np
              ENDIF
            ENDIF
          ELSE
            IF(D.LT.DMIN) THEN
              DMIN=D
              NMIN=np
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' The closest node to XP1 is # '',I4,'
     '    //''' at distance '',E11.4)') NMIN,DMIN
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF(ITYP10(1).GT.1) THEN !curvilinear coords
          WRITE(OP_STRING,'('' and the theta diff is '',E11.4)')
     '      D2MIN
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
CC$      call mp_unsetlock()
      ENDIF
      ie=0
      DO 40 noelem=1,NEELEM(0,1)
        ne=NEELEM(noelem,1)
        nb=NBJ(1,ne)
        DO 30 nn=1,NNT(nb)
          np=NPNE(nn,nb,ne)
          IF(np.NE.NMIN) GO TO 30
            IF(DMIN.LT.1.0D-6.AND. DABS(D2MIN).LT.TOL_XCOORD) THEN
              XI(1)=DBLE(nn)/2.0D0-DBLE(nn-1)/2.0D0
              XI(2)=DBLE(nn+1)/2.0D0-1.0D0
              GO TO 9998
            ENDIF
          ie=ie+1
          NEP(ie)=ne
          GO TO 40
 30     CONTINUE
 40   CONTINUE
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' Surrounding elements are:'',10I3)')
     '    (NEP(ne),ne=1,ie)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(ITYP10(1).GT.1.AND.DMIN.LT.1.D-06) THEN
C       assumed for central ar apex concindent nodes...
C       may need to search element node nn and use that to find
C       XI(2)=f(nn) ??
        XI(2)=1
      DO 480 nee=1,ie
        ne=NEP(nee)
        nb=NBJ(1,ne)
        NJT=NJ_LOC(NJL_GEOM,0,nr)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '    XA(1,1,ne),XE,XP,ERROR,*9999)
        XI(1)=0.0D0
C       DELTA(1)=0.10D0
        DELTA(1)=0.50D0
 445    CONTINUE
        nbb=NBJ(2,ne)
        XS(2)=PXI(IBT(1,1,nbb),IDO(1,1,0,nbb),INP(1,1,nbb),nbb,1,XI,
     '    XE(1,2))
        DIST= DABS(XP1(2)-XS(2))
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '(/'' El '',I2,'' Xi-coords: '','
     '      //'3F10.7,/'' Initial dist (theta) to XP1='',E11.4)')
     '      ne,(XI(I),I=1,3),DSQRT(DIST)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        XI(1)=XI(1)+DELTA(1)
        IF(XI(1).LE.1.0D0) GO TO 448
        XI(1)=XI(1)-2.0D0*DELTA(1)
        DELTA(1)=0.50D0*DELTA(1)
        GO TO 445
 448    CONTINUE
        nbb=NBJ(2,ne)
        XS(2)=PXI(IBT(1,1,nbb),IDO(1,1,0,nbb),INP(1,1,nbb),nbb,1,XI,
     '    XE(1,2))
        D= DABS(XP1(2)-XS(2))
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(28X,3F10.7,/'' Current dist to XP1='','
     '      //'E11.4)') (XI(I),I=1,3),DSQRT(D)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(D.LT.TOL_XCOORD) GOTO 85
        IF(D.LT.DIST) GOTO 445
        XI(1)=XI(1)-2.0D0*DELTA(1)
        IF(XI(1).LT.(-0.01d0*TOL_XCOORD)) XI(1)=XI(1)+DELTA(1)
        DELTA(1)=0.50D0*DELTA(1)
        IF(DELTA(1).LT.(0.1d0*TOL_XCOORD)) GOTO 480
        GOTO 445
 480  CONTINUE
      ENDIF

      DO 80 nee=1,ie
        ne=NEP(nee)
        nb=NBJ(1,ne)
        NJT=NJ_LOC(NJL_GEOM,0,nr)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '    XA(1,1,ne),XE,XP,ERROR,*9999)
        DO ni=1,NIT(nb)
          XI(ni)=0.0D0
C         DELTA(ni)=0.10D0
          DELTA(ni)=0.50D0
        ENDDO
 45     CONTINUE
        DO 60 ni=1,NIT(nb)
          IFLAG(ni)=0
          DO nj=1,NJT
            nbb=NBJ(nj,ne)
            XS(nj)=PXI(IBT(1,1,nbb),IDO(1,1,0,nbb),INP(1,1,nbb),nbb,1,
     '        XI,XE(1,nj))
          ENDDO
          CALL COORD(ITYP10(1),1,XS,YS,ERROR,*9999)
          DIST=0.0D0
          DO nj=1,NJT
            DIST=DIST+(YP1(nj)-YS(nj))**2
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' El '',I2,'' dir '',I1,'
     '        //''' Xi-coords: '',3F10.7,/'' Initial dist to XP1='','
     '        //'E11.4)') ne,ni,(XI(I),I=1,3),DSQRT(DIST)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
C 470      FORMAT(/' El ',I2,' dir ',I1,' Xi-coords: ',
C     '             3F10.7,/' Initial dist to XP1=',E11.4)
          XI(ni)=XI(ni)+DELTA(ni)
          IF(XI(ni).LE.1.0D0) GO TO 48
          XI(ni)=XI(ni)-2.0D0*DELTA(ni)
          DELTA(ni)=0.50D0*DELTA(ni)
          GO TO 60
 48       DO nj=1,NJT
            nbb=NBJ(nj,ne)
            XS(nj)=PXI(IBT(1,1,nbb),IDO(1,1,0,nbb),INP(1,1,nbb),nbb,1,
     '        XI,XE(1,nj))
          ENDDO
          CALL COORD(ITYP10(1),1,XS,YS,ERROR,*9999)
          D=0.0D0
          DO nj=1,NJT
            D=D+(YP1(nj)-YS(nj))**2
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(28X,3F10.7,/'' Current dist to XP1='','
     '        //'E11.4,'' IFLAG= '',3(I2,1X))')
     '        (XI(I),I=1,3),DSQRT(D),(IFLAG(I),I=1,NIT(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
C 520      FORMAT(28X,3F10.7,/' Current dist to XP1=',E11.4,' IFLAG= ',
C     '           3(I2,1X))
          IF(D.LT.TOL_XCOORD) GOTO 85
          IF(D.LT.DIST) GOTO 60
          XI(ni)=XI(ni)-2.0D0*DELTA(ni)
          IF(XI(ni).LT.-1.0D-8) XI(ni)=XI(ni)+DELTA(ni)
          DELTA(ni)=0.50D0*DELTA(ni)
          IF(IFLAG(ni).EQ.1) GOTO 60
C         IF(DELTA(ni).LT.(100.0d0*TOL_XCOORD)) IFLAG(ni)=1
          IF(DELTA(ni).LT.(0.01d0*TOL_XCOORD)) IFLAG(ni)=1
 60     CONTINUE
        ICHECK=0
        DO ni=1,NIT(nb)
           ICHECK=ICHECK+IFLAG(ni)
        ENDDO
C       If Icheck=nit(nb), then all directions have failed..
        IF(ICHECK.EQ.NIT(nb)) GOTO 80
        GOTO 45
 80   CONTINUE
      ne=0
 85   CONTINUE

 9998 CALL EXITS('XCOORD')
      RETURN
 9999 CALL ERRORS('XCOORD',ERROR)
      CALL EXITS('XCOORD')
      RETURN 1
      END

C CPB 25/10/94 This is not really used anymore ??

