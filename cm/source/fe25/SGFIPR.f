      SUBROUTINE SGFIPR(INDEX,IBT,IDO,INP,ISEG,ISFIPR,iw,ixi,NAN,NBJ,
     '  NELIST,NKJE,NPF,NPNE,NRE,NVJE,SE,XA,XE,XG,XIPOS,XP,YMAGN,
     '  CSEG,ERROR,*)

C#### Subroutine: SGFIPR
C###  Description:
C###    SGFIPR creates fibre angle profile segments ISFIPR. iw is 33.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISEG(*),ISFIPR,iw,ixi,NAN(NIM,NAM,NBFM),
     '  NBJ(NJM,NEM),NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XIPOS(*),XP(NKM,NVM,NJM,NPM),YMAGN
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,INDEX_OLD,
     '  N1XI,nb,NITB,ne,ni,nolist,nox,noxi,nr,
     '  NTX,NTXI,NTXMAX
      PARAMETER(NTXMAX=81)
      REAL*8 ETA,PTS(3,NTXMAX),X(2),XPROF(NTXMAX),
     '  YPROF(NTXMAX),YPROFM,Z(3),Z0(3)
      CHARACTER CHAR*7,CHAR1*1,TITLE*60

      CALL ENTERS('SGFIPR',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,*) 'NELIST(0)=',NELIST(0)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(CHAR1,'(I1)') ixi
      TITLE='Fibre Angle Profiles vs Xi_'//CHAR1
      NTX=NTXMAX
      NTXI=(NTX-1)/NELIST(0)
      NTX=NTXI*NELIST(0)+1
      nox=0
      XPROF(1)=0.0D0
      DO nolist=1,NELIST(0)
        ne=NELIST(nolist)
        nr=NRE(ne)
        IF(DOP) THEN
          WRITE(OP_STRING,*) 'ne =',ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        nb=NBJ(1,ne)
        IF(DOP) THEN
          WRITE(OP_STRING,*) 'nb =',nb
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        NITB=NIT(nb)
        IF(DOP) THEN
          WRITE(OP_STRING,*) 'NIT=',NITB
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        IF(nolist.EQ.1) THEN
          N1XI=0
        ELSE IF(nolist.GT.1) THEN
          N1XI=1
        ENDIF
        DO noxi=N1XI,NTXI
          nox=nox+1
          XIPOS(ixi)=DBLE(noxi)/DBLE(NTXI)
          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,XG,XIPOS,
     '      ERROR,*9999)
          IF(NELIST(0).EQ.1) THEN
            XPROF(nox)=XIPOS(ixi)
          ELSE IF(NELIST(0).GT.1) THEN
            IF(nox.EQ.1) THEN
              CALL XZ(ITYP10(nr),XG(1,1),Z0)
            ELSE IF(nox.GT.1) THEN
              CALL XZ(ITYP10(nr),XG(1,1),Z)
              XPROF(nox)=XPROF(nox-1)+DSQRT((Z(1)-Z0(1))**2
     '          +(Z(2)-Z0(2))**2+(Z(3)-Z0(3))**2)
              DO ni=1,NITB
                Z0(ni)=Z(ni)
              ENDDO
            ENDIF
          ENDIF
          nb=NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
          IF(DOP) THEN
            WRITE(OP_STRING,*)
     '        'nb=',nb,' NNT(nb)=',NNT(nb),'NAT(nb)=',NAT(nb)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          NR=NR
          IF(NNT(nb).GT.0) THEN
            IF(NJ_LOC(NJL_FIBR,0,nr).EQ.1) THEN
              ETA=XG(NJ_LOC(NJL_FIBR,1,nr),1)
            ELSE IF(NJ_LOC(NJL_FIBR,0,nr).EQ.2) THEN
            ENDIF
          ELSE IF(NAT(nb).GT.0) THEN
             ETA=XA(1,NJ_LOC(NJL_FIBR,1,nr),ne)
          ENDIF
          YPROF(nox)=ETA*180.0D0/PI
          IF(DOP) THEN
            WRITE(OP_STRING,*) 'X,Y=',XPROF(nox),YPROF(nox)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
      ENDDO

      IF(YMAGN.EQ.0.0D0) THEN
        YPROFM=0.0D0
        DO nox=1,NTX
          IF(DABS(YPROF(nox)).GT.YPROFM) THEN
            YPROFM=DABS(YPROF(nox))
          ENDIF
        ENDDO
      ELSE
        YPROFM=YMAGN
      ENDIF
      DO nox=1,NTX
        XPROF(nox)=XPROF(nox)/XPROF(NTX)
        YPROF(nox)=YPROF(nox)/YPROFM
      ENDDO

      CALL OPEN_SEGMENT(ISFIPR,ISEG,iw,'FIPR',INDEX,INDEX_OLD,
     '  iw,1,CSEG,ERROR,*9999)

      PTS(1,1)=0.0D0
      PTS(2,1)=0.0D0
      PTS(1,2)=1.0D0
      PTS(2,2)=0.0D0
      CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
      PTS(1,1)=0.0D0
      PTS(2,1)=-1.0D0
      PTS(1,2)=0.0D0
      PTS(2,2)=1.0D0
      CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
      PTS(2,1)=-0.95D0
      PTS(2,2)= 0.95D0
      DO nolist=1,NELIST(0)-1
        nox=nolist*(NTXI+1)+1
        PTS(1,1)=XPROF(nox)
        PTS(1,2)=XPROF(nox)
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
      ENDDO
      PTS(2,1)=0.0D0
      PTS(2,2)=0.03D0
      DO I=1,2*NELIST(0)
        PTS(1,1)=DBLE(I)/DBLE(2*NELIST(0))
        PTS(1,2)=PTS(1,1)
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
        IF(MOD(I,2).EQ.0) THEN
          WRITE(CHAR,'(F4.2)') X(1)
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
        ENDIF
      ENDDO
      PTS(1,2)=0.02D0
      DO i=1,5
        PTS(1,1)=0.0D0
        PTS(2,1)=DBLE(i-1)/2.0D0-1.0D0
        PTS(2,2)=PTS(2,1)
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
        WRITE(CHAR,'(F7.2)') DBLE(i)*YPROFM
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        PTS(1,1)=-0.01D0
        CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
      ENDDO
      PTS(1,1)=0.5D0
      PTS(2,1)=-1.0D0
      CALL STRING_TRIM(TITLE,IBEG,IEND)
      CALL TEXT(1,iw,TITLE(IBEG:IEND),PTS(1,1),ERROR,*9999)
      DO nox=1,NTX
        PTS(1,nox)=XPROF(nox)
        PTS(2,nox)=YPROF(nox)
      ENDDO
      CALL POLYLINE(INDEX,iw,NTX,PTS,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISFIPR,iw,ERROR,*9999)

      CALL EXITS('SGFIPR')
      RETURN
 9999 CALL ERRORS('SGFIPR',ERROR)
      CALL EXITS('SGFIPR')
      RETURN 1
      END


