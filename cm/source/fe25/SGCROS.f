      SUBROUTINE SGCROS(INDEX,IBT,IDO,INP,ISCROS,ISEG,iw,
     '  NAN,NBJ,NEELEM,NENP,NHE,nj,NKJE,NNB,NOCROS,NPF,NPNE,NRE,
     '  NVJE,NXI,nx,CSEG,DATA,
     '  PG,SE,SHEET_TYPE,XA,XE,XG,XP,ZDD,ZE,ZG,ZVAL,ERROR,*)

C#### Subroutine: SGCROS
C###  Description:
C###    SGCROS creates cross-section segment ISCROS(iw,nocros).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISCROS,ISEG(*),iw,NAN(NIM,NAM,NBFM),
     '  NBJ(NJM,NEM),NENP(NPM,0:NEPM,0:NRM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM),nj,NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),NOCROS,
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM),nx
      REAL*8 PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  ZDD(NJM),ZE(NSM,NHM),ZG(NHM,NUM),ZVAL
      CHARACTER CSEG(*)*(*),ERROR*(*),SHEET_TYPE*(*)
      LOGICAL DATA
!     Local Variables
      INTEGER ICOORD,INDEX_OLD,ixi,nb,ne,noelem,NICONT,nr,NROOTS
      REAL*8 ROOTS(2,12),TEXT_POINT(3),X(3,3),XDIFF,
     '  XICONT,XLM(3),XXDIF,XXMAX,XXMIN,YYDIF,YYMAX,YYMIN,ZZDIF,
     '  ZZMAX,ZZMIN
      CHARACTER STRG*17

      CALL ENTERS('SGCROS',*9999)
      CALL OPEN_SEGMENT(ISCROS,ISEG,iw,'CROS',INDEX,INDEX_OLD,
     '  NOCROS,1,CSEG,ERROR,*9999)

      IF(iw.LE.3) THEN
        XXMIN=DBLE(XMIN)+0.01D0*DBLE(XMAX-XMIN)
        XXMAX=DBLE(XMAX)-0.01D0*DBLE(XMAX-XMIN)
        YYMIN=DBLE(YMIN)+0.01D0*DBLE(YMAX-YMIN)
        YYMAX=DBLE(YMAX)-0.01D0*DBLE(YMAX-YMIN)
        ZZMIN=DBLE(ZMIN)+0.01D0*DBLE(ZMAX-ZMIN)
        ZZMAX=DBLE(ZMAX)-0.01D0*DBLE(ZMAX-ZMIN)
        XXDIF=XXMAX-XXMIN
        YYDIF=YYMAX-YYMIN
        ZZDIF=ZZMAX-ZZMIN
        TEXT_POINT(1)=XXMIN+0.01D0*XXDIF
        TEXT_POINT(2)=YYMIN+0.01D0*YYDIF
        TEXT_POINT(3)=ZZMIN+0.01D0*ZZDIF
      ELSE IF(iw.EQ.13) THEN !sheet angle plot
        XDIFF=DBLE(XMAX-XMIN)
        TEXT_POINT(1)=DBLE(XMAX)-0.05D0*XDIFF !1st coord is axial
        TEXT_POINT(2)=0.0D0             !2nd coord is radial
      ENDIF
      IF(iw.EQ.2) THEN       !nj=1
        ICOORD=1
C        STRG='X='//CFROMR(ZVAL,'(E10.3)')
        WRITE(STRG,'(''X='',E10.3)') ZVAL
      ELSE IF(iw.EQ.1) THEN  !nj=2
        ICOORD=1
C        STRG='Y='//CFROMR(ZVAL,'(E10.3)')
        WRITE(STRG,'(''Y='',E10.3)') ZVAL
      ELSE IF(iw.EQ.3) THEN  !nj=3
        ICOORD=1
C        STRG='Z='//CFROMR(ZVAL,'(E10.3)')
        WRITE(STRG,'(''Z='',E10.3)') ZVAL
      ELSE IF(iw.EQ.13) THEN !nj=2 or 3
        WRITE(OP_STRING,*) ' zval=',zval
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ICOORD=ITYP10(1)
C        STRG='Theta='//CFROMR(ZVAL*180.0D0/PI,'(F6.1)')
        WRITE(STRG,'(''Theta='',F6.1)') ZVAL*180.0D0/PI
      ENDIF
      CALL TEXT(1,iw,STRG,TEXT_POINT,ERROR,*9999)

      IF(DATA.AND.(iw.EQ.1.OR.iw.EQ.2)) THEN !iw=1,2
        X(1,1)=XXMIN
        X(2,1)=YYMIN
        X(3,1)=ZDD(3)
        X(1,2)=XXMAX
        X(2,2)=YYMAX
        X(3,2)=ZDD(3)
        CALL POLYLINE(INDEX,iw,2,X,ERROR,*9999)
        X(1,1)=ZDD(1)
        X(2,1)=ZDD(2)
        X(3,1)=ZZMIN
        X(1,2)=ZDD(1)
        X(2,2)=ZDD(2)
        X(3,2)=ZZMAX
        CALL POLYLINE(INDEX,iw,2,X,ERROR,*9999)
      ELSE IF(DATA.AND.iw.EQ.3) THEN  !iw=3
        X(1,1)=XXMIN
        X(2,1)=ZDD(2)
        X(1,2)=XXMAX
        X(2,2)=ZDD(2)
        CALL POLYLINE(INDEX,iw,2,X,ERROR,*9999)
        X(1,1)=ZDD(1)
        X(2,1)=YYMIN
        X(1,2)=ZDD(1)
        X(2,2)=YYMAX
        CALL POLYLINE(INDEX,iw,2,X,ERROR,*9999)
      ELSE IF(iw.EQ.13) THEN  !sheet plot axis
        XDIFF=DBLE(XMAX)-DBLE(XMIN)
        X(1,1)=DBLE(XMIN)+0.1D0*XDIFF !1st coord is axial
        X(2,1)=0.0D0            !2nd coord is radial
        X(1,2)=DBLE(XMAX)-0.1D0*XDIFF
        X(2,2)=0.0D0
        CALL POLYLINE(INDEX,iw,2,X,ERROR,*9999)
        X(1,1)=0.0D0            !1st coord is axial
        X(2,1)=0.0D0            !2nd coord is radial
        X(1,2)=0.0D0
        X(2,2)=0.7D0*XDIFF
        CALL POLYLINE(INDEX,iw,2,X,ERROR,*9999)
      ENDIF

      IF(iw.LE.3.OR.(iw.EQ.13.AND.SHEET_TYPE(1:6).EQ.'RADIAL')) THEN
        NICONT=3
        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
          nr=NRE(ne)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          nb=NBJ(nj,ne)
          DO ixi=0,1
            XICONT=DBLE(ixi)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' nj='',I1,'' nbh='',I1,'' iw='',I2,'
     '          //''' NICONT='',I1,'' ne='',I3,'' XICONT='',F3.1)')
     '          nj,nb,iw,NICONT,ne,XICONT
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
!           Note: nh & nbh replaced with nj & nbj 1-feb-1990
            CALL CROOTS('GEOMETRIC',IBT,IDO,INP,100,NAN,NBJ(1,ne),
     '        NBJ(1,ne),nj,NHE(ne),nj,ICOORD,NICONT,
     '        nr,NROOTS,nx,ROOTS,0.02D0,XICONT,PG,XE,XG,
     '        ZE,ZG,ZVAL,ERROR,*9999)
            CALL CONTOR(INDEX,'GEOMETRIC',IBT,IDO,
     '        INP,ICOORD,100,iw,NICONT,NAN,NBJ(1,ne),NBJ(1,ne),nj,
     '        NHE(ne),nj,nr,NROOTS,nx,
     '        ROOTS,1.0D-2,PG,XE,XG,XICONT,5.0D-2,XLM,
     '        ZE,ZG,ZVAL,ERROR,*9999)
          ENDDO
        ENDDO

      ELSE IF(iw.EQ.13.AND.SHEET_TYPE(1:6).EQ.'NORMAL') THEN
        CALL CROSSSECTION(INDEX,IBT,IDO,INP,iw,NBJ,NEELEM,NENP,nj,NKJE,
     '    NNB,NPF,NPNE,NRE,NVJE,NXI,SE,XA,XE,XP,ZVAL,ERROR,*9999)
      ENDIF

      CALL CLOSE_SEGMENT(ISCROS,iw,ERROR,*9999)

      CALL EXITS('SGCROS')
      RETURN
 9999 CALL ERRORS('SGCROS',ERROR)
      CALL EXITS('SGCROS')
      RETURN 1
      END


