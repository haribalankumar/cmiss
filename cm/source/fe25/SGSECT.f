      SUBROUTINE SGSECT(INDEX,IBT,IDO,INP,ISEG,ISSECT,iw,NBH,NBJ,
     '  NEELEM,NELEM,NHE,NHP,NKJE,NKH,NPF,NPNE,NPLIST,
     '  NPNODE,NRE,NTELEM,NTSECT,NVHP,NVJE,
     '  NYNE,NYNP,CSEG,SE,TYPE,XA,XD,XE,XI,XP,YP,ZA,ZD,ZP,ERROR,*)

C#### Subroutine: SGSECT
C###  Description:
C###    SGSECT creates segment ISSECT.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'press00.cmn'
      INCLUDE 'sect00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISEG(*),ISSECT,iw,NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELEM(*),NHE(NEM),
     '  NHP(NPM),NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPLIST(0:NPM),NRE(NEM),NTELEM,NTSECT,
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XD(*),XE(NSM,NJM),XI(*),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER I,IBEG,IEND,INDEX_OLD,
     '  N,N1CHORD,NB1,NB4,nc,nd,ne,ne1,ne2,ng,
     '  ng1,ng2,ni,nj,nochord,noelem,NOFIBR,nohist,nolist,nopts,np,nr,
     '  nrc,NSFIBR,NTCHORD,NTPTS,nv,nx,ny
      REAL*8 DLAMDA,ETA,FIBR(3,200),OFFSET,PTS(3,200),PTS2(3,2),PXI,
     '  RATIO,RLAMDA(21),SUM,T,X(3),X1,X2,XXE(4,2),YMAX
      CHARACTER CHAR*7

      CALL ENTERS('SGSECT',*9999)
      nc=1 !Temporary MPN 12-Nov-94
      nv=1 !Temporary
      CALL OPEN_SEGMENT(ISSECT,ISEG,iw,'SECT',INDEX,INDEX_OLD,
     '  NTSECT,1,CSEG,ERROR,*9999)

      IF(TYPE(1:8).EQ.'PRESSURE') THEN
        IF(NTSECT.EQ.1) THEN
          PTS(1,1)= 0.0D0
          PTS(2,1)=-0.5D0
          PTS(1,2)= 1.0D0
          PTS(2,2)= PTS(2,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          PTS(1,1)= 0.0D0
          PTS(2,1)=-1.0D0
          PTS(1,2)= 0.0D0
          PTS(2,2)= 1.0D0
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          PTS(1,2)=0.02D0
          DO I=1,5
            PTS(1,1)=0.0D0
            PTS(2,1)=DBLE(I-1)/2.0D0-1.0D0
            PTS(2,2)=PTS(2,1)
            CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
            WRITE(CHAR,'(F7.3)') (PTS(2,1)+0.5D0)*PRESS_MAX/FACTOR
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            PTS(1,1)=-0.01D0
            CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
          ENDDO
          DO I=1,2
            PTS(1,1)= 0.5D0*DBLE(I)
            PTS(2,1)=-0.5D0
            PTS(1,2)= PTS(1,1)
            PTS(2,2)= PTS(2,1)+0.03D0
            CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          ENDDO
          PTS(1,1)=-0.06D0
          PTS(2,1)= 0.8D0
          CALL TEXT(1,iw,'kPa',PTS(1,1),ERROR,*9999)
          PTS(1,1)=-0.5D0
          PTS(2,1)=-1.0D0
          CALL TEXT(1,iw,'Chord pressure distribution',PTS(1,1),
     '      ERROR,*9999)
        ENDIF
        NTCHORD=NPLIST(0)
        DO nochord=1,NTCHORD
          N1CHORD=NPLIST(nochord)
          ne1=2*INT(DBLE(N1CHORD-1)/3.d0+0.001d0)+1
          ne2=ne1+1
          ng1=3*(MOD(N1CHORD-1,3)+1)-2
          ng2=ng1+2
          IF(NTCHORD.GT.1) THEN
            RATIO=DBLE(nochord-1)/DBLE(NTCHORD-1)
          ELSE IF(NTCHORD.EQ.1) THEN
            RATIO=0.d0
          ENDIF
          OFFSET=-0.5d0+RATIO*0.5d0
          X1=0.5d0*RATIO
          X2=1.0d0-0.25d0*RATIO
          nopts=0
          DO ne=ne1,ne2
            DO ng=ng1,ng2
              nopts=nopts+1
              IF(PRESS_MAX.GT.1.0D-6) THEN
                PTS(2,nopts)=OFFSET+PRESS(ng,ne)*FACTOR/PRESS_MAX
              ELSE
                PTS(2,nopts)=OFFSET+PRESS(ng,ne)
              ENDIF
            ENDDO
          ENDDO
          NTPTS=nopts
          DO nopts=1,NTPTS
            PTS(1,nopts)=X1+(X2-X1)*DBLE(nopts-1)/DBLE(NTPTS-1)
          ENDDO
          CALL POLYLINE(INDEX,iw,NTPTS,PTS,ERROR,*9999)
          WRITE(CHAR,'(I2)') N1CHORD
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,4),ERROR,*9999)
        ENDDO
        IF(NTCHORD.GT.1) THEN
          PTS(1,1)= 0.0D0
          PTS(2,1)=-0.5D0
          PTS(1,2)= 0.5D0
          PTS(2,2)= 0.0D0
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          PTS(1,1)= 1.0D0
          PTS(2,1)=-0.5D0
          PTS(1,2)= 0.75D0
          PTS(2,2)= 0.0D0
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
        ENDIF

      ELSE IF(TYPE(1:4).EQ.'FILE') THEN
C **    Draw axes
        PTS(1,1)= 0.0D0
        PTS(2,1)= 0.0D0
        PTS(1,2)= 1.0D0
        PTS(2,2)= 0.0D0
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
        PTS(1,1)= 0.0D0
        PTS(2,1)=-1.0D0
        PTS(1,2)= 0.0D0
        PTS(2,2)= 1.0D0
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)

C **    Read history file
        YMAX=0.0D0
        nohist=0
        nx=1 !Temporary
        nr=1 !needs fixing
        nrc=2 !temporary
 200      nohist=nohist+1
          READ(9,'('' YP(ny,1) at t='',E11.4,'' :'')',END=400) T
          READ(9,'(10E13.5)') (YP(ny,1),ny=1,NYT(nrc,1,nx))
          READ(9,'('' YP(ny,8) at t='',E11.4,'' :'')') T
          READ(9,'(10E13.5)') (YP(ny,8),ny=1,NYT(nrc,1,nx))
          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          IF(nohist.EQ.1) THEN
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              PTS(1,nolist)=XP(1,nv,1,np)-XP(1,nv,1,NPLIST(1))
            ENDDO
            IF(DABS(PTS(1,NPLIST(0))).LE.1.0D-6) THEN
              DO nolist=1,NPLIST(0)
                np=NPLIST(nolist)
                PTS(1,nolist)=XP(1,nv,2,np)-XP(1,nv,2,NPLIST(1))
              ENDDO
            ENDIF
            IF(DABS(PTS(1,NPLIST(0))).GT.1.0D-6) THEN
              DO nolist=1,NPLIST(0)
                PTS(1,nolist)=PTS(1,nolist)/PTS(1,NPLIST(0))
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' PTS(1,'',I2,'')='',E11.3)')
     '              nolist,PTS(1,nolist)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                PTS2(1,1)=PTS(1,nolist)
                PTS2(2,1)=0.0D0
                PTS2(1,2)=PTS2(1,1)
                PTS2(2,2)=0.03D0
                CALL POLYLINE(INDEX,iw,2,PTS2,ERROR,*9999)
                np=NPLIST(nolist)
                WRITE(CHAR,'(I3)') np
                CALL STRING_TRIM(CHAR,IBEG,IEND)
                CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS2(1,1),ERROR,*9999)
              ENDDO
            ENDIF
          ENDIF
          IF(nohist.LE.2) THEN
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              IF(DABS(ZP(1,1,1,np,nc)).GT.YMAX)
     '          YMAX=DABS(ZP(1,1,1,np,nc))
            ENDDO
          ENDIF
          DO nolist=1,NPLIST(0)
            np=NPLIST(nolist)
            PTS(2,nolist)=ZP(1,1,1,np,nc)
          ENDDO
          IF(DABS(YMAX).GT.1.0D-6) THEN
            DO nolist=1,NPLIST(0)
              PTS(2,nolist)=PTS(2,nolist)*FACTOR/YMAX
            ENDDO
          ENDIF
          CALL POLYLINE(INDEX,iw,NPLIST(0),PTS,ERROR,*9999)
          GO TO 200
 400    CONTINUE

C **    Draw section plots
        PTS(1,2)=0.02D0
        DO I=1,5
          PTS(1,1)=0.0D0
          PTS(2,1)=DBLE(I-1)/2.0D0-1.0D0
          PTS(2,2)=PTS(2,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          WRITE(CHAR,'(F7.2)') PTS(2,1)*YMAX/FACTOR
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          PTS(1,1)=-0.01D0
          CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
        ENDDO
        WRITE(CHAR,'(I2)') NTSECT
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        PTS(1,1)=0.5D0
        PTS(2,1)=-1.0D0
        CALL TEXT(1,iw,'Solution at section '//CHAR(IBEG:IEND),
     '    PTS(1,1),ERROR,*9999)

      ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
        PTS(1,1)= 0.0D0
        PTS(2,1)= 0.0D0
        PTS(1,2)= 1.0D0
        PTS(2,2)= 0.0D0
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
        PTS(1,1)= 0.0D0
        PTS(2,1)=-PI/2.0D0
        PTS(1,2)= 0.0D0
        PTS(2,2)= PI/2.0D0
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
        DO I=1,2
          PTS(1,1)=DBLE(I)/2.0D0
          PTS(2,1)=0.0D0
          PTS(1,2)=PTS(1,1)
          PTS(2,2)=0.1D0
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          X(1)=DBLE(I)/2.0D0
          IF(I.EQ.1) THEN
            CHAR='endo'
          ELSE IF(I.EQ.2) THEN
            CHAR='epi'
          ENDIF
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          PTS(1,1)=DBLE(I-1)
          PTS(2,1)=-0.1D0
          CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
        ENDDO
        DO I=-1,1,2
          PTS(1,1)=0.0D0
          PTS(2,1)=PI/2.0D0*DBLE(I)
          PTS(1,2)=0.03D0
          PTS(2,2)=PTS(2,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          WRITE(CHAR,'(I3)') 90*I
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          PTS(1,1)=-0.1D0
          CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
        ENDDO
C ***   Plot fibre angles for elements in order of incr. radial coord.
        DO noelem=1,NTELEM
          ne=NELEM(noelem)
          nr=NRE(ne)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          NB1=NBJ(1,ne)
          NB4=NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
          NOFIBR=0
          DO I=1,21
            NSFIBR=NOFIBR
            XI(3)=DBLE(I-1)/20.0D0
C           Find Xi coords at theta (XD(3)) & mu (XD(2)) for bilinear basis
            XXE(1,1)=(1.0D0-XI(3))*XE(1,3)+XI(3)*XE(5,3) !is 1st theta on 2D plane
            XXE(2,1)=(1.0D0-XI(3))*XE(2,3)+XI(3)*XE(6,3) !is 2nd theta on 2D plane
            XXE(3,1)=(1.0D0-XI(3))*XE(3,3)+XI(3)*XE(7,3) !is 3rd theta on 2D plane
            XXE(4,1)=(1.0D0-XI(3))*XE(4,3)+XI(3)*XE(8,3) !is 4th theta on 2D plane
            XXE(1,2)=(1.0D0-XI(3))*XE(1,2)+XI(3)*XE(5,2) !is 1st mu on 2D plane
            XXE(2,2)=(1.0D0-XI(3))*XE(2,2)+XI(3)*XE(6,2) !is 2nd mu on 2D plane
            XXE(3,2)=(1.0D0-XI(3))*XE(3,2)+XI(3)*XE(7,2) !is 3rd mu on 2D plane
            XXE(4,2)=(1.0D0-XI(3))*XE(4,2)+XI(3)*XE(8,2) !is 4th mu on 2D plane
            CALL XYCOORD(XD(3),XD(2),XXE,XI(1),XI(2),ERROR,*9999)
            WRITE(OP_STRING,'('' Xi: '',3F7.3)') (XI(ni),ni=1,3)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C ***       Find lamda and abscissa for plot
            RLAMDA(I)=PXI(IBT(1,1,NB1),IDO(1,1,0,NB1),INP(1,1,NB1),
     '        NB1,1,XI,XE(1,1))
            PTS(1,I)=(DBLE(2-noelem)+XI(3))/DBLE(NTELEM)
C ***       If mu>0 plot data within DTHETA & DMU
            IF(I.GT.1.AND.I.LT.21.AND.DMU.GT.0.0D0) THEN
              XD(1)=RLAMDA(I)
              DLAMDA=(RLAMDA(I)-RLAMDA(I-1))/2.0D0
              WRITE(OP_STRING,'('' XD:'',3E12.3)') (XD(nj),nj=1,3)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nd=1,NDT
                CALL ZX(ITYP10(nr),ZD(1,nd),X) !X= fibre pos.n in prolate coords
                IF(X(1).LE.XD(1)+DLAMDA.AND.X(1).GT.XD(1)-DLAMDA  .AND
     '            .X(2).LE.XD(2)+DMU   .AND.X(2).GT.XD(2)-DMU     .AND
     '            .X(3).LE.XD(3)+DTHETA.AND.X(3).GT.XD(3)-DTHETA) THEN
                  WRITE(OP_STRING,
     '              '('' Found nd='',I5,'' at ZD:'',3E12.3,'
     '              //''' Fibre angle='',E12.3)') nd,(ZD(nj,nd),nj=1,4)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  IF(NOFIBR.LT.200) THEN !store fibre position on plot
                    NOFIBR=NOFIBR+1
                    FIBR(1,NOFIBR)=PTS(1,I)
                    FIBR(2,NOFIBR)=ZD(4,nd)
                  ENDIF
                ENDIF
              ENDDO
c             WRITE(IOOP,'('' dmu='',E12.3,'' dtheta='',E12.3)') DMU,DTHETA
              IF((NOFIBR-NSFIBR).GT.0) THEN
                IF(TYPE(6:9).EQ.'FILE') THEN
                  WRITE(9,'('' Fibre data at transmural position i= '','
     '              //'I2)') I
                  WRITE(9,'(1X,20I6)') (INT(FIBR(2,N)*180.0D0/PI),
     '              N=NSFIBR+1,NOFIBR)
                  SUM=0.0D0
                  DO N=NSFIBR+1,NOFIBR
                    SUM=SUM+FIBR(2,N)
                  ENDDO
                  WRITE(9,'('' Average fibre angle = '',I6)')
     '              INT(SUM/DBLE(NOFIBR-NSFIBR)*180.0D0/PI)
                ENDIF
              ELSE
                WRITE(9,'('' no data at transmural position    i= '','
     '            //'I2)') I
              ENDIF
            ENDIF
            IF(NJ_LOC(NJL_FIBR,0,nr).GT.0) THEN
              ETA=PXI(IBT(1,1,NB4),IDO(1,1,0,NB4),INP(1,1,NB4),NB4,1,XI,
     '          XE(1,NJ_LOC(NJL_FIBR,1,nr)))
            ELSE IF(NJ_LOC(NJL_FIBR,0,nr).EQ.2) THEN
            ENDIF
            PTS(2,I)=ETA
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' ne='',I3,'' xi(1)='',E11.3,'' xi(2)='',E11.3,'
     '          //''' xi(3)='',E11.3,'' x='',E11.4,'' y='',E11.3,'
     '          //''' lamda='',E11.3)')
     '          ne,XI(1),XI(2),XI(3),PTS(1,I),PTS(2,I),RLAMDA(I)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO

          CALL POLYLINE(INDEX,iw,21,PTS,ERROR,*9999)
          IF(NOFIBR.GT.0) THEN
            CALL POLYMARKER(INDEX,iw,NOFIBR,FIBR,ERROR,*9999)
          ENDIF
          WRITE(CHAR,'(I2)') NTSECT
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,10),ERROR,*9999)
          IF(TYPE(6:9).EQ.'FILE') THEN
            WRITE(9,'('' Element '',I3)') ne
            WRITE(9,'(1X,21F6.3)') (PTS(1,I),I=1,21)
            WRITE(9,'(1X,21F6.3)') (RLAMDA(I),I=1,21)
            WRITE(9,'(1X,21I6)') (INT(PTS(2,I)*180.0D0/PI),I=1,21)
          ENDIF
        ENDDO
      ENDIF

      CALL CLOSE_SEGMENT(ISSECT,iw,ERROR,*9999)

      CALL EXITS('SGSECT')
      RETURN
 9999 CALL ERRORS('SGSECT',ERROR)
      CALL EXITS('SGSECT')
      RETURN 1
      END


