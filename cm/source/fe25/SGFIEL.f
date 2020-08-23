      SUBROUTINE SGFIEL(IBT,IDO,INP,ISEG,ISFIEL,iw,MXI,NB_TYPE,NBJ,ne,
     '  NGAP,NLATNE,NQNE,NQNLAT,NQS,NQXI,CSEG,TYPE,XE,XQ,YG,YQ,ZE,
     '  ERROR,*)
      

C#### Subroutine: SGFIEL
C###  Description:
C###    SGFIEL creates field segment ISFIEL.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISFIEL,iw,MXI(2),NB_TYPE,NBJ(NJM),NGAP(NIM),
     '  NLATNE(NEQM+1),NQNE(NEQM,NQEM),NQNLAT(NEQM*NQEM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM)
      REAL*8 XE(NSM,NJM),XQ(NJM,NQM),YG(NIYGM,NGM),YQ(NYQM),ZE(NSM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER i,i1,i2,INDEX,INDEX_OLD,ip,j1,j2,
     '  n1,n2,n3,n4,ni1,ni2,nb,ne,ng,nj,no_points1,no_points2,nq,
     '  NQMAX,NQMAX2
      REAL*8 POINTS(3,40),PXI,X(3),XI(3),ZVAL,ZZP(4,400)

      CALL ENTERS('SGFIEL',*9999)
      CALL OPEN_SEGMENT(ISFIEL,ISEG,iw,'FIEL',1,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

      IF(iw.EQ.4) THEN
        PROJEC=MAP_PROJEC
        IF(PROJEC(1:2).EQ.'XI') THEN
          MXI1=MXI(1) !bottom left coords
          MXI2=MXI(2) !..for Xi map projection
        ENDIF
      ENDIF

      IF(NIT(NB_TYPE).EQ.1) THEN
        IF(TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:9).EQ.'DEPENDENT') THEN
          XI(1)=0.0D0
          DO nj=1,NJT !to find position of 1st point on line
            nb=NBJ(nj)
            POINTS(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XI,XE(1,nj))
          ENDDO
          DO i=2,20
            XI(1)=DBLE(i-1)/19.0D0
            DO nj=1,NJT !to find position of next pt on line
              nb=NBJ(nj)
              POINTS(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '          1,XI,XE(1,nj))
            ENDDO
            nb=NB_TYPE
            ZVAL=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,ZE)
            IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
            IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
            IF(COLOUR_WS) THEN
              INDEX=256-NINT((ZVAL-ZMINI)/ZDIFF*215.0D0)
            ELSE
              INDEX=8-NINT((ZVAL-ZMINI)/ZDIFF*3.0D0)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' I='',I3,'' Value='',E12.3,'' INDEX='',I3)')
     '          i,ZVAL,INDEX
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL POLYLINE(INDEX,iw,2,POINTS,ERROR,*9999)
            POINTS(1,1)=POINTS(1,2)
            POINTS(2,1)=POINTS(2,2)
          ENDDO

        ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN
          DO nj=1,3
            X(nj)=0.d0
          ENDDO
          DO ng=1,NGAP(1)-1
            XI(1)=DBLE(ng-1)/DBLE(NGAP(1)-1)
            XI(2)=DBLE(ng)/DBLE(NGAP(1)-1)
            DO nj=1,NJT !to find position of ends of line
              nb=NBJ(nj)
              X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XI(1),XE(1,nj))
            ENDDO
            CALL XZ(ITYP10(1),X,POINTS(1,1))
            DO nj=1,NJT !to find position of ends of line
              nb=NBJ(nj)
              X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XI(2),XE(1,nj))
            ENDDO
            CALL XZ(ITYP10(1),X,POINTS(1,2))
            ZVAL=(YG(1,ng)+YG(1,ng+1))/2.d0
            IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
            IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
            IF(COLOUR_WS) THEN
              INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
            ELSE
              INDEX=8-NINT((ZVAL-ZMINI)/ZDIFF*3.0D0)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Value='',E12.3,'' INDEX='',I3)')
     '          ZVAL,INDEX
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
C MLB 10/7/97 this call produces incorrect colours if INDEX drops
C below 41, this is reserved for particular colours.
C            CALL POLYLINE(INDEX,iw,2,POINTS,ERROR,*9999)
            CALL FILL_LINE(INDEX,iw,2,POINTS,ERROR,*9999)
          ENDDO !ng

C MLB 15/7/97 adding grid option
        ELSE IF(TYPE(1:4).EQ.'GRID') THEN
          NQMAX=NQXI(1,NQS(ne))-1
          DO nq=1,NQMAX
            DO nj=1,NJT
              POINTS(nj,1)=XQ(nj,NQNE(ne,nq))
              POINTS(nj,2)=XQ(nj,NQNE(ne,nq+1))
            ENDDO
            ZVAL=(YQ(NQNE(ne,nq))+YQ(NQNE(ne,nq+1)))/2.0d0
            IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
            IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
            IF(COLOUR_WS) THEN
              INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
            ELSE
              INDEX=8-NINT((ZVAL-ZMINI)/ZDIFF*3.0D0)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Value='',E12.3,'' INDEX='',I3)')
     '          ZVAL,INDEX
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FILL_LINE(INDEX,iw,2,POINTS,ERROR,*9999)
          ENDDO !nq

C*****  temp removed until have additional qualifier on 'de fie;s'
c       IF(JTYP14.EQ.0) THEN !draw shaded area either side of line element
!         Calculate fill area coordinates
c         DO I=1,20
c           XI(1)=DBLE(I-1)/19.0D0
c           DO nj=1,NJT !to find position on line
c             nb=NBJ(nj)
c             X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,XE(1,nj))
c           ENDDO
c           DO nj=1,NJT !to find Xi-deriv of line
c             nb=NBJ(nj)
c             DX(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,XE(1,nj))
c           ENDDO
c           DD=DSQRT(DX(1)**2+DX(2)**2)
c           nb=NB_TYPE
c           AREA=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,ZE)
c           IF(AREA.GT.0.0D0) THEN
c             RADIUS=DSQRT(AREA/PI)
c           ELSE
c             RADIUS=0.0D0
c           ENDIF
c           XINCR(1)=RADIUS*DX(2)/DD !is increment in x
c           XINCR(2)=RADIUS*DX(1)/DD !is increment in y
c           POINTS(1,I)   =X(1)+XINCR(1) !point on one side of line
c           POINTS(2,I)   =X(2)-XINCR(2)
c           POINTS(1,41-I)=X(1)-XINCR(1) !.. & point on other side of line
c           POINTS(2,41-I)=X(2)+XINCR(2)
c         ENDDO
c         IF(COLOUR_WS) THEN
c           INDEX=100
c         ELSE
c           INDEX=8
c         ENDIF
c         IF(DOP) WRITE(*,'('' Radius='',E12.3,'' INDEX='',I3)') RADIUS,INDEX
c         CALL FILL_AREA(INDEX,iw,40,POINTS,ERROR,*9999)
        ENDIF !type

      ELSE IF(NIT(NB_TYPE).GT.1) THEN !draw 20*20 mesh of fill areas
        DO i=1,3
          XI(i)=0.0d0
        ENDDO
        IF(iw.EQ.2) THEN
          ni1=2
          ni2=3
        ELSEIF(iw.EQ.3) THEN
          ni1=1
          ni2=3
        ELSE
          ni1=1
          ni2=2
        ENDIF
        IF(TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:9).EQ.'DEPENDENT') THEN
          IF(IBT(1,1,NB_TYPE).NE.3) THEN !Not Simplex
!           Calculate vertex coordinates ip=1..400
            DO i2=1,20
              XI(ni2)=DBLE(i2-1)/19.0D0
              DO i1=1,20
                XI(ni1)=DBLE(i1-1)/19.0D0
                ip=i1+20*(i2-1)
                IF(iw.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
                  X(ni1)=XI(ni1)
                  X(ni2)=XI(ni2)
                ELSE
                  DO nj=1,NJT
                    nb=NBJ(nj)
                    X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,1,XI,XE(1,nj))
                  ENDDO
                ENDIF
C CPB 1/11/92 Allow draw field for postscript
                IF(iw.LE.3.OR.iw.EQ.15) THEN !convert curvilinear to r.c. coords
                  CALL XZ(ITYP10(1),X,ZZP(1,ip))
                ELSE IF(iw.EQ.4) THEN !leave as curvi. coords (since call map4)
                  ZZP(1,ip)=X(1) !lamda for prolate
                  ZZP(2,ip)=X(2) !mu    for prolate
                  ZZP(3,ip)=X(3) !theta for prolate
                ENDIF
              ENDDO
            ENDDO
!           Calculate fill areas ie=1..361
            DO j2=1,19
              XI(ni2)=0.025d0+(0.95d0/18.0d0)*DBLE(j2)
              DO j1=1,19
                XI(ni1)=0.025d0+(0.95d0/18.0d0)*DBLE(j1)
                n1=j1+20*(j2-1) !is bottom left vertex
                n2=n1+1         !is bottom right
                n3=n2+20        !is top right
                n4=n1+20        !is top left
                DO nj=1,NJT
                  POINTS(nj,1)=ZZP(nj,n1)
                  POINTS(nj,2)=ZZP(nj,n2)
                  POINTS(nj,3)=ZZP(nj,n3)
                  POINTS(nj,4)=ZZP(nj,n4)
                ENDDO
                ZVAL=PXI(IBT(1,1,NB_TYPE),IDO(1,1,0,NB_TYPE),
     '            INP(1,1,NB_TYPE),NB_TYPE,1,XI,ZE)
                IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
                IF(COLOUR_WS) THEN
                  INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' Index='',I3)') ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
                  INDEX=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' Index='',I2)') ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                CALL FILL_AREA(INDEX,iw,4,POINTS,ERROR,*9999)
              ENDDO
            ENDDO
C CS 8/8/97 new draw field for simplex elements (Area coordinates).
C Divides simplex into 100 triangular areas and fills them.
          ELSE !Simplex
!           Calculate vertex coordinates
            ip=1
            DO i1=0,10
              XI(1)=DBLE(i1)/10.0d0
              DO i2=0,(10-i1)
                IF(i1.eq.10) THEN
                  XI(2)=0.0d0
                ELSE
                  XI(2)=DBLE(i2)/10.0d0
                ENDIF
                XI(3)=0.0d0
                IF(iw.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
                  X(1)=XI(1)
                  X(2)=XI(2)
                ELSE
                  DO nj=1,NJT
                    nb=NBJ(nj)
                    X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,1,XI,XE(1,nj))
                  ENDDO
                ENDIF
                ZVAL=PXI(IBT(1,1,NB_TYPE),IDO(1,1,0,NB_TYPE),
     '            INP(1,1,NB_TYPE),NB_TYPE,1,XI,ZE)
                ZZP(1,ip)=X(1)
                ZZP(2,ip)=X(2)
                ZZP(3,ip)=X(3)
                ZZP(4,ip)=ZVAL
                ip=ip+1
              ENDDO
            ENDDO
!           Calculate fill areas
            no_points1=0              !odd areas
            no_points2=0
            DO i1=0,9
              no_points1=no_points2
              DO i2=0,(9-i1)
                n1=(i2+1)+no_points1  !is bottom left vertex
                n2=n1+1               !is bottom right
                n3=n1+11-i1           !is top
                DO nj=1,NJT
                  POINTS(nj,1)=ZZP(nj,n1)
                  POINTS(nj,2)=ZZP(nj,n2)
                  POINTS(nj,3)=ZZP(nj,n3)
                ENDDO
                no_points2=no_points2+1
                ZVAL= ZZP(4,n3)
                IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
                IF(COLOUR_WS) THEN
                  INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' Index='',I3)') ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
                  INDEX=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' Index='',I2)') ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                CALL FILL_AREA(INDEX,iw,3,POINTS,ERROR,*9999)
              ENDDO
                no_points2=no_points2+1
            ENDDO
            no_points1=0                  !even areas
            no_points2=0
            DO i1=0,8
              no_points1=no_points2
              DO i2=0,(8-i1)
                n1=(i2+1)+no_points1+1    !is bottom vertex
                n2=n1+11-i1               !is top left
                n3=n1+11-i1-1             !is top right
                DO nj=1,NJT
                  POINTS(nj,1)=ZZP(nj,n1)
                  POINTS(nj,2)=ZZP(nj,n2)
                  POINTS(nj,3)=ZZP(nj,n3)
                ENDDO
                no_points2=no_points2+1
                ZVAL= ZZP(4,n3)
                IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
                IF(COLOUR_WS) THEN
                  INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' Index='',I3)') ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
                  INDEX=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '              '('' ZVAL='',E12.3,'' Index='',I2)') ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                CALL FILL_AREA(INDEX,iw,3,POINTS,ERROR,*9999)
              ENDDO
                no_points2=no_points2+2
            ENDDO
          ENDIF

        ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN
!         Calculate vertex coordinates ip=1..
          DO i2=1,NGAP(2)+1
            XI(2)=DBLE(i2-1)/DBLE(NGAP(2))
            DO i1=1,NGAP(1)+1
              XI(1)=DBLE(i1-1)/DBLE(NGAP(1))
              ip=i1+(NGAP(1)+1)*(i2-1)
              IF(iw.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
                X(1)=XI(1)
                X(2)=XI(2)
              ELSE
                DO nj=1,NJT
                  nb=NBJ(nj)
                  X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '              XI,XE(1,nj))
                ENDDO
              ENDIF
              IF(iw.LE.3) THEN !convert curvilinear to r.c. coords
                CALL XZ(ITYP10(1),X,ZZP(1,ip))
              ELSE IF(iw.EQ.4) THEN !leave as curvi. coords (since call map4)
                ZZP(1,ip)=X(1) !lamda for prolate
                ZZP(2,ip)=X(2) !mu    for prolate
                ZZP(3,ip)=X(3) !theta for prolate
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,
     '            '('' I2='',I2,'' I1='',I2,'' IP='',I3)') i2,i1,ip
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
!         Calculate fill areas
          ng=0
          DO j2=1,NGAP(2) !loops over Gauss pts in Xi(2) direction
            DO j1=1,NGAP(1)   !loops over Gauss pts in Xi(1) direction
              n1=j1+(NGAP(1)+1)*(j2-1) !is bottom left vertex
              n2=n1+1                  !is bottom right
              n3=n2+NGAP(1)+1          !is top right
              n4=n1+NGAP(1)+1          !is top left
              DO nj=1,NJT
                POINTS(nj,1)=ZZP(nj,n1)
                POINTS(nj,2)=ZZP(nj,n2)
                POINTS(nj,3)=ZZP(nj,n3)
                POINTS(nj,4)=ZZP(nj,n4)
              ENDDO
              ng=ng+1
              ZVAL=YG(1,ng)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' J2='',I2,'' J1='',I2,'' ng='',I3,'
     '            //''' N:'',4I3,'' ZVAL='',E12.4)')
     '            j2,j1,ng,n1,n2,n3,n4,ZVAL
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
              IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
              IF(COLOUR_WS) THEN
                INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ZVAL='',E12.3,'' Index='',I3)')
     '              ZVAL,INDEX
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ELSE
                INDEX=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ZVAL='',E12.3,'' Index='',I2)')
     '              ZVAL,INDEX
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
              CALL FILL_AREA(INDEX,iw,4,POINTS,ERROR,*9999)
            ENDDO
          ENDDO

C MLB 16/7/97 adding grid option
        ELSE IF(TYPE(1:4).EQ.'GRID') THEN
!         Calculate vertex coordinates ip=1..
          NQMAX=NQXI(1,NQS(ne))
          NQMAX2=NQXI(2,NQS(ne))
          DO j2=1,NQMAX2-1    !loops over Grid pts in Xi(2) direction
            DO j1=1,NQMAX-1   !loops over Grid pts in Xi(1) direction
              n1=j1+NQMAX*(j2-1)   !is bottom left vertex
              n2=n1+1              !is bottom right
              n3=n2+NQMAX          !is top right
              n4=n1+NQMAX          !is top left
              DO nj=1,NJT
C DAH adding in support for lattice method
                IF(USE_LAT.EQ.1) THEN
                  POINTS(nj,1)=XQ(nj,NQNLAT(NLATNE(ne)+n1-1))
                  POINTS(nj,2)=XQ(nj,NQNLAT(NLATNE(ne)+n2-1))
                  POINTS(nj,3)=XQ(nj,NQNLAT(NLATNE(ne)+n3-1))
                  POINTS(nj,4)=XQ(nj,NQNLAT(NLATNE(ne)+n4-1))
                ELSE
                  POINTS(nj,1)=XQ(nj,NQNE(ne,n1))
                  POINTS(nj,2)=XQ(nj,NQNE(ne,n2))
                  POINTS(nj,3)=XQ(nj,NQNE(ne,n3))
                  POINTS(nj,4)=XQ(nj,NQNE(ne,n4))
                ENDIF
              ENDDO
              IF(USE_LAT.EQ.1) THEN
                ZVAL=0.25d0*(YQ(NQNLAT(NLATNE(ne)+n1-1))+
     '            YQ(NQNLAT(NLATNE(ne)+n2-1))+
     '            YQ(NQNLAT(NLATNE(ne)+n3-1))
     '            +YQ(NQNLAT(NLATNE(ne)+n4-1)))
              ELSE
                ZVAL=0.25d0*(YQ(NQNE(ne,n1))+YQ(NQNE(ne,n2))+
     '            YQ(NQNE(ne,n3))+YQ(NQNE(ne,n4)))                  
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' J2='',I2,'' J1='',I2,'
     '            //''' N:'',4I3,'' ZVAL='',E12.4)')
     '            j2,j1,n1,n2,n3,n4,ZVAL
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
              IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
              IF(COLOUR_WS) THEN
                INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ZVAL='',E12.3,'' Index='',I3)')
     '              ZVAL,INDEX
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ELSE
                INDEX=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ZVAL='',E12.3,'' Index='',I2)')
     '              ZVAL,INDEX
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
              CALL FILL_AREA(INDEX,iw,4,POINTS,ERROR,*9999)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      CALL CLOSE_SEGMENT(ISFIEL,iw,ERROR,*9999)

      CALL EXITS('SGFIEL')
      RETURN
 9999 CALL ERRORS('SGFIEL',ERROR)
      CALL EXITS('SGFIEL')
      RETURN 1
      END


