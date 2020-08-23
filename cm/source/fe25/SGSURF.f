      SUBROUTINE SGSURF(INDEX,IBT,IDO,INP,ISEG,ISSURF,iw,
     '  NAN,NBH,NBJ,NELIST,NHE,NKHE,NKJE,NOSURF,
     '  NOXIPT,NPF,NPNE,NRE,NVHE,NVJE,NW,nx,
     '  CURVCORRECT,PG,RG,SE,XA,XE,XG,XI3,XP,ZA,ZE,ZG,ZP,
     '  COLOUR,DEFORM,STATIC,CSEG,ERROR,*)

C#### Subroutine: SGSURF
C###  Description:
C###    SGSURF creates element surface segment ISSURF.
C**** 22Aug89: Modified to interpolate surface points correctly in
C****          Lagrange/Hermite tensor-product elements, using ZEZW.
C****          Added DEFORM,NAN,NHE,NW,XG,ZG,ERROR to param list. ADMcC

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
      INCLUDE 'surf00.cmn'
      INCLUDE 'scal00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISEG(*),ISSURF,iw,NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NELIST(0:NEM),NHE(NEM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NOSURF,NOXIPT(*),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XI3,
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL COLOUR,DEFORM,STATIC
!     Local Variables
      INTEGER FACET_COLOURS(1000),i,i1,i2,ICOL,ie,II(3,1000),INDEX_OLD,
     '  ip,IPAT,j1,j2,k,na,n1,n1xipt,n2,n2xipt,n3,n4,nb,NBFF,NBG,nc,ne,
     '  NFACETS,nj,nk,nn,nolist,nr,ns,NSF,NSG,NSTB2,
     '  NVERTMX
      PARAMETER(NVERTMX=1000)
C      REAL FACET_NORMAL(3,1000)
      REAL*8 A(3),B(3),C(3),CC,DXI1,DXI2,
     '  DXIX(3,3),DXIXN(3,3),DXNXI(3,3),EG(3,3),PF1,PHI(3),PST(3),PXI,
     '  RM(3,3),TEMP,TEMP2,
     '  X(3),XI(3),XI1,XI2,XIP1,XIP2,
     '  Z(3),Z2(3),ZVAL,ZZP(3,NVERTMX)

      CALL ENTERS('SGSURF',*9999)
      nc=1 !Temporary MPN 12-Nov-94

      DO i=1,3
        DO k=1,3
          DXIX(i,k)=0.0d0
        ENDDO
      ENDDO

      CALL OPEN_SEGMENT(ISSURF,ISEG,iw,'SURF',INDEX,INDEX_OLD,
     '  NOSURF,1,CSEG,ERROR,*9999)

      DO nolist=1,NELIST(0) !is main element loop
        ne=NELIST(nolist)
        nr=NRE(ne)
        nb=NBJ(1,ne)
        IF(.NOT.STATIC) THEN
          nb=NBH(NH_LOC(1,nx),1,ne)
          CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '      ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
        ENDIF
        IF(DEFORM) THEN
          nb=NBH(NH_LOC(1,nx),1,ne)
          CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '      ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
        ELSE
          nb=NBJ(1,ne)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        ENDIF

        IF(SURFACE_TYPE(1:4).EQ.'DOTS') THEN
          nb=NBJ(1,ne)
          IF(NIT(nb).EQ.3) THEN
C ***       define nodes at interpolated position given by Xi_3
            NSTB2=NST(nb)/2
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              DO ns=1,NSTB2
                XE(ns,nj)=(1.0D0-XI3)*XE(ns,nj)+XI3*XE(ns+NSTB2,nj)
              ENDDO
            ENDDO
          ENDIF
          DXI1=1.0D0/DBLE(NOXIPT(1))
          DXI2=1.0D0/DBLE(NOXIPT(2))
          XIP1=-DXI1/2.0D0
          DO n1xipt=1,NOXIPT(1)
            XIP1=XIP1+DXI1
            XIP2=-DXI2/2.0D0
            DO n2xipt=1,NOXIPT(2)
              XIP2=XIP2+DXI2
              XI(1)=XIP1
              XI(2)=XIP2
              XI(3)=XI3
              IF(DOP) THEN
                WRITE(OP_STRING,'('' xi='',3F5.2)') (XI(i),i=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(DEFORM) THEN
                CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '            NJ_LOC(NJL_GEOM,0,NRE(ne)),
     '            NRE(ne),nx,DXIX,ZE,ZG,XI,ERROR,*9999)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  X(nj)=ZG(nj,1)
                ENDDO
              ELSE
                CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),XE,XG,XI,
     '            ERROR,*9999)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  X(nj)=XG(nj,1)
                ENDDO
              ENDIF
C             DO nj=1,NJT
C               nb=NBJ(nj,ne)
C               X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C    '            XE(1,nj))
C             ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' x='',3E12.3)') (X(i),i=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL XZ(ITYP10(1),X,Z)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' z='',3E12.3)') (Z(i),i=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL POLYMARKER(5,iw,1,Z,ERROR,*9999)
            ENDDO
          ENDDO

C LC 25/2/97 archived section :
C        I don't think anyone wants to use this, delete it? AAY 18-12-90

        ELSE IF(SURFACE_TYPE(1:7).EQ.'PATTERN') THEN
          XI(3)=XI3
          IF(.NOT.STATIC) THEN
            XI(NIT(NBJ(1,ne))+1)=TIME
          ENDIF

          IF(iw.LE.2) THEN
C ***       Calculate vertex coordinates ip=1..400
            DO i2=1,20
              XI(1)=DBLE(i2-1)/19.0D0
              DO i1=1,20
                XI(2)=DBLE(i1-1)/19.0D0
                ip=i1+20*(i2-1)
                DO nj=1,NJT
                  nb=NBJ(nj,ne)
                  X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '              XI,XE(1,nj))
                ENDDO
                CALL XZ(ITYP10(1),X,ZZP(1,ip)) !convert curvilinear to r.c. coords
                IF(.NOT.STATIC) THEN
                  DO nj=1,NJT
                    nb=NBH(nj,1,ne)
                    Z2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,1,XI,ZE(1,nj))
                    IF(KTYP58(nr).EQ.2) THEN !displacement
                      ZZP(nj,ip)=ZZP(nj,ip)+Z2(nj)
                    ELSE                !take motion from ze
                      ZZP(nj,ip)=Z2(nj)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
C ***       Calculate fill areas ie=1..361
            DO j2=1,19
              XI(1)=0.05D0*DBLE(j2)
              DO j1=1,19
                ie=j1+19*(j2-1)
                n1=j1+20*(j2-1) !is bottom left vertex
                n2=n1+1         !is bottom right
                n3=n2+20        !is top right
                n4=n1+20        !is top left
                XI(2)=0.05D0*DBLE(j1)
                nb=NBJ(NJ_LOC(NJL_FIEL,1,nr),ne)
                ZVAL=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,NJ_LOC(NJL_FIEL,1,nr)))
                IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                IF(COLOUR) THEN
                  INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' ICOL='',I3)')
     '                ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
                  INDEX=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' IPAT='',I2)')
     '                ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                CALL FILL_AREA(INDEX,iw,4,ZZP,ERROR,*9999)
              ENDDO
            ENDDO

          ELSE IF(iw.EQ.3.OR.iw.EQ.16) THEN
C ***       Calculate vertex coordinates
            DO i1=1,NOXIPT(1)+1
              XI(1)=DBLE(i1-1)/DBLE(NOXIPT(1))
              DO i2=1,NOXIPT(2)+1
                XI(2)=DBLE(i2-1)/DBLE(NOXIPT(2))
                ip=I1+(NOXIPT(1)+1)*(i2-1)
                DO nj=1,3
                  nb=NBJ(nj,ne)
                  X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XI,XE(1,nj))
                ENDDO
                CALL XZ(ITYP10(1),X,ZZP(1,ip))
                IF(.NOT.STATIC) THEN
                  DO nj=1,NJT
                    nb=NBH(nj,1,ne)
                    Z2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,1,XI,ZE(1,nj))
                    IF(KTYP58(nr).EQ.2) THEN !displacement
                      ZZP(nj,ip)=ZZP(nj,ip)+Z2(nj)
                    ELSE
                      ZZP(nj,ip)=Z2(nj)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

C ***       Calculate vertex list for each triangular facet
            NFACETS=NOXIPT(1)*NOXIPT(2)*2
C            NVERTICES=(NOXIPT(1)+1)*(NOXIPT(2)+1)
            ie=0
            DO i2=1,NOXIPT(2)
              DO i1=1,NOXIPT(1)
                n1=I1+(NOXIPT(1)+1)*(i2-1)
                n2=n1+1
                n3=n1+(NOXIPT(1)+1)
                n4=n3+1
                ie=ie+1
                II(1,ie)=n1-1
                II(2,ie)=n2-1
                II(3,ie)=n4-1
                IF(DOP) THEN
                  WRITE(OP_STRING,*)
     '              'II(1..3,ie)=',II(1,ie),II(2,ie),II(3,ie)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                ie=ie+1
                II(1,ie)=n1-1
                II(2,ie)=n4-1
                II(3,ie)=n3-1
                IF(DOP) THEN
                  WRITE(OP_STRING,*)
     '              ' II(1..3,ie)=',II(1,ie),II(2,ie),II(3,ie)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO

C ***       Calculate facet normals
            DO ie=1,NFACETS
              DO nj=1,3
                A(nj)=ZZP(nj,II(2,ie)+1)-ZZP(nj,II(1,ie)+1)
                B(nj)=ZZP(nj,II(3,ie)+1)-ZZP(nj,II(1,ie)+1)
              ENDDO
              C(1)=A(2)*B(3)-A(3)*B(2)
              C(2)=A(3)*B(1)-A(1)*B(3)
              C(3)=A(1)*B(2)-A(2)*B(1)
              CC=DSQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))
              IF(CC.GT.1.0D-5) THEN
                C(1)=C(1)/CC
                C(2)=C(2)/CC
                C(3)=C(3)/CC
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Normal: '',3F6.3)')
     '            C(1),C(2),C(3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
C              DO nj=1,3
C                FACET_NORMAL(nj,ie)=REAL(C(nj))
C              ENDDO
            ENDDO

C ***       Calculate facet colours for field intensity plotting only
            IF(VARIABLE_TYPE(1:5).EQ.'FIELD') THEN
              DXI1=1.0D0/DBLE(NOXIPT(1))
              DXI2=1.0D0/DBLE(NOXIPT(2))
              ie=0
              DO i2=1,NOXIPT(2)
                XI2=DBLE(i2-1)*DXI2
                DO i1=1,NOXIPT(1)
                  XI1=DBLE(i1-1)*DXI1
                  ie=ie+1 !is lower triangle
                  XI(1)=XI1+2.0D0/3.0D0*DXI1
                  XI(2)=XI2+1.0D0/3.0D0*DXI2
                  nb=NBJ(NJ_LOC(NJL_FIEL,1,nr),ne)
                  ZVAL=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XI,XE(1,NJ_LOC(NJL_FIEL,1,nr)))
                  IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                  IF(COLOUR) THEN
                    ICOL=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                    FACET_COLOURS(ie)=ICOL
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' IE='',I3,'' ZVAL='',E12.3,'
     '                  //''' ICOL='',I3,'' FACET_COLOURS:'',I4)')
     '                  ie,ZVAL,ICOL,FACET_COLOURS(ie)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE !greyscale
                    IPAT=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' IE='',I2,'' ZVAL='',E12.3,'
     '                  //''' IPAT='',I2)') ie,ZVAL,IPAT
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF
                  ie=ie+1 !is upper triangle
                  XI(1)=XI1+1.0D0/3.0D0*DXI1
                  XI(2)=XI2+2.0D0/3.0D0*DXI2
                  nb=NBJ(NJ_LOC(NJL_FIEL,1,nr),ne)
                  ZVAL=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XI,XE(1,NJ_LOC(NJL_FIEL,1,nr)))
                  IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                  IF(COLOUR) THEN
                    ICOL=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)+17
                    FACET_COLOURS(ie)=ICOL
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' IE='',I2,'' ZVAL='',E12.3,'
     '                  //''' ICOL='',I3,'' FACET_COLOURS:'',I4)')
     '                  ie,ZVAL,ICOL,FACET_COLOURS(ie)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE !greyscale
                    IPAT=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' IE='',I2,'' ZVAL='',E12.3,'
     '                  //''' IPAT='',I2)') ie,ZVAL,IPAT
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO

            ELSE IF(VARIABLE_TYPE(1:6).EQ.'STRAIN') THEN
              DXI1=1.0D0/DBLE(NOXIPT(1))
              DXI2=1.0D0/DBLE(NOXIPT(2))
              ie=0
              IF(.NOT.STATIC) THEN
C               see sgstra for comments on the following time intpl.
                DO nj=1,NJT
                  NBFF=NBH(nj,1,ne) !time basis
                  NBG=NBJ(nj,ne)    !corresponding spatial basis
                  NSG=0
                  DO nn=1,NNT(NBG)
                    DO nk=1,NKT(0,NBG)
                      NSG=NSG+1
                      TEMP=0.0D0
                      TEMP2=0.0D0
                      DO na=1,IBT(2,NIT(NBFF),NBFF) !#terms in series
                        NSF=(nn-1)*NKT(0,NBFF)+(na-1)*NKT(0,NBG)+nk
                        TEMP=TEMP+PF1(na,1,TIME)*ZE(NSF,nj)
                        TEMP2=TEMP2+PF1(na,1,TIME_REF)*ZE(NSF,nj)
                      ENDDO
                      ZE(NSG,nj)=TEMP+XE(NSG,nj)
                      XE(NSG,nj)=TEMP2+XE(NSG,nj)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
              DO i2=1,NOXIPT(2)
                XI2=DBLE(i2-1)*DXI2
                DO i1=1,NOXIPT(1)
                  XI1=DBLE(i1-1)*DXI1
                  DO k=1,2 !two triangles
                    ie=ie+1 !is triangle(facet) number
                    IF(k.EQ.1) THEN
                      XI(1)=XI1+2.0D0/3.0D0*DXI1
                      XI(2)=XI2+1.0D0/3.0D0*DXI2
                    ELSE IF(k.EQ.2) THEN
                      XI(1)=XI1+1.0D0/3.0D0*DXI1
                      XI(2)=XI2+2.0D0/3.0D0*DXI2
                    ENDIF
                    IF(STATIC) THEN
                      CALL ZEEX51(IBT,IDO,INP,NAN,NBH(1,1,ne),
     '                  NBJ(1,ne),0,NHE(ne),nr,nx,
     '                  DXIXN,DXNXI,EG,PG,PHI,PST,RG(1),RM,
     '                  XE,XG,XI,ZE,ZG,ERROR,*9999)
                    ELSE IF(.NOT.STATIC) THEN
                      CALL ZEEX51(IBT,IDO,INP,NAN,NBJ(1,ne),
     '                  NBJ(1,ne),0,NHE(ne),nr,nx,
     '                  DXIXN,DXNXI,EG,PG,PHI,PST,RG(1),RM,
     '                  XE,XG,XI,ZE,ZG,ERROR,*9999)
                    ENDIF
                    IF(VARIABLE_TYPE(7:8).EQ.'P1') THEN !1st ppl strain
                      ZVAL=PST(1)
                    ELSE IF(VARIABLE_TYPE(7:8).EQ.'P2') THEN !2nd pl stn
                      ZVAL=PST(2)
                    ELSE IF(VARIABLE_TYPE(7:8).EQ.'PA') THEN !ppl angle
                      ZVAL=PHI(1)
                    ELSE IF(VARIABLE_TYPE(7:9).EQ.'E11') THEN !circum
                      ZVAL=EG(1,1)
                    ELSE IF(VARIABLE_TYPE(7:9).EQ.'E22') THEN !long
                      ZVAL=EG(2,2)
                    ELSE IF(VARIABLE_TYPE(7:9).EQ.'E33') THEN !trans
                      ZVAL=EG(3,3)
                    ELSE IF(VARIABLE_TYPE(7:9).EQ.'E12') THEN !shear
                      ZVAL=EG(1,2)
                    ELSE IF(VARIABLE_TYPE(7:9).EQ.'E13') THEN !shear
                      ZVAL=EG(1,3)
                    ELSE IF(VARIABLE_TYPE(7:9).EQ.'E23') THEN !shear
                      ZVAL=EG(2,3)
                    ENDIF
                    IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                    IF(COLOUR) THEN
                      ICOL=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                      FACET_COLOURS(ie)=ICOL
                      IF(DOP) THEN
                        WRITE(OP_STRING,
     '                    '('' IE='',I3,'' ZVAL='',E12.3,'
     '                    //''' ICOL='',I3,'' FACET_COLOURS:'',I4)')
     '                    ie,ZVAL,ICOL,FACET_COLOURS(ie)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE !greyscale
                      IPAT=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                      IF(DOP) THEN
                        WRITE(OP_STRING,
     '                    '('' IE='',I2,'' ZVAL='',E12.3,'
     '                    //''' IPAT='',I2)') ie,ZVAL,IPAT
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

c           CALL SURFACE(INDEX_SURF,iw,NVERTICES,ZZP,FACET_COLOURS,
c    '        FACET_NORMAL,VERTEX_COLOURS,VERTEX_NORMAL,NFACETS,NVPF,
c    '        II,NVERTMX,ERROR,*9999)

          ELSE IF(iw.EQ.4) THEN
C ***       Calculate vertex coordinates ip=1..400
            DO i2=1,20
              XI(1)=DBLE(i2-1)/19.0D0
              DO i1=1,20
                XI(2)=DBLE(i1-1)/19.0D0
                ip=i1+20*(i2-1)
                DO nj=1,NJT
                  nb=NBJ(nj,ne)
                  X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '              XI,XE(1,nj))
                ENDDO
                ZZP(1,ip)=X(2) !mu for prolate
                ZZP(2,ip)=X(3) !theta for prolate
              ENDDO
            ENDDO
C ***       Calculate fill areas ie=1..361
            DO j2=1,19
              XI(1)=0.05D0*DBLE(j2)
              DO j1=1,19
                ie=j1+19*(j2-1)
                n1=j1+20*(j2-1) !is bottom left vertex
                n2=n1+1         !is bottom right
                n3=n2+20        !is top right
                n4=n1+20        !is top left
                XI(2)=0.05D0*DBLE(j1)
                nb=NBJ(NJ_LOC(NJL_FIEL,1,nr),ne)
                ZVAL=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,NJ_LOC(NJL_FIEL,1,nr)))
                IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
                IF(COLOUR) THEN
                  INDEX=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' ICOL='',I3)')
     '                ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
                  INDEX=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '('' ZVAL='',E12.3,'' IPAT='',I2)')
     '                ZVAL,INDEX
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                PROJEC=MAP_PROJEC
                CALL FILL_AREA(INDEX,iw,4,ZZP,ERROR,*9999)
              ENDDO !jq
            ENDDO !j2
          ENDIF
        ENDIF
      ENDDO !nolist (ne) main element loop

      CALL CLOSE_SEGMENT(ISSURF,iw,ERROR,*9999)

      CALL EXITS('SGSURF')
      RETURN
 9999 CALL ERRORS('SGSURF',ERROR)
      CALL EXITS('SGSURF')
      RETURN 1
      END
