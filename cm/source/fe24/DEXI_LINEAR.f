      SUBROUTINE DEXI_LINEAR(IBT,IDO,INP,LD,NBJ,
     '  NBH,nb,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,nj2,
     '  NKHE,NKJE,NKTB,NPF,NPNE,nolist,nr,NRE,ns1,ns2,ns3,ns4,
     '  NVHE,NVJE,NW,nx,START_ELEMENT,A1,A2,ALFA,B1,B2,
     '  BETA,C1,C2,CURVCORRECT,D1,D2,DELTA,DENOM1,DENOM2,
     '  DIFF,GAMA,SE,SQ,SQND,THETAMIN,THETAMAX,X0,X1,X2,
     '  XA,XD,XI,XE,XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,
     '  DEFORM,EXCLUDE,EXTRAPOLATE,FOUND,NEW,ORTHOG,
     '  SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,ERROR,*)

C#### Subroutine: DEXI_LINEAR
C###  Description:
C###    DEXI_LINEAR

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),nb,NBJ(NJM,NEM),NBH(NHM,NCM,NEM),nd,ND0,ND1,ne,
     '  NELIST(0:NEM),NHE(NEM,NXM),ni,NITB,nj1,nj2,
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NKTB,NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),nolist,nr,NRE(NEM),ns1,ns2,ns3,ns4,
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  START_ELEMENT
      REAL*8 A1,A2,ALFA,B1,B2,BETA,C1,C2,CURVCORRECT(2,2,NNM,NEM),D1,
     '  D2,DELTA,DENOM1,DENOM2,DIFF,GAMA,SE(NSM,NBFM,NEM),SQ(NDM),
     '  SQND,THETAMIN,THETAMAX,X0,X1,X2,XA(NAM,NJM,NEM),XD(3),
     '  XE(NSM,NJM),XI(3),XI_1,XI_2,XI_3,XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORM,EXCLUDE,EXTRAPOLATE,FOUND,NEW,ORTHOG,
     '  SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 PXI

      CALL ENTERS('DEXI_LINEAR',*9999)

      DO nd=ND0,ND1 !new AAY 23 March 95 limits set with ADD keyword
        IF(DOP) THEN
          WRITE(OP_STRING,'('' Data point '',I5)') nd
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        FOUND=.FALSE.

        IF(.NOT.NEW) THEN ! start Xi from previous values
          CALL DEXI_EXISTING(IBT,IDO,INP,LD,NBJ,
     '      NBH,nd,ne,NHE,ni,NITB,NJ1,NKHE,NKJE,NPF,NPNE,nr,
     '      NRE,NVHE,NVJE,NW,nx,START_ELEMENT,
     '      CURVCORRECT,SE,SQ,SQND,XA,XE,XI,XI_1,
     '      XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,FOUND,ORTHOG,SET_XI_1,
     '      SET_XI_2,SET_XI_3,SPECIFY,ERROR,*9999)
        ENDIF !not new

        IF(.NOT.FOUND) THEN !try and find the Xi point by
                            !looking at all the elements
C GMH 30/9/95 If specified, then notify that we have to search in
C             other elements
          IF(SPECIFY) THEN
            OP_STRING(1)=
     '        'Projection not found on specified '
     '        //'element - searching all elements'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          LD(nd)=0
          SQ(nd)=0.0D0
          nolist=0
          DO WHILE ((LD(nd).EQ.0).AND.(nolist.NE.NELIST(0)))
            nolist=nolist+1
            ne=NELIST(nolist)
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' Element '',I4)') ne
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            NITB=NIT(NBJ(NJ1,ne))
C GMH 30/10/95 Adding call to orthog projection procedure, accounting
C for deformed coordinates correctly (I think)
C GMH 30/10/95 Should NBJ be NBH for deformed?
            IF(DEFORM)THEN
c             Note that deformed coords --> XE
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '          NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '          NVHE(1,1,1,ne),NW(ne,1),nx,
     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '          XE,ZP,ERROR,*9999)
            ELSE
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),
     '          NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            ENDIF

            nb=NBJ(NJ1,ne) !nb for first geometric coord
            NKTB=NKT(0,nb)
            ns1=1
            ns2=1+NKTB
            ns3=1+2*NKTB
            ns4=1+3*NKTB
            A1=XE(ns2,NJ1)-XE(ns1,NJ1)
            B1=XE(ns3,NJ1)-XE(ns1,NJ1)
            C1=XE(ns1,NJ1)-XE(ns2,NJ1)-XE(ns3,NJ1)+
     '        XE(ns4,NJ1)
            nb=NBJ(NJ2,ne)
            NKTB=NKT(0,nb)
            ns1=1
            ns2=1+NKTB
            ns3=1+2*NKTB
            ns4=1+3*NKTB
            A2=XE(ns2,NJ2)-XE(ns1,NJ2)
            B2=XE(ns3,NJ2)-XE(ns1,NJ2)
            C2=XE(ns1,NJ2)-XE(ns2,NJ2)-XE(ns3,NJ2)+
     '        XE(ns4,NJ2)
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' Element '',I3)') ne
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' A1,B1,C1='',3E11.3)')
     '          A1,B1,C1
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' A2,B2,C2='',3E11.3)')
     '          A2,B2,C2
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            ALFA=A2*C1-A1*C2
            IF(LD(nd).EQ.0) THEN
              CALL ZX(ITYP10(1),ZD(1,nd),XD)
              EXCLUDE=.TRUE.
              IF(ITYP10(1).EQ.4) THEN !check if lies within element
                THETAMIN=DMIN1(XE(ns2,3),XE(ns4,3))
                THETAMAX=DMAX1(XE(ns1,3),XE(ns3,3))
                IF(XD(3).LT.THETAMAX.AND.
     '            XD(3).GT.THETAMIN) THEN
                  EXCLUDE=.FALSE.
                  D2=XD(NJ2)-XE(1,NJ2)
                ELSE IF((XD(3)+2.0D0*PI).LT.THETAMAX
     '              .AND.(XD(3)+2.0D0*PI).GT.THETAMIN) THEN
                  EXCLUDE=.FALSE.
                  D2=XD(NJ2)+2.0D0*PI-XE(1,NJ2)
                ELSE IF((XD(3)-2.0D0*PI).LT.THETAMAX
     '              .AND.(XD(3)-2.0D0*PI).GT.THETAMIN) THEN
                  EXCLUDE=.FALSE.
                  D2=XD(NJ2)-2.0D0*PI-XE(1,NJ2)
                ENDIF
              ELSE
                EXCLUDE=.FALSE.
                D2=XD(NJ2)-XE(1,NJ2)
              ENDIF
              IF(.NOT.EXCLUDE) THEN
                D1=XD(NJ1)-XE(1,NJ1)
                BETA=C2*D1-C1*D2+A2*B1-A1*B2
                GAMA=B2*D1-B1*D2
                DELTA=BETA*BETA-4.0D0*ALFA*GAMA
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' D1='',E11.3,'
     ,              //''' D2='',E11.3,'
     '              //''' ALFA='',E11.3,'' BETA='',E11.3,'
     '              //''' GAMA='',E11.3,'' DELTA='',E11.3'
     '              //')')D1,D2,ALFA,BETA,GAMA,DELTA
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(DABS(ALFA).GT.1.0D-6) THEN
                  IF(DELTA.GE.0.0D0) THEN
                    XI(1)=(-BETA+DSQRT(DELTA))/
     '                (2.0D0*ALFA)
                    IF(XI(1).LT.0.0D0.OR.
     '                XI(1).GT.1.0D0) THEN
                      XI(1)=(-BETA-DSQRT(DELTA))/
     '                  (2.0D0*ALFA)
                    ENDIF
                  ELSE
                     XI(1) = -1.0D0
                  ENDIF
                ELSE IF(DABS(ALFA).LE.1.0D-6) THEN
                  IF(DABS(BETA).GT.1.0D-6) THEN
                    XI(1)=-GAMA/BETA
                  ELSE
                    XI(1)=-1.0D0
                  ENDIF
                ENDIF
                IF(XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0) THEN
                  DENOM1=B1+C1*XI(1)
                  IF(DABS(DENOM1).GT.1.0D-6) THEN
                    XI(2)=(D1-A1*XI(1))/DENOM1
                  ELSE
                    DENOM2=B2+C2*XI(1)
                    IF(DABS(DENOM2).GT.1.0D-6) THEN
                      XI(2)=(D2-A2*XI(1))/DENOM2
                    ELSE
                      WRITE(OP_STRING,
     '                  '('' Xi(2) cannot be defined for '
     '                  //'nd='',I6,'' in element '',I5)')
     '                  nd,ne
                      CALL WRITES(IOFI,OP_STRING,ERROR,
     '                  *9999)
                    ENDIF
                  ENDIF
                  IF(XI(2).GE.0.0D0.AND
     '              .XI(2).LE.1.0D0) THEN
                    IF(NITB.EQ.2) THEN
                      IF(DOP) THEN
                        WRITE(OP_STRING,
     '                    '('' nd='',I4,'' Xi:'',3E11.3,'
     '                    //''' (linear calc.)'')')
     '                    nd,(XI(ni),ni=1,3)
                        CALL WRITES(IODI,OP_STRING,ERROR,
     '                    *9999)
                      ENDIF
                      nb=NBJ(NJ1,ne)
                      IF(NKT(0,nb).GT.1) THEN !update Xi(2)
                        XI(2)=0.0D0
                        X1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,1,XI,XE(1,NJ1))
                        XI(2)=1.0D0
                        X2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,1,XI,XE(1,NJ1))
                        IF(DABS(X2-X1).GT.1.0D-6) THEN !PJH 13APR91
                          XI(2)=(XD(NJ1)-X1)/(X2-X1)
                        ELSE
                          WRITE(OP_STRING,
     '                      '('' Warning!!!'',''X1=X2'')')
                          CALL WRITES(IOOP,OP_STRING,
     '                      ERROR,*9999)
                          XI(2)=1.0D0
                        ENDIF
                      ENDIF
                      XID(1,nd)=XI(1)
                      XID(2,nd)=XI(2)
                      LD(nd)=ne
                    ELSE IF(NITB.EQ.3) THEN
                      nb=NBJ(1,ne)
                      XI(3)=0.0D0
                      X0=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,1,XI,XE(1,1))
                      XI(3)=1.0D0
                      X1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,1,XI,XE(1,1))
                      DIFF=X1-X0
                      IF(DABS(DIFF).GT.1.0D-6) THEN
                        XI(3)=(XD(1)-X0)/DIFF
                      ELSE
                        XI(3)=0.0D0
                      ENDIF
                      IF(DOP) THEN
                        WRITE(OP_STRING,
     '                    '('' X0='',E11.3,'
     '                    //''' X1='',E11.3,'' XD(1)='','
     '                    //'E11.3,'' Xi(3)='',E11.3)')
     '                    X0,X1,XD(1),XI(3)
                        CALL WRITES(IODI,OP_STRING,ERROR,
     '                    *9999)
                      ENDIF
                      IF(XI(3).GT.0.0D0.AND.
     '                  XI(3).LE.1.0D0) THEN
                        LD(nd)=ne
                      ELSE IF(EXTRAPOLATE.AND.
     '                    XI(3).GT.1.0D0.AND.
     '                    XI(3).LE.1.1D0) THEN
                        WRITE(OP_STRING,
     '                    '('' TAGGING: Found nd='',I6,'
     '                    //''' with xi(3)='',F5.2)')
     '                    nd,XI(3)
                        CALL WRITES(IOOP,OP_STRING,
     '                    ERROR,*9999)
                        IF(LD(nd).EQ.0) THEN
                          LD(nd)=-ne !to tag for 'extrapolate' below
                        ENDIF
                      ENDIF
                      XID(1,nd)=XI(1)
                      XID(2,nd)=XI(2)
                      XID(3,nd)=XI(3)
!news                             AAY use set xi_3 option if required 23 March 95
                      IF(SET_XI_3) THEN
                        XI(3)=XI_3
                      ENDIF
!newe
                    ENDIF
                    IF(DOP) THEN
                      WRITE(OP_STRING,
     '                  '('' nd='',I4,'' Xi:'',3E11.3)')
     '                  nd,(XI(ni),ni=1,3)
                      CALL WRITES(IODI,OP_STRING,ERROR,
     '                  *9999)
                    ENDIF
                  ENDIF !0<xi(2)<1
                ENDIF !0<xi(1)<1
              ENDIF !.NOT.EXCLUDE
            ENDIF !LD(nd)=0

          ENDDO !nolist
        ENDIF ! not found
      ENDDO ! nd=1,NDT loop

      CALL EXITS('DEXI_LINEAR')
      RETURN
 9999 CALL ERRORS('DEXI_LINEAR',ERROR)
      CALL EXITS('DEXI_LINEAR')
      RETURN 1
      END

