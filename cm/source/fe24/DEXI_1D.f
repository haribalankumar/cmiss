      SUBROUTINE DEXI_1D(IBT,IDO,INP,LD,nb,NBJ,
     '  NBH,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,nj2,
     '  NKHE,NKJE,NKTB,NPF,NPNE,nolist,nr,NRE,ns1,ns2,NVHE,NVJE,
     '  NW,nx,START_ELEMENT,A1,A2,CURVCORRECT,D1,D2,DIST,
     '  SE,SLOPE,SQ,SQND,X1,X2,XD,XI,XI_1,XI_2,XI_3,
     '  XA,XE,XID,XP,ZA,ZD,ZP,DEFORM,FOUND,NEW,ORTHOG,
     '  SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,ERROR,*)

C#### Subroutine: DEXI_1D
C###  Description:
C###    DEXI_1D

      IMPLICIT NONE
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
     '  NPNE(NNM,NBFM,NEM),nolist,nr,NRE(NEM),ns1,ns2,
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  START_ELEMENT
      REAL*8 A1,A2,CURVCORRECT(2,2,NNM,NEM),D1,D2,DIST,
     '  SE(NSM,NBFM,NEM),SLOPE,SQ(NDM),SQND,XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XI_1,XI_2,XI_3,X1,X2,XD(3),XI(3),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORM,FOUND,NEW,ORTHOG,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY
      CHARACTER ERROR*(*)

      CALL ENTERS('DEXI_1D',*9999)

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

          DO nolist=1,NELIST(0)
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

            IF(ITYP10(1).EQ.2) THEN !cylindrical coords
              nb=NBJ(NJ2,ne) !is basis for theta
            ELSE
              nb=NBJ(NJ1,ne)
            ENDIF
            NKTB=NKT(0,nb)
            ns1=1
            ns2=1+NKTB
            A1=XE(ns2,NJ1)-XE(ns1,NJ1)
            IF(NJT.GT.2) THEN
              nb=NBJ(NJ2,ne)
              NKTB=NKT(0,nb)
              ns1=1
              ns2=1+NKTB
              A2=XE(ns2,NJ2)-XE(ns1,NJ2)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' Element '',I3)') ne
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' A1,A2='',2E11.3)')
     '          A1,A2
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL ZX(ITYP10(1),ZD(1,nd),XD)
            IF(NJT.EQ.2.AND.ITYP10(1).EQ.1) THEN !2D rect.cart.
              D1=XD(NJ1)-XE(ns1,NJ1)
              D2=XD(NJ2)-XE(ns1,NJ2)
              IF(DABS(A2).LT.1.0D-6) THEN
                XI(1)=D1/A1
              ELSE IF(DABS(A1).LT.1.0D-6) THEN
                XI(1)=D2/A2
              ELSE
                SLOPE=A2/A1
                XI(1)=(SLOPE*D2+D1)/(SLOPE*A2+A1)
              ENDIF
              IF(XI(1).LE.1.0D0.AND.XI(1).GE.0.0D0) THEN
                X1=XE(ns1,NJ1)*(1.0D0-XI(1))+XE(ns2,NJ1)*
     '            XI(1)
                X2=XE(ns1,NJ2)*(1.0D0-XI(1))+XE(ns2,NJ2)*
     '            XI(1)
                DIST=(X1-XD(NJ1))**2+(X2-XD(NJ2))**2
                IF(LD(nd).EQ.0.OR.DIST.LT.SQ(nd)) THEN
                  XID(1,nd)=XI(1)
                  LD(nd)=ne
                  SQ(nd)=DIST
                ENDIF
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,
     '            '('' ne='',I4,'' nd='',I4,'
     '            //''' Xi(1)='','//'E11.3,'' SQ(nd)='','
     '            //'E11.3)') ne,nd,XI(1),SQ(nd)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

            ELSE IF(NJT.EQ.2.AND.ITYP10(1).EQ.2) THEN !2D cyl.polar
              XI(1)=(XD(2)-XE(ns1,2))/
     '          (XE(ns2,2)-XE(ns1,2))
              IF(XI(1).LE.1.0D0.AND.XI(1).GE.0.0D0) THEN
                X1=XE(ns1,1)*(1.0D0-XI(1))+XE(ns2,1)*XI(1)
                DIST=(X1-XD(1))**2
                IF(LD(nd).EQ.0.OR.DIST.LT.SQ(nd)) THEN
                  XID(1,nd)=XI(1)
                  LD(nd)=ne
                  SQ(nd)=DIST
                ENDIF
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,
     '            '('' ne='',I4,'' nd='',I4,'' Xi(1)='','
     '            //'E11.3,'' SQ(nd)='',E11.3)')
     '            ne,nd,XI(1),SQ(nd)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ELSE IF(NJT.EQ.3.OR.ITYP10(1).NE.1) THEN
            ENDIF

          ENDDO !nolist
        ENDIF ! not found
      ENDDO ! nd=1,NDT loop

      CALL EXITS('DEXI_1D')
      RETURN
 9999 CALL ERRORS('DEXI_1D',ERROR)
      CALL EXITS('DEXI_1D')
      RETURN 1
      END


