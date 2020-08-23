      SUBROUTINE DEXI_ORTHOG(IBT,IDO,INP,LD,NBJ,
     '  NBH,nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,
     '  NKHE,NKJE,NPF,NPNE,nolist,nr,NRE,NVHE,NVJE,NW,nx,
     '  START_ELEMENT,CURVCORRECT,SE,SQ,SQND,XA,XE,XI,
     '  XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,
     '  NEW,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,
     '  ERROR,*)

C#### Subroutine: DEXI_ORTHOG
C###  Description:
C###    DEXI_ORTHOG

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),NBJ(NJM,NEM),NBH(NHM,NCM,NEM),nd,ND0,ND1,ne,
     '  NELIST(0:NEM),NHE(NEM,NXM),ni,NITB,nj1,NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),nolist,
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx,START_ELEMENT
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  SQ(NDM),SQND,XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XI_1,XI_2,XI_3,
     '  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORM,FOUND,NEW,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY
      CHARACTER ERROR*(*)

      CALL ENTERS('DEXI_ORTHOG',*9999)

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
     '      XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,FOUND,.TRUE.,SET_XI_1,
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

            DO ni=1,NITB !Start in the center of the element
              XI(ni)=0.5D0
            ENDDO
            IF(SET_XI_3) THEN
              XI(3)=XI_3
            ENDIF
            FOUND=.FALSE.
            CALL PROJ_ORTHOG(IBT,IDO,INP,NBJ(1,ne),SQND,
     '        XE,XI,ZD(1,nd),FOUND,ERROR,*9999)
            IF(FOUND) THEN
              IF(LD(nd).EQ.0.OR.SQND.LT.SQ(nd)) THEN
                DO ni=1,NITB !copy new xi values
                  XID(ni,nd)=XI(ni)
                ENDDO
                LD(nd)=ne
                SQ(nd)=SQND
              ENDIF
            ENDIF

          ENDDO !nolist
        ENDIF ! not found
      ENDDO ! nd=1,NDT loop

      CALL EXITS('DEXI_ORTHOG')
      RETURN
 9999 CALL ERRORS('DEXI_ORTHOG',ERROR)
      CALL EXITS('DEXI_ORTHOG')
      RETURN 1
      END


