      SUBROUTINE DEXI_EXISTING(IBT,IDO,INP,LD,NBJ,
     '  NBH,nd,ne,NHE,ni,NITB,NJ1,NKHE,NKJE,NPF,NPNE,nr,
     '  NRE,NVHE,NVJE,NW,nx,START_ELEMENT,
     '  CURVCORRECT,SE,SQ,SQND,XA,XE,XI,XI_1,
     '  XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,FOUND,ORTHOG,SET_XI_1,
     '  SET_XI_2,SET_XI_3,SPECIFY,ERROR,*)

C#### Subroutine: DEXI_EXISTING
C###  Description:
C###    DEXI_EXISTING

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),NBJ(NJM,NEM),NBH(NHM,NCM,NEM),nd,
     '  ne,NHE(NEM,NXM),ni,NITB,NJ1,NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,START_ELEMENT
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  SQ(NDM),SQND,XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XI_1,XI_2,XI_3,
     '  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORM,FOUND,ORTHOG,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY
      CHARACTER ERROR*(*)

      CALL ENTERS('DEXI_EXISTING',*9999)
C
C CPB 25/1/93 If the data projection has already been found previously,
C try and find the new projection in the same element
C
C GMH 30/9/95 Pick up the specified values
      IF(SPECIFY) THEN
!Set element number
        IF(START_ELEMENT.GT.0) THEN
          LD(nd)=START_ELEMENT
        ENDIF
      ENDIF
      IF(LD(nd).NE.0) THEN
        IF(ORTHOG) THEN ! use nonlinear Xi calculation
          ne=LD(nd) ! start at the previous element
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Element '',I5)') ne
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          NITB=NIT(NBJ(NJ1,ne))
C GMH 30/10/95 Adding call to orthog projection procedure, accounting
C for deformed coordinates correctly (I think)
C GMH 30/10/95 Should NBJ be NBH for deformed?
          IF(DEFORM)THEN
C           Note that deformed coords --> XE
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '        NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '        NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '        ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
          ELSE
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),
     '        NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          ENDIF

C!!! LKC 30-MAY-2002 Be aware that these maybe uninitialised
C if there are not previous values but we want to specify and element
C to search in (and default the xi's to 0.0)
C
          DO ni=1,NITB !start from previous values
            XI(ni)=XID(ni,nd)
          ENDDO
          IF(SPECIFY) THEN !Set any xi values
            IF(SET_XI_1) THEN
              XI(1)=XI_1
            ENDIF
            IF(SET_XI_2) THEN
              XI(2)=XI_2
            ENDIF
            IF(SET_XI_3) THEN
              XI(3)=XI_3
            ENDIF
          ENDIF
          CALL PROJ_ORTHOG(IBT,IDO,INP,NBJ(1,ne),SQND,
     '      XE,XI,ZD(1,nd),FOUND,ERROR,*9999)
          IF(FOUND) THEN
            DO ni=1,NITB !copy new xi values
              XID(ni,nd)=XI(ni)
            ENDDO
            LD(nd)=ne
            SQ(nd)=SQND
          ENDIF
        ENDIF !ld(nd).NE.0
      ENDIF !orthog

      CALL EXITS('DEXI_EXISTING')
      RETURN
 9999 CALL ERRORS('DEXI_EXISTING',ERROR)
      CALL EXITS('DEXI_EXISTING')
      RETURN 1
      END
