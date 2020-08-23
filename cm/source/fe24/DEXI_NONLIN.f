      SUBROUTINE DEXI_NONLIN(IBT,IDO,INP,LD,LDTEMP,NBJ,NBH,
     '  nd,ND0,ND1,ne,NELIST,NHE,ni,NITB,nj1,NKHE,NKJE,NPF,
     '  NPNE,nolist,nr,NRE,NVHE,NVJE,NW,nx,CURVCORRECT,
     '  SE,XA,XE,XI,XID,XIQ,XP,XQ,
     '  ZA,ZD,ZP,DEFORM,FOUND,GRPGRID,PASS2,
     '  SPECIFY,ERROR,*)

C#### Subroutine: DEXI_NONLIN
C###  Description:
C###    DEXI_NONLIN

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LDTEMP,NBJ(NJM,NEM),NBH(NHM,NCM,NEM),nd,ND0,ND1,ne,
     '  NELIST(0:NEM),NHE(NEM,NXM),ni,NITB,nj1,
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),nolist,nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XID(NIM,NDM),
     '  XIQ(NIM,NQM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORM,FOUND,GRPGRID,SPECIFY
      CHARACTER ERROR*(*)
      LOGICAL PASS2
!     Local Variables
      INTEGER ld_index

      CALL ENTERS('DEXI_NONLIN',*9999)

      DO nolist=1,NELIST(0)
        ne=NELIST(nolist)
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Element '',I4)') ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        FOUND=.FALSE.
        NITB=NIT(NBJ(NJ1,ne))
C GMH 30/10/95 Adding call to orthog projection procedure, accounting
C for deformed coordinates correctly (I think)
C GMH 30/10/95 Should NBJ be NBH for deformed?
        IF(DEFORM)THEN
c         Note that deformed coords --> XE
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '      NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '      NVHE(1,1,1,ne),NW(ne,1),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '      XE,ZP,ERROR,*9999)
        ELSE
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),
     '        NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        ENDIF

        DO nd=ND0,ND1 !new AAY 23 March 95 limits set with ADD keyword
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Data point '',I5)') nd
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          IF(GRPGRID) THEN
            ld_index=1
          ELSE
            ld_index=nd
          ENDIF

          IF((((LD(ld_index).EQ.ne).OR.(LD(ld_index).EQ.0))
     '      .AND..NOT.PASS2).OR.
     '      ((LD(ld_index).EQ.0).AND.PASS2)) THEN
            IF(LD(ld_index).NE.0) THEN !initialising XI
              IF(GRPGRID) THEN
                DO ni=1,NITB
                  XI(ni)=XIQ(ni,nd)
                ENDDO
              ELSE
                DO ni=1,NITB
                  XI(ni)=XID(ni,nd)
                ENDDO
              ENDIF
            ELSE
              DO ni=1,NITB
                XI(ni)=0.5d0
              ENDDO
            ENDIF

            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' LD(nd)='',I4)') LD(nd)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(/'' Initial XI(1..)='',3F12.4)')
     '          (XI(ni),ni=1,NITB)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

C GMH 30/9/95 If specified, then notify that we have to search in
C             other elements
            IF(SPECIFY) THEN
              OP_STRING(1)=
     '          'Projection not found on specified '
     '          //'element - searching all elements'
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(GRPGRID) THEN
              CALL DEXI_POINT(IBT,IDO,INP,LDTEMP,NBJ,
     '          ne,NITB,NRE(ne),0.d0,XE,XI,XIQ(1,nd),XQ(1,nd),
     '          .FALSE.,ERROR,*9999)
            ELSE
              CALL DEXI_POINT(IBT,IDO,INP,LDTEMP,NBJ,
     '          ne,NITB,NRE(ne),0.d0,XE,XI,XID(1,nd),ZD(1,nd),
     '          .FALSE.,ERROR,*9999)
            ENDIF

            LD(ld_index)=LDTEMP
          ELSE
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Not in this element or found '')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF ! LD
        ENDDO ! nd=1,NDT loop
      ENDDO !nolist

      IF(GRPGRID) THEN
        nd=nd-1
      ENDIF

      CALL EXITS('DEXI_NONLIN')
      RETURN
 9999 CALL ERRORS('DEXI_NONLIN',ERROR)
      CALL EXITS('DEXI_NONLIN')
      RETURN 1
      END


