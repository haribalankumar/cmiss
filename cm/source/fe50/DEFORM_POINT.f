      SUBROUTINE DEFORM_POINT(CURVCORRECT,IBT,IDO,INP,LD,NBH,NBJ,NEELEM,
     '  NHE,NKJE,NKHE,NPF,NPNE,NRE,NVHE,NVJE,NW,SE,SUCCESS,XA,XD,XE,XID,
     '  XP,ZA,ZDlocal,ZP,ERROR,*)

C##### SUBROUTINE: DEFORM_POINT
C####  DESCRIPTION:
C####    Takes a point defined in the undeformed coordinates and calculates that point
C####    in the deformed coordinates

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      LOGICAL SUCCESS
      CHARACTER ERROR*(*)
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),LD,
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM,NXM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XD(3),XE(NSM,NJM),XID(NIM),
     '  XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZDlocal(3)
!     Local Variables
      INTEGER nb,ne,nj,nr,nx,noelem,LDTEMP
      REAL*8 PXI,XI(3)

      CALL ENTERS('DEFORM_POINT',*9999)
      nr=1
      nx=1
      SUCCESS=.FALSE.
      LDTEMP=LD
      
      IF(LDTEMP.NE.0) THEN ! test the suggest solution
          DO nj=1,3
            XI(nj)=XID(nj)
          ENDDO
          CALL XPXE(NBJ(1,LDTEMP),NKJE(1,1,1,LDTEMP),
     '      NPF(1,1),NPNE(1,1,LDTEMP),NRE(LDTEMP),
     '      NVJE(1,1,1,LDTEMP),SE(1,1,LDTEMP),XA,
     '      XE,XP,ERROR,*9999)
          CALL DEXI_POINT(IBT,IDO,INP,LDTEMP,NBJ,
     '      LD,NIT(NBJ(1,LDTEMP)),NRE(LDTEMP),0.d0,XE,XI,XID,ZDlocal,
     '      .FALSE.,ERROR,*9999)
C           WRITE(OP_STRING,'(''XI: '',3E11.3)') XI(1),XI(2),XI(3)
C           CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO noelem=1,NEELEM(0,nr) ! test all other elements
        ne=NEELEM(noelem,nr)
        IF(LDTEMP.EQ.0.AND.ne.NE.LD) THEN
          DO nj=1,3
            XI(nj) = 0.5d0
          ENDDO
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '      NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '      NVJE(1,1,1,ne),SE(1,1,ne),XA,
     '      XE,XP,ERROR,*9999)
          CALL DEXI_POINT(IBT,IDO,INP,LDTEMP,NBJ,
     '      ne,NIT(NBJ(1,ne)),NRE(ne),0.d0,XE,XI,XID,ZDlocal,
     '      .FALSE.,ERROR,*9999)
        ENDIF
      ENDDO
C      WRITE(OP_STRING,'(''LDTEMP: '',I5)') LDTEMP
C      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      LD=LDTEMP
      IF(LDTEMP.NE.0) THEN
        CALL ZPZE(NBH(1,1,LDTEMP),1,NHE(LDTEMP,nx),
     '    NKHE(1,1,1,LDTEMP),NPF(1,1),NPNE(1,1,LDTEMP),nr,
     '    NVHE(1,1,1,LDTEMP),NW(LDTEMP,1,nx),nx,
     '    CURVCORRECT(1,1,1,LDTEMP),SE(1,1,LDTEMP),
     '    ZA(1,1,1,LDTEMP),XE,ZP,ERROR,*9999)
        DO nj=1,3
          nb=NBJ(nj,LDTEMP)
          XD(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '      INP(1,1,nb),nb,1,XID,XE(1,nj))
        ENDDO
        SUCCESS=.TRUE.
      ENDIF

      CALL EXITS('DEFORM_POINT')
      RETURN
 9999 CALL ERRORS('DEFORM_POINT',ERROR)
      CALL EXITS('DEFORM_POINT')
      RETURN 1
      END

