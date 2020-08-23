      SUBROUTINE CALC_LATT_XQ(IBT,IDO,INP,NAN,NBJ,NENQ,NKJE,
     '  NPF,NPNE,nr,NVJE,DXDXIQ2,SE,XA,XE,XIQ,
     '  XP,XQ,ERROR,*)

C#### Subroutine: CALC_LATT_XQ
C###  Description:
C###    CALC_LATT_XQ calculates the array XQ when the lattice
C###    grid scheme is used.      

      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'      
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),NENQ(0:8,NQM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 DXDXIQ2(3,3,NQM),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XIQ(NIM,NQM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)      

!     Local Variables
      INTEGER ne,nj,njj1,njj2,nq
      REAL*8 XQD(NJM,NUM),XQ_TEMP(3)
      
      CALL ENTERS('CALC_LATT_XQ',*9999)

      nr=1
      
      DO nq=1,NQT
        ne=NENQ(1,nq)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '    XA(1,1,ne),XE,XP,ERROR,*9999)
        CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),
     '    nr,XE,XQD,XIQ(1,nq),ERROR,*9999)
        DO nj=1,NJT
          XQ_TEMP(nj)=XQD(nj,1)
          DXDXIQ2(nj,1,nq)=XQD(nj,2)
          DXDXIQ2(nj,2,nq)=XQD(nj,4)
          IF(NJT.EQ.3) DXDXIQ2(nj,3,nq)=XQD(nj,7)
        ENDDO !nj        
        CALL XZ(ITYP10(nr),XQ_TEMP,XQ(1,nq))        
        njj1=NJL_FIBR !fibres
                        
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          XQ(nj,nq)=XQD(nj,1)
        ENDDO !nj-fibre
      ENDDO

      CALL EXITS('CALC_LATT_XQ')
      RETURN
 9999 CALL ERRORS('CALC_LATT_XQ',ERROR)
      CALL EXITS('CALC_LATT_XQ')
      RETURN 1
      END


