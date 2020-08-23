      SUBROUTINE UPDATA_VOLUME(IBT,IDO,INP,LD,LD_NP,NAN,NBH,NBJ,NELIST,
     &  NFF,NHE,njj_initial,njj_volume,NKHE,NKJE,NLL,NPF,NPNE,nr,NRE,
     &  NVHE,NVJE,NW,nx,CURVCORRECT,PG,RG,SE,WG,XA,XAB,XE,XG,XID,XP,ZA,
     &  ZD,ZE,ZG,ZP,DEFORMED,ERROR,*)

C####  Subroutine: UPDATA_VOLUME
C###   Description:
C###     UPDATA_VOLUME updates volume estimates at data points for 1.
C###     undeformed geometry by assigning mean volumes, or 2. deformed
C###     geometry by calculating sqrt(DET(AZL)) at the data point
C###     locations.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'geom00.cmn'
!     Parameter values
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  LD(NDM),LD_NP(NDM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     &  NBJ(NJM,NEM),NELIST(0:NEM),NFF(6,NEM),NHE(NEM),njj_initial,
     &  njj_volume,NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     &  NLL(12,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),nr,NRE(NEM),
     &  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     &  SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),XAB(NORM,NEM),
     &  XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     &  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZE(NSM,NHM),ZG(NHM,NUM),
     &  ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORMED
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb,nc,nd,ne,ne_last,ni,NITB,nj_initial,nj_volume,noelem,
     &  VOLTC(NBFM)
      REAL*8 AZ,AZU(3,3),AZL(3,3),DETM,dXidNu(3,3),ratio_volume_change,
     &  VOLAVG,VOL(NBFM),VOLT(NBFM),VOLTTL,XI(3)
!     Functions
      REAL*8 DET

      CALL ENTERS('UPDATA_VOLUME',*9999)

C***  Check validity and set field numbers
      CALL ASSERT(njj_initial.LE.NJ_LOC(NJL_FIEL,0,nr),
     &  '>> FROM field not defined',ERROR,*9999)
      CALL ASSERT(njj_volume.LE.NJ_LOC(NJL_FIEL,0,nr),
     &  '>> TO field not defined',ERROR,*9999)
      nj_initial=NJ_LOC(NJL_FIEL,njj_initial,nr)
      nj_volume=NJ_LOC(NJL_FIEL,njj_volume,nr)

      DO nb=1,NBFT
        VOLT(nb)=0.0d0
        VOLTC(nb)=0
      ENDDO
      VOLAVG=0.d0

C***  For each data point the mapped node is given by LD_NP(nd), the
C***  element that contains nd is stored in LD(nd). XI location of data
C***  points stored in XID(ni,nd).

      IF(DEFORMED)THEN !update using deformed coordinates
C***    Calculate total deformed volume, for estimate of error        
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBH(NH_LOC(1,nx),1,ne)
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     &      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     &      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,
     &      *9999)
          CALL OPELEMD(NBH(1,1,ne),NBJ(1,ne),ne,NFF(1,ne),NHE(ne),
     &      NKHE,NLL(1,ne),NPNE(1,1,ne),nr,NRE,NVHE(1,1,1,ne),
     &      nx,VOLTC,PG,RG,SE(1,1,ne),VOL,VOLT,WG,ZE,ZG,.TRUE.,ERROR,
     &      *9999)
        ENDDO !noelem

        nc=1 !temporary
        ne_last=0
        VOLTTL=0.d0
        DO nd=1,NDT
          ne=LD(nd) !element that 'hosts' point nd
          nb=NBH(NH_LOC(1,nx),1,ne)
          NITB=NIT(nb)
c          np=LD_NP(nd) !node that 'maps' to point nd
          DO ni=1,NIT(nb) !for the number of Xi coordinates in ne
            XI(ni)=XID(ni,nd)
          ENDDO !ni
          IF(ne_last.NE.ne)THEN !only if different host element from last nd
C***        Get undeformed element information
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &        NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C***        Get deformed element information
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     &        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     &        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,
     &        *9999)
          ENDIF !ne_last
          
C***      Interpolate XE at XI, return in XG        
          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,XG,XI,ERROR,*9999)
          CALL DXIDXM(NITB,nr,dXidNu,DETM,XG,'Fibre',ERROR,*9999)
C***      Interpolate dependent variables ZG and derivs wrt Nu (JP=1)
          CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH(1,nc,ne),NHE(ne),nr,nx,
     &      dXidNu,ZE,ZG,XI,ERROR,*9999)
C***      Evaluate the right Cauchy-Green tensor (AZL) at Xi        
          CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)
C***      Ratio of deformed:undeformed volume = sqrt(det(AZL))
          ratio_volume_change=DSQRT(DET(AZL))
C***      Calculate the deformed volume        
          ZD(nj_volume,nd)=ratio_volume_change*ZD(nj_initial,nd)

          VOLTTL=VOLTTL+ZD(nj_volume,nd)
          ne_last=ne
        ENDDO !nd

        WRITE(OP_STRING,'('' Error in volume='',F8.4,''%'')')
     &    (VOLTTL-VOLT(nb))/VOLT(nb)*100.d0
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
      ELSE !use average volume
        !calculate total volume, divide by number of data points
        VOLAVG=0.d0
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBH(NH_LOC(1,nx),1,ne)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &      NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          CALL OPELEM1(.FALSE.,IBT,NBJ(1,ne),ne,NFF(1,ne),nj_volume,
     &      NKJE(1,1,1,ne),NLL(1,ne),NPNE(1,1,ne),nr,NRE,NVJE(1,1,1,ne),
     &      VOLTC,PG,RG,SE(1,1,ne),VOL,VOLT,WG,XAB,XE,XG,.TRUE.,ERROR,
     &      *9999)
        ENDDO !noelem

        VOLAVG=VOLT(nb)/NDT
        DO nd=1,NDT
          ZD(nj_volume,nd)=VOLAVG
        ENDDO !nd
      ENDIF
      
      CALL EXITS('UPDATA_VOLUME')
      RETURN
 9999 CALL ERRORS('UPDATA_VOLUME',ERROR)
      CALL EXITS('UPDATA_VOLUME')
      RETURN 1
      END



