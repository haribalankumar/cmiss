      SUBROUTINE IPMOTI_LUNG(NBJ,NEELEM,NELIST,NELIST2,NENP,NJJ_COEFF,
     &  NJJ_FLOW,NJJ_PAO2,NJJ_PBO2,NJJ_PRESSURE,NJJ_RADIUS,NJJ_SOURCE,
     &  NKJ,NKJE,NORD,NPNE,NPNODE,nr,NVJE,NVJP,NXI,BBM,CE,FRC,MINPCNT,
     &  TLC,XAB,XP,BELOW,LUMPED_PARAMETER,SET_FRC,SET_TLC,SCALE,
     &  VMIN_STATE,ERROR,*)

C#### Subroutine: IPMOTI_LUNG
C###  Description:
C###    IPMOTI_LUNG inputs motion parameters for pulmonary problems in
C###    region nr.

C XAB(1) is percentage volume change
C XAB(2) is volume_dv
C XAB(3) is volume of element      
C XAB(4) is volume below the element

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),
     &  NELIST2(0:NEM),NENP(NPM,0:NEPM),NJJ_COEFF,NJJ_FLOW,NJJ_PAO2,
     &  NJJ_PBO2,NJJ_PRESSURE,NJJ_RADIUS,NJJ_SOURCE,NKJ(NJM,NPM),
     &  NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M),nr,NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),FRC,MINPCNT,TLC,XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM)
      LOGICAL BELOW,LUMPED_PARAMETER,SET_FRC,SET_TLC,VMIN_STATE
      CHARACTER SCALE*(10),ERROR*(*)
!     Local Variables
      INTEGER M,nb,ne,ne0,N_ELEM,NELIST_LOCAL(0:20000),NE_OLD(2000),
     &  NE_TEMP(2000),noelem,noelem2,np1,np2,NT_BNS,nv0,
     &  nv1,nv2
      REAL*8 dV_acinus,rad_scale,vol_acinus

      CALL ENTERS('IPMOTI_LUNG',*9999)

      nj_coeff=NJ_LOC(NJL_FIEL,NJJ_COEFF,nr)

C.....Store the volume below (and including) branch ne in CE(4,ne)  
      write(*,*) 'callling mesh volume 1'
      CALL MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,BBM,CE,XP,ERROR,
     &  *9999)

C.....dV(branch)=dV(acinus)*vol(branch)/vol(acinus)
      DO noelem=1,NELIST(0) !for each parent in the list
        ne0=NELIST(noelem)
        nb=NBJ(1,ne0)
        np2=NPNE(2,nb,ne0) !end node
        nv2=NVJE(2,nb,nj_flow,ne0)
        vol_acinus = CE(4,ne0)-CE(3,ne0)
        rad_scale=DSQRT(BBM(1,ne0)/vol_acinus)
        dV_acinus = XP(1,nv2,nj_flow,np2)
        NELIST_LOCAL(0)=0
        DO M=1,NXI(1,0,ne0)
          NE_OLD(M)=NXI(1,M,ne0)
          NELIST_LOCAL(0)=NELIST_LOCAL(0)+1
          NELIST_LOCAL(NELIST_LOCAL(0))=NXI(1,M,ne0)
        ENDDO
        NT_BNS=NXI(1,0,ne0)
        DO WHILE(NT_BNS.NE.0)
          N_ELEM=NT_BNS
          NT_BNS=0
          DO M=1,N_ELEM
            ne=NE_OLD(M)
            XAB(1,ne)=dV_acinus*CE(3,ne)/vol_acinus
            np1=NPNE(1,nb,ne) !start node
            nv1=NVJE(1,nb,nj_flow,ne)
            XP(1,nv1,nj_flow,np1)=XAB(1,ne) !initialise flow=dV
            DO noelem2=1,NXI(1,0,ne) !for each child
              NT_BNS=NT_BNS+1
              NE_TEMP(NT_BNS)=NXI(1,noelem2,ne)
              NELIST_LOCAL(0)=NELIST_LOCAL(0)+1
              NELIST_LOCAL(NELIST_LOCAL(0))=NXI(1,noelem2,ne)
            ENDDO
          ENDDO !M
          DO M=1,NT_BNS
            NE_OLD(M)=NE_TEMP(M)
          ENDDO !M
        ENDDO !while
        DO noelem2=NELIST_LOCAL(0),1,-1
          ne=NELIST_LOCAL(noelem2)
          np1=NPNE(1,nb,ne) !start node
          np2=NPNE(2,nb,ne) !end node
          nv1=NVJE(1,nb,nj_flow,ne)
          nv2=NVJE(2,nb,nj_flow,ne)
          XAB(2,ne)=XAB(1,ne) !initialise the dV_below
          XP(1,nv2,nj_flow,np2)=XP(1,nv1,nj_flow,np1)
          DO M=1,NXI(1,0,ne)
            nv0=NVJE(1,nb,nj_flow,NXI(1,M,ne)) !version in child
            IF(NORD(5,NXI(1,M,ne)).EQ.1)THEN
              XP(1,nv2,nj_flow,np2)=XP(1,nv2,nj_flow,np2)+2.d0*
     &          XP(1,nv0,nj_flow,np2)
              XAB(2,ne)=XAB(2,ne)+2.d0*XAB(2,NXI(1,M,ne))
            ELSE
              XP(1,nv2,nj_flow,np2)=XP(1,nv2,nj_flow,np2)+
     &          XP(1,nv0,nj_flow,np2)
              XAB(2,ne)=XAB(2,ne)+XAB(2,NXI(1,M,ne))
            ENDIF
          ENDDO
          XP(1,nv1,nj_flow,np1)=XP(1,nv2,nj_flow,np2)
! scale the size of the acinus
          nv1=NVJE(1,nb,nj_radius,ne)
          nv2=NVJE(2,nb,nj_radius,ne)
          XP(1,nv1,nj_radius,np1)=XP(1,nv1,nj_radius,np1)*rad_scale
          XP(1,nv2,nj_radius,np2)=XP(1,nv2,nj_radius,np2)*rad_scale
        ENDDO !local list of elements
        XAB(2,ne0)=XAB(1,ne0)
        DO M=1,NXI(1,0,ne0)
          IF(NORD(5,NXI(1,M,ne0)).EQ.1)THEN
            XAB(2,ne0)=XAB(2,ne0)+2.d0*XAB(2,NXI(1,M,ne0))
          ELSE
            XAB(2,ne0)=XAB(2,ne0)+XAB(2,NXI(1,M,ne0))
          ENDIF
        ENDDO
        np2=NPNE(2,nb,ne0) !end node
        
        BBM(1,ne0)=0.d0
      ENDDO !for each parent element

      write(*,*) 'callling mesh volume 2'
      CALL MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,BBM,CE,XP,ERROR,
     &  *9999)

      CALL EXITS('IPMOTI_LUNG')
      RETURN
 9999 CALL ERRORS('IPMOTI_LUNG',ERROR)
      CALL EXITS('IPMOTI_LUNG')
      RETURN 1
      END


