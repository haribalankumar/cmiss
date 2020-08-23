      SUBROUTINE AREAFIELD(IBT,IDO,INP,NBJ,NCENTRE_REG,NDDATA,NDP,            
     '  NEELEM,NKJE,NPF,NPOINTS_PER_ELEM,NPNE,NSURFACE_REG,NVJE,
     '  SE,RADIUS_OF_INTEREST,TOLERANCE,WD,XA,XE,XP,ZD,
     '  ERROR,*)

C#### Subroutine: AREAFIELD
C###  Description:
C###    AREAFIELD calculates the cross sectional area for a tubular mesh
C###    along the centreline specified
C**** Created by Peter Bier, Mar 2003      

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
C      INCLUDE 'cmiss$reference:call00.cmn'

!     parameters
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),      
     '  NBJ(NJM,NEM),NCENTRE_REG,NDDATA(0:NDM,0:NRM),NDP(NDM),
     '  NEELEM(0:NE_R_M,0:NRM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPOINTS_PER_ELEM,NPNE(NNM,NBFM,NEM),NSURFACE_REG,
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 RADIUS_OF_INTEREST,SE(NSM,NBFM,NEM),TOLERANCE,WD(NJM,NDM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER*(*) ERROR
      
!     local variables
      INTEGER nb,nd,ne,nj,noelem,nppelem
      REAL*8 SEC_AREA,XDERIV(3),XI(3),XPOINT(3)
!     functions
      REAL*8 PXI

      CALL ENTERS('AREAFIELD',*9999)

C     existing data will be clobbered
      NDT = 0
      NDDATA(0,0) = NDT
      NDDATA(0,NCENTRE_REG) = 0
      nd = NDT
            
      DO noelem = 1,NEELEM(0,NCENTRE_REG)
        ne = NEELEM(noelem,NCENTRE_REG)

C       calculate area at npoints per element
        DO nppelem = 1,NPOINTS_PER_ELEM
          XI(1) = 1.0d0 * (nppelem-1)/NPOINTS_PER_ELEM
          XI(2) = 0.0d0
          XI(3) = 0.0d0
          nb = NBJ(1,ne)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NCENTRE_REG,NVJE(1,1,1,ne),SE(1,1,ne),
     '      XA(1,1,ne),XE,XP,ERROR,*9999)                           
          
          DO nj = 1,NJT                
            XPOINT(nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,1,XI,XE(1,nj))
            XDERIV(nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,2,XI,XE(1,nj))          
          ENDDO ! nj
          
C         get point coords and derivative vector
          CALL CROSS_SECTION_AREA(IBT,IDO,INP,NBJ,NEELEM,NKJE,NPF,NPNE,
     '      NSURFACE_REG,NVJE,RADIUS_OF_INTEREST,SEC_AREA,SE,TOLERANCE,
     '      XA,XE,XDERIV,XP,XPOINT,ERROR,*9999)
          
C       write geometric data
          nd = nd + 1
          DO nj=1,NJT
            ZD(nj,nd) = XPOINT(nj)
            WD(nj,nd) = 1.0d0
          ENDDO !nj
C         write area field data
          ZD(NJT+1,nd) = SEC_AREA
          WD(NJT+1,nd) = 1.0d0

          NDDATA(nd,NCENTRE_REG) = nd
          NDP(nd)=nd

        ENDDO !NPOINTS_PER_ELEM
      ENDDO !NEELM

      NDDATA(0,0) = nd
      NDDATA(0,NCENTRE_REG) = nd
      NDT = NDDATA(0,0)

      CALL EXITS('AREAFIELD')
      RETURN
 9999 CALL ERRORS('AREAFIELD',ERROR)
      CALL EXITS('AREAFIELD')
      RETURN 1
      END
      

