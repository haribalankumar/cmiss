      SUBROUTINE CROSS_SECTION_AREA(IBT,IDO,INP,NBJ,NEELEM,NKJE,NPF,      

     '  NPNE,nsurface_reg,NVJE,RADIUS,SEC_AREA,SE,TOLERANCE,XA,XE,XNORM,      
     '  XP,XPOINT,ERROR,*)

C#### Fuction: CROSS_SECTION_AREA
C###  Description:
C###    returns the cross sectional area of the intersection of
C###    a plane in point normal form with the surface mesh in the      
C###    region specified. 
C**** Created by Peter Bier, April 2003
      
      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      
!     parameters
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),      
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),nsurface_reg,
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 RADIUS,SEC_AREA,SE(NSM,NBFM,NEM),TOLERANCE,XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XNORM(3),XP(NKM,NVM,NJM,NPM),XPOINT(3)
      CHARACTER ERROR*(*)
      
!     local variables
      INTEGER ne,noelem
      REAL*8 AREA_FROM_ELEM

      CALL ENTERS('CROSS_SECTION_AREA',*9999)

C     calculate contribution to area from each element
      SEC_AREA = 0.0d0
      noelem = 0
      
      DO WHILE (noelem .LT. NEELEM(0,nsurface_reg))
        
        noelem = noelem + 1 
        ne = NEELEM(noelem,nsurface_reg) ! element to try

C       load element information into XE array
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nsurface_reg,NVJE(1,1,1,ne),SE(1,1,ne),
     '    XA(1,1,ne),XE,XP,ERROR,*9999)                           

        CALL CROSS_SEC_AREA_FROM_ELEM(IBT,IDO,INP,NBJ,ne,AREA_FROM_ELEM,        
     '    RADIUS,TOLERANCE,XE,XNORM,XPOINT,ERROR,*9999)
        SEC_AREA = SEC_AREA + AREA_FROM_ELEM
        
      ENDDO

      WRITE(OP_STRING,*) 'Area of cross section is',SEC_AREA
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('CROSS_SECTION_AREA')
      RETURN
 9999 CALL ERRORS('CROSS_SECTION_AREA',ERROR)
      CALL EXITS('CROSS_SECTION_AREA')
      RETURN
      END

      
