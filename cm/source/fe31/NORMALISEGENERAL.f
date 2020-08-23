      SUBROUTINE NORMALISEGENERAL(NBJ,NEELEM,nj,NPNE,NXI,NVJE,XP,
     &  ERROR,*)
C#### Subroutine: NORMALISEGENERAL
C###  Description:
C###    NORMALISEGENERAL performs the same analysis as NORMALISEFLOWS
C###    but for any general field nj specified on the command line.
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn' 
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn' 
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung_nej00.cmn' 
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),nj,NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

      INTEGER nb,ne,ne0,noelem,np1,np2,nv0,nv1,nv2

      CALL ENTERS('NORMALISEGENERAL',*9999)

C  Calculate the flow proportions
      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        ne0=NXI(-1,1,ne) !parent element
        nb=NBJ(1,ne)
        np1=NPNE(1,nb,ne) !start node
        np2=NPNE(2,nb,ne) !end node
        nv1=NVJE(1,nb,nj,ne)
        nv2=NVJE(2,nb,nj,ne)
        IF(ne0.NE.0)THEN
          nv0=NVJE(2,nb,nj,ne0) !version at node in parent
          XP(1,nv1,nj,np1)=XP(1,nv1,nj,np1)/XP(1,1,nj,1)
          XP(1,nv2,nj,np2)=XP(1,nv2,nj,np2)/XP(1,1,nj,1)
        ELSE
          XP(1,nv1,nj,np1)=1.d0
          XP(1,nv2,nj,np2)=1.d0
        ENDIF
      ENDDO
      
      CALL EXITS('NORMALISEGENERAL')
      RETURN
 9999 CALL ERRORS('NORMALISEGENERAL',ERROR)
      CALL EXITS('NORMALISEGENERAL')
      RETURN 1
      END
