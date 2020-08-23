      SUBROUTINE CLOSEST_SURFACE_TO_NODE(IBT,IDO,INP,NBJ,ne_found,
     '  NEELEM,NKJE,np,NPF,NPNE,nr,NVJE,NXI,
     '  SE,X,XA,XE,XI,XP,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M),ne_found,
     &  NKJE(NKM,NNM,NJM,NEM),np,NPF(9,NFM),NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),X(3),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XI(3),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LDnp,IT,ITMAX,nb,nd,ne,neadj,nelast,neold,ni,nj,no_ne,
     &  noelem
      REAL*8 SQND,SQnp,TEMP
      LOGICAL FOUND
!     Functions
      REAL*8 PXI

      PARAMETER(ITMAX=20)

      CALL ENTERS('CLOSEST_SURFACE_TO_NODE',*9999)

      DO ni=1,3 !compare the distances to the centres of each element
        XI(ni)=0.5D0
      ENDDO
      LDnp=0
      SQnp=0.0D0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     &    nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,
     &    *9999)
        SQND=0.0d0
        DO nj=1,NJT
          nb=NBJ(nj,ne)
          TEMP=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '      XE(1,nj))-X(nj)
          SQND=SQND+TEMP*TEMP
        ENDDO
        IF(LDnp.EQ.0.OR.SQND.LT.SQnp) THEN
          LDnp=ne
          SQnp=SQND
        ENDIF
      ENDDO !no_ne
      nelast=0
      ne=LDnp
      DO ni=1,NJT
        XI(ni)=0.5d0
      ENDDO
      FOUND=.FALSE.
      IF(ne.NE.0) THEN
        IT=0
        DO WHILE(.NOT.FOUND.AND.IT.LT.ITMAX)
          IT=IT+1
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '      nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,
     &      *9999)
          FOUND=.TRUE. !find nearest point in element
          CALL PROJ_ORTHOG(IBT,IDO,INP,NBJ(1,ne),SQND,
     '      XE,XI,X,FOUND,ERROR,*9999)
          neold=nelast
          nelast=ne
          DO ni=1,NJT
            IF(XI(ni).EQ.0.0d0) THEN
              neadj=NXI(-ni,1,ne)
              IF(neadj.GT.0) THEN
                ne=neadj
                XI(ni)=1.0d0
              ENDIF
            ELSE IF(XI(ni).EQ.1.0d0) THEN
              neadj=NXI(ni,1,ne)
              IF(neadj.GT.0) THEN
                ne=neadj
                XI(ni)=0.0d0
              ENDIF
            ENDIF !in bounds
          ENDDO
          FOUND=ne.EQ.neold.OR.ne.EQ.nelast
        ENDDO ! until found
      ENDIF
      ne_found=ne

      CALL EXITS('CLOSEST_SURFACE_TO_NODE')
      RETURN
 9999 CALL ERRORS('CLOSEST_SURFACE_TO_NODE',ERROR)
      CALL EXITS('CLOSEST_SURFACE_TO_NODE')
      RETURN 1
      END



