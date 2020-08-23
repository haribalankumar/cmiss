      SUBROUTINE CHMESH_CALC_ARCLENGTHS(IBT,IDO,INP,NBJ,NEELEM,NGAP,
     '  NPNE,nr,DL,SE,TOTAL,WG,XIG,XP,Z,ERROR,*)

C#### Subroutine: CHMESH_CALC_ARCLENGTHS
C###  Description:
C###    CHMESH_CALC_ARCLENGTHS is J.Crocombes calc_arclengths
C###    subroutine.  Calculates the circumference of the mesh at a
C###    given Z value.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NGAP(NIM,NBM),
     '  NPNE(NNM,NBFM,NEM),nr
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM),TOTAL,WG(NGM,NBM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),Z
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 CHMESH_ARCLENGTH,INTEGRAL,XI(3),XI2PREV
      INTEGER nb,ne,ng,noelem,noelem2,none,max,TEMP,W,znb,zne

      CALL ENTERS('CHMESH_CALC_ARCLENGTHS',*9999)
      TOTAL=0.0d0
      XI(3)=0.d0
      none=0
      !This subroutine calculates the distance around the mesh at any
      !given height or z value.
      IF(NJT.GT.2) THEN
        !find element containing z value.
        DO noelem=1,NEELEM(0,nr)
          zne=0
          ne=NEELEM(noelem,nr)
          nb=NBJ(3,ne)
          IF(XP(1,1,3,NPNE(3,nb,ne)).LT.Z.AND.XP(1,1,3,NPNE(1,nb,ne))
     '       .GE.Z) THEN
            none=1
            zne=ne
            znb=nb
          ELSEIF(XP(1,1,3,NPNE(3,nb,ne)).GE.Z.AND.
     '       XP(1,1,3,NPNE(1,nb,ne)).LT.Z) THEN
            none=1
            zne=ne
            znb=nb
          ENDIF
          IF(zne.ne.0) THEN
            INTEGRAL=0.0d0
            XI(2)=0.5d0
            DO ng = 1,NGAP(2,nb)
              XI(1) = XIG(1,ng,1)
                XI2PREV=XI(2)
                CALL CHMESH_FIND_ZETA(IBT,IDO,INP,max,znb,zne,NPNE,SE,
     '            XP,Z,XI(1),XI(2),ERROR,*9999)
                IF(max.EQ.50.AND.XI2PREV.LT.0.5d0) THEN
                  TEMP = NPNE(1,znb,zne)
                  noelem2=1
                  zne= NEELEM(noelem2,nr)
                  DO WHILE (NPNE(3,znb,zne).NE.TEMP)
                    noelem2=noelem2+1
                    zne=NEELEM(noelem2,nr)
                    znb=NBJ(2,ne)
                  ENDDO
                  CALL CHMESH_FIND_ZETA(IBT,IDO,INP,max,znb,zne,NPNE,SE,
     '            XP,Z,XI(1),XI(2),ERROR,*9999)
                ELSEIF(max.EQ.50.AND.XI2PREV.GE.0.5d0) THEN
                  TEMP = NPNE(3,znb,zne)
                  noelem2=1
                  zne= NEELEM(noelem2,nr)
                  DO WHILE (NPNE(1,znb,zne).NE.TEMP)
                    noelem2=noelem2+1
                    zne=NEELEM(noelem2,nr)
                    znb=NBJ(2,ne)
                  ENDDO
                  CALL CHMESH_FIND_ZETA(IBT,IDO,INP,max,znb,zne,NPNE,SE,
     '            XP,Z,XI(1),XI(2),ERROR,*9999)
                ENDIF
              CALL ASSERT(max.NE.50,'>>Circumference not measured',
     '         ERROR,*9999)
              W=(ng-1)*3
              INTEGRAL = INTEGRAL + (WG(W+1,1)+WG(W+2,1)+WG(W+3,1))*
     '          CHMESH_ARCLENGTH(IBT,IDO,INP,znb,zne,NPNE,SE,XI,XP)
            ENDDO
            TOTAL=TOTAL+INTEGRAL
           ENDIF
        ENDDO !noelem
      ELSEIF(NJT.EQ.2) THEN
        CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
        zne=ne
        znb=NBJ(2,zne)
        INTEGRAL=DL(3,ne)
        TOTAL=TOTAL+INTEGRAL
      ENDIF
      CALL ASSERT(none.GT.0,'>>Invalid measurement request',
     '    ERROR,*9999)

      CALL EXITS('CHMESH_CALC_ARCLENGTHS')
      RETURN
 9999 CALL ERRORS('CHMESH_CALC_ARCLENGTHS',ERROR)
      CALL EXITS('CHMESH_CALC_ARCLENGTHS')
      RETURN 1
      END


