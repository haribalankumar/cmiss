      SUBROUTINE CHMESH_FIND_ZETA(IBT,IDO,INP,max,nb,ne,NPNE,SE,
     '  XP,Z,XI1,ZETA2,ERROR,*)

C#### Subroutine: CHMESH_FIND_ZETA
C###  Description:
C###    CHMESH_FIND_ZETA is J.Crocombes find_zeta subroutine.  Uses
C###    the Newton Raphson scheme and bicubic interpolation with
C###    constant zeta 1 to find the appropriate zeta 2 value.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  max,nb,ne,NPNE(NNM,NBFM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),Z,ZETA2
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nk,nn,ns
      REAL*8 PSI1,TEMP,FN,DFN,XI1,XI(3),ZERO_TOLERANCE

      CALL ENTERS('CHMESH_FIND_ZETA',*9999)
      ZERO_TOLERANCE=0.1d-6

      XI(2)=0.2d0
      TEMP=1.0d0
      max=0
      DO WHILE(DABS(XI(2)-TEMP).GT.ZERO_TOLERANCE.AND.max.LT.50)
        TEMP=XI(2)
        FN=0.d0
        DFN=0.d0
        XI(1)=XI1
        XI(3)=0.d0
        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            ns=(nn-1)*4+nk
            FN=FN+PSI1(IBT,IDO,INP,nb,1,nk,nn,XI)*
     '        XP(nk,1,3,NPNE(nn,nb,ne))*SE(ns,nb,ne)
            DFN=DFN+PSI1(IBT,IDO,INP,nb,4,nk,nn,XI)*
     '        XP(nk,1,3,NPNE(nn,nb,ne))*SE(ns,nb,ne)
          ENDDO !nk
        ENDDO !nn
        FN=FN-Z
        XI(2)=XI(2)-(FN/DFN)
        max=max+1
      ENDDO
      ZETA2=XI(2)

      CALL EXITS('CHMESH_FIND_ZETA')
      RETURN
 9999 CALL ERRORS('CHMESH_FIND_ZETA',ERROR)
      CALL EXITS('CHMESH_FIND_ZETA')
      RETURN 1
      END


