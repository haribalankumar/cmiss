      REAL*8 FUNCTION CHMESH_ARCLENGTH(IBT,IDO,INP,nb,ne,NPNE,SE,XI,XP)

C#### Function: CHMESH_ARCLENGTH
C###  Type: REAL*8
C###  Description:
C###    CHMESH_ARCLENGTH is one of J Crocombe's (temporary) change
C###    mesh functions.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  nb,ne,NPNE(NNM,NBFM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XI(3),XP(NKM,NVM,NJM,NPM)
!     Local Variables
      INTEGER ni,nj,nk,nn,ns,NU1(3)
      REAL*8 DERIV,PSI1,TEMP

      DATA NU1/2,4,7/

      ni=1
      TEMP = 0.0d0
      CHMESH_ARCLENGTH = 0.0d0

      DO nj=1,NJT
        DERIV=0.d0
        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            ns=(nn-1)*4+nk
            DERIV=DERIV+PSI1(IBT,IDO,INP,nb,NU1(ni),nk,nn,XI)*
     '        XP(nk,1,nj,NPNE(nn,nb,ne))*SE(ns,nb,ne)
          ENDDO !nk
        ENDDO !nn
        TEMP = TEMP + DERIV**2
      END DO !nj
      CHMESH_ARCLENGTH= DSQRT(TEMP)

      RETURN
      END


