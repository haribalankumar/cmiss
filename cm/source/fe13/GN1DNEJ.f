      SUBROUTINE GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,np,np_start,nr,NRE,NVJE,
     '  NVJP,SE,ERROR,*)

C#### Subroutine: GN1DNEJ
C###  Description:
C###    GN1DNEJ sets up NRE, NBJ, SE, NVJE, NKJE, NKJ, NVJP, NKJ, and
C###    NVJP for 1D linear elements. This routine is used to set up the
C###    arrays as each individual element in a 1D network is created.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'

!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),ne,NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),np,
     '  np_start,nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM)
      CHARACTER ERROR*(*)
      INTEGER nj,nj1,nj2,nk,nn,ns

      CALL ENTERS('GN1DNEJ',*9999)

      NRE(ne)=nr
      DO nj1=1,3 !NJL_GEOM,NJL_FIBR,NJL_FIEL
        DO nj2=1,NJ_LOC(nj1,0,nr)
          nj=NJ_LOC(nj1,nj2,nr)
          NBJ(nj,ne)=nb
        ENDDO !nj
      ENDDO !nj1
      DO ns=1,NST(nb)+NAT(nb)
        SE(ns,nb,ne)=1.0d0
      ENDDO !ns
      DO nn=1,NNT(nb)
        DO nj1=1,3
          DO nj2=1,NJ_LOC(nj1,0,nr)
            nj=NJ_LOC(nj1,nj2,nr)
            NVJE(nn,nb,nj,ne)=1
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=1
            ENDDO !nk
          ENDDO !nj
        ENDDO !nj1
      ENDDO !nn
      DO nj1=1,3
        DO nj2=1,NJ_LOC(nj1,0,nr)
          nj=NJ_LOC(nj1,nj2,nr)
          NKJ(nj,np)=1
          NVJP(nj,np)=1
          NKJ(nj,np_start)=1
          NVJP(nj,np_start)=1
        ENDDO !nj2
      ENDDO !nj1

      CALL EXITS('GN1DNEJ')
      RETURN
 9999 CALL ERRORS('GN1DNEJ',ERROR)
      CALL EXITS('GN1DNEJ')
      RETURN 1
      END

