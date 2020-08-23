      SUBROUTINE LINCAL_1D(NBJ,NEELEM,NKJE,nr,NRE,NVJE,SE,ERROR,*)

C#### Subroutine: LINCAL_1D
C###  Description:
C###    LINCAL_1D calculates global line parameters for 1D mesh.

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb,noelem,ne,nj,nk,nn,ns

      CALL ENTERS('LINCAL_1D',*9999)

c Unit scale factors; avoids time-consuming looping  in LINCAL
      nb=NBJ(1,NEELEM(1,nr))
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        NRE(ne)=nr
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          NBJ(nj,ne)=nb
        ENDDO !nj
        DO ns=1,NST(nb)+NAT(nb)
          SE(ns,nb,ne)=1.0d0
        ENDDO !ns
        DO nn=1,NNT(nb)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            NVJE(nn,nb,nj,ne)=1
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=1
            ENDDO !nk
          ENDDO !nj
        ENDDO !nn
      ENDDO

      CALL_LINE=.TRUE.

      CALL EXITS('LINCAL_1D')
      RETURN
 9999 CALL ERRORS('LINCAL_1D',ERROR)
      CALL EXITS('LINCAL_1D')
      RETURN 1
      END


