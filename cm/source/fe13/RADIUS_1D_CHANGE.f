      SUBROUTINE RADIUS_1D_CHANGE(NBJ,ne,NPNE,NVJE,radius,XP,ERROR,*)

C#### Subroutine: RADIUS_1D_CHANGE
C###  Description:
C###    RADIUS_1D_CHANGE inputs motion parameters for pulmonary problems in
C###    region nr.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),ne,NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 radius,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nn,npnn(2),nvnn(2)

      CALL ENTERS('RADIUS_1D_CHANGE',*9999)

      nb=NBJ(nj_radius,ne) !basis function for radii
      DO nn=1,NNT(nb) !versions for radii at nodes
        npnn(nn)=NPNE(nn,nb,ne)
        nvnn(nn)=NVJE(nn,nb,nj_radius,ne)
      ENDDO !nn
      XP(1,nvnn(1),nj_radius,npnn(1))=radius
      XP(1,nvnn(2),nj_radius,npnn(2))=radius

      CALL EXITS('RADIUS_1D_CHANGE')
      RETURN
 9999 CALL ERRORS('RADIUS_1D_CHANGE',ERROR)
      CALL EXITS('RADIUS_1D_CHANGE')
      RETURN 1
      END


