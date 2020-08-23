      SUBROUTINE MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '  A_VECTOR,B_VECTOR,C_VECTOR,XE,XG,XI,CALC_XG,ERROR,*)

C#### Subroutine: MAT_VEC_XI
C###  Description:
C###    MAT_VEC_XI calculates direction cosines of material vectors at
C###    XI in element ne.
C###    This routine assumes that XE contains undeformed
C###    element vertex coordinates, and microstructural material angles.

C     MPN rewritten 25-Apr-96

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBJ(NJM),nr
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  XE(NSM,NJM),XG(NJM,NUM),XI(3)
      CHARACTER ERROR*(*)
      LOGICAL CALC_XG
!     Local Variables
      INTEGER ni,NITB
      REAL*8 DXRCXI(3,3)

      CALL ENTERS('MAT_VEC_XI',*9999)

      NITB=NIT(NBJ(1))
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Xi coords: '',3D12.3)') (XI(ni),ni=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C     Interpolate midwall geometric vars XG and derivs wrt Xi
      IF(CALC_XG) THEN
        CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
      ENDIF

C     Compute derivatives of rc coords wrt Xi coords
      CALL DXRCDXI(NITB,nr,DXRCXI,XG,ERROR,*9999)

C     Calculate direction cosines
      CALL MAT_VEC(NITB,nr,A_VECTOR,B_VECTOR,C_VECTOR,DXRCXI,XG,
     '  ERROR,*9999)

C 21/2/97 LC archived section : old MPN 25-Apr-96
CC**** Note: XE must be passed in (calculated in previous call to XPXE).
CC**** A_VECTOR is fibre angle vector (in sheet)
CC**** B_VECTOR is sheet angle vector (in sheet orthog to fibres)
CC**** C_VECTOR   is normal to sheet
CC**** Alfa is fibre angle
CC**** Beta is imbrication angle
CC**** Gama is sheet angle
CC**** Delta is angle between h (wall normal) and v (normal to wall in
CC**** yz-plane)
CC**** Phi is angle between v and r
CC**** GAMA,PHI & DELTA are also returned.
CC**** Note: If Beta>0 A_VECTOR is calculated with Beta terms included.

      CALL EXITS('MAT_VEC_XI')
      RETURN
 9999 CALL ERRORS('MAT_VEC_XI',ERROR)
      CALL EXITS('MAT_VEC_XI')
      RETURN 1
      END


