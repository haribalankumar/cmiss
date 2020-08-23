      SUBROUTINE MESH_COFM(LD,nen,COFM,ZD,ERROR,*)

C#### Subroutine: MESH_COFM
C###  Description:
C###    MESH_COFM calculates the centre of mass of a collection of
C###    random points by averaging their coordinates.
C***  Created by Merryn Howatson Tawhai, February 1997

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INTEGER LD(NDM)
      REAL*8 COFM(3),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER DAT,nd,nj,nen,nsp

      CALL ENTERS('MESH_COFM',*9999)

      DAT=0
      DO nj=1,NJT
        COFM(nj)=0.0d0
      ENDDO !nj
      DO nd=1,NDT
        nsp=LD(nd) !the space # for the nd-th data point
        IF(nsp.EQ.nen)THEN
          DAT=DAT+1
          DO nj=1,NJT
            COFM(nj)=COFM(nj)+ZD(nj,nd)
          ENDDO !nj
        ENDIF
      ENDDO !nd
      IF(DAT.NE.0)THEN
        DO nj=1,NJT
          COFM(nj)=COFM(nj)/DAT !centre of mass
        ENDDO !nj
      ENDIF

      CALL EXITS('MESH_COFM')
      RETURN
 9999 CALL ERRORS('MESH_COFM',ERROR)
      CALL EXITS('MESH_COFM')
      RETURN 1
      END


