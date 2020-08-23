      SUBROUTINE CAP_PROJ2(NELIST,NPNE_ALV,centre,XP,ERROR,*)

C#### Subroutine: CAP_PROJ2
C###  Description:
C###    CAP_PROJ2 projects capillary mesh nodes back out to sphere
C###    surface so can join adjacent alveolar meshes together.

C*** Created by KSB, 28th May, 2002.

C*** It is not necessary to carry out this step for the 1st
C*** alveolus as there are no shared faces for the 1st one.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER NELIST(0:NEM),NPNE_ALV(0:NP_NE,NE_R_M)
      REAL*8 centre(NJT),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ne,nj,noelem,nonode,np
      REAL*8 length,radius,XP_A(3)

      CALL ENTERS('CAP_PROJ2',*9999)

      radius=1.d0 !unit sphere
      DO noelem=1,NELIST(0) !for each alveolar element/face
        ne=NELIST(noelem)
          DO nonode=1,NPNE_ALV(0,ne)
C... loops over all capillary nodes on the alveolar face
            np=NPNE_ALV(nonode,ne) !node #
C... project this node from alveolus onto surface of unit sphere
            length=0.d0
            DO nj=1,NJT !vector a (from np to centre)
              XP_A(nj)=XP(1,1,nj,np)-centre(nj)
              length=length+XP_A(nj)**2.d0
            ENDDO
            length=DSQRT(length) !length of vector a
            DO nj=1,NJT
              XP(1,1,nj,np)=centre(nj)+(radius/length)*XP_A(nj) !co-ords
            ENDDO !nj
          ENDDO !nonode
      ENDDO !noelem

      CALL EXITS('CAP_PROJ2')
      RETURN
 9999 CALL ERRORS('CAP_PROJ2',ERROR)
      CALL EXITS('CAP_PROJ2')
      RETURN 1
      END


