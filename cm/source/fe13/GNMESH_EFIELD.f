      SUBROUTINE GNMESH_EFIELD(NEELEM,nr,CONST_EFIELD,XAB,LPM_FIELD,
     &  EFIELD_VAL,ERROR,*)

C#### Subroutine: GNMESH_EFIELD
C###  Description:
C###  GNMESH_EFIELD sets up element field structure for pulmonary meshes. 

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER NEELEM(0:NE_R_M,0:NRM),nr
      REAL*8  CONST_EFIELD,XAB(NORM,NEM)
      LOGICAL LPM_FIELD,EFIELD_VAL
      CHARACTER ERROR*(*)  
! local variables
      INTEGER noelem,ne

      CALL ENTERS('GNMESH_EFIELD',*9999)


      IF(LPM_FIELD.AND.EFIELD_VAL)THEN
      NEJ_LOC(0,nr)=NEJ_LOC(0,nr)+1
      NEJ_LOC(0,0)=NEJ_LOC(0,0)+1
      NEJ_LOC(nej_cap,nr)=1
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            XAB(nej_cap,ne)=CONST_EFIELD
          ENDDO !neelem
      ENDIF

      CALL EXITS('GNMESH_EFIELD')
      RETURN
 9999 CALL ERRORS('GNMESH_EFIELD',ERROR)
      CALL EXITS('GNMESH_EFIELD')
      RETURN 1
      END

