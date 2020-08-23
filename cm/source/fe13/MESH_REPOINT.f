      SUBROUTINE MESH_REPOINT(LD,nb,N_ELM,NE_OLD,NE_REACTIVATE,
     &  NPNE,DISTANCE_LIMIT,XP,ZD,FIRST,ERROR,*)

C#### Subroutine: MESH_REPOINT
C###  Description:
C###  MESH_REPOINT reassigns data (seed) points to the closest ending of
C###  branches in the current generation.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

      !Parameter list
      INTEGER LD(NDM),nb,N_ELM,NE_OLD(NE_R_M),NE_REACTIVATE(0:NEM),
     &  NPNE(NNM,NBFM,NEM)
      REAL*8 DISTANCE_LIMIT,XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      LOGICAL FIRST
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER nd,ne,ne_min,noelem,np
      REAL*8 DIST,MIN_DIST
      CHARACTER CHAR*1
!     External functions
      INTEGER IDIGITS

      CALL ENTERS('MESH_REPOINT',*9999)

      NE_REACTIVATE(0)=0
      DO nd=1,NDT
        IF(LD(nd).NE.0)THEN
          MIN_DIST=1.d10
          IF(N_ELM.GT.NE_R_M) THEN
            WRITE(CHAR,'(I1)') IDIGITS(N_ELM)
            WRITE(ERROR,'(''>>Increase NE_R_M to '',I'//CHAR//')')
     '        N_ELM
            GO TO 9999
          ENDIF
          DO noelem=1,N_ELM
            ne=NE_OLD(noelem)
            np=NPNE(2,nb,ne)
            DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np))**2.d0+(ZD(2,nd)
     '        -XP(1,1,2,np))**2.d0+(ZD(3,nd)-XP(1,1,3,np))**2.d0)
            IF(DIST.LT.MIN_DIST)THEN
              ne_min=ne
              MIN_DIST=DIST
            ENDIF
          ENDDO
          IF(FIRST)THEN
            LD(nd)=ne_min
          ELSE
            IF(MIN_DIST.LT.DISTANCE_LIMIT)THEN !keep seed points
              LD(nd)=ne_min
            ELSE
c              LD(nd)=0 !too far from branch ends, so discard
c mht changing this so that seed point not discarded.
c              NE_REACTIVATE(0)=NE_REACTIVATE(0)+1
c              NE_REACTIVATE(NE_REACTIVATE(0))=LD(nd)
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      CALL EXITS('MESH_REPOINT')
      RETURN
 9999 CALL ERRORS('MESH_REPOINT',ERROR)
      CALL EXITS('MESH_REPOINT')
      RETURN 1
      END


