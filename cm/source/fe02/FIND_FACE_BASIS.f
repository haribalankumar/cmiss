      SUBROUTINE FIND_FACE_BASIS(IBT,nb,nbf,nef,NNF,ERROR,*)

C#### Subroutine:  FIND_FACE_BASIS
C###  Description:
C###    FIND_FACE_BASIS returns the face basis number nbf of the face
C###    that corresponds to local face nef of element ne which is
C###    interpolated with element basis nb. If no such basis exists
C###    FIND_FACE_BASIS returns 0.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),nb,nbf,nef,NNF(0:17,6,NBFM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER II(5),ni1,ni2,ni3,nicollapse,NUMCOLLAPSED
      LOGICAL FOUND,SECTOR

      DATA II/1,2,3,1,2/

      CALL ENTERS('FIND_FACE_BASIS',*9999)

      CALL ASSERT(nb.NE.0,'>>Element basis is zero',ERROR,*9999)
      CALL ASSERT(nef.NE.0,'>>Local element face is zero',ERROR,*9999)
      ni1=NNF(1,nef,nb)
      CALL ASSERT(ni1.NE.0,'>>Face normal is zero',ERROR,*9999)
      ni2=II(ni1+1)
      ni3=II(ni1+2)
      SECTOR=.FALSE.
      NUMCOLLAPSED=0
      IF(IBT(1,ni2,nb).EQ.5.OR.IBT(1,ni2,nb).EQ.6) THEN
        NUMCOLLAPSED=1
        nicollapse=ni2
        SECTOR=.TRUE.
      ENDIF
      IF(IBT(1,ni3,nb).EQ.5.OR.IBT(1,ni3,nb).EQ.6) THEN
        NUMCOLLAPSED=NUMCOLLAPSED+1
        nicollapse=ni3
        SECTOR=.TRUE.
      ENDIF
C     Try and find a global basis that matches the face basis.
      nbf=0
      FOUND=.FALSE.
      DO WHILE(.NOT.FOUND.AND.nbf.LE.NBT)
        nbf=nbf+1
C KAT 28Jan99: No longer demanding same scale factors as face scale
C              factors come from the element scale factors anyway.
C        IF(NIT(nbf).EQ.2.AND.(NBC(nbf).NE.5.AND.NBC(nbf).NE.6).AND
C     '    .(NKT(0,nbf).LT.2.OR.NNT(nbf).EQ.0.OR.NBI(nbf).EQ.NBI(nb)))
C     '    THEN
        IF(NIT(nbf).EQ.2.AND.(NBC(nbf).NE.5.AND.NBC(nbf).NE.6)) THEN
C         2D basis, not BEM, and correct scale factor if there are derivs.
          IF(SECTOR) THEN
C           Check the sector case. If we are on a face with a face
C           xi direction that has been collapsed (in the global
C           element) but the other xi direction in the face is not
C           the collapse direction then we need a special case.
            IF(IBT(3,nicollapse,nb).EQ.ni1) THEN
              IF(NUMCOLLAPSED.EQ.2) THEN
                IF(IBT(2,ni2,nb).EQ.2) THEN
                  IF(IBT(1,1,nbf).EQ.2.AND.
     '              IBT(2,1,nbf).EQ.1)
     '              FOUND=.TRUE.
                ELSE
                  IF(IBT(1,1,nbf).EQ.1.AND.
     '              IBT(2,1,nbf).EQ.IBT(2,ni2,nb))
     '              FOUND=.TRUE.
                ENDIF
                IF(IBT(2,ni3,nb).EQ.2) THEN
                  IF(IBT(1,1,nbf).EQ.2.AND.IBT(2,1,nbf).EQ.1) THEN
                    FOUND=.TRUE.
                  ELSE
                    FOUND=.FALSE.
                  ENDIF
                ELSE
                  IF(IBT(1,1,nbf).EQ.1.AND.
     '              IBT(2,1,nbf).EQ.IBT(2,ni2,nb)) THEN
                    FOUND=.TRUE.
                  ELSE
                    FOUND=.FALSE.
                  ENDIF
                ENDIF
              ELSE
                IF(nicollapse.EQ.ni2) THEN
                  IF(IBT(1,2,nbf).EQ.IBT(1,ni3,nb).AND.
     '              IBT(2,2,nbf).EQ.IBT(2,ni3,nb)) THEN
                    IF(IBT(2,ni2,nb).EQ.4) THEN
                      IF(IBT(1,1,nbf).EQ.2.AND.IBT(2,1,nbf).EQ.1)
     '                  FOUND=.TRUE.
                    ELSE
                      IF(IBT(1,1,nbf).EQ.1.AND.
     '                  IBT(2,1,nbf).EQ.IBT(2,ni2,nb))
     '                  FOUND=.TRUE.
                    ENDIF
                  ENDIF
                ELSE
                  IF(IBT(1,1,nbf).EQ.IBT(1,ni2,nb).AND.
     '              IBT(2,1,nbf).EQ.IBT(2,ni2,nb))
     '              THEN
                    IF(IBT(2,ni3,nb).EQ.4) THEN
                      IF(IBT(1,2,nbf).EQ.2.AND.IBT(2,2,nbf).EQ.1)
     '                  FOUND=.TRUE.
                    ELSE
                      IF(IBT(1,2,nbf).EQ.1.AND.
     '                  IBT(2,2,nbf).EQ.IBT(2,ni3,nb))
     '                  FOUND=.TRUE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ELSE
              IF(IBT(1,1,nbf).EQ.IBT(1,ni2,nb).AND.
     '          IBT(2,1,nbf).EQ.IBT(2,ni2,nb).AND.
     '          IBT(1,2,nbf).EQ.IBT(1,ni3,nb).AND.
     '          IBT(2,2,nbf).EQ.IBT(2,ni3,nb))
     '          FOUND=.TRUE.
            ENDIF
C         Check the standard case
          ELSE IF((IBT(1,1,nbf).EQ.IBT(1,ni2,nb)).AND.
     '        (IBT(1,2,nbf).EQ.IBT(1,ni3,nb)).AND.
     '        (IBT(2,1,nbf).EQ.IBT(2,ni2,nb)).AND.
     '        (IBT(2,2,nbf).EQ.IBT(2,ni3,nb))) THEN
            FOUND=.TRUE.
          ENDIF
        ENDIF
      ENDDO !nbf
      IF(.NOT.FOUND) nbf=0

      CALL EXITS('FIND_FACE_BASIS')
      RETURN
 9999 CALL ERRORS('FIND_FACE_BASIS',ERROR)
      CALL EXITS('FIND_FACE_BASIS')
      RETURN 1
      END


