      SUBROUTINE CHECK_POINT_INTERNAL(IBT,IDO,INP,NBJ,NDLIST,NEELEM,NEP,
     &  NHOST,NKJE,np,np1,NPF,NPNE,nr,NVJE,NXI,SE,XA,XE,XIP,XP,X,ZD,
     &  INTERNAL,UPPER,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  NBJ(NJM,NEM),NDLIST(0:NDM),NEELEM(0:NE_R_M,0:NRM),NEP(NPM),
     &  NHOST(NP_R_M),NKJE(NKM,NNM,NJM,NEM),np,np1,NPF(9,NFM),
     &  NPNE(NNM,NBFM,NEM),nr,NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     &  XIP(NIM,NPM),XP(NKM,NVM,NJM,NPM),X(3),ZD(NJM,NDM)
      LOGICAL INTERNAL,UPPER
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ne,nj,noelem,IT
      REAL*8 XI(3),Z_LEN
      LOGICAL FOUND,INELEM

      CALL ENTERS('CHECK_POINT_INTERNAL',*9999)

      FOUND=.FALSE.
      INTERNAL=.TRUE.

      IF(NHOST(np1).GT.0)THEN
        ne=NHOST(np1)
        DO nj=1,NJT
          XI(nj)=0.5d0
        ENDDO
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &    NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &    X,FOUND,INELEM,ERROR,*9999)
        IF(FOUND) NHOST(np)=ne
      ENDIF

      IF(.NOT.FOUND.AND.NXI(-1,0,NHOST(np1)).NE.0)THEN !check each neighbouring element
        ne=NXI(-1,1,NHOST(np1))
        DO nj=1,NJT
          XI(nj)=0.5d0
        ENDDO
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &    NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &    X,FOUND,INELEM,ERROR,*9999)
        IF(FOUND) NHOST(np)=ne
      ENDIF

      IF(.NOT.FOUND.AND.NXI(1,0,NHOST(np1)).NE.0)THEN !check each neighbouring element
        ne=NXI(1,1,NHOST(np1))
        DO nj=1,NJT
          XI(nj)=0.5d0
        ENDDO
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &    NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &    X,FOUND,INELEM,ERROR,*9999)
        IF(FOUND) NHOST(np)=ne
      ENDIF

      IF(.NOT.FOUND.AND.NXI(-2,0,NHOST(np1)).NE.0)THEN !check each neighbouring element
        ne=NXI(-2,1,NHOST(np1))
        DO nj=1,NJT
          XI(nj)=0.5d0
        ENDDO
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &    NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &    X,FOUND,INELEM,ERROR,*9999)
        IF(FOUND) NHOST(np)=ne
      ENDIF

      IF(.NOT.FOUND.AND.NXI(2,0,NHOST(np1)).NE.0)THEN !check each neighbouring element
        ne=NXI(2,1,NHOST(np1))
        DO nj=1,NJT
          XI(nj)=0.5d0
        ENDDO
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &    NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &    X,FOUND,INELEM,ERROR,*9999)
        IF(FOUND) NHOST(np)=ne
      ENDIF

      IF(.NOT.FOUND.AND.NXI(-3,0,NHOST(np1)).NE.0)THEN !check each neighbouring element
        ne=NXI(-3,1,NHOST(np1))
        DO nj=1,NJT
          XI(nj)=0.5d0
        ENDDO
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &    NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &    X,FOUND,INELEM,ERROR,*9999)
        IF(FOUND) NHOST(np)=ne
      ENDIF

      IF(.NOT.FOUND.AND.NXI(3,0,NHOST(np1)).NE.0)THEN !check each neighbouring element
        ne=NXI(3,1,NHOST(np1))
        DO nj=1,NJT
          XI(nj)=0.5d0
        ENDDO
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &    NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        INELEM=.FALSE.
        IT=0
        CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &    X,FOUND,INELEM,ERROR,*9999)
        IF(FOUND) NHOST(np)=ne
      ENDIF

      IF(.NOT.FOUND)THEN !if still not found, check all
        noelem=1
        FOUND=.FALSE.
        INTERNAL=.TRUE.
        DO WHILE(.NOT.FOUND)
          ne=NEELEM(noelem,nr)
          DO nj=1,NJT
            XI(nj)=0.5d0
          ENDDO
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     &      NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          INELEM=.FALSE.
          IT=0
          CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,LOOSE_TOL,XE,XI,
     &      X,FOUND,INELEM,ERROR,*9999)
          IF(.NOT.FOUND)THEN
            noelem=noelem+1
            IF(noelem.GT.NEELEM(0,nr))THEN
              FOUND=.TRUE.
              INTERNAL=.FALSE.
              NHOST(np)=0
            ENDIF
          ELSE
            NHOST(np)=ne
          ENDIF
        ENDDO !while (not found)
      ENDIF

      IF(INTERNAL)THEN
        DO nj=1,NJT
          XIP(nj,np)=XI(nj)
        ENDDO
        NEP(np)=ne

c        IF(.NOT.UPPER)THEN
c          CALL CHECK_CROSSING_FISSURE(NDLIST,X,ZD,INTERNAL,ERROR,*9999)
c        ENDIF
        
      ELSE
        IF(UPPER)THEN
          DO nj=1,NJT
            XIP(nj,np)=0.d0
          ENDDO
          NEP(np)=0
        ELSE
          CALL CLOSEST_SURFACE_TO_NODE(IBT,IDO,INP,NBJ,
     '      ne,NEELEM(0,nr),NKJE,np,NPF,NPNE,nr,NVJE,NXI,
     '      SE,X,XA,XE,XI,XP,ERROR,*9999)
          DO nj=1,NJT
            XIP(nj,np)=XI(nj)
          ENDDO
          NEP(np)=ne
c        write(*,*) 'check point',np,ne,XI(1),XI(2),XI(3)
        ENDIF
      ENDIF
      
      CALL EXITS('CHECK_POINT_INTERNAL')
      RETURN
 9999 CALL ERRORS('CHECK_POINT_INTERNAL',ERROR)
      CALL EXITS('CHECK_POINT_INTERNAL')
      RETURN 1
      END

