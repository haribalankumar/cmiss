      SUBROUTINE DIST(intscheme,nb,NLL,nnmin,NPNE,nr,NW,DL,
     '  MINDIST,XP,XPFP,ERROR,*)

C#### Subroutine: DIST
C###  Description:
C###    DIST determines the minimum distance (MINDIST) from the
C###    position XPFP and the local nodes of the element ne. The local
C###    node number corresponding to this minimum distance is returned
C###    in nnmin. It also returns the integration scheme to be used
C###    based on this minimum distance.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER intscheme,nb,NLL(12),nnmin,NPNE(NNM,NBFM),nr,NW
      REAL*8 DL(3,NLM),MINDIST,XP(NKM,NVM,NJM,NPM),XPFP(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nae,nbbt,nj,nl,nn
      REAL*8 AV_DIM,DISTANCE,SCALEDIST,SUM,XGC(3),XPC(3),XR_LOCAL(3)

      CALL ENTERS('DIST',*9999)

      DO nj=1,3
        XR_LOCAL(nj)=0.0D0
      ENDDO !nj
      SCALEDIST=-1.0d0
      MINDIST=RMAX
      IF(ITYP10(nr).EQ.1) THEN
        DO nn=1,NNT(nb)
          SUM=0.0d0
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            XR_LOCAL(nj)=XP(1,1,nj,NPNE(nn,nb))-XPFP(nj)
            SUM=SUM+XR_LOCAL(nj)*XR_LOCAL(nj)
          ENDDO !nj
          DISTANCE=DSQRT(SUM)
          IF(DISTANCE.LT.MINDIST) THEN
            MINDIST=DISTANCE
            NNMIN=nn
          ENDIF
        ENDDO !nn
      ELSE
C*** Transform to cartesian to find distance
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          XR_LOCAL(nj)=XPFP(nj)
        ENDDO
        CALL COORD(ITYP10(nr),1,XR_LOCAL,XPC,ERROR,*9999)
        DO nn=1,NNT(nb)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            XR_LOCAL(nj)=XP(1,1,nj,NPNE(nn,nb))
          ENDDO
          CALL COORD(ITYP10(nr),1,XR_LOCAL,XGC,ERROR,*9999)
          SUM=0.0d0
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            XR_LOCAL(nj)=XGC(nj)-XPC(nj)
            SUM=SUM+XR_LOCAL(nj)*XR_LOCAL(nj)
          ENDDO
          DISTANCE=DSQRT(SUM)
          IF(DISTANCE.LT.MINDIST) THEN
            MINDIST=DISTANCE
            NNMIN=nn
          ENDIF
        ENDDO !nn
      ENDIF
      IF(MINDIST.LE.RDELTA) THEN !Node within the element
        intscheme=1 !element splitting scheme
      ELSE
        IF(ADAPINT) THEN
C*** Scale distance by average element dimension.
          IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN !1d element
            nl=NLL(1)
            SCALEDIST=MINDIST/DL(3,nl) !Scale by len of element side
          ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN !2d element
C*** Average side length in xi1 and xi2 direction
            AV_DIM=0.0d0
            DO nae=1,4
              nl=NLL(nae)
              IF(nl.GT.0) AV_DIM=AV_DIM+DL(3,nl)
            ENDDO
            AV_DIM=AV_DIM/4.0d0
            CALL ASSERT(AV_DIM.GE.RDELTA,'>>Zero size element ?',
     '        ERROR,*9999)
            SCALEDIST=MINDIST/AV_DIM
          ENDIF
        ENDIF
        IF(ADAPINT.AND.SCALEDIST.GE.0.05d0.AND.SCALEDIST.LE.1.5d0) THEN
          intscheme=2 !adaptive integration
        ELSE
          nbbt=NBASEF(nb,0)
          IF(MINDIST.LE.DLIM(nbbt-3,NW)) THEN
            intscheme=3 !high order gauss scheme
          ELSE IF(MINDIST.LE.DLIM(nbbt-2,NW)) THEN !Medium scheme
            intscheme=4 !medium order gauss scheme
          ELSE
            intscheme=5 !low order gauss scheme
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('DIST')
      RETURN
9999  CALL ERRORS('DIST',ERROR)
      CALL EXITS('DIST')
      RETURN 1
      END


