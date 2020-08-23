      SUBROUTINE CALC_VORO_FACE2D(cnp,NBJ,NENP,np,
     '  NPNE,nr,VC,XNFV,XP,ZA,ADD,ERROR,*)

C#### Subroutine: CALC_VORO_FACE2D
C###  Description:
C###    Calculates the voronoi face area that is adjoined by
C###    nodes np and cnp. This area is also used to calculate
C###    the additional cell volume added to voronoi cell nvc.

      IMPLICIT NONE
      REAL*8 SMALL
      PARAMETER (SMALL=1.0D-12)
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
      INCLUDE 'voro01.inc'
!     Parameter list
      INTEGER cnp,NBJ(NJM,NEM),NENP(NPM,0:NEPM,0:NRM),np,
     '  NPNE(NNM,NBFM,NEM),nr
      REAL*8 VC,XNFV(-(NJM+1):NJM),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM)
      LOGICAL ADD
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER CNE(2),m,nj,njj,ne,nn,noelem,nb
      REAL*8 HEIGHT,DIST

      CALL ENTERS('CALC_VORO_FACE2D',*9999)

C     ..First compute the two elements common to nodes np and cnp
      m=1
      noelem=1
      DO WHILE(m.LE.2.AND.noelem.LE.NENP(np,0,nr))
        ne=NENP(np,noelem,nr)
        nb=NBJ(1,ne)
        DO nn=1,NNT(nb)
          IF(NPNE(nn,nb,ne).EQ.cnp) THEN
            CNE(m)=ne
            m=m+1
          ENDIF
        ENDDO
        noelem=noelem+1
      ENDDO

C     ..Calculate the unit normal, inverse distance between cnp and np
      XNFV(IDIST)=0.d0
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        nj=NJ_LOC(NJL_GEOM,njj,nr)
        XNFV(njj)=XP(1,1,nj,cnp)-XP(1,1,nj,np)
        XNFV(IDIST)=XNFV(IDIST)+XNFV(njj)**2
      ENDDO

      DIST=DSQRT(XNFV(IDIST))

C     ..Height = height of triangle for volume contribution
C              = half of distance between the two nodes
      HEIGHT=0.5d0*DIST

      XNFV(IDIST)=1.d0/DIST
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        XNFV(njj)=XNFV(njj)*XNFV(IDIST)
      ENDDO

C     ..Now have vertices and unit normal, can calculate the face area,
C       and the volume contribution
C       Face area = distance between two elements' circumcentres
      XNFV(FAREA)=0.d0
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        XNFV(FAREA)=XNFV(FAREA)+
     '    (ZA(1,njj,1,CNE(1))-ZA(1,njj,1,CNE(2)))**2
        XNFV(CENTINDX(njj))=
     '    0.5d0*(ZA(1,njj,1,CNE(1))+ZA(1,njj,1,CNE(2)))
      ENDDO
      XNFV(FAREA)=DSQRT(XNFV(FAREA))

      VC=VC+0.5d0*XNFV(FAREA)*HEIGHT

C     ..Only add to adjacency structure if area is large enough
      IF(DABS(XNFV(FAREA)).LT.SMALL) THEN
        ADD=.FALSE.
      ELSE
        ADD=.TRUE.
      ENDIF

      CALL EXITS('CALC_VORO_FACE2D')
      RETURN
 9999 CALL ERRORS('CALC_VORO_FACE2D',ERROR)
      CALL EXITS('CALC_VORO_FACE2D')
      RETURN 1
      END


