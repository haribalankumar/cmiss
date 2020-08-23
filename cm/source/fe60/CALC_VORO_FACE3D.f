      SUBROUTINE CALC_VORO_FACE3D(cnp,NBJ,NENFVC,NENP,nfv,np,NPNE,nr,
     '  NXI,VC,XNFV,XP,ZA,ADD,ERROR,*)

C#### Subroutine: CALC_VORO_FACE3D
C###  Description:
C###    Calculates the voronoi face area, volume, dist and normal

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
      INCLUDE 'voro01.inc'
      REAL*8 ONETHIRD,ONESIXTH,SMALL
      PARAMETER (ONETHIRD=0.33333333333333333d0,
     '  ONESIXTH=0.166666666666666d0,SMALL=1.0d-12)
!     Parameter list
      INTEGER cnp,NBJ(NJM,NEM),NENFVC(0:NFVCM,NFVM),
     '  NENP(NPM,0:NEPM,0:NRM),nfv,np,NPNE(NNM,NBFM,NEM),nr,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 VC,XNFV(-(NJM+1):NJM),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM)
      LOGICAL ADD
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER STARTELEM,ELEM,m,NEXTELEM,nj,njj,np1,
     '  NONAXNODE(2),NONAXLOCL(2),ne,nn,noelem,
     '  VERTTOTAL,VERT,nb
      REAL*8 CENTROID(3),DDOT,DOTPROD,FACEAREA,PROJECTION(3),
     '  VERTICES(3,60),XPROD(3)
      LOGICAL FIN,FIN2

      CALL ENTERS('CALC_VORO_FACE3D',*9999)

C     ..First compute an element common to np and cnp
      noelem=1
      FIN=.FALSE.
      DO WHILE(.NOT.FIN.AND.noelem.LE.NENP(np,0,nr))
        ne=NENP(np,noelem,nr)
        nb=NBJ(1,ne)
        nn=1
        DO WHILE(.NOT.FIN.AND.nn.LE.NNT(nb))
          IF(NPNE(nn,nb,ne).EQ.cnp) THEN
            STARTELEM=ne
            FIN=.TRUE.
          ELSE
            nn=nn+1
          ENDIF
        ENDDO
        noelem=noelem+1
      ENDDO

      NENFVC(0,nfv)=1 !# of ne contributing to face
      NENFVC(1,nfv)=STARTELEM

C     ..Now find non-axial nodes of that element
      m=1
      nn=1
      nb=NBJ(1,STARTELEM)
      DO WHILE(m.LE.2.AND.nn.LE.NNT(nb))
        np1=NPNE(nn,nb,STARTELEM)
        IF(np1.NE.cnp.AND.np1.NE.np) THEN
          NONAXNODE(m)=np1
          NONAXLOCL(m)=nn
          m=m+1
        ENDIF
        nn=nn+1
      ENDDO

C     ..Now compute the elements around the axis in circular order
      FIN=.FALSE.
      VERT=1
      ELEM=STARTELEM
      DO WHILE(.NOT.FIN)

C       ..Compute vertex from node to tetrahedral centre
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          VERTICES(njj,VERT)=ZA(1,njj,1,ELEM)-XP(1,1,nj,np)
        ENDDO

C       ..Get the next element
        NEXTELEM=NXI(0,NONAXLOCL(1),ELEM)
        IF(NEXTELEM.EQ.0) THEN
          WRITE(ERROR,'(''>>Found non-axial node opposite'//
     '      ' nothing: ne = '',I6,'' Local node = '',I1,'' np'//
     '      ' = '',I5,'' cnp = '',I5)') ELEM,NONAXLOCL(1),np,cnp
          GOTO 9999
        ENDIF

        IF(NEXTELEM.EQ.STARTELEM) THEN
          FIN=.TRUE.
        ELSE
          NONAXNODE(1)=NONAXNODE(2)
          nb=NBJ(1,NEXTELEM)
C         ..Record elements that contribute vertices to face nfv
          NENFVC(0,nfv)=NENFVC(0,nfv)+1
          NENFVC(NENFVC(0,nfv),nfv)=NEXTELEM

C         ..First find the local node of NONAXNODE(1)
          nn=1
          FIN2=.FALSE.
          DO WHILE(.NOT.FIN2.AND.nn.LE.NNT(nb))
            np1=NPNE(nn,nb,NEXTELEM)
            IF(np1.EQ.NONAXNODE(1)) THEN
              FIN2=.TRUE.
              NONAXLOCL(1)=nn
            ENDIF
            nn=nn+1
          ENDDO

C         ..Now compute the next axis
          nn=1
          FIN2=.FALSE.
          DO WHILE(.NOT.FIN2.AND.nn.LE.NNT(nb))
            np1=NPNE(nn,nb,NEXTELEM)
            IF(np1.NE.cnp.AND.np1.NE.np.AND.np1.NE.NONAXNODE(1)) THEN
              FIN2=.TRUE.
              NONAXNODE(2)=np1
              NONAXLOCL(2)=nn
            ENDIF
            nn=nn+1
          ENDDO
          ELEM=NEXTELEM
          VERT=VERT+1
        ENDIF
      ENDDO
      VERTTOTAL=VERT

C     ..Make the end VERTEX the start VERTEX to complete the loop
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        VERTICES(njj,VERT+1)=VERTICES(njj,1)
      ENDDO

C     ..Calculate the unit normal, distance
      XNFV(IDIST)=0.d0
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        nj=NJ_LOC(NJL_GEOM,njj,nr)
        XNFV(njj)=XP(1,1,nj,cnp)-XP(1,1,nj,np)
        XNFV(IDIST)=XNFV(IDIST)+XNFV(njj)**2
      ENDDO
      XNFV(IDIST)=1.d0/DSQRT(XNFV(IDIST))
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        XNFV(njj)=XNFV(njj)*XNFV(IDIST)
      ENDDO

C     ..Now have vertices and unit normal, can calculate the face area,
C       and the volume contribution
      FACEAREA=0.d0
      DOTPROD=DDOT(NJ_LOC(NJL_GEOM,0,nr),VERTICES(1,1),1,XNFV(1),1)
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        CENTROID(njj)=0.d0
        PROJECTION(njj)=DOTPROD*XNFV(njj)
      ENDDO
      DO VERT=1,VERTTOTAL
        CALL CROSS(VERTICES(1,VERT),VERTICES(1,VERT+1),XPROD)
        DOTPROD=DDOT(NJ_LOC(NJL_GEOM,0,nr),XPROD,1,XNFV(1),1)
        FACEAREA=FACEAREA+DOTPROD
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          CENTROID(njj)=CENTROID(njj)+DOTPROD*ONETHIRD
     '      *(VERTICES(njj,VERT)+VERTICES(njj,VERT+1)+PROJECTION(njj))
        ENDDO
      ENDDO
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        XNFV(CENTINDX(njj))=CENTROID(njj)
      ENDDO
      
      IF(DABS(FACEAREA).LT.SMALL) THEN
        ADD=.FALSE.
      ELSE
        ADD=.TRUE.
        XNFV(FAREA)=0.5d0*DABS(FACEAREA)
      ENDIF
      IF(FACEAREA.LT.0.d0) THEN
        VC=VC-ONESIXTH*DDOT(NJ_LOC(NJL_GEOM,0,nr),CENTROID,1,XNFV(1),1)
      ELSE
        VC=VC+ONESIXTH*DDOT(NJ_LOC(NJL_GEOM,0,nr),CENTROID,1,XNFV(1),1)
      ENDIF

C 9998 CALL EXITS('CALC_VORO_FACE3D')
C      RETURN
      CALL EXITS('CALC_VORO_FACE3D')
      RETURN
 9999 CALL ERRORS('CALC_VORO_FACE3D',ERROR)
      CALL EXITS('CALC_VORO_FACE3D')
      RETURN 1
      END


