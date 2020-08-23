      SUBROUTINE TRI_NE_NP(nb,NBJ,NEELEM,NKJ,NKJE,NPI,NPNE,nr,
     '  NRE,NVJE,NVJP,TRIANGLES,SE,VERTICES_XYZ,XP,ERROR,*)

C###  Subroutine: TRI_NE_NP
C###  Description:
C###    TRI_NE_NP sets up nodes and elements for Delaunay triangles,
C###    note that these are only useful for exporting to view.

C***  Created by: Kelly Burrowes, February 2002
C***  Last modified: June 2002

C*** The array NPI(np) is used to give the 1st I position in array
C*** VERTICES_XYZ(I) for np.

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NKJ(NJM,NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NPI(N_VERT_EST),NPNE(NNM,NBFM,NEM),
     '  nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  TRIANGLES(3*MAX_TRIANGLES)
      REAL*8 SE(NSM,NBFM,NEM),VERTICES_XYZ(3*N_VERT_EST),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      INTEGER I,ne,noelem,np,ntri,TRI(4)

      CALL ENTERS('TRI_NE_NP',*9999)

C... nodes already created now. KSB June 2002
      DO i=1,N_VERTICES
        np=NPI(I) !nodes already created, just changed co-ordinates
        XP(1,1,1,np)=VERTICES_XYZ(3*i-2)
        XP(1,1,2,np)=VERTICES_XYZ(3*i-1)
        XP(1,1,3,np)=VERTICES_XYZ(3*i)
      ENDDO !i
C... create 1D elements for delaunay triangulation
      noelem=NEELEM(0,nr)
      ne=NET(0)
      DO ntri=1,N_TRIANGLES
        DO I=1,3
          TRI(I)=TRIANGLES(3*(ntri-1)+I)
C... below puts in correct node # using NPI(I)=np
          TRI(I)=NPI(TRI(I)) !correct np #
        ENDDO
        TRI(4)=TRI(1)
        DO I=1,3
          ne=ne+1
          noelem=noelem+1
          CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',
     '      ERROR,*9999)
          CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
          NEELEM(noelem,nr)=ne
          NPNE(1,nb,ne)=TRI(I)
          NPNE(2,nb,ne)=TRI(I+1)
          CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(1,nb,ne),NPNE(2,nb,ne),
     &      nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
        ENDDO !I
      ENDDO !ntri

      NEELEM(0,nr)=noelem
      NEELEM(0,0)=ne
      NET(nr)=ne
      NET(0)=NET(nr)

      CALL EXITS('TRI_NE_NP')
      RETURN
 9999 CALL ERRORS('TRI_NE_NP',ERROR)
      CALL EXITS('TRI_NE_NP')
      RETURN 1
      END


