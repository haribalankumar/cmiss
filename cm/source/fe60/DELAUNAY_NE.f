      SUBROUTINE DELAUNAY_NE(nb,NBJ,NEELEM,NENP,NKJE,NPNE,NPNODE,nr,NRE,
     '  NVJE,ELEM,TETRA,SE,ERROR,*)

C#### Subroutine: DELAUNAY_NE
C###  Description:
C###    Writes tetrahedral element information to CMISS ipelem format.
CC JMB 16-NOV-2001

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  ELEM(LDELEM,4),TETRA(0:LDTETRA,4)
      REAL*8 SE(NSM,NBFM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,ne,NELEM,nj,nk,nn,noelem,nonode,np,ns,REF

      CALL ENTERS('DELAUNAY_NE',*9999)

      DO nonode=1,NPNODE(0,nr) !initialise
        np=NPNODE(nonode,nr)
        NENP(np,0,nr)=0
      ENDDO
      ! Write out only valid elements (tetrahedral)
      NELEM=0
      REF=TETRA(0,HEAD)
      DO WHILE(REF.NE.0)
        NELEM=NELEM+1
        REF=TETRA(REF,NEXT)
      ENDDO !while
      I=1
      noelem=0 !current region starts from 1
      ne=NET(0) !highest element number
      REF=TETRA(0,HEAD)
      DO WHILE(REF.NE.0)
        noelem=noelem+1
        ne=ne+1
        NEELEM(noelem,nr)=ne
        DO J=1,NPTS
c          WORK(J)=ELEM(REF,J)
          np=NPNODE(ELEM(REF,J),nr)
          CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
          CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
          NPNE(J,nb,ne)=np
          NENP(np,0,nr)=NENP(np,0,nr)+1
          NENP(np,NENP(np,0,nr),nr)=ne
        ENDDO !j
        NRE(ne)=nr
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          NBJ(nj,ne)=nb
        ENDDO !nj
        DO ns=1,NST(nb)+NAT(nb)
          SE(ns,nb,ne)=1.0d0
        ENDDO !ns
        DO nn=1,NNT(nb)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            NVJE(nn,nb,nj,ne)=1
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=1
            ENDDO !nk
          ENDDO !nj
        ENDDO !nn
        ! Increment counter
        REF=TETRA(REF,NEXT)
        I=I+1
      ENDDO !while(ref.ne.0)
      NEELEM(0,nr)=noelem
      NEELEM(0,0)=NEELEM(0,0)+noelem
      NET(nr)=ne
      NET(0)=ne

      CALL EXITS('DELAUNAY_NE')
      RETURN
 9999 CALL ERRORS('DELAUNAY_NE',ERROR)
      CALL EXITS('DELAUNAY_NE')
      RETURN 1
      END


