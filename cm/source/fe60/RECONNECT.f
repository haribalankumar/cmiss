      SUBROUTINE RECONNECT(NBJ,NEELEM,NELIST,NENP,NKJE,NPLIST,NPNE,
     '  NPNODE,nr,NRE,NVJE,NXI,SE,XP,ZA,ERROR,*)

C#### Subroutine: RECONNECT
C###  Description:
C###    Recalculates the delaunay mesh by executing a "flipping"
C###    procedure throughout the mesh.

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERROR_FLAG,ERROR_STRINGC(500),LENGTH,
     '  noelem,ne,nb,nn,nj,nk,ns,njj,np,nonode,
     '  OLDNUMELEMS,NEWNUMELEMS,ELEMSCREATED
      CHARACTER ERROR_STRINGF*100

      CALL ENTERS('RECONNECT',*9999)
C     CALL_RECONNECT=.TRUE.

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        CALL CIRCUM(NBJ,ne,NPNE,nr,NVJE,XP,ZA,ERROR,*9999)
      ENDDO
      
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
        CALL FLIP2D(NBJ,NEELEM,NPNE,nr,NVJE,NXI,XP,ZA,ERROR,*9999)
      ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
        ERROR_FLAG=0
        OLDNUMELEMS=NEELEM(0,nr)
C        CALL FLIP3D(NBJ,NEELEM,NPNE,nr,NVJE,NXI,XP,ZA,ERROR,*9999)
C         Obtain the nonode numbers foreach node number np
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            NPLIST(np)=nonode
          ENDDO
C         Obtain the noelem numbers foreach  number ne
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NELIST(ne)=noelem
          ENDDO
          CALL ALE_flip(NAM,NBJ(1,NEELEM(1,nr)),NBFM,NCM,NELIST,NEM,
     '      NEELEM(0,nr),%VAL(NEIM),NHM,NIM,NJM,NKM,NNM,NPLIST,NPM,NPNE,
     '      NPNODE(0,nr),NVJE,NVM,NXI,XP,ZA,ERROR_FLAG,
     '      ERROR_STRINGC)
        IF(ERROR_FLAG.NE.0) THEN
          CALL CSTRINGLEN(LENGTH,ERROR_STRINGC)
          CALL C2FSTRING(ERROR_STRINGC,LENGTH,ERROR_STRINGF)
          CALL ASSERT(.FALSE.,ERROR_STRINGF,ERROR,*9999)
        ENDIF
        NEWNUMELEMS=NEELEM(0,nr)
        ELEMSCREATED=NEWNUMELEMS-OLDNUMELEMS
        IF(ELEMSCREATED.NE.0) THEN
          WRITE(OP_STRING,'('' >>Warning: New range of elements'//
     '      ' is:'',I6,''   ->'',I6,I6)') NEELEM(1,nr),
     '        NEELEM(NEELEM(0,nr),nr),nr
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ASSERT(NEWNUMELEMS.LE.NEM,'>>Increase NEM',ERROR,*9999)
        CALL ASSERT(NEWNUMELEMS.LE.NE_R_M,'>>Increase NE_R_M',
     '    ERROR,*9999)
        NET(nr)=NET(nr)+ELEMSCREATED
        IF(ELEMSCREATED.GT.0) THEN
          DO noelem=OLDNUMELEMS+1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              nb=NBJ(nj,NEELEM(1,nr))
              NBJ(nj,ne)=nb
              DO nn=1,NNT(NBJ(1,ne))
                DO nk=1,NKT(nn,nb)
                  NKJE(nk,nn,nj,ne)=NKJE(nk,nn,nj,NEELEM(1,nr))
                ENDDO
              ENDDO
            ENDDO
            DO nb=1,NBT
              DO nn=1,NNT(nb)
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  NVJE(nn,nb,nj,ne)=1
                ENDDO
              ENDDO
            ENDDO
            NRE(ne)=nr
            DO ns=1,NST(NBJ(1,ne))+NAT(NBJ(1,ne))
              DO nb=1,NBT
                SE(ns,nb,ne)=SE(ns,nb,NEELEM(1,nr))
              ENDDO
            ENDDO
          ENDDO
        ELSEIF(ELEMSCREATED.LT.0) THEN
          DO noelem=NEELEM(0,nr)+1,OLDNUMELEMS
            ne=NEELEM(noelem,nr)
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              NBJ(nj,ne)=0
            ENDDO
            DO nn=1,NNT(NBJ(1,NEELEM(1,nr)))
              DO nb=1,NBT
                NPNE(nn,nb,ne)=0
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  NVJE(nn,nb,nj,ne)=0
                ENDDO
              ENDDO
            ENDDO
C KAT 2001-05-25: Not resetting as I don't expect this to be used.
C            DO nn=1,NNT(NBJ(1,NEELEM(1,nr)))
C              DO nb=1,NBT
C                DO nk=1,NKT(nn,nb)
C                  NKE(nk,nn,nb,ne)=0
C                ENDDO
C              ENDDO
C            ENDDO
            DO ns=1,NST(NBJ(1,NEELEM(1,nr)))+NAT(NBJ(1,NEELEM(1,nr)))
              DO nb=1,NBT
                SE(ns,nb,ne)=0.0D0
              ENDDO
            ENDDO
            DO nj=0,NEIM
              NXI(0,nj,ne)=0
            ENDDO
            DO nj=1,NJM
              ZA(1,nj,1,ne)=0.d0
            ENDDO
            NEELEM(noelem,nr)=0
            NRE(ne)=0
          ENDDO
        ENDIF
      ENDIF

      CALL CALC_NENP_VORO(NBJ,NEELEM,NENP,NPNE,NPNODE,nr,ERROR,*9999)

      CALL EXITS('RECONNECT')
      RETURN
 9999 CALL ERRORS('RECONNECT',ERROR)
      CALL EXITS('RECONNECT')
      RETURN 1
      END


