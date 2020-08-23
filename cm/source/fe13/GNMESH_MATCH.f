      SUBROUTINE GNMESH_MATCH(NBJ,NEELEM,NENP,NKJ,NKJE,NORD,
     &  NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,OFFSET,SE,XAB,XP,
     &  EFIELD,LADDER_MESH,REVERSE,TERMINAL_MESH,ERROR,*)     
      
C#### Subroutine: GNMESH_MATCH
C###  Description:
C###    GNMESH_MATCH creates a mesh to duplicate an existing 1D 
C###    branching mesh.
C###    This mesh can be  in the same direction as the original (i.e. a 
C###    direct  copy), or in the reverse direction to the original (i.e. 
C###    the creation of a venous tree that is to be connected to an
C###    arterial tree.
C###    The matching mesh can be connected in a 'ladder' structure - so
C###    with elements connecting every level of the structure.
C****   Created by ARC 2011-01-12
   
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NORD(5,NE_R_M),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 OFFSET(3),SE(NSM,NBFM,NEM),XAB(NORM,NEM),
     &   XP(NKM,NVM,NJM,NPM)
      LOGICAL EFIELD,LADDER_MESH,REVERSE,TERMINAL_MESH
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER cap_conns,cap_term,N,nb,ne,ne0,ne1,ne_m,nj,nj1,nj2,nk,nn,
     &  noelem,noelem0,nonode,nonode0,np,np1,np2,np0,np_m,
     &  NP_MAP(NP_R_M),ns
      
      CALL ENTERS('GNMESH_MATCH',*9999)

      CALL ASSERT(NVM.GE.2,'>> Increase NVM to 2 ',ERROR,*9999)

      nb=NBJ(1,NEELEM(1,nr))
      
      np0=NPT(0) !current highest global node #
      ne0=NET(0) !current highest global element #
      noelem0=NEELEM(0,nr) !current # of elements in nr
      nonode0=NPNODE(0,nr) !current # of nodes in nr


      DO nonode=1,NPNODE(0,nr)
        np=np0+nonode
        np_m=NPNODE(nonode,nr)
        NP_MAP(np_m)=np !maps new to old node numbering
        NPNODE(nonode0+nonode,nr)=np
        DO nj=1,NJT
          XP(1,1,nj,np)=XP(1,1,nj,np_m)+OFFSET(nj) !geometry
          XP(2,1,nj,np)=XP(2,1,nj,np_m) !direction
        ENDDO
        NP_INTERFACE(np,0)=1
        NP_INTERFACE(np,1)=nr
        NENP(np,0,nr)=0 !initialise
        DO nj1=1,3
          DO nj2=1,NJ_LOC(nj1,0,nr)
            nj=NJ_LOC(nj1,nj2,nr)
            NKJ(nj,np)=NKJ(nj,np_m)
            NVJP(nj,np)=NVJP(nj,np_m)
          ENDDO !nj2
        ENDDO !nj1
      ENDDO

      DO noelem=1,NEELEM(0,nr)
        ne=ne0+noelem
        ne_m=NEELEM(noelem,nr)
        NEELEM(noelem0+noelem,nr)=ne
        IF(.NOT.REVERSE)THEN
          NPNE(1,nb,ne)=NP_MAP(NPNE(1,nb,ne_m))
          NPNE(2,nb,ne)=NP_MAP(NPNE(2,nb,ne_m))
          NXI(1,0,ne)=NXI(1,0,ne_m)
          NXI(-1,0,ne)=NXI(-1,0,ne_m)
          DO n=1,NXI(1,0,ne)
            NXI(1,n,ne)=NXI(1,n,ne_m)+ne0
          ENDDO
          DO n=1,NXI(-1,0,ne)
            NXI(-1,n,ne)=NXI(-1,n,ne_m)+ne0
          ENDDO
        ELSE
          NPNE(2,nb,ne)=NP_MAP(NPNE(1,nb,ne_m))
          NPNE(1,nb,ne)=NP_MAP(NPNE(2,nb,ne_m))
          NXI(-1,0,ne)=NXI(1,0,ne_m)
          NXI(1,0,ne)=NXI(-1,0,ne_m)
          DO n=1,NXI(1,0,ne)
            NXI(-1,n,ne)=NXI(1,n,ne_m)+ne0
          ENDDO
          DO n=1,NXI(-1,0,ne)
            NXI(1,n,ne)=NXI(-1,n,ne_m)+ne0
          ENDDO
        ENDIF
        NRE(ne)=nr
        DO nj1=1,3 !NJL_GEOM,NJL_FIBR,NJL_FIEL
          DO nj2=1,NJ_LOC(nj1,0,nr)
            nj=NJ_LOC(nj1,nj2,nr)
            NBJ(nj,ne)=NBJ(nj,ne_m)
          ENDDO !nj
        ENDDO !nj1
        DO ns=1,NST(nb)+NAT(nb)
          SE(ns,nb,ne)=SE(ns,nb,ne_m)
        ENDDO !ns
        DO nn=1,NNT(nb)
          DO nj1=1,3
            DO nj2=1,NJ_LOC(nj1,0,nr)
              nj=NJ_LOC(nj1,nj2,nr)
              NVJE(nn,nb,nj,ne)=NVJE(nn,nb,nj,ne_m)
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,ne)=NKJE(nk,nn,nj,ne_m)
              ENDDO !nk
            ENDDO !nj
          ENDDO !nj1
        ENDDO !nn
        NENP(NPNE(1,nb,ne),0,nr)=NENP(NPNE(1,nb,ne),0,nr)+1
        NENP(NPNE(1,nb,ne),NENP(NPNE(1,nb,ne),0,nr),nr)=ne
        NENP(NPNE(2,nb,ne),0,nr)=NENP(NPNE(2,nb,ne),0,nr)+1
        NENP(NPNE(2,nb,ne),NENP(NPNE(2,nb,ne),0,nr),nr)=ne
        NORD(1,ne)=NORD(1,ne_m) !record generation number (same as before)
        NORD(2,ne)=NORD(2,ne_m) ! Horsfield order same as before
        NORD(3,ne)=NORD(3,ne_m) ! Strahler order same as before
        IF(EFIELD) XAB(nej_cap,ne)=XAB(nej_cap,ne_m)+2.d0 ! Identifies as new mesh
      ENDDO
      np0=np !current highest node
      ne1=ne !current highest element
      noelem0=NEELEM(0,nr)+noelem0
      IF(LADDER_MESH)THEN
        cap_conns=0
        cap_term=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NXI(1,0,ne).EQ.1)THEN
            cap_conns=cap_conns+1
            np1=NPNE(2,nb,ne)
            np2=NP_MAP(np1)
            noelem0=noelem0+1
            ne1=ne1+1
            NEELEM(noelem0,nr)=ne1
            NPNE(1,nb,ne1)=np1
            NPNE(2,nb,ne1)=np2
            NENP(np1,0,nr)=NENP(np1,0,nr)+1
            NENP(np1,NENP(np1,0,nr),nr)=ne1
            NENP(np2,0,nr)=NENP(np2,0,nr)+1
            NENP(np2,NENP(np2,0,nr),nr)=ne1
            NXI(-1,1,ne1)=ne
            NXI(1,1,ne1)=ne+ne0
            NORD(1,ne1)=NORD(1,ne) !gen no
            NORD(2,ne1)=NORD(2,ne) ! Horsfield order same as before
            NORD(3,ne1)=NORD(3,ne) ! Strahler order same as before
            IF(EFIELD) XAB(nej_cap,ne1)=XAB(nej_cap,ne)+1.d0 !dentifies as a connection between meshes
          ENDIF
          IF(NXI(1,0,ne).EQ.0)THEN
            cap_conns=cap_conns+1
            cap_term=cap_term+1
            np1=NPNE(2,nb,ne)
            np2=NP_MAP(np1)
            noelem0=noelem0+1
            ne1=ne1+1
            NEELEM(noelem0,nr)=ne1
            NPNE(1,nb,ne1)=np1
            NPNE(2,nb,ne1)=np2
            NENP(np1,0,nr)=NENP(np1,0,nr)+1
            NENP(np1,NENP(np1,0,nr),nr)=ne1
            NENP(np2,0,nr)=NENP(np2,0,nr)+1
            NENP(np2,NENP(np2,0,nr),nr)=ne1
            NXI(-1,1,ne1)=ne
            NXI(1,1,ne1)=ne+ne0
            NORD(1,ne1)=NORD(1,ne) !gen no
            NORD(2,ne1)=NORD(2,ne) ! Horsfield order same as before
            NORD(3,ne1)=NORD(3,ne) ! Strahler order same as before
            IF(EFIELD) XAB(nej_cap,ne1)=XAB(nej_cap,ne)+1.d0 !FIdentifies as a connection between meshes
          ENDIF
          NRE(ne1)=nr
          DO nj1=1,3 !NJL_GEOM,NJL_FIBR,NJL_FIEL
            DO nj2=1,NJ_LOC(nj1,0,nr)
              nj=NJ_LOC(nj1,nj2,nr)
              NBJ(nj,ne1)=nb
            ENDDO !nj
          ENDDO !nj1
          DO ns=1,NST(nb)+NAT(nb)
            SE(ns,nb,ne1)=1.d0
          ENDDO !ns
          DO nn=1,NNT(nb)
            DO nj1=1,3
              DO nj2=1,NJ_LOC(nj1,0,nr)
                nj=NJ_LOC(nj1,nj2,nr)
                NVJE(nn,nb,nj,ne1)=1
                DO nk=1,NKT(nn,nb)
                  NKJE(nk,nn,nj,ne1)=1
                ENDDO !nk
              ENDDO !nj
            ENDDO !nj1
          ENDDO !nn
        ENDDO
        WRITE(OP_STRING,'('' Number of connections '',i6)')cap_conns
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(i6,'' are at terminals '')')cap_term
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSEIF(TERMINAL_MESH)THEN
        cap_term=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NXI(1,0,ne).EQ.0)THEN
            cap_term=cap_term+1
            np1=NPNE(2,nb,ne)
            np2=NP_MAP(np1)
            noelem0=noelem0+1
            ne1=ne1+1
            NEELEM(noelem0,nr)=ne1
            NPNE(1,nb,ne1)=np1
            NPNE(2,nb,ne1)=np2
            NENP(np1,0,nr)=NENP(np1,0,nr)+1
            NENP(np1,NENP(np1,0,nr),nr)=ne1
            NENP(np2,0,nr)=NENP(np2,0,nr)+1
            NENP(np2,NENP(np2,0,nr),nr)=ne1
            NXI(-1,1,ne1)=ne
            NXI(1,1,ne1)=ne+ne0
            NORD(1,ne1)=NORD(1,ne) !gen no
            NORD(2,ne1)=NORD(2,ne) ! Horsfield order same as before
            NORD(3,ne1)=NORD(3,ne) ! Strahler order same as before
            IF(EFIELD) XAB(nej_cap,ne1)=XAB(nej_cap,ne)+1.d0 !FIdentifies as a connection between meshes
          ENDIF
          NRE(ne1)=nr
          DO nj1=1,3 !NJL_GEOM,NJL_FIBR,NJL_FIEL
            DO nj2=1,NJ_LOC(nj1,0,nr)
              nj=NJ_LOC(nj1,nj2,nr)
              NBJ(nj,ne1)=nb
            ENDDO !nj
          ENDDO !nj1
          DO ns=1,NST(nb)+NAT(nb)
            SE(ns,nb,ne1)=1.d0
          ENDDO !ns
          DO nn=1,NNT(nb)
            DO nj1=1,3
              DO nj2=1,NJ_LOC(nj1,0,nr)
                nj=NJ_LOC(nj1,nj2,nr)
                NVJE(nn,nb,nj,ne1)=1
                DO nk=1,NKT(nn,nb)
                  NKJE(nk,nn,nj,ne1)=1
                ENDDO !nk
              ENDDO !nj
            ENDDO !nj1
          ENDDO !nn
        ENDDO
        WRITE(OP_STRING,'('' Number of connections '',i6)')cap_term
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      
      NPT(nr)=np0
      NET(nr)=ne1
      NPT(0)=NPT(nr)
      NET(0)=NET(nr)
      NEELEM(0,nr)=noelem0
      NPNODE(0,nr)=NPNODE(0,nr)+nonode0
 
      CALL EXITS('GNMESH_MATCH')
      RETURN
 9999 CALL ERRORS('GNMESH_MATCH',ERROR)
      CALL EXITS('GNMESH_MATCH')
      RETURN 1
      END



