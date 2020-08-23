      SUBROUTINE GNMESH_MATCH_TRUNK(DISTAL,NBJ,NEELEM,NELIST2,NENP,NKJ,
     &         NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,NPLIST,NPLIST2,nr,NRE,
     &         NVJE,NVJP,NXI,OFFSET,SE,XAB,XP,EFIELD,LADDER_MESH,
     &         REVERSE,TERMINAL_MESH,ERROR,*)
C#### Subroutine: GNMESH_MATCH_TRUNK
C###  Description:
C###    GNMESH_MATCH_TRUNK creates a mesh to duplicate an existing 1D 
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
      INTEGER DISTAL,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST2(0:NEM),
     &  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NORD(5,NE_R_M),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),NPLIST(0:NPM),NPLIST2(0:NPM),nr,NRE(NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 OFFSET(3),SE(NSM,NBFM,NEM),XAB(NORM,NEM),
     &   XP(NKM,NVM,NJM,NPM)
      LOGICAL EFIELD,LADDER_MESH,REVERSE,TERMINAL_MESH
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER cap_conns,cap_term,ERR,N,nb,ne,ne0,ne1,ne_m,nj,nj1,nj2,nk,
     &  nn,
     &  noelem,noelem0,nonode,nonode0,not_copied,np,np1,np2,np0,np_m,
     &  NP_MAP(NP_R_M),ns
      
      CALL ENTERS('GNMESH_MATCH_TRUNK',*9999)

      nb=NBJ(1,NEELEM(1,nr))      
      np0=NPT(0) !current highest global node #
      ne0=NET(0) !current highest global element #
      noelem0=NEELEM(0,nr) !current # of elements in nr
      nonode0=NPNODE(0,nr) !current # of nodes in nr
C...  Initial set up - check enough nodes and elements and which ones to copy
      CALL ASSERT(NVM.GE.2,'>> Increase NVM to 2 ',ERROR,*9999)
      CALL ASSERT(NEM.GE.(NELIST2(0)+NET(0)),'>> Increase NEM',
     &  ERROR,*9999)

      not_copied=0
      !group all nodes in the element group to be copied
      NPLIST2(0)=0
      DO noelem=1,NELIST2(0) !loop through the element list
         ne=NELIST2(noelem)
         nb=NBJ(1,ne)
         DO nn=1,NNT(nb)
            NPLIST2(0)=NPLIST2(0)+1
           NPLIST2(NPLIST2(0))=NPNE(nn,nb,ne)
         ENDDO
      ENDDO
      CALL ILISTRMDUP(NPLIST2(0),NPLIST2(1),ERROR,*9999)
      IF((NPLIST2(0)+np0).GT.NPM) THEN
        CALL WRITE_CHAR(IOER,'Increase NPM to at least ',ERR)
        CALL WRITE_INT(IOER,NPLIST2(0),ERR)
        CALL WRITE_CHAR(IOER,NEWLINE,ERR)
      ENDIF

      NPLIST(0)=0 
      not_copied=0 
        DO nonode=1,NPLIST2(0)
          np_m=NPLIST2(nonode)
          np=NPLIST2(NPLIST2(0))+np_m
          NPLIST(np_m)=np !maps new to old node numbering for non-copied nodes that < than distal point
          NPLIST(np_m)=np !maps new to old node numbering for non-copied nodes that < than distal point
          NPLIST(0)=NPLIST(0)+1
          IF(np_m.LE.DISTAL)THEN      
            IF(np_m.eq.DISTAL)write(*,*) 'Copying ',NPLIST2(0)-NPLIST(0)
     &       ,' nodes and ',NELIST2(0),' elems'
            not_copied=not_copied+1
            IF(NXI(1,0,NENP(np,1,nr)).EQ.0)THEN !identifying terminals and 'attaching' them to new elts
!              NENP(np,0,nr)=NENP(np,0,nr)+1
!              DO np1=1,NENP(np,0,nr)
!                NENP(np,1,nr)=NENP(np_m,1,nr)+NELIST2(NELIST2(0))
!              ENDDO
              IF(NXI(1,0,NENP(np_m,1,nr)).NE.0)THEN
              NXI(1,0,NENP(np,1,nr))=NXI(1,0,NENP(np_m,1,nr)) 
              DO n=1,NXI(1,0,NENP(np,1,nr))
               NXI(1,n,NENP(np,1,nr))=NXI(1,n,NENP(np_m,1,nr))
              ENDDO
              ENDIF
           DO nj1=1,3
            DO nj2=1,NJ_LOC(nj1,0,nr)
              nj=NJ_LOC(nj1,nj2,nr)
              NKJ(nj,np)=NKJ(nj,np_m)
              NVJP(nj,np)=NVJP(nj,np_m)
            ENDDO !nj2
           ENDDO !nj1
            ENDIF
          ELSE !just do a standard copy
            NPNODE(np,nr)=np ! NEW NODE
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
          ENDIF       
        ENDDO

      DO noelem=1,NELIST2(0)
        ne=ne0+noelem
        ne_m=NELIST2(noelem)
        NEELEM(noelem0+noelem,nr)=ne
        IF(.NOT.REVERSE)THEN
          NPNE(1,nb,ne)=NPLIST(NPNE(1,nb,ne_m))
          NPNE(2,nb,ne)=NPLIST(NPNE(2,nb,ne_m))
          NXI(1,0,ne)=NXI(1,0,ne_m)
          NXI(-1,0,ne)=NXI(-1,0,ne_m)
          DO n=1,NXI(1,0,ne)
            NXI(1,n,ne)=NXI(1,n,ne_m)+NELIST2(NELIST2(0))
          ENDDO
          DO n=1,NXI(-1,0,ne)
            NXI(-1,n,ne)=NXI(-1,n,ne_m)+NELIST2(NELIST2(0))
          ENDDO
        ELSE
          NPNE(2,nb,ne)=NPLIST(NPNE(1,nb,ne_m))
          NPNE(1,nb,ne)=NPLIST(NPNE(2,nb,ne_m))
          NXI(-1,0,ne)=NXI(1,0,ne_m)
          NXI(1,0,ne)=NXI(-1,0,ne_m)
          DO n=1,NXI(1,0,ne)
            NXI(-1,n,ne)=NXI(1,n,ne_m)+NELIST2(NELIST2(0))
          ENDDO
          DO n=1,NXI(-1,0,ne)
            NXI(1,n,ne)=NXI(-1,n,ne_m)+NELIST2(NELIST2(0))
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
        IF(EFIELD) XAB(nej_cap,ne)=XAB(nej_cap,ne_m)+2.d0 ! Identifies as new mesh
      ENDDO

      np0=np !current highest node
      ne1=ne !current highest element
      noelem0=NEELEM(0,nr)+NELIST2(0)
      IF(TERMINAL_MESH)THEN
        cap_term=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(ne.LE.NELIST2(NELIST2(0)))THEN
            IF(NXI(1,0,ne).EQ.0)THEN
              cap_term=cap_term+1
              np1=NPNE(2,nb,ne)
              np2=np1+NPLIST2(NPLIST2(0))
              noelem0=noelem0+1
              ne1=ne1+1
              NEELEM(noelem0,nr)=ne1
              NPNE(1,nb,ne1)=np1
              NPNE(2,nb,ne1)=np2
              NENP(np1,0,nr)=NENP(np1,0,nr)+1
              NENP(np1,NENP(np1,0,nr),nr)=ne1
              NENP(np2,0,nr)=NENP(np2,0,nr)+1
              NENP(np2,NENP(np2,0,nr),nr)=ne1
              IF(NENP(np1,0,nr).GT.2)write(*,*)np1,NENP(np1,0,nr),ne,ne1
              IF(NENP(np2,0,nr).GT.2)write(*,*)np2,NENP(np2,0,nr),ne,ne1
 
              NXI(-1,1,ne1)=ne
              NXI(1,1,ne1)=ne+NELIST2(NELIST2(0))
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
          ENDIF
        ENDDO
        WRITE(OP_STRING,'('' Number of connections '',i6)')cap_term
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      
      NPT(nr)=np0
      NET(nr)=ne1
      NPT(0)=NPT(nr)
      NET(0)=NET(nr)
      NEELEM(0,nr)=noelem0
      NPNODE(0,nr)=NPNODE(0,nr)+NPLIST2(0)-not_copied
      write(*,*)'There are now',NEELEM(0,nr),' Elts and',NPNODE(0,nr)
     &  ,' Nodes'

 
      CALL EXITS('GNMESH_MATCH_TRUNK')
      RETURN
 9999 CALL ERRORS('GNMESH_MATCH_TRUNK',ERROR)
      CALL EXITS('GNMESH_MATCH_TRUNK')
      RETURN 1
      END



