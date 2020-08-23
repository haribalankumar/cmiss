      SUBROUTINE UPORDR_DYNAM(nb,NBJ,NEELEM,NELIST,NENP,NE_START,
     &  NE_TEMP,NE_TEMP_OLD,NKJ,NKJE,NORD,NORD_TEMP,NPLIST1,NPLIST2,
     &  NPNE,NPNE_TEMP,NPNODE,nr,NRE,NVJE,NVJP,NXI,OFFSET,SE,
     &  XP,XP_TEMP,ERROR,*)

C#### Subroutine: UPORDR_DYNAM
C###  Description:
C###    UPORDR_DYNAM reorders a 1D tree network so that elements and
C###    nodes increase with order.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'

!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NE_START,NE_TEMP(NEM),NE_TEMP_OLD(NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),
     &  NORD_TEMP(NE_R_M),NPLIST1(0:NPM),NPLIST2(0:NPM),NPNE(NNM,NBFM,
     &  NEM),NPNE_TEMP(2,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  OFFSET
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),XP_TEMP(NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ne,ne0,ne2,nj,njj,nk,nn,noelem,nonode,np,ns,N_TB,
     &  N_TB_OLD,nv,nz1,nz2,nzz
      LOGICAL SINGLE

      CALL ENTERS('UPORDR_DYNAM',*9999)

      NELIST(1)=NE_START !store start element
      NPLIST1(1)=NPNE(1,nb,NE_START) !store 1st node of start element
      NPLIST2(NPNE(1,nb,NE_START))=1 !store 1st node of start element
      nzz=0
C  Set up the list of correctly ordered elements and nodes.
      N_TB=1
      NE_TEMP_OLD(1)=NE_START
      DO WHILE(N_TB.NE.0) !while still some to reorder
        N_TB_OLD=N_TB
        N_TB=0
        DO nz1=1,N_TB_OLD !for parents in previous gen
          ne0=NE_TEMP_OLD(nz1)
          IF(DOP)THEN
            WRITE(OP_STRING,'('' Parent '',I6)') ne0
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          nzz=nzz+1 !increment counter
          CALL ASSERT(nzz.LE.NEM,'>>Increase NEM',ERROR,*9999)
          NELIST(nzz)=ne0 !store element #
          CALL ASSERT(nzz+1.LE.NPM,'>>Increase NPM',ERROR,*9999)
          NPLIST1(nzz+1)=NPNE(2,nb,ne0) !end node #
          NPLIST2(NPNE(2,nb,ne0))=nzz+1
          IF(DOP)THEN
            WRITE(OP_STRING,'('' nzz '',I4,''  NELIST(nzz) '',I4,'
     '      //'''  NPLIST1(nzz+1) '',I4,''  NPLIST2(NPNE) '',I4)') ne0,
     '        NELIST(nzz),NPLIST1(nzz+1),NPLIST2(NPNE(2,nb,ne0))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          SINGLE=.TRUE.
          DO WHILE(SINGLE)
            IF(NXI(1,0,ne0).EQ.1)THEN
              ne=NXI(1,1,ne0) !daughter element #
              nzz=nzz+1 !increment counter
              NELIST(nzz)=ne !store element #
              NPLIST1(nzz+1)=NPNE(2,nb,ne) !end node #
              NPLIST2(NPNE(2,nb,ne))=nzz+1
              ne0=ne
             IF(DOP)THEN
              WRITE(OP_STRING,
     '          '('' NXI=1, ne0_init '',I4,''  daughter ne '',I4,'
     '          //'''  ne0_new '',I4)')
     '          NELIST(nzz-1),ne,ne0
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             ENDIF
            ELSE IF(NXI(1,0,ne0).GE.2)THEN
              DO nz2=1,NXI(1,0,ne0) !for each daughter
                ne=NXI(1,nz2,ne0) !daughter element #
                N_TB=N_TB+1
                NE_TEMP(N_TB)=ne
              ENDDO !nz2
              SINGLE=.FALSE.
             IF(DOP)THEN
              WRITE(OP_STRING,'('' NXI=2'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             ENDIF
            ELSE IF(NXI(1,0,ne0).EQ.0)THEN
              SINGLE=.FALSE.
             IF(DOP)THEN
              WRITE(OP_STRING,'('' NXI=0'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             ENDIF
            ENDIF
          ENDDO !WHILE
        ENDDO !nz1
        DO nz1=1,N_TB
          NE_TEMP_OLD(nz1)=NE_TEMP(nz1)
        ENDDO !nz1
      ENDDO !WHILE
C  Put values into temporary arrays; move to correct array positions.
C     Assumes that new element and node numbering will start from 1.
      NEELEM(0,nr)=nzz
      DO noelem=1,NEELEM(0,nr)
        ne=NELIST(noelem) !for each reordered element
        IF(ne.NE.0)THEN
          NPNE_TEMP(1,noelem)=NPLIST2(NPNE(1,nb,ne))
          NPNE_TEMP(2,noelem)=NPLIST2(NPNE(2,nb,ne))
          NORD_TEMP(noelem)=NORD(5,ne)
        ENDIF
      ENDDO
      DO noelem=1,NEELEM(0,nr)
        ne=noelem+offset
        CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
        IF(ne.NE.0)THEN
          NPNE(1,nb,ne)=NPNE_TEMP(1,noelem)+offset
          NPNE(2,nb,ne)=NPNE_TEMP(2,noelem)+offset
c          NE_BBM(NEELEM(noelem,nr))=ne !map old to new
          NEELEM(noelem,nr)=ne
          NRE(ne)=nr
          NORD(5,ne)=NORD_TEMP(noelem)
          
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            NBJ(nj,ne)=nb
          ENDDO !nj
          DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            NBJ(nj,ne)=nb
          ENDDO !nj

          DO ns=1,NST(nb)+NAT(nb)
            SE(ns,nb,ne)=1.0d0
          ENDDO !ns

C         Set the number of versions at np in ne for both geometry and
C         field to be 1 (in NVJE); set a single value (no derivatives)
C         for geometry and field values (in NKJE). Note that after
C         CALC_NENP_1D is called NVJE is updated for the field versions.
          DO nn=1,NNT(nb)
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              NVJE(nn,nb,nj,ne)=1
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,ne)=1
              ENDDO !nk
            ENDDO !nj
            DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
              nj=NJ_LOC(NJL_FIEL,njj,nr)
              NVJE(nn,nb,nj,ne)=1
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,ne)=1
              ENDDO !nk
            ENDDO !nj
          ENDDO !nn
        ENDIF
      ENDDO !noelem

      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999) ! -> NENP
      CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*9999) ! -> NXI

C     New node numbering 
      DO nonode=1,NPNODE(0,nr)
        NPNODE(nonode,nr)=nonode+offset
      ENDDO !nonode
      
C     Update the versions of nodes for the fields.
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        IF(np.NE.0)THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            NKJ(nj,np+offset)=1
            NVJP(nj,np+offset)=1
          ENDDO !njj
          DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            NKJ(nj,np+offset)=1
c            NVJP(nj,np+offset)=NENP(np,0,nr)
            NVJP(nj,np+offset)=NENP(nonode,0,nr)
          ENDDO !njj
        ENDIF !np
      ENDDO !nonode
C  Sets up field arrays
      IF(nj_radius.NE.0.OR.nj_alveoli.NE.0) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            DO i=1,NENP(np,0,nr)
              ne2=NENP(np,i,nr)
              IF(ne2.EQ.ne)THEN
                IF(nj_radius.NE.0) NVJE(nn,nb,nj_radius,ne)=i
                !the version of node at nn'th position
                IF(nj_alveoli.NE.0) NVJE(nn,nb,nj_alveoli,ne)=i
                !the version of node at nn'th position
              ENDIF
            ENDDO !i
            IF(nj_radius.NE.0) NKJE(1,nn,nj_radius,ne)=1
            IF(nj_alveoli.NE.0) NKJE(1,nn,nj_alveoli,ne)=1
            IF(nj_alveoli.NE.0) nv=NVJE(nn,nb,nj_alveoli,ne)
            IF(nj_alveoli.NE.0) XP(1,nv,nj_alveoli,np)=1.d0 !default value
            SE(nn,nb,ne)=1.d0
          ENDDO !nn
        ENDDO !noelem
      ELSE
        WRITE(OP_STRING,
     &    '('' Fields not set up: define mesh field first '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF        
      DO nonode=1,NPNODE(0,nr)
        np=NPLIST1(nonode)
        IF(np.NE.0)THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DO nv=1,NVJP(nj,nonode)
              XP_TEMP(nv,nj,nonode)=XP(1,nv,nj,np)
            ENDDO !nv
          ENDDO !njj
          DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            DO nv=1,NVJP(nj,nonode)
              XP_TEMP(nv,nj,nonode)=XP(1,nv,nj,np)
            ENDDO !nv
          ENDDO !njj
        ENDIF
      ENDDO !nonode
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        IF(np.NE.0)THEN
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DO nv=1,NVJP(nj,np)
              XP(1,nv,nj,np+offset)=XP_TEMP(nv,nj,np)
            ENDDO !nv
          ENDDO !njj
          DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            DO nv=1,NVJP(nj,np)
              XP(1,nv,nj,np+offset)=XP_TEMP(nv,nj,np)
            ENDDO !nv
          ENDDO !njj
c          NPNODE(nonode,nr)=np+offset
        ENDIF
      ENDDO !nonode

      WRITE(OP_STRING,'('' Write out nodes, elements, and fields'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('UPORDR_DYNAM')
      RETURN
 9999 CALL ERRORS('UPORDR_DYNAM',ERROR)
      CALL EXITS('UPORDR_DYNAM')
      RETURN 1
      END


