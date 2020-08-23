      SUBROUTINE IPMESH2(MIN_ORDER,NBJ,NEELEM,NELIST,NELIST2,
     &  NENP,NKJE,NKJ,NLL,NNB,NORD,NP_ATTACH,NP_INTERFACE,NPL,
     &  NPLIST,NPNE,NPNODE,nr,NRE,Nrefine,nr_host,NSTORE,NVJE,NVJP,NXI,
     &  TERMINAL_ORDER,BBM,CE,RSTORE,SE,Spread,XAB,XP,XSTORE,ADD_SUPER,
     &  ARTERIES,LUNG_TOTAL,ORDER,REDUCE,VEINS,ERROR,*)
      
C#### Subroutine: IPMESH2
C###  Description:
C###    IPMESH2 defines mesh parameters for fractal trees.

C****   NORD(1,ne) Generation
C****   NORD(2,ne) Order (Horsfield)
C****   NORD(3,ne) Order (Strahler)
C****   NORD(4,ne) Order (Diameter-defined Strahler)
C****   NORD(5,ne) Indicator of branch ending for symmetric tree

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER MIN_ORDER,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NELIST2(0:NEM),
     &  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     &  NLL(12,NEM),NNB(4,4,4,NBFM),NORD(5,NE_R_M),
     &  NP_ATTACH,NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     &  NRE(NEM),Nrefine,nr_host,NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),TERMINAL_ORDER
      REAL*8 BBM(2,NEM),CE(NMM,NEM),SE(NSM,NBFM,NEM),Spread,
     &  XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),XSTORE(NE_R_M,3)
      LOGICAL ADD_SUPER,ARTERIES,LUNG_TOTAL,ORDER,REDUCE,VEINS
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER gen,i,IBEG,ICHAR,IEND,INFO,INLET,j,MAX_STRAHLER,nb,
     &  nb_airway,nb_cap,nb_circ,nb_host,N_BOUNDARY_EST,ne,n_generation,
     &  nhx,ni,nj,nn,noelem,nogrel,NOQUES,np,np1,np2,nr0,nr_feed,
     &  NSTORE(NE_R_M,3),nv1,nv2,OUTLET,NORDER
      INTEGER*4 NBJ_TEMP_PTR,NEELEM_TEMP_PTR,NKJ_TEMP_PTR,NKJE_TEMP_PTR,
     &  NORD_TEMP_PTR,NPNE_TEMP_PTR,NPNODE_TEMP_PTR,NRE_TEMP_PTR,
     &  NVJE_TEMP_PTR,NVJP_TEMP_PTR,XP_TEMP_PTR
      REAL*8 a_A,ANG,C(MAX_ALVEOLI,NJT),centre(NJT),centre2(NJT),
     &  depth,depth_sum,DOT_PROD,length,length1,length2,length3,
     &  RSTORE(NE_R_M,3),SA,SA_sum,SA_total,SA_vol,start(3),V1(NJT),
     &  V2(NJT),V3(NJT),V4(NJT),vol,vol_total
      CHARACTER CHAR1*5,CHAR2*2,CHAR3*26,CHAR6*6,CHAR7*6,LABEL*30,
     &  MESH_TYPE*30
      LOGICAL FILEIP,GROUP,REGULAR

      CALL ENTERS('IPMESH2',*9999)

      CIRCULATION=.FALSE.
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      NTT_GEN=0
      IF(ARTERIES.OR.VEINS)THEN
        CIRCULATION=.TRUE.
        AIRWAY=.FALSE.
      ELSE IF(LUNG_TOTAL)THEN
        NBJ_TEMP_PTR=0
        NEELEM_TEMP_PTR=0
        NKJ_TEMP_PTR=0
        NKJE_TEMP_PTR=0
        NORD_TEMP_PTR=0
        NPNE_TEMP_PTR=0
        NPNODE_TEMP_PTR=0
        NRE_TEMP_PTR=0
        NVJE_TEMP_PTR=0
        NVJP_TEMP_PTR=0
        XP_TEMP_PTR=0
        CALL ALLOCATE_MEMORY(NJM*NEM,0,INTTYPE,NBJ_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NE_R_M,0,INTTYPE,NEELEM_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NJM*NPM,0,INTTYPE,NKJ_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NKM*NNM*NJM*NEM,0,INTTYPE,NKJE_TEMP_PTR,
     &    .TRUE.,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NE_R_M,0,INTTYPE,NORD_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NNM*NBFM*NEM,0,INTTYPE,NPNE_TEMP_PTR,
     &    .TRUE.,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NP_R_M,0,INTTYPE,NPNODE_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NEM,0,INTTYPE,NRE_TEMP_PTR,.TRUE.,ERROR,
     &    *9999)
        CALL ALLOCATE_MEMORY(NNM*NBFM*NJM*NEM,0,INTTYPE,NVJE_TEMP_PTR,
     &    .TRUE.,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NJM*NPM,0,INTTYPE,NVJP_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NKM*NVM*NJM*NPM,0,DPTYPE,XP_TEMP_PTR,
     &    .TRUE.,ERROR,*9999)

        CALL IPMESH2_1(NBJ,%VAL(NBJ_TEMP_PTR),NEELEM,
     &    %VAL(NEELEM_TEMP_PTR),NELIST2,NENP,NKJ,%VAL(NKJ_TEMP_PTR),
     &    NKJE,%VAL(NKJE_TEMP_PTR),NORD,%VAL(NORD_TEMP_PTR),NP_ATTACH,
     &    NP_INTERFACE,NPLIST,NPNE,%VAL(NPNE_TEMP_PTR),NPNODE,
     &    %VAL(NPNODE_TEMP_PTR),nr,NRE,%VAL(NRE_TEMP_PTR),NVJE,
     &    %VAL(NVJE_TEMP_PTR),NVJP,%VAL(NVJP_TEMP_PTR),NXI,SE,XP,
     &    %VAL(XP_TEMP_PTR),ERROR,*9999)
          CALL FREE_MEMORY(NBJ_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NEELEM_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NKJ_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NKJE_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NORD_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NPNE_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NPNODE_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NVJE_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(NVJP_TEMP_PTR,ERROR,*9999)
          CALL FREE_MEMORY(XP_TEMP_PTR,ERROR,*9999)
          IF(NPT(nr).GT.NPT(0)) NPT(0)=NPT(nr)
          AIRWAY=.FALSE.
          CIRCULATION=.FALSE.
      ELSE
        FORMAT='('' The mesh is for the [1]:'''//
     '    '/''   (1) Airways'''//
     '    '/''   (2) Circulation'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=1
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3)THEN
          IF(IDATA(1).EQ.1)THEN
            AIRWAY=.TRUE.
          ELSE IF(IDATA(1).EQ.2)THEN
            AIRWAY=.FALSE.
            CIRCULATION=.TRUE.
          ENDIF
        ENDIF
      ENDIF

      IF(AIRWAY)THEN !airways
C***  AIRWAY MODEL
        FORMAT='('' The airway tree is [1]:'''//
     '    '/''   (1) symmetric'''//
     '    '/''   (2) multi-branching'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=LTYP_C(1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) LTYP_C(1)=IDATA(1)
        IF(LTYP_C(1).LE.3)THEN
          FORMAT='($,'' Enter the basis function [1]: '',I3)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) nb_airway=IDATA(1)
        ENDIF
        IF(LTYP_C(1).EQ.1)THEN !symmetric conducting airways
C***    SYMMETRIC CONDUCTING AIRWAY MODEL
          MESH_TYPE='SYMMETRIC'
          FORMAT=
     &      '($,'' Enter the number of airway generations [1]:  '',I3)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=NTT_GEN
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,30,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NTT_GEN=IDATA(1)

          DO n_generation=1,NTT_GEN
            WRITE(CHAR2,'(I2)') n_generation
            FORMAT='($,'' Enter the number of elements in '//
     &        'generation '//CHAR2//' [1]:  '',I6)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=B_BRANCH_GEN(n_generation)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) B_BRANCH_GEN(n_generation)=IDATA(1)
              
            FORMAT='($,'' Enter the length of each element in '//
     &        'generation '//CHAR2//'[1.0]: '',D12.4)'
            RDEFLT(1)=1.d0
            IF(IOTYPE.EQ.3) RDATA(1)=L_BRANCH_GEN(n_generation)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) L_BRANCH_GEN(n_generation)=RDATA(1)
          ENDDO !n_generation

        ELSE IF(LTYP_C(1).EQ.2)THEN !multi-branching
c          NTB=0
C***    Multi-branching respiratory airway model
          FORMAT='('' Enter multi-branching model type [1]:'''//
     '      '/''   (1) Define mean branch dimensions'''//
     '      '/''   (2) Define asymmetry parameter'''//
     '      '/''   (3) Read in from file'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=LTYP_R(2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) LTYP_R(2)=IDATA(1)
          IF(LTYP_R(2).EQ.1)THEN
            MESH_TYPE='MULTI_BRANCH'
            FORMAT='($,'' Enter number of generations [9]: '',I2)'
            IDEFLT(1)=9
            IF(IOTYPE.EQ.3) IDATA(1)=9
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NT_GEN=IDATA(1)

            DO gen=1,NT_GEN
              WRITE(CHAR2,'(I2)') gen
              FORMAT='($,'' Enter branching ratio mean and SD for '//
     '          'generation '//CHAR2//' [2,0]: '',2(D12.4))'
              RDEFLT(1)=2.0d0
              RDEFLT(1)=0.0d0
              IF(IOTYPE.EQ.3)THEN
                RDATA(1)=MNB(gen)
                RDATA(2)=SDNB(gen)
              ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,3.0d0,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3)THEN
                MNB(gen)=RDATA(1)
                SDNB(gen)=RDATA(2)
              ENDIF
              FORMAT='($,'' Enter the length of each element in '//
     &          'generation '//CHAR2//'[1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=L_BRANCH_GEN(gen)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) L_BRANCH_GEN(gen)=RDATA(1)
            ENDDO !gen
          ELSE IF(LTYP_R(2).EQ.3)THEN !read in
            MESH_TYPE='MULTI_READ'
            FORMAT='($,'' Enter the number of respiratory '//
     '        'airways [1]: '',I6)'
            IDEFLT(1)=1
c            N_RESP=IDEFLT(1)
            IF(IOTYPE.EQ.3) IDATA(1)=N_RESP
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NE_R_M,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) N_RESP=IDATA(1)
c            IF(IOTYPE.EQ.3)THEN
c              WRITE(STRING,'(''Enter parent element, generation, and '//
c     &          'symmetry (1/0) for:'')')
c              CALL WRITES(IFILE,STRING,ERROR,*9999)
c            ENDIF

            IF(IOTYPE.EQ.3)THEN
              ne=NEELEM(1,nr) !must be the first element
              nb=NBJ(1,ne)
              np=NPNE(1,nb,ne)
              DO nj=1,NJT
                start(nj)=XP(1,1,nj,np)
              ENDDO !nj
            ENDIF

            DO noelem=1,N_RESP
              WRITE(CHAR6,'(I6)') noelem !+759 !local acinus element #
              WRITE(CHAR7,'(I6)') noelem-1 !+759
              FORMAT='($,'' Element '//CHAR6//' '//
     '          '['//CHAR7//',1,0]: '',3(I6)))'
              IDEFLT(1)=noelem-1
              IDEFLT(2)=1
              IDEFLT(3)=0
              NSTORE(noelem,1)=IDEFLT(1) 
              NSTORE(noelem,2)=IDEFLT(2)
              NSTORE(noelem,3)=IDEFLT(3)
              IF(IOTYPE.EQ.3)THEN
                ne=NEELEM(noelem,nr)
                IDATA(1)=NXI(-1,1,ne)
                IDATA(2)=NORD(1,ne)
                IDATA(3)=NORD(5,ne)
              ENDIF
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,NE_R_M,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3)THEN
                NSTORE(noelem,1)=IDATA(1) !-759 !local parent
                NSTORE(noelem,2)=IDATA(2) !-16 !local generation
                NSTORE(noelem,3)=IDATA(3) !symmetry
              ENDIF
              
            ENDDO !noelem
            DO noelem=1,N_RESP !for each acinus element
              WRITE(CHAR6,'(I6)') noelem !+759 !local acinus element #
              FORMAT='($,'' Element '//CHAR6//' '//
     '          '[1.0,1.0,1.0]: '',3(D12.4))'
              IF(IOTYPE.EQ.3)THEN
                ne=NEELEM(noelem,nr)
                nb=NBJ(1,ne)
                np1=NPNE(1,nb,ne)
                np2=NPNE(2,nb,ne)
                nv1=NVJE(1,nb,1,ne)
                nv2=NVJE(2,nb,1,ne)
                length=0.d0
                DO nj=1,NJT
                  length=length+(XP(1,nv1,nj,np1)-XP(1,nv2,nj,np2))**2
                ENDDO
                nv1=NVJE(1,nb,nj_radius,ne)
                nv2=NVJE(2,nb,nj_radius,ne)
                a_A=0.5d0*(XP(1,nv1,nj_alveoli,np1)+XP(1,nv2,nj_alveoli,
     &            np2))
                RDEFLT(1)=DSQRT(length)
                RDEFLT(2)=0.5d0*(XP(1,nv1,nj_radius,np1)+XP(1,nv2,
     &            nj_radius,np2))
                RDEFLT(3)=a_A

                nb=NBJ(1,ne)
                np=NPNE(2,nb,ne)
              ELSE
                DO I=1,3
                  RDEFLT(I)=1.d0
                ENDDO !I
              ENDIF
              DO I=1,3
                RDATA(I)=RDEFLT(I)
              ENDDO !I

              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,NE_R_M,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3)THEN
                RSTORE(noelem,1)=RDATA(1) !length
                RSTORE(noelem,2)=RDATA(2) !radius
                RSTORE(noelem,3)=RDATA(3) !a/A
              ENDIF
            ENDDO !noelem

            DO noelem=1,N_RESP !for each acinus element
              WRITE(CHAR6,'(I6)') noelem !+759 !local acinus element #
              FORMAT='($,'' Element '//CHAR6//' '//
     '          '[1.0,1.0,1.0]: '',3(F12.4))'
              IF(IOTYPE.EQ.3)THEN
                ne=NEELEM(noelem,nr)
                nb=NBJ(1,ne)
                np=NPNE(2,nb,ne)
                DO nj=1,NJT
                  RDEFLT(nj)=XP(1,1,nj,np)-start(nj)
                ENDDO !nj
              ELSE
                DO I=1,3
                  RDEFLT(I)=0.d0
                ENDDO !I
              ENDIF
              DO I=1,3
                RDATA(I)=RDEFLT(I)
              ENDDO !I

              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,NE_R_M,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3)THEN
                XSTORE(noelem,1)=RDATA(1)
                XSTORE(noelem,2)=RDATA(2)
                XSTORE(noelem,3)=RDATA(3)
              ENDIF
            ENDDO !noelem
            
          ENDIF !LTYP_R(2).EQ.3
          B_ANGLE_XY(2)=90.d0*PI/180.0d0

        ENDIF
      ELSE IF(CIRCULATION)THEN !circulation
        IF(VEINS)THEN
          CALL GNVEIN(NBJ,NEELEM,NENP,NKJ,NKJE,NP_INTERFACE,NPLIST,
     &      NPNE,NPNODE,nr_host,NRE,nr,NVJE,NVJP,NXI,CE,SE,XP,ERROR,
     &      *9999)
          IF(DOP) THEN
            WRITE(OP_STRING,'(''Done GNVEIN'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(ARTERIES)THEN
C         CALL ASSERT()CHECK AIRWAY MESH ALREADY DEFINED
          CALL GNARTRY(MIN_ORDER,NBJ(1,NEELEM(1,nr_host)),NBJ,NEELEM,
     '      NENP,NKJ,NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,1,NRE,2,NVJE,
     '      NVJP,NXI,CE,SE,XP,ADD_SUPER,ERROR,
     '      *9999)
          IF(DOP) THEN
            WRITE(OP_STRING,'(''Done GNARTRY'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE

C KSBs 7-MAY-2001 Enter circulation model
C### Pulmonary Circulation stuff still in experimental stage, this
C### includes GENCIRC & GENPCAP subroutines following.

        FORMAT='($,'' Enter basis function for circulation [3]: '',I1)'
        IDEFLT(1)=3
        IF(IOTYPE.EQ.3) IDATA(1)=3
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) nb_circ=IDATA(1)
        DO nhx=1,2
          WRITE(CHAR3,'(A)') OPT_VESSEL(nhx)
          FORMAT='('' Enter model for '//CHAR3(1:15)//' [1]:'''//
     '      '/''   (0) None'''//
     '      '/''   (1) Strahler ordering model'''//
     '      '/''   (2) Horsfield ordering model'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) CIRC_MODEL(nhx)=IDATA(1)
        ENDDO !nhx
        DO nhx=1,2
          IF(CIRC_MODEL(nhx).NE.0) THEN
            WRITE(OP_STRING,'('' Warning, nr not correct'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            nr=nr+1
            WRITE(OP_STRING,'(A)') OPT_VESSEL(nhx)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !CIRC_MODEL(nhx).NE.0
          IF(CIRC_MODEL(nhx).EQ.1) THEN
            FORMAT='($,'' Enter highest strahler order [11]: '',I2)'
            IDEFLT(1)=11
            IF(IOTYPE.EQ.3) IDATA(1)=11
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        30,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) MAX_STRAHLER=IDATA(1)
            FORMAT='($,'' Solve down to Strahler order # [7]: '',I2)'
            IDEFLT(1)=7
            IF(IOTYPE.EQ.3) IDATA(1)=7
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        11,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NORDER=IDATA(1)
            DO nr0=MAX_STRAHLER,1,-1
              IF(nhx.EQ.1) IDEFLT(1)=N_STRAHLER_ORD(nr0)
              IF(nhx.EQ.2) IDEFLT(1)=N_STRAHLER_ORDVEN(nr0)
              WRITE(CHAR1,'(I5)') IDEFLT(1)
              WRITE(CHAR2,'(I2)') nr0
              FORMAT='($,'' Enter # branches ['//CHAR1(1:5)
     '          //'] in order ['//CHAR2(1:2)//']: '',I5)'
              IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
            ENDDO !nr0
            FORMAT='($,'' Enter angle for viewing [25]: '''//
     '        ',D13.5)'
            RDEFLT(1)=25d0
            IF(IOTYPE.EQ.3) RDATA(1)=25d0
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        50,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANG=RDATA(1)
          ENDIF !CIRC_MODEL(nhx).EQ.1
          IF(CIRC_MODEL(nhx).NE.0) THEN
            CALL GENCIRC(MAX_STRAHLER,nb_circ,NBJ,NEELEM,NPNE,NPNODE,nr,
     '        NRE,NXI,NORDER,ANG,CE,nhx,XP,ERROR,*9999)
          ENDIF !CIRC_MODEL(nhx).NE.0
        ENDDO !nhx=1,2
C KSB(start) Sept 2001, capillary option
        FORMAT='('' Capillary bed inlet is: [1]:'''//
     '    '/''   (1) None'''//
     '    '/''   (2) Arteriole created in gencirc'''//
     '    '/''   (3) Create a new arteriole '''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=1
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,LDATA,
     '    LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) INLET=IDATA(1)
        FORMAT='('' Capillary bed outlet is: [1]:'''//
     '    '/''   (1) None'''//
     '    '/''   (2) Venule created in gencirc'''//
     '    '/''   (3) Create a new venule '''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=1
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,LDATA,
     '    LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) OUTLET=IDATA(1)
        FORMAT='('' Capillary mesh is: [1]:'''//
     '    '/''   (1) None  '''//
     '    '/''   (2) Create test mesh on uniform mesh'''//
     '    '/''   (3) Create voronoi mesh on unit sphere'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=1
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,LDATA,
     '    LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) CAP_MESH=IDATA(1)
        FORMAT='($,'' Enter host basis function # [4]: '',I1)'
        IDEFLT(1)=4
        IF(IOTYPE.EQ.3) IDATA(1)=4
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) nb_host=IDATA(1)
        FORMAT='($,'' Enter host region number [1]: '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=1
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NRM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) nr_host=IDATA(1)
        FORMAT='($,'' Enter feed vessel region number [1]: '',I1)'
        IDEFLT(1)=5
        IF(IOTYPE.EQ.3) IDATA(1)=5
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NRM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) nr_feed=IDATA(1)
        IF(CAP_MESH.NE.1)THEN
          FORMAT='($,'' Enter basis function for capillaries [3]: '''
     '      //', I1)'
          IDEFLT(1)=3
          IF(IOTYPE.EQ.3) IDATA(1)=3
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) nb_cap=IDATA(1)
C put in an assert to check the basis function is 1D
          IF(CAP_MESH.EQ.3)THEN
            FORMAT='('' Are the seed points [1]:'''//
     '        '/''   (1) Randomly generated'''//
     '        '/''   (2) Regular '''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     &        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3)THEN
              IF(IDATA(1).EQ.1)THEN
                REGULAR=.FALSE.
              ELSE
                REGULAR=.TRUE.
              ENDIF
            ENDIF
            FORMAT='($,'' Enter angle limit for inlet '//
     '        'orifice [90 deg]: '',D12.4)'
            RDEFLT(1)=90.d0
            IF(IOTYPE.EQ.3) RDATA(1)=90.d0
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,135.d0,INFO,
     &        ERROR,*9999)
            IF(IOTYPE.NE.3) low_limit=RDATA(1)
            IF(.NOT.REGULAR)THEN
              FORMAT='($,'' Enter # random points [100]: '',I5)'
              IDEFLT(1)=100
              IF(IOTYPE.EQ.3) IDATA(1)=100
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          6,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IOTYPE.NE.3) N_VERT_EST=IDATA(1)
            ELSE !IF REGULAR POINTS
              FORMAT='($,'' Enter # of refinements [3]: '',I5)'
              IDEFLT(1)=3
              IF(IOTYPE.EQ.3) IDATA(1)=3
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,4,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IOTYPE.NE.3) N_REFINE=IDATA(1)
C... N_VERT_EST is N_VERTICES estimate used for array sizing
C... initial N_VERTICES * MAX_TRIANGLES DEFINED FOR ARRAY SIZING
C...  correct values calculated in SPHERE_POINTS
              IF(N_REFINE.EQ.0) N_VERT_EST=12
              IF(N_REFINE.EQ.1) N_VERT_EST=72
              IF(N_REFINE.EQ.2) N_VERT_EST=312
              IF(N_REFINE.EQ.3) N_VERT_EST=1272
            ENDIF !.NOT.REGULAR
C... Setting up recommended default values for N_BOUNDARY
C... fix to include # of vertices if want to use random mesh
            IF(low_limit.LT.40.d0.OR.low_limit.GT.135.d0) THEN
              IF(N_REFINE.EQ.0) N_BOUNDARY_EST=8
              IF(N_REFINE.EQ.1) N_BOUNDARY_EST=12
              IF(N_REFINE.EQ.2) N_BOUNDARY_EST=20
              IF(N_REFINE.EQ.3) N_BOUNDARY_EST=30
            ELSE IF(low_limit.GE.40.d0.AND.low_limit.LE.135.d0) THEN
              IF(N_REFINE.EQ.0) N_BOUNDARY_EST=12
              IF(N_REFINE.EQ.1) N_BOUNDARY_EST=16
              IF(N_REFINE.EQ.2) N_BOUNDARY_EST=30
              IF(N_REFINE.EQ.3) N_BOUNDARY_EST=45
            ELSE
              N_BOUNDARY_EST=50
            ENDIF
            IF(low_limit.EQ.0.d0) N_BOUNDARY_EST=0
            IDEFLT(1)=N_BOUNDARY_EST
            WRITE(CHAR1,'(I5)') IDEFLT(1)
            FORMAT='($,'' Enter # boundary points ['//CHAR1//']: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=N_BOUNDARY_EST
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) N_BOUNDARY=IDATA(1)
C... NB/ insufficient # of boundary pts can cause problems in mesh
            IF(N_BOUNDARY.LT.N_BOUNDARY_EST) THEN
              WRITE(OP_STRING,'('' WARNING:insufficient # of boundary'//
     '          ' points can cause problems in mesh'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            N_VERT_EST=N_VERT_EST+2*N_BOUNDARY
            MAX_TRIANGLES=N_VERT_EST*4 !redefined in sphere_points
            FORMAT='('' The capillary mesh is projected onto [1]:'''//
     '        '/''   (1) Test alveolar mesh'''//
     '        '/''   (2) Voronoi alveolar mesh'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3)THEN
              IF(IDATA(1).EQ.1)THEN
                VORONOI_ALV=.FALSE.
              ELSE IF(IDATA(1).EQ.2)THEN
                VORONOI_ALV=.TRUE.
              ENDIF
            ENDIF
            N_ALVEOLI=0 !initialise
            NETALV=0 !initialise total # of alveoli elements
            vol_total=0.d0 !initialise total alveolar volume
            SA_vol=0.d0
            SA_sum=0.d0
            depth_sum=0.d0
 802        FORMAT='($,'' Enter alveolar element group name '
     '        //' [EXIT]: '',A32)'
            CDEFLT(1)='EXIT'
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(CDATA(1).NE.'EXIT') THEN !not default exit
              N_ALVEOLI=N_ALVEOLI+1
              CALL CUPPER(CDATA(1),CDATA(1))
              CALL STRING_TRIM(CDATA(1),IBEG,IEND)
              ALVEOLI_NAMES(N_ALVEOLI)=CDATA(1)(IBEG:IEND)
              DO i=1,N_ALVEOLI-1
                CALL ASSERT(ALVEOLI_NAMES(i).NE.
     '            ALVEOLI_NAMES(N_ALVEOLI),
     '            '>>Element group names must be unique',ERROR,*9999)
              ENDDO !i
              GROUP=.FALSE.
              DO nogrel=1,NTGREL !checks this is an element group
                CALL CUPPER(LAGREL(nogrel),LABEL)
                CALL STRING_TRIM(LABEL,IBEG,IEND)
                IF(ALVEOLI_NAMES(N_ALVEOLI).EQ.LABEL(IBEG:IEND))
     '            GROUP=.TRUE.
              ENDDO !nogrel
              CALL ASSERT(GROUP,
     '          '>>This is not an existing element group',ERROR,*9999)
              CDATA(2)='ELEMENTS' !puts element group list in NELIST
              CALL PARSILG(NELIST,NEM,CDATA(2),CDATA(1),ERROR,*9999)
              NETALV=NETALV+NELIST(0) !total ne's in all alveolar groups
C... calculate default values for alveolar centre
              DO nj=1,NJT 
                centre(nj)=0.d0 !initialise 
                centre2(nj)=0.d0 
              ENDDO
              DO noelem=1,NELIST(0)
                ne=NELIST(noelem)
                DO nn=2,NNT(nb_host) !1st node repeated in 2D elem
                  np=NPNE(nn,nb_host,ne) !node #'s of ne
                  DO nj=1,NJT
                    centre(nj)=centre(nj)+XP(1,1,nj,np)/
     '                (DBLE(NNT(nb_host))-1.d0)
                  ENDDO !nj
                ENDDO !nn
              ENDDO !noelem
              DO nj=1,NJT
                centre(nj)=centre(nj)/NELIST(0) !alveolar centre
              ENDDO
C... to determine which nodes are at entrance of alveoli
              IF(VORONOI_ALV) THEN !calculates centre & entrance plane
                CALL ALVPARA(nb_host,NELIST,NLL,NPL,NPLIST,NPNE,NXI,
     '            centre,centre2,XP,ERROR,*9999)
              ENDIF 
              DO nj=1,NJT
                RDEFLT(1)=(centre(nj)+centre2(nj))/2.d0
                !default value for alv centre
                WRITE(CHAR3,'(D12.4)') RDEFLT(1)
                WRITE(CHAR2,'(I1)') nj
                CALL STRING_TRIM(CHAR2,IBEG,IEND)
                FORMAT='($,'' Enter X('//CHAR2(IBEG:IEND)//') co-ord '
     '            //' for alveolar centre ['//CHAR3//']: '',D8.4)  '
                IF(IOTYPE.EQ.3) RDATA(1)=centre(nj)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) C(N_ALVEOLI,nj)=RDATA(1)
              ENDDO !nj
C... To calculate surface area and volume of alveolus
              IF(.NOT.VORONOI_ALV) THEN
                DO nj=1,NJT
                  centre2(nj)=centre(nj)
                ENDDO
              ENDIF
              vol=0.d0
              SA_total=0.d0
              depth=1.d-6 !initialise
              DO noelem=1,NELIST(0)
                ne=NELIST(noelem)
                length1=0.d0
                length2=0.d0
                length3=0.d0
C... Calculating alveolar volume and (maximum) alveolar depth             
                DO nj=1,NJT !use centre2 b/c this is on entrance plane
                  V1(nj)=XP(1,1,nj,NPNE(2,nb_host,ne))-centre2(nj)
                  length1=length1+V1(nj)**2.d0
                  V2(nj)=XP(1,1,nj,NPNE(3,nb_host,ne))-centre2(nj)
                  length2=V2(nj)**2.d0
                  V3(nj)=XP(1,1,nj,NPNE(4,nb_host,ne))-centre2(nj)
                  length3=V3(nj)**2.d0
                ENDDO !nj
                length1=DSQRT(length1)
                length2=DSQRT(length2)
                length3=DSQRT(length3)
                IF(length1.GT.depth) THEN
                  depth=length1
C                  np_max=NPNE(2,nb_host,ne)
                ENDIF
                IF(length2.GT.depth) THEN
                  depth=length2
C                  np_max=NPNE(3,nb_host,ne)
                ENDIF
                IF(length3.GT.depth) THEN
                  depth=length3
C                  np_max=NPNE(4,nb_host,ne)
                ENDIF  
                CALL CROSS(V2,V3,V4) !volume=1/6*(a.(bxc))
                vol=vol+(1.d0/6.d0)*DABS(DOT_PROD(V1,V4)) !alveolar volume
C... For tetrahedral elements, sum volume of 2 triangles
                IF(NPNE(1,nb_host,ne).NE.NPNE(2,nb_host,ne)) THEN
                  DO nj=1,NJT 
                    V1(nj)=XP(1,1,nj,NPNE(2,nb_host,ne))-centre2(nj)
                    V2(nj)=XP(1,1,nj,NPNE(3,nb_host,ne))-centre2(nj)
                    V3(nj)=XP(1 ,1,nj,NPNE(1,nb_host,ne))-centre2(nj)
                  ENDDO !nj                  
                  CALL CROSS(V2,V3,V4)
                 vol=vol+(1.d0/6.d0)*DABS(DOT_PROD(V1,V4)) 
                ENDIF
C... Calculating alveolar surface area                
                DO nj=1,NJT
                  V1(nj)=XP(1,1,nj,NPNE(3,nb_host,ne))-
     '              XP(1,1,nj,NPNE(2,nb_host,ne))
                  V2(nj)=XP(1,1,nj,NPNE(4,nb_host,ne))-
     '              XP(1,1,nj,NPNE(2,nb_host,ne))
                ENDDO !nj
                CALL CROSS(V1,V2,V4)
                SA=0.d0
                DO nj=1,NJT
                  SA=SA+V4(nj)**2.d0
                ENDDO
                SA=DSQRT(SA)/2.d0 !surface area
                SA_total=SA_total+SA !sums alveolar surface area
C... For tetrahedral elements, sum surface area of 2 triangles
                IF(NPNE(1,nb_host,ne).NE.NPNE(2,nb_host,ne)) THEN
                  DO nj=1,NJT
                    V1(nj)=XP(1,1,nj,NPNE(1,nb_host,ne))-
     '                XP(1,1,nj,NPNE(2,nb_host,ne))
                  ENDDO !nj
                  CALL CROSS(V1,V2,V4)
                  SA=0.d0
                  DO nj=1,NJT
                    SA=SA+V4(nj)**2.d0
                  ENDDO
                  SA=DSQRT(SA)/2.d0 !surface area
                  SA_total=SA_total+SA !sums alveolar surface area
                ENDIF
              ENDDO !noelem
              WRITE(OP_STRING,'('' Group name:'',A,'
     '          //''' Volume='',D12.4,'
     '          //''' Surface Area (SA)='',D12.4,'
     &          //''' Max depth='',D12.4)')
     '          CDATA(1),vol,SA_total,depth
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' SA-Vol ratio='',D12.4)')
     '          SA_total/(vol**(2.d0/3.d0)) !from Weibel:1963
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              vol_total=vol_total+vol !sums alveolar volumes
              SA_vol=SA_vol+SA_total/(vol**(2.d0/3.d0)) !sums SA/vol ratio
              SA_sum=SA_sum+SA_total
              depth_sum=depth_sum+depth
              GOTO 802 !continue until EXIT
            ENDIF !CDATA
C... determine scale_factor from cmiss units to alveolar size
C...        from Weibel:1963 (pg 68) average alveolar volume .011224mm
C...        **3, smaller value=0.00544d0 
            vol=vol_total/N_ALVEOLI !average alveolar volume
            SA_vol=SA_vol/N_ALVEOLI !average surface area to volume ratio
            SA_sum=SA_sum/N_ALVEOLI !average surface area
            depth_sum=depth_sum/N_ALVEOLI
            WRITE(OP_STRING,'('' Average SA-Vol ratio='',D12.4)')
     '        SA_vol
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     &        '('' Average SA='',D12.4,''Average Vol'',D12.4)') SA_sum,
     &        vol
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Average (maximum) depth='',D12.4)')
     '        depth_sum
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C           scale_factor=(0.011224d0/vol)**(1.d0/3.d0)
C           scale_factor=(0.00544d0/vol)**(1.d0/3.d0)
            scale_factor=(0.00744d0/vol)**(1.d0/3.d0)
            CALL ASSERT(N_ALVEOLI.NE.0,
     '        '>>No alveolar element groups defined',ERROR,*9999)
          ENDIF !CAP_MESH.EQ.3
          FIRST_PROJECT=.TRUE. !for 1st time through projection code
          CALL GENPCAP(INLET,nb_cap,nb_host,NBJ,NEELEM,NELIST,NELIST2,
     '      NENP,NKJ,NKJE,NP_INTERFACE,NPLIST,NPNE,NPNODE,NRE,nr,
     '      nr_feed,nr_host,NVJE,NVJP,NXI,OUTLET,C,SE,XP,REGULAR,
     '      ERROR,*9999)
        ENDIF !CAP_MESH.NE.1
C KSB(end) MAY 2002
        CALL CLOSEF(IOFILE2,ERROR,*9999)
        CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      ENDIF

      ENDIF !AIRWAY.OR.CIRCULATION

      IF(NET(0).EQ.0)THEN !no elements defined yet, initialise NXI,NENP
        DO ni=-NIM,NIM
          DO nj=0,NEIM
            DO ne=1,NEM
              NXI(ni,nj,ne)=0
            ENDDO !ne
          ENDDO !i
        ENDDO !ni
        DO np=1,NPM
          NENP(np,0,nr)=0
        ENDDO !np
      ENDIF
      IF(IOTYPE.NE.3)THEN !not writing out the file
        IF(AIRWAY) THEN
          DO i=-NIM,NIM !initialise NXI for elements
            DO j=0,NEIM ! that will be created
              NXI(i,j,0)=0
              DO noelem=NET(0)+1,NEM
                NXI(i,j,noelem)=0
              ENDDO !noelem
            ENDDO !j
          ENDDO !i
          CALL GNMESHR(nb_airway,NBJ,NEELEM,NELIST2,NENP,NKJ,NKJE,NORD,
     &      NP_INTERFACE,NPNE,NPNODE,nr,NRE,Nrefine,NSTORE,NVJE,
     &      NVJP,NXI,CE,RSTORE,SE,Spread,XP,XSTORE,MESH_TYPE,
     &      ERROR,*9999)
          IF(REDUCE)THEN
            CALL MESH_REDUCE(NBJ,NEELEM,NENP,NNB,NORD,NP_INTERFACE,NPNE,
     '        NPNODE,nr,NXI,TERMINAL_ORDER,ORDER,ERROR,*9999)
          ENDIF
        ENDIF !AIRWAY
      ENDIF

      CALL EXITS('IPMESH2')
      RETURN
 9999 CALL ERRORS('IPMESH2',ERROR)
      CALL EXITS('IPMESH2')
      RETURN 1
      END


