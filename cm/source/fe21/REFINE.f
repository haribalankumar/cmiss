      SUBROUTINE REFINE(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NEL,NELIST,NELIST3,NENP,NFF,NFFACE,NHE,NHP,NKEF,NKH,NKHE,NKJ,
     '  NKJE,NLF,NLL,NLLINE,NNB,NNF,NNL,NONY,NPF,NP_INTERFACE,NPL,
     '  NPLIST,NPNE,NPNF,NPNODE,NPNY,NRE,NRLIST,NUNK,NVHE,NVHP,NVJE,
     '  NVJF,NVJL,NVJP,NW,NWP,NXI,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '  CE,CONY,CYNO,DF,DL,PG,RG,SE,
     '  SF,SP,WG,XA,XE,XG,XP,YP,ZA,FIX,STRING,ERROR,*)
C LKC 5-NOV-97 unused ,NPLIST1

C#### Subroutine: REFINE
C###  Description:
C###    REFINE subdivides current mesh (geometry refinement only)=.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'mesh00.cmn'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NELIST3(0:NEM),NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),NRLIST(0:NRM),
     '  NUNK(NKM,NJM,NPM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),NVJP(NJM,NPM),NW(NEM,3,NXM),
     '  NWP(NPM,2),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CE(NMM,NEM,NXM),CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),DF(NFM),DL(3,NLM),
     '  PG(NSM,NUM,NGM,NBM),
     '  RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),SP(NKM,NBFM,NPM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER I,IBEG,IBEG4,IDRN,IEND,IEND4,IFROMC,MATCH,N3CO,nb,NB1,NB2,
     '  ne,NE2,NELIST_OLD,NELIST2(0:2),NET_OLD,ni,NILIST(10),nj,nl,
     '  nnne2,nn1,nn2,noelem,noelem2,nolist,none,noneli,nonr1,nonr2,
     '  notimes,NP1,NP2,NPSTART,nr,NR1,NR2,NR_REFINE,NRLIST2(0:2),
     '  NTLIST,NTTIMES,num_ne_ref,num_times,NUMTIMES,nx,XI1,XI2
      INTEGER*4 WORK_PTR
      REAL*8 MIN_X,PXI,RFROMC,X,XI,XI_SAMPLE(3)
      CHARACTER CHAR4*4
      LOGICAL ALL_REGIONS,BEM_MATRIX,CBBREV,EXIT,FOUND_NE,NE_INTERFACE,
     '  NOCROSS,NOINTERFACE,ONE_D,SECTOR,SIMPLEX

C LKC 5-NOV-97 unused INTEGER NPLIST1(0:NPM),

      CALL ENTERS('REFINE',*9999)

      IN_REFINE=.TRUE.

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR4(1:4),'(I4)') NPT(0)+1
        CALL STRING_TRIM(CHAR4,IBEG4,IEND4)

C---------------------------------------------------------------------

C#### Command: FEM refine
C###  Parameter:      <xi DIRECTION#s[1]>
C###    Defines the direction in which the element is refined
C###  Parameter:      <element (#s/all)[all]>
C###    Specifies the elements to be refined
C###  Parameter:      <at XI#[0.5]>
C###    Specifies the position of the refinement along the xi-axis
C###  Parameter:      <times #[1]>
C###    Specifies the number of refinements to make
C###  Parameter:      <ntimes #[1]>
C###    Specifies the number of cumulative refinements to make
C###  Parameter:      <from NODE#[1]>
C###    Specifies the node number assigned to the next node created by
C###    the refinement. By default the node number will be the next
C###    number after the largest in the file.
C###  Parameter:      <(cross/nocross)[cross]>
C###    Specifies whether cross derivatives are to be calculated
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:      <(interface/nointerface)[nointerface]>
C###  Description:
C###    Refine specified elements in specified xi-direction at
C###    specified xi location by specified number of times from
C###    specified node number.
C###  See-Example: 113
C###  Parameter:      <1dimension>
C###  Description:
C###    Refine for a region with 1-D basis functions only.

        OP_STRING(1)=STRING(1:IEND)//' <xi DIRECTION#s[1]>'
        OP_STRING(2)=BLANK(1:13)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:13)//'<at XI#[0.5]>'
        OP_STRING(4)=BLANK(1:13)//'<times #[1]>'
        OP_STRING(5)=BLANK(1:13)//'<ntimes #[1]>'
        OP_STRING(6)=BLANK(1:13)//'<from NODE#['
     '    //CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(7)=BLANK(1:13)//'<(cross/nocross)[cross]>'
        OP_STRING(8)=BLANK(1:13)//'<region (#s/all)[1]>'
        OP_STRING(9)=BLANK(1:13)//'<(interface/nointerface)'
     '    //'[nointerface]>'
        OP_STRING(10)=BLANK(1:13)//'<1dimension>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','REFINE',ERROR,*9999)
      ELSE

        CALL ASSERT(CALL_BASE,'Define a basis first',ERROR,*9999)
        CALL ASSERT(CALL_ELEM,'Define elements first',ERROR,*9999)

        nx=1 !temporary

        IF(ISC_GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GKM*USE_SPARSE,0,INTTYPE,
     '      ISC_GK_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GK_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GKM*USE_SPARSE,0,INTTYPE,
     '      ISR_GK_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISC_GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISC_GQM*USE_SPARSE,0,INTTYPE,
     '      ISC_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF
        IF(ISR_GQ_PTR(nx).EQ.0) THEN
          CALL ALLOCATE_MEMORY(NISR_GQM*USE_SPARSE,0,INTTYPE,
     '      ISR_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),10,NTLIST,NILIST,ERROR,*9999)
        ELSE IF(CALL_DEREFI) THEN
          NTLIST=1
          NILIST(1)=REFINE_DIR
        ELSE
          NTLIST=1
          NILIST(1)=1
        ENDIF
        IF(CBBREV(CO,'1DIMENSION',3,noco+1,NTCO,N3CO)) THEN
          ONE_D=.TRUE.
        ELSE
          ONE_D=.FALSE.
        ENDIF

        NPLIST(0)=0

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(CALL_DEREFI) THEN
          nr=1 !temp
          NELIST(0)=0
          DO none=1,NEELEM(0,nr)
            ne=NEELEM(none,nr)
            nj=NJ_LOC(NJL_FIEL,1,nr)
            nb=NBJ(nj,ne)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            IF(REFINE_FIELD_SAMPLE_POINT.EQ.1) THEN
              Xi_SAMPLE(1)=0.5d0
              Xi_SAMPLE(2)=0.5d0
              IF(NIT(nb).EQ.3) THEN
                Xi_SAMPLE(3)=0.5d0
              ELSE
                Xi_SAMPLE(3)=0.0d0
              ENDIF
              X=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,1,Xi_SAMPLE,XE(1,nj))
            ELSE IF(REFINE_FIELD_SAMPLE_POINT.EQ.2) THEN
              MIN_X=RMAX
C KAT 11Dec98: XI3_DIR not used
C              IF(NIT(nb).EQ.3) THEN
C                XI3_DIR=1
C              ELSE
C                XI3_DIR=0
C              ENDIF
              DO XI1=0,1
                DO XI2=0,1
C                  DO XI3=0,XI3_DIR
                    Xi_SAMPLE(1)=DBLE(XI1)
                    Xi_SAMPLE(2)=DBLE(XI2)
C                    Xi_SAMPLE(3)=DBLE(XI3)
                    Xi_SAMPLE(3)=0.0d0
                    X=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,1,Xi_SAMPLE,XE(1,nj))
                    IF(X.LT.MIN_X) MIN_X=X
                  ENDDO !XI1
                ENDDO !XI2
C              ENDDO !XI3
            ENDIF
            IF(X.LE.REFINE_FIELD_VALUE) THEN
              NELIST(0)=NELIST(0)+1
              NELIST(NELIST(0))=ne
            ENDIF
          ENDDO !none
        ELSE
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)

C LKC 1-JUL-1999 Check elements exist
          CALL ASSERT(NELIST(0).GT.0,
     '      '>> No Elements to refine',ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
          XI=RFROMC(CO(N3CO+1))
        ELSE IF(CALL_DEREFI) THEN
          XI=REFINE_XI
        ELSE
          XI=0.5d0
        ENDIF

        IF(CBBREV(CO,'TIMES',1,noco+1,NTCO,N3CO)) THEN
          NTTIMES=IFROMC(CO(N3CO+1))
        ELSE
          NTTIMES=1
        ENDIF

        IF(CBBREV(CO,'NTIMES',1,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(NTTIMES.EQ.1,'>>TIMES must be equal to 1',
     '      ERROR,*9999)
          NUMTIMES=IFROMC(CO(N3CO+1))
        ELSE
          NUMTIMES=1
        ENDIF

        IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
          NPSTART=IFROMC(CO(N3CO+1))
        ELSE
          NPSTART=NPT(0)+1
        ENDIF

        IF(CBBREV(CO,'NOCROSS',3,noco+1,NTCO,N3CO)) THEN
C KAT 14Dec98: The process of determining cross derivatives has been
C         changed so I think it is now stable.  I have therefore
C         changed the default to cross.
          NOCROSS=.TRUE.
!         !Added AJP 14-6-93.  If nocross is true then a refined
!         !hermite element will have the cross derivatives set to
!         !zero at the new nodes.  This is to keep the refinement
!         !stable (see CMISS example ??? to illustrate this).
        ELSE
          NOCROSS=.FALSE.
        ENDIF

        IF(CBBREV(CO,'INTERFACE',2,noco+1,NTCO,N3CO)) THEN
          NOINTERFACE=.FALSE.
!         !Added AJP 12-7-96.
!         !If NOINTERFACE is .true. then no attempt is made to
!         !provide a consistent coupling between elements across
!         !an interface between different regions.  If NOINTERFACE is
!         !.false. then an ATTEMPT is made to ensure that when an
!         !element has been refined and is on an interface between 2
!         !regions the elements on the other side of the interface
!         !in the other region are also refined.  NOTE: this does not
!         !work very well.  It was designed (and should work) in the
!         !situation where a BE region is refined and is surrounded by
!         !an FE region (and the FE region has elements only 1 layer
!         !thick).
        ELSE
          NOINTERFACE=.TRUE.
        ENDIF

        CALL ASSERT(NET(0)+NTTIMES*NTLIST*NELIST(0).LE.NE_R_M,
     '    '>>NE_R_M too small',ERROR,*9999)

        DO num_times=1,NUMTIMES !will only do next loop once
          num_ne_ref=0
          DO notimes=1,NTTIMES !Loop over #times
            DO nolist=1,NTLIST !Loop over #Xi-directions
              NET_OLD=NET(0) !Total number of elements before refining
              IDRN=NILIST(nolist)
              DO noneli=1,NELIST(0) !Loop over elements
                ne=NELIST(noneli)
                num_ne_ref=num_ne_ref+1
                NELIST3(num_ne_ref)=ne
                NR_REFINE=NRE(ne) !is region# for current element
                NE_INTERFACE=.FALSE.
                IF(.NOT.NOINTERFACE) THEN
 !Check if ne is an interface element.
                  NB1=NBJ(1,ne)
                  nn1=1
                  EXIT=.FALSE.
                  DO WHILE((nn1.LE.NNT(NB1)).AND.(.NOT.EXIT))
                    NP1=NPNE(nn1,NB1,ne)
 !Check all nodes on element ne to see if they are
 !shared by the same region(s).
                    IF(NP_INTERFACE(NP1,0).EQ.1)THEN
                      EXIT=.TRUE.
                    ENDIF
                    IF(.NOT.EXIT)THEN
                      IF(nn1.GT.1)THEN
                        nn2=nn1-1 !Check previous node regions
                        NP2=NPNE(nn2,NB1,ne)
                        DO nonr1=1,NP_INTERFACE(NP1,0)
 !Check all regions
                          NR1=NP_INTERFACE(NP1,nonr1)
                          nonr2=1
                          NR2=NP_INTERFACE(NP2,nonr2)
                          DO WHILE((nr1.ne.NR2).AND.(.NOT.EXIT))
                            nonr2=nonr2+1
                            IF(nonr2.GT.NP_INTERFACE(NP2,0))THEN
                              EXIT=.TRUE.
                            ELSE
                              NR2=NP_INTERFACE(NP2,nonr2)
                            ENDIF
                          ENDDO
                        ENDDO
                      ENDIF
                      IF((nn1.EQ.NNT(NB1)).AND.(.NOT.EXIT))THEN
                        NE_INTERFACE=.TRUE.
                      ENDIF
                    ENDIF
                    nn1=nn1+1
                  ENDDO
                ENDIF !nointerface
                NELIST2(0)=0
                NRLIST2(0)=0
                IF(NE_INTERFACE)THEN
 !                 !Find other element number(s) containing the same
 !                 !interface nodes as element ne.
                  NP1=NPNE(1,NB1,ne)
                  DO nonr1=1,NP_INTERFACE(NP1,0)
                    NR2=NP_INTERFACE(NP1,nonr1)
 !   New AJP 20-1-94
                    IF(nr2.ne.NR_REFINE)THEN
                      FOUND_NE=.FALSE.
                      noelem2=0
                      DO WHILE(.NOT.FOUND_NE)
                        noelem2=noelem2+1
                        CALL ASSERT(noelem2.LE.NEELEM(0,NR2),
     '                    '   >>Error in finding element',ERROR,*9999)
                        NE2=NEELEM(noelem2,NR2)
                        NB2=NBJ(1,NE2)
                        MATCH=0
                        DO nn2=1,NNT(NB2)
                          NP2=NPNE(nn2,NB2,NE2)
                          IF(NP1.EQ.NP2)THEN
 !Interface node NP1 is in element NE2.
 !Check other nodes in ne are also in NE2
                            MATCH=1
                            DO nn1=2,NNT(NB1)
 !NP1 is node at node 1 on ne
                              DO nnne2=1,NNT(NB2)
                                IF(NPNE(nnne2,NB2,NE2).EQ.
     '                            NPNE(nn1,NB1,ne))MATCH=MATCH+1
                              ENDDO
                            ENDDO
                          ENDIF
                        ENDDO
                        IF(MATCH.EQ.NNT(NB1))THEN
                          FOUND_NE=.TRUE.
                          NELIST2(0)=NELIST2(0)+1
                          NELIST2(NELIST2(0))=NE2
                          NRLIST2(NELIST2(0))=NR2
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDDO
                ENDIF !Found all other elements  - NELIST2

                IF(DOP)THEN
                  IF(NE_INTERFACE)THEN
                    WRITE(OP_STRING,*)' ne=',ne
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,*)' Number of interface elements=',
     '                NELIST2(0)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,*)' Interface elements : ',
     '                (NELIST2(I),I=1,NELIST2(0))
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,*)' Region numbers : ',
     '                (NRLIST2(I),I=1,NELIST2(0))
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF

c cpb 27/11/96 Adding sectors
 ! Decide if element ne is a special hermite simplex or not.
                nl=NLL(2,ne) !Second line number of element ne
                SECTOR=.FALSE.
                DO ni=1,NIT(NBJ(1,ne))
                  IF(IBT(1,ni,NBJ(1,ne)).EQ.5.OR.
     '              IBT(1,ni,NBJ(1,ne)).EQ.6) SECTOR=.TRUE.
                ENDDO !ni
                IF(SECTOR) THEN
                  CALL DIVS2(IBT,IDO,IDRN,INP,NBH,NBJ,ne,NEELEM,
     '              NHE(1,nx),NKJE,NKJ,NPF,NPLIST,NPNE,NPNODE,NPSTART,
     '              NR_REFINE,NRE,NUNK,NVJE,NVJP,NW(1,1,nx),nx,
     '              CE(1,1,nx),SE,SP,XA,XE,XI,XP,NOCROSS,ERROR,*9999)
                ELSE IF(nl.GT.0) THEN
                  IF((NPL(1,1,nl).EQ.6).OR.(NPL(1,1,nl).EQ.7)) THEN
 !                   !A special hermite simplex
 !                   !Refine current element ne
                    CALL DIVS1(IBT,IDO,IDRN,INP,NBH,NBJ,
     '                ne,NEELEM,NHE(1,nx),NKJE,NKJ,NPF,NPNE,NPNODE,
     '                NPSTART,NR_REFINE,NRE,NVJE,NVJP,NW(1,1,nx),nx,
     '                CE(1,1,nx),SE,XA,XE,XI,XP,NOCROSS,ERROR,*9999)
                  ELSE
 !                   !Refine current element ne
                    CALL DIVH1(IBT,IDO,IDRN,INP,NBH,NBJ,ne,NEELEM,
     '                NHE(1,nx),NKJE,NKJ,NPF,NPLIST,NPNE,NPNODE,
     '                NPSTART,NR_REFINE,NRE,NUNK,NVJE,NVJP,NW(1,1,nx),
     '                nx,CE(1,1,nx),SE,SP,XA,XE,XI,XP,NOCROSS,ERROR,
     '                *9999)
                  ENDIF
                ELSE
 !                 !Refine current element ne
                  CALL DIVH1(IBT,IDO,IDRN,INP,NBH,NBJ,ne,NEELEM,
     '              NHE(1,nx),NKJE,NKJ,NPF,NPLIST,NPNE,NPNODE,NPSTART,
     '              NR_REFINE,NRE,NUNK,NVJE,NVJP,NW(1,1,nx),nx,
     '              CE(1,1,nx),SE,SP,XA,XE,XI,XP,NOCROSS,ERROR,*9999)
                ENDIF
                num_ne_ref=num_ne_ref+1
                NELIST3(num_ne_ref)=NET(0) !record new element #
 !Refine elements in other regions sharing same nodes
                DO noelem=1,NELIST2(0)
                  ne2=NELIST2(noelem)
                  nr2=NRLIST2(noelem)
                  nl=NLL(2,NE2) !Second line number of element ne
                  SECTOR=.FALSE.
                  SIMPLEX=.FALSE.
                  IF(NIT(NBJ(1,ne2)).GT.1) THEN
                    DO ni=1,NIT(NBJ(1,ne2))
                      IF(IBT(1,ni,NBJ(1,ne2)).EQ.5.OR.
     '                  IBT(1,ni,NBJ(1,ne2)).EQ.6)SECTOR=.TRUE.
                    ENDDO !ni
                    IF((nl.GT.0).AND.((NPL(1,1,nl).EQ.6).OR.
     '                (NPL(1,1,nl).EQ.7))) SIMPLEX=.TRUE.
                  ENDIF
                  IF(SECTOR) THEN
                    CALL DIVS2(IBT,IDO,IDRN,INP,NBH,NBJ,ne2,NEELEM,
     '                NHE(1,nx),NKJE,NKJ,NPF,NPLIST,NPNE,NPNODE,
     '                NPSTART,nr2,NRE,NUNK,NVJE,NVJP,NW(1,1,nx),nx,
     '                CE(1,1,nx),SE,SP,XA,XE,XI,XP,NOCROSS,ERROR,
     '                *9999)
                  ELSE IF(SIMPLEX)THEN
                    CALL DIVS1(IBT,IDO,IDRN,INP,NBH,NBJ,
     '                ne2,NEELEM,NHE(1,nx),NKJE,NKJ,NPF,NPNE,
     '                NPNODE,NPSTART,nr2,NRE,NVJE,NVJP,NW(1,1,nx),nx,
     '                CE(1,1,nx),SE,XA,XE,XI,XP,NOCROSS,ERROR,*9999)
                  ELSE
                    CALL DIVH1(IBT,IDO,IDRN,INP,NBH,NBJ,ne2,NEELEM,
     '                NHE(1,nx),NKJE,NKJ,NPF,NPLIST,NPNE,NPNODE,
     '                NPSTART,nr2,NRE,NUNK,NVJE,NVJP,NW(1,1,nx),nx,
     '                CE(1,1,nx),SE,SP,XA,XE,XI,XP,NOCROSS,ERROR,
     '                *9999)
                  ENDIF
                  num_ne_ref=num_ne_ref+1
                  NELIST3(num_ne_ref)=NET(0) !record new element #
                ENDDO
 ! Update interface info
                CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
              ENDDO !elements

              CALL ASSERT(NET(0).LE.NEM,'>>Increase NEM',ERROR,*9999)
              CALL ASSERT(NPNODE(0,0).LE.NPM,'>>Increase NPM',
     '          ERROR,*9999)
              DO nr=1,NRT
                CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase NP_R_M',
     '            ERROR,*9999)
              ENDDO !nr

c cpb 20/7/95 NENP needs to be calculated before the lines are
c calculated in LINSEG
              IF(ONE_D)THEN
                !for 1-d basis function over entire region
                !should use something better than "NBJ(1,1)"
                CALL CALC_NENP_1D(NBJ(1,1),NEELEM,NENP,NPNE,NRLIST(1),
     '            ERROR,*9999)
              ELSE
                CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)
              ENDIF

              NB1=NBJ(1,NEELEM(1,NR_REFINE))
              IF(NBI(NB1).EQ.2.OR.NBI(NB1).EQ.3) THEN
 ! scale factors based on specified element or global derivs
C ***           Transfer SE to DL
                CALL LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '            NLLINE,NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*9999)
                CALL SEDL(IBT,IDO,NB1,NEELEM,NLL,NNL,NPL,DL,SE,
     '            ERROR,*9999)
              ENDIF

C             C CPB Adding average arc-length scaling so must determine
C             NXI before SE
              IF(ONE_D)THEN
                !for 1-d basis function over entire region
                !should use something better than "NBJ(1,1)"
                CALL NENXI_1D(NBJ(1,1),NEELEM,NENP,NPNE,NRLIST(1),NXI,
     '            ERROR,*9999)
                CALL LINCAL_1D(NBJ,NEELEM,NKJE,NRLIST(1),NRE,NVJE,
     '            SE,ERROR,*9999)
              ELSE
                CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
                CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '            NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '            DL,SE,XP,ERROR,*9999)
                CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '            NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,
     '            NRE,NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,
     '            XE,XG,XP,ERROR,*9999)
              ENDIF

              CALL ASSERT(NET(1).LE.NEM,'>>NEM too small',ERROR,*9999)
              CALL ASSERT(NFT   .LE.NFM,'>>NFM too small',ERROR,*9999)
              CALL ASSERT(NLT   .LE.NLM,'>>NLM too small',ERROR,*9999)
              CALL ASSERT(NPT(1).LE.NPM,'>>NPM too small',ERROR,*9999)

C Need to update NELIST so that new elements are refined in other xi
C directions (if required)
              NELIST_OLD=NELIST(0)
              NELIST(0)=NELIST(0)+NET(0)-NET_OLD
              DO noneli=NELIST_OLD+1,NELIST(0)
                NELIST(noneli)=noneli !new elems are in sequntl order
              ENDDO
            ENDDO !End of xi directions
          ENDDO !End of Ntimes loop
          NELIST(0)=num_ne_ref
          DO num_ne_ref=1,NELIST(0)
            NELIST(num_ne_ref)=NELIST3(num_ne_ref)
          ENDDO !num_ne_ref
        ENDDO !end of num_times loop

C cpb 5/3/96 Check if any dependent varaible information to refine
C SAB 23/01/01 Don't do this for a fit problem cause the routines are
C       not set up to handle a fit.
        IF(CALL_EQUA.OR.CALL_OPTI) THEN
C         New AJP 23/3/95 Setting up dependent variable info
          BEM_MATRIX=.FALSE.
          DO nr=1,NRT

            IF(ITYP4(nr,nx).EQ.2) BEM_MATRIX=.TRUE.

C*** Calculate NHP
            CALL CALC_NHP(NBH,NEELEM,NHE(1,nx),NHP(1,0,nx),NPNE,
     '        NPNODE,nr,nx,ERROR,*9999)

C*** Calculate NVHP,NVHE
            CALL CALC_VERSIONS_DEP(NBH,NBJ,NEELEM,NHE(1,nx),
     '        NHP(1,nr,nx),NPNE,NPNODE,nr,NVHE,NVHP,NVJE,NVJP,
     '        nx,ERROR,*9999)

C*** Calculate NKH
            IF(IOTYPE.NE.3) THEN
              IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN
C               Don't do for cardiac cellular modelling
              ELSE
                CALL CALC_NKH(NBH,NEELEM,NHP(1,nr,nx),NKH,
     '            NPNE,NPNODE,nr,NW(1,1,nx),nx,ERROR,*9999)
              ENDIF
            ENDIF

C*** Calculate ny's etc
            IF(IOTYPE.NE.3) THEN
              IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN
C               Don't do for cardiac cellular modelling
              ELSE
                CALL CALC_NY_MAPS_DEP(NBH,NEELEM,
     '            NHP(1,0,nx),NKH,
     '            NP_INTERFACE,NPNODE,NPNY(0,1,0,nx),nr,NVHP,nx,
     '            NYNE,NYNP,NYNR(0,0,1,0,nx),ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO !nr

          IF(ITYP4(NRLIST(1),nx).GT.0) THEN !A solve option has been set

C CPB 6/11/95 Temporary work array allocation
            IF(KTYP24.NE.0) THEN
              WORK_PTR=0
              CALL ALLOCATE_MEMORY(NYT(1,1,nx)*NYT(2,1,nx),1,CHARTYPE,
     '          WORK_PTR,MEM_INIT,ERROR,*9999)
            ENDIF
            CALL CALC_SPARSE_SOLVE(NISC_GKM,NISR_GKM,
     '        %VAL(ISC_GK_PTR(nx)),%VAL(ISR_GK_PTR(nx)),
     '        LGE,NYT(1,1,nx),NYT(2,1,nx),NBH,1,NEELEM,NHE,
     '        NPNE,NPNY,NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZ_GK_M,KTYP24,
     '        %VAL(WORK_PTR),FIX(1,1,nx),.TRUE.,ERROR,*9999)
            IF(KTYP24.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
            IF(BEM_MATRIX) THEN !BEM
C CPB 6/11/95 Temporary work array allocation
              IF(KTYP24.NE.0) THEN
                WORK_PTR=0
                CALL ALLOCATE_MEMORY(NYT(1,2,nx)*NYT(2,2,nx),1,
     '            CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
              ENDIF
              CALL CALC_SPARSE_SOLVE(NISC_GQM,NISR_GQM,
     '          %VAL(ISC_GQ_PTR(nx)),%VAL(ISR_GQ_PTR(nx)),
     '          LGE,NYT(1,2,nx),NYT(2,2,nx),NBH,2,NEELEM,NHE,
     '          NPNE,NPNY,NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZ_GQ_M,
     '          KTYP24,%VAL(WORK_PTR),FIX(1,1,nx),.TRUE.,ERROR,*9999)
              IF(KTYP24.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
            ENDIF
          ENDIF

C End new AJP 23/3/95

        ENDIF

C KAT 18Oct00: This is now done in REFINE_SETNODE
C       C ***     Reduce angles to 0 to 2*pi
C       DO nr=1,NRT
C       IF(ITYP10(nr).GE.2) THEN
C       DO nonode=1,NPNODE(0,nr)
C       np=NPNODE(nonode,nr)
C       C GMH 8/1/97 Update cmgui link
C       CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C       IF(XP(1,1,2,np).GT.2.d0*PI) XP(1,1,2,np)=
C       '                                      XP(1,1,2,np)-2.d0*PI
C       IF(ITYP10(nr).GE.3) THEN
C       IF(XP(1,1,3,np).GT.2.d0*PI) XP(1,1,3,np)=
C       '                                        XP(1,1,3,np)-2.d0*PI
C       ENDIF
C       ENDDO
C       ENDIF
C       ENDDO

C***      new CS 15/4/98
C removing for the moment - should put a flag option on refine for this
        IF(JTYP2C.EQ.1) THEN
C         DO nr=1,NRT
          nr=1
          CALL HANGING_NODE_DETECT(IBT,IDO,IDRN,INP,NBJ,NEELEM,
     '      NELIST,NENP,NKJE,NPF,NPLIST,NPNE,NPNODE,nr,NVJE,
     '      NVJP,NWP,SE,SP,XA,XE,XP,ERROR,*9999)
C         ENDDO
        ENDIF

        CALL_INIT=.FALSE. !Initial conditions need to be set at new
        CALL_SOLV=.FALSE. !nodes

      ENDIF

      IN_REFINE=.FALSE.

      CALL EXITS('REFINE')
      RETURN
 9999 CALL ERRORS('REFINE',ERROR)
      CALL EXITS('REFINE')
      RETURN 1
      END


