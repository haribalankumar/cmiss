      SUBROUTINE UPORDR(NBJ,NEELEM,NELIST,NENP,NKJ,NKJE,NORD,NPLIST1,
     '  NPLIST2,NPNE,NPNODE,NRE,NRLIST,NVJE,NVJP,NXI,SE,XP,STRING,
     '  ERROR,*)

C#### Subroutine: UPORDR
C###  Description:
C###    UPORDR reorders a 1D tree network so that elements and nodes
C###    increase with order.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NPLIST1(0:NPM),NPLIST2(0:NPM),NORD(5,NE_R_M),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRE(NEM),
     '  NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,nb,NE_START,nr,N3CO,OFFSET
      INTEGER*4 NE_TEMP_PTR,NE_TEMP_OLD_PTR,NORD_TEMP_PTR,NPNE_TEMP_PTR,
     &  XP_TEMP_PTR
      LOGICAL ALL_REGIONS,CBBREV

      CALL ENTERS('UPORDR',*9999)

      IF(DOP)THEN
        WRITE(OP_STRING,'('' Reached call to UPORDR_DYNAM'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update ordering
C###  Parameter:      <start_element NE_START[1]>
C###  Specify the start element for reordering.
C###  Parameter:      <offset #[0]>
C###    Specify an offset to add to the nodes and elements.
C###  Parameter:      <region #[1]>
C###    Specify the region numbers to update.
C###  Description: updates 1D tree network ordering to give increasing
C###    node and element numbers with 'branch' order. The reordered tree
C###    has nodes and elements each starting from 1. Any fields in XP
C###    are also updated at the same time; these should be written out
C###    for using with the reordered tree.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<start_element NE_START[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE
        IF(CBBREV(CO,'START_ELEMENT',1,noco+1,NTCO,N3CO)) THEN
          NE_START=IFROMC(CO(N3CO+1))
        ELSE
          NE_START=1
        ENDIF
        CALL ASSERT(NXI(-1,0,NE_START).EQ.0,
     '    '>>Start element is not the tree stem',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL ASSERT(NRLIST(0).EQ.1,
     '    '>> Only implemented for 1 region',ERROR,*9999)
        nr=NRLIST(1)
        nb=NBJ(1,NE_START) !basis fn. for start element
        IF(DOP)THEN
          WRITE(OP_STRING,'('' Reached call to UPORDR_DYNAM'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(CBBREV(CO,'OFFSET',3,noco+1,NTCO,N3CO)) THEN
          OFFSET=IFROMC(CO(N3CO+1))
        ELSE
          OFFSET=0
        ENDIF

        NE_TEMP_PTR=0 !NE_TEMP(NEM)
        NE_TEMP_OLD_PTR=0 !NE_TEMP_OLD(NEM)
        NORD_TEMP_PTR=0
        NPNE_TEMP_PTR=0 !NPNE_TEMP(2,NEM)
        XP_TEMP_PTR=0 !XP_TEMP(NVM,NJM,NPM)
        CALL ALLOCATE_MEMORY(NEM,1,INTTYPE,NE_TEMP_PTR,.TRUE.,ERROR,
     &    *9999)
        CALL ALLOCATE_MEMORY(NEM,1,INTTYPE,NE_TEMP_OLD_PTR,.TRUE.,ERROR,
     &    *9999)
        CALL ALLOCATE_MEMORY(NE_R_M,1,INTTYPE,NORD_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(2*NEM,1,INTTYPE,NPNE_TEMP_PTR,.TRUE.,ERROR,
     &    *9999)
        CALL ALLOCATE_MEMORY(NVM*NJM*NPM,1,DPTYPE,XP_TEMP_PTR,.TRUE.,
     &    ERROR,*9999)

        CALL UPORDR_DYNAM(nb,NBJ,NEELEM,NELIST,NENP,NE_START,
     &    %VAL(NE_TEMP_PTR),%VAL(NE_TEMP_OLD_PTR),NKJ,NKJE,NORD,
     &    %VAL(NORD_TEMP_PTR),NPLIST1,NPLIST2,NPNE,%VAL(NPNE_TEMP_PTR),
     &    NPNODE,nr,NRE,NVJE,NVJP,NXI,OFFSET,SE,
     &    XP,%VAL(XP_TEMP_PTR),ERROR,*9999)

        CALL FREE_MEMORY(NE_TEMP_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NE_TEMP_OLD_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NORD_TEMP_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NPNE_TEMP_PTR,ERROR,*9999)
        CALL FREE_MEMORY(XP_TEMP_PTR,ERROR,*9999)

      ENDIF

      CALL EXITS('UPORDR')
      RETURN
 9999 CALL ERRORS('UPORDR',ERROR)
      CALL EXITS('UPORDR')
      RETURN 1
      END



