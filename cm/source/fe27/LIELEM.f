      SUBROUTINE LIELEM(IBT,NBH,NBJ,NEELEM,NELIST,
     '  NFF,NHE,NHP,NKH,NKHE,NKJE,NLL,NP_INTERFACE,
     '  NPF,NPLIST,NPNE,NPNODE,NRE,NRLIST,NVHE,NVHP,
     '  NVJE,NW,NXLIST,NYNE,NYNP,VOLTC,
     '  CURVCORRECT,PG,RG,SE,
     '  VOL,VOLT,WG,XA,XAB,XE,XG,XP,YP,
     '  ZA,ZE,ZG,ZP,STRING,ERROR,*)

C#### Subroutine: LIELEM
C###  Description:
C###    LIELEM lists element parameters NPNE(nn,nb,ne) and XE(ns,nj).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NFF(6,NEM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NLL(12,NEM),NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),VOLTC(NBFM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),
     '  VOL(NIM,NBFM),VOLT(NIM,NBFM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nh_update,nhx,nj,nj_volume,
     '  noelem,nolist,NP_average,nr,nr_interface,nx,nxc
      CHARACTER FILE*100,TYPE*15,TYPE1*12
      LOGICAL ALL_REGIONS,AVGENODE,CBBREV,OPFILE,OUTPUT_TOT_ONLY,
     '  ELEM_VOL

      CALL ENTERS('LIELEM',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list elements<;FILENAME> [undeformed]
C###  Description:
C###    Lists specified undeformed element(s) to screen
C###    or file FILENAME.opelem if qualifier present.
C###  Parameter:    <elements (GROUP/#s/all)[all]>
C###    Specify either element group, element numbers or all elements
C###    to be listed
C###  Parameter:    <(total/all)[all]>
C###    Specifying total will list all the element numbers and the
C###    total area of all elements. The all parameter will list each
C###    element, along with with global node, face and line information.
C###  Parameter:    <region (#s/all)[all]>
C###    Specify the region number of the elements to list, or all
C###    regions.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> [undeformed]'
        OP_STRING(2)=BLANK(1:15)//'<elements (GROUP/#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(total/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list elements<;FILENAME> deformed
C###  Description:
C###    Lists specified deformed element(s) to screen
C###    or file FILENAME.opelem if qualifier present.
C###  Parameter:    <elements (GROUP/#s/all)[all]>
C###    Specify either element group, element numbers or all elements
C###    to be listed
C###  Parameter:    <(total/all)[all]>
C###    Specifying total will list all the element numbers and the
C###    total area of all elements. The all parameter will list each
C###    element, along with with global node, face and line information.
C###  Parameter:    <average NODE_NUM#[1]
C###    The 'average NODE_NUM' option is used to calculate the average
C###    position of the given list of deformed nodes, and then
C###    change the deformed position for node NODE_NUM to that average
C###    for use in the computation of deformed element volumes.
C###    The node list and NODE_NUM must belong to the specified region.
C###  Parameter:       <in NH#[1]>
C###    Specifies the dependent variable (nh) number to use
C###  Parameter:       <node (GROUP/#s/all)[all]>>
C###    Specifies the node group, numbers or all nodes to use
C###  Parameter:    <region (#s/all)[all]>
C###    Specifies the region of the elements to list
C###  Parameter:    <class #[1]>
C###    Specfies the class number (of solve type) of the elements to
C###    list

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> deformed'
        OP_STRING(2)=BLANK(1:15)//'<elements (GROUP/#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(total/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<average NODE_NUM#[1]'
        OP_STRING(5)=BLANK(1:15)//'   <in NH#[1]>'
        OP_STRING(6)=BLANK(1:15)//'   <node (#s/GROUP/all)[all]>>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list elements<;FILENAME> groups
C###  Description:
C###    Lists all element groups to screen
C###    or file FILENAME.opelem if qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> groups'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list elements<;FILENAME> interface
C###  Parameter:    <nr_interface (#)[1]>
C###    Specify the region which the nodes should interface with
C###  Parameter:    <region (#s/all)[all]>
C###    Specify the region number of the elements to list, or all
C###    regions.
C###  Description:
C###    Lists all element which are defined by nodes which
C###    interface the regions specified by nr_interface and region.
C###    This is useful for exporting surfaces of bem meshs where
C###    multiple surfaces are defined for each region.


        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> interface'
        OP_STRING(2)=BLANK(1:15)//'<nr_interface (#)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list elements<;FILENAME> cavity_volume
C###  Description:
C###    Lists deformed elements specifically for
C###    uncoupled cavity problems.
C###  Parameter:    <region (#s/all)[all]>
C###    Specifies the region of the elements to list

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> cavity_volume'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIELEM',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opelem','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(CBBREV(CO,'ELEMENTS',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ELSE !list all elements in region list
          NELIST(0)=0
          DO nolist=1,NRLIST(0)
            nr=NRLIST(nolist)
            DO noelem=1,NEELEM(0,nr)
              NELIST(NELIST(0)+noelem)=NEELEM(noelem,nr)
            ENDDO !noelem
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO !nr
        ENDIF
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

C!!!    Temporary until people get used to command
        IF(CBBREV(CO,'NUMBER',1,noco+1,noco+3,N3CO)) THEN
          WRITE(OP_STRING,'(''>>Repeat command replacing the '
     '      //' option NUMBER with ELEMENT (and the #(s)'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GO TO 9999
        ENDIF

        IF(CBBREV(CO,'AVERAGE',1,noco+1,NTCO,N3CO)) THEN
          NP_average=IFROMC(CO(N3CO+1))
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          CALL ASSERT(NPLIST(0).GT.0,'>>No nodes in list!',ERROR,*9999)
          AVGENODE=.TRUE.
          IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '        ERROR,*9999)
            nhx=IFROMC(CO(N3CO+1))
            nh_update=NH_LOC(nhx,nx)
            CALL ASSERT(nh_update.GT.0,'>>nh # not defined',
     '        ERROR,*9999)
          ELSE
            nh_update=1
          ENDIF
        ELSE
          AVGENODE=.FALSE.
        ENDIF

        nr_interface=-1 !setting to an invalid number
        IF(CBBREV(CO,'GROUPS',1,noco+1,NTCO,N3CO)) THEN
          TYPE1='GROUPS'

C LKC 1-JUL-1999 List elem interface
        ELSEIF(CBBREV(CO,'INTERFACE',3,noco+1,NTCO,N3CO)) THEN
          TYPE1='INTERFACE'
          IF(CBBREV(CO,'NR_INTERFACE',3,noco+1,NTCO,N3CO)) THEN
            nr_interface=IFROMC(CO(N3CO+1))
          ELSE
            nr_interface=1
          ENDIF

C LKC 25-SEP-1999 Change the way interface works -
C   specify a specific region
C          CALL ASSERT(nr_interface.LE.3,
C     '      '>> NUM_REGIONS must be <= 3',ERROR,*9999)
C          CALL ASSERT(nr_interface.GE.0,
C     '      '>> NUM_REGIONS must be >= 0',ERROR,*9999)

          CALL ASSERT(NEELEM(0,nr_interface).GE.0,
     '      '>> No elements defined in this region',ERROR,*9999)


        ELSE
          TYPE1='ALL_ELEMENTS'
        ENDIF

        IF(CBBREV(CO,'DEFORMED',1,noco+1,NTCO,N3CO)) THEN
          TYPE='DEFORMED'
          CALL ASSERT(CALL_INIT,'>>No initial conditions defined',
     '      ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE
          nx=1
          TYPE='UNDEFORMED'
        ENDIF

        IF(CBBREV(CO,'CAVITY_VOLUME',1,noco+1,NTCO,N3CO)) THEN
          TYPE='DEFORMED_CAVITY'
          CALL ASSERT(CALL_INIT,'>>No initial conditions defined',
     '      ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'TOTAL_VOLUME',1,noco+1,NTCO,N3CO)) THEN
          OUTPUT_TOT_ONLY=.TRUE.
        ELSE
          OUTPUT_TOT_ONLY=.FALSE.
        ENDIF

        IF(CBBREV(CO,'VOLUME_FIELD',3,noco+1,NTCO,N3CO))THEN
          ELEM_VOL=.TRUE.
          nj=IFROMC(CO(N3CO+1))
          nj_volume=NEJ_LOC(nj,nr)
          write(*,*)'nj ',nj,' nej_loc ',NEJ_LOC(nj,nr),nr
          CALL ASSERT(nj_volume.GT.0,'>>Non-existent field',ERROR,*9999)
        ELSE
          ELEM_VOL=.FALSE.
        ENDIF

        IF(TYPE1(1:6).EQ.'GROUPS') THEN
          CALL OPELEMG(ERROR,*9999)

C old MPN 6Mar96
C        ELSE IF(TYPE1(1:11).EQ.'ONE_ELEMENT') THEN
C          nr=NRE(ne) !AJP 17-4-93
C          CALL ASSERT(ne.gt.0.and.NE.LE.NET(nr),
C     '      '>>Element number out of range',ERROR,*9999)
C          IF(TYPE(1:10).EQ.'UNDEFORMED')THEN
C            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,ne),
C     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
C         '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C            CALL OPELEM1(NBJ(1,ne),ne,NFF(1,ne),
C     '        NKE(1,1,1,ne),NLL(1,ne),NPNE(1,1,ne),nr,
C     '        NVJE(1,1,1,ne),PG,RG,SE(1,1,ne),VOL,VOLT,
C     '        WG,XE,XG,OUTPUT_TOT_ONLY,ERROR,*9999)
C          ELSE IF(TYPE(1:8).EQ.'DEFORMED')THEN
C            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
C     '        NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
C     '        ERROR,*9999)
C            CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),
C     '        NPF(1,ne),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),
C     '         nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
C     '        ERROR,*9999)
C            CALL OPELEMD(NBH(1,1,ne),NBJ(1,ne),
C     '        ne,NFF(1,ne),NHE(ne,nx),NKE(1,1,1,ne),
C     '        NLL(1,ne),NPNE(1,1,ne),
C     '        nr,NVHE(1,1,1,ne),nx,
C     '        PG,RG,SE(1,1,ne),VOL,VOLT,
C     '        WG,ZE,ZG,OUTPUT_TOT_ONLY,ERROR,*9999)
C          ENDIF
C        ELSE IF(TYPE1(1:12).EQ.'ALL_ELEMENTS') THEN


C LKC 1-JUL-1999 new interface option
        ELSEIF(TYPE1(1:9).EQ.'INTERFACE') THEN

          CALL OPELEM_INTERFACE(NBJ,NEELEM,NELIST,NP_INTERFACE,
     '      NPNE,nr_interface,NRLIST,ERROR,*9999)

        ELSE
C Adding nr loop to fix MPN bug AJP 14Apr97
CC new MPN 8Apr97: regions now handled within opelem
C          DO nolist=1,NRLIST(0)
C            nr=NRLIST(nolist)
            CALL OPELEM(IBT,NBJ,NBH,NEELEM,NELIST,NFF,NHE(1,nx),
     '        NHP(1,0,nx),nh_update,NKH,NKHE,NKJE,NLL,NP_average,NPF,
     '        NPLIST,NPNE,NPNODE,NRE,NVHE,NVHP,NVJE,NW(1,1,nx),nx,
     '        NYNE,NYNP,VOLTC,ELEM_VOL,nj_volume,
     '        CURVCORRECT,PG,RG,SE,VOL,VOLT,WG,XA,XAB,XE,XG,
     '        XP,YP(1,1,nx),ZA,ZE,ZG,ZP,
     '        AVGENODE,OUTPUT_TOT_ONLY,TYPE,ERROR,*9999)
C          ENDDO !nr
C old MPN 8Apr97: regions handled incorrectly
CC GMH 3/9/95 Must initialise VOLT before calls
C          DO ni=1,NIM
C            DO nb=1,NBFM
C              VOLT(ni,nb)=0.0D0
C            ENDDO !nb
C          ENDDO !ni
C          DO nolist=1,NRLIST(0)
C            nr=NRLIST(nolist)
C            IF(TYPE(1:8).EQ.'DEFORMED')THEN
C              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
C     '          NKH(1,1,1,nr),
C     '          NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
C     '          ERROR,*9999)
C              IF(AVGENODE) THEN
CC               Calc the average nj_update value for nodes in NPLIST
C                SUM=0.0d0
C                DO nonode=1,NPLIST(0)
C                  np=NPLIST(nonode)
C                  SUM=SUM+ZP(1,1,nj_update,np,1)
C                ENDDO !np
C                ZP(1,1,nj_update,NP_average,1)=SUM/DBLE(NPLIST(0))
C              ENDIF
C            ENDIF
C            CALL OPELEM(IBT,NBJ,NBH,NEELEM,NELIST,NFF,NHE(1,nx),
C     '        NKE,NLL,NPF,NPNE,nr,NVHE,NVJE,NW(1,1,nx),nx,
C     '        PG,RG,SE,VOL,VOLT,WG,XA,XE,XG,XP,ZA,ZE,ZG,ZP,
C     '        OUTPUT_TOT_ONLY,TYPE,ERROR,*9999)
C          ENDDO
C end old
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIELEM')
      RETURN
 9999 CALL ERRORS('LIELEM',ERROR)
      CALL EXITS('LIELEM')
      RETURN 1
      END


