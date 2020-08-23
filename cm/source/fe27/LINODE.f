      SUBROUTINE LINODE(NBH,NEELEM,NHE,NHP,NKH,NKJ,NP_INTERFACE,NPNODE,
     '  NPLIST,NRLIST,NVHP,NVJP,NWP,NXLIST,NYNE,NYNP,XP,YP,YP1,ZA,ZP,
     '  STRING,ERROR,*)

C#### Subroutine: LINODE
C###  Description:
C###    LINODE lists nodal coordinates and deformed position or
C###    displacements or reactions.
C###    AJP 29-3-94.  Also lists gradient of solution at nodes.
C###    JWF 02-3-02   Also lists contact forces at nodes.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NP_INTERFACE(0:NPM,0:3),NPNODE(0:NP_R_M,0:NRM),NPLIST(0:NPM),
     '  NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NWP(NPM,2),NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IY,N3CO,nc,nj,nk,
     '  nolist,nonode,np,nr,nv,nx,nxc,ny,ny1,niy,nxx
      REAL*8 GRADPHI(3),X(3),Z(3),YP1(NYM,NIYM,NXM)
      CHARACTER FILE*100,TYPE*14
      LOGICAL ALL_REGIONS,CBBREV,EXTEND,INTERFACE,OPFILE,PNODES

      CALL ENTERS('LINODE',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list nodes<;FILENAME>
C###  Description:
C###    List the coordinate and derivative values of each node to the
C###    screen or to file FILENAME.opnode if qualifier is present.
C###  Parameter:      <(solution/displacement/reaction/flux/cartesian
C###  /gradient/contact forces)>
C###    Specify another piece of nodal information to list in addition
C###    to geometry.
C###  Parameter:      <nodes (all/GROUP/#s)[all]>
C###    Specify the nodes to list. The 'all' command prompts all
C###    currently defined nodes to be listed.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to list. The "all" keyword indicates
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    For specify the class number (of solve type) for the additional
C###    nodal information.

      OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
      OP_STRING(2)=BLANK(1:15)
     '    //'<(solution/displacement/reaction/flux/cartesian/gradient/co
     'ntact_forces)>'
      OP_STRING(3)=BLANK(1:15)//'<nodes (all/GROUP/#s)[all]>'
      OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
      OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
      CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list nodes<;FILENAME> interface
C###  Description:
C###    List the regions that in which each nodes is defined to the
C###    screen or to file FILENAME.opnode if qualifier is present.
C###  Parameter:      <nodes (all/GROUP/#s)[all]>
C###    Specify the nodes to list. The 'all' command prompts all
C###    currently defined nodes to be listed.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to list. The "all" keyword indicates
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> interface'
        OP_STRING(2)=BLANK(1:15)//'<nodes (all/GROUP/#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (all/#s)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list nodes<;FILENAME> groups
C###  Description:
C###    List node groups to the screen or to FILENAME.opnode if
C###    qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> groups'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LINODE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opnode','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          PNODES=.TRUE.
        ELSE !list all nodes
          PNODES=.FALSE.
C MLB 28/5/97 moved inside main region loop
C          NPLIST(0)=0
C          DO nolist=1,NRLIST(0)
C            nr=NRLIST(nolist)
C            DO nonode=1,NPNODE(0,nr)
C              NPLIST(NPLIST(0)+nonode)=NPNODE(nonode,nr)
C            ENDDO !nonode
C            NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
C          ENDDO !nr
        ENDIF

C!!!    Temporary until people get used to command
        IF(CBBREV(CO,'NUMBER',1,noco+1,noco+3,N3CO)) THEN
          WRITE(OP_STRING,'(''>>Repeat command replacing the '
     '      //' option NUMBER with NODE (and the #(s)'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GO TO 9999
        ENDIF

        nc=1   !defaults
        iy=1
        IF(CBBREV(CO,'CARTESIAN',2,noco+1,NTCO,N3CO)) THEN
          TYPE='Cartesian'
        ELSE IF(CBBREV(CO,'DISPLACEMENT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='Displacement'
        ELSE IF(CBBREV(CO,'FLUX',1,noco+1,NTCO,N3CO)) THEN
          TYPE='Flux'
          nc=2
        ELSE IF(CBBREV(CO,'GRADIENT',3,noco+1,NTCO,N3CO)) THEN
          TYPE='Gradient'
        ELSE IF(CBBREV(CO,'GROUPS',3,noco+1,NTCO,N3CO)) THEN
          TYPE='GROUPS'
        ELSE IF(CBBREV(CO,'REACTION',1,noco+1,NTCO,N3CO)) THEN
          TYPE='Reaction'
        ELSE IF(CBBREV(CO,'SOLUTION',1,noco+1,NTCO,N3CO)) THEN
          TYPE='Solution'
        ELSE IF(CBBREV(CO,'CONTACT_FORCES',2,noco+1,NTCO,N3CO)) THEN
          TYPE='Contact_forces'
        ELSE
          TYPE='Geometry'
        ENDIF

        IF(CBBREV(CO,'INTERFACE',2,noco+1,NTCO,N3CO)) THEN
          INTERFACE=.TRUE.
        ELSE
          INTERFACE=.FALSE.
        ENDIF

        IF(TYPE(1:6).EQ.'GROUPS') THEN
          CALL OPNODEG(ERROR,*9999)

        ELSE !write out nodal variables

          IF(TYPE(1:8).NE.'Geometry') THEN
            CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
            nxc=NXLIST(1)
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9999)
            IF(TYPE(1:8).EQ.'Reaction') THEN
              IF(ITYP6(NRLIST(1),nx).EQ.2) THEN !nonlinear problems
                nc=1
                IY=4
              ELSE
                nc=2
                IY=1
              ENDIF
            ENDIF
            IF(TYPE(1:14).EQ.'Contact_forces') THEN
              IF(ITYP6(NRLIST(1),nx).EQ.2) THEN !nonlinear problems
                nc=1
                IY=6
              ELSE
                nc=2
                IY=1
              ENDIF
            ENDIF
          ENDIF
C*** JHC 20/2/08 For reaction or contact forces, ny for YP(ny,4 or 6) is not the same as ny for YP(ny,1) for nonlinear problem
C        especially for multi-region problem. YP1 is temporary array, copy of YP but using appropriate ny index 
          IF((TYPE(1:8).EQ.'Reaction'.OR.
     &      TYPE(1:14).EQ.'Contact_forces').AND.nc.EQ.1) THEN ! contact force or reaction
            DO ny=1,NYM
              DO niy=1,NIYM
                DO nxx=1,NXM
                  YP1(ny,niy,nxx)=0.0d0
                ENDDO
              ENDDO
            ENDDO
            DO nolist=1,NRLIST(0)
              nr=NRLIST(nolist)
              IF(.NOT.PNODES) THEN
                NPLIST(0)=0
                DO nonode=1,NPNODE(0,nr)
                  NPLIST(NPLIST(0)+nonode)=NPNODE(nonode,nr)
                ENDDO !nonode
                NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
              ENDIF
              DO nonode=1,NPLIST(0)
                np=NPLIST(nonode)
                DO nj=1,NJT
                  DO nv=1,NVJP(nj,np)
                    DO nk=1,NKJ(nj,np)
                      ny=NYNP(nk,nv,nj,np,1,nc,nr)
                      ny1=NYNP(nk,nv,nj,np,0,nc,nr)
                      DO nxx=1,NXM
                        YP1(ny1,IY,nxx)=YP(ny,IY,nxx)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C*** END JHC
          DO nolist=1,NRLIST(0)
            nr=NRLIST(nolist)
            IF(.NOT.PNODES) THEN
              NPLIST(0)=0
              DO nonode=1,NPNODE(0,nr)
                NPLIST(NPLIST(0)+nonode)=NPNODE(nonode,nr)
              ENDDO !nonode
              NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
            ENDIF
            IF(TYPE(1:8).NE.'Geometry') THEN
              EXTEND=.TRUE.
              IF(TYPE(1:8).EQ.'Gradient') THEN
                DO nonode=1,NPLIST(0)
                  np=NPLIST(nonode)
C                 CALL GRADPHI_N(NKH(1,1,1,nr),NP_INTERFACE,
C    '              np,nr,NYNP,
C    '              GRADPHI,XG,XN,XP,YP(1,1,nx),ERROR,*9999)
                  DO nj=1,NJT
                    DO nv=1,NVJP(nj,np)
                      ZP(1,nv,nj,np,nc)=GRADPHI(nj)
                    ENDDO !nv
                  ENDDO !nj
                ENDDO !nonode (np)
              ELSE
C*** JHC 20/2/08 Rather than passing YP, use YP1 for reaction or contact forces for nonlinear problem
                IF((TYPE(1:8).EQ.'Reaction'.OR.
     &            TYPE(1:14).EQ.'Contact_forces').AND.nc.EQ.1) THEN ! contact force or reaction
                  CALL YPZP(iy,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '              NKH(1,1,1,nr),NPNODE,
     '              nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP1(1,1,nx),ZA,ZP,
     '              ERROR,*9999)
                ELSE
                  CALL YPZP(iy,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '              NKH(1,1,1,nr),NPNODE,
     '              nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '              ERROR,*9999)
                ENDIF
C*** END JHC
C CPB 23/3/96 Only subtract XP if solution is not already a displ.
               IF(TYPE(1:12).EQ.'Displacement'.AND.
     '           ITYP6(nr,nx).EQ.2.AND.KTYP58(nr).EQ.1) THEN
                  DO nonode=1,NPLIST(0)
                    np=NPLIST(nonode)
                    DO nj=1,NJT
                      DO nv=1,NVJP(nj,np)
                        DO nk=1,NKJ(nj,np)
                          ZP(nk,nv,nj,np,nc)=ZP(nk,nv,nj,np,nc)-
     '                      XP(nk,nv,nj,np)
                        ENDDO !nk
                      ENDDO !nv
                    ENDDO !nj
                  ENDDO !nonode (np)
                ENDIF
              ENDIF !type
            ELSE
              EXTEND=.FALSE.
            ENDIF !type

            IF(TYPE(1:9).EQ.'Cartesian') THEN
              EXTEND=.TRUE.
              DO nonode=1,NPLIST(0)
                np=NPLIST(nonode)
                DO nj=1,NJT
                  X(nj)=XP(1,1,nj,np)
                ENDDO
                CALL XZ(ITYP10(nr),X,Z)
                DO nj=1,NJT
                  ZP(1,1,nj,np,nc)=Z(nj)
                  DO nk=2,NKJ(nj,np)
                    ZP(nk,1,nj,np,nc)=0.0d0
                  ENDDO
                ENDDO
              ENDDO
C old MPN 6Mar96
C              IF(CBBREV(CO,'NUMBER',1,noco+1,noco+3,N3CO)) THEN
C                np=IFROMC(CO(N3CO+1))
C                CALL OPNODE1(NKJ(1,np),NKJ(1,np),np,
C     '            NP_INTERFACE,nr,NVJP,
C     '            XP(1,1,1,np),ZP(1,1,1,np,nc),
C     '            EXTEND,INTERFACE,ERROR,*9999)
C              ELSE
                CALL OPNODE(nc,NHP,
     '            NKH(1,1,1,nr),NKJ,NP_INTERFACE,
     '            NPLIST,NPNODE,nr,
     '            NVJP,NWP,nx,XP,ZP,
     '            EXTEND,INTERFACE,TYPE,ERROR,*9999)
C              ENDIF

            ELSE
              !write variables to screen or file
C old MPN 6Mar96
C              IF(CBBREV(CO,'NUMBER',1,noco+1,noco+3,N3CO)) THEN
C          np=IFROMC(CO(N3CO+1))
C          CALL OPNODE1(NKH(1,np,nc,nr),NKJ(1,np),np,
C     '            NP_INTERFACE,nr,NVJP,
C     '      XP(1,1,1,np),ZP(1,1,1,np,nc),
C     '      EXTEND,INTERFACE,ERROR,*9999)
C              ELSE
              CALL OPNODE(nc,NHP,NKH(1,1,1,nr),NKJ,NP_INTERFACE,
     '          NPLIST,NPNODE,nr,
     '          NVJP,NWP,nx,XP,ZP,
     '          EXTEND,INTERFACE,TYPE,ERROR,*9999)
C              ENDIF
            ENDIF !type
          ENDDO !nr
        ENDIF !groups/nodes

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LINODE')
      RETURN
 9999 CALL ERRORS('LINODE',ERROR)
      CALL EXITS('LINODE')
      RETURN 1
      END


