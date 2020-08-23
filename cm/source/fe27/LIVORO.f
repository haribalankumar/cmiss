      SUBROUTINE LIVORO(NFVC,NODENVC,NODENVCB,NPNODE,NVCB,NVCNODE,
     '  VC,XNFV,XP,STRING,ERROR,*)

C#### Subroutine: LIVORO
C###  Description:
C###    LIVORO outputs the Voronoi mesh

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),NODENVCB(NVCBM),
     '  NPNODE(0:NP_R_M,0:NRM),NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M)
      REAL*8 VC(0:NVCM),XNFV(-(NJM+1):NJM,NFVM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,nr,N3CO,nvc,nfvl,np,nj,nfv,
     '  nonode,cnonode,MAX_LOCFACE,bnvc,bnp,bnonode,i,cnvc
      REAL*8 TOTAREA,TOTDIST,AVEAREA,AVEDIST,SURFAREA,VOL_DIFF_TOT,
     '  VOL_DIFF_AVG,NVC_VOL,VOL_DIFF_MAX,VOL_DIFF
      CHARACTER FILE*100
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LIVORO',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list voronoi;filename
C###  Description:
C###    Lists Voronoi cell attributes
C###  Parameter:      <total> <region #[1]>
C###    total lists the total Voronoi statistics
C###    region lists Voronoi region #

        OP_STRING(1)=BLANK(1:15)//'total <region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        OP_STRING(2)=BLANK(1:15)//'wall/inlet/outlet/free/driving'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opvoro','NEW',
     '      'SEQUEN','FORMATTED',160,ERROR,*9999)
        ENDIF
        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
          WRITE(*,*) 'Defaulting to region 1.  Is this the correct' //
     '      ' region?'
        ENDIF
        IF(CBBREV(CO,'WALL',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_INIT,'Need to define BCs.',ERROR,*9999)
          WRITE(OP_STRING,'('' Listing Wall Boundary Connecs'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO bnvc=1,NVBT
            bnonode=NODENVCB(bnvc)
            bnp=NPNODE(bnonode,nr)
            IF(NVCB(BCTYPE,bnvc).EQ.WALL) THEN
              WRITE(OP_STRING,'(/,'' Boundary cell: '',I6,'' bnp '//
     '          '= '',I6,'' nonode = '',I6)') bnvc,bnp,bnonode
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Connects to:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO i=1,NVCB(0,bnvc)
                nonode=NVCB(i,bnvc)
                np=NPNODE(nonode,nr)
                nvc=NVCNODE(MAP,nonode)
C             ..Find the connection
                DO nfvl=1,NFVC(1,0,nvc)
                  IF(NFVC(1,nfvl,nvc).EQ.bnonode) THEN
                    nfv=NFVC(2,nfvl,nvc)
                    GOTO 10
                  ENDIF
                ENDDO
 10             CONTINUE
                WRITE(OP_STRING,'('' nvc = '',I6,'' np '//
     '            '= '',I6,'' Area = '',D12.4,'' D'//
     '            'ist = '',D12.4,'' Normal = '',3D12.4)')
     '            nvc,np,XNFV(FAREA,nfv),1.0D0/XNFV(IDIST,nfv),
     '            (XNFV(nj,nfv),nj=1,NJ_LOC(NJL_GEOM,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ELSEIF(CBBREV(CO,'INLET',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_INIT,'Need to define BCs.',ERROR,*9999)
          WRITE(OP_STRING,'('' Listing Inlet Boundary Connecs'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO bnvc=1,NVBT
            bnonode=NODENVCB(bnvc)
            bnp=NPNODE(bnonode,nr)
            IF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN
              WRITE(OP_STRING,'(/,'' Boundary cell: '',I6,'' bnp '//
     '          '= '',I6,'' nonode = '',I6)') bnvc,bnp,bnonode
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Connects to:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO i=1,NVCB(0,bnvc)
                nonode=NVCB(i,bnvc)
                np=NPNODE(nonode,nr)
                nvc=NVCNODE(MAP,nonode)
C             ..Find the connection
                DO nfvl=1,NFVC(1,0,nvc)
                  IF(NFVC(1,nfvl,nvc).EQ.bnonode) THEN
                    nfv=NFVC(2,nfvl,nvc)
                    GOTO 20
                  ENDIF
                ENDDO
 20             CONTINUE
                WRITE(OP_STRING,'('' nvc = '',I6,'' np '//
     '            '= '',I6,'' Area = '',D12.4,'' D'//
     '            'ist = '',D12.4,'' Normal = '',3D12.4)')
     '            nvc,np,XNFV(FAREA,nfv),1.0D0/XNFV(IDIST,nfv),
     '            (XNFV(nj,nfv),nj=1,NJ_LOC(NJL_GEOM,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO

        ELSEIF(CBBREV(CO,'OUTLET',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_INIT,'Need to define BCs.',ERROR,*9999)
          WRITE(OP_STRING,'('' Listing Outlet Boundary Connecs'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO bnvc=1,NVBT
            bnonode=NODENVCB(bnvc)
            bnp=NPNODE(bnonode,nr)
            IF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN
              WRITE(OP_STRING,'(/,'' Boundary cell: '',I6,'' bnp '//
     '          '= '',I6,'' nonode = '',I6)') bnvc,bnp,bnonode
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Connects to:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO i=1,NVCB(0,bnvc)
                nonode=NVCB(i,bnvc)
                np=NPNODE(nonode,nr)
                nvc=NVCNODE(MAP,nonode)
C             ..Find the connection
                DO nfvl=1,NFVC(1,0,nvc)
                  IF(NFVC(1,nfvl,nvc).EQ.bnonode) THEN
                    nfv=NFVC(2,nfvl,nvc)
                    GOTO 30
                  ENDIF
                ENDDO
 30             CONTINUE
                WRITE(OP_STRING,'('' nvc = '',I6,'' np '//
     '            '= '',I6,'' Area = '',D12.4,'' D'//
     '            'ist = '',D12.4,'' Normal = '',3D12.4)')
     '            nvc,np,XNFV(FAREA,nfv),1.0D0/XNFV(IDIST,nfv),
     '            (XNFV(nj,nfv),nj=1,NJ_LOC(NJL_GEOM,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ELSEIF(CBBREV(CO,'FREE',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_INIT,'Need to define BCs.',ERROR,*9999)
          WRITE(OP_STRING,'('' Listing Freeslip Boundary Connecs'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO bnvc=1,NVBT
            bnonode=NODENVCB(bnvc)
            bnp=NPNODE(bnonode,nr)
            IF(NVCB(BCTYPE,bnvc).EQ.FREESLIP) THEN
              WRITE(OP_STRING,'(/,'' Boundary cell: '',I6,'' bnp '//
     '          '= '',I6,'' nonode = '',I6)') bnvc,bnp,bnonode
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Connects to:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO i=1,NVCB(0,bnvc)
                nonode=NVCB(i,bnvc)
                np=NPNODE(nonode,nr)
                nvc=NVCNODE(MAP,nonode)
C             ..Find the connection
                DO nfvl=1,NFVC(1,0,nvc)
                  IF(NFVC(1,nfvl,nvc).EQ.bnonode) THEN
                    nfv=NFVC(2,nfvl,nvc)
                    GOTO 40
                  ENDIF
                ENDDO
 40             CONTINUE
                WRITE(OP_STRING,'('' nvc = '',I6,'' np '//
     '            '= '',I6,'' Area = '',D12.4,'' D'//
     '            'ist = '',D12.4,'' Normal = '',3D12.4)')
     '            nvc,np,XNFV(FAREA,nfv),1.0D0/XNFV(IDIST,nfv),
     '            (XNFV(nj,nfv),nj=1,NJ_LOC(NJL_GEOM,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ELSEIF(CBBREV(CO,'DRIVING',3,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALL_INIT,'Need to define BCs.',ERROR,*9999)
          WRITE(OP_STRING,'('' Listing Driving Boundary Connecs'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO bnvc=1,NVBT
            bnonode=NODENVCB(bnvc)
            bnp=NPNODE(bnonode,nr)
            IF(NVCB(BCTYPE,bnvc).EQ.DRIVING) THEN
              WRITE(OP_STRING,'(/,'' Boundary cell: '',I6,'' bnp '//
     '          '= '',I6,'' nonode = '',I6)') bnvc,bnp,bnonode
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Connects to:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO i=1,NVCB(0,bnvc)
                nonode=NVCB(i,bnvc)
                np=NPNODE(nonode,nr)
                nvc=NVCNODE(MAP,nonode)
C             ..Find the connection
                DO nfvl=1,NFVC(1,0,nvc)
                  IF(NFVC(1,nfvl,nvc).EQ.bnonode) THEN
                    nfv=NFVC(2,nfvl,nvc)
                    GOTO 50
                  ENDIF
                ENDDO
 50             CONTINUE
                WRITE(OP_STRING,'('' nvc = '',I6,'' np '//
     '            '= '',I6,'' Area = '',D12.4,'' D'//
     '            'ist = '',D12.4,'' Normal = '',3D12.4)')
     '            nvc,np,XNFV(FAREA,nfv),1.0D0/XNFV(IDIST,nfv),
     '            (XNFV(nj,nfv),nj=1,NJ_LOC(NJL_GEOM,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO
        ELSEIF(CBBREV(CO,'TOTAL',1,noco+1,NTCO,N3CO)) THEN
          WRITE(OP_STRING,'('' Listing Voronoi cell total:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/,'' Total no. of cells..... '',I14)') NVCT
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Total Volume........... '',D14.4)')
     '      VC(0)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average cell volume.... '',D14.4)')
     '      VC(0)/DBLE(NVCT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/,'' Connectivity Statistics:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/,'' Total no. of connections.... '',I14)')
     '      NFVT
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average no. of connections.. '',
     '      D14.4)') DBLE(NFVT)/DBLE(NVCT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          TOTAREA=0.d0
          TOTDIST=0.d0
          DO nfv=1,NFVT
            TOTAREA=TOTAREA+XNFV(FAREA,nfv)
            TOTDIST=TOTDIST+1.d0/XNFV(IDIST,nfv)
          ENDDO

C         ..Adding mesh quality quantification
          VOL_DIFF_AVG=0.d0
          VOL_DIFF_MAX=0.d0
          DO nvc=1,NVCT
            NVC_VOL=VC(nvc)
            VOL_DIFF_TOT=0.d0
            DO nfvl=1,NFVC(1,0,nvc)
              cnonode=NFVC(1,nfvl,nvc)
              IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
                cnvc=NVCNODE(MAP,cnonode)
                VOL_DIFF=DABS(NVC_VOL-VC(cnvc))
                VOL_DIFF_TOT=VOL_DIFF_TOT+VOL_DIFF
                VOL_DIFF=VOL_DIFF/NVC_VOL
                IF(VOL_DIFF.GT.VOL_DIFF_MAX) VOL_DIFF_MAX=VOL_DIFF
              ENDIF
            ENDDO
            VOL_DIFF_TOT=VOL_DIFF_TOT/NVC_VOL
            VOL_DIFF_AVG=VOL_DIFF_AVG+VOL_DIFF_TOT/DBLE(NFVC(1,0,nvc))
          ENDDO
          VOL_DIFF_AVG=VOL_DIFF_AVG/(2.d0*DBLE(NVCT))


C         .. Adding external surface area
          SURFAREA=0.d0
          DO nvc=1,NVCT
            DO nfvl=1,NFVC(1,0,nvc)
              cnonode=NFVC(1,nfvl,nvc)
              nfv=NFVC(2,nfvl,nvc)
              IF(NVCNODE(TYPE,cnonode).EQ.BOUNDARY) THEN
                SURFAREA=SURFAREA+XNFV(FAREA,nfv)
              ENDIF
            ENDDO
          ENDDO
          TOTAREA=TOTAREA/2.d0
          TOTDIST=TOTDIST/2.d0
          AVEAREA=TOTAREA/DBLE(NFVT)
          AVEDIST=TOTDIST/DBLE(NFVT)

C         .. Max number of local connections ..
          MAX_LOCFACE=0
          DO nvc=1,NVCT
            IF(NFVC(1,0,nvc).GT.MAX_LOCFACE) MAX_LOCFACE=NFVC(1,0,nvc)
          ENDDO

          WRITE(OP_STRING,'('' Max no. of local faces...... '',I14)')
     '      MAX_LOCFACE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Total area of all cells..... '',D14.4)')
     '      TOTAREA
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Total external area......... '',D14.4)')
     '      SURFAREA
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average area of all cells... '',D14.4)')
     '      AVEAREA
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average dist between cells.. '',D14.4,
     '      /$)')
     '      AVEDIST
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average volume difference... '',D14.4)')
     '      VOL_DIFF_AVG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Maximum volume difference... '',D14.4)')
     '      VOL_DIFF_MAX
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE
          CALL OPVORO(NFVC,NODENVC,NPNODE,nr,VC,XNFV,XP,ERROR,*9999)
c          WRITE(OP_STRING,'('' Listing Voronoi cells:'')')
c          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c          DO nvc=1,NVCT
c            nonode=NODENVC(nvc)
c            np=NPNODE(nonode,nr)
c            CALL ASSERT( np.NE.0, '>> Global node not found. Check ' //
c     '        ' that the correct region is specified',ERROR,*9999)
c            WRITE(OP_STRING,'(/$,''Voronoi cell:'',I6)') nvc
c            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            WRITE(OP_STRING,'(''Node:        '',
c     '        I6,'' (nonode = '',I6,'')'')') np,nonode
c            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            WRITE(OP_STRING,'(''Coordinates:    '',3D12.4)')
c     '        (XP(1,1,nj,np),nj=1,NJ_LOC(NJL_GEOM,0,nr))
c            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            WRITE(OP_STRING,'(/$,''Adjacency Molecule:'')')
c            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            WRITE(OP_STRING,'(''==================='')')
c            CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C TVK 06/04/99 Changed output to accomdate 2D and 3D

c            IF(NJT.EQ.2) THEN
c              WRITE(OP_STRING,'(''  Node  Face        Area'//
c     '        '        Dist  Centroid         '//
c     '        '       Norm'')')
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              DO nfvl=1,NFVC(1,0,nvc)
c                WRITE(OP_STRING,
c     '            '(I6,I6,D12.4,D12.4,D12.4,D12.4,2D12.4)')
c     '            NFVC(1,nfvl,nvc),NFVC(2,nfvl,nvc),
c     '            XNFV(FAREA,NFVC(2,nfvl,nvc)),
c     '            1.0D0/XNFV(IDIST,NFVC(2,nfvl,nvc)),
c     '            XNFV(-2,NFVC(2,nfvl,nvc)),
c     '            XNFV(-3,NFVC(2,nfvl,nvc)),
c     '            (XNFV(nj,NFVC(2,nfvl,nvc)),nj=1,NJ_LOC(NJL_GEOM,0,nr))
c                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              ENDDO
c            ELSE
c              WRITE(OP_STRING,'(''  Node  Face        Area'//
c     '        '        Dist  Centroid                     '//
c     '        '       Norm'')')
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              DO nfvl=1,NFVC(1,0,nvc)
c                WRITE(OP_STRING,
c     '            '(I6,I6,D12.4,D12.4,D12.4,D12.4,D12.4,3D12.4)')
c     '            NFVC(1,nfvl,nvc),NFVC(2,nfvl,nvc),
c     '            XNFV(FAREA,NFVC(2,nfvl,nvc)),
c     '            1.0D0/XNFV(IDIST,NFVC(2,nfvl,nvc)),
c     '            XNFV(-2,NFVC(2,nfvl,nvc)),
c     '            XNFV(-3,NFVC(2,nfvl,nvc)),
c     '            XNFV(-4,NFVC(2,nfvl,nvc)),
c     '            (XNFV(nj,NFVC(2,nfvl,nvc)),nj=1,NJ_LOC(NJL_GEOM,0,nr))
c                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              ENDDO
c            ENDIF !NJT
c            WRITE(OP_STRING,'(/$,''Volume:'',D12.4)') VC(nvc)
c            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c          ENDDO
c          WRITE(OP_STRING,'(/$,''Total Volume:'',D12.4)') VC(0)
c          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(OPFILE) THEN
            CALL CLOSEF(IOFI,ERROR,*9999)
            IOFI=IOOP
          ENDIF
        ENDIF
      ENDIF
      CALL EXITS('LIVORO')
      RETURN
 9999 CALL ERRORS('LIVORO',ERROR)
      CALL EXITS('LIVORO')
      RETURN 1
      END


