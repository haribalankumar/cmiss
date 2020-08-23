      SUBROUTINE IPINI6(NHP,NKH,NODENVCB,NPNODE,nr,
     '  NVCB,NVCNODE,NVHP,nx,NYNP,YP,ALL_REGIONS,ERROR,*)

C#### Subroutine: IPINI6
C###  Description:
C###    IPINI6 inputs initial conditions and boundary conditions for
C###    FE60 problems.

C**** YP(ny,1) contains essential b.c.s
C**** YP(ny,3)    "     initial solution.
C**** YP(ny,8) contains solution at the "old" time step
C****

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NHP(NPM),NKH(NHM,NPM,NCM),NODENVCB(NVCBM), !NODENVC(NVCM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),
     '  NVHP(NHM,NPM,NCM),nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL ALL_REGIONS
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES,nonode,np,bnvc,nk,nv,nh,nhx,
     '  nonode2,np2,n,np_to_find,ny,bnonode,bnp,ib,bny,i
      CHARACTER CHARnonode*5,CHARnp*5,CHARbnvc*5,
     '  CHAR1*1,CHAR2*12
      LOGICAL FILEIP

      CALL ENTERS('IPINI6',*9999)

C     ..Early exit if possible
      CALL ASSERT(USE_VORONOI.EQ.1,'>>Set USE_VORONOI = 1',ERROR,*9999)

C     ..Initialise Variables
      ICHAR=999
      NOQUES=0
      FILEIP=.FALSE.

C     Need these asserts to ensure we can store old solution
      CALL ASSERT(NIYM.GE.12,'Increase NIYM >= 12',ERROR,*9999)

      IF(IOTYPE.ne.3) THEN
C***    Need to only initialise those belonging to the
C***    current region.
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhx=1,NHP(np)
            nh=nh_loc(nhx,nx)
            DO nv=1,NVHP(nh,np,1)
              DO nk=1,MAX(NKH(nh,np,1)-KTYP93(1,nr),1)
                ny=NYNP(nk,nv,nh,np,0,1,nr)
                YP(ny,1)=0.d0
                YP(ny,3)=0.d0
                YP(ny,8)=0.d0
                YP(ny,12)=0.d0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF(ITYP1(nr,nx).EQ.6) THEN !Fluid dynamics

        DO bnvc=1,NVBT

          nonode=NODENVCB(bnvc)
          np=NPNODE(nonode,nr)
          WRITE(CHARnonode,'(I5)') nonode
          WRITE(CHARnp,'(I5)') np
          WRITE(CHARbnvc,'(I5)') bnvc
          FORMAT='(/$,'' Enter the node numbers (IB nodes)'//
     '      ' adjoining B'//
     '      ' node '//CHARbnvc//' (np = '//CHARnp//'): '',I6)'

          IF(IOTYPE.EQ.3) THEN
            IDATA(0)=NVCB(0,bnvc)
            DO n=1,NVCB(0,bnvc)
              IDATA(n)=NVCB(n,bnvc)
            ENDDO
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '      FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '      IDATA,IONE,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '      RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CALL ASSERT(IDATA(0).NE.0,'>>A B node must map onto an IB'//
     '        ' node',ERROR,*9999)
            NVCB(0,bnvc)=IDATA(0)
            IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
              CALL ASSERT(NVCB(0,bnvc).LE.2,'>>Cant match a B'//
     '          ' node onto more than two IB nodes',ERROR,*9999)
            ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
              CALL ASSERT(NVCB(0,bnvc).LE.2,'>>Cant match a B'//
     '          ' node onto more than three IB nodes',ERROR,*9999)
            ENDIF
            DO n=1,IDATA(0)
              np_to_find=IDATA(n)
              DO nonode2=1,NPNODE(0,nr)
                np2=NPNODE(nonode2,nr)
                IF(np2.EQ.np_to_find) THEN
                  CALL ASSERT(NVCNODE(TYPE,nonode2).EQ.INTBOUN,
     '              '>>A B node must map onto an IB node',ERROR,*9999)
                  NVCB(n,bnvc)=nonode2
                  GOTO 102
                ENDIF
              ENDDO !nonode
 102          CONTINUE
            ENDDO !n
          ENDIF !iotype.NE.3

          IDEFLT(1)=0
          FORMAT='('' Specify the boundary condition type '
     '      //'[0]: '''//
     '      '/''   (0) Wall'''//
     '      '/''   (1) Inlet'''//
     '      '/''   (2) Outlet'''//
     '      '/''   (3) Free'''//
     '      '/''   (4) Driving'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NVCB(BCTYPE,bnvc)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NVCB(BCTYPE,bnvc)=IDATA(1)

          IF(NVCB(BCTYPE,bnvc).EQ.WALL.OR.NVCB(BCTYPE,bnvc).EQ.INLET)
     '      THEN
            DO nhx=1,nh_loc(0,nx)-1
              ny=NYNP(1,1,nh_loc(nhx,nx),np,0,1,nr)
              RDEFLT(1)=0.d0
              WRITE(CHAR1,'(I1)') nhx
              WRITE(CHAR2,'(E12.5)') RDEFLT(1)
              FORMAT='($,'' Velocity component '//CHAR1(1:1)//' is'//
     '          ' ['//CHAR2(1:12)//']: '',G25.17)'
              IF(IOTYPE.EQ.3) RDATA(1)=YP(ny,3)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '          FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '          IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '          RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                YP(ny,1)=RDATA(1)
                YP(ny,3)=RDATA(1)
                YP(ny,8)=RDATA(1)
                YP(ny,12)=RDATA(1)
              ENDIF
            ENDDO
          ENDIF
        ENDDO !bnvc
      ENDIF !ityp1(nr,nx).eq.6

!     Adding fixed pressure node - rgb 190598
      FORMAT='(/$,'' Is the solution domain confined ? [N] '',A)'
      ADEFLT(1)='N'
      IF(IOTYPE.EQ.3) THEN
        IF(CONFINED) THEN
          ADATA(1)='Y'
        ELSE
          ADATA(1)='N'
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF((ADATA(1).EQ.'Y').OR.(ADATA(1).EQ.'y')) THEN
          CONFINED=.TRUE.
        ELSE
          CONFINED=.FALSE.
        ENDIF
      ENDIF
      IF(CONFINED) THEN
        IDEFLT(1)=1
        FORMAT='($,'' Enter the node number that you wish to'
     '    //' have the pressure fixed at [1]: '',I6)'
        IF(IOTYPE.EQ.3) IDATA(1)=FIXDNODE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     '    NPNODE(0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '    ERROR,*9999)
        IF(IOTYPE.NE.3) FIXDNODE=IDATA(1)
      ENDIF

!     Adding duct flow outlet boundary condition
      IF(.NOT.CONFINED) THEN
        FORMAT='(/$,'' Do you wish to use duct flow outlets ? [N] '',A)'
        ADEFLT(1)='N'
        IF(IOTYPE.EQ.3) THEN
          IF(DUCTFLOW) THEN
            ADATA(1)='Y'
            ADEFLT(1)='Y'
          ELSE
            ADATA(1)='N'
            ADEFLT(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF((ADATA(1).EQ.'Y').OR.(ADATA(1).EQ.'y')) THEN
            DUCTFLOW=.TRUE.
          ELSE
            DUCTFLOW=.FALSE.
          ENDIF
        ENDIF
        IF(DUCTFLOW) THEN
          NUM_OUTLETS=1
          IDEFLT(1)=NUM_OUTLETS
          IF(IOTYPE.EQ.3) IDATA(1)=NUM_OUTLETS
          FORMAT='($,'' The number of outlets is [1]: '',I1)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '      FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '      IDATA,IDEFLT,1,10,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '      RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NUM_OUTLETS=IDATA(1)
          CALL ASSERT(NUM_OUTLETS.LE.9,'>>Can only handle up to 9'//
     '      ' outlets at present',ERROR,*9999)
          DO i=1,NUM_OUTLETS
            WRITE(CHAR1,'(I1)') i
            IDEFLT(1)=OUTLET_FIXDNODES(i)
            IF(IOTYPE.EQ.3) IDATA(1)=OUTLET_FIXDNODES(i)
            FORMAT='($,'' Enter the fixed pressure node '//
     '        'for outlet '//CHAR1(1:1)//': '',I5)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '        FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '        IDATA,IONE,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '        RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) OUTLET_FIXDNODES(i)=IDATA(1)
          ENDDO !i
        ENDIF
      ENDIF

      DO bnvc=1,NVBT
        IF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN
          bnonode=NODENVCB(bnvc)
          bnp=NPNODE(bnonode,nr)
          DO ib=1,NVCB(0,bnvc)
            nonode=NVCB(ib,bnvc)
            np=NPNODE(nonode,nr)
            DO nhx=1,nh_loc(0,nx)-1
              bny=NYNP(1,1,nh_loc(nhx,nx),bnp,0,1,nr)
              ny=NYNP(1,1,nh_loc(nhx,nx),np,0,1,nr)
              YP(ny,1)=YP(bny,3)
              YP(ny,8)=YP(bny,3)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      IF(FILEIP.AND..NOT.ALL_REGIONS) CALL CLOSEF(IFILE,ERROR,*9999)

      CALL EXITS('IPINI6')
      RETURN
 9999 CALL ERRORS('IPINI6',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPINI6')
      RETURN 1
      END


