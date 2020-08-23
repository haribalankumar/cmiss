      SUBROUTINE CALC_NY_MAPS_DEP(NBH,NEELEM,NHP,NKH,
     '  NP_INTERFACE,NPNODE,NPNY,nr,NVHP,nx,NYNE,NYNP,NYNR,
     '  ERROR,*)

C#### Subroutine: CALC_NY_MAPS_DEP
C###  Description:
C###    CALC_NY_MAPS_DEP calculates the mapping arrays NYNE, NYNP and
C###    their inverse NPNY for dependent variables.

C#### Variable: NPNY(0:6,ny,0:nrc,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NY_MAPS_DEP
C###  Description:
C###    <HTML>
C###    <P>
C###    NPNY(0:6,ny,0,nx) for problem nx is the reverse mapping to NYNP
C###    and NYNE for global ny variables only.  For nrc=1 the array is
C###    not set up, and for nrc=2 it is only set up for matrix nc=1 and
C###    so should probably not be used.</P>
C###    <P>
C###    The NPNY(0,ny,nrc) position identifies whether
C###    the ny (for nrc) is node based (=1) or element based (=2).</P>
C###    <PRE>
C###    If NPNY(0,ny,nrc)=1 then:
C###      NPNY(1,ny,nrc)=nk
C###      NPNY(2,ny,nrc)=nv
C###      NPNY(3,ny,nrc)=nh
C###      NPNY(4,ny,nrc)=np
C###      NPNY(5,ny,nrc)=nc
C###      NPNY(6,ny,nrc)=nr
C###    If NPNY(0,ny,nrc)=2 then:
C###      NPNY(1,ny,nrc)=na
C###      NPNY(2,ny,nrc)=nh
C###      NPNY(3,ny,nrc)=nc
C###      NPNY(4,ny,nrc)=ne
C###      NPNY(5,ny,nrc)=nr
C###    </PRE> </HTML>

C#### Variable: nrc
C###  Type: INTEGER
C###  Set_up: CALC_NY_MAPS_DEP
C###  Description:
C###    nrc indicates whether an associated ny value is a row or column
C###    number.  A value of 1 indicates an equation (or row) number.
C###    Where a value of 0 is possible (NPNY,NYNP,NYNE,NYNP), 0 indicates
C###    the global variable (mesh dof) number, and 2 the column number
C###    local to the matrix nc.  Where only values of 1 and 2 are possible
C###    (NYNO, NONY), 2 indicates the global variable number.

C#### Variable: NYNE(na,nh,0:nrc,nc,ne)
C###  Type: INTEGER
C###  Set_up: CALC_NY_MAPS_DEP
C###  Description:
C###    NYNE(na,nh,nrc,nc,ne) for problem nx is the mapping between
C###    the auxillary parameters na,nh,nc,ne and mesh dof ny for
C###    equations/rows (nrc=1), variables (nrc=0), and columns (nrc=2).

C#### Variable: NYNP(nk,nv,nh,np,0:nrc,nc,nr)
C###  Type: INTEGER
C###  Set_up: CALC_NY_MAPS_DEP
C###  Description:
C###    NYNP(nk,nv,nh,np,0:nrc,nc,nr) for problem nx is the mapping
C###    between the mesh parameters nk,nv,nh,np,nc,nr and mesh dof ny.
C###    nrc=1/2 returns a local row[equation]/column[variable] number
C###    for the matrix nc. nrc=0 gives the global
C###    mesh dof (ny or variable number) for region nr.

C#### Variable: NYNR(0:ny,0:nrc,nc,0:nr,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NY_MAPS_DEP
C###  Description:
C###    NYNR(0:ny,0:nrc,nc,0:nr,nx) is the list of nys for a region.
C###    NYNR(0:ny,0,nc,0:nr,nx) gives the global list of nys (dofs)
C###    for a region.  NYNR(0:ny,1,nc,0:nr,nx) gives the list of local
C###    [=global in all cases] rows for region nr and matrix nc.
C###    NYNR(0:ny,2,nc,0:nr,nx) gives the list of local column numbers
C###    for matrix nc in region nr. NYNR(0,nrc,nc,nr,nx) is the number
C###    in the list and NYNR(1..,nrc,nc,nr,nx) is the list of nys.
C###    nr=0 gives the mappings for the entire coupled problem.
C###    NOTE: This is not set up until after the solution mapping
C###    arrays have been set up.

C#### Variable: NYT(nrc,nc,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NY_MAPS_DEP
C###  Description:
C###    NYT(nrc,nc,nx) is the total number of ny values for rows (nrc=1)
C###    and columns (nrc=2) for matrix nc and problem type nx.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHP(NPM,0:NRM),
     '  NKH(NHM,NPM,NCM,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,
     '  NVHP(NHM,NPM,NCM,0:NRM),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,na,nb,nc,ne,nh,nhx,nk,NK_TOT,
     '  noelem,nonode,np,nrr,nrc,nv,ny,ny_start(0:2,3)
      CHARACTER STRING*255
      LOGICAL CALC_NY

      CALL ENTERS('CALC_NY_MAPS_DEP',*9999)

C***  Find the starting ny for the current region

      CALL ASSERT(NRCM.GE.2,'>>Increase NRCM to be >= 2',
     '  ERROR,*9999)

      DO nrc=0,2
        DO nc=1,NCT(nr,nx)
          ny_start(nrc,nc)=0
        ENDDO
        DO nrr=1,nr-1
          DO nc=1,NCT(nrr,nx) !LHS,RHS and GD variables
            DO nonode=1,NPNODE(0,nrr)
              np=NPNODE(nonode,nrr)
              DO nhx=1,NHP(np,nrr)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc,nrr)

c cpb 16/5/95 Removing (if necessary) cross derivatives from the ny
c calculation

C ??? cpb 26/2/94 Why is this done for nrc <= 1 here but only nrc = 1
C ??? below.
                  IF(nrc.LE.1) THEN
!CC                  IF(nrc.EQ.1) THEN
                    NK_TOT=MAX(NKH(nh,np,1,nrr)-KTYP93(1,nrr),
     '                NKH(nh,np,2,nrr)-KTYP93(2,nrr),1)
                  ELSE
                    NK_TOT=MAX(NKH(nh,np,nc,nrr)-KTYP93(nc,nrr),1)
                  ENDIF
                  DO nk=1,NK_TOT
                    ny_start(nrc,nc)=ny_start(nrc,nc)+1
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDDO !nonode (np)
            DO noelem=1,NEELEM(0,nrr)
              ne=NEELEM(noelem,nrr)
              DO nhx=1,NH_LOC(0,nx)
                nh=NH_LOC(nhx,nx)
C ??? cpb 26/2/94 Why is this done for nrc <= 1 here but only nrc = 1
C ??? below.
                IF(nrc.LE.1) THEN
!CC                IF(nrc.EQ.1) THEN
                  nb=NBH(nh,1,ne) !Use nc=1 to determine # of eqn's
                ELSE
                  nb=NBH(nh,nc,ne)
                ENDIF
! AJP 19/7/96 Why should nb be set up for dep var nh in another region?
! Hence the new "IF" statement
                IF(nb.GT.0) THEN
                  DO na=1,NAT(nb)
                    ny_start(nrc,nc)=ny_start(nrc,nc)+1
                  ENDDO !na
                ENDIF
              ENDDO !nh
            ENDDO !noelem (ne)
          ENDDO !nc
        ENDDO !nrr
      ENDDO !nrc

C***  Initialise mapping arrays above current region

      DO nrc=0,2
        IF(nrc.EQ.0) THEN
          ny=0
          DO nc=1,NCT(nr,nx)  !LHS,RHS and GD variables
            ny=ny+ny_start(nrc,nc)
          ENDDO
        ENDIF
        DO nrr=nr,NRT
          DO nc=1,NCT(nrr,nx) !LHS,RHS and GD variables
            IF(nrc.NE.0) ny=ny_start(nrc,nc)
            DO nonode=1,NPNODE(0,nrr)
              np=NPNODE(nonode,nrr)
C              IF(NP_INTERFACE(np,0).GT.0) THEN
C                INTERFACE=.TRUE.
C              ELSE
C                INTERFACE=.FALSE.
C              ENDIF
              DO nhx=1,NHP(np,nrr)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc,nrr)
c cpb 16/5/95 Removing (if necessary) cross derivatives from the ny
c calculation
                  IF(nrc.EQ.1) THEN
                    NK_TOT=MAX(NKH(nh,np,1,nrr)-KTYP93(1,nrr),
     '                NKH(nh,np,2,nrr)-KTYP93(2,nrr),1)
                  ELSE
                    NK_TOT=MAX(NKH(nh,np,nc,nrr)-KTYP93(nc,nrr),1)
                  ENDIF
                  DO nk=1,NK_TOT
                    NYNP(nk,nv,nh,np,nrc,nc,nr)=0
                    ny=ny+1
C                    IF(INTERFACE) THEN
C                      IF(nrr.EQ.NP_INTERFACE(np,1)) THEN
C                        NYNP(nk,nv,nh,np,nrc,nc)=0
C                      ENDIF
C                    ELSE
C                      NYNP(nk,nv,nh,np,nrc,nc)=0
C                    ENDIF
                    IF(ny.LE.NYM) THEN
                      DO i=0,6
                        NPNY(i,ny,nrc)=0
                      ENDDO
                    ENDIF !nym
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDDO !nonode (np)

            DO noelem=1,NEELEM(0,nrr)
              ne=NEELEM(noelem,nrr)
              DO nhx=1,NH_LOC(0,nx)
                nh=NH_LOC(nhx,nx)
                IF(nrc.EQ.1) THEN
                  nb=NBH(nh,1,ne) !Use nc=1 to determine # of eqn's
                ELSE
                  nb=NBH(nh,nc,ne)
                ENDIF
! AJP 19/7/96 Why should nb be set up for dep var nh in another region?
! Hence the new "IF" statement
                IF(nb.GT.0) THEN
                  DO na=1,NAT(nb)
                    ny=ny+1
                    NYNE(na,nh,nrc,nc,ne)=0
                    IF(ny.LE.NYM)THEN
C JPP 4 July 2003
C CMISS bombs out in IOHIST if we don't initialize the 6'th row!
C                      DO i=0,5
                      DO i=0,6
                        NPNY(i,ny,nrc)=0
                      ENDDO
                    ENDIF !nym
                  ENDDO !na
                ENDIF !check on nb
              ENDDO !nh
            ENDDO !noelem (ne)
            NYNR(0,nrc,nc,nrr)=0
          ENDDO !nc
        ENDDO !nrr
      ENDDO !nrc
      CALL ASSERT(ny.LE.NYM,'>>Increase NYM',ERROR,*9999)

C***  Set up mapping arrays for current region.
      DO nrc=0,2
        IF(nrc.EQ.0) THEN
          ny=0
          DO nc=1,NCT(nr,nx) !LHS and RHS and GD variables
            ny=ny+ny_start(nrc,nc)
          ENDDO
        ENDIF
        DO nc=1,NCT(nr,nx) !LHS and RHS and GD variables
          IF(nrc.NE.0) ny=ny_start(nrc,nc)
          NYNR(0,nrc,nc,nr)=0
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            CALC_NY=.TRUE.
c            IF(NP_INTERFACE(np,0).GT.1) THEN !node is on an interface
c             INTERFACE=.TRUE.
cC*** This is problem/equation dependent in different regions. If there
cC*** are multiple ny's for a given node on the interface then each
cC*** problem will need to work out if another ny does not need to be
cC*** set up. The default is that both ny's will be set up for nc=1,2.
cC*** If not the problem will need to set CALC_NY to false.
c            ELSE
c             INTERFACE=.FALSE.
c            ENDIF
            IF(CALC_NY) THEN
              DO nhx=1,NHP(np,nr)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc,nr)
c cpb 16/5/95 Removing (if necessary) cross derivatives from the ny
c calculation
                  IF(nrc.EQ.1) THEN
                    NK_TOT=MAX(NKH(nh,np,1,nr)-KTYP93(1,nr),
     '                NKH(nh,np,2,nr)-KTYP93(2,nr),1)
                  ELSE
                    NK_TOT=MAX(NKH(nh,np,nc,nr)-KTYP93(nc,nr),1)
                  ENDIF
                  DO nk=1,NK_TOT
                    ny=ny+1
                    NYNR(0,nrc,nc,nr)=NYNR(0,nrc,nc,nr)+1
                    IF(nrc.NE.0) THEN
                      IF(ny.GT.NYT(nrc,nc,nx)) NYT(nrc,nc,nx)=ny
                    ENDIF
                    IF(NYNR(0,nrc,nc,nr).LE.NY_R_M)
     '                NYNR(NYNR(0,nrc,nc,nr),nrc,nc,nr)=ny
                    IF(ny.LE.NYM) THEN
                      NYNP(nk,nv,nh,np,nrc,nc,nr)=ny
C cpb 15/6/95 adding nr to NYNP
C                      IF(INTERFACE) THEN
C                        IF(nr.EQ.NP_INTERFACE(np,1)) THEN
CC!!! As the NYNP array is not specific for a region then at an
CC!!! interface there is multiple ny's for a given nk, nv, nh and np.
CC!!! By definition only the ny in the master (lowest region #) is
CC!!! stored in NYNP.
C                          NYNP(nk,nv,nh,np,nrc,nc)=ny
C                        ENDIF
C                      ELSE
C                        NYNP(nk,nv,nh,np,nrc,nc)=ny
C                      ENDIF
C cpb 1/3/95 Adding nrc to NPNY
                      IF(nrc.NE.1.OR.nc.EQ.1) THEN
C                       Don't set NPNY for nrc=1(rows), nc<>1(LHS) cases
                        NPNY(0,ny,nrc)=1 !mesh dof is nodal based
                        NPNY(1,ny,nrc)=nk
                        NPNY(2,ny,nrc)=nv
                        NPNY(3,ny,nrc)=nh
                        NPNY(4,ny,nrc)=np
                        NPNY(5,ny,nrc)=nc
                        NPNY(6,ny,nrc)=nr
                      ENDIF
                    ENDIF !ny<NYM
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDIF !calc_ny
          ENDDO !nonode (np)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              IF(nrc.EQ.1) THEN
                nb=NBH(nh,1,ne) !Use nc=1 to determine # of eqn's
              ELSE
                nb=NBH(nh,nc,ne)
              ENDIF
              DO na=1,NAT(nb)
                ny=ny+1
                NYNR(0,nrc,nc,nr)=NYNR(0,nrc,nc,nr)+1
                IF(nrc.NE.0) THEN
                  IF(ny.GT.NYT(nrc,nc,nx)) NYT(nrc,nc,nx)=ny
                ENDIF
                IF(NYNR(0,nrc,nc,nr).LE.NY_R_M)
     '            NYNR(NYNR(0,nrc,nc,nr),nrc,nc,nr)=ny
                IF(ny.LE.NYM) THEN
                  NYNE(na,nh,nrc,nc,ne)=ny
C cpb 1/3/95 Adding nrc to NPNY
                  IF(nrc.NE.1.OR.nc.EQ.1) THEN
C                   Don't set NPNY for nrc=1(rows), nc<>1(LHS) cases
                    NPNY(0,ny,nrc)=2 !mesh dof is element based
                    NPNY(1,ny,nrc)=na
                    NPNY(2,ny,nrc)=nh
                    NPNY(3,ny,nrc)=nc
                    NPNY(4,ny,nrc)=ne
                    NPNY(5,ny,nrc)=nr
                  ENDIF
                ENDIF
              ENDDO !na
            ENDDO !nh
          ENDDO !noelem (ne)
          CALL ASSERT(NYNR(0,nrc,nc,nr).LE.NY_R_M,'>>Increase NY_R_M',
     '      ERROR,*9999)
          IF(nrc.EQ.1) THEN
            CALL ASSERT(ny.LE.NYROWM,'>>Increase NYROWM',ERROR,*9999)
          ENDIF
          CALL ASSERT(NYNR(0,nrc,nc,nr).LE.NYM,'>>Increase NYM',
     '      ERROR,*9999)
          WRITE(STRING,'(''>>Increase NYM. Try '''
     &      //',I6)') ny
          CALL ASSERT(ny.LE.NYM,STRING,ERROR,*9999)
        ENDDO !nc
      ENDDO !nrc

      CALL EXITS('CALC_NY_MAPS_DEP')
      RETURN
 9999 CALL ERRORS('CALC_NY_MAPS_DEP',ERROR)
      CALL EXITS('CALC_NY_MAPS_DEP')
      RETURN 1
      END

