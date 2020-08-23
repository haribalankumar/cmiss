      SUBROUTINE CALC_VERSIONS_DEP(NBH,NBJ,NEELEM,
     '  NHE,NHP,NPNE,NPNODE,nr,NVHE,NVHP,NVJE,NVJP,nx,ERROR,*)

C#### Subroutine: CALC_VERSIONS_DEP
C###  Description:
C###    CALC_VERSIONS_DEP calculates the version for dependent
C###    variables stored in NVHP and NVHE arrays. Handles some
C###    problem specific version mappings; others are prompted
C###    for in IPEQUA.

C!!!  NOTE: this may need to change to IPVERSIONS_DEP (or sth) if
C!!!        version prompting is required in this routine

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM),NHP(NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nb1,nb2,nbj1,nc,ne,ne1,ne2,nh,nhx,nj,njj,njj1,nj_max,
     '  nn,no_aero,noelem,nonode,no_wake,no_wake_node,np

      CALL ENTERS('CALC_VERSIONS_DEP',*9999)

C     Initialise NVHP for nh's in this nx in current region
      DO nc=1,NCT(nr,nx)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C MPN 3Mar96 DO nhx=1,NH_LOC(0,nx)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            NVHP(nh,np,nc,nr)=0
          ENDDO !nh
        ENDDO !np
      ENDDO !nc

      IF(ITYP1(nr,nx).EQ.3) THEN !FE30 problems
C       Set default version arrays NVHP(nh,np,nc,nr) & NVHE(nn,nb,nh,ne)
        DO nc=1,NCT(nr,nx) !to set default to 1
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
C MPN 3Mar96 DO nhx=1,NH_LOC(0,nx)
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              NVHP(nh,np,nc,nr)=1
            ENDDO !nh
          ENDDO !nonode
        ENDDO !nc
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                NVHE(nn,nb,nh,ne)=1
              ENDDO !nn
            ENDDO !nb
          ENDDO !nh
        ENDDO !noelem

C       Setup special cases for version arrays
        IF(ITYP5(nr,nx).EQ.1.AND   !static analysis
     '    .ITYP2(nr,nx).EQ.3.AND   !Laplace equation
     '    .ITYP3(nr,nx).EQ.3) THEN !aerofoil analysis
C         Wake nodes - 2 versions for potential
          DO no_wake_node=1,NP_WAKE(0,1)
            NVHP(1,NP_WAKE(no_wake_node,1),1,nr)=2
            NVHP(1,NP_WAKE(no_wake_node,1),2,nr)=2
          ENDDO !no_wake_node

C         Specify node versions used in wake elements
          DO no_wake=1,NL_WAKE(0,1) !wake lines
            ne1=NE_WAKE(no_wake,1) !element# above wake line
            ne2=NE_WAKE(no_wake,2) !element# below wake line
            nh=NH_LOC(1,nx)
            nb1=NBH(nh,1,ne1) !nc=1
            nb2=NBH(nh,1,ne2) !nc=1
            NVHE(1,nb1,nh,ne1)=2 !nn=1 for element above wake line
            NVHE(2,nb1,nh,ne1)=2 !nn=2  "     "      "     "    "
            NVHE(3,nb2,nh,ne2)=1 !nn=3 for element below wake line
            NVHE(4,nb2,nh,ne2)=1 !nn=4  "     "      "     "    "
          ENDDO !no_wake

C         Specify node versions used in trailing edge aerofoil elements
          DO no_aero=1,NE_AERO(0,1) !aero lines
            ne1=NE_AERO(no_aero,1) !element# above aerofoil
            ne2=NE_AERO(no_aero,2) !element# below aerofoil
            nh=Nh_LOC(1,nx)
            nb1=NBH(nh,1,ne1) !nc=1
            nb2=NBH(nh,1,ne2) !nc=1
            IF(NPNE(2,nb1,ne1).EQ.NP_aero_TE1) THEN !2nd node is TE
              NVHE(2,nb1,nh,ne1)=2 !nn=2 for element above aerofoil
            ENDIF
            IF(NPNE(4,nb2,ne2).EQ.NP_aero_TE2) THEN !4th node is TE
              NVHE(4,nb2,nh,ne2)=1 !nn=4 for element below aerofoil
            ENDIF
          ENDDO !no_wake
        ELSE IF(JTYP2A.EQ.1) THEN
          IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3.
     '      AND.NJ_LOC(NJL_FIEL,0,nr).GE.2) THEN !eikonal equation
C           Use field 2 to determine the dependent variable versions.
            nj=NJ_LOC(NJL_FIEL,2,nr)
            DO nc=1,NCT(nr,nx)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(NVJP(nj,np).GT.1) THEN
                  DO nhx=1,NHP(np)
                    nh=NH_LOC(nhx,nx)
                    NVHP(nh,np,nc,nr)=NVJP(nj,np)
                  ENDDO !nh
                ENDIF
              ENDDO !nonode
            ENDDO !nc
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nbj1=NBJ(nj,ne)
              DO nhx=1,NHE(ne)
                nh=NH_LOC(nhx,nx)
                nb=NBH(nh,1,ne) !assuming nc=2 use same basis as nc=1
                DO nn=1,NNT(nb)
                  np=NPNE(nn,nb,ne)
                  IF(NVHP(nh,np,1,nr).GT.1) THEN
                    NVHE(nn,nb,nh,ne)=NVJE(nn,nbj1,nj,ne)
                  ENDIF
                ENDDO !nn
              ENDDO !nh
            ENDDO !noelem
          ELSE
C           Base dependent variable version on the geometric variable with
C           the largest number of versions.
            DO nc=1,NCT(nr,nx)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np)
                  nh=NH_LOC(nhx,nx)
                  DO njj1=NJL_GEOM,NJL_FIBR,NJL_FIBR-NJL_GEOM !geometry,fibre
                    DO njj=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj,nr)
                      IF(NVJP(nj,np).GT.NVHP(nh,np,nc,nr)
     '                  ) NVHP(nh,np,nc,nr)=NVJP(nj,np)
                    ENDDO !njj
                  ENDDO !njj1
                ENDDO !nh
              ENDDO !nonode
            ENDDO !nc
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nhx=1,NHE(ne)
                nh=NH_LOC(nhx,nx)
                nb=NBH(nh,1,ne) !assuming nc=2 use same basis as nc=1
                DO nn=1,NNT(nb)
                  np=NPNE(nn,nb,ne)
                  IF(NVHP(nh,np,1,nr).GT.1) THEN
                    njj1=NJL_GEOM !geometry
                    nj_max=NJ_LOC(njj1,1,nr)
                    DO njj=2,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj,nr)
                      IF(NVJP(nj,np).GT.NVJP(nj_max,np)) nj_max=nj
                    ENDDO
                    njj1=NJL_FIBR !fibre
                    DO njj=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj,nr)
                      IF(NBJ(nj,ne).NE.0.
     '                  AND.NVJP(nj,np).GT.NVJP(nj_max,np)) nj_max=nj
                    ENDDO
                    nbj1=NBJ(nj_max,ne)
                    NVHE(nn,nb,nh,ne)=NVJE(nn,nbj1,nj_max,ne)
                  ENDIF
                ENDDO !nn
              ENDDO !nh
            ENDDO !noelem
          ENDIF !problem type
        ENDIF !aerofoil/coincident version mapping

      ELSE IF(ITYP1(nr,nx).EQ.5) THEN  !FE50 problems
C!!!    NOTE the dependent variable versions for geometry default must
C!!!    be the same as independent variable versions.
C!!!    This may need to become problem dependent.
        DO nc=1,NCM !to set default to 1
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
C MPN 3Mar96  DO nhx=1,NH_LOC(0,nx)
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              IF(nhx.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
                nj=NJ_LOC(NJL_GEOM,nhx,nr)
                NVHP(nh,np,nc,nr)=NVJP(nj,np)
              ELSE
                NVHP(nh,np,nc,nr)=1
              ENDIF
            ENDDO !nh
          ENDDO !nonode
        ENDDO !nc
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
! GBS 04-04-1996            nbh1=NBH(nh,1,ne)
! Modified to initialise for all bases
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                IF(nhx.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
                  nj=NJ_LOC(NJL_GEOM,nhx,nr)
                  nbj1=NBJ(nj,ne)
                  NVHE(nn,nb,nh,ne)=NVJE(nn,nbj1,nj,ne)
                ELSE
                  NVHE(nn,nb,nh,ne)=1
                ENDIF
              ENDDO !nn
            ENDDO !nb
          ENDDO !nh
        ENDDO !noelem

      ELSE IF(ITYP1(nr,nx).EQ.6) THEN !FE60 problems
C!! The problem type in fe60 is to have geometric nodes spaced slightly
C!! apart from the actual Voronoi boundary. However, because of the
C!! finite volume staggered grid, the dependent variables are stored
C!! at the nodes, not the faces. Therefore, since these nodes placed
C!! at the boundary are versions 2,3 for external/internal nodes in the
C!! geometric sense, we need to put dependent variables at versions
C!! 1,2 for external,internal nodes.

C       Set up NVHP
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            IF(nhx.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
              nj=NJ_LOC(NJL_GEOM,nhx,nr)
              NVHP(nh,np,1,nr)=NVJP(nj,np)
            ELSE ! Pressure - more dependent variables than geometric
              NVHP(nh,np,1,nr)=NVJP(1,np)
            ENDIF
          ENDDO !nh
        ENDDO !nonode
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                IF(nhx.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
                  nj=NJ_LOC(NJL_GEOM,nhx,nr)
                  nbj1=NBJ(nj,ne)
                  NVHE(nn,nb,nh,ne)=NVJE(nn,nbj1,nj,ne)
                ELSE ! Pressure - more dependent variables than geomtric
                  NVHE(nn,nb,nh,ne)=NVJE(nn,nbj1,1,ne)
                ENDIF
              ENDDO !nn
            ENDDO !nb
          ENDDO !nh
        ENDDO !noelem

      ELSE !other problem types
C       Set default version arrays NVHP(nh,np,nc,nr) & NVHE(nn,nb,nh,ne)
        DO nc=1,NCM !to set default to 1
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
C MPN 3Mar96 DO nhx=1,NH_LOC(0,nx)
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              NVHP(nh,np,nc,nr)=1
            ENDDO !nh
          ENDDO !nonode
        ENDDO !nc
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            DO nb=1,NBFT
              DO nn=1,NNT(nb)
                NVHE(nn,nb,nh,ne)=1
              ENDDO !nn
            ENDDO !nb
          ENDDO !nh
        ENDDO !noelem
      ENDIF

C     Find max NVHP for zero'th nr location
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        DO nc=1,NCT(nr,nx)
          DO nhx=1,NH_LOC(0,nx)
            nh=NH_LOC(nhx,nx)
            IF(NVHP(nh,np,nc,nr).GT.NVHP(nh,np,nc,0))
     '        NVHP(nh,np,nc,0)=NVHP(nh,np,nc,nr)
          ENDDO !nh
        ENDDO !nc
      ENDDO !nonode (np)

      CALL EXITS('CALC_VERSIONS_DEP')
      RETURN
 9999 CALL ERRORS('CALC_VERSIONS_DEP',ERROR)
      CALL EXITS('CALC_VERSIONS_DEP')
      RETURN 1
      END


