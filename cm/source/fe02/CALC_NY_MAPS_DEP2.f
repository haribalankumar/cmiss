      SUBROUTINE CALC_NY_MAPS_DEP2(NBH,NEELEM,NHP,NKH,
     '  NPNE,NPNODE,NPNY,nr,NVHP,nx,NYNE,NYNP,NYNR,
     '  ERROR,*)

C#### Subroutine: CALC_NY_MAPS_DEP2
C###  Description:
C###    CALC_NY_MAPS_DEP2 calculates the mapping arrays NYNE, NYNP and
C###    their inverse NPNY for dependent variables.
C###  This variation of CALC_NY_MAPS_DEP sets up NYNP & NYNE consecutively,
C###  i.e. ny's defined for the nodes of an element then the element itself.
C### NB/ NP_INTERFACE stuff taken out for now may be needed later
C###  KSB 24 October 2001

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHP(NPM,0:NRM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,
     '  NVHP(NHM,NPM,NCM,0:NRM),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,na,nb,nc,ne,nh,nhx,nk,NK_TOT,nn,
     '  noelem,nonode,np,nrr,nrc,nv,ny,ny_start(0:2,3)
      LOGICAL CALC_NY
      CHARACTER CHAR*1
!     External function
      INTEGER IDIGITS

      CALL ENTERS('CALC_NY_MAPS_DEP2',*9999)

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
                  IF(nrc.LE.1) THEN
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
                IF(nrc.LE.1) THEN
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
                      DO i=0,5
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
      IF(ny.GT.NYM) THEN
        WRITE(CHAR,'(I1)') IDIGITS(ny)
        WRITE(ERROR,'(''>>Increase NYM to '',I'//CHAR//')') ny
        GO TO 9999
      ENDIF

C***  Set up mapping arrays for current region.
C###  Altered from here includes nested loop creating NYNP & NYNE
C###  consecutively
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
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO nn=1,NNT(nb)
              np=NPNE(nn,nb,ne)
              CALC_NY=.TRUE.
C## KSB 10/01: Have removed INTERFACE stuff from this subroutine for now
              IF(CALC_NY) THEN
                DO nhx=1,NHP(np,nr)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc,nr)
c cpb 16/5/95 Removing (if necessary) cross derivatives from the ny
c calculation
                    IF(nrc.EQ.1) THEN
                      NK_TOT=MAX(NKH(nh,np,1,nr)-KTYP93(1,nr),
     '                  NKH(nh,np,2,nr)-KTYP93(2,nr),1)
                    ELSE
                      NK_TOT=MAX(NKH(nh,np,nc,nr)-KTYP93(nc,nr),1)
                    ENDIF
                    DO nk=1,NK_TOT
                      IF(NYNP(nk,nv,nh,np,nrc,nc,nr).EQ.0)THEN
                      !Checks if this node already been done
                        ny=ny+1
                        NYNR(0,nrc,nc,nr)=NYNR(0,nrc,nc,nr)+1
                        IF(nrc.NE.0) THEN
                          IF(ny.GT.NYT(nrc,nc,nx)) NYT(nrc,nc,nx)=ny
                        ENDIF
                        IF(NYNR(0,nrc,nc,nr).LE.NY_R_M)
     '                    NYNR(NYNR(0,nrc,nc,nr),nrc,nc,nr)=ny
                        IF(ny.LE.NYM) THEN
                          NYNP(nk,nv,nh,np,nrc,nc,nr)=ny
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
                          ENDIF !nrc.NE.1
                        ENDIF !ny.LE.NYM
                      ENDIF !NYNP.EQ.0
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nhx
              ENDIF !calc_ny
            ENDDO !nn (np)
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
C                 Don't set NPNY for nrc=1(rows), nc<>1(LHS) cases
                    NPNY(0,ny,nrc)=2 !mesh dof is element based
                    NPNY(1,ny,nrc)=na
                    NPNY(2,ny,nrc)=nh
                    NPNY(3,ny,nrc)=nc
                    NPNY(4,ny,nrc)=ne
                    NPNY(5,ny,nrc)=nr
                  ENDIF
                ENDIF
              ENDDO !na
            ENDDO !nhx
          ENDDO !noelem (ne)
          IF(NYNR(0,nrc,nc,nr).GT.NY_R_M) THEN
            WRITE(CHAR,'(I1)') IDIGITS(NYNR(0,nrc,nc,nr))
            WRITE(ERROR,'(''>>Increase NY_R_M to '',I'//CHAR//')')
     &        NYNR(0,nrc,nc,nr) 
            GO TO 9999
          ENDIF
          IF(nrc.EQ.1) THEN
            IF(ny.GT.NYROWM) THEN
              WRITE(CHAR,'(I1)') IDIGITS(ny)
              WRITE(ERROR,'(''>>Increase NYROWM to '',I'//CHAR//')') ny
              GO TO 9999
            ENDIF
          ENDIF
          IF(NYNR(0,nrc,nc,nr).GT.NYM) THEN
            WRITE(CHAR,'(I1)') IDIGITS(NYNR(0,nrc,nc,nr))
            WRITE(ERROR,'(''>>Increase NYM to '',I'//CHAR//')')
     &        NYNR(0,nrc,nc,nr)
            GO TO 9999
          ENDIF
          IF(ny.GT.NYM) THEN
            WRITE(CHAR,'(I1)') IDIGITS(ny)
            WRITE(ERROR,'(''>>Increase NYM to '',I'//CHAR//')') ny
            GO TO 9999
          ENDIF
        ENDDO !nc
      ENDDO !nrc

      CALL EXITS('CALC_NY_MAPS_DEP2')
      RETURN
 9999 CALL ERRORS('CALC_NY_MAPS_DEP2',ERROR)
      CALL EXITS('CALC_NY_MAPS_DEP2')
      RETURN 1
      END


