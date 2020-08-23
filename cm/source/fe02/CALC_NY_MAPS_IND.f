      SUBROUTINE CALC_NY_MAPS_IND(NKH,NP_INTERFACE,NPNODE,NPNY,nr,NVHP,
     '  nx,NYNP,NYNR,ERROR,*)

C#### Subroutine: CALC_NY_MAPS_IND
C###  Description:
C###    CALC_NY_MAPS_IND calculates the mapping arrays NYNP and NPNY
C###    for the `independent variables' ie. fitting.
C###  See-Also: CALC_NY_MAPS_DEP

C**** To keep this routine consistent with the structures used in
C**** CALC_NY_MAPS_DEP only the nc=1 variables are considered.
C**** See the routine CALC_NY_MAPS_DEP for a definition of the
C**** variables.

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NKH(NHM,NPM,NCM,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,
     '  NVHP(NHM,NPM,NCM,0:NRM),nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,nh,nhj,nhx,njj,nk,nonode,np,nrc,nrr,nv,ny,
     '  ny_start(0:2)
      LOGICAL CALC_NY

      CALL ENTERS('CALC_NY_MAPS_IND',*9999)

C***  Find the starting ny for the current region

      CALL ASSERT(NRCM.GE.2,'>>Increase NRCM to be >= 2',
     '  ERROR,*9999)

      DO nrc=0,2
        ny_start(nrc)=0
        DO nrr=1,nr-1
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO njj=1,NUM_FIT(0)
              DO nhj=1,NUM_FIT(njj)
                nhx=NLH_FIT(nhj,3,njj)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,nrr)
                  DO nk=1,NKH(nh,np,1,nrr)
                    ny_start(nrc)=ny_start(nrc)+1
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhj
            ENDDO !njj
          ENDDO !nonode (np)
          IF(nrc.NE.0) NYT(nrc,1,nx)=ny_start(nrc)
        ENDDO !nrr
      ENDDO !nrc

C***  Initialise mapping arrays above current region

      DO nrc=0,2
        ny=ny_start(nrc)
        DO nrr=nr,NRT
          DO nonode=1,NPNODE(0,nrr)
            np=NPNODE(nonode,nrr)
            DO njj=1,NUM_FIT(0)
              DO nhj=1,NUM_FIT(njj)
                nhx=NLH_FIT(nhj,3,njj)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,nrr)
                  DO nk=1,NKH(nh,np,1,nrr)
                    ny=ny+1
                    NYNP(nk,nv,nh,np,nrc,1,nr)=0
                    IF(ny.LE.NYM) THEN
                      DO i=0,6
                        NPNY(i,ny,nrc)=0
                      ENDDO
                    ENDIF !nym
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhj
            ENDDO !njj
          ENDDO !nonode (np)
          NYNR(0,nrc,1,nrr)=0
        ENDDO !nrr
      ENDDO !nrc

      CALL ASSERT(ny.LE.NYM,'>>Increase NYM',ERROR,*9999)

C***  Set up mapping arrays for current region.

      DO nrc=0,2
        ny=ny_start(nrc)
        NYNR(0,nrc,1,nr)=0
C KAT 16Nov98: Moving this
C        IF(nrc.NE.0) NYT(nrc,1,nx)=ny_start(nrc)
C        DO njj=1,NUM_FIT(0)
        IF(nrc.NE.0) ny=ny_start(nrc)
        DO njj=1,NUM_FIT(0)
          DO nhj=1,NUM_FIT(njj)
            nhx=NLH_FIT(nhj,3,njj)
            nh=NH_LOC(nhx,nx)
C GMH 22/10/95 Removed on the advice of Chris
C            IF(nrc.NE.0) ny=0 !reset the ny counter back to zero
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              CALC_NY=.TRUE.
              IF(NP_INTERFACE(np,0).GT.1) THEN !node is on an interface
C               Do nothing
              ENDIF
              IF(CALC_NY) THEN
                DO nv=1,NVHP(nh,np,1,nr)
                  DO nk=1,NKH(nh,np,1,nr)
                    ny=ny+1
C                    IF(nrc.NE.0) NYT(nrc,1,nx)=NYT(nrc,1,nx)+1
                    NYNR(0,nrc,1,nr)=NYNR(0,nrc,1,nr)+1
                    IF(NYNR(0,nrc,1,nr).LE.NY_R_M)
     '                NYNR(NYNR(0,nrc,1,nr),nrc,1,nr)=ny
                    IF(ny.LE.NYM) THEN
                      NYNP(nk,nv,nh,np,nrc,1,nr)=ny
                      NPNY(0,ny,nrc)=1 !mesh dof is node based
                      NPNY(1,ny,nrc)=nk
                      NPNY(2,ny,nrc)=nv
                      NPNY(3,ny,nrc)=nh
                      NPNY(4,ny,nrc)=np
                      NPNY(5,ny,nrc)=1
                      NPNY(6,ny,nrc)=nr
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDIF !calc_ny
            ENDDO !nonode
          ENDDO !nhj
        ENDDO !njj
        CALL ASSERT(NYNR(0,nrc,1,nr).LE.NY_R_M,'>>Increase NY_R_M',
     '    ERROR,*9999)
        IF(nrc.EQ.1.AND.KTYP8.NE.6) THEN
          CALL ASSERT(ny.LE.NYROWM,'>>Increase NYROWM',ERROR,*9999)
        ENDIF
        IF(nrc.NE.0) THEN
          CALL ASSERT(ny.LE.NYM,'>>Increase NYM',ERROR,*9999)
          NYT(nrc,1,nx)=ny
        ENDIF
      ENDDO !nrc

      CALL EXITS('CALC_NY_MAPS_IND')
      RETURN
 9999 CALL ERRORS('CALC_NY_MAPS_IND',ERROR)
      CALL EXITS('CALC_NY_MAPS_IND')
      RETURN 1
      END


