      SUBROUTINE CALC_NY_MAPS_IND_OPT(NJH,NKJ,
     '  NPNODE,NPNY,nr,NVJP,nx,NYNP,NYNR,ERROR,*)

C#### Subroutine: CALC_NY_MAPS_IND_OPT
C###  Description:
C###    CALC_NY_MAPS_IND_OPT calculates the mapping arrays NYNP and NPNY
C###    for the `independent variables' in fibre optimisation probs
C###  See-Also: CALC_NY_MAPS_IND
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
      INTEGER NJH(3),NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),nr,NVJP(NJM,NPM),
     '  nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,nh,nhx,nj,nk,nonode,np,nrc,nrr,nv,ny

      CALL ENTERS('CALC_NY_MAPS_IND_OPT',*9999)

C***  Find the starting ny for the current region

      CALL ASSERT(NRCM.GE.2,'>>Increase NRCM to be >= 2',
     '  ERROR,*9999)

C***  Initialise mapping arrays

      ny=0
      DO nrc=0,2
        DO nrr=nr,NRT
          DO nonode=1,NPNODE(0,nrr)
            np=NPNODE(nonode,nrr)
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              nj=NJH(nhx)
              DO nv=1,NVJP(nj,np)
                DO nk=1,NKJ(nj,np)
                  ny=ny+1
                  NYNP(nk,nv,nh,np,nrc,1,nr)=0
                  IF(ny.LE.NYM) THEN
                    DO i=0,6
                      NPNY(i,ny,nrc,nx)=0
                    ENDDO
                  ENDIF !nym
                ENDDO !nk
              ENDDO !nv
            ENDDO !nhx
          ENDDO !nonode (np)
          NYNR(0,nrc,1,nrr,nx)=0
        ENDDO !nrr
      ENDDO !nrc

      CALL ASSERT(ny.LE.NYM,'>>Increase NYM',ERROR,*9999)

C***  Set up mapping arrays for current region.

      ny=0
      DO nrc=0,2
        DO nrr=nr,NRT
          DO nonode=1,NPNODE(0,nrr)
            np=NPNODE(nonode,nrr)
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              nj=NJH(nhx)
              DO nv=1,NVJP(nj,np)
                DO nk=1,NKJ(nj,np)
                  ny=ny+1
                  NYNR(0,nrc,1,nrr,nx)=NYNR(0,nrc,1,nrr,nx)+1
                  IF(NYNR(0,nrc,1,nrr,nx).LE.NY_R_M)
     '              NYNR(NYNR(0,nrc,1,nrr,nx),nrc,1,nrr,nx)=ny
                  IF(ny.LE.NYM) THEN
                    NYNP(nk,nv,nh,np,nrc,1,nrr)=ny
                    NPNY(0,ny,nrc,nx)=1 !mesh dof is node based
                    NPNY(1,ny,nrc,nx)=nk
                    NPNY(2,ny,nrc,nx)=nv
                    NPNY(3,ny,nrc,nx)=nj
                    NPNY(4,ny,nrc,nx)=np
                    NPNY(5,ny,nrc,nx)=1
                    NPNY(6,ny,nrc,nx)=nrr
                  ENDIF
                ENDDO !nk
              ENDDO !nv
            ENDDO !nhx
          ENDDO !nonode (np)
        ENDDO !nrr
      ENDDO !nrc

      CALL ASSERT(NYNR(0,0,1,nr,nx).LE.NY_R_M,'>>Increase NY_R_M',
     '  ERROR,*9999)
      IF(nrc.EQ.1.AND.KTYP8.NE.6) THEN
        CALL ASSERT(ny.LE.NYROWM,'>>Increase NYROWM',ERROR,*9999)
      ENDIF
      IF(nrc.NE.0) THEN
        CALL ASSERT(ny.LE.NYM,'>>Increase NYM',ERROR,*9999)
        NYT(1,1,nx)=ny
      ENDIF

      CALL EXITS('CALC_NY_MAPS_IND_OPT')
      RETURN
 9999 CALL ERRORS('CALC_NY_MAPS_IND_OPT',ERROR)
      CALL EXITS('CALC_NY_MAPS_IND_OPT')
      RETURN 1
      END


