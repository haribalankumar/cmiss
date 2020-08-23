      SUBROUTINE UPDELA(IBT,NBJ,NEELEM,NELIST,NENP,NKJE,NPLIST,NPNE,
     '  NPNODE,NRE,NVJE,NVJP,NXI,SE,XP,ZA,STRING,ERROR,*)

C#### Subroutine: UPDELA
C###  Description:
C###    Makes triangulation Delaunay

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter list
      INTEGER IBT(3,NIM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER noelem,ne,njj,IBEG,IEND,nb,IFROMC,N3CO,nonode,
     '  np,nn,nr,nj,ERROR_FLAG,ERROR_STRINGC(500),LENGTH
      CHARACTER ERROR_STRINGF*100
      LOGICAL CBBREV

      CALL ENTERS('UPDELA',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update delaunay
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Updates the delaunay triangulation

        OP_STRING(1)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE
        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

C       Do all the assert stuff
        CALL ASSERT(NEELEM(0,nr).GT.1,
     '    '>>Need more than 1 element to be able to flip',
     '    ERROR,*9999)
        CALL ASSERT(NHM.GE.(NJ_LOC(NJL_GEOM,0,nr)+1),
     '    '>>Increase NHM',ERROR,*9999)
        IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
          DO noelem=1,NEELEM(0,nr)-1
            ne=NEELEM(noelem,nr)
            CALL ASSERT(ne.EQ.(NEELEM(noelem+1,nr)-1),
     '        '>>3d flipping has to have sequential elements',
     '        ERROR,*9999)
          ENDDO
        ENDIF
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          CALL ASSERT(NNT(nb).EQ.(NJ_LOC(NJL_GEOM,0,nr)+1).AND.
     '      IBT(1,1,nb).EQ.3,
     '      '>>Need to specify linear simplex basis functions',
     '      ERROR,*9999)
        ENDDO
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            CALL ASSERT(NVJP(nj,np).LE.2,'>>Can only have one '//
     '        'version for Delaunay triangulations',ERROR,*9999)
          ENDDO
        ENDDO

C       End of asserts
        CALL CALC_NENP_VORO(NBJ,NEELEM,NENP,NPNE,NPNODE,nr,ERROR,*9999)
        CALL NENXI_VORO(NBJ,NEELEM,NENP,NPNE,NPNODE,nr,NXI,ERROR,*9999)

        ERROR_FLAG=0
        IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN

C         Obtain the nonode numbers foreach node number np
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            NPLIST(np)=nonode
          ENDDO
C         Obtain the noelem numbers foreach  number ne
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NELIST(ne)=noelem
          ENDDO
          CALL ALE_initialise(NAM,NBJ(1,NEELEM(1,nr)),NBFM,NCM,NELIST,
     '      NEM,NEELEM(0,nr),%VAL(NEIM),NHM,NIM,NJM,NKM,NNM,NPLIST,NPM,
     '      NPNE,NPNODE(0,nr),NVJE,NVM,NXI,XP,ZA,ZERO_TOL,ERROR_FLAG,
     '      ERROR_STRINGC)
        ENDIF
        IF(ERROR_FLAG.NE.0) THEN
          CALL CSTRINGLEN(LENGTH,ERROR_STRINGC)
          CALL C2FSTRING(ERROR_STRINGC,LENGTH,ERROR_STRINGF)
          CALL ASSERT(.FALSE.,ERROR_STRINGF,ERROR,*9999)
        ENDIF

        CALL RECONNECT(NBJ,NEELEM,NELIST,NENP,NKJE,NPLIST,NPNE,NPNODE,
     '    nr,NRE,NVJE,NXI,SE,XP,ZA,ERROR,*9999)

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
          WRITE(OP_STRING,'('' Simplex/Voronoi element'//
     '      ' listing:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            WRITE(OP_STRING,'('' Element  Node Version   Opp'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            DO nn=1,NNT(nb)
              WRITE(OP_STRING,'(''  '',I6,I6,''  '',I6,I6)') ne,
     '          NPNE(nn,nb,ne),NVJE(nn,nb,1,ne),
     '          NXI(0,nn,ne)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nn
          ENDDO !noelem
          WRITE(OP_STRING,'('' Element adjacency'//
     '      ' listing:'')')
          WRITE(OP_STRING,'(''  Node Adjacent Elements'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO noelem=1,NENP(np,0,nr)
              ne=NENP(np,noelem,nr)
              WRITE(OP_STRING,'(I6,I6)') np,ne
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !noelem
          ENDDO !nonode
CC$      call mp_unsetlock()
        ENDIF !dop
      ENDIF! CO(noco)

      CALL EXITS('UPDELA')
      RETURN
 9999 CALL ERRORS('UPDELA',ERROR)
      CALL EXITS('UPDELA')
      RETURN 1
      END


C      SUBROUTINE UPELEM(IBT,IDO,INP,NBJ,NEELEM,NENQ,NKE,NPF,NPNE,
C     '  NQLIST,NRLIST,NVJE,NWQ,NXLIST,NYNP,SE,XA,XE,XIQ,XP,XQ,STRING,
C     '  FIX,FIXQ,ERROR,*)
C
CC#### Subroutine: UPELEM
CC###  Description:
CC###    Updates element variables
CC**** Written by Martin Buist 4-June-1999
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:b00.cmn'
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:cbdi10.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:grid00.cmn'
C      INCLUDE 'cmiss$reference:mxch.inc'
C      INCLUDE 'cmiss$reference:loc00.cmn'
C      INCLUDE 'cmiss$reference:loc00.inc'
C      INCLUDE 'cmiss$reference:time02.cmn'
C
C!     Parameter List
C      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
C     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),
C     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
C     '  NQLIST(0:NQM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),
C     '  NWQ(8,0:NQM,NAM),NXLIST(0:NXM),
C     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
C      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XIQ(NIM,NQM),
C     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM)
C      CHARACTER ERROR*(*),STRING*(MXCH)
C      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXQ(NYQM,NIYFIXM,NXM)
C!     Local variables
C      INTEGER IBEG,IEND,IT,ITMAX,nb,nb_blood,nb_torso,nc,ne_bem,ne_t,
C     '  nh,ni,NITB,nj,nn,NODES(27),noelemt,np,nr,nr_blood,nr_grid,
C     '  nr_torso,nq,nqq,nv,nx,nx_blood,nxc,nx_torso,nx_upd,ny,N3CO
C      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
C      REAL*8 DISTANCE,DQ,DQMAX,XI(3),XP1(3)
C      LOGICAL ALL_REGIONS,BLOOD,CBBREV,ERROR_FLAG,POTENTIAL,
C     '  FLUX
CCC*s
CC      LOGICAL SPECIAL
CCC*f
C
C      CALL ENTERS('UPELEM',*9999)
C
C 1    IF(CO(noco+1).EQ.'?') THEN
C        CALL STRING_TRIM(STRING,IBEG,IEND)
C
CC---------------------------------------------------------------------
C
CC#### Command: FEM update elements
CC###  Parameter:      <region #s[1,2,3]>
CC###    Three regions are required, region 1 is the region in which
CC###    grid points are defined. Region 2 is the ventricle region and
CC###    region 3 is the torso cavity region.
CC###  Parameter:      <class #s[1,2,3]>
CC###    Three classes are required, class 1 is the extracellular
CC###    update class, class 2 is the ventricle solution class and
CC###    class 3 is the torso cavity-body surface solution class.
CC###  Parameter:      <grids name[none]>
CC###    This is an optional parameter which allows a subset of
CC###    grid points to be updated. Only points in the grid group
CC###    will be affected.
CC###  Description:
CC###    Generate a boundary element number and xi location
CC###    for external grid points in coupled FD/BEM problems
C
C        OP_STRING(1)=STRING(1:IEND)//' '
C        OP_STRING(2)=BLANK(1:15)//'<region #s[1,2,3]>'
C        OP_STRING(3)=BLANK(1:15)//'<class #s[1,2,3]>'
C        OP_STRING(4)=BLANK(1:15)//'<grids name[none]>'
C        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C
CC---------------------------------------------------------------------
C
C      ELSE IF(CO(noco+1).EQ.'??') THEN
C        CALL DOCUM('fe23','doc','UPELEM',ERROR,*9999)
C      ELSE
C
C        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
C        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
C
C        IF(CBBREV(CO,'REGION',3,noco+1,NTCO,n3co)) THEN
C          CALL ASSERT(NRLIST(0).EQ.3,' >>Must define 3 regions',
C     '      ERROR,*9999)
C          nr_grid=NRLIST(1)
C          nr_blood=NRLIST(2)
C          nr_torso=NRLIST(3)
C        ELSE
C          nr_grid=1
C          nr_blood=2
C          nr_torso=3
C        ENDIF
C
C        IF(CBBREV(CO,'CLASS',2,noco+1,NTCO,n3co)) THEN
C          CALL ASSERT(NXLIST(0).EQ.3,' >>Must define 3 classes',
C     '      ERROR,*9999)
C          nxc=NXLIST(1)
C          CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
C          CALL ASSERT(nx_upd.GT.0,
C     '      '>>No nx defined for this solve class',ERROR,*9999)
C          nxc=NXLIST(2)
C          CALL NX_LOC(NX_INQUIRE,nxc,nx_blood,NX_SOLVE,ERROR,*9999)
CC may be zero
CC            CALL ASSERT(nx_blood.GT.0,
CC     '        '>>No nx defined for this solve class',ERROR,*9999)
C          nxc=NXLIST(3)
C          CALL NX_LOC(NX_INQUIRE,nxc,nx_torso,NX_SOLVE,ERROR,*9999)
CC may be zero
CC            CALL ASSERT(nx_torso.GT.0,
CC     '        '>>No nx defined for this solve class',ERROR,*9999)
C        ELSE
C          nxc=1
C          CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
C          CALL ASSERT(nx_upd.GT.0,
C     '      '>>No nx defined for this solve class',ERROR,*9999)
C          nxc=2
C          CALL NX_LOC(NX_INQUIRE,nxc,nx_blood,NX_SOLVE,ERROR,*9999)
C          CALL ASSERT(nx_blood.GT.0,
C     '      '>>No nx defined for this solve class',ERROR,*9999)
C          nxc=3
C          CALL NX_LOC(NX_INQUIRE,nxc,nx_torso,NX_SOLVE,ERROR,*9999)
C          CALL ASSERT(nx_torso.GT.0,
C     '      '>>No nx defined for this solve class',ERROR,*9999)
C        ENDIF
C
C        !Initialise
C        DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
C          FIXQ(nq,1,nx_upd)=.FALSE.
C          FIXQ(nq,2,nx_upd)=.FALSE.
C          FIXQ(nq,3,nx_upd)=.FALSE.
C        ENDDO !nq
C
C        IF(CBBREV(CO,'GRIDS',2,noco+1,NTCO,n3co)) THEN
C          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
C            FIXQ(nq,3,nx_upd)=.TRUE.
C          ENDDO !nq
C          CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)
C          DO nqq=1,NQLIST(0)
C            nq=NQLIST(nqq)
C            FIXQ(nq,3,nx_upd)=.FALSE.
C          ENDDO !nqq
C        ELSE
C          NQLIST(0)=0
C          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
C            IF(nr_blood.EQ.0) THEN
C              IF(NWQ(1,nq,1).LT.nq) THEN
C                NQLIST(0)=NQLIST(0)+1
C                NQLIST(NQLIST(0))=nq
C              ELSE
C                FIXQ(nq,3,nx_upd)=.TRUE.
C              ENDIF
C            ELSE IF(nr_torso.EQ.0) THEN
C              IF(NWQ(1,nq,1).GT.nq) THEN
C                NQLIST(0)=NQLIST(0)+1
C                NQLIST(NQLIST(0))=nq
C              ENDIF
C            ELSE
C              NQLIST(0)=NQLIST(0)+1
C              NQLIST(NQLIST(0))=nq
C            ENDIF
C          ENDDO !nq
C        ENDIF !grid_group
C
C        !init
C        ERROR_FLAG=.FALSE.
CC          ITMAX=10
C        ITMAX=50
C        DO nq=1,NQT
CC MLB change for march8_coupled 7/4/00 to flux bc.
CC          IF(NWQ(1,nq,1).GT.0) FIXQ(nq,1,nx_upd)=.TRUE.
C          IF(NWQ(1,nq,1).GT.0) FIXQ(nq,2,nx_upd)=.TRUE.
C        ENDDO
C
C        CALL CPU_TIMER(CPU_USER,TIME_START)
C
CC$OMP PARALLEL DO
CC$&   PRIVATE(BLOOD,DISTANCE,DQ,DQMAX,FLUX,IT,nb,nb_blood,nb_torso,nc,
CC$&     ne_bem,ne_t,nh,ni,NITB,nj,nn,NODES,noelemt,np,nq,nqq,nr,nv,nx,
CC$&     ny,POTENTIAL,XE,XI,XP1)
CC$&   SHARED(ERROR_FLAG,FIX,FIXQ,IBT,IDO,INP,ITMAX,NBJ,NEELEM,NENQ,NKE,
CC$&     NPF,NPNE,NQLIST,NQR,nr_blood,nr_grid,nr_torso,NVJE,NWQ,nx_upd,
CC$&     NYNP,SE,XA,XIQ,XP,XQ)
CC          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
C        DO nqq=1,NQLIST(0)
C          nq=NQLIST(nqq)
C          IF(.NOT.ERROR_FLAG) THEN
C
C            !initialise to -ve so we can check later
C            DO ni=1,NIM
C              XIQ(ni,nq)=-1.0d0
C            ENDDO
C            IF(NWQ(1,nq,1).NE.0) THEN !boundary point
C
C              !calculate ne_bem
C              DQMAX=RMAX
C              ne_bem=0
C              BLOOD=.FALSE.
C
C              !check elements from the torso region
C              IF(nr_torso.GT.0) THEN
C                DO noelemt=1,NEELEM(0,nr_torso)
C                  ne_t=NEELEM(noelemt,nr_torso)
C                  nb_torso=NBJ(1,ne_t)
C                  DQ=0.0d0
C                  DO nn=1,NNT(nb_torso)
C                    np=NPNE(nn,nb_torso,ne_t)
C                    DO nj=1,NJT
C                      DQ=DQ+(XQ(nj,nq)-XP(1,1,nj,np))**2.0d0
C                    ENDDO
C                  ENDDO
C                  DQ=DSQRT(DQ)
C                  IF(DQ.LT.DQMAX) THEN
C                    DQMAX=DQ
C                    ne_bem=ne_t
C                  ENDIF
C                ENDDO
C              ENDIF
C
C              !check elements from the blood region
C              IF(nr_blood.GT.0) THEN
C                DO noelemt=1,NEELEM(0,nr_blood)
C                  ne_t=NEELEM(noelemt,nr_blood)
C                  nb_blood=NBJ(1,ne_t)
C                  DQ=0.0d0
C                  DO nn=1,NNT(nb_blood)
C                    np=NPNE(nn,nb_blood,ne_t)
C                    DO nj=1,NJT
C                      DQ=DQ+(XQ(nj,nq)-XP(1,1,nj,np))**2.0d0
C                    ENDDO
C                  ENDDO
C                  DQ=DSQRT(DQ)
C                  IF(DQ.LT.DQMAX) THEN
C                    DQMAX=DQ
C                    ne_bem=ne_t
C                    BLOOD=.TRUE.
C                  ENDIF
C                ENDDO
C              ENDIF
C
C              CALL ASSERT(ne_bem.GT.0,'>>No torso element found',
C     '          ERROR,*100)
C
C              !store element in NENQ
C              CALL ASSERT(NENQ(0,nq).LE.7,
C     '          '>>Increase NENQ allocation',ERROR,*100)
C              NENQ(0,nq)=NENQ(0,nq)+1
C              NENQ(NENQ(0,nq),nq)=ne_bem
C
C              !calculate xi within ne_bem
C              IF(BLOOD) THEN
C                nb=nb_blood
C                nr=nr_blood
C                nx=nx_blood
C              ELSE
C                nb=nb_torso
C                nr=nr_torso
C                nx=nx_torso
C              ENDIF
C
C              nb=NBJ(1,ne_bem)
C              NITB=NIT(nb)
C              DO ni=1,NITB !initialising XI
C                XI(ni)=0.5d0
C              ENDDO
C              DO nj=1,NJT !initialising XP1
C                XP1(nj)=XQ(nj,nq)
C              ENDDO
C
C              CALL XPXE(NBJ(1,ne_bem),NKE(1,1,1,ne_bem),NPF(1,1),
C     '          NPNE(1,1,ne_bem),nr,NVJE(1,1,1,ne_bem),SE(1,1,ne_bem),
C     '          XA(1,1,ne_bem),XE,XP,ERROR,*100)
C              CALL CLOS11(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne_bem),
C     '          DISTANCE,XE,XI(1),XP1,.TRUE.,ERROR,*100)
C
C              !store xi in XIQ
C              DO ni=1,NITB
CC*** Try with this on and off!
C                IF(XI(ni).LT.0.01d0) XI(ni)=0.0d0
C                IF(XI(ni).GT.0.99d0) XI(ni)=1.0d0
C                XIQ(ni,nq)=XI(ni)
C              ENDDO
C
C              !calculate the FIXQ array for nx_upd
C              DO ni=1,NIYFIXM
C                FIXQ(nq,ni,nx_upd)=.FALSE.
C              ENDDO
C              POTENTIAL=.FALSE.
C              FLUX=.FALSE.
C
CCC*s
CC              SPECIAL=.FALSE.
CCC*f
C              DO nn=1,NNT(nb)
C                NODES(nn)=0
C                nv=1
C                nh=NH_LOC(1,nx)
C                np=NPNE(nn,nb,ne_bem)
CCC*s
CC                IF((np.EQ.6).OR.(np.EQ.8)) SPECIAL=.TRUE.
CCC*f
C                nc=1
C                ny=NYNP(1,nv,nh,np,0,nc,nr)
C                IF(ny.GT.0) THEN
C                  IF(FIX(ny,1,nx)) THEN
C                    NODES(nn)=1
C                  ELSE
C                    nc=2
C                    ny=NYNP(1,nv,nh,np,0,nc,nr)
C                    IF(ny.GT.0) NODES(nn)=2
C                  ENDIF
C                ENDIF
C              ENDDO
C              IF((NODES(1).EQ.1).AND.(NODES(2).EQ.1)) POTENTIAL=.TRUE.
C              IF((NODES(1).EQ.2).AND.(NODES(2).EQ.2)) FLUX=.TRUE.
C              IF((NODES(1).EQ.1).AND.(NODES(2).EQ.2)) THEN
C                IF(XI(1).LE.0.5d0) THEN
C                  POTENTIAL=.TRUE.
C                ELSE
C                  FLUX=.TRUE.
C                ENDIF
C              ENDIF
C              IF((NODES(1).EQ.2).AND.(NODES(2).EQ.1)) THEN
C                IF(XI(1).LE.0.5d0) THEN
C                  FLUX=.TRUE.
C                ELSE
C                  POTENTIAL=.TRUE.
C                ENDIF
C              ENDIF
C              IF(BLOOD) FIXQ(nq,3,nx_upd)=.TRUE.
C
CCC*s
CC              IF(SPECIAL) FLUX=.TRUE.
CCC*f
C
CCC new check to see if excluded node is in the current element
CC              DO nn=1,NNT(nb)
CC                np=NPNE(nn,nb,ne_bem)
CC                DO nn2=1,CPLST(0,1)
CC                  IF(np.EQ.CPLST(nn2,1)) FLUX=.TRUE.
CC                ENDDO !nn2
CC              ENDDO !nn
CC end new check
C
C              IF(POTENTIAL) THEN
C                FIXQ(nq,2,nx_upd)=.TRUE.
C              ELSEIF(FLUX) THEN
C                FIXQ(nq,1,nx_upd)=.TRUE.
C              ELSE
C                CALL ASSERT(ne_bem.GT.0,'>>No ny value found',
C     '            ERROR,*100)
C              ENDIF
C
C            ENDIF !external
C            GO TO 102
C            !this statement is designed to be skipped if no error
C            !occurs. However if a error occurs within a subroutine
C            !the alternate return jumps to line 100 to set the flag
C 100        CONTINUE
CC$OMP CRITICAL(UPELEM_1)
C            ERROR_FLAG=.TRUE.
C            WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
C     '        //'results may be unreliable'')')
C            CALL WRITES(IODI,OP_STRING,ERROR,*101)
C 101        CONTINUE
CC$OMP END CRITICAL(UPELEM_1)
C 102        CONTINUE
C          ENDIF !not error_flag
C        ENDDO
CC$OMP END PARALLEL DO
C        UP_NENQ=.FALSE.
C        CALL CPU_TIMER(CPU_USER,TIME_STOP)
C        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
CC$OMP CRITICAL(UPELEM_2)
C        WRITE(OP_STRING,'(1X,''Time to update NENQ'',F6.2''s cpu'')')
C     '    ELAPSED_TIME
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(UPELEM_2)
C
C        IF(DOP) THEN
CC$OMP CRITICAL(UPELEM_3)
C          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
C            IF(NWQ(1,nq,1).NE.0) THEN !boundary point
C              WRITE(OP_STRING,'('' nq,ne,xi'',2I8,F12.6)')
C     '          nq,NENQ(NENQ(0,nq),nq),XIQ(1,nq)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDDO
CC$OMP END CRITICAL(UPELEM_3)
C        ENDIF
C
C      ENDIF
C
C      CALL EXITS('UPELEM')
C      RETURN
C 9999 CALL ERRORS('UPELEM',ERROR)
C      CALL EXITS('UPELEM')
C      RETURN 1
C      END


