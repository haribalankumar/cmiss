      SUBROUTINE UPELEM(IBT,IDO,INP,NBJ,NEELEM,NENQ,NKJE,NPF,NPNE,
     '  NQLIST,NRLIST,NVJE,NWQ,NXLIST,NYNP,SE,XA,XE,XIQ,XP,XQ,STRING,
     '  FIX,FIXQ,ERROR,*)

C#### Subroutine: UPELEM
C###  Description:
C###    Updates element variables
C**** Written by Martin Buist 4-June-1999

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NQLIST(0:NQM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NWQ(8,0:NQM,NAM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XIQ(NIM,NQM),
     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXQ(NYQM,NIYFIXM,NXM)
!     Local variables
      INTEGER IBEG,IEND,IT,ITMAX,nb,nb_torso,nc,ne_bem,ne_t,
     '  nh,ni,NITB,nj,nn,NODES(27),noelemt,np,nr,nr_blood,nr_grid,nrr,
     '  nr_torso,nq,nqq,nv,nx,nxc,NX_TORSO_LIST(9),
     '  nx_upd,nxx,ny,N3CO
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 DISTANCE,DQ,DQMAX,XI(3),XP1(3)
      LOGICAL ALL_REGIONS,BLOOD,CBBREV,ERROR_FLAG,POTENTIAL,
     '  FLUX

      CALL ENTERS('UPELEM',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update elements
C###  Parameter:      <region #s[1,2,3]>
C###    Three regions are required, region 1 is the region in which
C###    grid points are defined. Region 2 is the ventricle region and
C###    region 3 is the torso cavity region.
C###  Parameter:      <class #s[1,2,3]>
C###    Three classes are required, class 1 is the extracellular
C###    update class, class 2 is the ventricle solution class and
C###    class 3 is the torso cavity-body surface solution class.
C###  Parameter:      <grids name[none]>
C###    This is an optional parameter which allows a subset of
C###    grid points to be updated. Only points in the grid group
C###    will be affected.
C###  Description:
C###    Generate a boundary element number and xi location
C###    for external grid points in coupled FD/BEM problems

        OP_STRING(1)=STRING(1:IEND)//' '
        OP_STRING(2)=BLANK(1:15)//'<region #s[1,2,3]>'
        OP_STRING(3)=BLANK(1:15)//'<class #s[1,2,3]>'
        OP_STRING(4)=BLANK(1:15)//'<grids name[none]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPELEM',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'REGION',3,noco+1,NTCO,n3co)) THEN
          CALL ASSERT(NRLIST(0).LE.9,' >>Maximum 9 regions',
     '      ERROR,*9999)
          nr_grid=NRLIST(1)
        ELSE
          CALL ASSERT(.FALSE.,
     '      '>>Region parameter is compulsory',ERROR,*9999)
        ENDIF
        nr_blood=NRLIST(2)
        nr_torso=NRLIST(3)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)

        IF(CBBREV(CO,'CLASS',2,noco+1,NTCO,n3co)) THEN
          CALL ASSERT(NXLIST(0).LE.9,' >>Maximum 9 classes',
     '      ERROR,*9999)

          !Grid class
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_upd.GT.0,
     '      '>>No nx defined for this solve class',ERROR,*9999)

          DO nxx=1,9
            NX_TORSO_LIST(nxx)=0
          ENDDO

          !Bem classes
          DO nxx=2,NXLIST(0)
            nxc=NXLIST(nxx)
            CALL NX_LOC(NX_INQUIRE,nxc,NX_TORSO_LIST(nxx),
     '        NX_SOLVE,ERROR,*9999)
            CALL ASSERT(NX_TORSO_LIST(nxx).GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ENDDO !nxx
        ELSE
          CALL ASSERT(.FALSE.,
     '      '>>Class parameter is compulsory',ERROR,*9999)
        ENDIF

        CALL ASSERT(NXLIST(0).EQ.NRLIST(0),
     '    '>>Must have equal #s of classes and regions',ERROR,*9999)

C        CALL ASSERT(NENQ(0,nq).LE.7,
C     '    '>>Increase NENQ allocation',ERROR,*100)
        !Initialise
        DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
          FIXQ(nq,1,nx_upd)=.FALSE.
          FIXQ(nq,2,nx_upd)=.FALSE.
          FIXQ(nq,3,nx_upd)=.FALSE.
          NENQ(0,nq)=NENQ(0,nq)+1
          NENQ(NENQ(0,nq),nq)=0
        ENDDO !nq

        IF(CBBREV(CO,'GRIDS',2,noco+1,NTCO,n3co)) THEN
          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
            FIXQ(nq,3,nx_upd)=.TRUE.
          ENDDO !nq
          CALL PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*9999)
          DO nqq=1,NQLIST(0)
            nq=NQLIST(nqq)
            FIXQ(nq,3,nx_upd)=.FALSE.
          ENDDO !nqq
        ELSE
          NQLIST(0)=0
          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
            IF(nr_blood.EQ.0) THEN
              IF(NWQ(1,nq,1).LT.nq) THEN
                NQLIST(0)=NQLIST(0)+1
                NQLIST(NQLIST(0))=nq
              ELSE
                FIXQ(nq,3,nx_upd)=.TRUE.
              ENDIF
            ELSE IF(nr_torso.EQ.0) THEN
              IF(NWQ(1,nq,1).GT.nq) THEN
                NQLIST(0)=NQLIST(0)+1
                NQLIST(NQLIST(0))=nq
              ENDIF
            ELSE
              NQLIST(0)=NQLIST(0)+1
              NQLIST(NQLIST(0))=nq
            ENDIF
          ENDDO !nq
        ENDIF !grid_group

        !init
        ERROR_FLAG=.FALSE.
        ITMAX=50
        DO nq=1,NQT
          IF(NWQ(1,nq,1).GT.0) FIXQ(nq,2,nx_upd)=.TRUE.
        ENDDO

        CALL CPU_TIMER(CPU_USER,TIME_START)

        DO nqq=1,NQLIST(0)
          nq=NQLIST(nqq)
          IF(.NOT.ERROR_FLAG) THEN

            !initialise to -ve so we can check later
            DO ni=1,NIM
              XIQ(ni,nq)=-1.0d0
            ENDDO
            IF(NWQ(1,nq,1).NE.0) THEN !boundary point

              DQMAX=RMAX
              ne_bem=0
              nr=0
              nx=0

              !check elements from the torso regions
              DO nrr=2,NRLIST(0)
                nr_torso=NRLIST(nrr)
                IF(nr_torso.GT.0) THEN
                  DO noelemt=1,NEELEM(0,nr_torso)
                    ne_t=NEELEM(noelemt,nr_torso)
                    nb_torso=NBJ(1,ne_t)
                    DQ=0.0d0
                    DO nn=1,NNT(nb_torso)
                      np=NPNE(nn,nb_torso,ne_t)
                      DO nj=1,NJT
                        DQ=DQ+(XQ(nj,nq)-XP(1,1,nj,np))**2.0d0
                      ENDDO
                    ENDDO
                    DQ=DSQRT(DQ)
                    IF(DQ.LT.DQMAX) THEN
                      DQMAX=DQ
                      ne_bem=ne_t
                      nr=nr_torso
                      nx=NX_TORSO_LIST(nrr)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              CALL ASSERT(ne_bem.GT.0,'>>No torso element found',
     '          ERROR,*100)

              BLOOD=.FALSE.
              !look for blood region
              IF(NWQ(1,nq,1).GT.nq) BLOOD=.TRUE.

              !store element in NENQ
C              CALL ASSERT(NENQ(0,nq).LE.7,
C     '          '>>Increase NENQ allocation',ERROR,*100)
C              NENQ(0,nq)=NENQ(0,nq)+1
              NENQ(NENQ(0,nq),nq)=ne_bem

              !Calculate xi position
              nb=NBJ(1,ne_bem)
              NITB=NIT(nb)
              DO ni=1,NITB !initialising XI
                XI(ni)=0.5d0
              ENDDO
              DO nj=1,NJT !initialising XP1
                XP1(nj)=XQ(nj,nq)
              ENDDO

              CALL XPXE(NBJ(1,ne_bem),NKJE(1,1,1,ne_bem),NPF(1,1),
     '          NPNE(1,1,ne_bem),nr,NVJE(1,1,1,ne_bem),SE(1,1,ne_bem),
     '          XA(1,1,ne_bem),XE,XP,ERROR,*100)
              CALL CLOS11(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne_bem),
     '          DISTANCE,XE,XI(1),XP1,.TRUE.,ERROR,*100)

              !store xi in XIQ
              DO ni=1,NITB
                IF(XI(ni).LT.0.01d0) XI(ni)=0.0d0
                IF(XI(ni).GT.0.99d0) XI(ni)=1.0d0
                XIQ(ni,nq)=XI(ni)
              ENDDO

              !calculate the FIXQ array for nx_upd
              DO ni=1,NIYFIXM
                FIXQ(nq,ni,nx_upd)=.FALSE.
              ENDDO
              POTENTIAL=.FALSE.
              FLUX=.FALSE.

              DO nn=1,NNT(nb)
                NODES(nn)=0
                nv=1
                nh=NH_LOC(1,nx)
                np=NPNE(nn,nb,ne_bem)
                nc=1
                ny=NYNP(1,nv,nh,np,0,nc,nr)
                IF(ny.GT.0) THEN
                  IF(FIX(ny,1,nx)) THEN
                    NODES(nn)=1
                  ELSE
                    nc=2
                    ny=NYNP(1,nv,nh,np,0,nc,nr)
                    IF(ny.GT.0) NODES(nn)=2
                  ENDIF
                ENDIF
              ENDDO
              IF((NODES(1).EQ.1).AND.(NODES(2).EQ.1)) POTENTIAL=.TRUE.
              IF((NODES(1).EQ.2).AND.(NODES(2).EQ.2)) FLUX=.TRUE.
              IF((NODES(1).EQ.1).AND.(NODES(2).EQ.2)) THEN
                IF(XI(1).LE.0.5d0) THEN
                  POTENTIAL=.TRUE.
                ELSE
                  FLUX=.TRUE.
                ENDIF
              ENDIF
              IF((NODES(1).EQ.2).AND.(NODES(2).EQ.1)) THEN
                IF(XI(1).LE.0.5d0) THEN
                  FLUX=.TRUE.
                ELSE
                  POTENTIAL=.TRUE.
                ENDIF
              ENDIF
              IF(BLOOD) FIXQ(nq,3,nx_upd)=.TRUE.

              IF(POTENTIAL) THEN
                FIXQ(nq,2,nx_upd)=.TRUE.
              ELSEIF(FLUX) THEN
                FIXQ(nq,1,nx_upd)=.TRUE.
              ELSE
                CALL ASSERT(ne_bem.GT.0,'>>No ny value found',
     '            ERROR,*100)
              ENDIF

            ENDIF !external
            GO TO 102
            !this statement is designed to be skipped if no error
            !occurs. However if a error occurs within a subroutine
            !the alternate return jumps to line 100 to set the flag
 100        CONTINUE
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '        //'results may be unreliable'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
 102        CONTINUE
          ENDIF !not error_flag
        ENDDO
        UP_NENQ=.FALSE.
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        WRITE(OP_STRING,'(1X,''Time to update NENQ'',F6.2,''s cpu'')')
     '    ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        IF(DOP) THEN
          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
            IF(NWQ(1,nq,1).NE.0) THEN !boundary point
              WRITE(OP_STRING,'('' nq,ne,xi'',2I8,F12.6)')
     '          nq,NENQ(NENQ(0,nq),nq),XIQ(1,nq)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF

      ENDIF

      CALL EXITS('UPELEM')
      RETURN
 9999 CALL ERRORS('UPELEM',ERROR)
      CALL EXITS('UPELEM')
      RETURN 1
      END

