      SUBROUTINE CHKSOL(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,NBH,NBHF,
     '  NBJ,NBJF,NDIPOLES,NEELEM,NFF,NHE,NHP,
     '  NKEF,NKH,NKHE,NKJE,NNF,NPF,NP_INTERFACE,NPNE,NPNF,NPNODE,
     '  NRLIST,NVHE,NVHF,NVJE,NVJF,NVHP,NW,NXLIST,NYNE,NYNP,CE,
     '  CURVCORRECT,DIPOLE_CEN,DIPOLE_DIR,PG,RG,SE,SF,WG,XA,XE,XG,XP,
     '  YP,ZA,ZE,ZP,STRING,ERROR,*)

C#### Subroutine: CHKSOL
C###  Description:
C###    CHKSOL outputs solutions and compares to analytic solution.
C###    For the relative error measure, REL_SCALE is used as a rough
C###    measure of the 'scale' of the solution (ie if voltages are
C###    approx 1kV, then errors of 1V are small.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),
     '  IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NDIPOLES(NRM,NXM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NNF(0:17,6,NBFM),NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHF(NNM,NBFM,NHM),NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),
     '  NVHP(NHM,NPM,NCM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM)
      REAL*8 CE(NMM,NEM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IBEG1,IEND,IEND1,IVAL,iy,N3CO,nb,nc,ne,nef,neff,nf,
     '  ng,nh,nhx,nk,NKHF(NKM,NNM,NHM),NKJF(NKM,NNM,NJM),nn,noelem,
     '  no_nrlist,nonode,np,npp,nr,ns,NTOTAL_ABS,NTOTAL_PER,NTOTAL_REL,
     '  NEXACT,NUM_NC,nv,nx,nxc,ny
      REAL*8 ABSOLUTE,ABSEXACT,COND,DIPOLESTRENGTH,DXIX(3,3),
     '  EXACT,INT,INTVALUE(4),
     '  INTVALUESQ(4),GL(3,3),GU(3,3),MIN_VAL,MAX_VAL,
     '  PERCENT,RAD,REL_SCALE,RMS_ERROR_ABS,
     '  RMS_ERROR_PER,RMS_ERROR_REL,RMSEXACT,RDM_ERROR,RDM_EXACT,
     '  RDM_SUM,RELATIVE,SUM,SUM1,PHI,THETA,W,XX(3),ZZ(3)
      CHARACTER FILE*(MXCH)
      LOGICAL ALL_REGIONS,ANALDIFF,CALCANAL,CBBREV,CONT,CROSSDERIV,
     '  FACE,HERMITE_Q,HERMITE_U,INTERFACE,MATCH,NORMAL,OPFILE,
     '  OUTPUT,STORE,STORE_ABSOLUTE,STORE_PERCENTAGE,STORE_RELATIVE,
     '  SPHERE,SUMMARY,RDM

      CALL ENTERS('CHKSOL',*9999)

      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------
C#### Command: FEM check solution<;FILENAME>
C###  Description:
C###    Check the solution versus an analytic solution.
C###  Parameter: <store absolute/percentage/relative>[absolute]
C###  Description:
C###    Store the absolute/percentage/relative error in YP(ny,5)
C###  Parameter: <sphere>
C###  Description:
C###    Assume circular or spherical geometry in order to calculate
C###    integral errors.
C###  Parameter: <interface>
C###  Description:
C###    List the errors for only those nodes on an interface.
C###  Parameter: <rdm>
C###  Description:
C###    Print out the rdm error
C###  Parameter: <summary>
C###  Description:
C###    Does not print out individual node comparisons,
C###     only summaries
C###  Parameter: <region #s/ALL>[1]
C###  Description:
C###    Specify the regions to check the solution for.
C###  Parameter: <class #>[1]
C###  Description:
C###    Specify the class to check the solution for.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<store absolute/percentage/relative>'
     '    //'[absolute]'
        OP_STRING(3)=BLANK(1:15)//'<sphere>'
        OP_STRING(4)=BLANK(1:15)//'<interface>'
        OP_STRING(5)=BLANK(1:15)//'<rdm>'
        OP_STRING(5)=BLANK(1:15)//'<region #s/ALL>[1]'
        OP_STRING(6)=BLANK(1:15)//'<class #>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C----------------------------------------------------------------------

C#### Command: FEM check solution<;FILENAME>
C###  Parameter: <forward/inverse [inverse]>
C###  Parameter: <infile1 FILENAME[$current]>
C###  Parameter: <infile2 FILENAME[$current_new]>
C###  Parameter: <(ascii/binary)[ascii]>
C###  Parameter: <class #>[1]
C###  Description:
C###    LKC 29-SEP-1999 THIS DOES NOT APPEAR TO BE IMPLEMENTED YET.
C###    Checks the accuracy of a forward/inverse calculation on two
C###    history files - infile1 and infile2.  A typical usage would
C###    involve applying a forward or inverse matrix to some data in a
C###    history file (infile2 say, which would contain results on both
C###    inner and outer surfaces) and then comparing the resulting
C###    timeseries on one of the surfaces (stored in infile1) to it.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<forward/inverse [inverse]>'
        OP_STRING(3)=BLANK(1:15)//' <infile1 FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//' <infile2 FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(5)=BLANK(1:15)//' <(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<class #>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C----------------------------------------------------------------------

      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','CHKSOL',ERROR,*9999)
      ELSE
        IF(NTCOQU(NOCO).GT.0) THEN
          FILE=COQU(NOCO,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opchks','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,NOCO,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '   ERROR,*9999)


C GMH 15/11/95 Adding store option to YP
        IF(CBBREV(CO,'STORE',1,NOCO+1,NTCO,N3CO)) THEN
          STORE=.TRUE.
          STORE_ABSOLUTE=.FALSE.
          STORE_PERCENTAGE=.FALSE.
          STORE_RELATIVE=.FALSE.
          IF(CBBREV(CO,'ABSOLUTE',2,NOCO+1,NTCO,N3CO)) THEN
            STORE_ABSOLUTE=.TRUE.
          ELSE IF(CBBREV(CO,'PERCENTAGE',2,NOCO+1,NTCO,N3CO)) THEN
            STORE_PERCENTAGE=.TRUE.
          ELSE IF(CBBREV(CO,'RELATIVE',2,NOCO+1,NTCO,N3CO)) THEN
            STORE_RELATIVE=.TRUE.
          ELSE
            STORE_ABSOLUTE=.FALSE.
          ENDIF
        ELSE
          STORE=.FALSE.
        ENDIF

C CPB 8/11/96 Adding sphere option
        IF(CBBREV(CO,'SPHERE',1,NOCO+1,NTCO,N3CO)) THEN
          SPHERE=.TRUE.
        ELSE
          SPHERE=.FALSE.
        ENDIF

C CPB 1/2/96 Adding interface option
        IF(CBBREV(CO,'INTERFACE',1,NOCO+1,NTCO,N3CO)) THEN
          INTERFACE=.TRUE.
        ELSE
          INTERFACE=.FALSE.
        ENDIF

C CPB 5/7/00 Adding RDM option
        IF(CBBREV(CO,'RDM',1,NOCO+1,NTCO,N3CO)) THEN
          RDM=.TRUE.
        ELSE
          RDM=.FALSE.
        ENDIF


C LKC 29-SEP-1999 The history command is not implemented yet -
C  but is used in some examples incorrectly.
        IF(CBBREV(CO,'INFILE1',7,NOCO+1,NTCO,N3CO)) THEN
          ERROR='History comparison not implemented'
          GOTO 9999
        ENDIF

C LKC 31-JUL-2000 Adding the summary option
        SUMMARY=.FALSE.
        IF(CBBREV(CO,'SUMMARY',1,NOCO+1,NTCO,N3CO)) THEN
          SUMMARY=.TRUE.
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          WRITE(OP_STRING,'(/'' Region '',I1)') nr
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(.NOT.(APPLY_TRANSFER.OR.APPLY_INVERSE)) THEN
            CALL ASSERT(ANAL_CHOICE(nr).GT.0,
     '        '>>Generate analytic b.c. expressions first',ERROR,*9999)
          ENDIF

          IF(.NOT.SUMMARY) THEN
            WRITE(OP_STRING,'(/'' Dependent variable '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          nc=1 ! dependent variable
          FORMAT='('' Node #'',3X,''Numerical'','
     '      //'4X,''Analytic'',2X,''% error'','
     '      //'1X,''Absolute error'',1X,''Relative error'')'
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          RMS_ERROR_ABS=0.0d0
          RMS_ERROR_PER=0.0d0
          RMS_ERROR_REL=0.0d0
          RDM_SUM=0.0d0
          NTOTAL_ABS=0
          NTOTAL_PER=0
          NTOTAL_REL=0
          ABSEXACT=0.0d0
          RMSEXACT=0.0d0
          RDM_EXACT=0.0d0
          NEXACT=0
          MIN_VAL=RMAX
          MAX_VAL=-RMAX
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '        NP_INTERFACE(np,1).EQ.nr)) THEN
              ny=NYNP(1,1,NH_LOC(1,nx),np,0,nc,nr)
              ABSEXACT=ABSEXACT+DABS(YP(ny,7,nx))
              RDM_EXACT=RDM_EXACT+YP(ny,7,nx)**2
              RMSEXACT=RMSEXACT+YP(ny,7,nx)**2
              IF(YP(ny,7,nx).LT.MIN_VAL) MIN_VAL=YP(ny,7,nx)
              IF(YP(ny,7,nx).GT.MAX_VAL) MAX_VAL=YP(ny,7,nx)
              NEXACT=NEXACT+1
            ENDIF
          ENDDO !nonode
          IF(NEXACT.EQ.0) THEN
            ABSEXACT=0.0d0
            RMSEXACT=-1.0d0 !Impossible amount
            REL_SCALE=1.0d0
          ELSE
            ABSEXACT=ABSEXACT/DBLE(NEXACT)
            RMSEXACT=DSQRT(RMSEXACT/DBLE(NEXACT))
            REL_SCALE=MAX_VAL-MIN_VAL
            IF(REL_SCALE.LT.ZERO_TOL) REL_SCALE=1.0d0 !use actual diff
          ENDIF
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '        NP_INTERFACE(np,1).EQ.nr)) THEN
              ny=NYNP(1,1,NH_LOC(1,nx),np,0,nc,nr)
              EXACT=YP(ny,7,nx) !Exact expressions calcu in IPANA# and
C                                stored in YP(ny,7,nx).
              ABSOLUTE=DABS(EXACT-YP(ny,1,nx))
C cpb 26/1/96 GMH over-ruled
CC GMH 9/11/95 Change to error/(1+abs)
C            PERCENT=DABS((EXACT-YP(ny,1,nx))/
C     '        (1.0D0+DABS(EXACT)))*100.0d0
              IF(DABS(EXACT).LE.0.001d0) THEN
                PERCENT=(EXACT-YP(ny,1,nx))/(1.0d0+DABS(EXACT))*100.0d0
              ELSE
                PERCENT=(EXACT-YP(ny,1,nx))/DABS(EXACT)*100.0d0
              ENDIF
              RELATIVE=ABSOLUTE/REL_SCALE
              NTOTAL_ABS=NTOTAL_ABS+1
              NTOTAL_PER=NTOTAL_PER+1
              NTOTAL_REL=NTOTAL_REL+1
              RMS_ERROR_ABS=RMS_ERROR_ABS+ABSOLUTE**2
              RMS_ERROR_PER=RMS_ERROR_PER+PERCENT**2
              RMS_ERROR_REL=RMS_ERROR_REL+RELATIVE**2
              RDM_SUM=RDM_SUM+ABSOLUTE**2

              IF(.NOT.SUMMARY) THEN
                FORMAT=
     '            '(2X,I5,1X,D11.4,1X,D11.4,1X,F8.2,1X,D14.7,1X,D14.7)'
                WRITE(OP_STRING,FORMAT) np,YP(ny,1,nx),EXACT,PERCENT,
     '            ABSOLUTE,RELATIVE
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF


C GMH 15/11/95 Adding store option to YP
              IF(STORE) THEN
                IF(STORE_ABSOLUTE) THEN
                  YP(ny,5,nx)=ABSOLUTE
                ELSE IF(STORE_PERCENTAGE) THEN
                  YP(ny,5,nx)=PERCENT
                ELSE IF(STORE_RELATIVE) THEN
                  YP(ny,5,nx)=RELATIVE
                ENDIF
              ENDIF
            ENDIF
          ENDDO !nonode (np)
          IF(NTOTAL_ABS.GT.0) THEN
            RMS_ERROR_ABS=DSQRT(RMS_ERROR_ABS/DBLE(NTOTAL_ABS))
          ELSE
            RMS_ERROR_ABS=-1.0d0 !A default value
          ENDIF
          IF(NTOTAL_PER.GT.0) THEN
            RMS_ERROR_PER=DSQRT(RMS_ERROR_PER/DBLE(NTOTAL_PER))
          ELSE
            RMS_ERROR_PER=-1.0d0 !A default value
          ENDIF
          IF(NTOTAL_REL.GT.0) THEN
            RMS_ERROR_REL=DSQRT(RMS_ERROR_REL/DBLE(NTOTAL_REL))
          ELSE
            RMS_ERROR_REL=-1.0d0 !A default value
          ENDIF
          FORMAT='('' RMS error (Percent)  = '',D11.4,'' ('','
     '      //'F8.2,''%)'')'
          WRITE(OP_STRING,FORMAT) RMS_ERROR_PER,RMS_ERROR_PER
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='('' RMS error (Absolute) = '',D11.4)'
          WRITE(OP_STRING,FORMAT) RMS_ERROR_ABS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='('' RMS error (Relative) = '',D11.4,'' ('','
     '      //'F8.2,''%)'')'
          WRITE(OP_STRING,FORMAT) RMS_ERROR_REL,RMS_ERROR_REL*100.0d0
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C CPB 5/7/00 Adding RDM
          IF(RDM) THEN
            IF(DABS(RDM_EXACT).GT.ZERO_TOL) THEN
              RDM_ERROR=DSQRT(RDM_SUM/RDM_EXACT)
            ELSE
              RDM_ERROR=-1.0d0
            ENDIF
            FORMAT='('' RDM error            = '',D11.4,'' ('','
     '        //'F8.2,''%)'')'
            WRITE(OP_STRING,FORMAT) RDM_ERROR,RDM_ERROR*100.0d0
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
C cpb 2/2/96 This is not the best way to find hermite_u!!!!
          HERMITE_U=.FALSE.
          nb=NBASEF(NBH(NH_LOC(1,nx),2,NEELEM(1,nr)),1) !Basis fn for q
          IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
            IF(IBT(1,1,nb).EQ.2) THEN !Hermite
              HERMITE_U=.TRUE.
            ENDIF
          ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
            IF(IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2) THEN
              HERMITE_U=.TRUE.
            ELSE IF(IBT(1,1,nb).EQ.1.AND.IBT(1,2,nb).EQ.2) THEN
              HERMITE_U=.TRUE.
            ELSE IF(IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.1) THEN
              HERMITE_U=.TRUE.
            ELSE IF(IBT(1,1,nb).EQ.3) THEN !Simplex
              IF(IBT(2,1,nb).EQ.4) THEN !Special hermite simplex
                HERMITE_U=.TRUE.
              ENDIF
            ELSE IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6.OR.
     '          IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6) THEN !Sector
              IF(IBT(2,1,nb).EQ.4.OR.IBT(2,2,nb).EQ.4) THEN
                HERMITE_U=.TRUE.
              ENDIF
            ENDIF
          ENDIF
C          NORMAL=(CALL_MESH.AND.JTYP14.EQ.2.AND.NJT.EQ.2).OR.HERMITE_U
          NORMAL=.TRUE.
          IF(NORMAL) THEN
            IF(.NOT.(APPLY_TRANSFER.OR.APPLY_INVERSE)) THEN
C***          Don't need to look at the normal derivative if transfer
C***          matrices used.
              RMS_ERROR_ABS=0.0d0
              RMS_ERROR_PER=0.0d0
              RMS_ERROR_REL=0.0d0
              RDM_SUM=0.0d0
              NTOTAL_ABS=0
              NTOTAL_PER=0
              NTOTAL_REL=0
              IF(ITYP4(nr,nx).EQ.2) THEN !BE region
                WRITE(OP_STRING,'(/'' Normal Derivative '')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                nc=2 ! normal derivative for bem region
                FORMAT='('' Node #'',3X,''Numerical'','
     '            //'4X,''Analytic'',2X,''% error'','
     '            //'1X,''Absolute error'',1X,''Relative error'')'
                WRITE(OP_STRING,FORMAT)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ABSEXACT=0.0d0
                RMSEXACT=0.0d0
                RDM_EXACT=0.0d0
                NEXACT=0
                MIN_VAL=RMAX
                MAX_VAL=-RMAX
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '              NP_INTERFACE(np,1).EQ.nr)) THEN
C cpb 1/2/96 Output normal derivative at all nodes
C                    IF(((NJ_LOC(NJL_GEOM,0,nr).EQ.2).AND.
C     '                (NKH(NH_LOC(1,nx),np,nc)-KTYP93(nc,nr).GE.
C     '                2)).OR.(NJ_LOC(NJL_GEOM,0,nr).EQ.2.AND.
C     '              NKJ(1,np).GT.2).OR.
C     '                ((NJ_LOC(NJL_GEOM,0,nr).EQ.3).AND.
C     '                (NKH(NH_LOC(1,nx),np,nc)-KTYP93(nc,nr).GE.3)))
C     '                THEN
                    ny=NYNP(1,1,NH_LOC(1,nx),np,0,2,nr)
                    ABSEXACT=ABSEXACT+DABS(YP(ny,7,nx))
                    RMSEXACT=RMSEXACT+YP(ny,7,nx)**2
                    RDM_EXACT=RDM_EXACT+YP(ny,7,nx)**2
                    IF(YP(ny,7,nx).LT.MIN_VAL) MIN_VAL=YP(ny,7,nx)
                    IF(YP(ny,7,nx).GT.MAX_VAL) MAX_VAL=YP(ny,7,nx)
                    NEXACT=NEXACT+1
C                    ENDIF
                  ENDIF
                ENDDO !nonode
                IF(NEXACT.EQ.0) THEN
                  ABSEXACT=0.0d0
                  RMSEXACT=-1.0d0 !Impossible amount
                  REL_SCALE=1.0d0
                ELSE
                  ABSEXACT=ABSEXACT/DBLE(NEXACT)
                  RMSEXACT=DSQRT(RMSEXACT/DBLE(NEXACT))
                  REL_SCALE=MAX_VAL-MIN_VAL
                  IF(REL_SCALE.LT.ZERO_TOL) REL_SCALE=1.0d0
                ENDIF
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '              NP_INTERFACE(np,1).EQ.nr)) THEN
C cpb 1/2/96 Output normal derivative at all nodes
C                    OUTPUT=.FALSE.
C                    IF((NJ_LOC(NJL_GEOM,0,nr).EQ.2).AND.
C     '                (NKH(NH_LOC(1,nx),np,nc)-KTYP93(nc,nr).GE.2)) THEN
C                      OUTPUT=.TRUE.
C                    ELSEIF((NJ_LOC(NJL_GEOM,0,nr).EQ.2).AND.
C     '                  (NKJ(1,np)-KTYP93(nc,nr).GE.2)) THEN
C                      OUTPUT=.TRUE.
C     ELSE IF((NJ_LOC(NJL_GEOM,0,nr).EQ.3).AND.(NKH(NH_LOC(1,nx),np,nc)
C                      -KTYP93(nc,nr).GE.3)) THEN
C                      OUTPUT=.TRUE.
C                    ELSE IF(CALL_MESH.AND.JTYP14.EQ.2.AND.NJT.EQ.2) THEN
C                      OUTPUT=.TRUE.
C                    ENDIF
                    OUTPUT=.TRUE.
                    IF(OUTPUT) THEN
                      ny=NYNP(1,1,NH_LOC(1,nx),np,0,nc,nr)
                      EXACT=YP(ny,7,nx)
C cpb 26/1/96 Changing the way %'s are calculated
C                    PERCENT=DABS((EXACT-YP(ny,1,nx))/
C     '                (1.0D0+DABS(EXACT)))*100.0d0
                      IF(DABS(EXACT).LE.0.001d0) THEN
                        PERCENT=(EXACT-YP(ny,1,nx))/(1.0d0+DABS(EXACT))*
     '                    100.0d0
                      ELSE
                        PERCENT=(EXACT-YP(ny,1,nx))/DABS(EXACT)*100.0d0
                      ENDIF
                      ABSOLUTE=DABS(EXACT-YP(ny,1,nx))
                      RELATIVE=ABSOLUTE/REL_SCALE
                      NTOTAL_ABS=NTOTAL_ABS+1
                      NTOTAL_PER=NTOTAL_PER+1
                      NTOTAL_REL=NTOTAL_REL+1
                      RMS_ERROR_ABS=RMS_ERROR_ABS+ABSOLUTE**2
                      RMS_ERROR_PER=RMS_ERROR_PER+PERCENT**2
                      RMS_ERROR_REL=RMS_ERROR_REL+RELATIVE**2
                      RDM_SUM=RDM_SUM+ABSOLUTE**2
                      FORMAT='(2X,I5,1X,D11.4,1X,D11.4,1X,F8.2,1X,'
     '                  //'D14.7,1X,D14.7)'
                      WRITE(OP_STRING,FORMAT) np,YP(ny,1,nx),EXACT,
     '                  PERCENT,ABSOLUTE,RELATIVE
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C GMH 15/11/95 Adding store option to YP
                      IF(STORE) THEN
                        IF(STORE_ABSOLUTE) THEN
                          YP(ny,5,nx)=ABSOLUTE
                        ELSE IF(STORE_PERCENTAGE) THEN
                          YP(ny,5,nx)=PERCENT
                        ELSE IF(STORE_RELATIVE) THEN
                          YP(ny,5,nx)=RELATIVE
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nonode (np)
                IF(NTOTAL_ABS.GT.0) THEN
                  RMS_ERROR_ABS=DSQRT(RMS_ERROR_ABS/DBLE(NTOTAL_ABS))
                ELSE
                  RMS_ERROR_ABS=-1.0d0 !A default value
                ENDIF
                IF(NTOTAL_PER.GT.0) THEN
                  RMS_ERROR_PER=DSQRT(RMS_ERROR_PER/DBLE(NTOTAL_PER))
                ELSE
                  RMS_ERROR_PER=-1.0d0 !A default value
                ENDIF
                IF(NTOTAL_REL.GT.0) THEN
                  RMS_ERROR_REL=DSQRT(RMS_ERROR_REL/DBLE(NTOTAL_REL))
                ELSE
                  RMS_ERROR_REL=-1.0d0 !A default value
                ENDIF
                FORMAT='('' RMS error (Percent)  = '',D11.4,'' ('','
     '            //'F8.2,''%)'')'
                WRITE(OP_STRING,FORMAT) RMS_ERROR_PER,RMS_ERROR_PER
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                FORMAT='('' RMS error (Absolute) = '',D11.4)'
                WRITE(OP_STRING,FORMAT) RMS_ERROR_ABS
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                FORMAT='('' RMS error (Relative) = '',D11.4,'' ('','
     '            //'F8.2,''%)'')'
                WRITE(OP_STRING,FORMAT) RMS_ERROR_REL,RMS_ERROR_REL*
     '            100.0d0
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C CPB 5/7/00 Adding RDM
                IF(RDM) THEN
                  IF(DABS(RDM_EXACT).GT.ZERO_TOL) THEN
                    RDM_ERROR=DSQRT(RDM_SUM/RDM_EXACT)
                  ELSE
                    RDM_ERROR=-1.0d0
                  ENDIF
                  FORMAT='('' RDM error            = '',D11.4,'' ('','
     '              //'F8.2,''%)'')'
                  WRITE(OP_STRING,FORMAT) RDM_ERROR,RDM_ERROR*100.0d0
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF !Not transfer matrices
            ENDIF !Only check normal derivative for BE routines
          ENDIF !normal
          IF(HERMITE_U) THEN
            RMS_ERROR_ABS=0.0d0
            RMS_ERROR_PER=0.0d0
            RMS_ERROR_REL=0.0d0
            RDM_SUM=0.0d0
            NTOTAL_ABS=0
            NTOTAL_PER=0
            NTOTAL_REL=0
            WRITE(OP_STRING,'(/'' Arc Length Derivative '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            nc=1 ! dependent variable
            FORMAT='('' Node #'',3X,''Numerical'','
     '        //'4X,''Analytic'',2X,''% error'','
     '        //'1X,''Absolute error'',1X,''Relative error'')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ABSEXACT=0.0d0
            RMSEXACT=0.0d0
            RDM_EXACT=0.0d0
            NEXACT=0
            MIN_VAL=RMAX
            MAX_VAL=-RMAX
            CROSSDERIV=.FALSE.
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '          NP_INTERFACE(np,1).EQ.nr)) THEN
                IF(((NJ_LOC(NJL_GEOM,0,nr).EQ.2).AND.
     '            (NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr).GE.
     '            2)).OR.((NJ_LOC(NJL_GEOM,0,nr).EQ.3).AND.
     '            (NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr).GE.3))) THEN
                  DO nk=2,NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr)
                    IF(nk.NE.4) THEN
                      ny=NYNP(nk,1,NH_LOC(1,nx),np,0,1,nr)
                      ABSEXACT=ABSEXACT+DABS(YP(ny,7,nx))
                      RMSEXACT=RMSEXACT+YP(ny,7,nx)**2
                      RDM_EXACT=RDM_EXACT+YP(ny,7,nx)**2
                      IF(YP(ny,7,nx).LT.MIN_VAL) MIN_VAL=YP(ny,7,nx)
                      IF(YP(ny,7,nx).GT.MAX_VAL) MAX_VAL=YP(ny,7,nx)
                      NEXACT=NEXACT+1
                    ELSE
                      CROSSDERIV=.TRUE.
                    ENDIF
                  ENDDO !nk
                ENDIF
              ENDIF
            ENDDO !nonode
            IF(NEXACT.EQ.0) THEN
              ABSEXACT=0.0d0
              RMSEXACT=-1.0d0 !Impossible amount
              REL_SCALE=1.0d0
            ELSE
              ABSEXACT=ABSEXACT/DBLE(NEXACT)
              RMSEXACT=DSQRT(RMSEXACT/DBLE(NEXACT))
              REL_SCALE=MAX_VAL-MIN_VAL
              IF(REL_SCALE.LT.ZERO_TOL) REL_SCALE=1.0d0 !use actual diff
            ENDIF
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '          NP_INTERFACE(np,1).EQ.nr)) THEN
                OUTPUT=.FALSE.
C               Only output results in situations where tangents
C               can be found.
                IF((NJ_LOC(NJL_GEOM,0,nr).EQ.2).AND.
     '            (NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr).GE.2)) THEN
                  OUTPUT=.TRUE.
                ELSE IF((NJ_LOC(NJL_GEOM,0,nr).EQ.3).AND.
     '              (NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr).GE.3)) THEN
                  OUTPUT=.TRUE.
                ENDIF
                IF(OUTPUT) THEN
                  DO nk=2,NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr)
                    IF(nk.NE.4) THEN
                      ny=NYNP(nk,1,NH_LOC(1,nx),np,0,1,nr)
                      EXACT=YP(ny,7,nx)
                      ABSOLUTE=DABS(EXACT-YP(ny,1,nx))
C cpb 26/1/96 Changing the way %'s are calculated
C                    PERCENT=DABS((EXACT-YP(ny,1,nx))/
C     '                (1.0D0+DABS(EXACT)))*100.0d0
                      IF(DABS(EXACT).LE.0.001d0) THEN
                        PERCENT=(EXACT-YP(ny,1,nx))/(1.0d0+DABS(EXACT))*
     '                    100.0d0
                      ELSE
                        PERCENT=(EXACT-YP(ny,1,nx))/DABS(EXACT)*100.0d0
                      ENDIF
                      RELATIVE=ABSOLUTE/REL_SCALE
                      NTOTAL_ABS=NTOTAL_ABS+1
                      NTOTAL_PER=NTOTAL_PER+1
                      NTOTAL_REL=NTOTAL_REL+1
                      RMS_ERROR_ABS=RMS_ERROR_ABS+ABSOLUTE**2
                      RMS_ERROR_PER=RMS_ERROR_PER+PERCENT**2
                      RMS_ERROR_REL=RMS_ERROR_REL+RELATIVE**2
                      RDM_SUM=RDM_SUM+ABSOLUTE**2
                      FORMAT='(2X,I5,1X,D11.4,1X,D11.4,1X,F8.2,1X,'
     '                  //'D14.7,1X,D14.7)'
                      WRITE(OP_STRING,FORMAT) np,YP(ny,1,nx),EXACT,
     '                  PERCENT,ABSOLUTE,RELATIVE
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C GMH 15/11/95 Adding store option to YP
                      IF(STORE) THEN
                        IF(STORE_ABSOLUTE) THEN
                          YP(ny,5,nx)=ABSOLUTE
                        ELSE IF(STORE_PERCENTAGE) THEN
                          YP(ny,5,nx)=PERCENT
                        ELSE IF(STORE_RELATIVE) THEN
                          YP(ny,5,nx)=RELATIVE
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO ! nk
                ENDIF
              ENDIF
            ENDDO !End np
            IF(NTOTAL_ABS.GT.0) THEN
              RMS_ERROR_ABS=DSQRT(RMS_ERROR_ABS/DBLE(NTOTAL_ABS))
            ELSE
              RMS_ERROR_ABS=-1.0d0 !A default value
            ENDIF
            IF(NTOTAL_PER.GT.0) THEN
              RMS_ERROR_PER=DSQRT(RMS_ERROR_PER/DBLE(NTOTAL_PER))
            ELSE
              RMS_ERROR_PER=-1.0d0 !A default value
            ENDIF
            IF(NTOTAL_REL.GT.0) THEN
              RMS_ERROR_REL=DSQRT(RMS_ERROR_REL/DBLE(NTOTAL_REL))
            ELSE
              RMS_ERROR_REL=-1.0d0 !A default value
            ENDIF
            FORMAT='('' RMS error (Percent)  = '',D11.4,'' ('','
     '        //'F8.2,''%)'')'
            WRITE(OP_STRING,FORMAT) RMS_ERROR_PER,RMS_ERROR_PER
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='('' RMS error (Absolute) = '',D11.4)'
            WRITE(OP_STRING,FORMAT) RMS_ERROR_ABS
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='('' RMS error (Relative) = '',D11.4,'' ('','
     '        //'F8.2,''%)'')'
            WRITE(OP_STRING,FORMAT) RMS_ERROR_REL,RMS_ERROR_REL*100.0d0
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C CPB 5/7/00 Adding RDM
            IF(RDM) THEN
              IF(DABS(RDM_EXACT).GT.ZERO_TOL) THEN
                RDM_ERROR=DSQRT(RDM_SUM/RDM_EXACT)
              ELSE
                RDM_ERROR=-1.0d0
              ENDIF
              FORMAT='('' RDM error            = '',D11.4,'' ('','
     '          //'F8.2,''%)'')'
              WRITE(OP_STRING,FORMAT) RDM_ERROR,RDM_ERROR*100.0d0
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(CROSSDERIV) THEN
              RMS_ERROR_ABS=0.0d0
              RMS_ERROR_PER=0.0d0
              RMS_ERROR_REL=0.0d0
              RDM_SUM=0.0d0
              NTOTAL_ABS=0
              NTOTAL_PER=0
              NTOTAL_REL=0
              WRITE(OP_STRING,'(/'' Arc Length Cross Derivative '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              nc=1 ! dependent variable
              FORMAT='('' Node #'',3X,''Numerical'','
     '          //'4X,''Analytic'',2X,''% error'','
     '          //'1X,''Absolute error'',1X,''Relative error'')'
              WRITE(OP_STRING,FORMAT)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ABSEXACT=0.0d0
              RMSEXACT=0.0d0
              RDM_EXACT=0.0d0
              NEXACT=0
              MIN_VAL=RMAX
              MAX_VAL=-RMAX
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '            NP_INTERFACE(np,1).EQ.nr)) THEN
                  IF((NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr)).EQ.4) THEN
                    ny=NYNP(4,1,NH_LOC(1,nx),np,0,1,nr)
                    ABSEXACT=ABSEXACT+DABS(YP(ny,7,nx))
                    RMSEXACT=RMSEXACT+YP(ny,7,nx)**2
                    RDM_EXACT=RDM_EXACT+YP(ny,7,nx)**2
                    IF(YP(ny,7,nx).LT.MIN_VAL) MIN_VAL=YP(ny,7,nx)
                    IF(YP(ny,7,nx).GT.MAX_VAL) MAX_VAL=YP(ny,7,nx)
                    NEXACT=NEXACT+1
                  ENDIF
                ENDIF
              ENDDO !nonode
              IF(NEXACT.EQ.0) THEN
                ABSEXACT=0.0d0
                RMSEXACT=-1.0d0 !Impossible amount
                REL_SCALE=1.0d0
              ELSE
                ABSEXACT=ABSEXACT/DBLE(NEXACT)
                RMSEXACT=DSQRT(RMSEXACT/DBLE(NEXACT))
                REL_SCALE=MAX_VAL-MIN_VAL
                IF(REL_SCALE.LT.ZERO_TOL) REL_SCALE=1.0d0
              ENDIF
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '            NP_INTERFACE(np,1).EQ.nr)) THEN
                  IF((NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr)).EQ.4) THEN
                    ny=NYNP(4,1,NH_LOC(1,nx),np,0,1,nr)
                    EXACT=YP(ny,7,nx)
                    ABSOLUTE=DABS(EXACT-YP(ny,1,nx))
                    IF(DABS(EXACT).LE.0.001d0) THEN
                      PERCENT=(EXACT-YP(ny,1,nx))/(1.0d0+DABS(EXACT))*
     '                  100.0d0
                    ELSE
                      PERCENT=(EXACT-YP(ny,1,nx))/DABS(EXACT)*100.0d0
                    ENDIF
                    RELATIVE=ABSOLUTE/REL_SCALE
                    NTOTAL_ABS=NTOTAL_ABS+1
                    NTOTAL_PER=NTOTAL_PER+1
                    NTOTAL_REL=NTOTAL_REL+1
                    RMS_ERROR_ABS=RMS_ERROR_ABS+ABSOLUTE**2
                    RMS_ERROR_PER=RMS_ERROR_PER+PERCENT**2
                    RMS_ERROR_REL=RMS_ERROR_REL+RELATIVE**2
                    RDM_SUM=RDM_SUM+ABSOLUTE**2
                    FORMAT='(2X,I5,1X,D11.4,1X,D11.4,1X,F8.2,1X,'
     '                //'D14.7,1X,D14.7)'
                    WRITE(OP_STRING,FORMAT) np,YP(ny,1,nx),EXACT,
     '                PERCENT,ABSOLUTE,RELATIVE
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C GMH 15/11/95 Adding store option to YP
                    IF(STORE) THEN
                      IF(STORE_ABSOLUTE) THEN
                        YP(ny,5,nx)=ABSOLUTE
                      ELSE IF(STORE_PERCENTAGE) THEN
                        YP(ny,5,nx)=PERCENT
                      ELSE IF(STORE_RELATIVE) THEN
                        YP(ny,5,nx)=RELATIVE
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !End np
              IF(NTOTAL_ABS.GT.0) THEN
                RMS_ERROR_ABS=DSQRT(RMS_ERROR_ABS/DBLE(NTOTAL_ABS))
              ELSE
                RMS_ERROR_ABS=-1.0d0 !A default value
              ENDIF
              IF(NTOTAL_PER.GT.0) THEN
                RMS_ERROR_PER=DSQRT(RMS_ERROR_PER/DBLE(NTOTAL_PER))
              ELSE
                RMS_ERROR_PER=-1.0d0 !A default value
              ENDIF
              IF(NTOTAL_REL.GT.0) THEN
                RMS_ERROR_REL=DSQRT(RMS_ERROR_REL/DBLE(NTOTAL_REL))
              ELSE
                RMS_ERROR_REL=-1.0d0 !A default value
              ENDIF
              FORMAT='('' RMS error (Percent)  = '',D11.4,'' ('','
     '          //'F8.2,''%)'')'
              WRITE(OP_STRING,FORMAT) RMS_ERROR_PER,RMS_ERROR_PER
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              FORMAT='('' RMS error (Absolute) = '',D11.4)'
              WRITE(OP_STRING,FORMAT) RMS_ERROR_ABS
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              FORMAT='('' RMS error (Relative) = '',D11.4,'' ('','
     '          //'F8.2,''%)'')'
              WRITE(OP_STRING,FORMAT) RMS_ERROR_REL,RMS_ERROR_REL*
     '          100.0d0
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C CPB 5/7/00 Adding RDM
              IF(RDM) THEN
                IF(DABS(RDM_EXACT).GT.ZERO_TOL) THEN
                  RDM_ERROR=DSQRT(RDM_SUM/RDM_EXACT)
                ELSE
                  RDM_ERROR=-1.0d0
                ENDIF
                FORMAT='('' RDM error            = '',D11.4,'' ('','
     '            //'F8.2,''%)'')'
                WRITE(OP_STRING,FORMAT) RDM_ERROR,RDM_ERROR*100.0d0
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            IF(.NOT.(APPLY_TRANSFER.OR.APPLY_INVERSE)) THEN
C***          Don't need to look at the normal derivative if transfer
C***          matrices used.
              RMS_ERROR_ABS=0.0d0
              RMS_ERROR_PER=0.0d0
              RMS_ERROR_REL=0.0d0
              RDM_SUM=0.0d0
              NTOTAL_ABS=0
              NTOTAL_PER=0
              NTOTAL_REL=0
              HERMITE_Q=.FALSE.
              nb=NBASEF(NBH(NH_LOC(1,nx),2,NEELEM(1,nr)),1) !Basis fn for q
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
                IF(IBT(1,1,nb).EQ.2) THEN !Hermite
                  HERMITE_Q=.TRUE.
                ENDIF
              ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
                IF(IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2) THEN
                  HERMITE_Q=.TRUE.
                ELSE IF(IBT(1,1,nb).EQ.1.AND.IBT(1,2,nb).EQ.2) THEN
                  HERMITE_Q=.TRUE.
                ELSE IF(IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.1) THEN
                  HERMITE_Q=.TRUE.
                ELSE IF(IBT(1,1,nb).EQ.3) THEN !Simplex
                  IF(IBT(2,1,nb).EQ.4) THEN !Special hermite simplex
                    HERMITE_Q=.TRUE.
                  ENDIF
                ELSE IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6.OR.
     '              IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6) THEN !Sector
                  IF(IBT(2,1,nb).EQ.4.OR.IBT(2,2,nb).EQ.4) THEN
                    HERMITE_Q=.TRUE.
                  ENDIF
                ENDIF
              ENDIF
              IF(ITYP4(nr,nx).EQ.2.AND.HERMITE_Q) THEN !BE region
C*** Assumes that the domain is a circle in 2d or a sphere 3d.
C*** NOTE: For a region between 2 spheres or circles only the magnitude
C*** of the solution can be given.  It is not known if a node lies on
C*** the inner or outer circle/sphere and so the analytic solution
C*** (which requires this information) is determined!only up to the
C*** sign. For this reason only magnitudes are compared for a single
C*** region problem.
                WRITE(OP_STRING,'(/'' Arc length derivatives of'
     '            //' the normal derivative '')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                nc=2 ! normal derivative for BEM problems
                IF(.NOT.BEMCURVATURECORRECTION) THEN
                  WRITE(OP_STRING,'('' >>Warning: VALID FOR CIRCULAR '
     '              //'OR SPHERICAL DOMAINS ONLY'')')
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
C                IF(NRT.EQ.1) THEN
C                  WRITE(OP_STRING,'('' MAGNITUDE ONLY  '')')
C                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                ENDIF
                FORMAT='('' Node #'',3X,''Numerical'','
     '            //'4X,''Analytic'',2X,''% error'','
     '            //'1X,''Absolute error'',1X,''Relative error'')'
                WRITE(OP_STRING,FORMAT)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ABSEXACT=0.0d0
                RMSEXACT=0.0d0
                RDM_EXACT=0.0d0
                NEXACT=0
                MIN_VAL=RMAX
                MAX_VAL=-RMAX
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '              NP_INTERFACE(np,1).EQ.nr)) THEN
                    IF(((NJ_LOC(NJL_GEOM,0,nr).EQ.2).AND.
     '                (NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr).GE.
     '                2)).OR.((NJ_LOC(NJL_GEOM,0,nr).EQ.3).AND.
     '                (NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr).GE.3))) THEN
                      DO nk=2,NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr)
                        ny=NYNP(nk,1,NH_LOC(1,nx),np,0,2,nr)
                        ABSEXACT=ABSEXACT+DABS(YP(ny,7,nx))
                        RMSEXACT=RMSEXACT+YP(ny,7,nx)**2
                        RDM_EXACT=RDM_EXACT+YP(ny,7,nx)**2
                        IF(YP(ny,7,nx).LT.MIN_VAL) MIN_VAL=YP(ny,7,nx)
                        IF(YP(ny,7,nx).GT.MAX_VAL) MAX_VAL=YP(ny,7,nx)
                        NEXACT=NEXACT+1
                      ENDDO !nk
                    ENDIF
                  ENDIF
                ENDDO !nonode (np)
                IF(NEXACT.EQ.0) THEN
                  ABSEXACT=0.0d0
                  RMSEXACT=-1.0d0 !Impossible value
                  REL_SCALE=1.0d0
                ELSE
                  ABSEXACT=ABSEXACT/DBLE(NEXACT)
                  RMSEXACT=DSQRT(RMSEXACT/DBLE(NEXACT))
                  REL_SCALE=MAX_VAL-MIN_VAL
                  IF(REL_SCALE.LT.ZERO_TOL) REL_SCALE=1.0d0
                ENDIF
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  IF((.NOT.INTERFACE).OR.(INTERFACE.AND.
     '              NP_INTERFACE(np,1).EQ.nr)) THEN
C                   Only output results in situations where tangents
C                   can be found.
                    OUTPUT=.FALSE.
                    IF((NJ_LOC(NJL_GEOM,0,nr).EQ.2).AND.
     '                (NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr).GE.2)) THEN
                      OUTPUT=.TRUE.
                    ELSE IF((NJ_LOC(NJL_GEOM,0,nr).EQ.3).AND.
     '                  (NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr).GE.3)) THEN
                      OUTPUT=.TRUE.
                    ENDIF
                    IF(OUTPUT) THEN
                      DO nk=2,NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr)
                        ny=NYNP(nk,1,NH_LOC(1,nx),np,0,2,nr)
                        EXACT=YP(ny,7,nx)
C LC 26/2/97 archiving section :
C   Signs are now adjusted for multiple regions
C               except line :
C
                        ABSOLUTE=DABS(EXACT-YP(ny,1,nx))

                        IF(DABS(EXACT).LE.0.001d0) THEN
                          PERCENT=(EXACT-YP(ny,1,nx))/
     '                      (1.0d0+DABS(EXACT))*100.0d0
                        ELSE
                          PERCENT=(EXACT-YP(ny,1,nx))/
     '                      DABS(EXACT)*100.0d0
                        ENDIF
C                      ENDIF
                        RELATIVE=ABSOLUTE/REL_SCALE
                        NTOTAL_ABS=NTOTAL_ABS+1
                        NTOTAL_PER=NTOTAL_PER+1
                        NTOTAL_REL=NTOTAL_REL+1
                        RMS_ERROR_ABS=RMS_ERROR_ABS+ABSOLUTE**2
                        RMS_ERROR_PER=RMS_ERROR_PER+PERCENT**2
                        RMS_ERROR_REL=RMS_ERROR_REL+RELATIVE**2
                        RDM_SUM=RDM_SUM+ABSOLUTE**2
                        FORMAT='(2X,I5,1X,D11.4,1X,D11.4,1X,F8.2,'
     '                    //'1X,D14.7,1X,D14.7)'
                        WRITE(OP_STRING,FORMAT) np,YP(ny,1,nx),EXACT,
     '                    PERCENT,ABSOLUTE,RELATIVE
                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C GMH 15/11/95 Adding store option to YP
                        IF(STORE) THEN
                          IF(STORE_ABSOLUTE) THEN
                            YP(ny,5,nx)=ABSOLUTE
                          ELSE IF(STORE_PERCENTAGE) THEN
                            YP(ny,5,nx)=PERCENT
                          ELSE IF(STORE_RELATIVE) THEN
                            YP(ny,5,nx)=RELATIVE
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDIF
                ENDDO !End np
                IF(NTOTAL_ABS.GT.0) THEN
                  RMS_ERROR_ABS=DSQRT(RMS_ERROR_ABS/DBLE(NTOTAL_ABS))
                ELSE
                  RMS_ERROR_ABS=-1.0d0 !A default value
                ENDIF
                IF(NTOTAL_PER.GT.0) THEN
                  RMS_ERROR_PER=DSQRT(RMS_ERROR_PER/DBLE(NTOTAL_PER))
                ELSE
                  RMS_ERROR_PER=-1.0d0 !A default value
                ENDIF
                IF(NTOTAL_REL.GT.0) THEN
                  RMS_ERROR_REL=DSQRT(RMS_ERROR_REL/DBLE(NTOTAL_REL))
                ELSE
                  RMS_ERROR_REL=-1.0d0 !A default value
                ENDIF
                FORMAT='('' RMS error (Percent)  = '',D11.4,'' ('','
     '            //'F8.2,''%)'')'
                WRITE(OP_STRING,FORMAT) RMS_ERROR_PER,RMS_ERROR_PER
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                FORMAT='('' RMS error (Absolute) = '',D11.4)'
                WRITE(OP_STRING,FORMAT) RMS_ERROR_ABS
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                FORMAT='('' RMS error (Relative) = '',D11.4,'' ('','
     '            //'F8.2,''%)'')'
                WRITE(OP_STRING,FORMAT) RMS_ERROR_REL,RMS_ERROR_REL*
     '            100.0d0
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C CPB 5/7/00 Adding RDM
                IF(RDM) THEN
                  IF(DABS(RDM_EXACT).GT.ZERO_TOL) THEN
                    RDM_ERROR=DSQRT(RDM_SUM/RDM_EXACT)
                  ELSE
                    RDM_ERROR=-1.0d0
                  ENDIF
                  FORMAT='('' RDM error            = '',D11.4,'' ('','
     '              //'F8.2,''%)'')'
                  WRITE(OP_STRING,FORMAT) RDM_ERROR,RDM_ERROR*100.0d0
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF !Only check are length derivatives of
                    !normal derivative for BE routines
            ENDIF !Not transfer matrices
          ENDIF !End of HERMITE_U loop

C cpb 11/1/96 Adding integral error estimate
          IF(ITYP4(nr,nx).EQ.2.AND.
     '      .NOT.(APPLY_TRANSFER.OR.APPLY_INVERSE)) THEN
            NUM_NC=2
          ELSE
            NUM_NC=1
          ENDIF
          DO nc=1,NUM_NC
            IF(nc.EQ.1) THEN
              WRITE(OP_STRING,'(/,'' Dependent variable integral '
     '          //'error:'')')
            ELSE IF(nc.EQ.2) THEN
              WRITE(OP_STRING,'(/,'' Flux integral error:'')')
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='(10X,''Numerical'',4X,''Analytic'',2X,''% error'','
     '        //'1X,''Absolute error'',1X,''Relative error'')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C cpb 3/11/96 Adding analytic integral error for a centric dipole
C for a single circle/sphere
C cpb 8/11/96 Adding analytic integral error for x^2+2xy-y^2 on a
C circle
            IF(SPHERE) THEN
              IF((NJT.EQ.2.AND.ANAL_CHOICE(nr).EQ.5).OR.
     '          (NJT.EQ.3.AND.ANAL_CHOICE(nr).EQ.12)) THEN
                INTVALUE(2)=0.0d0
                COND=CE(1,NEELEM(1,nr),nx)
                np=NPNODE(1,nr)
                IF(nc.EQ.1) THEN
                  IF(NJT.EQ.2) THEN
                    RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2)
                    DIPOLESTRENGTH=ANAL_A**2+ANAL_B**2
                    IF(DABS(COND).LE.ZERO_TOL) THEN
                      INTVALUESQ(2)=DIPOLESTRENGTH/(PI*RAD)
                    ELSE
                      INTVALUESQ(2)=DIPOLESTRENGTH/(PI*COND*COND*RAD)
                    ENDIF
                  ELSE IF(NJT.EQ.3) THEN
                    RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2+
     '                XP(1,1,3,np)**2)
                    DIPOLESTRENGTH=ANAL_A**2+ANAL_B**2+ANAL_C**2
                    IF(DABS(COND).LE.ZERO_TOL) THEN
                      INTVALUESQ(2)=3.0d0*DIPOLESTRENGTH/
     '                  (4.0d0*PI*RAD**2)
                    ELSE
                      INTVALUESQ(2)=3.0d0*DIPOLESTRENGTH/
     '                  (4.0d0*PI*COND*COND*RAD**2)
                    ENDIF
                  ENDIF
                ELSE
                  INTVALUESQ(2)=0.0d0
                ENDIF
              ELSE IF(NJT.EQ.3.AND.ANAL_CHOICE(nr).EQ.15) THEN
                INTVALUE(2)=0.0d0
                COND=CE(1,NEELEM(1,nr),nx)
                np=NPNODE(1,nr)
                IF(nc.EQ.1) THEN
                  IVAL=1
                ELSE
                  IVAL=4
                ENDIF
                SUM=0.0d0
                SUM1=0.0d0
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(INTERFACE.AND.NIT(NBJ(1,ne)).EQ.3) THEN
                    nef=0
                    DO neff=1,NFE(NBJ(1,ne))
                      MATCH=.TRUE.
                      DO nn=1,NNF(0,neff,NBJ(1,ne))
                        npp=NPNE(NNF(1+nn,neff,NBJ(1,ne)),NBJ(1,ne),ne)
                        IF(NP_INTERFACE(npp,1).NE.nr) MATCH=.FALSE.
                      ENDDO !nn
                      IF(MATCH) THEN
                        nef=neff
                      ENDIF
                    ENDDO !nef
                    CALL ASSERT(nef.NE.0,'>>Could not find face',
     '                ERROR,*9999)
                    nf=NFF(nef,ne)
                    FACE=.TRUE.
                    CONT=.TRUE.
                  ELSE
                    np=NPNE(1,NBJ(1,ne),ne)
                    CONT=(.NOT.INTERFACE).OR.(INTERFACE.AND.
     '                NP_INTERFACE(np,1).EQ.nr)
                    FACE=.FALSE.
                  ENDIF
                  IF(CONT) THEN
                    IF(FACE) THEN
                      CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
     '                  NBJF(1,nf),nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,
     '                  NPNE(1,1,ne),NPNF,nr,NVJE(1,1,1,ne),NVJF,
     '                  SE(1,1,ne),SF,ERROR,*9999)
                      CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,NVJF,
     '                  SF,XA(1,1,1),XE,XP,ERROR,*9999)
                      nb=NBJF(1,nf)
                      np=NPNF(1,nb)
                    ELSE
                      CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '                  NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '                  XA(1,1,ne),XE,XP,ERROR,*9999)
                      nb=NBJ(1,ne)
                      np=NPNE(1,nb,ne)
                    ENDIF
                    RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2+
     '                XP(1,1,3,np)**2)
                    DO ng=1,NGT(nb)
                      IF(FACE) THEN
                        CALL XEXG(NBJF(1,nf),ng,nr,PG,XE,XG,ERROR,*9999)
                        CALL XGMG(0,NIT(nb),NBJF(1,nf),nr,DXIX,GL,
     '                    GU,RG(ng),XG,ERROR,*9999)
                      ELSE
                        CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,ERROR,*9999)
                        CALL XGMG(0,NIT(nb),NBJ(1,ne),nr,DXIX,GL,
     '                    GU,RG(ng),XG,ERROR,*9999)
                      ENDIF
                      ZZ(1)=XG(1,1)
                      ZZ(2)=XG(2,1)
                      ZZ(3)=XG(3,1)
                      CALL ZX(3,ZZ,XX)
                      XX(1)=RAD
                      CALL XZ(3,XX,ZZ)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME(1,1,nx),
     '                  DIPOLE_DIR_NTIME(1,1,nx),IVAL,NDIPOLES(1,nx),
     '                  np,NP_INTERFACE,nr,nx,CE,DIPOLE_CEN(1,0,1,1,nx),
     '                  DIPOLE_DIR(1,0,1,1,nx),XP,W,ZZ(1),ZZ(2),ZZ(3),
     '                  ERROR,*9999)
                      SUM=SUM+W*RG(ng)*WG(ng,nb)
                      SUM1=SUM1+W**2*RG(ng)*WG(ng,nb)
                    ENDDO !ng
                  ENDIF
                ENDDO
                INTVALUE(2)=SUM
                INTVALUESQ(2)=SUM1
                CALCANAL=.TRUE.
              ELSE IF(NJT.EQ.2) THEN
                IF(ANAL_CHOICE(nr).GE.1.AND.ANAL_CHOICE(nr).LE.3) THEN
                  CALCANAL=.TRUE.
                  INTVALUE(2)=0.0d0
                  np=NPNODE(1,nr)
                  RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2)
                  IF(ANAL_CHOICE(nr).EQ.1) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*2.0d0*PI*RAD**3
                    ELSE
                      INTVALUESQ(2)=ANAL_K*2.0d0*PI*RAD
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*PI*RAD**5
                    ELSE
                      INTVALUESQ(2)=ANAL_K*2.0d0*PI*RAD**3
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*2.0d0*PI*RAD**5
                    ELSE
                      INTVALUESQ(2)=ANAL_K*8.0d0*PI*RAD**3
                    ENDIF
                  ENDIF
                ELSE
                  CALCANAL=.FALSE.
                ENDIF
              ELSE IF(NJT.EQ.3) THEN
                IF(ANAL_CHOICE(nr).GE.1.AND.ANAL_CHOICE(nr).LE.11) THEN
                  CALCANAL=.TRUE.
                  INTVALUE(2)=0.0d0
                  np=NPNODE(1,nr)
                  RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2+
     '              XP(1,1,3,np)**2)
                  IF(ANAL_CHOICE(nr).EQ.1) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*8.0d0/3.0d0*PI*RAD**4
                    ELSE
                      INTVALUESQ(2)=ANAL_K*8.0d0/3.0d0*PI*RAD**2
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*4.0d0/3.0d0*PI*RAD**4
                    ELSE
                      INTVALUESQ(2)=ANAL_K*4.0d0/3.0d0*PI*RAD**2
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*16.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*64.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*16.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*64.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*16.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*64.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*48.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*192.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*48.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*192.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*48.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*192.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*32.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*128.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*32.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*128.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN
                    IF(nc.EQ.1) THEN
                      INTVALUESQ(2)=ANAL_K*32.0d0/15.0d0*PI*RAD**6
                    ELSE
                      INTVALUESQ(2)=ANAL_K*128.0d0/15.0d0*PI*RAD**4
                    ENDIF
                  ENDIF
                ELSE
                  CALCANAL=.FALSE.
                ENDIF
              ENDIF
            ELSE
              CALCANAL=.FALSE.
            ENDIF

c cpb 4/1/97 Adding integral of difffernce
c cpb 28/8/97 Adding normalised integral of difference squared (NID*)
            DO i=1,4
C             Numerical (i=1), anal (i=2), diff (i=3), norm diff (i=4)
              IF(.NOT.(i.EQ.2.AND.CALCANAL)) THEN
                IF(i.EQ.2) THEN
                  iy=7
                ELSE
                  iy=1
                ENDIF
                INTVALUE(i)=0.0d0
                INTVALUESQ(i)=0.0d0
                CALL YPZP(iy,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH,
     '            NPNODE,nr,NVHP,nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,
     '            *9999)
                ANALDIFF=.FALSE.
                IF(i.EQ.3.OR.i.EQ.4) THEN
                  IF(SPHERE.AND.((NJT.EQ.2.AND.(ANAL_CHOICE(nr).EQ.3.OR.
     '              ANAL_CHOICE(nr).EQ.5)).OR.(NJT.EQ.3.AND.
     '              (ANAL_CHOICE(nr).EQ.6.OR.ANAL_CHOICE(nr).EQ.12))))
     '              THEN
                    ANALDIFF=.TRUE.
                    np=NPNODE(1,nr)
                    IF(NJT.EQ.2) THEN
                      RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2)
                      IF(ANAL_CHOICE(nr).EQ.5.AND.
     '                  ITYP3(nr,nx).EQ.2) THEN
                        COND=CE(1,NEELEM(1,nr),nx)
                      ELSE
                        COND=1.0d0
                      ENDIF
                    ELSE
                      RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2+
     '                  XP(1,1,3,np)**2)
                      IF(ANAL_CHOICE(nr).EQ.12.AND.
     '                  ITYP3(nr,nx).EQ.2) THEN
                        COND=CE(1,NEELEM(1,nr),nx)
                      ELSE
                        COND=1.0d0
                      ENDIF
                    ENDIF
                  ELSE IF(SPHERE.AND.NJT.EQ.3.AND.
     '                ANAL_CHOICE(nr).EQ.15) THEN
                    RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2+
     '                XP(1,1,3,np)**2)
                    ANALDIFF=.TRUE.
                  ELSE
                    DO nonode=1,NPNODE(0,nr)
                      np=NPNODE(nonode,nr)
                      DO nhx=1,NHP(np,nr,nx)
                        nh=NH_LOC(nhx,nx)
                        DO nv=1,NVHP(nh,np,nc)
                          DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                            ny=NYNP(nk,nv,nh,np,0,nc,nr)
                            IF(i.EQ.3) THEN
                              ZP(nk,nv,nh,np,nc)=YP(ny,7,nx)-
     '                          ZP(nk,nv,nh,np,nc)
                            ELSE
                              IF(INTVALUESQ(1).GT.ZERO_TOL.AND.
     '                          INTVALUESQ(2).GT.ZERO_TOL) THEN
                                ZP(nk,nv,nh,np,nc)=YP(ny,7,nx)/
     '                            DSQRT(INTVALUESQ(2))-
     '                            ZP(nk,nv,nh,np,nc)/
     '                            DSQRT(INTVALUESQ(1))
                              ELSE
                                ZP(nk,nv,nh,np,nc)=0.0d0
                              ENDIF
                            ENDIF
                          ENDDO !nk
                        ENDDO !nv
                      ENDDO !nh
                    ENDDO !nonode (np)
                  ENDIF
                ENDIF
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(INTERFACE.AND.NIT(NBJ(1,ne)).EQ.3) THEN
                    nef=0
                    DO neff=1,NFE(NBJ(1,ne))
                      MATCH=.TRUE.
                      DO nn=1,NNF(0,neff,NBJ(1,ne))
                        npp=NPNE(NNF(1+nn,neff,NBJ(1,ne)),NBJ(1,ne),ne)
                        IF(NP_INTERFACE(npp,1).NE.nr) MATCH=.FALSE.
                      ENDDO !nn
                      IF(MATCH) THEN
                        nef=neff
                      ENDIF
                    ENDDO !nef
                    CALL ASSERT(nef.NE.0,'>>Could not find face',
     '                ERROR,*9999)
                    nf=NFF(nef,ne)
                    FACE=.TRUE.
                    CONT=.TRUE.
                  ELSE
                    np=NPNE(1,NBJ(1,ne),ne)
                    CONT=(.NOT.INTERFACE).OR.(INTERFACE.AND.
     '                NP_INTERFACE(np,1).EQ.nr)
                    FACE=.FALSE.
                  ENDIF
                  IF(CONT) THEN
                    IF(FACE) THEN
                      CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
     '                  NBJF(1,nf),nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,
     '                  NPNE(1,1,ne),NPNF,nr,NVJE(1,1,1,ne),NVJF,
     '                  SE(1,1,ne),SF,ERROR,*9999)
                      CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,NVJF,
     '                  SF,XA(1,1,1),XE,XP,ERROR,*9999)
                      CALL CALC_FACE_INFORMATION_DEP(NBH(1,nc,ne),
     '                  NBHF(1,nc,nf),nef,NHE(ne,nx),
     '                  NKHE(1,1,1,ne),NKEF,NKHF,NNF,NPNE(1,1,ne),
     '                  NPNF,NVHE(1,1,1,ne),NVHF,nx,SE(1,1,ne),SF,
     '                  ERROR,*9999)
                      CALL ZPZE(NBHF(1,1,nf),nc,NHE(ne,nx),
     '                  NKHF,NPF(1,nf),NPNF,nr,NVHF,NW(ne,1,nx),nx,
     '                  CURVCORRECT(1,1,1,ne),SF,ZA(1,1,1,ne),
     '                  ZE,ZP,ERROR,*9999)
                      nb=NBHF(NH_LOC(1,nx),nc,nf)
                      np=NPNF(1,nb)
                    ELSE
                      CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '                  NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '                  XA(1,1,ne),XE,XP,ERROR,*9999)
                      CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),
     '                  NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '                  NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '                  CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '                  ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                      nb=NBH(NH_LOC(1,nx),nc,ne)
                      np=NPNE(1,nb,ne)
                    ENDIF
                    IF(NJT.EQ.2) THEN
                      RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2)
                    ELSE
                      RAD=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2+
     '                  XP(1,1,3,np)**2)
                    ENDIF
                    SUM=0.0d0
                    SUM1=0.0d0
                    DO ng=1,NGT(nb)
                      IF(FACE) THEN
                        CALL XEXG(NBJF(1,nf),ng,nr,PG,XE,XG,ERROR,*9999)
                        CALL XGMG(0,NIT(nb),NBJF(1,nf),nr,DXIX,GL,
     '                    GU,RG(ng),XG,ERROR,*9999)
                      ELSE
                        CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,ERROR,*9999)
                        CALL XGMG(0,NIT(nb),NBJ(1,ne),nr,DXIX,GL,
     '                    GU,RG(ng),XG,ERROR,*9999)
                      ENDIF
                      IF(ANALDIFF) THEN
                        IF(NJT.EQ.2) THEN
                          ZZ(1)=XG(1,1)
                          ZZ(2)=XG(2,1)
                          CALL ZX(2,ZZ,XX)
                          THETA=XX(2)
                          IF(ANAL_CHOICE(nr).EQ.3) THEN
                            XX(1)=RAD
                            CALL XZ(2,XX,ZZ)
                            IF(nc.EQ.1) THEN
                              EXACT=ANAL_K*(ZZ(1)**2+2.0d0*ZZ(1)*ZZ(2)-
     '                          ZZ(2)**2)
                            ELSE
                              EXACT=2.0d0*ANAL_K*(ZZ(1)**2+2.0d0*ZZ(1)*
     '                          ZZ(2)-ZZ(2)**2)/RAD
                            ENDIF
                          ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN
                            IF(nc.EQ.1) THEN
                              EXACT=-(ANAL_A*DCOS(THETA)+ANAL_B*
     '                          DSIN(THETA))/(PI*COND*RAD)
                            ELSE
                              EXACT=0.0d0
                            ENDIF
                          ENDIF
                        ELSE
                          ZZ(1)=XG(1,1)
                          ZZ(2)=XG(2,1)
                          ZZ(3)=XG(3,1)
                          CALL ZX(3,ZZ,XX)
                          THETA=XX(2)
                          PHI=XX(3)
                          IF(ANAL_CHOICE(nr).EQ.6) THEN
                            XX(1)=RAD
                            CALL XZ(3,XX,ZZ)
                            IF(nc.EQ.1) THEN
                              EXACT=ANAL_K*(ZZ(1)**2+ZZ(2)**2-2.0d0*
     '                          ZZ(3)**2)
                            ELSE
                              EXACT=2.0d0*ANAL_K*(ZZ(1)**2+ZZ(2)**2-
     '                          2.0d0*ZZ(3)**2)/RAD
                            ENDIF
                          ELSE IF(ANAL_CHOICE(nr).EQ.12) THEN
                            IF(nc.EQ.1) THEN
                              EXACT=3.0d0/(4.0d0*PI*COND*RAD**2)*
     '                          (ANAL_A*DCOS(THETA)*DCOS(PHI)+
     '                          ANAL_B*DSIN(THETA)*DCOS(PHI)+
     '                          ANAL_C*DSIN(PHI))
                            ELSE
                              EXACT=0.0d0
                            ENDIF
                          ELSE IF(ANAL_CHOICE(nr).EQ.15) THEN
                            XX(1)=RAD
                            CALL XZ(3,XX,ZZ)
                            CALL DIPOLE_EVALUATE(
     '                        DIPOLE_CEN_NTIME(1,1,nx),
     '                        DIPOLE_DIR_NTIME(1,1,nx),IVAL,
     '                        NDIPOLES(1,nx),np,NP_INTERFACE,
     '                        nr,nx,CE,DIPOLE_CEN(1,0,1,1,nx),
     '                        DIPOLE_DIR(1,0,1,1,nx),XP,W,ZZ(1),
     '                        ZZ(2),ZZ(3),ERROR,*9999)
                            EXACT=W
                          ENDIF
                        ENDIF
                        INT=0.0d0
                        DO ns=1,NST(nb)
                          INT=INT+PG(ns,1,ng,nb)*ZE(ns,NH_LOC(1,nx))
                        ENDDO !ns
                        IF(i.EQ.3) THEN
                          SUM=SUM+(EXACT-INT)*RG(ng)*WG(ng,nb)
                          SUM1=SUM1+(EXACT-INT)**2*RG(ng)*WG(ng,nb)
                        ELSE
                          IF(INTVALUESQ(1).GT.ZERO_TOL.AND.
     '                      INTVALUESQ(2).GT.ZERO_TOL) THEN
                            SUM=SUM+(EXACT/DSQRT(INTVALUESQ(2))-
     '                        INT/DSQRT(INTVALUESQ(1)))*RG(ng)*WG(ng,nb)
                            SUM1=SUM1+(EXACT/DSQRT(INTVALUESQ(2))-
     '                        INT/DSQRT(INTVALUESQ(1)))**2*RG(ng)*
     '                        WG(ng,nb)
                          ELSE
                            SUM=0.0d0
                            SUM1=0.0d0
                          ENDIF
                        ENDIF
                      ELSE
                        INT=0.0d0
                        DO ns=1,NST(nb)
                          INT=INT+PG(ns,1,ng,nb)*ZE(ns,NH_LOC(1,nx))
                        ENDDO !ns
                        SUM=SUM+INT*RG(ng)*WG(ng,nb)
                        SUM1=SUM1+INT**2*RG(ng)*WG(ng,nb)
                      ENDIF
                    ENDDO !ng
                    INTVALUE(i)=INTVALUE(i)+SUM
                    INTVALUESQ(i)=INTVALUESQ(i)+SUM1
                  ENDIF
                ENDDO !ne
              ENDIF
            ENDDO !i (numerical and analtyic)
            EXACT=INTVALUE(2)
            ABSOLUTE=DABS(EXACT-INTVALUE(1))
C cpb 26/1/96 Changing the way %'s are calculated
C            PERCENT=DABS((EXACT-YP(ny,1,nx))/
C     '        (1.0D0+DABS(EXACT)))*100.0d0
            IF(DABS(EXACT).LE.0.001d0) THEN
              PERCENT=(EXACT-INTVALUE(1))/(1.0d0+DABS(EXACT))*100.0d0
            ELSE
              PERCENT=(EXACT-INTVALUE(1))/DABS(EXACT)*100.0d0
            ENDIF
            RELATIVE=ABSOLUTE/(1.0d0+DABS(EXACT))
            FORMAT='(2X,''Intgl'',1X,D11.4,1X,D11.4,1X,F8.2,1X,'
     '        //'D14.7,1X,D14.7)'
            WRITE(OP_STRING,FORMAT) INTVALUE(1),EXACT,PERCENT,
     '        ABSOLUTE,RELATIVE
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            EXACT=INTVALUESQ(2)
            ABSOLUTE=DABS(EXACT-INTVALUESQ(1))
            IF(DABS(EXACT).LE.0.001d0) THEN
              PERCENT=(EXACT-INTVALUESQ(1))/(1.0d0+DABS(EXACT))*100.0d0
            ELSE
              PERCENT=(EXACT-INTVALUESQ(1))/DABS(EXACT)*100.0d0
            ENDIF
            RELATIVE=ABSOLUTE/(1.0d0+DABS(EXACT))
            FORMAT='(2X,''Int^2'',1X,D11.4,1X,D11.4,1X,F8.2,1X,'
     '        //'D14.7,1X,D14.7)'
            WRITE(OP_STRING,FORMAT) INTVALUESQ(1),EXACT,PERCENT,
     '        ABSOLUTE,RELATIVE
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='(10X,''Numerical'',4X,'
     '        //'''     NID                    NID*'')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(DABS(INTVALUE(2)).GT.ZERO_TOL) THEN
              FORMAT='(2X,''Diff.'',1X,D11.4,1X,D11.4,1X,''('',F8.2,'
     '          //'''%)'')'
              WRITE(OP_STRING,FORMAT) INTVALUE(3),
     '          INTVALUE(3)/INTVALUE(2),INTVALUE(3)/INTVALUE(2)*100.0d0
            ELSE
              FORMAT='(2X,''Diff.'',1X,D11.4,1X,'
     '          //'''   Infinity ( Infinity)'')'
              WRITE(OP_STRING,FORMAT) INTVALUE(3)
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(INTVALUESQ(1).GT.ZERO_TOL.AND.
     '        INTVALUESQ(2).GT.ZERO_TOL) THEN
              FORMAT='(2X,''Dif^2'',1X,D11.4,1X,D11.4,1X,''('',F8.2,'
     '          //'''%)'',1X,D11.4,1X,''('',F8.2,''%)'')'
               WRITE(OP_STRING,FORMAT) INTVALUESQ(3),
     '          DSQRT(INTVALUESQ(3)/INTVALUESQ(2)),
     '          DSQRT(INTVALUESQ(3)/INTVALUESQ(2))*100.0d0,
     '          DSQRT(INTVALUESQ(4)),DSQRT(INTVALUESQ(4))*100.0d0
            ELSE IF(INTVALUESQ(2).GT.ZERO_TOL) THEN
              FORMAT='(2X,''Dif^2'',1X,D11.4,1X,D11.4,1X,''('',F8.2,'
     '          //'''%)'',1X,''   Infinity ( Infinity)'')'
               WRITE(OP_STRING,FORMAT) INTVALUESQ(3),
     '          DSQRT(INTVALUESQ(3)/INTVALUESQ(2)),
     '          DSQRT(INTVALUESQ(3)/INTVALUESQ(2))*100.0d0
            ELSE
              FORMAT='(2X,''Dif^2'',1X,D11.4,1X,''   Infinity'',1X,'
     '          //'''( Infinity)    Infinity ( Infinity)'')'
              WRITE(OP_STRING,FORMAT) INTVALUESQ(3)
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nc
        ENDDO !End nr

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('CHKSOL')
      RETURN
 9999 CALL ERRORS('CHKSOL',ERROR)
      CALL EXITS('CHKSOL')
      RETURN 1
      END


