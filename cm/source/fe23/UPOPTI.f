      SUBROUTINE UPOPTI(CE,CELL_RCQS_VALUE,CP,DIPOLE_CEN_NTIME,
     &  DIPOLE_DIR_NTIME,
     '  ILPIN,ISIZE_PHI,NDIPOLES,NEELEM,NMNO,NONY,
     '  NP1OPT,NP2OPT,NP3OPT,NP_INTERFACE,NPNY,NRLIST,NXLIST,
     '  NYNO,NYNP,NYNR,PAOPTY,
     '  CQ,DIPOLE_CEN,DIPOLE_DIR,PAOPTI,PBOPTI,XP,YP,STRING,ERROR,*)

C#### Subroutine: UPOPTI
C###  Description:
C###    UPOPTI updates finite element variables from optimisation
C###    parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'aero00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER ISIZE_PHI(2),DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),ILPIN(NMM,NRM,NXM),
     '  NDIPOLES(NRM),NEELEM(0:NE_R_M,0:NRM),
     '  NMNO(1:2,0:NOPM,NXM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NP1OPT(NOPM),NP2OPT(NOPM),NP3OPT(NOPM),
     '  NP_INTERFACE(0:NPM,0:3),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NXLIST(0:NXM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),NRLIST(0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),PAOPTY(NOPM)
      REAL*8 CE(NMM,NEM,NXM),CELL_RCQS_VALUE(NQRM,NQVM),
     &  CP(NMM,NPM,NXM),CQ(NMM,NQM,NXM),
     '  DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  PAOPTI(*),PBOPTI(*),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER ERR,IBEG,IEND,N3CO,ndipole,ne,nj,nh,nk,nm,no,
     '  no_interface,no_nynr,noelem,noopti,nonr,noy,np,NP1,NP2,
     '  nq,nr,nt,
     '  nx_opt,nx_sol,nxc,SSTART,SEND,ny,nyo
      REAL*8 SLOPE,TSTART,TEND,X_aero,Y_aero
      LOGICAL ALL_REGIONS,CBBREV,DOING_OPT

!     Functions
      INTEGER CALC_SAMPLE_FROM_TIME
      REAL*8 RFROMC

      CALL ENTERS('UPOPTI',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update optimisation
C###  Parameter:      <(optimise/solve)[optimise]>
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <nr (#s/all)[1]>
C###    Specify the region number(s)
C###  Parameter:      <tstart (#)[0.0]>
C###    Specify the start time for time dependent problems
C###  Parameter:      <tend (#)[0.0]>
C###    Specify the end time for time dependent problems
C###  Description: updates finite element variables from optimisation
C###    parameters.

        OP_STRING(1)=STRING(1:IEND)//' <(optimise/solve)[optimise]>'
        OP_STRING(2)=BLANK(1:15)//' <class #[1]>'
        OP_STRING(3)=BLANK(1:15)//' <region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//' <tstart #[0.0]>'
        OP_STRING(5)=BLANK(1:15)//' <tend #[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPOPTI',ERROR,*9999)
      ELSE

C LKC 2-MAY-2002 remove temporary nr=1
C        nr=1 !temporary
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'SOLVE',2,noco+1,NTCO,N3CO)) THEN
          DOING_OPT=.FALSE.
        ELSE
          DOING_OPT=.TRUE.
        ENDIF

        IF(CBBREV(CO,'TSTART',3,noco+1,NTCO,N3CO)) THEN
          TSTART=RFROMC(CO(N3CO+1))
        ELSE
          TSTART=0.D0
        ENDIF

        IF(CBBREV(CO,'TEND',3,noco+1,NTCO,N3CO)) THEN
          TEND=RFROMC(CO(N3CO+1))
        ELSE
          TEND=0.D0
        ENDIF

        ERR=0
        SSTART=CALC_SAMPLE_FROM_TIME(TSTART,ERR,ERROR)
        IF(ERR.NE.0) GOTO 9999
        SEND=CALC_SAMPLE_FROM_TIME(TEND,ERR,ERROR)
        IF(ERR.NE.0) GOTO 9999

        IF(DOING_OPT) THEN
          IF(KTYP26.EQ.1) THEN !material parameters
            CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
            CALL ASSERT(nx_opt.GT.0,
     '        '>>No nx defined for this optimisation class',
     '        ERROR,*9999)
            CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx_sol.GT.0,
     '        '>>No nx defined for this solve class',
     '        ERROR,*9999)
            IF(ITYP5(nr,nx_sol).EQ.2.AND.ITYP19(nr,nx_sol).EQ.1.AND.
     '        ITYP2(nr,nx_sol).EQ.9) THEN
              !activation model
              DO noopti=1,NMNO(1,0,nx_opt)
                nm=NMNO(1,noopti,nx_opt)
                DO nq=1,NQT
                  CQ(nm,nq,nx_opt)=PAOPTI(noopti)
                ENDDO
c              IF(DOP) THEN
                WRITE(OP_STRING,'('' CQ('',I1,'',nq,nx_opt)='',D10.4)')
     '            nm,PAOPTI(noopti)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              ENDIF
              ENDDO
            ELSE
C KFA 2003/01/14 Added general material update from opti (for fe prob)
              DO noopti=1,NTOPTI
                nm=NMNO(1,noopti,nx_opt)
C news
C HS, 14/5/04, added loop for constitutivelaw via CellMl envrionment
                IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN !grid coupling
                  CELL_RCQS_VALUE(nm,1)=PAOPTI(noopti)
                ELSE
C newe
                  IF(ILPIN(nm,nr,nx_sol).EQ.1) THEN
                    DO noelem=1,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      CE(nm,ne,nx_sol)=PAOPTI(noopti)
                    ENDDO
C                    IF(DOP) THEN
                      WRITE(OP_STRING,
     '                  '('' CE('',I2,'',.,'',I3,'')='',D10.4)')
     '                  nm,nx_sol,PAOPTI(noopti)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    ENDIF
                  ELSE IF(ILPIN(nm,nr,nx_sol).EQ.2) THEN
                    CE(nm,NMNO(2,noopti,nx_opt),nx_sol)=PAOPTI(noopti)
                    IF(DOP) THEN
                      WRITE(OP_STRING,
     '                  '('' CE('',I2,'','',I3,'','',I3,'')='',D10.4)')
     '                  nm,NMNO(2,noopti,nx_opt),nx_sol,PAOPTI(noopti)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE IF(ILPIN(nm,nr,nx_sol).EQ.3) THEN
                    CP(nm,NMNO(2,noopti,nx_opt),nx_sol)=PAOPTI(noopti)
                  ELSE
                    ERROR='Material Parameters need to be constant, '
     '             //'or piece-wise constant'
                    GOTO 9999
                  ENDIF !ILPIN
                ENDIF !KTYP54 & KTYP3B
              ENDDO !NTOPTI
            ENDIF

          ELSE IF(KTYP26.EQ.2) THEN !geometric parameters
 ! Update nodal coordinates from optimisation params
            IF(KTYP27.EQ.4) THEN !Fluid interface condition
              DO noopti=1,NTOPTI
                NP1=NP1OPT(noopti)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
                XP(1,1,2,NP1)=PAOPTI(noopti)/PBOPTI(noopti)*
     '            XP(1,1,2,NP1)
                PBOPTI(noopti)=PAOPTI(noopti)
              ENDDO

            ELSE IF(KTYP27.EQ.7) THEN !Aero wake pres diff & sail stress
              IF(N_OPTI(1).GE.NL_WAKE(0,1)) THEN
 !             Update wake node y-coords
                DO noopti=1,NL_WAKE(0,1)
                  NP1=NP1OPT(noopti)
                  NP2=NP2OPT(noopti)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
                  XP(1,1,2,NP1)=PAOPTI(noopti)
                  XP(1,1,2,NP2)=XP(1,1,2,NP1)
                  PBOPTI(noopti)=PAOPTI(noopti)
                ENDDO
              ENDIF

              IF(N_OPTI(2).GT.0) THEN !sail stress parameters included
                CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
                CALL ASSERT(nx_opt.GT.0,
     '            '>>No nx defined for this optimisation class',
     '            ERROR,*9999)
 !             Update sail node y-coords
                nr=2 !is assumed sail stress region
                DO noopti=N_OPTI(1)+1,NTOPTI
c cpb this should probably be done through nyno and npny not np#opt
                  np=NP1OPT(noopti)
                  nh=NP2OPT(noopti)
                  nk=NP3OPT(noopti)
                  ny=NYNP(nk,1,nh,np,0,1,nr)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
 !               Update stress field deformed node positions
                  YP(ny,4,nx_opt)=PAOPTI(noopti)
 !               Update flow field undeformed node positions
                  DO no_interface=1,NP_INTERFACE(0,1)
                    IF(NP_INTERFACE(no_interface,3).eq.np) THEN
                      XP(nk,1,nh,NP_INTERFACE(no_interface,1))=
     '                  PAOPTI(noopti)
                      XP(nk,1,nh,NP_INTERFACE(no_interface,2))=
     '                  PAOPTI(noopti)
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

            ELSE IF(KTYP27.EQ.8) THEN !Aerofoil lift & wake press. diff.
 !           Update wake node y-coords (if wake optimisation included)
              DO noopti=1,N_OPTI(1)
                NP1=NP1OPT(noopti)
                NP2=NP2OPT(noopti)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
                XP(1,1,2,NP1)=PAOPTI(noopti)
                XP(1,1,2,NP2)=XP(1,1,2,NP1)
                PBOPTI(noopti)=PAOPTI(noopti)
              ENDDO

 !           Calc. slope of aerofoil centreline
              SLOPE=(XP(1,1,2,NP_aero_TE1)-XP(1,1,2,NP_aero_LE))
     '          /(XP(1,1,1,NP_aero_TE1)-XP(1,1,1,NP_aero_LE))
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' Slope='',E13.5)') SLOPE
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF

 !           Update aerofoil node y-coords
              DO noopti=N_OPTI(1)+1,NTOPTI
                NP1=NP1OPT(noopti) !upper surface node
                NP2=NP2OPT(noopti) !lower surface node
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
 !             Calc. new upper surface node position
                XP(1,1,2,NP1)=PAOPTI(noopti) !new y value
 !             Calc. new lower surface node position
                IF(DABS(SLOPE).GT.1.0d-6) THEN
                  Y_aero=(XP(1,1,1,NP1)-XP(1,1,1,NP_aero_LE)
     '              +SLOPE*XP(1,1,2,NP1)
     '              +XP(1,1,2,NP_aero_LE)/SLOPE)/(SLOPE+1.0d0/SLOPE)
                  X_aero=XP(1,1,1,NP_aero_LE)
     '              +(Y_aero-XP(1,1,2,NP_aero_LE))/SLOPE
                  XP(1,1,1,NP2)=2.0d0*X_aero-XP(1,1,1,NP1)
                  XP(1,1,2,NP2)=2.0d0*Y_aero-XP(1,1,2,NP1)
                ELSE
                  XP(1,1,2,NP2)=-XP(1,1,2,NP1) !symmet about horiz axis
                ENDIF
                PBOPTI(noopti)=PAOPTI(noopti)
              ENDDO
            ELSEIF(KTYP27.EQ.12) THEN !activation times
              CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
              CALL ASSERT(nx_opt.GT.0,
     '          '>>No nx defined for this optimisation class',
     '          ERROR,*9999)
              CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx_sol.GT.0,
     '          '>>No nx defined for this solve class',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '          ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(1).GT.0,'>>Evaluate PHI first',
     '          ERROR,*9999)
              CALL ASSERT(NTOPTI.GT.0,'>>Zero optimisation variables?',
     '          ERROR,*9999)
              CALL ASSERT(CALL_OPTI,'>>Define optimisation first',
     '          ERROR,*9999)
C             Copy in the current estimates of YP (activation times)
C             from PAOPTI
              nr=TRSF_NR_FIRST
              DO noopti=1,NTOPTI
                IF(PAOPTY(noopti).EQ.1) THEN
                  DO nyo=1,NYNO(0,noopti,2,nr,nx_opt)
                    ny=NYNO(nyo,noopti,2,nr,nx_opt)
                    nk=NPNY(1,ny,0,nx_sol)
C                    nv=NPNY(2,ny,0,nx_sol)
C                    nj=NPNY(3,ny,0,nx_sol)
                    np=NPNY(4,ny,0,nx_sol)
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    YP(ny,1,nx_sol)=PAOPTI(noopti)
                  ENDDO !nyo
                ELSE IF(PAOPTY(noopti).EQ.2) THEN
                  TRSF_ACTN_WAVE_JUMP=PAOPTI(noopti)
                ELSE IF(PAOPTY(noopti).EQ.3) THEN
                  TRSF_ACTN_WAVE_WIDTH=PAOPTI(noopti)
                ELSE
                  ERROR='>>Invalid PAOPTY for activation optimisation'
                  GOTO 9999
                ENDIF
              ENDDO


C*** Dipole position and orientation optimisation

            ELSEIF(KTYP27.EQ.13) THEN !dipole

C LKC 7-JUL-2005 Should really be using GET_DIPOLE as the dipole
C   storage doesn't really have anything to do with the MFI/SAMPLE
C              
C              ERR=0
C              SSTART=CALC_SAMPLE_FROM_TIME(TSTART,ERR,ERROR)
C              IF(ERR.NE.0) GOTO 9999
C              SEND=CALC_SAMPLE_FROM_TIME(TEND,ERR,ERROR)
C              IF(ERR.NE.0) GOTO 9999

              CALL ASSERT(SSTART.EQ.SEND,
     '          '>> Only implemented for 1 time step',ERROR,*9999)
              CALL ASSERT(NRLIST(0).EQ.1,
     '          '>> Only implemented for 1 region',ERROR,*9999)
              CALL ASSERT(NDIPOLES(nr).GT.0,
     '          '>> No dipoles in optimisation region',ERROR,*9999)
              CALL ASSERT(NDIPOLES(nr).EQ.1,
     '          '>> Only implemented for 1 dipole',ERROR,*9999)

              DO nonr=1,NRLIST(0)
                nr=NRLIST(nonr)
                DO ndipole=1,NDIPOLES(nr)

                  IF(DIPOLE_CEN_NTIME(ndipole,nr).EQ.0) THEN
                    SSTART=1
                    SEND=1
                  ENDIF !time dependent

                  DO nt=SSTART,SEND
                    DO nj=1,NJT ! dipole position/orientation
C LKC 14-JUL-2005 I think the -1 is to compensate the way GET_DIPOLE calculates the
C  time sample
                      DIPOLE_CEN(nj,nt-1,ndipole,nr)=PAOPTI(nj)
C                      DIPOLE_CEN(nj,nt,ndipole,nr)=PAOPTI(nj)
                    ENDDO !nj

C                    IF(DOP) THEN
C                      WRITE(OP_STRING,'('' UPOPTI center '',I3,3F15.9)')
C     '                  nt,(PAOPTI(nj),nj=1,3)
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    ENDIF

                  ENDDO !nt

                  IF(DIPOLE_DIR_NTIME(ndipole,nr).EQ.0) THEN
                    SSTART=1
                    SEND=1
                  ENDIF
                  DO nt=SSTART,SEND
                    DO nj=1,NJT ! dipole position/orientation
C LKC 14-JUL-2005 I think the -1 is to compensate the way GET_DIPOLE calculates the
C  time sample
                      DIPOLE_DIR(nj,nt-1,ndipole,nr)=PAOPTI(nj+NJT)
C                      DIPOLE_DIR(nj,nt,ndipole,nr)=PAOPTI(nj+NJT)
                    ENDDO !nj
                    
C                   IF(DOP) THEN
C                      WRITE(OP_STRING,'('' UPOPTI orient '',I3,3F15.9)')
C     '                  nt,(PAOPTI(nj+NJT),nj=1,3)
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    ENDIF

                  ENDDO !nt
                ENDDO !ndipole
              ENDDO !nr

            ELSE
              ERROR='Unknown KTYP27'
              GOTO 9999
            ENDIF !ktyp27
          ENDIF !ktyp26


        ELSE !doing a solve using E04UPF for FE50 problems
          CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_sol.GT.0,
     '      '>>No nx defined for this solve class',ERROR,*9999)
          nr=1 !temporary MPN 17/9/93
          DO no_nynr=1,NYNR(0,0,1,nr,nx_sol) !loop over global variables
            ny=NYNR(no_nynr,0,1,nr,nx_sol) !is global variable number
            IF(NPNY(0,ny,0,nx_sol).EQ.1) THEN
              np=NPNY(4,ny,0,nx_sol)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF
            DO noy=1,NONY(0,ny,2,nr,nx_sol)
              no=NONY(noy,ny,2,nr,nx_sol)
              YP(ny,1,nx_sol)=PBOPTI(no) !PBOPTI set up in FUNCT4
            ENDDO !noy (no)
          ENDDO !no_nynr (ny)

        ENDIF !DOING_OPT
      ENDIF

      CALL EXITS('UPOPTI')
      RETURN
 9999 CALL ERRORS('UPOPTI',ERROR)
      CALL EXITS('UPOPTI')
      RETURN 1
      END



