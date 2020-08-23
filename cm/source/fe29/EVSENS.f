      SUBROUTINE EVSENS(ISTATE,IUSER,NEELEM,NMNO,NPNODE,NVJP,NYNP,
     '  PAOPTI,PE,PF,PMIN,PMAX,R,RESID,RESJAC,USER,XC,XP,YP,STRING,
     '  FIX,FIXP,ERROR,*)

C#### Subroutine: EVSENS
C###  Description:
C###    EVSENS evaluates sensitivity of optimisation parameters
C###    wrt geometry pressure bcs or force bcs.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER ISTATE(*),IUSER(*),NEELEM(0:NE_R_M,0:NRM),
     '  NMNO(1:2,0:NOPM,NXM),NPNODE(0:NP_R_M,0:NRM),NVJP(NJM,NPM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 PAOPTI(*),PE(2,NEM),PF(2,NEM),PMIN(*),PMAX(*),R(NOPM,*),
     '  RESID(*),RESJAC(NREM,*),USER(*),XC(*),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXP(2,NEM)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IFROMC,IIY,IY,L_SENSIT,L_XC2,N3CO,
     '  noopti,nr,nx
      INTEGER*4 SENSIT_PTR,XC2_PTR
      CHARACTER CTEMP*30,FILE*(MXCH),TYPE*8
      LOGICAL ABBREV,CBBREV,CONSTR,DEFORMED,OUTPTPRN,OPFILE


      CALL ENTERS('EVSENS',*9999)

      nx=1 !Temporary cpb 22/11/94
      SENSIT_PTR=0
      XC2_PTR=0

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate sensitivity<;FILENAME>
C###  Parameter:        wrt nodal_parameters
C###    Evaluate wrt nodal parameters
C###  Parameter:        deformed (default)
C###    Evaluate wrt deformed coordinates.
C###  Parameter:        <(constrained/all)[constrained]>
C###    Constrains dofs
C###  Parameter:        <iy=IY[4]>
C###    Specify a varible number
C###  Parameter:        full_output
C###    Outputs further information from the optimiser
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Evauluate sensitivity of optimisation parameters.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'wrt nodal_parameters'
        OP_STRING(3)=BLANK(1:15)//'deformed (default)'
        OP_STRING(4)=BLANK(1:15)//'<(constrained/all)[constrained]>'
        OP_STRING(5)=BLANK(1:15)//'<iy=IY[4]>'
        OP_STRING(6)=BLANK(1:15)//'full_output'
        OP_STRING(7)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate sensitivity reference
C###  Parameter:        <(constrained/all)[all]>
C###    Contrain dofs
C###  Parameter:        full_output
C###    Outputs further information from the optimiser
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Evaluate a sensitivity wrt of the solution
C###    to reference state.

        OP_STRING(1)=STRING(1:IEND)//' reference'
        OP_STRING(2)=BLANK(1:15)//'<(constrained/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'full_output'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate sensitivity wrt pressure_bcs
C###  Parameter:        full_output
C###    Outputs further information from the optimiser
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Evaluates the sensitity of the solution
C###    wrt pressure boundary conditions.

        OP_STRING(1)=STRING(1:IEND)//' wrt pressure_bcs'
        OP_STRING(2)=BLANK(1:15)//'full_output'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVSENS',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opsens','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL ASSERT(CALL_OPTI,'>>Optimisation parameters not defined',
     '    ERROR,*9999)

        TYPE=' '
        IF(CBBREV(CO,'WRT',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'NODAL_PARAMETERS',1)) THEN
            TYPE='NODAL'
            IF(CBBREV(CO,'REFERENCE',1,noco+1,NTCO,N3CO)) THEN
              DEFORMED=.FALSE.
              IF(CBBREV(CO,'CONSTRAINED',1,noco+1,NTCO,N3CO)) THEN
                CONSTR=.TRUE.
              ELSE  !all dof
                CONSTR=.FALSE.
              ENDIF

            ELSE  !vary YP
              DEFORMED=.TRUE.
              IF(CBBREV(CO,'ALL',1,noco+1,NTCO,N3CO)) THEN
                CONSTR=.FALSE.
              ELSE  !constrained dof
                CONSTR=.TRUE.
              ENDIF
              IF(CBBREV(CO,'IY',1,noco+1,NTCO,N3CO)) THEN
                IY=IFROMC(CO(N3CO+1))
              ELSE
                IY=4
              ENDIF
              IF(IY.EQ.4) THEN
                IIY=3  !FIX(3) applies to YP(4)
              ELSE
                IIY=IY
              ENDIF
            ENDIF

          ELSE IF(ABBREV(CO(N3CO+1),'PRESSURE_BCS',1)) THEN
            TYPE='PRESSURE'
            DEFORMED=.FALSE.
          ENDIF
        ENDIF
        IF(TYPE(1:1).EQ.' ') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
          CALL STRING_TRIM(CO(noco),IBEG1,IEND1)
          CTEMP=CO(noco)(IBEG1:IEND1)
          CALL STRING_TRIM(CTEMP,IBEG1,IEND1)
          STRING=STRING(1:IEND)//' '//CTEMP(IBEG1:IEND1)
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(CBBREV(CO,'FULL_OUTPUT',2,noco+1,NTCO,N3CO)) THEN
          OUTPTPRN=.TRUE.
        ELSE
          OUTPTPRN=.FALSE.
        ENDIF

        IF(CBBREV(CO,'REGION',3,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Optimising for unperturbed '
     '      //'parameters'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C ***   Optimise mat params for unperturbed nodal params/press bcs
        CALL NAGMINA(ISTATE,PAOPTI,PMIN,PMAX,R,RESID,RESJAC,XC,
     '    IUSER,USER,OUTPTPRN,ERROR,*9999)

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Material param # '',5(I2,10X)/,'
     '      //'(20X,5(I2,10X)))') (NMNO(1,noopti,nx),noopti=1,NTOPTI)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Solution    '',5D12.4/,(13X,5D12.4))')
     '      (XC(noopti),noopti=1,NTOPTI)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

        ! sensitivity wrt nodal parameters or pressure bcs
        IF(TYPE(1:5).EQ.'NODAL'.OR.TYPE(1:8).EQ.'PRESSURE') THEN
          L_SENSIT=NOPM*NJM*NPM
          L_XC2=NOPM
          CALL ALLOCATE_MEMORY(L_SENSIT,1,DPTYPE,SENSIT_PTR,MEM_INIT,
     '      ERROR,*9999)
          CALL ALLOCATE_MEMORY(L_XC2,1,DPTYPE,XC2_PTR,MEM_INIT,
     '      ERROR,*9999)

          CALL EVSENS_CALC(ISTATE,IUSER,IIY,IY,NEELEM,NMNO,NPNODE,nr,
     '      NVJP,NYNP,PAOPTI,PE,PF,PMIN,PMAX,R,RESID,RESJAC,
     '      %VAL(SENSIT_PTR),USER,XC,%VAL(XC2_PTR),XP,YP,TYPE,CONSTR,
     '      DEFORMED,FIX,FIXP,OPFILE,OUTPTPRN,ERROR,*9999)

          CALL FREE_MEMORY(XC2_PTR,ERROR,*9999)
          CALL FREE_MEMORY(SENSIT_PTR,ERROR,*9999)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVSENS')
      RETURN
 9999 CALL ERRORS('EVSENS',ERROR)
      CALL EXITS('EVSENS')
      RETURN 1
      END


