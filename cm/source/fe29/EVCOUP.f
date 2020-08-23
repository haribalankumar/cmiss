      SUBROUTINE EVCOUP(NPNODE,NRLIST,NYNP,STRING,XP,YP,ERROR,*)
C#### Subroutine: EVCOUP
C###  Description:
C###    EVCOUP evaluates a coupling model.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
!      LOGICAL
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,COUP_TYPE,COUP_SUBTYPE,COUP_INSTANCE,
     '  COUPNO,N3CO,coup_no,nh,nj,nonode,nonrlist,noq,np,nr,nx,ny
      REAL*8 RVALUE(99)
      LOGICAL ALL_REGIONS,CBBREV,FOUND

      CALL ENTERS('EVCOUP',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate coupling model
C###  Parameter:        <type #[1]>
C###    Specify the coupling model type.
C###  Parameter:        <subtype #[1]>
C###    Specify the coupling model subtype.
C###  Parameter:        <instance #[1]>
C###    Specify the coupling model instance.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    evaluates a coupling model

        OP_STRING(1)=BLANK(1:15)//'<type #[1]>'
        OP_STRING(2)=BLANK(1:15)//'<subtype #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<instance #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
      ELSE
        IF(CBBREV(CO,'TYPE',2,noco+1,NTCO,N3CO)) THEN
          COUP_TYPE=IFROMC(CO(N3CO+1))
        ELSE
          COUP_TYPE=1
        ENDIF
        IF(CBBREV(CO,'SUBTYPE',2,noco+1,NTCO,N3CO)) THEN
          COUP_SUBTYPE=IFROMC(CO(N3CO+1))
        ELSE
          COUP_SUBTYPE=1
        ENDIF
        IF(CBBREV(CO,'INSTANCE',2,noco+1,NTCO,N3CO)) THEN
          COUP_INSTANCE=IFROMC(CO(N3CO+1))
        ELSE
          COUP_INSTANCE=1
        ENDIF
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
      ENDIF

      FOUND=.FALSE.
      IF(DOP)THEN
        write(*,*) 'Coupling to evaluate:',COUP_TYPE,COUP_SUBTYPE,
     '    COUP_INSTANCE
        write(*,*) 'Coupling List :'
        DO coup_no=1,COUP_CL(0,1,0)
          write(*,*) COUP_CL(coup_no,1,0),COUP_CL(coup_no,2,0),
     '      COUP_CL(coup_no,3,0)
        ENDDO
      ENDIF
      DO coup_no=1,COUP_CL(0,1,0)
        IF((COUP_CL(coup_no,1,0).EQ.COUP_TYPE).AND.
     '    (COUP_CL(coup_no,2,0).EQ.COUP_SUBTYPE).AND.
     '    (COUP_CL(coup_no,3,0).EQ.COUP_INSTANCE))THEN
          CALL ASSERT(.NOT.FOUND,'>>Duplicate coupling model found',
     '      ERROR,*9999)
          FOUND=.TRUE.
          COUPNO=coup_no
        ENDIF
      ENDDO
      CALL ASSERT(FOUND,'>>Coupling model not defined',
     '  ERROR,*9999)
      DO nonrlist=1,NRLIST(0)
        nr=NRLIST(nonrlist)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          !get rvalues
          DO noq=1,COUP_CL(COUPNO,5,0)
            IF(COUP_CL(COUPNO,3,noq).EQ.1)THEN
              IF(COUP_CL(COUPNO,4,noq).EQ.1)THEN
                RVALUE(noq)=CCR(COUP_CL(COUPNO,5,noq))
              ELSEIF(COUP_CL(COUPNO,4,noq).EQ.2)THEN
                nj=NJ_LOC(NJL_FIEL,COUP_CL(COUPNO,5,noq),nr)
                RVALUE(noq)=XP(1,1,nj,np)
              ELSEIF(COUP_CL(COUPNO,4,noq).EQ.3)THEN
                nx=COUP_CL(COUPNO,6,noq)
                nh=NH_LOC(COUP_CL(COUPNO,7,noq),nx)
                ny=NYNP(1,1,nh,np,0,1,nr)
                RVALUE(noq)=YP(ny,1,nx)
              ENDIF
            ENDIF
          ENDDO

          !call coupling model equation
          CALL COUP_EVALUATE(COUP_TYPE,COUP_SUBTYPE,RVALUE,
     '      ERROR,*9999)
          !update output/s
          DO noq=1,COUP_CL(COUPNO,5,0)
            IF(COUP_CL(COUPNO,3,noq).EQ.0)THEN
              nj=NJ_LOC(NJL_FIEL,COUP_CL(COUPNO,5,noq),nr)
              XP(1,1,nj,np)=RVALUE(noq)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('EVCOUP')
      RETURN
 9999 CALL ERRORS('EVCOUP',ERROR)
      CALL EXITS('EVCOUP')
      RETURN 1
      END


