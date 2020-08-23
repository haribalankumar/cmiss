      SUBROUTINE ADDFIEL(NKH,NKJ,NPNODE,NRLIST,NVHP,NVJP,XP,
     '  STRING,ERROR,*)

C#### Subroutine: ADDFIEL
C###  Description:
C###    ADDFIEL updates a specified field from a linear combination
C###    of specified fields i.e. field3=a*field1+b*field2.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
C MHT 24-03-00 common blocks not needed
C      INCLUDE 'cmiss$reference:b00.cmn'
C      INCLUDE 'cmiss$reference:call00.cmn'
C      INCLUDE 'cmiss$reference:cbfe01.cmn'
C      INCLUDE 'cmiss$reference:ktyp90.cmn'
!     Parameter List
      INTEGER NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,FIELD1,FIELD2,FIELD3,nc,nh,nj1,nj2,
     '  nj3,nk,no_nrlist,nonode,np,nr,nv
      REAL*8 A_COEFF,B_COEFF,RFROMC
      LOGICAL ABSOLUTE,ALL_REGIONS,CBBREV

      CALL ENTERS('ADDFIEL',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM add field
C###  Parameter:      <field1 FIELD_VAR#[1]>
C###    The number of the first field
C###  Parameter:      <A #[1.0]>
C###    The weighting coefficient for the first field
C###  Parameter:      <field2 FIELD_VAR#[2]>
C###    The number of the second field
C###  Parameter:      <B #[1.0]>
C###    The weighting coefficient for the second field
C###  Parameter:      <field3 FIELD_VAR#[3]>
C###    The field to update
C###  Parameter:      <absolute>
C###    The field to update is the updated with the absolute value
C###    of the linear combination of specified fields
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Updates a specified field from a linear combination
C###    of specified fields i.e. field3=a*field1+b*field2.

C**** Created by Carey Stevens January 1999

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<field1 FIELD_VAR#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<A #[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<field2 FIELD_VAR#[2]>'
        OP_STRING(5)=BLANK(1:15)//'<B #[1.0]>'
        OP_STRING(6)=BLANK(1:15)//'<field3 FIELD_VAR#[3]>'
        OP_STRING(7)=BLANK(1:15)//'<absolute>'
        OP_STRING(8)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','ADDFIEL',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'FIELD1',6,noco+1,NTCO,N3CO)) THEN
          IF(N3CO.LT.NTCO) THEN
            FIELD1=IFROMC(CO(N3CO+1))
          ELSE
            ERROR='>>No field number supplied'
            GOTO 9999
          ENDIF
        ELSE
          FIELD1=1
        ENDIF

        IF(CBBREV(CO,'A',1,noco+1,NTCO,N3CO)) THEN
          IF(N3CO.LT.NTCO) THEN
            A_COEFF=RFROMC(CO(N3CO+1))
          ELSE
            ERROR='>>No a_coefficient value supplied'
            GOTO 9999
          ENDIF
        ELSE
          A_COEFF=1.0d0
        ENDIF

        IF(CBBREV(CO,'FIELD2',6,noco+1,NTCO,N3CO)) THEN
          IF(N3CO.LT.NTCO) THEN
            FIELD2=IFROMC(CO(N3CO+1))
          ELSE
            ERROR='>>No field number supplied'
            GOTO 9999
          ENDIF
        ELSE
          FIELD2=2
        ENDIF

        IF(CBBREV(CO,'B',1,noco+1,NTCO,N3CO)) THEN
          IF(N3CO.LT.NTCO) THEN
            B_COEFF=RFROMC(CO(N3CO+1))
          ELSE
            ERROR='>>No b_coefficient value supplied'
            GOTO 9999
          ENDIF
        ELSE
          B_COEFF=1.0d0
        ENDIF

        IF(CBBREV(CO,'FIELD3',6,noco+1,NTCO,N3CO)) THEN
          IF(N3CO.LT.NTCO) THEN
            FIELD3=IFROMC(CO(N3CO+1))
          ELSE
            ERROR='>>No field number supplied'
            GOTO 9999
          ENDIF
        ELSE
          FIELD3=3
        ENDIF

        IF(CBBREV(CO,'ABSOLUTE',6,noco+1,NTCO,N3CO)) THEN
          ABSOLUTE=.TRUE.
        ELSE
          ABSOLUTE=.FALSE.
        ENDIF

        nh=1 !should be specified on command line
        nc=1 !should be specified on command line
        DO no_nrlist=1,NRLIST(0) !loop over regions
          nr=NRLIST(no_nrlist)
          CALL ASSERT(FIELD1.GT.0.AND.FIELD1.LE.NJ_LOC(NJL_FIEL,0,nr),
     '      '>>Field variable number 1 out of range',ERROR,*9999)
          CALL ASSERT(FIELD2.GT.0.AND.FIELD2.LE.NJ_LOC(NJL_FIEL,0,nr),
     '      '>>Field variable number 2 out of range',ERROR,*9999)
          CALL ASSERT(FIELD3.GT.0.AND.FIELD3.LE.NJ_LOC(NJL_FIEL,0,nr),
     '      '>>Field variable number 3 out of range',ERROR,*9999)
          nj1=NJ_LOC(NJL_FIEL,FIELD1,nr)
          nj2=NJ_LOC(NJL_FIEL,FIELD2,nr)
          nj3=NJ_LOC(NJL_FIEL,FIELD3,nr)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            NVJP(nj3,np)=NVHP(nh,np,nc,nr)
            NKJ(nj3,np)=NKH(nh,np,nc,nr)
            DO nv=1,NVJP(nj3,np)
              DO nk=1,NKJ(nj3,np)
                XP(nk,nv,nj3,np)=A_COEFF*XP(nk,nv,nj1,np)+
     '            B_COEFF*XP(nk,nv,nj2,np)
                IF(ABSOLUTE) XP(nk,nv,nj3,np)=DABS(XP(nk,nv,nj3,np))
              ENDDO !nk
            ENDDO !nv
          ENDDO !nonode
        ENDDO !no_nrlist
      ENDIF

      CALL EXITS('ADDFIEL')
      RETURN
 9999 CALL ERRORS('ADDFIEL',ERROR)
      CALL EXITS('ADDFIEL')
      RETURN 1
      END


