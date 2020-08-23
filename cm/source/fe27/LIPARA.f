      SUBROUTINE LIPARA(NEELEM,NENP,NFFACE,NLLINE,NPNODE,NXLIST,NYNR,
     '  STRING,ERROR,*)

C#### Subroutine: LIPARA
C###  Description:
C###    LIPARA lists CMISS parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NFFACE(0:NF_R_M,NRM),NLLINE(0:NL_R_M,0:NRM),
     '  NPNODE(0:NP_R_M,0:NRM),NXLIST(0:NXM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nx,nxc
      CHARACTER CLOWER*3,TYPE*10,PARAM*3,FILE*100
      LOGICAL ABBREV,CBBREV,NX_DEFINED,OPFILE

      CALL ENTERS('LIPARA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list parameters<;FILENAME>
C###  Parameter:    <(n*/all)[all]>
C###   Indicates which parameters to list. The n* option provides
C###   for the listing of a particular category of parameters,
C###   substitute the letter of the category to be listed for the *.
C###   eg. To list element parameters use ne
C###  Parameter:    <(dimensions/itype/jtype/ktype/nh_loc/nj_loc/nx_list
C###  Parameter:       /calls/offsets/use/versions
C###  Parameter:       /dipoles/timevars)[dimensions]>
C###   Specify the parameter type to list.
C###  Parameter:    <using (fit/optimisation/solve)[solve]>
C###   Specify which problem type to list the parameters for
C###  Parameter:    <class #[1]>
C###    Specify the class number to list the parameters for.
C###  Parameter:    <force_nx #>
C###    Specify the nx number to list the parameters for
C###  Description:
C###    List parameters lists the parameter set that cmiss is using.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<(n*/all)[all]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<(dimensions/itype/jtype/ktype/nh_loc/nj_loc/nx_list'
        OP_STRING(4)=BLANK(1:17)
     '    //'/calls/offsets/use/versions/dipoles)[dimensions]>'
        OP_STRING(5)=BLANK(1:15)
     '    //'<using (fit/optimisation/solve)[solve]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<force_nx #>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIPARA',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.oppara','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'ITYPE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='ITYPE'
        ELSE IF(CBBREV(CO,'JTYPE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='JTYPE'
        ELSE IF(CBBREV(CO,'KTYPE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='KTYPE'
        ELSE IF(CBBREV(CO,'NH_LOC',3,noco+1,NTCO,N3CO)) THEN
          TYPE='NH_LOC'
        ELSE IF(CBBREV(CO,'NJ_LOC',3,noco+1,NTCO,N3CO)) THEN
          TYPE='NJ_LOC'
        ELSE IF(CBBREV(CO,'NX_LIST',3,noco+1,NTCO,N3CO)) THEN
          TYPE='NX_LIST'
        ELSE IF(CBBREV(CO,'CALLS',2,noco+1,NTCO,N3CO)) THEN
          TYPE='CALLS'
        ELSE IF(CBBREV(CO,'OFFSETS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='OFFSETS'
        ELSE IF(CBBREV(CO,'USE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='USE'
        ELSE IF(CBBREV(CO,'VERSIONS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='VERSIONS'
        ELSE IF(CLOWER(CO(noco+1)(1:1)).EQ.'n') THEN
          TYPE='DIMENSIONS'
          PARAM=CLOWER(CO(noco+1)(1:2))
        ELSE IF(CBBREV(CO,'DIMENSIONS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='DIMENSIONS'
          PARAM='all'
        ELSE IF(CBBREV(CO,'DIPOLES',2,noco+1,NTCO,N3CO)) THEN
          TYPE='DIPOLES'
        ELSE IF(CBBREV(CO,'TIMEVARS',4,noco+1,NTCO,N3CO)) THEN
          TYPE='TIMEVARS'
        ELSE
          TYPE='DIMENSIONS'
          PARAM='all'
        ENDIF

        NX_DEFINED=.FALSE.
        IF (CBBREV(CO,'FORCE_NX',8,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(N3CO+1.LE.NTCO,
     '      '>>must enter value for nx',ERROR,*9999)
          nx=IFROMC(CO(N3CO+1))
          CALL ASSERT((nx.GT.0.AND.nx.LE.NXM),
     '      '>>nx outside permissible limits',ERROR,*9999)
          NX_DEFINED=.TRUE.
        ENDIF

        IF(TYPE(1:7).NE.'NX_LIST'.AND.
     '    (CALL_EQUA.OR.CALL_FIT.OR.CALL_OPTI)) THEN
          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'SOLVE',2)) THEN
               CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ELSE IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
               CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this fit class',ERROR,*9999)
            ELSE IF(ABBREV(CO(N3CO+1),'OPTIMISE',2)) THEN
               CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this optimisation class',
     '          ERROR,*9999)
            ELSE
               CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ENDIF
          ELSE
             CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ENDIF
          NX_DEFINED=.TRUE.
        ELSE
          IF(TYPE(1:5).EQ.'ITYPE') THEN
            nx=1 !is needed for itype variables
          ENDIF
        ENDIF

        CALL OPPARA(TYPE,PARAM,NEELEM,NENP,NFFACE,NLLINE,NPNODE,nx,NYNR,
     '    NX_DEFINED,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF

      ENDIF

      CALL EXITS('LIPARA')
      RETURN
 9999 CALL ERRORS('LIPARA',ERROR)
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9998)
        IOFI=IOOP
      ENDIF
 9998 CALL EXITS('LIPARA')
      RETURN 1
      END


