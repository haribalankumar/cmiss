      SUBROUTINE DEFIEL(NKJ,NEELEM,NENP,NPNODE,NRLIST,NVJP,XAB,XP,
     &  STRING,ERROR,*)

C#### Subroutine: DEFIEL
C###  Description:
C###    DEFIEL defines a field variable.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVJP(NJM,NPM)
      REAL*8 XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IFROMC,IPFILE,IWK(6),N3CO,njj,NJJ_END,NJJ_START,
     '  no_nrlist,nr,NTIW,NUM_FIELD_DEFAULT,ny
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,DEPEND,DERIVS,ELEM_FIELD,FILIO,
     '  FIRST_TIME,FLOW_NODAL,GENER,HYPOXIA,MECH,MOUSE,STRETCH,TREE,
     &  VALUES,VERSIONS

      CALL ENTERS('DEFIEL',*9999)
      
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define field;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define a scalar field. This command establishes a field
C###    variable to be carried at the nodes.  It is essential to
C###    define a field prior to fitting the field values to data.
C###    Field values are read from or written to the file FILENAME
C###    (with extension .ipfiel) in the directory specified by PATH.
C###  Parameter:      <ny=#>
C###    If ny=NUMBER is
C###    specified, the YP(NUMBER) array may be written out to a file.
C###  Parameter:      <field_variables=#>
C###    Field_variables is used to set the default number for a defualt
C###    input. i.e. this is only relevant when ;d; is used.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<ny=#>'
        OP_STRING(3)=BLANK(1:15)//'<field_variables #[3]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEFIEL',ERROR,*9999)
      ELSE
        IPFILE=3 !is input file version number on 10-Mar-1997
c        CALL ASSERT(CALL_NODE,'>>Define nodes first',ERROR,*9999)
c        CALL ASSERT(CALL_ELFD,'>>Define field elements first',
c       '    ERROR,*9999)

        CALL PARSE_QUALIFIERS('DLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(MOUSE) CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'ny',1,noco+1,NTCO,N3CO)) THEN
          DEPEND=.TRUE.
          ny=IFROMC(CO(N3CO+1))
          CALL ASSERT(NHM.EQ.NJM,'>>NHM.NE.NJM',ERROR,*9999)
          CALL ASSERT(ny.LE.16,'>>ny.gt.16',ERROR,*9999)
          ERROR='>> Check code for dependent field'
          GOTO 9999
        ELSE
          DEPEND=.FALSE.
        ENDIF

C       MHT 09Apr03 so that number of field variables can be specified
C       on command line for default input.
        IF(CBBREV(CO,'FIELD_VARIABLES',2,noco+1,NTCO,N3CO)) THEN
          NUM_FIELD_DEFAULT=IFROMC(CO(N3CO+1))
        ELSE
          NUM_FIELD_DEFAULT=NJT !previous default
        ENDIF

C       Specify whether the default number of versions is the
C       same as the versions for geometry        
        IF(CBBREV(CO,'TREE',3,noco+1,NTCO,N3CO)) THEN
          TREE=.TRUE.
        ELSE
          TREE=.FALSE.
        ENDIF
        IF(CBBREV(CO,'DERIVATIVES',3,noco+1,NTCO,N3CO)) THEN
          DERIVS=.TRUE.
        ELSE
          DERIVS=.FALSE.
        ENDIF
        IF(CBBREV(CO,'VERSIONS',3,noco+1,NTCO,N3CO)) THEN
          VERSIONS=.TRUE.
        ELSE
          VERSIONS=.FALSE.
        ENDIF
        IF(CBBREV(CO,'VALUES',3,noco+1,NTCO,N3CO)) THEN
          VALUES=.TRUE.
        ELSE
          VALUES=.FALSE.
        ENDIF
        IF(CBBREV(CO,'START',2,noco+1,NTCO,N3CO)) THEN
          NJJ_START=IFROMC(CO(N3CO+1))
        ELSE
          NJJ_START=1
        ENDIF
        IF(CBBREV(CO,'END',2,noco+1,NTCO,N3CO)) THEN
          NJJ_END=IFROMC(CO(N3CO+1))
        ELSE
          nr=NRLIST(1)
          NJJ_END=NJ_LOC(NJL_FIEL,0,nr)
        ENDIF
        IF(CBBREV(CO,'ELEMENTS',4,noco+1,NTCO,N3CO))THEN
          ELEM_FIELD=.TRUE.
        ELSE
          ELEM_FIELD=.FALSE.
        ENDIF

        MECH=.FALSE.
        IF(CBBREV(CO,'PLEURALS',4,noco+1,NTCO,N3CO)) MECH=.TRUE.
        STRETCH=.FALSE.
        IF(CBBREV(CO,'STRETCH',4,noco+1,NTCO,N3CO)) STRETCH=.TRUE.
        FLOW_NODAL=.FALSE.
        IF(CBBREV(CO,'FLOW_NODAL',4,noco+1,NTCO,N3CO)) FLOW_NODAL=.TRUE.
        HYPOXIA=.FALSE.
        IF(CBBREV(CO,'HYPOXIA',4,noco+1,NTCO,N3CO)) HYPOXIA=.TRUE.


        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'fiel',
     '        STATUS,ERR,ERROR,*9999)
! This needs tidying up ...PJH 21May93
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(DEPEND) THEN
C                nx=1 !temporary
C                JTYP11=NHP(NPNODE(1,nr),nr,nx)
C                IF(IOTYPE.EQ.4) THEN
C                  CALL YPZP(ny,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
C     '              NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
C     '              YP(1,1,nx),ZA,ZP,ERROR,*9999)
C                ENDIF
C                CALL IPFIEL(NBJ,NEELEM,NKJ,NPNODE,
C     '            nr,NVJE,NVJP,XP,ERROR,*9999)
C                IF(IOTYPE.EQ.2) THEN
C                  CALL ZPYP(ny,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
C     '              NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
C     '              YP(1,1,nx),ZA,ZP,ERROR,*9999)
C                ENDIF
C                JTYP11=0
              ELSE
                CALL IPFIEL(NEELEM,NENP,NJJ_END,
     &            NJJ_START,NKJ,NPNODE,nr,NUM_FIELD_DEFAULT,NVJP,XAB,XP,
     &            DERIVS,ELEM_FIELD,FLOW_NODAL,MECH,HYPOXIA,STRETCH,
     &            TREE,
     &            VALUES,VERSIONS,
     &            ERROR,*9999)


                  IF(CBBREV(CO,'FLOW_FIELD',4,noco+1,NTCO,N3CO))
     '            nej_flow=NEJ_LOC(IFROMC(CO(N3CO+1)),nr)
                  IF(CBBREV(CO,'RESISTANCE_FIELD',4,noco+1,NTCO,N3CO))
     '              nej_resis=NEJ_LOC(IFROMC(CO(N3CO+1)),nr)
                  IF(CBBREV(CO,'RADIUS',4,noco+1,NTCO,N3CO))
     '              nej_strtrad=NEJ_LOC(IFROMC(CO(N3CO+1)),nr)
                  IF(CBBREV(CO,'VELOCITY_FIELD',4,noco+1,NTCO,N3CO))
     '              nej_vel=NEJ_LOC(IFROMC(CO(N3CO+1)),nr)
               ENDIF !depend
            ENDDO !no_nrlist
            CALL CLOSEF(IFILE,ERROR,*9999)
          ENDDO !while first_time or backup

        ELSE IF(MOUSE) THEN
        ENDIF !FILIO

C... KSB 01/2008 To define nj_radius field pointer on command line
        IF(.NOT.ADD.AND..NOT.ELEM_FIELD) THEN
          IF(CBBREV(CO,'NJ_RADIUS',4,noco+1,NTCO,n3co)) THEN
C AJS 11/2010 Setting up nj_radius properly 
!             nj_radius=IFROMC(CO(N3CO+1)) 
            njj=IFROMC(CO(N3CO+1)) 
            nj_radius=NJ_LOC(NJL_FIEL,njj,nr)
          ELSE
!             nj_radius=4
            njj=1
            nj_radius=NJ_LOC(NJL_FIEL,njj,nr)
          ENDIF
        ENDIF


        CALL_FIEL=.TRUE.
      ENDIF

      CALL EXITS('DEFIEL')
      RETURN
 9999 CALL ERRORS('DEFIEL',ERROR)
      CALL EXITS('DEFIEL')
      RETURN 1
      END


