      SUBROUTINE DEFINE_LABEL(STRING,ERROR,*)

C#### Subroutine: DEFINE_LABEL
C###  Description:
C###    DEFINE_LABEL reads in a maths file, like a cellml, and generates
C###    the code required for the maths to be used in cmiss

C Author: Duane Malcolm
C Created: 11 March 2004

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      CHARACTER STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nj,njj,nr,nxc,nx,VERBOSE
      LOGICAL CBBREV


      CALL ENTERS('DEFINE_LABEL',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define label
C###  Parameter:      NAME
C###    Specifies the label
C###  Parameter:      <geometry [x/y/z]>
C###    Specifies the geometry to label
C###  Parameter:      <field field#>
C###    Specifies the field to label
C###  Parameter:      <equation class#>
C###    Specifies the equation to label
C###  Parameter:      <verbose verbosity#>
C###    Specifies the verbosity level
C###  Parameter:    <region #[1]>
C###    Specifies the region
C###  Description:
C###    This command defines labels to objects in cmiss.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'NAME'
        OP_STRING(3)=BLANK(1:15)//'<geometry [x/y/z]>'
        OP_STRING(4)=BLANK(1:15)//'<field field#>'
        OP_STRING(5)=BLANK(1:15)//'<equation class#>'
        OP_STRING(6)=BLANK(1:15)//'<verbose verbosity#>'
        OP_STRING(7)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','DEFINE_LABEL',ERROR,*9999)
      ELSE
        
        IF(CBBREV(CO,'REGION',3,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF
        CALL ASSERT(nrt.ge.nr,' >>Region not defined',ERROR,*9999)
        CALL ASSERT(nr.GT.0,' >>Region number 0 specified',ERROR,*9999)

        IF(CBBREV(CO,'VERBOSE',3,noco+1,NTCO,N3CO)) THEN
          VERBOSE=IFROMC(CO(N3CO+1))
        ELSE
          VERBOSE=0
        ENDIF
        
        IF(CBBREV(CO,'GEOMETRY',4,noco+1,NTCO,N3CO)) THEN
          
          IF(CO(N3CO+1).EQ.'X') THEN
            njj=1
          ELSEIF(CO(N3CO+1).EQ.'Y') THEN
            njj=2
          ELSEIF(CO(N3CO+1).EQ.'Z') THEN
            njj=3
          ELSE
            CALL ASSERT(.TRUE.,' >> No geometry defined',ERROR,*9999)
          ENDIF
          
          
          
          CALL ASSERT(njj.LE.NJ_LOC(NJL_GEOM,0,nr),
     &      'Field number does not exist',ERROR,*9999)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          CALL SETUP_FIELD(CO(noco+1),nj,ERROR,*9999)
        ELSEIF(CBBREV(CO,'FIELD',4,noco+1,NTCO,N3CO)) THEN
          njj=IFROMC(CO(N3CO+1))
          CALL ASSERT(njj.LE.NJ_LOC(NJL_FIEL,0,nr),
     &      'Field number does not exist',ERROR,*9999)
          nj=NJ_LOC(NJL_FIEL,njj,nr)
          CALL SETUP_FIELD(CO(noco+1),nj,ERROR,*9999)
        ELSEIF(CBBREV(CO,'EQUATION',4,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,
     &      'Equation class number does not exist',ERROR,*9999)
          CALL ADD_EQUATION(CO(noco+1),ERROR,*9999)
          CALL SET_EQUATION_CLASS(CO(noco+1),nxc,ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL SET_EQUATION_NX(CO(noco+1),nx,ERROR,*9999)
          CALL SET_EQUATION_VERBOSE(CO(noco+1),VERBOSE,ERROR,*9999)
          CALL SETUP_EQUATION(CO(noco+1),nr,ERROR,*9999)
        ELSE
          ERROR='>> Object type to label not recognised'
          GOTO 9999
        ENDIF

      ENDIF

      CALL EXITS('DEFINE_LABEL')
      RETURN
 9999 CALL ERRORS('DEFINE_LABEL',ERROR)
      CALL EXITS('DEFINE_LABEL')
      RETURN 1
      END


