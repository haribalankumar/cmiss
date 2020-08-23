      SUBROUTINE DEFIBR(NKJ,NPNODE,
     '  NRLIST,NVJP,XP,STRING,ERROR,*)

C#### Subroutine: DEFIBR
C###  Description:
C###    DEFIBRE defines element fibre orientations.  Constant vectors
C###    of Xi-coordinate length DXIF are drawn on plane Xi(3)=XIF at
C###    increments of DXI1 and DXI2 in Xi(1) and Xi(2) direction.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVJP(NJM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,
     '  IEND2,IPFILE,n3co,no_nrlist,nr
      CHARACTER FILE*(MXCH),STATUS*3,TYPE*9
      LOGICAL ABBREV,ALL_REGIONS,CALCU,CBBREV,FILIO,
     '  FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEFIBR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define fibre;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    The fibre field is read from or written to a file
C###    FILENAME.ipfibr.  This command defines a fibre field, the
C###    values being specified at the nodes.  The fibre field can be
C###    carried either as an angle eta (in degrees or radians), or
C###    preferably as a two variable field defined as cos(2*eta),
C###    sin(2*eta).
C###  Parameter:      <from geometry>
C###    Specifies the default number of versions to match the first
C###    geometric variable
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The all value
C###    specifies all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<from geometry>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEFIBR',ERROR,*9999)
      ELSE
C        IPFILE=2 !is input file version number on 10-Mar-1997
C!!! CS Changed the default behavior of sheet angles, see
C!!! MAT_VEC_ROTATE
        IPFILE=3 !is input file version number on 29-Jan-2001 CS

        CALL ASSERT(CALL_NODE,'>>Define nodes first',ERROR,*9999)

        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,n3co)) THEN
          IF(ABBREV(CO(n3co+1),'GEOMETRY',1)) THEN
            TYPE='GEOMETRY'
          ENDIF
        ENDIF

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'fibr',
     '        STATUS,ERR,ERROR,*9999)
            IF(ERR.EQ.-1) THEN
              CALL ASSERT(.FALSE.,'>>ERROR: File update/s necessary'
     '          ,ERROR,*9999)
            ENDIF
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL IPFIBR(NKJ,NPNODE,nr,NVJP,XP,TYPE,ERROR,*9999)
            ENDDO
            CALL_FIBR=.TRUE.
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('DEFIBR')
      RETURN
 9999 CALL ERRORS('DEFIBR',ERROR)
      CALL EXITS('DEFIBR')
      RETURN 1
      END


