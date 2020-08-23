      SUBROUTINE DESING(STRING,ERROR,*)

C#### Subroutine: DESING
C###  Description:
C###    <HTML>
C###    DESING defines the node and element number at which a physical
C###    singularity occurs.  Examples are at the corner of surface
C###    cavities (governed by the modified Helmholtz equation) and
C###    stress singularities in some corners. This information is stored
C###    in the NPSING and NESING arrays:
C###    <UL>
C###    <LI>NPSING(i) gives the node number np of the ith singularity
C###    <LI>NPSING(0) stores the total number of NPSING entries
C###    <LI>NESING(i) gives the element number in which the singularity
C###        occurs
C###    <LI>NESING(0) stores the total number of NESING entries
C###    </UL>
C###    Currently only one node and one element per singularity.
C###    If a singularity occurs at a corner (which is usual) then it
C###    is possible to have a physical singularity in more than one
C###    element.
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'sing00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DESING',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define singularity;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines the node and element number at which a physical
C###    singularity occurs for boundary element problems. The
C###    singularity properties are read from or written to the file
C###    FILENAME (with extension .ipsing) in the directory specified
C###    by PATH.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DESING',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(NPT(1).GT.0,'>>no nodes are defined',
     '    ERROR,*9999)
        CALL ASSERT(NET(1).GT.0,'>>no elements are defined',
     '    ERROR,*9999)

C For the case when some singularity positions have already been defined
        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'sing',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPSING(ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
        IF(DOP) THEN
          DO i=0,NPSING(0)
            WRITE(OP_STRING,*)' NPSING(',I,')=',NPSING(i)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)' NESING(',I,')=',NESING(i)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

      ENDIF

      CALL EXITS('DESING')
      RETURN
 9999 CALL ERRORS('DESING',ERROR)
      CALL EXITS('DESING')
      RETURN 1
      END


