      SUBROUTINE LIXI(NEP,NPNODE,NRLIST,LD,XID,XIP,STRING,ERROR,*)

C#### Subroutine: LIXI
C###  Description:
C###    LIXI lists Xi coordinates of data (XID) or nodes (XIP)
C###    and their associated element coordinates, LD for data points
C###     and NEP for nodes.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER LD(NDM),NEP(NPM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM)
      REAL*8 XID(NIM,NDM),XIP(NIM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,no_nrlist,nr
      CHARACTER FILE*100,XI_TYPE*4
      LOGICAL ALL_REGIONS,OPFILE,CBBREV

      CALL ENTERS('LIXI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list xi<;FILENAME>
C###  Description:
C###    Lists Xi coordinates of nodes or data to screen or
C###    FILENAME.opxi in the directory
C###    specified by PATH with $current specifing the current default
C###    file if qualifier FILENAME is present.
C###  Parameter:      <(data/node)[node]>
C###    Specifies weither it is the xi position of nodes or data that
C###    is being listed.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=STRING(1:IEND)//' (data/node)[data]'

C LKC 25-MAY-1998 XID has no region associated with it
C GBS 26-Oct-2000 Replaced: wanting to list for nodes in a region
C###  Parameter:    <region (#s/all)[1]>
        OP_STRING(3)=BLANK(1:15)//'  <region (#s/all[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIXI',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opxi','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

C***    Set xi type
        IF(CBBREV(CO,'NODE',3,noco+1,NTCO,N3CO)) THEN
          XI_TYPE='NODE'
        ELSE
          XI_TYPE='DATA'
        ENDIF


C LKC 25-MAY-1998 This is nothing to do with xi -> rewrite
C GBS 26-Oct-2000 Replaced, with new OPXI call

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
C          CALL OPXI(NPNODE,nr,XP,ERROR,*9999)
          CALL OPXI(NEP,NPNODE(0,nr),XI_TYPE,LD,XID,XIP,ERROR,*9999)
        ENDDO
C        CALL OPXI(NEP,NPNODE,XI_TYPE,LD,XID,XIP,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIXI')
      RETURN
 9999 CALL ERRORS('LIXI',ERROR)
      CALL EXITS('LIXI')
      RETURN 1
      END

