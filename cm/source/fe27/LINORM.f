      SUBROUTINE LINORM(NEELEM,NRLIST,NW,STRING,ERROR,*)

C#### Subroutine: LINORM
C###  Description:
C###    LINORM lists bem normals which will be standard unless
C###    they have been changed by the define normal command.
C***  Created by Martin Buist April 1998

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NRLIST(0:NRM),NW(NEM,3,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nx
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,OPFILE

      CALL ENTERS('LINORM',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C-------------------------------------------------------------------

C#### Command: FEM list normals<;FILENAME>
C###  Description:
C###    Lists normal reversals (for BEM problems) to screen or to a
C###    file where a filename is specified.
C###  Parameter:   <region (all/#s)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:IEND)//'<region (all/#s)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C-------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LINORM',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opnorm','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        nx=1 ! MPN 14Jun2000  may need generalising?

        CALL OPNORM(NEELEM,NRLIST,NW(1,1,nx),ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LINORM')
      RETURN
 9999 CALL ERRORS('LINORM',ERROR)
      CALL EXITS('LINORM')
      RETURN 1
      END


