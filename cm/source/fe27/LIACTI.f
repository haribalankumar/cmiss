      SUBROUTINE LIACTI(NBJ,NEELEM,NRLIST,FEXT,STRING,ERROR,*)

C#### Subroutine: LIACTI
C###  Description:
C###    LIACTI lists active muscle model parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NRLIST(0:NRM)
      REAL*8 FEXT(NIFEXTM,NGM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,no_nrlist,nr
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,CBBREV,OPFILE,FULL_OUTPUT

      CALL ENTERS('LIACTI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list active<;FILENAME>
C###  Description:
C###    List active muscle model parameters to screen or FILENAME.opacti.
C###  Parameter:    <full>
C###    The FULL option also outputs the FEXT array. See: FEXT
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(3)=BLANK(1:15)//'<full>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIACTI',ERROR,*9999)
      ELSE
        CALL ASSERT(CALL_ACTI,'>>Active muscle model parameters '
     '    //'not defined',ERROR,*9999)
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opacti','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
C!!! Currently not needed since neither the active parameters or the
C!!! FEXT array are class dependent
C        CALL PARSE_CLASS(noco,NTCO,nxc,CO,ERROR,*9999)
C         CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C     '    ERROR,*9999)

        IF(CBBREV(CO,'FULL',2,noco+1,NTCO,N3CO)) THEN
          FULL_OUTPUT=.TRUE.
        ELSE
          FULL_OUTPUT=.FALSE.
        ENDIF

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL OPACTI(NBJ,NEELEM,nr,FULL_OUTPUT,FEXT,ERROR,*9999)
        ENDDO !no_nrlist (nr)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIACTI')
      RETURN
 9999 CALL ERRORS('LIACTI',ERROR)
      CALL EXITS('LIACTI')
      RETURN 1
      END


