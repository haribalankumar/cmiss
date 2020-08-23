      SUBROUTINE LIBASE(IBT,IDO,INP,NAN,NGAP,NKEF,NNF,NNL,
     '  PG,XIG,STRING,ERROR,*)

C#### Subroutine: LIBASE
C###  Description:
C###    LIBASE lists basis function arrays.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NGAP(NIM,NBM),NKEF(0:4,16,6,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM)
      REAL*8 PG(NSM,NUM,NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nb,NFAM,nu
      CHARACTER FILE*100
      LOGICAL CBBREV,OPFILE

      CALL ENTERS('LIBASE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list bases<;FILENAME>
C###  Parameter:    <number (#/all)[all]>
C###    Lists the specified bases only.
C###  Parameter:    <full>
C###    Includes additional information in the listing.
C###  Parameter:    <derivative #[1]>
C###    List only the information associated with derivative #.
C###  Description:
C###    Lists basis functions to screen or file FILENAME.opbase if
C###    qualifier present.  A more complete listing is given if the
C###    parameter 'full' is included.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<number (#/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<full>'
        OP_STRING(4)=BLANK(1:15)//'<derivative #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list bases
C###  Parameter:    <family (#/all)[all]>
C###    Lists only the specified family.
C###  Description:
C###    Lists the defined bases. Can be used to list individual
C###    families of basis functions.

        OP_STRING(1)=STRING(1:IEND)//'<family (#/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIBASE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opbase','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'FULL',2,noco+1,NTCO,N3CO)) THEN
          JTYP8=1
          IF(CBBREV(CO,'DERIVATIVE',1,noco+1,NTCO,N3CO)) THEN
            nu=IFROMC(CO(N3CO+1))
          ELSE
            nu=1
          ENDIF
C CPB 15/7/92 ADDED BASIS FAMILY CODE
        ELSE IF(CBBREV(CO,'FAMILY',2,noco+1,NTCO,N3CO)) THEN
          JTYP8=2
        ELSE
          JTYP8=0
        ENDIF

        NFAM=0
        IF(jtyp8.NE.2.AND.CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
          nb=IFROMC(CO(N3CO+1))
          CALL ASSERT(nb.gt.0.and.nb.LE.NBT,
     '      '>>Basis number out of range',ERROR,*9999)
          CALL OPBASE1(IBT,IDO,INP,NAN,nb,NFAM,NGAP,NKEF,NNF,NNL,nu,
     '      PG,XIG,ERROR,*9999)
C CPB 15/7/92 ADDED BASIS FAMILY CODE
        ELSE IF(JTYP8.EQ.2) THEN
          IF(NTCO.GT.N3CO) NFAM=IFROMC(CO(N3CO+1))
          CALL ASSERT(NFAM.gt.0.and.NFAM.LE.NBFT,
     '      '>>Family number out of range',ERROR,*9999)
          CALL OPBASE(IBT,IDO,INP,NAN,NFAM,NGAP,NKEF,NNF,NNL,nu,PG,XIG,
     '      ERROR,*9999)
        ELSE
          CALL OPBASE(IBT,IDO,INP,NAN,NFAM,NGAP,NKEF,NNF,NNL,nu,PG,XIG,
     '      ERROR,*9999)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIBASE')
      RETURN
 9999 CALL ERRORS('LIBASE',ERROR)
      CALL EXITS('LIBASE')
      RETURN 1
      END


