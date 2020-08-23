      SUBROUTINE LIFACE(NBJ,NBJF,NLF,NNF,NPF,NPNE,NPNF,NRLIST,DF,SF,
     '  STRING,ERROR,*)

C#### Subroutine: LIFACE
C###  Description:
C###    LIFACE lists face parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBJF(NJM,NFM),NLF(4,NFM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NRLIST(0:NRM)
      REAL*8 DF(NFM),SF(NSM,NBFM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nf
      CHARACTER FILE*100,TYPE1*12
      LOGICAL ALL_REGIONS,CBBREV,OPFILE

      CALL ENTERS('LIFACE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list faces<;FILENAME>
C###  Description:
C###    Lists faces to screen or file FILENAME.opface if qualifier
C###    present.
C###  Parameter:    <number (FACE#/all)[all]>
C###    Specifies the face number(s) or all of the faces to list
C###  Parameter:    <region (#s/all)[1]>
C###    Specifies the region of the faces to list

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<number (FACE#/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list faces<;FILENAME> groups
C###  Description:
C###    Lists face groups to the screen or to FILENAME.opface if
C###    qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> groups'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIFACE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opface','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'GROUPS',1,noco+1,NTCO,N3CO)) THEN
          TYPE1='GROUPS'
        ELSE IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
          TYPE1='ONE_FACE'
          nf=IFROMC(CO(N3CO+1))
        ELSE
          TYPE1='ALL_FACES'
        ENDIF

        IF(TYPE1(1:6).EQ.'GROUPS') THEN
          CALL OPFACEG(ERROR,*9999)

        ELSE IF(TYPE1(1:8).EQ.'ONE_FACE') THEN
          CALL OPFACE1(NBJ,NBJF(1,nf),nf,NLF(1,nf),NNF,NPF(1,nf),NPNE,
     '      NPNF,DF(nf),SF,ERROR,*9999)

        ELSE IF(TYPE1(1:9).EQ.'ALL_FACES') THEN
          CALL OPFACE(NBJ,NBJF,NLF,NNF,NPF,NPNE,NPNF,DF,SF,ERROR,*9999)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIFACE')
      RETURN
 9999 CALL ERRORS('LIFACE',ERROR)
      CALL EXITS('LIFACE')
      RETURN 1
      END


