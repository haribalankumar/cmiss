      SUBROUTINE LIMODA(LIST,NRLIST,NXLIST,EIGVAL,EIGVEC,STRING,ERROR,*)

C#### Subroutine: LIMODA
C###  Description:
C###    LIMODA lists modal values (eigenvalues and eigenvectors).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'eige00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER LIST(0:NLISTM),NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 EIGVAL(NTM,2),EIGVEC(NOM,NTM,2)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,no_nt,nr,nt,nx,nxc
      CHARACTER FILE*100
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,OPFILE,VALUE,VECTOR

      CALL ENTERS('LIMODA',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list modal<;FILENAME>
C###  Parameter:    <mode (#s/all)[all]>
C###    Specifies which modes are required to be listed. The 'all'
C###    command prompts all the modes obtained to be listed.
C###  Parameter:    <(both/value/vector)[both]>
C###    Specifies whether the eigenvalue or eigenvector for each mode
C###    is to be outputted. 'Both' can also be specified.
C###  Parameter:    <region #[1]>
C###    Specify the element file region numbers to be defined.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Lists modal values i.e. eigenvalues and eigenvectors.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<mode (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(both/value/vector)[both]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIMODA',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opmoda','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        CALL ASSERT(NLISTM.GE.NTM,'>>Increase NLISTM to be ge NTM',
     '    ERROR,*9999)
        IF(CBBREV(CO,'MODE',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO),'ALL',2)) THEN
            LIST(0)=NUMEIGEN
            DO nt=1,NUMEIGEN
              LIST(nt)=nt
            ENDDO !nt
          ELSE
            CALL PARSIL(CO(N3CO+1),NLISTM,LIST(0),LIST(1),ERROR,
     '        *9999)
            DO no_nt=1,LIST(0)
              nt=LIST(no_nt)
C              CALL ASSERT(nt.GT.0.AND.nr.LE.NUMEIGEN,
C     '          '>>Invalid mode number',ERROR,*9999)
              CALL ASSERT(nt.GT.0.AND.nt.LE.NUMEIGEN,
     '          '>>Invalid mode number',ERROR,*9999)
            ENDDO !no_nt (nt)
          ENDIF
        ELSE
          LIST(0)=NUMEIGEN
          DO nt=1,NUMEIGEN
            LIST(nt)=nt
          ENDDO !nt
        ENDIF
        VALUE=.TRUE.
        VECTOR=.TRUE.
        IF(CBBREV(CO,'VALUE',2,noco+1,NTCO,N3CO)) THEN
          VECTOR=.FALSE.
        ELSE IF(CBBREV(CO,'VECTOR',2,noco+1,NTCO,N3CO)) THEN
          VALUE=.FALSE.
        ENDIF

        nr=NRLIST(1)
        CALL ASSERT(ITYP5(nr,nx).EQ.3,'>>Modal analysis not defined',
     '    ERROR,*9999)

        CALL OPMODA(LIST,nr,nx,EIGVAL,EIGVEC,VALUE,VECTOR,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF

      ENDIF

      CALL EXITS('LIMODA')
      RETURN
 9999 CALL ERRORS('LIMODA',ERROR)
      CALL EXITS('LIMODA')
      RETURN 1
      END


