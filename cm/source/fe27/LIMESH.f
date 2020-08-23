      SUBROUTINE LIMESH(NBJ,NEELEM,NELIST,NENP,NORD,NPLIST,NPNE,NPNODE,
     &  NRLIST,NVJE,NXI,NXLIST,NYNP,CE,CP,XP,YP,STRING,ERROR,*)

C#### Subroutine: LIMESH
C###  Description:
C###    LIMESH list specialized mesh parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NORD(5,NE_R_M),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     &  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NXLIST(0:NXM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CP(NMM,NPM,NXM),XP(NKM,NVM,NJM,NPM),
     &  YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nr,nx,N3CO
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,CBBREV,OPFILE,LISOLN

      CALL ENTERS('LIMESH',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list mesh<;FILENAME>
C###  Description:
C###    List specialised mesh information. Mesh information is
C###    written to the file FILENAME (with extension .opmesh).
C###  Parameter:    <for REGION[1]>
C###    Specify which mesh region is to be listed.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<for REGION[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list mesh<;FILENAME> solution
C###  Description:
C###    List specialised mesh information. Mesh solution information is
C###    written to the file FILENAME (with extension .opmesh).
C###  Parameter:    <nodes #s,all [all]>
C###    Specify which nodes are to be listed.
C###  Parameter:    <region # [1]>
C###    Specify which mesh region is to be listed.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME solution>'
        OP_STRING(2)=BLANK(1:15)//'<nodes #s/all [all]>'
        OP_STRING(3)=BLANK(1:15)//'<region # [1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIMESH',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opmesh','NEW',
     '      'SEQUEN','FORMATTED',160,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

C MHT 25-02-03 parse regions to specify output region        
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     &    ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &    ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nx=NXLIST(1)
        
        IF(CBBREV(CO,'SOLUTION',3,noco+1,NTCO,N3CO)) THEN
          LISOLN=.TRUE.
        ELSE
          LISOLN=.FALSE.
        ENDIF
        IF(CBBREV(CO,'FIELDS',2,noco+1,NTCO,N3CO)) THEN
          WRITE(OP_STRING,
     &      '('' RADIUS     defined in field number:'',I6)')
     &      NJ_TYPE(nj_radius,2)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     &      '('' ALVEOLI    defined in field number:'',I6)')
     &      NJ_TYPE(nj_alveoli,2)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE
          nr=NRLIST(1)
          CALL OPMESH(NBJ,NEELEM(0,nr),NELIST,NENP(1,0,nr),NORD,NPLIST,
     &      NPNE,nr,NVJE,nx,NXI,NYNP(1,1,1,1,0,1,nr),CE(1,1,nx),
     &      CP(1,1,1),XP,YP(1,1,nx),LISOLN,ERROR,*9999)
        ENDIF
        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIMESH')
      RETURN
 9999 CALL ERRORS('LIMESH',ERROR)
      CALL EXITS('LIMESH')
      RETURN 1
      END


