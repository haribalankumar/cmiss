      SUBROUTINE OPVOLU(NBJ,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,NW,
     '  PG,RG,SE,WG,XA,XE,XG,XN,XP,ERROR,*)

C#### Subroutine: OPVOLU
C###  Description:
C###    OPVOLU outputs the volume enclosed by a selected group of
C###    elements. Currently only implemented for closed groups of
C###    2d elements in 3d space.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
C MHT 24-03-00 common blocks not used
C      INCLUDE 'cmiss$reference:file00.cmn'
C      INCLUDE 'cmiss$reference:inout00.cmn'
C      INCLUDE 'cmiss$reference:ityp00.cmn'
C      INCLUDE 'cmiss$reference:jtyp00.cmn'
C      INCLUDE 'cmiss$reference:loc00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3)
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  VOL,WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nr
      REAL*8 DUMMY(1)

      CALL ENTERS('OPVOLU',*9999)
      DUMMY(1)=0.0d0
      nr=NRLIST(1) !Assumes all elements in the same region.
      CALL VOLUME(NBJ,NELIST,NKJE,NPF,NPNE,nr,NVJE,NW,
     '  PG,RG,SE,VOL,WG,XA,XE,XG,XN,XP,ERROR,*9999)
      WRITE(OP_STRING,'(/'' Volume enclosed by elements is : '',E12.6)')
     '  VOL
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' Element region is: '',I2)') nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' #Elements in region '',I1,'' is '','
     '  //'I5,'', Elements are:'')') nr,NELIST(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      CALL WRITE_LONG(INTTYPE,1,1,IOFI,NELIST(0),10,10,
     '  NELIST(1),DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)


      CALL EXITS('OPVOLU')
      RETURN
 9999 CALL ERRORS('OPVOLU',ERROR)
      CALL EXITS('OPVOLU')
      RETURN 1
      END


