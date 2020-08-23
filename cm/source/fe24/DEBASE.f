      SUBROUTINE DEBASE(IBT,IDO,INP,NAN,NDET,NGAP,NKB,NNB,NSB,
     '  DET,PG,WG,XIG,STRING,ERROR,*)

C#### Subroutine: DEBASE
C###  Description:
C###    DEBASE defines basis parameters.

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
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NDET(NBFM,0:NNM),NGAP(NIM,NBM),
     '  NKB(2,2,2,NNM,NBFM),NNB(4,4,4,NBFM),NSB(NKM,NNM,NBFM)
      REAL*8 DET(NBFM,0:NNM,NGM,6),PG(NSM,NUM,NGM,NBM),
     '  WG(NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,nb
c     INTEGER ID_DEVICE,INPUT_STATUS,noch
c     REAL R4DATA(2)
      CHARACTER FILE*(MXCH),STATUS*3
c     CHARACTER SDATA*10
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,
     '  GENER,MOUSE

      CALL ENTERS('DEBASE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define bases;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Basis functions are read from or written to a file
C###    FILENAME.ipbase.  This command defines the basis functions to
C###    be used in defining both the geometric variables and the
C###    dependent variables of the problem.  The following types of
C###    basis function may be defined:
C###    (0) Auxiliary basis only
C###    (1) Lagrange/Hermite tensor product
C###    (2) Simplex/Serendipity/Lagrange
C###    (3) B-spline tensor product
C###    (4) Blending function
C###    (5) Chebyshev polynomial
C###    (6) Singular
C###    (7) Non-conforming
C###    Types 1..7 may also be defined with an auxiliary basis.  For
C###    each Xi direction the number of Gauss points (up to 7) is
C###    entered, determining the accuracy of Gauss-Legendre quadrature
C###    in that direction - n Gauss points integrate a (2n-1)th degree
C###    polynomial.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C       OP_STRING(1)=STRING(1:IEND)//';m'
C       CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEBASE',ERROR,*9999)
      ELSE
        IPFILE=2 !is input file version number on 14-Oct-95
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'base',
     '        STATUS,ERR,ERROR,*9999)
            CALL IPBASE(IBT,IDO,INP,NAN,0,NDET,NGAP,NKB,NNB,NSB,
     '        DET,PG,WG,XIG,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            CALL ASSERT(NBT.LE.NBM,'>>NBM too small',ERROR,*9999)
            CALL ASSERT(NBFT.LE.NBFM,'>>NBFM too small',ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
          DO nb=1,NBFT
C LKC 8-APR-1999 splitting the assert up
C            CALL ASSERT(NIT(nb).LE.NIM.
C     '        AND.NKT(0,nb).LE.NKM.AND.NST(nb).LE.NSM,
C     '        '>>NIM or NKM or NSM too small',ERROR,*9999)
            CALL ASSERT(NIT(nb).LE.NIM,'>> Increase NIM',ERROR,*9999)
            CALL ASSERT(NKT(0,nb).LE.NKM,'>> Increase NKM',ERROR,*9999)
            CALL ASSERT(NST(nb).LE.NSM,'>> Increase NSM',ERROR,*9999)
          ENDDO

          DO nb=1,NBT
            CALL ASSERT(NGT(nb).LE.NGM,'>>NGM too small',ERROR,*9999)
          ENDDO
        ENDIF
        CALL_BASE=.TRUE.
      ENDIF

      CALL EXITS('DEBASE')
      RETURN
 9999 CALL ERRORS('DEBASE',ERROR)
      CLOSE(UNIT=IFILE)
      CALL EXITS('DEBASE')
      RETURN 1
      END


