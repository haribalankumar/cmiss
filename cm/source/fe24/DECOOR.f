      SUBROUTINE DECOOR(NRLIST,STRING,ERROR,*)

C#### Subroutine: DECOOR
C###  Description:
C###    DECOOR sets number of coordinates (NJT) and coordinate
C###    type (ITYP10(1)):

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NRLIST(0:NRM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,
     '  IEND2,IL(2),IPFILE,nj,njj1,
     '  njj2,no_nrlist,nr,NTIL,N3CO
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ABBREV,ALL_REGIONS,CALCU,CBBREV,FILIO,FIRST_TIME,GENER,
     '  MOUSE

      CALL ENTERS('DECOOR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define coordinates;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    This command establishes the coordinate system
C###    to be used in the definition of the problem.
C###    The coordinate type is read from or written to a file
C###    FILENAME.ipcoor in the directory specified by PATH with $current
C###    specifing the current default file. The options that
C###    are important at the moment are the type of coordinates
C###    (rectangular cartesian, cylindrical polar, spherical polar,
C###    prolate spheroidal, oblate spheroidal), the number of
C###    coordinates (2 or 3 dimension), and whether the geometry is
C###    unsymmetric, cylindrically symmetric, or spherically symmetric.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define coordinates NUMBER,TYPE
C###  Description:
C###    This provides a shorthand way of defining the coordinate type.
C###    This method doesn't save the defined coordinates.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//' NUMBER,TYPE'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECOOR',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS(' DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)

C MPN 21Apr97 fixed bugs
C CPB 19/4/97 Defaulting to all regions
C        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN
            NRLIST(0)=NRT
            DO nr=1,NRT
              NRLIST(nr)=nr
            ENDDO
          ELSE
            CALL PARSIL(CO(N3CO+1),NRM,NRLIST(0),NRLIST(1),ERROR,*9999)
          ENDIF
        ELSE
          NRLIST(0)=NRT
          DO nr=1,NRT
            NRLIST(nr)=nr
          ENDDO
        ENDIF
        ALL_REGIONS=NRLIST(0).GT.1

C done below
C        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'coor',
     '        STATUS,ERR,ERROR,*9999)
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

              CALL IPCOOR(nr,ERROR,*9999)

            ENDDO
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO

        ELSE IF(.NOT.FILIO) THEN
          CALL PARSIL(CO(noco+1),2,NTIL,IL,ERROR,*9999)
          IF(NTIL.EQ.2) THEN
            NJT=IL(1)
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              ITYP10(nr)=IL(2)
              DO nj=1,NJT
                NJ_LOC(NJL_GEOM,nj,nr)=NJ
                NJ_TYPE(nj,1)=NJL_GEOM
                NJ_TYPE(nj,2)=NJ
              ENDDO
              NJ_LOC(NJL_GEOM,0,nr)=NJT
              IF(NJT.GT.NJ_LOC(NJL_GEOM,0,0)) NJ_LOC(NJL_GEOM,0,0)=NJT
              NJ_LOC(0,0,nr)=0
              DO njj1=1,3
                DO njj2=1,NJ_LOC(njj1,0,nr)
                  nj=NJ_LOC(njj1,njj2,nr)
                  IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=NJ
                ENDDO
              ENDDO
              IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0))
     '          NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
            ENDDO
          ENDIF
        ENDIF

        CALL ASSERT(NJT.LE.NJM,'>>NJM too small',ERROR,*9999)
      ENDIF

      CALL EXITS('DECOOR')
      RETURN
 9999 CALL ERRORS('DECOOR',ERROR)
      CALL EXITS('DECOOR')
      RETURN 1
      END


