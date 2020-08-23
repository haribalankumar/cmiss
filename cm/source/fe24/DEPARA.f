      SUBROUTINE DEPARA(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,ISIZE_MFI,
     '  ISIZE_TBH,NDIPOLES,NDLT,NEELEM,NEL,NENP,NFFACE,NLLINE,
     '  NONY,NPNODE,NQET,NQXI,NVHP,NVJP,NYNO,STRING,YP,ERROR,*)

C#### Subroutine: DEPARA
C###  Description:
C###    DEPARA defines array dimension parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),
     '  ISIZE_MFI(3,NSSM),ISIZE_TBH(2),NDLT(NEM),
     '  NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),NFFACE(0:NF_R_M,NRM),
     '  NLLINE(0:NL_R_M,0:NRM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),NQXI(0:NIM,NQSCM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM)
      REAL*8 YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,
     '  iw,MAXNITW,N3CO
      CHARACTER STATUS*3,FILE*(MXCH)
      LOGICAL ALL_REGIONS,BACKUP,CALCU,CBBREV,FILIO,FIRST_TIME,
     '  GENER,MINIMAL,MOUSE,REALLOCATE

      CALL ENTERS('DEPARA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define parameters;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines array/problem sizes/parameters. Parameters can be read
C###    from or written to the file FILENAME.ippara, with $current
C###    specifing the current default file.
C###  Parameter:      <minimal>
C###    Using the parameter minimal,after fully defining a problem,
C###    will write out the minimum parameter set required to solve
C###    the problem.
C###  Parameter:      <reallocate/no_reallocate>
C###    Specifies wiether arrays are reallocated after defining a new
C###    set of parameters. The default is to reallocate.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<minimal>'
        OP_STRING(3)=BLANK(1:15)//'<reallocate'
     '    //'/no_reallocate[reallocate]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEPARA',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS('LPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(CBBREV(CO,'MINIMAL',1,noco+1,NTCO,N3CO)) THEN
          MINIMAL=.TRUE.
        ELSE
          MINIMAL=.FALSE.
        ENDIF
        IF(CBBREV(CO,'NO_REALLOCATE',1,noco+1,NTCO,N3CO)) THEN
          REALLOCATE=.FALSE.
        ELSE
          REALLOCATE=.TRUE.
        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'para',
     '        STATUS,ERR,ERROR,*9999)
            CALL IPPARA(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '        ISIZE_MFI,ISIZE_TBH,MINIMAL,NDIPOLES,NDLT,NEELEM,NEL,NENP,
     '        NFFACE,NLLINE,NONY,NPNODE,NQET,NQXI,NVHP,NVJP,NYNO,YP,
     '        ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2

C reallocates after defining a new set of parameters

            IF (REALLOCATE.AND.(IOTYPE.NE.3)) THEN
              REALLOCATE_FEM=.TRUE.
              REALLOCATE_SYNTAX=.TRUE.
              MAXNITW=6 ! max # windows that can be parsed
              DO iw=1,MAXNITW
                IF(IWKG(iw).EQ.1) THEN
                  IWKG(iw)=0
                  CALL CLOSE_WS(iw,ERROR,*9999)
                ENDIF
              ENDDO !iw
            ENDIF
          ENDDO
        ENDIF !filio
      ENDIF

      CALL EXITS('DEPARA')
      RETURN
 9999 CALL ERRORS('DEPARA',ERROR)
      CALL EXITS('DEPARA')
      RETURN 1
      END


C Subroutine: DEPOTE -- defines potentials.
C    *** Archived 4 November 1998 ***


