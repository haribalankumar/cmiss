      SUBROUTINE DRISOC(ISEG,ISISOC,ITHRES,
     '  NBJ,NEELEM,NGAP,NKJE,NPF,NPNE,NRE,NVJE,
     '  SE,THRES,XA,XE,XP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRISOC
C###  Description:
C###    DRISOC draws isochrones in segments ISISOC(iw).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'colo00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks001.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISISOC(NWM),ITHRES(3,NGM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NGAP(NIM,NBM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),THRES(3,NGM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,ICOL,IEND,IEND1,IFILL,INDEX,INDEX_POLYLINE,
     '  IPAT,iw,IWK(6),LINE,LINE0,LINE1,
     '  N3CO,nb,ne,NEE,ng,noelem,noisoc,noiw,NTIW
      CHARACTER CHAR*22,FILE*255
      LOGICAL CBBREV,FILEIP,FOUND

      CALL ENTERS('DRISOC',*9999)
C GMH 2/9/95 Unused      nc=1 ! Temporary MPN 12-Nov-94
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw isochrones<;FILENAME>
C###  Description:
C###    Draw isochrones.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRISOC',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(NTCOQU(noco).EQ.1) THEN
          FILEIP=.TRUE.
          FILE=COQU(noco,1)
        ELSE
          FILEIP=.FALSE.
        ENDIF

C ***   Find range of excitation variable
        ZMINI=-2.0D0
        ZMAXI= 2.0D0
        ZDIFF=ZMAXI-ZMINI
        IF(DABS(ZDIFF).LT.1.0D-6) ZDIFF=1.0D0
        IF(DOP) THEN
          WRITE(OP_STRING,'('' ZMINI='',E12.3,'' ZMAXI='',E12.3)')
     '      ZMINI,ZMAXI
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(NINDICES.EQ.2) THEN      !monochrome
          COLOUR_WS=.FALSE.
          DO noelem=1,NEELEM(0,1)
            ne=NEELEM(noelem,1)
            nb=NBJ(1,ne)
            DO ng=1,NGT(nb)
              IFILL=ng+(ne-1)*NGT(nb) !is fill area index
              IPAT=16-NINT((THRES(2,ng,ne)-ZMINI)/ZDIFF*15.0D0)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' IFILL='',I4,'' IPAT='',I4)')
     '            IFILL,IPAT
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
c             CALL SET_FILL_REP(1,IFILL,'PATTERN',12+(17-IPAT),1,ERROR,
c    '          *9999)
            ENDDO
          ENDDO
        ELSE IF(NINDICES.GT.2) THEN !colour
          COLOUR_WS=.TRUE.
          DO noelem=1,NEELEM(0,1)
            ne=NEELEM(noelem,1)
            nb=NBJ(1,ne)
            DO ng=1,NGT(nb)
              IFILL=ng+(ne-1)*NGT(nb) !is fill area index
              ICOL=249-NINT((THRES(2,ng,ne)-ZMINI)/ZDIFF*232.0D0)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' IFILL='',I4,'' ICOL='',I4)')
     '            IFILL,ICOL
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
c             CALL SET_FILL_REP(1,IFILL,'SOLID',1,ICOL,ERROR,*9999)
            ENDDO
          ENDDO
        ENDIF

        IF(FILEIP) THEN   !!!!Note needs changing for type=area
          IOTYPE=2
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.history','OLD',
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          LINE0=1
          LINE1=10000
          DO noisoc=1,NTISOC
C            CHAR=' Isochrone number '//CFROMI(noisoc,'(I4)')
            WRITE(CHAR,'('' Isochrone number '',I4)') noisoc
            CALL FIND(IFILE,CHAR,LINE0,LINE1,LINE,FOUND,ERROR,*9999)
            LINE0=LINE
            IF(FOUND) THEN
              DO noelem=1,NEELEM(0,1)
                ne=NEELEM(noelem,1)
                nb=NBJ(1,ne)
                LINE=LINE+1
                READ(UNIT=IFILE,FMT='('' ne='',I4,'' ITHRES: '',27I1)',
     '            REC=LIne) NEE,(ITHRES(1,ng,ne),ng=1,NGT(nb))
C GMH 2/9/95 Unused - is NEE used anywhere?
                LINE=LINE+1
                IF(NIT(nb).EQ.2) THEN
                  READ(UNIT=IFILE,FMT='(9E12.3)',REC=LIne)
     '              (THRES(1,ng,ne),ng=1,NGT(nb))
                ELSE IF(NIT(nb).EQ.3) THEN
                  READ(UNIT=IFILE,FMT='(9E12.3)',REC=LIne)
     '              (THRES(1,ng,ne),ng=1,9)
                  LINE=LINE+1
                  READ(UNIT=IFILE,FMT='(9E12.3)',REC=LIne)
     '              (THRES(1,ng,ne),ng=10,18)
                  LINE=LINE+1
                  READ(UNIT=IFILE,FMT='(9E12.3)',REC=LIne)
     '              (THRES(1,ng,ne),ng=19,27)
                ENDIF
              ENDDO
              DO noiw=1,NTIW
                IW=IWK(noiw)
C CPB 28/3/94 Changed workstation update mode
                CALL ACWK(iw,1,ERROR,*9999)
C                CALL ACWK(iw,0,ERROR,*9999)
                CALL SGISOC(INDEX,ISEG,ISISOC(iw),iw,NBJ,
     '            NEELEM,NGAP,NKJE,NPF,NPNE,NRE,NVJE,
     '            SE,THRES,XA,XE,XP,CSEG,ERROR,*9999)
                CALL DAWK(iw,1,ERROR,*9999)
C                CALL DAWK(iw,0,ERROR,*9999)
              ENDDO
            ELSE
              WRITE(OP_STRING,*) CHAR//' not found'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
          CALL CLOSEF(IFILE,ERROR,*9999)

        ELSE IF(.NOT.FILEIP) THEN
          DO noiw=1,NTIW
            IW=IWK(noiw)
C CPB 28/3/94 Changed workstation update mode
            CALL ACWK(iw,1,ERROR,*9999)
C            CALL ACWK(iw,0,ERROR,*9999)
            CALL SGISOC(INDEX,ISEG,ISISOC(iw),iw,NBJ,
     '        NEELEM,NGAP,NKJE,NPF,NPNE,NRE,NVJE,
     '        SE,THRES,XA,XE,XP,CSEG,ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
C            CALL DAWK(iw,0,ERROR,*9999)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('DRISOC')
      RETURN
 9999 CALL ERRORS('DRISOC',ERROR)
      CALL EXITS('DRISOC')
      RETURN 1
      END

