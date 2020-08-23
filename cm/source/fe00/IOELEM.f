      SUBROUTINE IOELEM(FILEID,NBJ,NEELEM,NKJE,NPNE,NRE,NRLIST,NVJE,
     '  COMMAND,FILEFORMAT,FILENAME,SUBCOMMAND,ENDFILE,ERROR,*)

C#### Subroutine: IOELEM
C###  Description:
C###    Handles i/o of binary elements

C###  Comment: Binary Element File Format
C###  Description:
C###    <HTML>
C###     Soon elements done in binary.....
C###    <H3>
C###    <B>Note:</B></H3>
C###
C###    <H4>
C###    The binary element format is one of the first CMISS binary
C###    files that is using the version 2 file type. This means that
C###    the data is stored in a hierarchial manner, consistent with
C###    the CMISS framework of regions etc.</H4>
C###
C###    <H3>
C###    Level 0</H3>
C###    This is an integer that indicates the number of tags at level 1
C###    This should be equal to the number of regions.
C###    <BLOCKQUOTE>
C###    <H3>
C###    Level 1</H3>
C###    This tag is the regional data tag. The number of sub-tags is
C###    equal to 2.
C###    <BLOCKQUOTE>
C###    <H3>
C###    Level 2</H3>
C###    There are two tags in level 2. They are:
C###    <BLOCKQUOTE>
C###    <LI>
C###    Region data tag</LI>
C###
C###    <LI>
C###    Element data tag</LI>
C###
C###    <H4>
C###    <I>Region data tag</I></H4>
C###    The number of sub-tags should be equal to 0, as this is a leaf
C###    tag. The data is set out as follows:
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer for the region number</LI>
C###
C###    <LI>
C###    An integer for the number of variable types
C###    (Check, overwrite if new type)</LI>
C###
C###    <LI>
C###    Foreach variable type</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer for the number of variables in the type
C###    (Check, overwrite if
C###    new type)</LI>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###
C###    <H4>
C###    <I>Element data tag</I></H4>
C###    The number of sub-tags should be equal to 0, as this is a
C###     leaf tag. The data is set out as follows:
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer to specify the number of elements in the region</LI>
C###
C###    <LI>
C###    Foreach element in the region</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer for the element number</LI>
C###
C###    <LI>
C###    Foreach variable type</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    Foreach variable of the variable type</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer for the basis type number</LI>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###
C###    <LI>
C###    An integer for the number of bases defined for the element</LI>
C###
C###    <LI>
C###    Foreach basis</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer for the basis function number
C###    (Check, crash if 1&lt;=nb&lt;=nbt)</LI>
C###
C###    <LI>
C###    An integer for the number of local nodes in the basis</LI>
C###
C###    <LI>
C###    Foreach local node of the basis</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer to specify the global node number</LI>
C###
C###    <LI>
C###    Foreach variable type</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    Foreach variable of that variable type</LI>
C###
C###    <BLOCKQUOTE>
C###    <LI>
C###    An integer to specify the global version number</LI>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###    </BLOCKQUOTE>
C###
C###    <UL>
C###    <UL>
C###    <UL>&nbsp;</UL>
C###    </UL>
C###    </UL>
C###
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'binf00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER FILEID,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM)
      CHARACTER COMMAND*(*),ERROR*(*),FILEFORMAT*(*),FILENAME*(*),
     '  SUBCOMMAND*(*)
      LOGICAL ENDFILE
!     Local variables
      INTEGER CERROR(50),ERR,VERSION(3),FILETYPE,nonr,nr,nr1(1),
     '  NUM_VAR_TYPES,njtype,NJTYPES(3),no_njtype,njj,noelem,ne,nb,
     '  NUM_BAS_FUNCS,BAS_FUNCS(100),nonb,np,nv,nn,nj,i,nk,TEMP(1)
      LOGICAL ISBINFILEOPEN

      CALL ENTERS('IOELEM',*9999)

      IF(COMMAND(1:4).EQ.'READ') THEN
        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          ERROR='>>Ascii not implemented in ioelem yet'
          GOTO 9999

C LKC 27-SEP-1999 Unused at the moment
C          IF(SUBCOMMAND(1:4).EQ.'OPEN') THEN
C          ELSE IF(SUBCOMMAND(1:5).EQ.'RESET') THEN
C          ELSE
C            ERROR='>>Invalid subcommand'
C            GOTO 9999
C          ENDIF

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
          IF(SUBCOMMAND(1:4).EQ.'OPEN') THEN
            FILETYPE=5
            VERSION(1)=1
            VERSION(2)=0
            VERSION(3)=0
            NUMBERTAGS=1 !Just header tag for the moment.
            CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '        'READ','ele',FILENAME,ERROR,*9999)

          ELSEIF(SUBCOMMAND(1:8).EQ.'ELEMDATA') THEN
            IF(ISBINFILEOPEN(FILEID)) THEN

              CALL BINREADFILE(FILEID,INTTYPE,1,NRLIST(0),REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C             ..Initialise highest element number in all regions
              NET(0)=0

              DO nonr=1,NRLIST(0)

C               ..Read in the region data = subtag 1

C               ..The region number
                CALL BINREADFILE(FILEID,INTTYPE,1,nr1,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                NRLIST(nonr)=nr1(1)
                nr=NRLIST(nonr)

C               ..Read in the number of variable types
                CALL BINREADFILE(FILEID,INTTYPE,1,TEMP,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
                NUM_VAR_TYPES=TEMP(1)

C               ..Read in the variable types
                DO no_njtype=1,NUM_VAR_TYPES
                  CALL BINREADFILE(FILEID,INTTYPE,1,NJTYPES(no_njtype),
     '              REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '              ERROR,*9999)
                ENDDO
                DO no_njtype=1,NUM_VAR_TYPES
                  njtype=NJTYPES(no_njtype)
                  CALL BINREADFILE(FILEID,INTTYPE,1,
     '              NJ_LOC(njtype,0,nr),REAL4DATA,REAL8DATA,CHARDATA,
     '              LOGDATA,SINTDATA,ERROR,*9999)
                ENDDO

C               ..Read in the element data = subtag 2

C               ..Read in the number of elements
                CALL BINREADFILE(FILEID,INTTYPE,1,NEELEM(0,nr),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)

                DO noelem=1,NEELEM(0,nr)

C                 ..Read in the element number
                  CALL BINREADFILE(FILEID,INTTYPE,1,NEELEM(noelem,nr),
     '              REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '              ERROR,*9999)
                  ne=NEELEM(noelem,nr)
                  NRE(ne)=nr

                  DO no_njtype=1,NUM_VAR_TYPES
                    njtype=NJTYPES(no_njtype)
                    DO njj=1,NJ_LOC(njtype,0,nr)
                      nj=NJ_LOC(njtype,njj,nr)

C                     ..Read in the basis function number for
C                       element ne, variable nj
                      CALL BINREADFILE(FILEID,INTTYPE,1,NBJ(nj,ne),
     '                  REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                  ERROR,*9999)

C                     ..Element derivatives
                      nb=NBJ(nj,ne)
                      DO nn=1,NNT(nb)
                        DO nk=1,NKT(nn,nb)
                          NKJE(nk,nn,nj,ne)=nk
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO

C                 ..Read in the number of basis functions for
C                   element ne
                  CALL BINREADFILE(FILEID,INTTYPE,1,TEMP,
     '              REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '              ERROR,*9999)
                  NUM_BAS_FUNCS=TEMP(1)

                  CALL ASSERT(NUM_BAS_FUNCS.LE.NBT,'>> More basis '//
     '              'functions than there is set up',ERROR,*9999)

                  DO nonb=1,NUM_BAS_FUNCS

C                   ..Read in the basis function number
                    CALL BINREADFILE(FILEID,INTTYPE,1,TEMP,REAL4DATA,
     '                REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                    nb=TEMP(1)

C                   ..Read in the number of local nodes for this
C                     basis
                    CALL BINREADFILE(FILEID,INTTYPE,1,NNT(nb),
     '                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                ERROR,*9999)

                    DO nn=1,NNT(nb)

C                     ..Read in the global node number
                      CALL BINREADFILE(FILEID,INTTYPE,1,NPNE(nn,nb,ne),
     '                  REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                  ERROR,*9999)

C KAT 23Feb01: now handled by NKJE above
CC                     ..Element derivatives
C                      DO nk=1,NKT(nn,nb)
C                        NKE(nk,nn,nb,ne)=nk
C                      ENDDO

                      DO no_njtype=1,NUM_VAR_TYPES
                        njtype=NJTYPES(no_njtype)
                        DO njj=1,NJ_LOC(njtype,0,nr)
                          nj=NJ_LOC(njtype,njj,nr)

C                         ..Read in the global version number
                          CALL BINREADFILE(FILEID,INTTYPE,1,
     '                      NVJE(nn,nb,nj,ne),REAL4DATA,REAL8DATA,
     '                      CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                        ENDDO
                      ENDDO
                    ENDDO !nn
                  ENDDO !nonb
                ENDDO !ne
C               ..Set up NET
                NET(nr)=0
                DO noelem=1,NEELEM(0,nr)
                  IF(NEELEM(noelem,nr).GT.NET(nr))
     '              NET(nr)=NEELEM(noelem,nr)
                  IF(NEELEM(noelem,nr).GT.NET(0))
     '              NET(0)=NEELEM(noelem,nr)
                ENDDO
              ENDDO ! nonr
            ENDIF
          ELSE IF(SUBCOMMAND(1:5).EQ.'RESET') THEN
          ELSE
            ERROR='>>Invalid subcommand'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF
        ENDFILE=.FALSE.
      ELSE IF(COMMAND(1:5).EQ.'WRITE') THEN
        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          IF(SUBCOMMAND(1:4).EQ.'OPEN') THEN
            ERROR='>>Ascii not implemented in ioelem yet'
            GOTO 9999
          ELSE
            ERROR='>>Not implemented yet'
            GOTO 9999
          ENDIF
        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
          IF(SUBCOMMAND(1:4).EQ.'OPEN') THEN

C           ..Open the binary element file for writing

            FILETYPE=5 !Binary element file
            VERSION(1)=1
            VERSION(2)=0
            VERSION(3)=0
            CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '        'WRITE','ele',FILENAME,ERROR,*9999)

          ELSEIF(SUBCOMMAND(1:8).EQ.'ELEMDATA') THEN
            IF(ISBINFILEOPEN(FILEID)) THEN

C             ..Write out the number of regions

              CALL BINWRITEFILE(FILEID,INTTYPE,1,NRLIST(0),REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

              DO nonr=1,NRLIST(0)
                nr=NRLIST(nonr)

C               ..Write out the region data = subtag 1

C               ..The region number
                nr1(1)=nr
                CALL BINWRITEFILE(FILEID,INTTYPE,1,nr1,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C               ..Compute the number of variable types and which
C                 of these types are defined, ie geom/fibre/field
                NUM_VAR_TYPES=0
                DO njtype=1,3
                  IF(NJ_LOC(njtype,0,nr).GE.1) THEN
                    NUM_VAR_TYPES=NUM_VAR_TYPES+1
                    NJTYPES(NUM_VAR_TYPES)=njtype
                  ENDIF
                ENDDO

C               ..Write out the number of variable types
                TEMP(1)=NUM_VAR_TYPES
                CALL BINWRITEFILE(FILEID,INTTYPE,1,TEMP,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)

C               ..Write out the variable types
                DO no_njtype=1,NUM_VAR_TYPES
                  njtype=NJTYPES(no_njtype)
                  TEMP(1)=njtype
                  CALL BINWRITEFILE(FILEID,INTTYPE,1,TEMP,REAL4DATA,
     '              REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                ENDDO
                DO no_njtype=1,NUM_VAR_TYPES
                  njtype=NJTYPES(no_njtype)
                  CALL BINWRITEFILE(FILEID,INTTYPE,1,
     '              NJ_LOC(njtype,0,nr),REAL4DATA,REAL8DATA,CHARDATA,
     '              LOGDATA,SINTDATA,ERROR,*9999)
                ENDDO

C               ..Write out the element data = subtag 2

C               ..Write out the number of elements
                CALL BINWRITEFILE(FILEID,INTTYPE,1,NEELEM(0,nr),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)

                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)

C                 ..Write out the element number
                  TEMP(1)=ne
                  CALL BINWRITEFILE(FILEID,INTTYPE,1,TEMP,
     '              REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '              ERROR,*9999)

                  DO no_njtype=1,NUM_VAR_TYPES
                    njtype=NJTYPES(no_njtype)
                    DO njj=1,NJ_LOC(njtype,0,nr)
                      nj=NJ_LOC(njtype,njj,nr)

C                     ..Write out the basis function number for
C                       element ne, variable nj
                      CALL BINWRITEFILE(FILEID,INTTYPE,1,NBJ(nj,ne),
     '                  REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                  ERROR,*9999)
                    ENDDO
                  ENDDO

C                 ..Compute the number of bases for element ne
C                   and the basis function numbers
                  NUM_BAS_FUNCS=0
                  DO i=1,NUM_BAS_FUNCS
                    BAS_FUNCS(i)=0
                  ENDDO

                  DO no_njtype=1,NUM_VAR_TYPES
                    njtype=NJTYPES(no_njtype)
                    DO njj=1,NJ_LOC(njtype,0,nr)
                      nj=NJ_LOC(njtype,njj,nr)
                      nb=NBJ(nj,ne)
                      DO i=1,NUM_BAS_FUNCS
                        IF(BAS_FUNCS(i).EQ.nb) GOTO 200
                      ENDDO
                      NUM_BAS_FUNCS=NUM_BAS_FUNCS+1
                      BAS_FUNCS(NUM_BAS_FUNCS)=nb
 200                  CONTINUE
                    ENDDO
                  ENDDO

C                 ..Write out the number of basis functions for
C                   element ne
                  TEMP(1)=NUM_BAS_FUNCS
                  CALL BINWRITEFILE(FILEID,INTTYPE,1,TEMP,
     '              REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '              ERROR,*9999)

                  DO nonb=1,NUM_BAS_FUNCS
                    nb=BAS_FUNCS(nonb)

C                   ..Write out the basis function number
                    TEMP(1)=nb
                    CALL BINWRITEFILE(FILEID,INTTYPE,1,TEMP,REAL4DATA,
     '                REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C                   ..Write out the number of local nodes for this
C                     basis
                    CALL BINWRITEFILE(FILEID,INTTYPE,1,NNT(nb),
     '                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                ERROR,*9999)

                    DO nn=1,NNT(nb)

                      np=NPNE(nn,nb,ne)

C                     ..Write out the global node number
                      TEMP(1)=np
                      CALL BINWRITEFILE(FILEID,INTTYPE,1,TEMP,
     '                  REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                  ERROR,*9999)

                      DO no_njtype=1,NUM_VAR_TYPES
                        njtype=NJTYPES(no_njtype)
                        DO njj=1,NJ_LOC(njtype,0,nr)
                          nj=NJ_LOC(njtype,njj,nr)

                          nv=NVJE(nn,nb,nj,ne)

C                         ..Write out the global version number
                          TEMP(1)=nv
                          CALL BINWRITEFILE(FILEID,INTTYPE,1,TEMP,
     '                      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,
     '                      SINTDATA,ERROR,*9999)
                        ENDDO
                      ENDDO
                    ENDDO !nn
                  ENDDO !nonb
                ENDDO !ne
              ENDDO ! nonr
            ENDIF
          ELSE
            ERROR='>>Invalid subcommand'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF
      ELSE IF(COMMAND(1:5).EQ.'CLOSE') THEN
        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          ERROR='>>Ascii not implemented in ioelem yet'
          GOTO 9999

C LKC 27-SEP-1999 Unused at the moment
C          CALL CLOSEF(FILEID,ERROR,*9999)

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
          IF(ISBINFILEOPEN(FILEID))
     '      CALL BINCLOSEFILE(FILEID,ERROR,*9999)
        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF
      ELSE
        ERROR='>>Invalid command'
        GOTO 9999
      ENDIF

      CALL EXITS('IOELEM')
      RETURN
 9999 CALL ERRORS('IOELEM',ERROR)
      CALL EXITS('IOELEM')
      IF(ISBINFILEOPEN(FILEID)) CALL BINARYCLOSEFILE(FILEID,ERR,CERROR)
      RETURN 1
      END


