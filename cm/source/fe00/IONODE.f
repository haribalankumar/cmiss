      SUBROUTINE IONODE(FILEID,NKJ,NPNODE,NRLIST,NVJP,XP,
     '  COMMAND,FILEFORMAT,FILENAME,SUBCOMMAND,ERROR,*)

C#### Subroutine: IONODE
C###  Description:
C###    <HTML>
C###    Input and output of Nodal information. Will read/write an
C###    ipnode or or irnode file if ASCII format or a binnode file in
C###    BINARY format.
C###
C###    Usage:
C###    <UL>
C###      COMMAND - READ,WRITE,CLOSE
C###      FILEFORMAT - ASCII,BINARY
C###      SUBCOMMAND - OPEN,RESET,NODE_DATA,REGIONDATA,TAGDATA
C###    </UL>
C###    </HTML>


C###  Comment: Binary Node File Format
C###  Description:
C###    <HTML>
C###
C###    <H4>LEVEL 0 HEAD TAG: Node Header Information</H4>
C###    <UL>
C###      <LI> INTEGER - num regions, also indicates the number of subtags
C###    </UL>
C###
C###
C###    <H4>LEVEL 1 BRANCH TAG: </H4>
C###    <UL>
C###      <LI> INTEGER indicating 2 subtags, region and node data
C###    </UL>
C###
C###
C###    <H4>LEVEL 2a LEAF TAG: Region Data</H4>
C###    <UL>
C###      <LI> INTEGER region number
C###      <LI> INTEGER number of global coordinates
C###      <LI> INTEGER number of variables for each NJ_LOC (3)
C###      <LI> INTEGER number of nodes
C###      <LI> INTEGER number of derivatives
C###    </UL>
C###
C###
C###   <H4>LEVEL 2b LEAF TAG: Node Data</H4>
C###   <UL>
C###     <LI> for each node - INTEGER node number
C###     <UL>
C###       <LI> for each nj - INTEGER number of versions
C###       <UL>
C###         <LI> for each nv
C###         <UL>
C###           <LI> for each nk - REAL*8 information
C###         </UL>
C###       </UL>
C###     </UL>
C###   </UL>
C###
C###   </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'binf00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'

!     Parameter List
      INTEGER FILEID,NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVJP(NJM,NPM)
      CHARACTER COMMAND*(*),ERROR*(*),FILEFORMAT*(*),FILENAME*(*),
     '  SUBCOMMAND*(*)
      REAL*8 XP(NKM,NVM,NJM,NPM)

!     Local variables
      INTEGER CERROR(50),CERRLEN,ERR,FILETYPE,nj,nk,np,nr,nv,nonode,
     '  VERSION(3)
      LOGICAL ISBINFILEOPEN

! New Local Variables
      INTEGER SUBTAGS(1),njj,nonj_list,nonrlist
      REAL*8 BOB(1)


      CALL ENTERS('IONODE',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Command : '',A,'
     '    //'/'' Subcommand : '',A,'
     '    //'/'' Fileformat : '',A,'
     '    //'/'' Filename : '',A,'
     '    //'/'' FileID : '',I5)') COMMAND,SUBCOMMAND,FILEFORMAT,
     '    FILENAME,FILEID
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF


      IF(COMMAND(1:4).EQ.'READ') THEN

        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          CALL ASSERT(FILEFORMAT(1:5).EQ.'ASCII',
     '      '>>Not Implemented for ASCII',ERROR,*9999)

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

          IF(SUBCOMMAND(1:4).EQ.'OPEN') THEN
C***        Open the signal file
            FILETYPE=4        !Binary node file
            VERSION(1)=1  !expecting this version Sep-98
            VERSION(2)=0
            VERSION(3)=0
            NUMBERTAGS=1      !Just the tag header for now
            CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '        'READ','nod',FILENAME,ERROR,*9999)

C***        Write the tag header
C***        Write number of signal description bytes
C***        Write the signal description


          ELSE IF(SUBCOMMAND(1:5).EQ.'RESET') THEN

C           Reset file to beginning of node data
            CALL BINSETFILE(FILEID,0,ERROR,*9999)
            CALL BINSKIPCMHEADER(FILEID,3,ERROR,*9999)

C***        Read node header tag
            CALL READ_BIN_TAG_HEADER(FILEID,ERROR,*9999)
            CALL ASSERT(TAGINDEX.EQ.1,
     '        '>>Node information header not found',ERROR,*9999)
            CALL BINSKIPFILE(FILEID,NUMTAGBYTES,ERROR,*9999)

C***        Read level 1  tag


          ELSE IF(SUBCOMMAND(1:9).EQ.'NODE_DATA') THEN


C Check that nj_locs do not exist

            CALL ASSERT(NJ_LOC(NJL_FIBR,0,0).EQ.0,
     '        '>> Existing Fibres Info',ERROR,*9999)
            CALL ASSERT(NJ_LOC(NJL_FIEL,0,0).EQ.0,
     '        '>> Existing Field Info',ERROR,*9999)

C*** read - LEVEL 0 TAG
            CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '        NRLIST(0),REAL4DATA,REAL8DATA,CHARDATA,
     '        LOGDATA,SINTDATA,ERROR,*9999)


            DO nonrlist=1,NRLIST(0)
C*** read - LEVEL 1 TAG
              CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '          SUBTAGS(1),REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,
     '          SINTDATA(1),ERROR,*9999)
              CALL ASSERT(SUBTAGS(1).EQ.2,
     '          '>> SUBTAG <> 2',ERROR,*9999)

C*** read - LEVEL 2a TAG
              CALL BINREADFILE(IOFILE1,INTTYPE,1, !nr
     '          NRLIST(nonrlist),REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              nr=NRLIST(nonrlist)

              CALL BINREADFILE(IOFILE1,INTTYPE,1, ! coordinate
     '          INTDATA,REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

              CALL ASSERT(ITYP10(nr).EQ.INTDATA(1),
     '          '>> Wrong coordinate system in binfile',ERROR,*9999)

              IF(ITYP10(nr).EQ.4) THEN !focus if prolate
                CALL BINREADFILE(IOFILE1,DPTYPE,1, !#focus
     '            INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                FOCUS=REAL8DATA(1)
              ENDIF

              CALL BINREADFILE(IOFILE1,INTTYPE,1, !#nodes
     '          NPNODE(0,nr),REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C read NJ_LOC stuff
              DO njj=1,3 !geometry  #, fibre, field
                CALL BINREADFILE(IOFILE1,INTTYPE,1, !tot nj
     '            NJ_LOC(njj,0,nr),REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                DO nonj_list=1,NJ_LOC(njj,0,nr)
                  CALL BINREADFILE(IOFILE1,INTTYPE,1, !the nj
     '              NJ_LOC(njj,nonj_list,nr),REAL4DATA,
     '              REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                  CALL ASSERT(NJ_LOC(njj,nonj_list,nr).GT.0,
     '              '>> Invalid NJ_LOC',ERROR,*9999)
C Setup njtype
                NJ_TYPE(NJ_LOC(njj,nonj_list,nr),1)=njj
                NJ_TYPE(NJ_LOC(njj,nonj_list,nr),2)=nonj_list

                ENDDO !nonj_list
              ENDDO !njj

C*** read - LEVEL 2b TAG
              DO nonode=1,NPNODE(0,nr)
                CALL BINREADFILE(IOFILE1,INTTYPE,1, !np
     '            NPNODE(nonode,nr),REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                np=NPNODE(nonode,nr)

C Need to clear out NJ_LOC ?
                DO njj=1,3 !geometry  #, fibre, field

C                  CALL BINREADFILE(IOFILE1,INTTYPE,1, !#nj
C     '              NJ_LOC(njj,0,nr),REAL4DATA,
C     '              REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                  DO nonj_list=1,NJ_LOC(njj,0,nr)

C                  CALL BINREADFILE(IOFILE1,INTTYPE,1, !the nj
C     '              NJ_LOC(njj,nonj_list,nr),REAL4DATA,
C     '              REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                    nj=NJ_LOC(njj,nonj_list,nr)

                    CALL BINREADFILE(IOFILE1,INTTYPE,1, !#nk
     '                NKJ(nj,np),REAL4DATA,REAL8DATA,
     '                CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                    DO nk=1,NKJ(nj,np)
                      CALL BINREADFILE(IOFILE1,INTTYPE,1, !#nv
     '                  NVJP(nj,np),REAL4DATA,REAL8DATA,
     '                  CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                      DO nv=1,NVJP(nj,np)
                        CALL BINREADFILE(IOFILE1,DPTYPE,1, !#info
     '                    INTDATA,REAL4DATA,
     '                    bob(1),CHARDATA,LOGDATA,
     '                    SINTDATA,ERROR,*9999)
                        XP(nk,nv,nj,np)=bob(1)

                      ENDDO !nv
                    ENDDO !nk
                  ENDDO !nj
                ENDDO !njj
              ENDDO !np
            ENDDO !nr

          ELSE
            ERROR='>>Invalid subcommand'
            GOTO 9999
          ENDIF !SUBCOMMAND
        ELSE
          ERROR='>>Invalid file format, READ'
          GOTO 9999
        ENDIF !FILEFORMAT

C*** write

      ELSE IF(COMMAND(1:5).EQ.'WRITE') THEN
        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          CALL ASSERT(FILEFORMAT(1:5).EQ.'ASCII',
     '      '>>Not Implemented for ASCII',ERROR,*9999)


        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

          IF(SUBCOMMAND(1:4).EQ.'OPEN') THEN

C***        Open the signal file for writing
            FILETYPE=4       !Binary node file
            VERSION(1)=1 !Current version Sep-98
            VERSION(2)=0
            VERSION(3)=0
            NUMBERTAGS=1     !Just header tag for the moment.

            CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '        'WRITE','nod',FILENAME,ERROR,*9999)

          ELSE IF(SUBCOMMAND(1:9).EQ.'NODE_DATA') THEN

C*** write - LEVEL 0 TAG
            CALL BINWRITEFILE(IOFILE1,INTTYPE,1,
     '        NRLIST(0),REAL4DATA,REAL8DATA,CHARDATA,
     '        LOGDATA,SINTDATA,ERROR,*9999)

            SUBTAGS(1)=2
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)

C*** write - LEVEL 1 TAG
              CALL BINWRITEFILE(IOFILE1,INTTYPE,1,
     '          SUBTAGS(1),REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,
     '          SINTDATA,ERROR,*9999)

C*** write LEVEL 2a TAG
              CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !nr
     '          NRLIST(nr),REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

              CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !coord sys
     '          ITYP10(nr),REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

              IF(ITYP10(nr).EQ.4) THEN !focus if prolate
                REAL8DATA(1)=FOCUS
                CALL BINWRITEFILE(IOFILE1,DPTYPE,1, !#focus
     '            INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

              ENDIF

              CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !#nodes
     '          NPNODE(0,NRLIST(nr)),REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C write NJ_LOC stuff
              DO njj=1,3 !geometry, fibre, field
                CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !tot nj
     '            NJ_LOC(njj,0,nr),REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                DO nonj_list=1,NJ_LOC(njj,0,nr)
                  CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !tot nj
     '              NJ_LOC(njj,nonj_list,nr),REAL4DATA,
     '              REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                ENDDO !nonj_list
              ENDDO !njj


C*** write LEVEL 2b TAG
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !np
     '            NPNODE(nonode,nr),REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                DO njj=1,3 !geometry, fibre, field

C                      CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !#nj
C     '                  NJ_LOC(njj,0,nr),REAL4DATA,
C     '                  REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                  DO nonj_list=1,NJ_LOC(njj,0,nr)
                    nj=NJ_LOC(njj,nonj_list,nr)

C                        CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !the nj
C     '                    nj,REAL4DATA,REAL8DATA,
C     '                    CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                    CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !#nk
     '                NKJ(nj,np),REAL4DATA,REAL8DATA,
     '                CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                    DO nk=1,NKJ(nj,np)
                      CALL BINWRITEFILE(IOFILE1,INTTYPE,1, !#nv
     '                  NVJP(nj,np),REAL4DATA,REAL8DATA,
     '                  CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

                      DO nv=1,NVJP(nj,np)
                        bob(1)=XP(nk,nv,nj,np)
                        CALL BINWRITEFILE(IOFILE1,DPTYPE,1, !#info
     '                    INTDATA,REAL4DATA,
     '                    bob,CHARDATA,LOGDATA,
     '                    SINTDATA,ERROR,*9999)

C note : remove bob's

                      ENDDO !nv
                    ENDDO !nk
                  ENDDO !nj
                ENDDO !njj
              ENDDO !np
            ENDDO ! nr



C need to take into account reading 1 out of a group of regions?


          ELSE IF(SUBCOMMAND(1:11).EQ.'REGION_DATA') THEN

          ELSE
            ERROR='>>Invalid subcommand'
            GOTO 9999
          ENDIF !SUBCOMMAND
        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF !FILEFORMAT


      ELSE IF(COMMAND(1:5).EQ.'CLOSE') THEN

C*** CLOSE-ASCII
        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          CALL ASSERT(FILEFORMAT(1:5).EQ.'ASCII',
     '      '>>Not Implemented',ERROR,*9999)

C*** CLOSE-BINARY
        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
          CALL BINARYCLOSEFILE(FILEID,ERR,CERROR)
C          CALL BINCLOSEFILE(FILEID,CERROR,*9999)
          IF(ERR.NE.0) THEN
            CALL CSTRINGLEN(CERRLEN,CERROR)
            CALL C2FSTRING(CERROR,CERRLEN,ERROR)
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>Invalid file format'
          GOTO 9999
        ENDIF !FILEFORMAT
      ELSE
        ERROR='>>Invalid command'
        GOTO 9999
      ENDIF !COMMAND

      CALL EXITS('IONODE')
      RETURN
 9999 CALL ERRORS('IONODE',ERROR)
      CALL EXITS('IONODE')
      IF(ISBINFILEOPEN(FILEID)) CALL BINARYCLOSEFILE(FILEID,ERR,CERROR)
      RETURN 1
      END



