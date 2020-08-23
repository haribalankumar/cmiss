      SUBROUTINE OPEN_BIN_FILE(FILEID,TYPE,NUMTAGS,VERSION,COMMAND,
     '  EXTENSION,FILE,ERROR,*)

C#### Subroutine: OPEN_BIN_FILE
C###  Description:
C###    OPEN_BIN_FILE opens binary file called file.bin*** (where
C###    *** is the extension) specified by FILEID and reads/writes
C###    (specified by command) file type, version number, heading
C###    and number of tags. If the file is being opened for reading
C###    fileversion is set to the version in the file. If the file
C###    is opened for writing FILEVERSION is set to VERSION.

C#### Comment: BINARY FILE FORMAT
C###  Description:
C###    <HTML>
C###    Each CMISS binary file is organised as follows:
C###    <UL>
C###      <LI> Header section 
C###        <UL>
C###          <LI> An identity header section
C###          <LI> A machine header section
C###          <LI> A file header section
C###        </UL>
C###      <LI> Data section made up of tagged data 
C###    </UL>
C###
C###    <H2> Header Section </H2>
C###
C###    There are three sections in the header:
C###    <UL>
C###      <LI> an identity header section,
C###      <LI> a machine header section,
C###      <LI> a file header section.
C###    </UL>
C###
C###    <H3> Identity Header </H3>
C###
C###    The format of the identity header section is:
C###    <UL>
C###      <LI>
C###        1 identity byte, currently set to the value 0x07, intended
C###        to indicate that the file is a CMISS binary file.
C###      <LI>
C###        1 byte to indicate the revision of the binary file format
C###        (currently 0x02).
C###    </UL>
C###
C###    <H3> Machine Header </H3>
C###
C###    The format of the machine header section for the current
C###    binary file format revision (0x02) is:
C###    <UL>
C###      <LI>
C###        1 byte to indicate the number (in standard base two) of
C###        following bytes in the machine header section, which is
C###        21 (= 0x15).
C###      <LI>
C###        1 byte to indicate the machine type which created the
C###        file (0x01 - DEC Alpha, 0x02 - SGI, 0x03 - IBM,
C###        0x04 - Cray). Not used - always set to MACH_ANY now to aid
C###        example testing (DPN 20 January 2003)
C###      <LI>
C###        1 byte to indicate the operating system type of the machine
C###        that created the file (redundant - for future development).
C###      <LI>
C###        1 byte to indicate the endian ordering of the numbers in
C###        in the file (0x01 - big endian, 0x02 - little endian).
C###        Currently all numbers are converted to big endian format.
C###      <LI>
C###        1 byte to indicate the format of characters (0x01 -
C###        Ascii, 0x02 - Unicode).
C###      <LI>
C###        1 byte to indicate the format of integers (0x01 -
C###        twos complement, 0x02 - Signed magnitude).
C###      <LI>
C###        1 byte to indicate the format of unsigned integers (0x01 -
C###        base two).
C###      <LI>
C###        1 byte to indicate the format of standard single precision
C###        numbers (0x01 - IEEE single standard).
C###      <LI>
C###        1 byte to indicate the format of standard double precision
C###        numbers (0x01 - IEEE double standard).
C###      <LI>
C###        1 byte to indicate the number of bytes a character type
C###        takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes an integer type
C###        takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes a short integer type
C###        takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes a long integer type
C###        takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes an unsigned integer
C###        type takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes an unsigned  short
C###        integer type takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes an unsigned long
C###        integer type takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes a single precision
C###        type takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes a double precision
C###        type takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes a logical
C###        type takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes a single precision
C###        complex type takes up in the file.
C###      <LI>
C###        1 byte to indicate the number of bytes a double precision
C###        complex type takes up in the file.
C###    </UL>
C###
C###    The format for revision 0x01 is
C###    <UL>
C###    <LI>1 byte to indicate the number (in standard base two) of
C###        following bytes in the machine header section. This is
C###        11 for this revision.
C###    <LI>1 byte to indicate the machine type which created the
C###        file (0x01 - DEC Alpha, 0x02 - SGI, 0x03 - IBM,
C###        0x04 - Cray).
C###    <LI>1 byte to indicate the operating system type of the machine
C###        which created the file (redundant).
C###    <LI>1 byte to indicate the endian ordering of the numbers in
C###        in the file (0x01 - big endian, 0x02 - little endian).
C###        Currently all numbers are converted to big endian format.
C###    <LI>1 byte to indicate the format of standard single precision
C###        numbers (0x01 - IEEE 754 single standard).
C###    <LI>1 byte to indicate the format of standard double precision
C###        numbers (0x01 - IEEE 754 double standard).
C###    <LI>1 byte to indicate the number of bytes a character type
C###        takes up in the file.
C###    <LI>1 byte to indicate the number of bytes a integer type
C###        takes up in the file.
C###    <LI>1 byte to indicate the number of bytes a short integer type
C###        takes up in the file.
C###    <LI>1 byte to indicate the number of bytes a single precision
C###        type takes up in the file.
C###    <LI>1 byte to indicate the number of bytes a double precision
C###        type takes up in the file.
C###    <LI>1 byte to indicate the number of bytes a logical
C###        type takes up in the file.
C###    </UL>
C###
C###    The format for revision 0x00 is:
C###    <UL>
C###      <LI> 2 dummy bytes (to be skipped).
C###    </UL>
C###
C###    <H3> File Header </H3>
C###
C###    The format of the file header section for the current
C###    binary file format revision (0x02) is:
C###    <UL>
C###    <LI>Integer to specify the file type. Current binary
C###        file types are: 1 - Binary matrix file, 2 - Binary time
C###        series file, 3 - Binary signal file, 4 - Binary node file??,
C###        5 - Binary element file??
C###    <LI>Three integers to specify the version of the file in the
C###        form xxx.xxx.xxx
C###    <LI>Integer to specify the number of bytes in the heading.
C###    <LI>The heading (as a string of character bytes).
C###    <LI>Integer to specify how many tags are in the top level of
C###        the binary file.
C###    </UL>
C###
C###    If the binary file format number is below 0x02 then the format
C###    of the file header section is:
C###    <UL>
C###    <LI>Integer to specify the file type. Current binary
C###        file types are: 1 - Binary matrix file, 2 - Binary time
C###        series file, 3 - Binary signal file.
C###    <LI>Single precision number to specify the version of the file.
C###    <LI>Integer to specify the number of bytes in the heading.
C###    <LI>The heading (as a string of character bytes).
C###    <LI>Integer to specify how many `tags' of data are in the file.
C###    </UL>
C###
C###    <H2> Data Section </H2>
C###
C###    The rest of the file is made up of a hierarchy of `tags' which
C###    contain the actual data.  These tags are connected in a tree like
C###    structure giving two types of tags:  if a tag has further tags
C###    below it in the hierarchy it is termed a &quot;node tag&quot;;
C###    otherwise it is termed a &quot;leaf tag&quot;.  For the current
C###    binary file format revision (0x02), the format of each tag is:
C###    <UL>
C###    <LI>Integer to specify the type of tag.
C###    <LI>Integer to specify the number of bytes in the tag heading.
C###    <LI>The tag heading (as a string of character bytes).
C###    <LI>An integer to specify the number of tags below this tag.
C###        If this number is > 0 the tag is a node tag.
C###        If this number is = 0 the tag is a leaf tag.
C###    <LI>If the tag is a leaf tag
C###      <UL>
C###      <LI>
C###        Integer to specify the number of bytes in the tag
C###        (excluding the tag header information). NOTE: this restricts
C###        tag sizes to be under 2 GB. Future revisions of the binary
C###        format will change this quantity to a long integer but this
C###        will have to wait until Intel based machines can handle
C###        64-bit integers.
C###      <LI>The tag data.
C###      </UL>
C###    </UL>
C###    
C###    If the binary file format number is below 0x02, then there is
C###    only one level of tags, so all tags are leaf tags.
C###    The format of each tag is:
C###    <UL>
C###    <LI>Integer to specify the type of tag.
C###    <LI>Integer to specify the number of bytes in the tag
C###        heading.
C###    <LI>The tag heading (as a string of character bytes).
C###    <LI>Integer to specify the number of bytes in the tag
C###        (excluding the tag header information).
C###    <LI>The tag data.
C###    </UL>
C###    Note: that if any of the byte numbers are unknown they will have
C###    the value of 0xFF.
C###    </HTML>
C###  See-Also: BINARY MATRIX FILE FORMAT, BINARY HISTORY FILE FORMAT,
C###    BINARY SIGNAL FILE FORMAT, BINARY NODE FILE FORMAT,
C###    BIANRY ELEMENT FILE FORMAT

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'binf00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER FILEID,TYPE,NUMTAGS,VERSION(3)
      CHARACTER COMMAND*(*),EXTENSION*(*),FILE*(*),ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),ERR,HEADSZE(1),IBEG,IBEG1,IEND,IEND1,
     '  FILETYPE(1),NMTAGS(1),NUMMACHHEADERBYTES,
     '  TYPE2(1)
      REAL FVERSION(1)
      CHARACTER FILENAME*(MXCH)
      LOGICAL ISBINFILEOPEN

      CALL ENTERS('OPEN_BIN_FILE',*9999)

      CALL STRING_TRIM(FILE,IBEG,IEND)
      CALL STRING_TRIM(EXTENSION,IBEG1,IEND1)
      FILENAME=FILE(IBEG:IEND)//'.bin'//EXTENSION(IBEG1:IEND1)


C LKC 18-APR-2000 Not sure why we want to do this.
C  DOCDIR is setup in CHECKF but this is usally not called for these
C  files. The examples path eg. $example should be used before each file.
C
CC rgb 25/10/1999 append file name to example pathname
C      IF(DOCDIR) THEN
C        CALL STRING_TRIM(EXAMPLES_DIR,IBEG2,IEND2)
C        FILENAME=EXAMPLES_DIR(IBEG2:IEND2)//FILENAME
C      ENDIF




      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' FileID : '',I5,/'' Command : '',A,'
     '    //'/'' Filename : '',A)') FILEID,COMMAND,FILENAME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      CALL BINOPENFILE(FILEID,COMMAND,FILENAME,ERROR,*9999)

      IF(COMMAND(1:4).EQ.'READ') THEN

C*** Read identity header section
        CALL BINREADFILE(FILEID,CHARTYPE,2,INTDATA,REAL4DATA,
     '    REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        CALL ASSERT(CHARDATA(1:1).EQ.CHAR(7),
     '    '>>Not a CMISS binary file',ERROR,*9999)
        FILEBINVERTYPE(FILEID)=CHARDATA(2:2)
C*** Read machine header section
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' File identity header information:'','
     '      //'/'' Binary file header format version : '',I1)')
     '      ICHAR(FILEBINVERTYPE(FILEID))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(FILEBINVERTYPE(FILEID).EQ.CHAR(0)) THEN
C         Identity format 0 - has 2 extra dummy bytes. Skip them
          WRITE(OP_STRING,'('' >>Warning: Old binary file identity '
     '      //'format found. Please update file'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL BINREADFILE(FILEID,CHARTYPE,2,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        ELSE IF(FILEBINVERTYPE(FILEID).EQ.CHAR(1)) THEN
C         Identity format 1 - has machine header section.
          WRITE(OP_STRING,'('' >>Warning: Old binary file identity '
     '      //'format found. Please update file'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL BINREADFILE(FILEID,CHARTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NUMMACHHEADERBYTES=ICHAR(CHARDATA(1:1))
          CALL BINREADFILE(FILEID,CHARTYPE,NUMMACHHEADERBYTES,
     '      INTDATA,REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '      ERROR,*9999)
          IF(NUMMACHHEADERBYTES.EQ.11) THEN
            FILEMACHTYPE(FILEID)=CHARDATA(1:1)
            FILEOSTYPE(FILEID)=CHARDATA(2:2)
            FILEENDIANTYPE(FILEID)=CHARDATA(3:3)
            FILESPFORMTYPE(FILEID)=CHARDATA(4:4)
            FILEDPFORMTYPE(FILEID)=CHARDATA(5:5)
            FILECHARSIZE(FILEID)=CHARDATA(6:6)
            FILEINTSIZE(FILEID)=CHARDATA(7:7)
            FILESINTSIZE(FILEID)=CHARDATA(8:8)
            FILESPSIZE(FILEID)=CHARDATA(9:9)
            FILEDPSIZE(FILEID)=CHARDATA(10:10)
            FILELOGSIZE(FILEID)=CHARDATA(11:11)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' File machine header information:'','
     '          //'/,'' Machine type : '',I3,'
     '          //'/,'' Operating system type : '',I3,'
     '          //'/,'' Endian type : '',I3,'
     '          //'/,'' Single precision format type : '',I3,'
     '          //'/,'' Double precision format type : '',I3,'
     '          //'/,'' Character size : '',I3,'
     '          //'/,'' Integer size : '',I3,'
     '          //'/,'' Short integer size : '',I3,'
     '          //'/,'' Single precision size : '',I3,'
     '          //'/,'' Double precision size : '',I3,'
     '          //'/,'' Logical size : '',I3)')
     '          ICHAR(FILEMACHTYPE(FILEID)),ICHAR(FILEOSTYPE(FILEID)),
     '          ICHAR(FILEENDIANTYPE(FILEID)),
     '          ICHAR(FILESPFORMTYPE(FILEID)),
     '          ICHAR(FILEDPFORMTYPE(FILEID)),
     '          ICHAR(FILECHARSIZE(FILEID)),ICHAR(FILEINTSIZE(FILEID)),
     '          ICHAR(FILESINTSIZE(FILEID)),ICHAR(FILESPSIZE(FILEID)),
     '          ICHAR(FILEDPSIZE(FILEID)),ICHAR(FILELOGSIZE(FILEID))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ELSE
            ERROR='>>Invalid number of machine header bytes'
            GOTO 9999
          ENDIF
        ELSE IF(FILEBINVERTYPE(FILEID).EQ.CHAR(2)) THEN
C         Identity format 2 - has expanded machine header section and
C         subtags
          CALL BINREADFILE(FILEID,CHARTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NUMMACHHEADERBYTES=ICHAR(CHARDATA(1:1))
          CALL BINREADFILE(FILEID,CHARTYPE,NUMMACHHEADERBYTES,
     '      INTDATA,REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '      ERROR,*9999)
          IF(NUMMACHHEADERBYTES.EQ.16) THEN
            FILEMACHTYPE(FILEID)=CHARDATA(1:1)
            FILEOSTYPE(FILEID)=CHARDATA(2:2)
            FILEENDIANTYPE(FILEID)=CHARDATA(3:3)
            FILECHARFORMTYPE(FILEID)=CHARDATA(4:4)
            FILEINTFORMTYPE(FILEID)=CHARDATA(5:5)
            FILESPFORMTYPE(FILEID)=CHARDATA(6:6)
            FILEDPFORMTYPE(FILEID)=CHARDATA(7:7)
            FILECHARSIZE(FILEID)=CHARDATA(8:8)
            FILEINTSIZE(FILEID)=CHARDATA(9:9)
            FILESINTSIZE(FILEID)=CHARDATA(10:10)
            FILELINTSIZE(FILEID)=CHARDATA(11:11)
            FILESPSIZE(FILEID)=CHARDATA(12:12)
            FILEDPSIZE(FILEID)=CHARDATA(13:13)
            FILELOGSIZE(FILEID)=CHARDATA(14:14)
            FILESPCSIZE(FILEID)=CHARDATA(15:15)
            FILEDPCSIZE(FILEID)=CHARDATA(16:16)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' File machine header information:'','
     '          //'/,'' Machine type : '',I3,'
     '          //'/,'' Operating system type : '',I3,'
     '          //'/,'' Endian type : '',I3,'
     '          //'/,'' Character format type : '',I3,'
     '          //'/,'' Integer format type : '',I3,'
     '          //'/,'' Single precision format type : '',I3,'
     '          //'/,'' Double precision format type : '',I3,'
     '          //'/,'' Character size : '',I3,'
     '          //'/,'' Integer size : '',I3,'
     '          //'/,'' Short integer size : '',I3,'
     '          //'/,'' Long integer size : '',I3,'
     '          //'/,'' Single precision size : '',I3,'
     '          //'/,'' Double precision size : '',I3,'
     '          //'/,'' Logical size : '',I3,'
     '          //'/,'' Single precision complex size : '',I3,'
     '          //'/,'' Double precision complex size : '',I3)')
     '          ICHAR(FILEMACHTYPE(FILEID)),ICHAR(FILEOSTYPE(FILEID)),
     '          ICHAR(FILEENDIANTYPE(FILEID)),FILECHARFORMTYPE(FILEID),
     '          FILEINTFORMTYPE(FILEID),
     '          ICHAR(FILESPFORMTYPE(FILEID)),
     '          ICHAR(FILEDPFORMTYPE(FILEID)),
     '          ICHAR(FILECHARSIZE(FILEID)),
     '          ICHAR(FILEINTSIZE(FILEID)),ICHAR(FILESINTSIZE(FILEID)),
     '          ICHAR(FILELINTSIZE(FILEID)),ICHAR(FILESPSIZE(FILEID)),
     '          ICHAR(FILEDPSIZE(FILEID)),ICHAR(FILELOGSIZE(FILEID)),
     '          ICHAR(FILESPCSIZE(FILEID)),
     '          ICHAR(FILEDPCSIZE(FILEID))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ELSE
            ERROR='>>Invalid number of machine header bytes'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>Unknown binary file header identity format'
          GOTO 9999
        ENDIF

C*** Read file header section
C*** Read file type, version and header

        CALL BINREADFILE(FILEID,INTTYPE,1,FILETYPE,REAL4DATA,
     '    REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        IF(FILETYPE(1).NE.TYPE) THEN
          WRITE(OP_STRING,'('' >>Warning: File has a different '
     '      //'file type'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''     File type is          '',I3)')
     '      FILETYPE(1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''     Expected file type is '',I3)') TYPE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(FILEBINVERTYPE(FILEID).EQ.CHAR(0).OR.
     '    FILEBINVERTYPE(FILEID).EQ.CHAR(1)) THEN
          CALL BINREADFILE(FILEID,SPTYPE,1,INTDATA,FVERSION,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          BINFILEVERSION(1,FILEID)=INT(FVERSION(1))
          FVERSION(1)=(FVERSION(1)-REAL(BINFILEVERSION(1,FILEID)))*10.0
          BINFILEVERSION(2,FILEID)=INT(FVERSION(1))
          FVERSION(1)=(FVERSION(1)-REAL(BINFILEVERSION(2,FILEID)))*10.0
          BINFILEVERSION(3,FILEID)=INT(FVERSION(1))
          IF(BINFILEVERSION(1,FILEID).NE.VERSION(1)) THEN
            WRITE(OP_STRING,'('' >>Warning: File has a different '
     '        //'version number'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C LKC 27-APR-2000 Redoing this section a bit better
C
C            WRITE(OP_STRING,'(''     File version is          '',I2,'
C     '        //'''.'',I2,''.'',I2)') BINFILEVERSION(1,FILEID),
C     '        BINFILEVERSION(2,FILEID),BINFILEVERSION(3,FILEID)
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            WRITE(OP_STRING,'(''     Expected file version is '',I2,'
C     '        //'''.'',I2,''.'',I2)') VERSION(1),VERSION(2),VERSION(3)



            IEND=0
            CALL APPENDC(IEND,'     File version is          '
     '        ,OP_STRING(1))
            CALL APPENDI(IEND,BINFILEVERSION(1,FILEID),OP_STRING(1))
            CALL APPENDC(IEND,'.',OP_STRING(1))
            CALL APPENDI(IEND,BINFILEVERSION(2,FILEID),OP_STRING(1))
            CALL APPENDC(IEND,'.',OP_STRING(1))
            CALL APPENDI(IEND,BINFILEVERSION(3,FILEID),OP_STRING(1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


            IEND=0
            CALL APPENDC(IEND,' Expected version is          '
     '        ,OP_STRING(1))
            CALL APPENDI(IEND,VERSION(1),OP_STRING(1))
            CALL APPENDC(IEND,'.',OP_STRING(1))
            CALL APPENDI(IEND,VERSION(2),OP_STRING(1))
            CALL APPENDC(IEND,'.',OP_STRING(1))
            CALL APPENDI(IEND,VERSION(3),OP_STRING(1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)



            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE
          CALL BINREADFILE(FILEID,INTTYPE,3,BINFILEVERSION(1,FILEID),
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          IF(BINFILEVERSION(1,FILEID).NE.VERSION(1).OR.
     '       BINFILEVERSION(2,FILEID).NE.VERSION(2).OR.
     '       BINFILEVERSION(3,FILEID).NE.VERSION(3)) THEN
            WRITE(OP_STRING,'('' >>Warning: File has a different '
     '        //'version number'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''     File version is          '',I2,'
     '        //'''.'',I2,''.'',I2)') BINFILEVERSION(1,FILEID),
     '        BINFILEVERSION(2,FILEID),BINFILEVERSION(3,FILEID)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''     Expected file version is '',I2,'
     '        //'''.'',I2,''.'',I2)') VERSION(1),VERSION(2),VERSION(3)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

C***    Read Heading
        CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '    REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        HEADSIZE=INTDATA(1)
        CALL BINREADFILE(FILEID,CHARTYPE,HEADSIZE,INTDATA,
     '    REAL4DATA,REAL8DATA,HEADING,LOGDATA,SINTDATA,ERROR,*9999)
        IF(HEADSIZE.EQ.0) THEN
          WRITE(OP_STRING,'('' File heading: '')')
        ELSE
          WRITE(OP_STRING,'('' File heading: '',A)')
     '      HEADING(1:HEADSIZE)
        ENDIF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C***    Read num_tags
        CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,
     '    REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        NUMTAGS=INTDATA(1)

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' File header information:'','
     '      //'/'' File Type              : '',I5,'
     '      //'/'' File Version           : '',I2,''.'',I2,''.'',I2,'
     '      //'/'' Number of header bytes : '',I5)')
     '      FILETYPE(1),BINFILEVERSION(1,FILEID),
     '      BINFILEVERSION(2,FILEID),BINFILEVERSION(3,FILEID),
     '      HEADSIZE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

      ELSE IF(COMMAND(1:5).EQ.'WRITE') THEN

CPB 10/3/99 Now binary file revision 0x02

C*** Write identity header section
        CHARDATA(1:1)=CHAR(7)
        CHARDATA(2:2)=CHAR(2)
        CALL BINWRITEFILE(FILEID,CHARTYPE,2,INTDATA,REAL4DATA,REAL8DATA,
     '    CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        FILEBINVERTYPE(FILEID)=CHAR(2)

C*** Write the machine header section
        CHARDATA(1:1)=CHAR(16)
c DPN 20 January 2003 - Changing the machine type bit to be consistent
c       across all platforms since it doesn't seem to be used and causes
c       problems with multi-platform testing of binary files.
        !CHARDATA(2:2)=MACHTYPE
        CHARDATA(2:2)=CHAR(MACH_ANY)
        CHARDATA(3:3)=OSTYPE
c cpb 10/22/95 All numbers are in big endian for the moment
C        CHARDATA(4:4)=ENDIANTYPE
        CHARDATA(4:4)=CHAR(MACH_BIGENDIAN)
        CHARDATA(5:5)=CHARFORMTYPE
        CHARDATA(6:6)=INTFORMTYPE
        CHARDATA(7:7)=SPFORMTYPE
        CHARDATA(8:8)=DPFORMTYPE
        CHARDATA(9:9)=CHAR(CHARSIZE)
        CHARDATA(10:10)=CHAR(INTSIZE)
        CHARDATA(11:11)=CHAR(SINTSIZE)
        CHARDATA(12:12)=CHAR(LINTSIZE)
        CHARDATA(13:13)=CHAR(SPSIZE)
        CHARDATA(14:14)=CHAR(DPSIZE)
        CHARDATA(15:15)=CHAR(LOGSIZE)
        CHARDATA(16:16)=CHAR(SPCSIZE)
        CHARDATA(17:17)=CHAR(DPCSIZE)

C cpb 16/3/97 Initialise the filetypes for writing
        FILEMACHTYPE(FILEID)=MACHTYPE
        FILEOSTYPE(FILEID)=OSTYPE
        FILEENDIANTYPE(FILEID)=CHAR(MACH_BIGENDIAN)
        FILECHARFORMTYPE(FILEID)=CHARFORMTYPE
        FILEINTFORMTYPE(FILEID)=INTFORMTYPE
        FILESPFORMTYPE(FILEID)=SPFORMTYPE
        FILEDPFORMTYPE(FILEID)=DPFORMTYPE
        FILECHARSIZE(FILEID)=CHAR(CHARSIZE)
        FILEINTSIZE(FILEID)=CHAR(INTSIZE)
        FILESINTSIZE(FILEID)=CHAR(SINTSIZE)
        FILELINTSIZE(FILEID)=CHAR(SINTSIZE)
        FILESPSIZE(FILEID)=CHAR(SPSIZE)
        FILEDPSIZE(FILEID)=CHAR(DPSIZE)
        FILELOGSIZE(FILEID)=CHAR(LOGSIZE)
        FILESPCSIZE(FILEID)=CHAR(SPCSIZE)
        FILEDPCSIZE(FILEID)=CHAR(DPCSIZE)
        CALL BINWRITEFILE(FILEID,CHARTYPE,17,INTDATA,REAL4DATA,
     '    REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C*** Write the file header section

C*** Write file type, version and header

        TYPE2(1)=TYPE
        CALL BINWRITEFILE(FILEID,INTTYPE,1,TYPE2,REAL4DATA,REAL8DATA,
     '    CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

        BINFILEVERSION(1,FILEID)=VERSION(1)
        BINFILEVERSION(2,FILEID)=VERSION(2)
        BINFILEVERSION(3,FILEID)=VERSION(3)

        CALL BINWRITEFILE(FILEID,INTTYPE,3,BINFILEVERSION(1,FILEID),
     '    REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

        CALL STRING_TRIM(HEADING,IBEG,IEND)
        IF(IBEG.GT.IEND) THEN
          HEADSIZE=0
        ELSE IF(IBEG.EQ.IEND.AND.HEADING(IBEG:IEND).EQ.' ') THEN
          HEADSIZE=0
        ELSE
          HEADSIZE=IEND-IBEG+1
        ENDIF
        HEADSZE(1)=HEADSIZE
        CALL BINWRITEFILE(FILEID,INTTYPE,1,HEADSZE,REAL4DATA,
     '    REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

        CALL BINWRITEFILE(FILEID,CHARTYPE,HEADSIZE,INTDATA,REAL4DATA,
     '    REAL8DATA,HEADING,LOGDATA,SINTDATA,ERROR,*9999)
        NMTAGS(1)=NUMTAGS
        CALL BINWRITEFILE(FILEID,INTTYPE,1,NMTAGS,REAL4DATA,REAL8DATA,
     '    CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

      ENDIF ! read/write

      CALL EXITS('OPEN_BIN_FILE')
      RETURN
 9999 CALL ERRORS('OPEN_BIN_FILE',ERROR)
      CALL EXITS('OPEN_BIN_FILE')
      IF(ISBINFILEOPEN(FILEID)) CALL BINARYCLOSEFILE(FILEID,ERR,CERROR)
      RETURN 1
      END


