      SUBROUTINE GINPUT(IPTYPE,IUNIT,NDAT,A_DATA,A_DEFLT,C_DATA,C_DEFLT,
     '  CLINE,I_DATA,I_DEFLT,IMIN,IMAX,L_DATA,L_DEFLT,R_DATA,R_DEFLT,
     &  RMIN,RMAX,FILEIP,FILEIP_FIRST,INFO,ERROR,*,*)

C#### Subroutine: GINPUT
C###  Description:
C###    GINPUT reads answer,char,int,log or real*8 values from IUNIT.
C**** IUNIT has a logical record length of IRECL.
C**** Variables separated by blanks or commas are read into
C**** *DATA(nValues) nValues=1..
C**** Limiting values for an integer read are defined by IMIN and IMAX.
C**** Limiting values for a real*8 read are defined by RMIN and RMAX.
C**** If an end of file or carriage return is detected while attempting
C**** to read *DATA(nValues) it is assigned the value *DEFLT(nValues).
C**** No range check is made on the default variables.
C****
C**** If help request is made (?) a help file (based on INFO) is opened.
C**** If a quit request is made (q) the program is stopped.
C**** If an edit request is made (e) control returns to
C**** the command line and the file is closed.
C**** If restart request made (r) control returns to
C**** the command line and the file is deleted.
C**** Typing 'a' gives 'apply' return (apply previous val to subsequent)
C**** Typing 'd' gives a 'default' return
C**** If a backup request is made (b unless reading char in which case
C**** back is required) the second alternate return path is
C**** used to repeat the previous command.
C****
C**** If a file request is made (f) the program prompts for
C**** file name (FILE07) and FILEIP is set .true. & data is read from
C**** FILE07.iphist during a time-dep iteration in subroutine MARCH1, or
C**** FILE07.ippres during a time-dep loading in subroutine NONLIN.
C**** If a file request is made (g) the program prompts for
C**** file name (FILE07) (1st time only) & the data is read from this
C**** file immediately.
C****
C**** If a ? is typed, DOCUM subroutine is called with filename DOCFILE
C**** and labels consisting of DOCLABEL appended with the integer INFO.
C****
C**** If an error is detected or a help request is made ERROR is
C**** returned with a diagnostic message.
C**** NOTE: The Q edit descriptor in the read format is nonstandard!

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'docu00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER I_DATA(0:*),I_DEFLT(*),IMAX,IMIN,INFO,IPTYPE,IUNIT,NDAT
      REAL*8 R_DATA(*),R_DEFLT(*),RMAX,RMIN
      CHARACTER A_DATA(*)*(*),A_DEFLT(*)*(*),C_DATA(*)*(*),
     &  C_DEFLT(*)*(*),CLINE*(*),ERROR*(*)
      LOGICAL FILEIP,FILEIP_FIRST,L_DATA(*),L_DEFLT(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IBEG1,IBEG2,ic,ICOL,IEND,IEND1,IEND2,
     '  ILISTMBR,index,IOSTAT,
     '  ISTART,max_dim,n,NCHAR,nd,nValues,
     '  noelem,nogrel,nogrfa,nogrga,nogrgr,nogrli,nogrno,NTIL,num_points
      CHARACTER CHAR*5,CHAR10*10,CHARTEMP*20,CHDATA*16,C1*4,C4*4
      LOGICAL FIRST,FOUND_GROUP,IVALID,RVALID

      CALL ENTERS('GINPUT',*9999)
      nValues=0 !keeps track of #values read in
      FIRST=.TRUE.

!news ajp 3 11 95
C      IF(IPTYPE.EQ.3) THEN !handling different lengths of strings
C        max_dim=2999 !check with inout00.cmn if changed
C      ELSE
C        max_dim=2999 !check with inout00.cmn if changed
C      ENDIF
!newe

! rgb 27-09-1999 more complete - still need to check with inout00.cmn
      IF(IPTYPE.EQ.IPMESS) THEN
        max_dim=1
      ELSEIF(IPTYPE.EQ.IPANSW) THEN ! answer yes/no
        max_dim=4
      ELSEIF(IPTYPE.EQ.IPCHAR) THEN ! Character
        max_dim=20
      ELSEIF(IPTYPE.EQ.IPINTE) THEN ! Integer
        max_dim=9000
      ELSEIF(IPTYPE.EQ.IPLOGI) THEN ! Logical
        max_dim=4
      ELSEIF(IPTYPE.EQ.IPREAL) THEN ! Real
        max_dim=99
      ELSE
        ERROR='Incorrect input type'
        GOTO 9999
      ENDIF
! rgb end

      DO WHILE(FIRST.OR.nValues.LT.NDAT) !allows data spread over several lines
        FIRST=.FALSE.

C!!! GINOUT expects CLINE to contain the whole string but it will only
C!!! do so if all data is on one line.
        CALL READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*9999)

        IF(IOSTAT.EQ.0.AND.NCHAR.GT.0) THEN !not default
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' NCHAR='',I3)') NCHAR
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        !Deal with each data value in turn
          ICOL=1
          DO nd=nValues+1,max_dim !deal with remaining data values
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' nd='',I3,'' ICOL='',I3,'' NCHAR='','
     '          //'I3)') nd,ICOL,NCHAR
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(ICOL.GT.NCHAR) GOTO 4
            ISTART=0
            DO ic=ICOL,NCHAR !examine characters left in CLINE
              IF((CLINE(ic:ic).EQ.' ').OR.(CLINE(ic:ic).EQ.',')) THEN
                IF(ISTART.GT.0) THEN
                  IEND=ic-1
                  GOTO 2
                ENDIF
              ELSE !character is not blank or ','
                IF(ISTART.EQ.0) THEN
                  ISTART=ic
                ENDIF
              ENDIF
            ENDDO !ic
            IEND=NCHAR
    2       IF(ISTART.GT.0) THEN
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' istart='',I4,'' iend='',I4)')
     '            ISTART,IEND
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              ICOL=IEND+2
              CHDATA=CLINE(ISTART:IEND)
            !Check for ?
              IF(CHDATA(1:1).EQ.'?') THEN
                WRITE(UNIT=CHAR,FMT='(I5)') INFO
                CALL STRING_TRIM(DOCLABEL,IBEG1,IEND1)
                CALL STRING_TRIM(CHAR,IBEG2,IEND2)
                CALL DOCUM(DOCFILE,'doc',' '//DOCLABEL(IBEG1:IEND1)//' '
     '            //CHAR(IBEG2:IEND2),ERROR,*9997)
 9997           IF(ERROR.EQ.' ') THEN
                  ERROR=' ?'
                  GO TO 9998
                ELSE IF(ERROR(1:14).EQ.'File not found') THEN
                  WRITE(OP_STRING,*) '>>Help not available'
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR=' ?'
                  GO TO 9998
                ELSE
                  GO TO 9999
                ENDIF
              ENDIF !chdata=?
            !Deal with character input
              IF(IPTYPE.EQ.IPCHAR) THEN !Type Character
                CALL CUPPER(CHDATA(1:4),C4)
C                IF(CUPPER(CHDATA(1:4)).EQ.'BACK') GO TO 9995
                IF(C4.EQ.'BACK') GO TO 9995
                nValues=nValues+1
                C_DATA(nValues)=CLINE(ISTART:IEND)
                IF(nValues.EQ.NDAT) THEN
                  ERROR=' '
                  GO TO 9998
                ENDIF
            !Check for one character control characters
              ELSE IF(ISTART.EQ.1.AND.IEND.EQ.1) THEN !one charac only
                CALL CUPPER(CHDATA(1:1),C1)
                IF(C1.EQ.'B') THEN
                  GO TO 9995
C KAT 6Sep00: This is not implemented
C                ELSE IF(C1.EQ.'A') THEN
C                  ERROR='Apply'
C                  GO TO 9001
                ELSE IF(C1.EQ.'D') THEN
                  ERROR='Default'
                  GO TO 9001
                ELSE IF(C1.EQ.'E') THEN
                  ERROR='Edit'
                  GO TO 9001
                ELSE IF(C1.EQ.'F') THEN
                  IF(.NOT.FILEIP) THEN
                    FILEIP=.TRUE.
                    FILEIP_FIRST=.TRUE.
                    WRITE(IUNIT,'($,'' Enter FILENAME: '')')
                    READ(IUNIT,'(A)') FILE07
                  ELSE
                    FILEIP_FIRST=.FALSE.
                  ENDIF
C                ELSE IF(CUPPER(CHDATA(1:1)).EQ.'G') THEN
                ELSE IF(C1.EQ.'G') THEN
                  IF(.NOT.FILEIP) THEN
                    FILEIP=.TRUE.
                    FILEIP_FIRST=.TRUE.
                    WRITE(IUNIT,'($,'' Enter FILENAME(.data): '')')
                    READ(IUNIT,'(A)') FILE07
                    CALL STRING_TRIM(FILE07,IBEG,IEND)
                    CALL OPENF(7,'DISK',FILE07(IBEG:IEND)//'.data',
     '                'OLD','SEQUEN','FORMATTED',132,ERROR,*9999)
                  ELSE
                    FILEIP_FIRST=.FALSE.
                  ENDIF
                  nValues=nValues+1
                  IF(IPTYPE.EQ.IPANSW) THEN
                    READ(7,*) A_DATA(nValues)
                  ELSE IF(IPTYPE.EQ.IPCHAR) THEN
                    READ(7,*) C_DATA(nValues)
                  ELSE IF(IPTYPE.EQ.IPINTE) THEN
                    READ(7,*) I_DATA(nValues)
                  ELSE IF(IPTYPE.EQ.IPLOGI) THEN
                    READ(7,*) L_DATA(nValues)
                  ELSE IF(IPTYPE.EQ.IPREAL) THEN
                    READ(7,*) R_DATA(nValues)
                  ENDIF
                  ERROR=' '
                  GO TO 9998
C                ELSE IF(CUPPER(CHDATA(1:1)).EQ.'Q') THEN
                ELSE IF(C1.EQ.'Q') THEN
                  CLOSE(UNIT=IUNIT)
                  STOP
C                ELSE IF(CUPPER(CHDATA(1:1)).EQ.'R') THEN
                ELSE IF(C1.EQ.'R') THEN
                  ERROR='Restart'
                  GO TO 9001
                ENDIF !cupper(chdata(1:1))=a,b,d,e,f,g,q or r
              ENDIF !character i/p or one control character only

            !Interpret character string for given IPTYPE
              IF(IPTYPE.EQ.IPANSW) THEN !Type Yes/No
                CALL CUPPER(CHDATA(1:1),C1)
C                IF(CUPPER(CHDATA(1:1)).EQ.'Y') THEN
                IF(C1.EQ.'Y') THEN
                  nValues=nValues+1
                  A_DATA(nValues)='Y'
                  IF(nValues.EQ.NDAT) THEN
                    ERROR=' '
                    GO TO 9998
                  ENDIF
C                ELSE IF(CUPPER(CHDATA(1:1)).EQ.'N') THEN
                ELSE IF(C1.EQ.'N') THEN
                  nValues=nValues+1
                  A_DATA(nValues)='N'
                  IF(nValues.EQ.NDAT) THEN
                    ERROR=' '
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Answer data range error'
                  GO TO 9998
                ENDIF

              ELSE IF(IPTYPE.EQ.IPINTE) THEN !Type Integer
C KAT 2005-12-24: This went bad w length > 16
C                 LENGTH=IEND-ISTART+1
C                 CHDATA=' '
C                 CHDATA(17-LENGTH:16)=CLINE(ISTART:IEND) !I16 format
C                 IF(IVALID(CHDATA)) THEN !valid integer data
C                   nValues=nValues+1
C                   READ(UNIT=CHDATA,FMT='(I16)',IOSTAT=IOSTAT)
C      '              I_DATA(nValues)
                IF(IVALID(CLINE(ISTART:IEND))) THEN !valid integer data
                  nValues=nValues+1
                  READ(UNIT=CLINE(ISTART:IEND),FMT=*,IOSTAT=IOSTAT)
     '              I_DATA(nValues)
                  write(*,*) I_DATA(nValues)
                  IF(IOSTAT.EQ.0) THEN
                    IF((I_DATA(nValues).GE.IMIN).AND.
     '                 (I_DATA(nValues).LE.IMAX)) THEN
C Potential bug???
C MPN 23Feb2000: Should really check against IOIMX (see 'inout00.cmn')
C                here, but this has not been 'include'd above.
C                Instead I_DATA has been passed through parameter lists.
C                Should IOIMX also be passed through also???
C KAT 25/2/00:        The size of I_DATA should depend on the C_DATA.  If C_DATA=
C                     'GRIDS' then the maximum dimension should be NGM.
C                     Perhaps it would be easiest to have a parameter.
                      CALL ASSERT(nValues.LE.9000,
     '                  '>>Too many integer values in list',
     '                  ERROR,*9999)
                      I_DATA(0)=nValues !is #integers read in
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' I_DATA(0)='',I3,'
     '                    //''' I_DATA('',I2,'')='',I5)')
     '                    I_DATA(0),nValues,I_DATA(nValues)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                    ELSE
                      ERROR=' Integer data range error'
                      GO TO 9998
                    ENDIF
                  ELSE
                    ERROR=' Integer data read error'
                    GO TO 9998
                  ENDIF

              !check for list in form FIRST..LAST:INCR
                ELSE IF(CLOCAT('..',CLINE(ISTART:IEND)).GT.0) THEN
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' CLINE(ISTART:IEND) is'',A)')
     '                CHDATA(1:16)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF !dop
                  CALL PARSIL(CLINE(ISTART:IEND),max_dim-nValues,NTIL,
     '              I_DATA(nValues+1),ERROR,*9999)
                  nValues=nValues+NTIL
                  I_DATA(0)=nValues
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' NTIL='',I4,'' nValues='',I4)')
     '                NTIL,nValues
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' I_DATA:'',10I5,/:(10I5))')
     '                (I_DATA(n),n=1,I_DATA(0))
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF !dop

              !check for group name
                ELSE
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' CLINE is '',A)')
     '                CLINE(istart:iend)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' C_DATA(1)='',A)')
     '                C_DATA(1)(1:10)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF !dop
                  FOUND_GROUP=.FALSE.
                  IF(C_DATA(1)(1:5).EQ.'NODES') THEN
                    !put group nodes into I_DATA(n)
                    DO nogrno=1,NTGRNO
                      IF(CLINE(istart:iend).EQ.
     '                  LAGRNO(nogrno)(1:iend-istart+1)) THEN
C TVK 06July1999: Dynamic Groups
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,'('' Group name is '',A)')
     '                      LAGRNO(nogrno)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                          WRITE(OP_STRING,'('' Group nodes are: '','
C     '                      //'10I3)')
C     '                      (LIGRNO(n,nogrno),n=1,LIGRNO(0,nogrno))
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF !dop
C                        DO n=1,LIGRNO(0,nogrno)
C                          I_DATA(nValues+n)=LIGRNO(n,nogrno)
C                        ENDDO
C                        nValues=nValues+LIGRNO(0,nogrno)
                        CALL ILIST_COPY(NLIGRNO(nogrno),
     '                    %VAL(LIGRNO_PTR(nogrno)),I_DATA(nValues+1))
                        nValues=nValues+NLIGRNO(nogrno)
                        I_DATA(0)=nValues
                        FOUND_GROUP=.TRUE.
                      ENDIF !cline
                    ENDDO !nogrno

                  ELSE IF(C_DATA(1)(1:8).EQ.'ELEMENTS') THEN
                    !put group elements into I_DATA(n)
                    DO nogrel=1,NTGREL
                      IF(CLINE(istart:iend).EQ.
     '                  LAGREL(nogrel)(1:iend-istart+1)) THEN
C TVK 14July1999: Dynamic Groups
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,'('' Group name is '',A)')
     '                      LAGREL(nogrel)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                          WRITE(OP_STRING,'('' Group elements are: '','
C     '                      //'10I3)')
C     '                      (LIGREL(n,nogrel),n=1,LIGREL(0,nogrel))
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF !dop
C                        DO n=1,LIGREL(0,nogrel)
C                          I_DATA(nValues+n)=LIGREL(n,nogrel)
C                        ENDDO
C                        nValues=nValues+LIGREL(0,nogrel)
                        CALL ILIST_COPY(NLIGREL(nogrel),
     '                    %VAL(LIGREL_PTR(nogrel)),I_DATA(nValues+1))
                        nValues=nValues+NLIGREL(nogrel)
                        I_DATA(0)=nValues
                        FOUND_GROUP=.TRUE.
                      ENDIF !cline
                    ENDDO !nogrel

                  ELSE IF(C_DATA(1)(1:5).EQ.'FACES') THEN
                    !put group faces into I_DATA(n)
                    DO nogrfa=1,NTGRFA
                      IF(CLINE(istart:iend).EQ.
     '                  LAGRFA(nogrfa)(1:iend-istart+1)) THEN
C TVK 09August1999: Dynamic Groups
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,'('' Group name is '',A)')
     '                      LAGRFA(nogrfa)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                          WRITE(OP_STRING,'('' Group faces are: '','
C     '                      //'10I3)')
C     '                      (LIGRFA(n,nogrfa),n=1,LIGRFA(0,nogrfa))
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF !dop
C                        DO n=1,LIGRFA(0,nogrfa)
C                          I_DATA(nValues+n)=LIGRFA(n,nogrfa)
C                        ENDDO
C                        nValues=nValues+LIGRFA(0,nogrfa)
                        CALL ILIST_COPY(NLIGRFA(nogrfa),
     '                    %VAL(LIGRFA_PTR(nogrfa)),I_DATA(nValues+1))
                        nValues=nValues+NLIGRFA(nogrfa)
                        I_DATA(0)=nValues
                        FOUND_GROUP=.TRUE.
                      ENDIF !cline
                    ENDDO !nogrfa

                  ELSE IF(C_DATA(1)(1:5).EQ.'LINES') THEN
                    !put group lines into I_DATA(n)
                    DO nogrli=1,NTGRLI
                      IF(CLINE(istart:iend).EQ.
     '                  LAGRLI(nogrli)(1:iend-istart+1)) THEN
C TVK 29July1999: Dynamic Groups
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,'('' Group name is '',A)')
     '                      LAGRLI(nogrli)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                          WRITE(OP_STRING,'('' Group lines are: '','
C     '                      //'10I3)')
C     '                      (LIGRLI(n,nogrli),n=1,LIGRLI(0,nogrli))
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF !dop
C                        DO n=1,LIGRLI(0,nogrli)
C                          I_DATA(nValues+n)=LIGRLI(n,nogrli)
C                        ENDDO
C                        nValues=nValues+LIGRLI(0,nogrli)
                        CALL ILIST_COPY(NLIGRLI(nogrli),
     '                    %VAL(LIGRLI_PTR(nogrli)),I_DATA(nValues+1))
                        nValues=nValues+NLIGRLI(nogrli)
                        I_DATA(0)=nValues
                        FOUND_GROUP=.TRUE.
                      ENDIF !cline
                    ENDDO !nogrli

                  ELSE IF(C_DATA(1)(1:5).EQ.'GRIDS') THEN
                    !put group grid points into I_DATA(n)
                    DO nogrgr=1,NTGRGR
                      IF(CLINE(istart:iend).EQ.
     '                  LAGRGR(nogrgr)(1:iend-istart+1)) THEN
C KAT 11May99: dynamic groups
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,'('' Group name is '',A)')
     '                      LAGRGR(nogrgr)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                          WRITE(OP_STRING,'('' Group grid pts are: '','
C     '                      //'10I3)')
C     '                      (LIGRGR(n,nogrgr),n=1,LIGRGR(0,nogrgr))
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF !dop
C                        DO n=1,LIGRGR(0,nogrgr)
C                          I_DATA(nValues+n)=LIGRGR(n,nogrgr)
C                        ENDDO
C                        nValues=nValues+LIGRGR(0,nogrgr)
C Potential bug???
C MPN 23Feb2000: Should really check against IOIMX (see 'inout00.cmn')
C                here, but this has not been 'include'd above.
C                Instead I_DATA has been passed through parameter lists.
C                Should IOIMX also be passed through also???
C KAT 25/2/00:          Should probably check something like
C                       nValues+NLIGRGR(nogrgr) <= NGM.
CC new MPN 23Feb2000: check for list overrun
C                        CALL ASSERT(NLIGRGR(nogrgr).LE.9000,
C     '                    '>>Too many integer values in list',
C     '                    ERROR,*9999)
CC new end
                        CALL ILIST_COPY(NLIGRGR(nogrgr),
     '                    %VAL(LIGRGR_PTR(nogrgr)),I_DATA(nValues+1))
                        nValues=nValues+NLIGRGR(nogrgr)
                        I_DATA(0)=nValues
                        FOUND_GROUP=.TRUE.
                      ENDIF !cline
                    ENDDO !nogrgr

                  ELSE IF(C_DATA(1)(1:5).EQ.'GAUSS') THEN
                    !put group gauss points into I_DATA(n)
                    DO nogrga=1,NTGRGA
                      IF(CLINE(istart:iend).EQ.
     '                  LAGRGA(nogrga)(1:iend-istart+1)) THEN
                        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                        call mp_setlock()
                          WRITE(OP_STRING,'('' Group name is '',A)')
     '                      LAGRGA(nogrga)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                        call mp_unsetlock()
                        ENDIF !dop
                        I_DATA(0)=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),1)
                        index=2
                        DO noelem=1,I_DATA(0)
                          I_DATA(index-1)=
     '                      ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index)
                          I_DATA(index)=
     '                      ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index+1)
                          DO num_points=1,I_DATA(index)
                          I_DATA(index+num_points)=
     '                      ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),
     '                      index+1+num_points)
                          ENDDO !num_points
                          index=index+2+num_points
                        ENDDO !noelem
                        nValues=nValues+NLIGRGA(nogrga)
                        FOUND_GROUP=.TRUE.
                      ENDIF !cline
                    ENDDO !nogrga
                  ENDIF !cdata
                  IF(.NOT.FOUND_GROUP) THEN
                    ERROR=' No group of that name'
                    GO TO 9998
                  ENDIF
                ENDIF !valid integer data or group name

              ELSE IF(IPTYPE.EQ.IPLOGI) THEN !Type Logical
                CALL CUPPER(CHDATA(1:1),C1)
C                IF(CUPPER(CHDATA(1:1)).EQ.'T') THEN
                IF(C1.EQ.'T') THEN
                  nValues=nValues+1
                  L_DATA(nValues)=.TRUE.
                  IF(nValues.EQ.NDAT) THEN
                    ERROR=' '
                    GO TO 9998
                  ENDIF
C                ELSE IF(CUPPER(CHDATA(1:1)).EQ.'F') THEN
                ELSE IF(C1.EQ.'F') THEN
                  nValues=nValues+1
                  L_DATA(nValues)=.FALSE.
                  IF(nValues.EQ.NDAT) THEN
                    ERROR=' '
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Logical data range error'
                  GO TO 9998
                ENDIF

              ELSE IF(IPTYPE.EQ.IPREAL) THEN !Type Real*8
C KAT 31/3/00: This only worked for up to 16 characters
C                LENGTH=IEND-ISTART+1
C                CHDATA=' '
C                CHDATA(17-LENGTH:16)=CLINE(ISTART:IEND)
C                IF(RVALID(CHDATA)) THEN
C                  nValues=nValues+1
C                  READ(UNIT=CHDATA,FMT='(G16.0)',IOSTAT=IOSTAT)
C     '              R_DATA(nValues)
                IF(RVALID(CLINE(ISTART:IEND))) THEN
                  nValues=nValues+1
                  READ(UNIT=CLINE(ISTART:IEND),FMT=*,IOSTAT=IOSTAT)
     '              R_DATA(nValues)
                  IF(IOSTAT.EQ.0) THEN
                    IF((R_DATA(nValues).GE.RMIN).AND.
     '                 (R_DATA(nValues).LE.RMAX)) THEN
                      IF(nValues.EQ.NDAT) THEN
                        ERROR=' '
                        GO TO 9998
                      ENDIF
                    ELSE
                      ERROR=' Real data range error'
                      GO TO 9998
                    ENDIF
                  ELSE
                    ERROR=' Real data read error'
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Real data type error'
                  GO TO 9998
                ENDIF
              ENDIF !iptype

            ELSE  !istart=0
              GOTO 4
            ENDIF !istart>0
          ENDDO !nd (i.e. index for list of data values in CLINE)

        ELSE IF(IOSTAT.EQ.0.AND.NCHAR.EQ.0) THEN !handle default
          CLINE=' '
          DO nd=1,NDAT
            !define CLINE as character string containing defaults
            CALL STRING_TRIM(CLINE,IBEG,IEND)
            IF(IPTYPE.EQ.IPANSW) THEN      !input is answer(yes/no)
              A_DATA(nd)=A_DEFLT(nd)
              CLINE=CLINE(IBEG:IEND)//A_DEFLT(nd)
            ELSE IF(IPTYPE.EQ.IPCHAR) THEN !input is character
              C_DATA(nd)=C_DEFLT(nd)
C new MPN 31-Jul-95: delimit strings by commas
              CLINE=CLINE(IBEG:IEND)//','//C_DEFLT(nd)
C old         CLINE=CLINE(IBEG:IEND)//C_DEFLT(nd)
            ELSE IF(IPTYPE.EQ.IPINTE) THEN !input is integer
              I_DATA(nd)=I_DEFLT(nd)
              WRITE(UNIT=CHAR10,FMT='(I10)') I_DEFLT(nd)
              CALL STRING_TRIM(CHAR10,IBEG2,IEND2)
              CLINE=CLINE(IBEG:IEND)//' '//CHAR10(IBEG2:IEND2)
            ELSE IF(IPTYPE.EQ.IPLOGI) THEN !input is logical
              L_DATA(nd)=L_DEFLT(nd)
              IF(L_DEFLT(nd)) THEN
                CLINE=CLINE(IBEG:IEND)//'T'
              ELSE
                CLINE=CLINE(IBEG:IEND)//'F'
              ENDIF
            ELSE IF(IPTYPE.EQ.IPREAL) THEN !input is real
              R_DATA(nd)=R_DEFLT(nd)
              WRITE(CHARTEMP,'(D13.5)') R_DEFLT(nd)
              CLINE=CLINE(IBEG:IEND)//CHARTEMP
            ENDIF
          ENDDO !nd
          nValues=nValues+NDAT
          IF(IPTYPE.EQ.IPCHAR) THEN !remove leading comma for character ip
            CALL STRING_TRIM(CLINE,IBEG,IEND)
            IF(CLINE(IBEG:IBEG).EQ.',') CLINE=CLINE(IBEG+1:IEND)
          ENDIF

        ELSE
          REWIND(UNIT=IUNIT)
          ERROR=' Character data read error'
          GO TO 9998
        ENDIF
    4   CONTINUE
      ENDDO !while nValues<NDAT

      ERROR=' '
 9998 CALL EXITS('GINPUT')
      RETURN

 9995 ERROR=' '
      CALL EXITS('GINPUT')
      RETURN 2

 9999 CALL ERRORS('GINPUT',ERROR)
 9001 CALL EXITS('GINPUT')
      RETURN 1
      END


