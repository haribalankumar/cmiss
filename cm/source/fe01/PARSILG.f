      SUBROUTINE PARSILG(ILIST,ILISTM,CDATA,CLINE,ERROR,*)

C#### Subroutine: PARSILG
C###  Description:
C###    PARSILG handles input of a sequence of integer data including
C###    group names in the character string CLINE.
C###    Set CDATA to 'NODES','ELEMENTS','FACES','GRIDS','GAUSS'
C###    or 'LINES' to pick up appropriate group name.
C###    Values are output in ILIST(n),n=1,ILIST(0).    
C***    Rewritten MPN 4-4-96

C***    Note that there is no longer any checking for duplicates here
C***    (was not implemented for non-group lists anyway). Now this
C***    should be done in the calling routine through ILISTRMDUP. Also
C***    note the setting of CDATA to ' ' at the end of the routine has
C***    been removed.
C***    Re-rewritten MLB 19-2-2004

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'grou00.cmn'

!     Parameter List
      INTEGER ILISTM,ILIST(0:ILISTM)
      CHARACTER CDATA*(*),CLINE*(*),ERROR*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG1,IBEG2,IBEG3,IEND1,IEND2,IEND3,ILIST_START,
     &  IPOS1,IPOS2,N1GR,nogr,NUMINLIST
      CHARACTER CHAR1*(GR_MAXNAME),CHARAC*1,LABEL*(GR_MAXNAME),TAB*1

      CALL ENTERS('PARSILG',*9999)

      TAB=CHAR(9) !ASCII HT (horizontal TAB)
      CALL STRING_TRIM(CLINE,IBEG1,IEND1)
      !Check for null string
      IF(IEND1.EQ.0) THEN !no data
        ILIST(0)=0
        ILIST(1)=0
        GOTO 100
      ENDIF
      IPOS1=IBEG1
      ILIST(0)=0 !current total number of integers in list
      DO WHILE(IPOS1.LE.IEND1)
        CHARAC=CLINE(IPOS1:IPOS1)
        IF((CHARAC.GE.'0'.AND.CHARAC.LE.'9').OR. !a digit
     &    CHARAC.EQ.'-'.OR. !a minus sign
     &    CHARAC.EQ.'+'.OR. !a plus sign
     &    CHARAC.EQ.' ') THEN !a white space
          !A list of #'s follows
          IPOS2=IPOS1
          CHARAC=CLINE(IPOS2:IPOS2)
          DO WHILE(((CHARAC.GE.'0'.AND.CHARAC.LE.'9').OR. !a digit
     &      CHARAC.EQ.'.'.OR. !a period
     &      CHARAC.EQ.','.OR. !a comma
     &      CHARAC.EQ.':'.OR. !a colon
     &      CHARAC.EQ.'-'.OR. !a minus sign
     &      CHARAC.EQ.'+'.OR. !a plus sign
     &      CHARAC.EQ.TAB.OR. !a horizontal tab
     &      CHARAC.EQ.' ').AND. !a white space
     &      IPOS2.LT.IEND1)
            IPOS2=IPOS2+1
            CHARAC=CLINE(IPOS2:IPOS2)
          ENDDO
          IF(IPOS2.NE.IEND1) IPOS2=IPOS2-2
          !Parse and add current list of integers to list
          CALL PARSIL(CLINE(IPOS1:IPOS2),ILISTM-ILIST(0),
     &      NUMINLIST,ILIST(ILIST(0)+1),ERROR,*9999)
          ILIST(0)=ILIST(0)+NUMINLIST

        ELSE !possibly a group name entered
          IPOS2=IPOS1+CLOCAT(',',CLINE(IPOS1:IEND1))
          IF(IPOS2.EQ.IPOS1) THEN
            IPOS2=IEND1
          ELSE
            IPOS2=IPOS2-2
          ENDIF
          IF(IPOS2-IPOS1+1.GT.LEN(CHAR1)) THEN
            ! name greater than length of char1: can't use cupper,
            ! but it doesn't match any group name because they are shorter.
            N1GR=-1
          ELSE
            CALL CUPPER(CLINE(IPOS1:IPOS2),CHAR1)
            CALL STRING_TRIM(CHAR1,IBEG2,IEND2)
            N1GR=0
          ENDIF

          IF(CDATA(1:5).EQ.'NODES') THEN
            !Find which group has been requested
            nogr=0
            DO WHILE(nogr.LT.NTGRNO.AND.N1GR.EQ.0)
              nogr=nogr+1
              CALL CUPPER(LAGRNO(nogr),LABEL)
              CALL STRING_TRIM(LABEL,IBEG3,IEND3)
              IF(CHAR1(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) N1GR=nogr
            ENDDO !while(nogr.LT.NTGRNO.AND.N1GR.EQ.0)

            IF(N1GR.LE.0) THEN
              !No group found
              ERROR='>>"'//CLINE(IPOS1:IPOS2)//'" is not a node group'
              GO TO 9999
            ELSE
              !Check ILIST is big enough
              CALL ASSERT(ILIST(0)+NLIGRNO(N1GR).LE.ILISTM,
     &          '>>ERROR: Node list exceeds storage',
     &          ERROR,*9999)

              !Copy the current group into ILIST
              ILIST_START=ILIST(0)+1
              CALL ILIST_COPY(NLIGRNO(N1GR),%VAL(LIGRNO_PTR(N1GR)),
     &          ILIST(ILIST_START))
              ILIST(0)=ILIST(0)+NLIGRNO(N1GR)
            ENDIF

          ELSE IF(CDATA(1:8).EQ.'ELEMENTS') THEN
            !Find which group has been requested
            nogr=0
            DO WHILE(nogr.LT.NTGREL.AND.N1GR.EQ.0)
              nogr=nogr+1
              CALL CUPPER(LAGREL(nogr),LABEL)
              CALL STRING_TRIM(LABEL,IBEG3,IEND3)
              IF(CHAR1(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) N1GR=nogr
            ENDDO !while(nogr.LT.NTGREL.AND.N1GR.EQ.0)

            IF(N1GR.LE.0) THEN
              !No group found
              ERROR='>>"'//CLINE(IPOS1:IPOS2)
     &          //'" is not an element group'
              GO TO 9999
            ELSE
              !Check ILIST is big enough
              CALL ASSERT(ILIST(0)+NLIGREL(N1GR).LE.ILISTM,
     &          '>>ERROR: Node list exceeds storage',
     &          ERROR,*9999)

              !Copy the current group into ILIST
              ILIST_START=ILIST(0)+1
              CALL ILIST_COPY(NLIGREL(N1GR),%VAL(LIGREL_PTR(N1GR)),
     &          ILIST(ILIST_START))
              ILIST(0)=ILIST(0)+NLIGREL(N1GR)
            ENDIF

          ELSE IF(CDATA(1:5).EQ.'FACES') THEN
            !Find which group has been requested
            nogr=0
            DO WHILE(nogr.LT.NTGRFA.AND.N1GR.EQ.0)
              nogr=nogr+1
              CALL CUPPER(LAGRFA(nogr),LABEL)
              CALL STRING_TRIM(LABEL,IBEG3,IEND3)
              IF(CHAR1(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) N1GR=nogr
            ENDDO !while(nogr.LT.NTGRFA.AND.N1GR.EQ.0)

            IF(N1GR.LE.0) THEN
              !No group found
              ERROR='>>"'//CLINE(IPOS1:IPOS2)//'" is not a face group'
              GO TO 9999
            ELSE
              !Check ILIST is big enough
              CALL ASSERT(ILIST(0)+NLIGRFA(N1GR).LE.ILISTM,
     &          '>>ERROR: Node list exceeds storage',
     &          ERROR,*9999)

              !Copy the current group into ILIST
              ILIST_START=ILIST(0)+1
              CALL ILIST_COPY(NLIGRFA(N1GR),%VAL(LIGRFA_PTR(N1GR)),
     &          ILIST(ILIST_START))
              ILIST(0)=ILIST(0)+NLIGRFA(N1GR)
            ENDIF

          ELSE IF(CDATA(1:5).EQ.'GAUSS') THEN
            !Find which group has been requested
            nogr=0
            DO WHILE(nogr.LT.NTGRGA.AND.N1GR.EQ.0)
              nogr=nogr+1
              CALL CUPPER(LAGRGA(nogr),LABEL)
              CALL STRING_TRIM(LABEL,IBEG3,IEND3)
              IF(CHAR1(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) N1GR=nogr
            ENDDO !while(nogr.LT.NTGRGA.AND.N1GR.EQ.0)

            IF(N1GR.LE.0) THEN
              !No group found
              ERROR='>>"'//CLINE(IPOS1:IPOS2)//'" is not a gauss group'
              GO TO 9999
            ELSE
              !Check ILIST is big enough
              CALL ASSERT(ILIST(0)+NLIGRGA(N1GR).LE.ILISTM,
     &          '>>ERROR: Gauss list exceeds storage',
     &          ERROR,*9999)

              !Copy the current group into ILIST
              ILIST_START=ILIST(0)+1
              CALL ILIST_COPY(NLIGRGA(N1GR),%VAL(LIGRGA_PTR(N1GR)),
     &          ILIST(ILIST_START))
              ILIST(0)=ILIST(0)+NLIGRGA(N1GR)
            ENDIF

          ELSE IF(CDATA(1:5).EQ.'GRIDS') THEN
            !Find which group has been requested
            nogr=0
            DO WHILE(nogr.LT.NTGRGR.AND.N1GR.EQ.0)
              nogr=nogr+1
              CALL CUPPER(LAGRGR(nogr),LABEL)
              CALL STRING_TRIM(LABEL,IBEG3,IEND3)
              IF(CHAR1(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) N1GR=nogr
            ENDDO !while(nogr.LT.NTGRGR.AND.N1GR.EQ.0)

            IF(N1GR.LE.0) THEN 
              !No group found
              ERROR='>>"'//CLINE(IPOS1:IPOS2)
     &          //'" is not a grid point group'
              GO TO 9999
            ELSE
              !Check ILIST is big enough
              CALL ASSERT(ILIST(0)+NLIGRGR(N1GR).LE.ILISTM,
     &          '>>ERROR: Grid point list exceeds storage',
     &          ERROR,*9999)

              !Copy the current group into ILIST
              ILIST_START=ILIST(0)+1
              CALL ILIST_COPY(NLIGRGR(N1GR),%VAL(LIGRGR_PTR(N1GR)),
     &          ILIST(ILIST_START))
              ILIST(0)=ILIST(0)+NLIGRGR(N1GR)
            ENDIF

          ELSE IF(CDATA(1:5).EQ.'LINES') THEN
            !Find which group has been requested
            nogr=0
            DO WHILE(nogr.LT.NTGRLI.AND.N1GR.EQ.0)
              nogr=nogr+1
              CALL CUPPER(LAGRLI(nogr),LABEL)
              CALL STRING_TRIM(LABEL,IBEG3,IEND3)
              IF(CHAR1(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) N1GR=nogr
            ENDDO !while(nogr.LT.NTGRLI.AND.N1GR.EQ.0)

            IF(N1GR.LE.0) THEN
              !No group found
              ERROR='>>"'//CLINE(IPOS1:IPOS2)//'" is not a line group'
              GO TO 9999
            ELSE
              !Check ILIST is big enough
              CALL ASSERT(ILIST(0)+NLIGRLI(N1GR).LE.ILISTM,
     &          '>>ERROR: Grid point list exceeds storage',
     &          ERROR,*9999)

              !Copy the current group into ILIST
              ILIST_START=ILIST(0)+1
              CALL ILIST_COPY(NLIGRLI(N1GR),%VAL(LIGRLI_PTR(N1GR)),
     &          ILIST(ILIST_START))
              ILIST(0)=ILIST(0)+NLIGRLI(N1GR)
            ENDIF

          ELSE IF(CDATA.EQ.' ') THEN
            ERROR='Data not integer: '//CLINE(IPOS1:IEND1)
            GOTO 9999

          ELSE
            CALL STRING_TRIM(CDATA,IBEG1,IEND1)
            WRITE(ERROR,'(''>>Invalid Option '',A)') CDATA(IBEG1:IEND1)
            GOTO 9999
          ENDIF !cdata=NODES/ELEMENTS/FACES/GRIDS/LINES
        ENDIF
        IPOS1=IPOS2+2
      ENDDO !while(IPOS1.LE.IEND)

 100  CONTINUE

      CALL EXITS('PARSILG')
      RETURN
 9999 CALL ERRORS('PARSILG',ERROR)
      CALL EXITS('PARSILG')
      RETURN 1
      END
