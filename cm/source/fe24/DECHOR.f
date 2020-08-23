      SUBROUTINE DECHOR(NKJ,NLCHOR,NPL,DL,XP,STRING,ERROR,*)

C#### Subroutine: DECHOR
C###  Description:
C###    DECHOR defines chords.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NKJ(NJM,NPM),NLCHOR(0:10,NRM),NPL(5,0:3,NLM)
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER IBEG,IEND,INSTAT,NCOY,NCOZ,NTITLV(1),nora,NTLV,NTRA
      REAL*8 RA(100),XWC1,YWC1,ZPOS
      CHARACTER X1STR*132,X2STR*132,XSTR1*132,XSTR2*132,ZCHAR*1
      LOGICAL ABBREV

      CALL ENTERS('DECHOR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define chord
C###  Description:
C###    Define sail chords shapes to specify the shape of a sail.
C###    See Fiona McPheat's thesis on Finite Element Analysis of Yacht Sails.
C###  Parameter:      <from X_1#[0]> <to X_2#[0]>
C###    Specify the positions in the forward to aft direction.
C###  Parameter:      <at z=Z_LIST[locate]>
C###    Specify positions up the sail.
C###  Parameter:      <with y=Y_1#;Y_2#;..[0;0;..]>
C###    Specify the sail shape as a displacement from the chord.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<from X_1#[0]> <to X_2#[0]>'
        OP_STRING(3)=BLANK(1:15)//'<at z=Z_LIST[locate]>'
        OP_STRING(4)=BLANK(1:15)//'<with y=Y_1#;Y_2#;..[0;0;..]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECHOR',ERROR,*9999)
      ELSE IF(NJT.EQ.3) THEN
        CALL ACWK(1,0,ERROR,*9999)
        IF(ABBREV(CO(noco+1),'FROM',1)) THEN
          CALL STRING_TRIM(CO(noco+2),IBEG,IEND)
          XSTR1=CO(noco+2)(IBEG:IEND)
        ELSE
          XSTR1='0'
        ENDIF
        IF(ABBREV(CO(noco+1),'TO',1)) THEN
          CALL STRING_TRIM(CO(noco+2),IBEG,IEND)
          XSTR2=CO(noco+2)(IBEG:IEND)
        ELSE IF(ABBREV(CO(noco+3),'TO',1)) THEN
          CALL STRING_TRIM(CO(noco+4),IBEG,IEND)
          XSTR2=CO(noco+4)(IBEG:IEND)
        ELSE
          XSTR2='0'
        ENDIF
        IF(ABBREV(CO(noco+3),'AT',1)) THEN
          NCOZ=noco+4
        ELSE IF(ABBREV(CO(noco+5),'AT',1)) THEN
          NCOZ=noco+6
        ELSE
          NCOZ=0
        ENDIF
        IF(ABBREV(CO(noco+3),'WITH',1)) THEN
          NCOY=noco+4
        ELSE IF(ABBREV(CO(noco+5),'WITH',1)) THEN
          NCOY=noco+6
        ELSE IF(ABBREV(CO(noco+6),'WITH',1)) THEN
          NCOY=noco+7
        ELSE IF(ABBREV(CO(noco+8),'WITH',1)) THEN
          NCOY=noco+9
        ELSE
          NCOY=0
        ENDIF
        IF(NCOZ.GT.0) THEN
          CALL STRING_TRIM(CO(NCOZ),IBEG,IEND)
          ZCHAR=CO(NCOZ)(IBEG:IEND)
          CALL PARSTR(CO(NCOZ+1),1,NTLV,NTITLV,100,RA,ERROR,*9999)
          NTRA=NTITLV(1)
          IF(DOP) THEN
            WRITE(OP_STRING,'(1X,A)') 'ntra= ',NTRA
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DO nora=1,NTRA
            ZPOS=RA(nora)
            IF(DOP) THEN
              WRITE(OP_STRING,'(1X,A)') 'ZPOS= ',ZPOS
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            X1STR=XSTR1
            X2STR=XSTR2
            CALL CHORD(NCOY,NKJ,NLCHOR,NPL,DL,
     '        XP,ZPOS,STRING,X1STR,X2STR,ZCHAR,ERROR,*9999)
          ENDDO
        ELSE IF(NCOZ.EQ.0) THEN
          CALL LOCATOR(INSTAT,0.0D0,XWC1,0.0D0,YWC1,
     '      ERROR,*9999)
          DO WHILE(INSTAT.EQ.1)
            ZPOS=YWC1
            X1STR=XSTR1
            X2STR=XSTR2
            CALL CHORD(NCOY,NKJ,NLCHOR,NPL,DL,
     '        XP,ZPOS,STRING,X1STR,X2STR,ZCHAR,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.0D0,XWC1,0.0D0,YWC1,
     '        ERROR,*9999)
          ENDDO
        ENDIF
        CALL DAWK(1,0,ERROR,*9999)
      ENDIF

C MHT 04-05-01 no path to 9998
C  9998 CALL EXITS('DECHOR')
C      RETURN

      CALL EXITS('DECHOR')
      RETURN
 9999 CALL ERRORS('DECHOR',ERROR)
      CALL EXITS('DECHOR')
      RETURN 1
      END


