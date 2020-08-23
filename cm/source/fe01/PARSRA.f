      SUBROUTINE PARSRA(STRING,NXLV,NTLV,NTITLV,NXRA,RA,ERROR,*)

C#### Subroutine: PARSRA
C###  Description:
C###    PARSRA breaks a nested list of items found in STRING
C###    into separate REAL*8 numbers stored in the REAL*8 array RA.
C**** The maximum level of nesting is returned as NTLV.
C**** The size of each dimension is carried in the first NTLV
C**** components of the integer array NTITLV.
C**** Each nesting is delimited by left and right square brackets
C**** while each item is separated by commas.

      IMPLICIT NONE
!     Parameter List
      INTEGER NTLV,NXLV,NXRA,NTITLV(NXLV)
      REAL*8 RA(NXRA)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER MXLV
      PARAMETER (MXLV=16)
      INTEGER IBEG,IEND,n1ch,n2ch,NOITLV(MXLV),nolv,nora,NTRA
      CHARACTER CHAR

      CALL ENTERS('PARSRA',*9999)
      CALL STRING_TRIM(STRING,IBEG,IEND)
      DO nolv=1,NXLV
        NTITLV(nolv)=-1
      ENDDO
      NTLV=-1
      nolv=0
      NTRA=0
      DO n2ch=IBEG,IEND+1
        CHAR=STRING(n2ch:n2ch)
        IF(CHAR.EQ.'[') THEN
          nolv=nolv+1
          IF((nolv.GT.NXLV).OR.(nolv.GT.MXLV)) THEN
            ERROR=' Maximum array nesting exceeded'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=0
          n1ch=n2ch
        ELSE IF(CHAR.EQ.']') THEN
          IF(NTLV.LT.0) THEN
            NTLV=nolv
          ENDIF
          IF(nolv.EQ.NTLV) THEN
            CALL PARSRL(STRING(n1ch+1:n2ch-1),NXRA-NTRA,nora,RA(NTRA+1),
     '        ERROR,*9999)
            NOITLV(nolv)=nora
            NTRA=NTRA+nora
          ENDIF
          IF(NTITLV(nolv).LT.0) THEN
            NTITLV(nolv)=NOITLV(nolv)
          ELSE IF(NTITLV(nolv).NE.NOITLV(nolv)) THEN
            ERROR=' Inconsistent number of items in array'
            GOTO 9999
          ENDIF
          nolv=nolv-1
          IF(nolv.LT.0) THEN
            ERROR=' Inconsistent nesting of array'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=NOITLV(nolv)+1
        ENDIF
      ENDDO
      IF(nolv.EQ.0) THEN
        NTITLV(1)=NOITLV(1)
      ELSE
        ERROR=' Inconsistent nesting of array'
        GOTO 9999
      ENDIF

      CALL EXITS('PARSRA')
      RETURN
 9999 CALL ERRORS('PARSRA',ERROR)
      CALL EXITS('PARSRA')
      RETURN 1
      END


