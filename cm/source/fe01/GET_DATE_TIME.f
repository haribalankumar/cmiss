      SUBROUTINE GET_DATE_TIME(DATE,IDATEFMT,ERROR,*)

C#### Subroutine: GET_DATE_TIME
C###  Description:
C###    <HTML> <PRE>
C###    GET_DATE_TIME uses runtime library call to return date and time
C###    with IDATEFMT equal to:
C###    1  -  returns date in DEC format "DD-MMM-YYYY HH:MM:SS.HH"
C###          ie 3:45 in the afternoon on 27th of November 1990
C###          becomes "27-NOV-1990 15:45:00.00". (Default).
C###    2  -  returns date as all numbers in IGES format "YYMMDD.HHMMSS"
C###          ie the same date is returned as "901127.154500"
C###    3  -  returns a 7 or 8 digit random number based on the time
C###          as "HHSSMMHH" working from hundredths back to hours
C###    </PRE> </HTML>

      IMPLICIT NONE
!     Parameter List
      INTEGER IDATEFMT
      CHARACTER DATE*(*),ERROR*(*)
!     Local Variables
c1      INTEGER MONTHNUM
      CHARACTER MONTH*(3),TEMPDATE*(30)

      CALL ENTERS('GET_DATE_TIME',*9999)

!old  CALL LIB$DATE_TIME(TEMPDATE)
      TEMPDATE='01-JAN-1992 00:00:00.00'
      IF (IDATEFMT.EQ.2) THEN
        MONTH=TEMPDATE(4:6)
        IF (MONTH.EQ.'JAN') THEN
          MONTH='01'
        ELSE IF (MONTH.EQ.'FEB') THEN
          MONTH='02'
        ELSE IF (MONTH.EQ.'MAR') THEN
          MONTH='03'
        ELSE IF (MONTH.EQ.'APR') THEN
          MONTH='04'
        ELSE IF (MONTH.EQ.'MAY') THEN
          MONTH='05'
        ELSE IF (MONTH.EQ.'JUN') THEN
          MONTH='06'
        ELSE IF (MONTH.EQ.'JUL') THEN
          MONTH='07'
        ELSE IF (MONTH.EQ.'AUG') THEN
          MONTH='08'
        ELSE IF (MONTH.EQ.'SEP') THEN
          MONTH='09'
        ELSE IF (MONTH.EQ.'OCT') THEN
          MONTH='10'
        ELSE IF (MONTH.EQ.'NOV') THEN
          MONTH='11'
        ELSE IF (MONTH.EQ.'DEC') THEN
          MONTH='12'
        ENDIF
        DATE=TEMPDATE(10:11)//MONTH(1:2)//TEMPDATE(1:2)//'.'//
     '       TEMPDATE(13:14)//TEMPDATE(16:17)//TEMPDATE(19:20)
      ELSE IF (IDATEFMT.EQ.3) THEN
        DATE=TEMPDATE(22:23)//TEMPDATE(19:20)//TEMPDATE(16:17)
     '    //TEMPDATE(13:14)
      ELSE
        DATE=TEMPDATE
      ENDIF

      CALL EXITS('GET_DATE_TIME')
      RETURN
 9999 CALL ERRORS('GET_DATE_TIME',ERROR)
      CALL EXITS('GET_DATE_TIME')
      RETURN 1
      END


