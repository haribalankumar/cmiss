      SUBROUTINE SGPLOTXYLABEL(IBEG,IEND,IPTYPE,RDATA,FORMAT_STR,
     '  STRING,ERROR,*)

C#### Subroutine SGPLOTXYLABEL
C###  Description:
C###    Auxiliary routine for SGPLOTXY. Create the axes-scales.
CC JMB 10-APR-2000

      IMPLICIT NONE

!     Parameter List
      INTEGER IBEG,IEND,IPTYPE
      REAL*8 RDATA
      CHARACTER FORMAT_STR,STRING*(*),ERROR*(*)

!     Local Variables
      INTEGER i,IFEND

      CALL ENTERS('SGPLOTXYLABEL',*9999)

      IF(FORMAT_STR.EQ.'E') THEN
        IF(IPTYPE.EQ.1) THEN
          ! Log-scale
          WRITE(STRING,'(I3)') INT(RDATA)
          CALL STRING_TRIM(STRING,IBEG,IEND)
          IF(INT(RDATA).GE.0) THEN
            STRING='e+'//STRING(IBEG:IEND)
          ELSE
            STRING='e'//STRING(IBEG:IEND)
          ENDIF
          CALL STRING_TRIM(STRING,IBEG,IEND)
        ELSEIF(IPTYPE.EQ.2) THEN
          WRITE(STRING,'(E10.4)') RDATA
          CALL STRING_TRIM(STRING,IBEG,IEND)
         IFEND=IEND-4
          DO i=IFEND,IBEG,-1
            IF(STRING(i:i).EQ.'0') THEN
               IFEND=IFEND-1
            ELSE
              GOTO 1
            ENDIF
          ENDDO
 1        STRING=STRING(IBEG:IFEND)//STRING(IEND-3:IEND)
          CALL STRING_TRIM(STRING,IBEG,IEND)
        ENDIF
      ELSEIF(FORMAT_STR.EQ.'F') THEN
        WRITE(STRING,'(F10.6)') RDATA
        CALL STRING_TRIM(STRING,IBEG,IEND)
        DO i=IEND,IBEG,-1
          IF(STRING(i:i).EQ.'0') THEN
            IEND=IEND-1
          ELSEIF(STRING(i:i).EQ.'.') THEN
            IEND=IEND-1
            GOTO 2
          ELSE
            GOTO 2
          ENDIF
        ENDDO
      ELSEIF(FORMAT_STR.EQ.'I') THEN
        WRITE(STRING,'(I4)') INT(RDATA)
        CALL STRING_TRIM(STRING,IBEG,IEND)
      ENDIF

 2    CALL EXITS('SGPLOTXYLABEL')
      RETURN
 9999 CALL ERRORS('SGPLOTXYLABEL',ERROR)
      CALL EXITS('SGPLOTXYLABEL')
      RETURN 1
      END

