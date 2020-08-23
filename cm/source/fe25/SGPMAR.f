      SUBROUTINE SGPMAR(INDEX,ISEG,ISPMAR,iw,NTPTS,CSEG,XLIST,YLIST,
     '  ERROR,*)

C#### Subroutine: SGPMAR
C###  Description:
C###    SGPMAR creates polymarker segment ISPMAR.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISPMAR,iw,NTPTS
      REAL*8 XLIST(*),YLIST(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX_OLD,MAXPTS,N,nopt
      PARAMETER (MAXPTS=100)  !note limit of 100 points
      REAL*8 PTS(3,MAXPTS)

      CALL ENTERS('SGPMAR',*9999)
      CALL OPEN_SEGMENT(ISPMAR,ISEG,iw,'PMAR',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' X_list: '',10E11.3)')
     '    (XLIST(N),N=1,NTPTS)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Y_list: '',10E11.3)')
     '    (YLIST(N),N=1,NTPTS)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(NTPTS.GT.MAXPTS) THEN
        ERROR='Too many points in SGPLIN'
        GOTO 9999
      ENDIF
      DO nopt=1,NTPTS
        PTS(1,nopt)=XLIST(nopt)
        PTS(2,nopt)=YLIST(nopt)
      ENDDO
      CALL POLYMARKER(INDEX,iw,NTPTS,PTS,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISPMAR,iw,ERROR,*9999)

      CALL EXITS('SGPMAR')
      RETURN
 9999 CALL ERRORS('SGPMAR',ERROR)
      CALL EXITS('SGPMAR')
      RETURN 1
      END


