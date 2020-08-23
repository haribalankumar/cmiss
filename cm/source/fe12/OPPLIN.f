      SUBROUTINE OPPLIN(ERROR,*)

C#### Subroutine: OPPLIN
C###  Description:
C###    OPPLIN outputs polyline parameters.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'plin00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,no_plin,no_point,NT_POINT
      CHARACTER TITLE(2)*16

      DATA TITLE/'Piecewise linear',
     '           'Bezier curve    '/

      CALL ENTERS('OPPLIN',*9999)
      WRITE(OP_STRING,'('' Number of polylines defined is '',I2)')
     '  NT_PLIN
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO no_plin=1,NT_PLIN
        WRITE(OP_STRING,'(''('',I2,'') Type = '',A,'
     '    //'''   Number of sections = '',I2,'' Segment #'',I3)')
     '    no_plin,TITLE(INDEX_PLIN_TYPE(no_plin)),
     '    NT_PLIN_SECTIONS(no_plin),
     '    NOPLIN_INDEX(no_plin)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(INDEX_PLIN_TYPE(no_plin).EQ.1) THEN      !piecewise linear
          NT_POINT=NT_PLIN_SECTIONS(no_plin)+1
        ELSE IF(INDEX_PLIN_TYPE(no_plin).EQ.2) THEN !Bezier curve
          NT_POINT=3*NT_PLIN_SECTIONS(no_plin)+1
        ENDIF
        DO no_point=1,NT_POINT
          WRITE(OP_STRING,'(''     Rect. cart. coords of point '',I2,'
     '      //''':  '',3E12.3)') no_point,
     '      (PLIN_DATA(nj,no_point,no_plin),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      CALL EXITS('OPPLIN')
      RETURN
 9999 CALL ERRORS('OPPLIN',ERROR)
      CALL EXITS('OPPLIN')
      RETURN 1
      END


