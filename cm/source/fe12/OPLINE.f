      SUBROUTINE OPLINE(NEL,NLLINE,NPL,DL,ERROR,*)

C#### Subroutine: OPLINE
C###  Description:
C###    OPLINE outputs line segments.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NEL(0:NELM,NLM),NLLINE(0:NL_R_M,0:NRM),NPL(5,0:3,NLM)
      REAL*8 DL(3,NLM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ILT,l,nl,nj,nr

      CALL ENTERS('OPLINE',*9999)

      DO nr=1,NRT
        WRITE(OP_STRING,'(/''  Number of lines in region '',I1,'
     '    //''' is '',I6)') nr,NLLINE(0,nr)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      IF(NLT.GT.0) THEN
        WRITE(OP_STRING,'(/1X,'' Global line segment data : The '','
     '    //'''total no. of line segments = '',I5/,3X,''Line  Xi  '','
     '    //'''Basis types'',2X,''Deriv.1'',5X,''Deriv.2'',5X,'
     '    //'''Length'',6X,''Nodes and derivs   '',3X,''Elements'')')
     '    NLT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nl=1,NLT
          IF(NPL(1,1,nl).EQ.2.OR.NPL(1,2,nl).EQ.2
     '      .OR.NPL(1,3,nl).EQ.2) THEN
            ILT=3
            FORMAT=
     '        '(2X,I5,3X,I1,3X,3I2,2X,3(1X,E11.4),1X,3I5, 8X,10I4)'
          ELSE IF(NPL(1,1,nl).GE.3.OR.NPL(1,2,nl).GE.3
     '      .OR.NPL(1,3,nl).GE.3) THEN
            ILT=4
            IF(JTYP2B.EQ.1) THEN
              IF(NPL(4,0,nl).GT.0) THEN
                FORMAT='(2X,I5,3X,I1,3X,3I2,2X,3(1X,E11.4),1X,4I5,3X,'
     '            //'2I4,'' - coupled to line '',I4)'
              ELSE IF(NPL(4,0,nl).LT.0) THEN
                FORMAT='(2X,I5,3X,I1,3X,3I2,2X,3(1X,E11.4),1X,4I5,3X,'
     '            //'2I4,'' - mapped to line  '',I4)'
              ELSE
                FORMAT='(2X,I5,3X,I1,3X,3I2,2X,3(1X,E11.4),1X,4I5,3X,'
     '            //'10I4)'
              ENDIF
            ELSE
              FORMAT=
     '          '(2X,I5,3X,I1,3X,3I2,2X,3(1X,E11.4),1X,4I5,3X,10I4)'
            ENDIF
          ELSE
            ILT=2
            FORMAT='(2X,I5,3X,I1,3X,3I2,2X,3(1X,E11.4),1X,2I5,13X,10I4)'
          ENDIF

          IF(JTYP2B.EQ.1) THEN
            IF(NPL(4,0,nl).NE.0) THEN
              WRITE(OP_STRING,FORMAT) nl,NPL(1,0,nl),(NPL(1,nj,nl),
     '          nj=1,3),(DL(i,nl),i=1,3),(NPL(1+l,1,nl),l=1,ILT),
     '          (NEL(l,nl),l=1,NEL(0,nl)),ABS(NPL(4,0,nl))
            ELSE
              WRITE(OP_STRING,FORMAT) nl,NPL(1,0,nl),(NPL(1,nj,nl),
     '          nj=1,3),(DL(i,nl),i=1,3),(NPL(1+l,1,nl),l=1,ILT),
     '          (NEL(l,nl),l=1,NEL(0,nl))
            ENDIF
          ELSE
            WRITE(OP_STRING,FORMAT) nl,NPL(1,0,nl),
     '        (NPL(1,nj,nl),nj=1,3),
     '        (DL(i,nl),i=1,3),(NPL(1+l,1,nl),l=1,ILT),
     '        (NEL(l,nl),l=1,NEL(0,nl))
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('OPLINE')
      RETURN
 9999 CALL ERRORS('OPLINE',ERROR)
      CALL EXITS('OPLINE')
      RETURN 1
      END


