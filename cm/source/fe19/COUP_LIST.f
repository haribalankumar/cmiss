      SUBROUTINE COUP_LIST(FSTRING,MAXOPT,COUP_TYPE,ERROR,*)

C#### Subroutine: COUP_LIST
C###  Description:
C###    COUP_LIST is the list of coupling models
C**** Written by Duane Malcolm, 26 August 2002

      IMPLICIT NONE


      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER MAXOPT,COUP_TYPE
!      REAL*8
      CHARACTER FSTRING*1024,ERROR*(*)
!     Local Variables
!      INTEGER
!      REAL*8

      CALL ENTERS('COUP_LIST',*9999)

      IF(COUP_TYPE.EQ.0)THEN
        MAXOPT=7 ! number of COUP_TYPE
      ELSEIF(COUP_TYPE.EQ.1)THEN
        FSTRING='(/'' Specify the CellML model #s [EXIT]:'''//
     '    '/$,''    '',I3,14(1X,I3))'
        MAXOPT=0
      ELSEIF(COUP_TYPE.EQ.2)THEN
        FSTRING='(/'' Specify the membrane channel #s [EXIT]:'''//
     '    '/''   (1) INa, sodium channels                     '''//
     '    '/''   (2) IK, potassium channels                   '''//
     '    '/''   (3) ICl, chloride channels                   '''//
     '    '/$,''    '',I3,14(1X,I3))'
        MAXOPT=3
      ELSEIF(COUP_TYPE.EQ.3)THEN
        FSTRING='(/'' Specify the membrane pump #s [EXIT]:'''//
     '    '/''   (1) IKNa, sodium/potassium ATPase pump    '''//
     '    '/$,''    '',I3,14(1X,I3))'
        MAXOPT=1
      ELSEIF(COUP_TYPE.EQ.4)THEN
        FSTRING='(/'' Specify the membrane exchangers #s [EXIT]:'''//
     '    '/$,''    '',I3,14(1X,I3))'
        MAXOPT=0
      ELSEIF(COUP_TYPE.EQ.5)THEN
        FSTRING='(/'' Specify the reaction pathway #s [EXIT]:'''//
     '    '/$,''    '',I3,14(1X,I3))'
        MAXOPT=0
      ELSEIF(COUP_TYPE.EQ.6)THEN
        FSTRING='(/'' Specify the physical process #s [EXIT]:'''//
     '    '/''   (1) Osmotic Pressure    '''//
     '    '/$,''    '',I3,14(1X,I3))'
        MAXOPT=1
      ELSEIF(COUP_TYPE.EQ.7)THEN
        FSTRING='(/'' Specify the polynomial function #s [EXIT]:'''//
     '    '/''   (1) linear y=ax+b           '''//
     '    '/''   (2) quadratic y=ax^2+bx+c   '''//
     '    '/''   (3) cubic y=ax^3+bx^2+cx+d  '''//
     '    '/$,''    '',I3,14(1X,I3))'
        MAXOPT=3
      ENDIF

      CALL EXITS('COUP_LIST')
      RETURN
 9999 CALL ERRORS('COUP_LIST',ERROR)
      CALL EXITS('COUP_LIST')
      RETURN 1
      END


