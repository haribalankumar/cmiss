C*************FE00
C Subroutine WRITES  write string
C*************FE01
C Function    ABBREV    true if 1st argument is abbreviation of 2nd
C Function    ABBRV     true if 1st argument is abbreviation of 2nd
C Subroutine ALLOCATE_MEMORYFrees (if necessary) and allocates mem
C Subroutine ASSERT  return error if logical expression false
C Function    CFROMI    character string from integer variable
C Function    CLOCAT    returns position of string within string
C Subroutine CROSS   return the vector cross product of two vectors
C Subroutine  CUPPER    converts string to upper case
C Function    DET       determinant of 3*3 matrix
C Function    DXZ       gradient of curvilinear coords wrt cartesian
C Function    DZX       gradient of cartesian coords wrt curvilinear
C Function    DZXX      2nd deriv of cartesian coords wrt curvilinear
C Subroutine ERRORS  appends new error message to existing string
C Subroutine EXITS   traces exit from a subprogram
C Subroutine FREE_MEMORYfrees dynamically allocated memory.
C Subroutine FREETIMERfrees a timer handle
C Subroutine GET_MEMORYdynamically allocates memory
C Subroutine GETTIMERreturns a handle for use with timer
C Function    IFROMC    integer from character
C Subroutine INVERT  invert a 2*2 or 3*3 matrix
C Function    LOWCAS    true if argument lowercase
C Subroutine NORMALISEnormalises a vector
C Subroutine ROT_COORDSYSreturns orthog rotation of coord system
C Subroutine ROTATIONreturns orthogonal rotation of input vector
C Function    TIMER     calls a 'c' routine to get system times
C Subroutine TRAN    transform a 2nd order 3*3 tensor
C Subroutine TRIM    find 1st & last non-blank char.s in a string
C Function    VTIME     returns the current time in seconds
C*************FE02
C Function    GETNYR    returns the corresponding rhs ny for a lhs ny
C Subroutine CPCG    transfer material params to Gauss point array
C Subroutine CPXI    interpolate material params at Xi position
C Subroutine DEFMGRADRCeval defm gradient tensor wrt rc coords
C Subroutine DLZJDX  eval covar derivs of deformed theta coords
C Subroutine DXIDNU  calc undef derivs of Xi wrt Nu coords & inv.
C Subroutine DXIDZN  calc def derivs of Xi wrt nu coords & inverse
C Subroutine FIBRE_REF_VECScalc dir.n of fibre ref'nce axes at ng
C Subroutine FIBRE_REF_VECS_DEFcalc def fibre ref'nce axes at ng
C Subroutine MAT_VEC_DEFcalc dir.n cosines of deformed mat axes
C Subroutine MAT_VEC_NGcalc dir.n cosines of undef mat axes at ng
C Subroutine MAT_VEC_ROTATErotates fib ref vecs into matl vectors
C Subroutine MAT_VEC_XIcalc dir.n cosines of undef mat axes at xi
C Subroutine MELGE   calc #element variables for solution
C Subroutine TOFFEL  eval cpts of Christoffel symbol of 2nd kind
C Subroutine XEXG    Gauss point array from elem node array
C Subroutine XEXW    interp independ. var. array XE at Xi position
C Subroutine XGMG    metric tensor arrays from Gauss pt array
C Subroutine ZEZG    eval Gauss point array from element array
C Subroutine ZEZW    interp depend. var. array ZE at Xi position
C Subroutine ZGMG    eval cpts of metric tensor in deformed state
C*************FE05
C Function    PAF       eval 1D Fourier basis functions
C Function    PAL0      eval 1D Legendre aux basis fn  (constant)
C Function    PAL1      eval 1D Legendre aux basis fn  (linear)
C Function    PAL2      eval 1D Legendre aux basis fns (order>=2)
C Function    PAP4      eval 1D Pressure aux quartic (hat) bas fn
C Function    PF1       eval Fourier series basis function
C Function    PFXI      interpolate nodal array XE at XI
C Function    PGX       eval 1st deriv of basis fn wrt Xj
C Function    PH2       eval 1D quadratic Hermite basis function
C Function    PH3       eval 1D cubic Hermite basis function
C Function    PL1       eval 1D linear Lagrange basis function
C Function    PL2       eval 1D quadratic Lagrange basis function
C Function    PL3       eval 1D cubic Lagrange basis funtion
C Function    PSI1      eval tensor prod. Lagrange & Hermite basis fns
C Function    PSI8      eval auxiliary lagrange basis functions
C*************FE07
C Subroutine ZEES     calculate element tangent stiffness matrix
C*************FE10_SGI
C Subroutine ENTERS       traces entry to a subprogram
C Subroutine GET_COMMAND_LINE gets the command line arguments
C*************FE50
C Function    AG         integrand of virtual work equation
C Function    D_AG       deriv of integrand in virtual work equation
C Subroutine D_ENERGY calc 2nd derivs of strain energy function
C Subroutine D_PFRE_NEderiv of pres constr resid wrt mat/geom pars
C Subroutine D_PFRF   deriv of press loading resids
C Subroutine D_ZERE50 deriv of elem resids wrt mat/geom params
C Subroutine D_ZGTG53 deriv of stress cpts (3D) wrt mat/geom vars
C Subroutine D_ZGTG54 deriv of stress cpts (mem) wrt mat/geom vars
C Subroutine ENERGY   calc derivatives of strain energy function
C Subroutine PFRE_NE  calc stress constr resid terms for elem vars
C Subroutine PFRE_NP  calc stress constr resid terms for node vars
C Subroutine PFRF     calc pressure loading residual terms
C Subroutine USER51   user strain energy fn of principal invars
C Subroutine USER52   user strain energy fn of principal extns
C Subroutine USER53   user strain energy fn of fibre/trans strains
C Subroutine ZERE50   calc element residual
C Subroutine ZERE55   press-vol eles for ventric cavity loading
C Subroutine ZETX50   calc stress cpts wrt fibre/ref coords
C Subroutine ZGTG51   eval stress tensor for pl stress problems
C Subroutine ZGTG52   eval stress tensor for pl strain problems
C Subroutine ZGTG53   eval stress tensor for 3D problems
C Subroutine ZGTG54   eval cpts of stress tensor for membrane
C Subroutine ZGTG55   eval cpts of stress tensor for string
C Subroutine ZGTG5A   eval actively developed fibre stress
C*************FZ30
C Subroutine ZERE30
C*************FZ40
C Subroutine ZERE40


      SUBROUTINE WRITES(id,OP_STRING,ERROR,*)

C#### Subroutine: WRITES
C###  Description:
C###    WRITES writes string to output identified by id as follows:
C**** id=2 (=IOOP) listing
C****    3 (=IODI) diagnostics
C****    4 (=IOTR) trace
C****    5 (=IOER) errors
C****    6 (=IOH1) ? help
C****    7 (=IOH2) ?? help
C****    8 (=IOH3) ??? help
C****    9
C**** id>9 is used for file output
C**** Note: WRITES cannot call Enters or Exits (since would get into
C**** an infinite loop).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:disp00.cmn'
      INCLUDE 'cmiss$reference:fsklib.inc'
!     Parameter List
      INTEGER id
      CHARACTER ERROR*(*),OP_STRING(*)*(*)
!     Local Variables
      INTEGER i,IBEG,IEND(100),j,NUM_BLANKS,NUM_RECORDS
      INTEGER DATA_TYPE
      PARAMETER(DATA_TYPE = 0)  !Cmiss Output
      INTEGER CLEN,INTSTR(1024)

C *** Calculate number of records in OP_STRING
      IF(OS_TYPE(1:3).EQ.'VMS') THEN
        i=1
        DO WHILE (OP_STRING(i)(1:1).NE.CHAR(0).AND.I.LT.100)
          i=i+1
        ENDDO
        NUM_RECORDS=i-1
        IF(i.EQ.100) NUM_RECORDS=100
      ELSE IF(OS_TYPE(1:4).EQ.'UNIX') THEN
        i=1
        NUM_BLANKS=0
        CALL STRING_TRIM(OP_STRING(i),IBEG,IEND(i))
        DO WHILE (i.LT.100.AND.NUM_BLANKS.LT.2)
          i=i+1
          CALL STRING_TRIM(OP_STRING(i),IBEG,IEND(i))
          IF(IEND(i).EQ.1) THEN
            NUM_BLANKS=NUM_BLANKS+1
          ELSE
            NUM_BLANKS=0
          ENDIF
        ENDDO
        NUM_RECORDS=i-NUM_BLANKS
        IF(i.EQ.100) NUM_RECORDS=100
      ENDIF

      DO i=1,NUM_RECORDS
        CALL STRING_TRIM(OP_STRING(i),IBEG,IEND(i))
        IEND(i)=IEND(i)+1
      ENDDO

      IF(id.LE.9) THEN !not file output
        IF(USE_SOCKET) THEN !transfer string to frontend thru socket
          DO i=1,NUM_RECORDS
            IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1)
     '        GOTO 9999
            IF(FSKWRITE(id,SK_LONG_INT,1,CONNID2).EQ.-1)
     '        GOTO 9999
            CLEN=FSKLEN(OP_STRING(i))
            CALL FSKF2C(OP_STRING(i),CLEN,INTSTR)
            IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1)
     '        GOTO 9999
            IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1)
     '        GOTO 9999
          ENDDO
        ELSE IF(.NOT.USE_SOCKET) THEN !write out directly
          DO i=1,NUM_RECORDS
C CPB 8/2/93 increasing string length before cut-off to 132
            IF(IEND(i).LE.132) THEN
              WRITE(*,'(A)') OP_STRING(i)(1:IEND(i))
            ELSE IF(IEND(i).GT.132.AND.IEND(i).LE.255) THEN
              WRITE(*,'(A)') OP_STRING(i)(1:132)
              WRITE(*,'(A)') OP_STRING(i)(133:IEND(i))
            ELSE
              WRITE(*,'(A)') OP_STRING(i)(1:132)
              WRITE(*,'(A)') OP_STRING(i)(133:255)
            ENDIF
          ENDDO
        ENDIF

      ELSE !file output
        DO i=1,NUM_RECORDS
          WRITE(id,'(A)') OP_STRING(i)(1:IEND(i))
        ENDDO

      ENDIF

C*** Output strings to echo output file if required

      IF(ECHO_OUTPUT) THEN
        DO i=1,NUM_RECORDS
          IF(IEND(i).LE.132) THEN
            WRITE(IOOUT,'(A)') OP_STRING(i)(1:IEND(i))
          ELSE IF(IEND(i).GT.132.AND.IEND(i).LE.255) THEN
            WRITE(IOOUT,'(A)') OP_STRING(i)(1:132)
            WRITE(IOOUT,'(A)') OP_STRING(i)(133:IEND(i))
          ELSE
            WRITE(IOOUT,'(A)') OP_STRING(i)(1:132)
            WRITE(IOOUT,'(A)') OP_STRING(i)(133:255)
          ENDIF
        ENDDO
      ENDIF

C *** Reset OP_STRING
      IF(OS_TYPE(1:3).EQ.'VMS') THEN
        DO i=1,NUM_RECORDS
          OP_STRING(i)(1:1)=CHAR(0)
        ENDDO
      ELSE IF(OS_TYPE(1:4).EQ.'UNIX') THEN
        DO i=1,NUM_RECORDS
          DO j=1,MXCH
            OP_STRING(i)(j:j)=' '
          ENDDO
        ENDDO
      ENDIF

      RETURN
 9999 CALL ERRORS('WRITES',ERROR)
      RETURN 1
      END


      LOGICAL FUNCTION ABBREV(SHORT,LONG,MNCH)

C#### Function: ABBREV
C###  Type: LOGICAL
C###  Description:
C###    ABBREV returns .TRUE. if the character string SHORT is an
C###    abbreviation of LONG.  SHORT must be at least MNCH characters
C###    long.  As a side effect, if ABBREV is .TRUE., SHORT is
C###    overwritten with LONG.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER MNCH
      CHARACTER LONG*(*),SHORT*(*)
!     Local Variables
      CHARACTER BUFFER*(MXCH),ERROR*(500)
      LOGICAL ABBRV

C      BUFFER=CUPPER(SHORT)
      CALL CUPPER(SHORT,BUFFER,ERROR,*9999)
      ABBREV=ABBRV(BUFFER,LONG,MNCH)
      IF(ABBREV) THEN
        SHORT=LONG
      ENDIF
      RETURN
 9999 WRITE(OP_STRING,'('' >>ERROR: RETURN FROM CUPPER'',A)') ERROR
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
 9998 RETURN
      END


      LOGICAL FUNCTION ABBRV(SHORT,LONG,MNCH)

C#### Function: ABBRV
C###  Type: LOGICAL
C###  Description:
C###    ABBRV returns .TRUE. if the character string SHORT is an
C###    abbreviation of LONG.  SHORT must be at least MNCH
C###    characters long.

      IMPLICIT NONE
!     Parameter List
      INTEGER MNCH
      CHARACTER LONG*(*),SHORT*(*)
!     Local Variables
      INTEGER noch,NXCH

      ABBRV=.FALSE.
      NXCH=MIN(LEN(LONG),LEN(SHORT))
      DO noch=MNCH,NXCH
        IF(SHORT.EQ.LONG(:noch)) THEN
          ABBRV=.TRUE.
          RETURN
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE ALLOCATE_MEMORY(N,ITEM_TYPE,ADDR_PTR,INITIALISE,
     '  ERROR,*)

C#### Subroutine: ALLOCATE_MEMORY
C###  Description:
C###    ALLOCATE_MEMORY dynamically allocates N items of type ITEM_TYPE
C###    in memory. If the address pointer is not zero ALLOCATE_MEMORY
C###    will free the memory first before allocating. ITEM_TYPES are
C###    1/2/3/4/5/6/7/8 for INTEGER,REAL,REAL*8,CHARACTER,LOGICAL,
C###    INTEGER*2,COMPLEX and COMPLEX*16. Constants are setup in
C###    mach00.cmn for these types. If INITIALISE is true the variables
C###    are initialised to a special value (currently only for
C###    INTEGER's, REAL's and REAL*8's).
C###  See-Also: GET_MEMORY,FREE_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER N,ITEM_TYPE
      INTEGER*4 ADDR_PTR
      LOGICAL INITIALISE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('ALLOCATE_MEMORY',*9999)


      IF(ADDR_PTR.NE.0) CALL FREE_MEMORY(ADDR_PTR,ERROR,*9999)

      CALL GET_MEMORY(N,ITEM_TYPE,ADDR_PTR,INITIALISE,ERROR,*9999)

      CALL EXITS('ALLOCATE_MEMORY')
      RETURN
 9999 CALL ERRORS('ALLOCATE_MEMORY',ERROR)
      CALL EXITS('ALLOCATE_MEMORY')
      RETURN 1
      END


      SUBROUTINE ASSERT(LOEXPR,ERRMSG,ERROR,*)

C#### Subroutine: ASSERT
C###  Description:
C###    ASSERT assigns the error message ERRMSG to ERROR if the
C###    logical expression LOEXPR is false, and control is returned
C###    via the alternative return argument.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERRMSG*(*),ERROR*(*)
      LOGICAL LOEXPR

      CALL ENTERS('ASSERT',*9999)
      IF(LOEXPR) THEN
        ERROR=' '
      ELSE
        ERROR=ERRMSG
        GO TO 9999
      ENDIF

      CALL EXITS('ASSERT')
      RETURN
 9999 CALL EXITS('ASSERT')
      RETURN 1
      END


      CHARACTER*(*) FUNCTION CFROMI(IDATA,FORMAT)

C#### Function: CFROMI
C###  Type: CHARACTER
C###  Description:
C###    Converts integer variable IDATA to character string CFROMI
C###    as specified by the character string FORMAT.

      IMPLICIT NONE
!     Parameter List
      INTEGER IDATA
      CHARACTER FORMAT*(*)
!     Local Variables
      CHARACTER CDATA*500

      WRITE(UNIT=CDATA,FMT=FORMAT) IDATA
      CFROMI=CDATA

      RETURN
      END


      INTEGER FUNCTION CLOCAT(SUBSTR,STRING)

C#### Function: CLOCAT
C###  Type: INTEGER
C###  Description:
C###    CLOCAT returns the position of the first occurrence of SUBSTR
C###    within STRING.  If SUBSTR is not found within STRING, CLOCAT
C###    returns zero.

      IMPLICIT NONE
!     Parameter List
      CHARACTER SUBSTR*(*),STRING*(*)
!     Local Variables
      INTEGER l0,l1,LEN,LENSUB,LENSTR

      LENSUB=LEN(SUBSTR)
      LENSTR=LEN(STRING)
      l0=LENSUB-1
      CLOCAT=0
      DO l1=1,(LENSTR-LENSUB)
        IF(SUBSTR.EQ.STRING(l1:l0+l1)) THEN
          CLOCAT=l1
          RETURN
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE CROSS(A,B,C)

C#### Subroutine: CROSS
C###  Description:
C###    CROSS returns the vector cross product of A*B in C.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A(3),B(3),C(3)

C     CALL ENTERS('CROSS',*9999)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)

C     CALL EXITS('CROSS')
      RETURN
      END


C old MLB 19-9-96 changed to subroutine
C      CHARACTER*(*) FUNCTION CUPPER(STRING)
C
CC#### Function: CUPPER
CC###  Type: CHARACTER
CC###  Description:
CC###    CUPPER converts the character string STRING to upper case.
C
C      IMPLICIT NONE
C!     Parameter List
C      CHARACTER STRING*(*)
C!     Local Variables
C      INTEGER i,OFFSET
C      LOGICAL LOWCAS
C
C      OFFSET=ICHAR('A')-ICHAR('a')
C      CUPPER=STRING
C      DO i=1,LEN(STRING)
C        IF(LOWCAS(STRING(i:i))) THEN
C          CUPPER(i:i)=CHAR(ICHAR(STRING(i:i))+OFFSET)
C        ENDIF
C      ENDDO
C      RETURN
C      END


      SUBROUTINE CUPPER(INSTRING,OUTSTRING,ERROR,*)

C#### Subroutine: CUPPER
C###  Description:
C###    CUPPER converts a string to upper case.
C###    It has two arguments (+ error args), the first is the
C###    input string and the second is the output string.

      IMPLICIT NONE
!     Parameter List
      CHARACTER INSTRING*(*),OUTSTRING*(*),ERROR*(*)
!     Local Variables
      INTEGER i,OFFSET
      LOGICAL LOWCAS

      CALL ENTERS('CUPPER',*9999)

      OFFSET=ICHAR('A')-ICHAR('a')
      CALL ASSERT(LEN(INSTRING).LE.LEN(OUTSTRING),'>>Input string'
     '  //'longer than ouput string in CUPPER',ERROR,*9999)
      OUTSTRING=INSTRING
      DO i=1,LEN(INSTRING)
        IF(LOWCAS(INSTRING(i:i))) THEN
          OUTSTRING(i:i)=CHAR(ICHAR(INSTRING(i:i))+OFFSET)
        ENDIF
      ENDDO
      CALL EXITS('CUPPER')
      RETURN
 9999 CALL ERRORS('CUPPER',ERROR)
      CALL EXITS('CUPPER')
      RETURN 1
      END


      REAL*8 FUNCTION DET(A)

C#### Function: DET
C###  Type: REAL*8
C###  Description:
C###    DET evaluates determinant of 3*3 matrix A.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A(3,3)

      DET=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
     '   +A(1,2)*(A(2,3)*A(3,1)-A(3,3)*A(2,1))
     '   +A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

      RETURN
      END


      REAL*8 FUNCTION DXZ(ICOORD,i,k,X)

C#### Function: DXZ
C###  Type: REAL*8
C###  Description:
C###    DXZ calculates DX(I)/DZ(K) at X, where Z(K) are rectangular
C###    Cartesian and X(I) are curvilinear coordinates defined by
C###    ICOORD.
C###    ICOORD=1 for rectangular cartesian coordinates;
C###           2 for cylindrical polar coordinates;
C###           3 for spherical polar coordinates;
C###           4 for prolate spheriodal coordinates; and
C###           5 for oblate spheroidal coordinates.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b14.cmn'
!     Parameter List
      INTEGER ICOORD,i,k
      REAL*8 RD,X(3)

      GO TO (100,200,300,400,500),ICOORD
        DXZ=0.0D0
        RETURN
  100   DXZ=0.0D0
        IF(i.EQ.k) DXZ=1.0D0
        RETURN
  200   GO TO (210,220,230),i
  210     GO TO (211,212,213),k
  211       DXZ=DCOS(X(2))
            RETURN
  212       DXZ=DSIN(X(2))
            RETURN
  213       DXZ=0.0D0
            RETURN
  220     GO TO (221,222,223),k
  221       DXZ=-DSIN(X(2))/X(1)
            RETURN
  222       DXZ=DCOS(X(2))/X(1)
            RETURN
  223       DXZ=0.0D0
            RETURN
  230     GO TO (231,232,233),k
  231       DXZ=0.0D0
            RETURN
  232       DXZ=0.0D0
            RETURN
  233       DXZ=1.0D0
            RETURN
  300   GO TO (310,320,330),i
  310     GO TO (311,312,313),k
  311       DXZ=DCOS(X(2))*DCOS(X(3))
            RETURN
  312       DXZ=DSIN(X(2))*DCOS(X(3))
            RETURN
  313       DXZ=DSIN(X(3))
            RETURN
  320     GO TO (321,322,323),k
  321       DXZ=-DSIN(X(2))/(X(1)*DCOS(X(3)))
            RETURN
  322       DXZ=DCOS(X(2))/(X(1)*DCOS(X(3)))
            RETURN
  323       DXZ=0.0D0
            RETURN
  330     GO TO (331,332,333),k
  331       DXZ=-DCOS(X(2))*DSIN(X(3))/X(1)
            RETURN
  332       DXZ=-DSIN(X(2))*DSIN(X(3))/X(1)
            RETURN
  333       DXZ=DCOS(X(3))/X(1)
            RETURN
  400   RD=FOCUS*(DCOSH(X(1))*DCOSH(X(1))-DCOS(X(2))*DCOS(X(2)))
        GO TO (410,420,430),i
  410     GO TO (411,412,413),k
  411       DXZ=DSINH(X(1))*DCOS(X(2))/RD
            RETURN
  412       DXZ=DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))/RD
            RETURN
  413       DXZ=DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))/RD
            RETURN
  420     GO TO (421,422,423),k
  421       DXZ=-DCOSH(X(1))*DSIN(X(2))/RD
            RETURN
  422       DXZ=DSINH(X(1))*DCOS(X(2))*DCOS(X(3))/RD
            RETURN
  423       DXZ=DSINH(X(1))*DCOS(X(2))*DSIN(X(3))/RD
            RETURN
  430     GO TO (431,432,433),k
  431       DXZ=0.0D0
            RETURN
  432       DXZ=-DSIN(X(3))/(FOCUS*DSINH(X(1))*DSIN(X(2)))
            RETURN
  433       DXZ=DCOS(X(3))/(FOCUS*DSINH(X(1))*DSIN(X(2)))
            RETURN
  500   GO TO (510,520,530),i
  510     GO TO (511,512,513),k
  511       DXZ=0.0D0
            RETURN
  512       DXZ=0.0D0
            RETURN
  513       DXZ=0.0D0
            RETURN
  520     GO TO (521,522,523),k
  521       DXZ=0.0D0
            RETURN
  522       DXZ=0.0D0
            RETURN
  523       DXZ=0.0D0
            RETURN
  530     GO TO (531,532,533),k
  531       DXZ=0.0D0
            RETURN
  532       DXZ=0.0D0
            RETURN
  533       DXZ=0.0D0
            RETURN
      END


      REAL*8 FUNCTION DZX(ICOORD,i,j,X)

C#### Function: DZX
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    DZX calculates DZ(i)/DX(j) at X, where Z(i) are rectangular
C###    Cartesian and X(j) are curvilinear coordinates defined by
C###    ICOORD.
C###    <PRE>
C###    ICOORD=1 for rectangular cartesian coordinates;
C###           2 for cylindrical polar coordinates;
C###           3 for spherical polar coordinates;
C###           4 for prolate spheriodal coordinates; and
C###           5 for oblate spheroidal coordinates.
C###    </PRE></HTML>

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b14.cmn'
!     Parameter List
      INTEGER ICOORD,i,j
      REAL*8  X(3)

      GO TO (100,200,300,400,500),ICOORD
        DZX=0.0D0
        RETURN
  100   DZX=0.0D0
        IF(i.EQ.j) DZX=1.0D0
        RETURN
  200   GO TO (210,220,230),i
  210     GO TO (211,212,213),j
  211       DZX=DCOS(X(2))
            RETURN
  212       DZX=-X(1)*DSIN(X(2))
            RETURN
  213       DZX=0.0D0
            RETURN
  220     GO TO (221,222,223),j
  221       DZX=DSIN(X(2))
            RETURN
  222       DZX=X(1)*DCOS(X(2))
            RETURN
  223       DZX=0.0D0
            RETURN
  230     GO TO (231,232,233),j
  231       DZX=0.0D0
            RETURN
  232       DZX=0.0D0
            RETURN
  233       DZX=1.0D0
            RETURN
  300   GO TO (310,320,330),i
  310     GO TO (311,312,313),j
  311       DZX=DCOS(X(2))*DCOS(X(3))
            RETURN
  312       DZX=-X(1)*DSIN(X(2))*DCOS(X(3))
            RETURN
  313       DZX=-X(1)*DCOS(X(2))*DSIN(X(3))
            RETURN
  320     GO TO (321,322,323),j
  321       DZX=DSIN(X(2))*DCOS(X(3))
            RETURN
  322       DZX=X(1)*DCOS(X(2))*DCOS(X(3))
            RETURN
  323       DZX=-X(1)*DSIN(X(2))*DSIN(X(3))
            RETURN
  330     GO TO (331,332,333),j
  331       DZX=DSIN(X(3))
            RETURN
  332       DZX=0.0D0
            RETURN
  333       DZX=X(1)*DCOS(X(3))
            RETURN
  400   GO TO (410,420,430),i
  410     GO TO (411,412,413),j
  411       DZX=FOCUS*DSINH(X(1))*DCOS(X(2))
            RETURN
  412       DZX=-FOCUS*DCOSH(X(1))*DSIN(X(2))
            RETURN
  413       DZX=0.0D0
            RETURN
  420     GO TO (421,422,423),j
  421       DZX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            RETURN
  422       DZX=FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            RETURN
  423       DZX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            RETURN
  430     GO TO (431,432,433),j
  431       DZX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            RETURN
  432       DZX=FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            RETURN
  433       DZX=FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            RETURN
  500   GO TO (510,520,530),i
  510     GO TO (511,512,513),j
  511       DZX= FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            RETURN
  512       DZX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            RETURN
  513       DZX=-FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))
            RETURN
  520     GO TO (521,522,523),j
  521       DZX= FOCUS*DCOSH(X(1))*DSIN(X(2))
            RETURN
  522       DZX= FOCUS*DSINH(X(1))*DCOS(X(2))
            RETURN
  523       DZX= 0.0D0
            RETURN
  530     GO TO (531,532,533),j
  531       DZX= FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            RETURN
  532       DZX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            RETURN
  533       DZX= FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
            RETURN
      END


      REAL*8 FUNCTION DZXX(ICOORD,i,j,k,X)

C#### Function: DZXX
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    DZXX calculates D2Z(i)/[DX(j).DX(k)] at X, where Z(i) are
C###    rectangular Cartesian and X(k) are curvilinear coordinates
C###    defined by ICOORD.
C###    <PRE>
C###    ICOORD=1 for rectangular cartesian coordinates;
C###           2 for cylindrical polar coordinates;
C###           3 for spherical polar coordinates;
C###           4 for prolate spheriodal coordinates; and
C###           5 for oblate spheroidal coordinates.
C###    </PRE></HTML>

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER ICOORD,i,j,k
      REAL*8  X(3)
      CHARACTER ERROR*10

      DZXX=0.0d0
      IF(ICOORD.EQ.1) THEN !rectangular cartesian coords
C       second derivs all zero

      ELSE IF(ICOORD.EQ.2) THEN !cylindrical polar coords
        IF(k.EQ.1) THEN
          IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.2) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-X(1)*DCOS(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DSIN(X(2))
            ENDIF
          ENDIF
        ENDIF

      ELSE IF(ICOORD.EQ.3) THEN !spherical polar coords
        IF(k.EQ.1) THEN
          IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.1) THEN
              DZXX=-DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=DCOS(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.2) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-X(1)*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DSIN(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.1) THEN
              DZXX=X(1)*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.3) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=X(1)*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.1) THEN
              DZXX=-X(1)*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=-X(1)*DSIN(X(3))
            ENDIF
          ENDIF
        ENDIF

      ELSE IF(ICOORD.EQ.4) THEN !prolate spheroidal coords
        IF(k.EQ.1) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.2) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-FOCUS*DCOSH(X(1))*DCOS(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.3) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            ENDIF
          ENDIF
        ENDIF

      ELSE IF(ICOORD.EQ.5) THEN !oblate spheroidal coords
        WRITE(OP_STRING,'('' 2nd derivs not implemented for oblate '
     '    //'spheroidal coords'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


      SUBROUTINE ERRORS(NAME,ERROR)

C#### Subroutine: ERRORS
C###  Description:
C###    ERRORS is the standard error handling routine.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*),NAME*(*)
!     Local Variables
      INTEGER IBEG,IEND

C     CALL ENTERS('ERRORS',*9999)
      CALL STRING_TRIM(ERROR,IBEG,IEND)
      ERROR=ERROR(IBEG:IEND)//'>'//NAME

C     CALL EXITS('ERRORS')
      RETURN
      END


      SUBROUTINE EXITS(NAME)

C#### Subroutine: EXITS
C###  Description:
C###    EXITS traces exit from a subprogram recording the level of
C###    nesting, the exit time, and writing the subprogram name to a
C###    trace file.  If diagnostics are on (DIAGNO=.TRUE.) DOP is
C###    switched off if .NOT.ALLSUB & NAME=SUBNAM.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      CHARACTER NAME*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,n1lv,n1sb
      REAL*8 TM,TMEX,TMTL,VTIME
      CHARACTER CFROMI*10,COLV*10,ERROR*10

      IF(DIAGNO) THEN
        CALL STRING_TRIM(NAME,IBEG1,IEND1)
        IF(ALLSUB) THEN
          WRITE(OP_STRING,'('' *** Exit  '',A)') NAME(IBEG1:IEND1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE IF(.NOT.ALLSUB) THEN
          IF(DOP_STACK(NUM_STACK)) THEN
            WRITE(OP_STRING,'('' *** Exit  '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          NUM_STACK=NUM_STACK-1
          DOP=DOP_STACK(NUM_STACK)
        ENDIF
      ENDIF

      IF(TR01) THEN
        IF(NAME.EQ.TRSB) THEN
          TR02=.FALSE.
        ENDIF
        TR01=.FALSE.
        IF((nolv.GT.0).AND.(nolv.LE.NXLV)) THEN
          TM=VTIME()
          TMEX=TM-TMST
          TMTL=TMEX-TMEN(nolv)
          TMEL(nolv)=TMTL-TMEL(nolv)
          nosb=NOSBLV(nolv)
          IF(nosb.LE.NXSB) THEN
            IF(SB(nosb).NE.NAME) THEN
              WRITE(OP_STRING,'(T32,''!Error in EXIT: '''//
     '          ',A,'' entered, '',A,'' exited.'')')
     '          SB(nosb),NAME
              CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
              RETURN
            ENDIF
          ENDIF
          TMELSM(nosb)=TMELSM(nosb)+TMEL(nolv)
          TMTLSM(nosb)=TMTLSM(nosb)+TMTL
          IF(TR02) THEN
            COLV=CFROMI(nolv,'(I10)')
            CALL STRING_TRIM(COLV,IBEG,IEND)
            WRITE(OP_STRING,'(1X,3(F10.3),'//COLV(IBEG:IEND)//
     '        '('' <''),A)')
     '        TMEX,TMTL,TMEL(nolv),NAME
            CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
          ENDIF
          DO n1lv=1,nolv-1
            TMEL(n1lv)=TMEL(n1lv)+TMEL(nolv)
          ENDDO
        ENDIF
        nolv=nolv-1
        IF(nolv.EQ.0) THEN
          WRITE(OP_STRING,'(/''                 Average times:''/'//
     '      '''     Calls:    Total:   Actual:  Subprogram:'')')
          CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
          DO n1sb=1,NTSB
            WRITE(OP_STRING,'(1X,I10,2(F10.3),2X,A)')
     '        NOSBSM(n1sb),TMTLSM(n1sb)/NOSBSM(n1sb),
     '        TMELSM(n1sb)/NOSBSM(n1sb),SB(n1sb)
            CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF
        IF(NAME.EQ.TRSB) THEN
          TR02=.FALSE.
        ENDIF
        TR01=.TRUE.
      ENDIF

 9999 RETURN
      END


      SUBROUTINE FREE_MEMORY(ADDR_PTR,ERROR,*)

C#### Subroutine: FREE_MEMORY
C###  Description:
C###    FREE_MEMORY deallocates memory previous allocated with
C###    GET_MEMORY.
C###  See-Also: ALLOCATE_MEMORY,GET_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER*4 ADDR_PTR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('FREE_MEMORY',*9999)

      ERR=0
      CALL FREEMEMORY(ADDR_PTR,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        ERROR(CERRLEN:)=' '
        GOTO 9999
      ENDIF

      CALL EXITS('FREE_MEMORY')
      RETURN
 9999 CALL ERRORS('FREE_MEMORY',ERROR)
      CALL EXITS('FREE_MEMORY')
      RETURN 1
      END


      SUBROUTINE FREETIMER(IHANDLE,ERROR,*)

C#### Subroutine: FREETIMER
C###  Description:
C###    Frees a timer allocated via GETTIMER
C###  See-Also: GETTIMER

      IMPLICIT NONE
!     Parameter List
      INTEGER IHANDLE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('FREETIMER',*9999)

      ERR=0
      CALL CFREETIMER(IHANDLE,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        ERROR(CERRLEN:)=' '
        GOTO 9999
      ENDIF

      CALL EXITS('FREETIMER')
      RETURN
 9999 CALL ERRORS('FREETIMER',ERROR)
      CALL EXITS('FREETIMER')
      RETURN 1
      END


      SUBROUTINE GET_MEMORY(N,ITEM_TYPE,ADDR_PTR,INITIALISE,ERROR,*)

C#### Subroutine: GET_MEMORY
C###  Description:
C###    GET_MEMORY dynamically allocates N items of type ITEM_TYPE
C###    in memory. ITEM_TYPES are 1/2/3/4/5/6/7/8 for INTEGER,REAL,
C###    REAL*8,CHARACTER,LOGICAL,INTEGER*2,COMPLEX and COMPLEX*16.
C###    Constants are setup in mach00.cmn for these types. If
C###    INITIALISE is true the variables are initialised to a
C###    special value (currently only for INTEGER's, REAL's and
C###    REAL*8's).
C###  See-Also: ALLOCATE_MEMORY,FREE_MEMORY

      IMPLICIT NONE
!     Parameter List
      INTEGER N,ITEM_TYPE
      INTEGER*4 ADDR_PTR
      LOGICAL INITIALISE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR,INIT

      CALL ENTERS('GET_MEMORY',*9999)

      ERR=0
      IF(INITIALISE) THEN
        INIT=1
      ELSE
        INIT=0
      ENDIF
      CALL MALLOCMEMORY(N,ITEM_TYPE,INIT,ADDR_PTR,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        GOTO 9999
      ENDIF

      CALL EXITS('GET_MEMORY')
      RETURN
 9999 CALL ERRORS('GET_MEMORY',ERROR)
      CALL EXITS('GET_MEMORY')
      RETURN 1
      END


      SUBROUTINE GETTIMER(IHANDLE,ERROR,*)

C#### Subroutine: GETTIMER
C###  Description:
C###    Returns a new timer handle for use with timer

      IMPLICIT NONE
!     Parameter List
      INTEGER IHANDLE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('GETTIMER',*9999)

      ERR=0
      CALL CGETTIMER(IHANDLE,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        ERROR(CERRLEN:)=' '
        GOTO 9999
      ENDIF

      CALL EXITS('GETTIMER')
      RETURN
 9999 CALL ERRORS('GETTIMER',ERROR)
      CALL EXITS('GETTIMER')
      RETURN 1
      END


      INTEGER FUNCTION IFROMC(CDATA)

C#### Function: IFROMC
C###  Type: INTEGER
C###  Description:
C###    IFROMC converts character string CDATA to an integer variable.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      CHARACTER CDATA*(*)
!     Local Variables
      INTEGER IDATA
      CHARACTER ERROR*10

C     CALL ENTERS('IFROMC',*9999)
      READ(UNIT=CDATA,FMT=*,ERR=9998) IDATA
      IFROMC=IDATA

C     CALL EXITS('IFROMC')
      RETURN
 9998 IFROMC=0
      WRITE(OP_STRING,'('' IFROMC cannot convert '',A,'
     '  //''' to an integer'')') CDATA
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 9999 RETURN
      END


      SUBROUTINE INVERT(n,A,B,AA)

C#### Subroutine: INVERT
C###  Description:
C###    INVERT returns the inverse of matrix A as B and det(A) as AA.
C###    Matrix A may be no larger than 3*3 (N=3). Note that in
C###    both cases A and B are dimensioned to A(3,3) and B(3,3).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER n
      REAL*8 A(3,3),AA,B(3,3)
!     Local Variables
      CHARACTER ERROR*10
      INTEGER i,j,M(5)
      REAL*8 DET,ZERO_TOLERANCE

      PARAMETER(ZERO_TOLERANCE=1.0d-15)

      DATA M/1,2,3,1,2/

C     CALL ENTERS('INVERT',*9999)

      IF(N.EQ.1) THEN !1*1 matrix (for completeness - do not delete)
        AA=A(1,1)
        IF(DABS(AA).GT.ZERO_TOLERANCE) THEN
          B(1,1)=1.0d0/AA
        ELSE
          WRITE(OP_STRING,'('' >>Warning: A(1,1) is zero in INVERT'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(N.EQ.2) THEN !2*2 matrix
        AA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
        IF(DABS(AA).GT.ZERO_TOLERANCE) THEN
          B(1,1)= A(2,2)/AA
          B(2,1)=-A(2,1)/AA
          B(1,2)=-A(1,2)/AA
          B(2,2)= A(1,1)/AA
        ELSE
          WRITE(OP_STRING,
     '      '('' >>Warning: Zero determinant in 2*2 INVERT'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(N.EQ.3) THEN !3*3 matrix
        AA=DET(A)
        IF(DABS(AA).GT.ZERO_TOLERANCE) THEN
          DO i=1,3
            DO j=1,3
              B(i,j)=(A(M(j+1),M(i+1))*A(M(j+2),M(i+2))
     '               -A(M(j+2),M(i+1))*A(M(j+1),M(i+2)))/AA
            ENDDO !j
          ENDDO !i
        ELSE
          WRITE(OP_STRING,
     '      '('' >>Warning: Zero determinant in 3*3 INVERT'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE
        WRITE(OP_STRING,
     '    '('' >>Warning: Matrix larger than 3x3 - cannot INVERT!'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C     CALL EXITS('INVERT')

 9999 RETURN
      END


      LOGICAL FUNCTION LOWCAS(CHAR)

C#### Function: LOWCAS
C###  Type: LOGICAL
C###  Description:
C###    LOWCAS returns .TRUE. if the character CHAR is lower
C###    case alphabetic.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CHAR

      IF(LGE(CHAR,'a').AND.LLE(CHAR,'z')) THEN
        LOWCAS=.TRUE.
      ELSE
        LOWCAS=.FALSE.
      ENDIF
      RETURN
      END


      SUBROUTINE NORMALISE(NUMCMPTS,VECTOR,ERROR,*)

C#### Subroutine: NORMALISE
C###  Description:
C###    NORMALISE divides the components of VECTOR by it's length

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER NUMCMPTS
      REAL*8 VECTOR(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni
      REAL*8 SQR_VECTOR_LENGTH,ZERO_TOL

      DATA ZERO_TOL /1.0D-8/

      CALL ENTERS('NORMALISE',*9999)

      SQR_VECTOR_LENGTH=0.0d0
      DO ni=1,NUMCMPTS
        SQR_VECTOR_LENGTH=SQR_VECTOR_LENGTH+VECTOR(ni)*VECTOR(ni)
      ENDDO !ni
      IF(SQR_VECTOR_LENGTH.GT.ZERO_TOL) THEN
        DO ni=1,NUMCMPTS
          VECTOR(ni)=VECTOR(ni)/DSQRT(SQR_VECTOR_LENGTH)
        ENDDO !ni
      ELSE
        WRITE(OP_STRING,
     '    '('' >>WARNING: Cannot normalise a zero length vector'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('NORMALISE')
      RETURN
 9999 CALL ERRORS('NORMALISE',ERROR)
      CALL EXITS('NORMALISE')
      RETURN 1
      END


      SUBROUTINE ROT_COORDSYS(id,ANGLE,COORD_MATRIX,ERROR,*)

C#### Subroutine: ROT_COORDSYS
C###  Description:
C###    ROT_COORDSYS returns 3D orthogonal rotation of input
C###    COORD_MATRIX about ID-axis. NOTE this corresponds to a
C###    post-multiplication of COORD_MATRIX by the transpose of
C###    the orthogonal rotation matrix

      IMPLICIT NONE
!     Parameter List
      INTEGER id
      REAL*8 ANGLE,COORD_MATRIX(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER id1,id2
      REAL*8 a(3),b(3),COS_ANGLE,SIN_ANGLE

      CALL ENTERS('ROT_COORDSYS',*9999)
      IF(id.EQ.1) THEN !rotate about 1st axis
        id1=2
        id2=3
      ELSE IF(id.EQ.2) THEN !rotate about 2nd axis
        id1=1
        id2=3
      ELSE IF(id.EQ.3) THEN !rotate about 3rd axis
        id1=1
        id2=2
      ENDIF
      COS_ANGLE=DCOS(ANGLE)
      SIN_ANGLE=DSIN(ANGLE)
      a(1)= COS_ANGLE*COORD_MATRIX(1,id1)+SIN_ANGLE*COORD_MATRIX(1,id2)
      a(2)= COS_ANGLE*COORD_MATRIX(2,id1)+SIN_ANGLE*COORD_MATRIX(2,id2)
      a(3)= COS_ANGLE*COORD_MATRIX(3,id1)+SIN_ANGLE*COORD_MATRIX(3,id2)
      b(1)=-SIN_ANGLE*COORD_MATRIX(1,id1)+COS_ANGLE*COORD_MATRIX(1,id2)
      b(2)=-SIN_ANGLE*COORD_MATRIX(2,id1)+COS_ANGLE*COORD_MATRIX(2,id2)
      b(3)=-SIN_ANGLE*COORD_MATRIX(3,id1)+COS_ANGLE*COORD_MATRIX(3,id2)
      COORD_MATRIX(1,id1)=a(1)
      COORD_MATRIX(2,id1)=a(2)
      COORD_MATRIX(3,id1)=a(3)
      COORD_MATRIX(1,id2)=b(1)
      COORD_MATRIX(2,id2)=b(2)
      COORD_MATRIX(3,id2)=b(3)

      CALL EXITS('ROT_COORDSYS')
      RETURN
 9999 CALL ERRORS('ROT_COORDSYS',ERROR)
      CALL EXITS('ROT_COORDSYS')
      RETURN 1
      END


      SUBROUTINE ROTATION(id,ANGLE,VECTOR,ERROR,*)

C#### Subroutine: ROTATION
C###  Description:
C###    ROTATION returns 3D orthogonal rotation of input vector
C###    about ID-axis.

      IMPLICIT NONE
!     Parameter List
      INTEGER id
      REAL*8 ANGLE,VECTOR(3)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 COS_ANGLE,SIN_ANGLE,V(3)

      CALL ENTERS('ROTATION',*9999)
      COS_ANGLE=DCOS(ANGLE)
      SIN_ANGLE=DSIN(ANGLE)

      IF(id.EQ.1) THEN      !rotate about 1st axis
        V(1)=VECTOR(1)
        V(2)=COS_ANGLE*VECTOR(2)-SIN_ANGLE*VECTOR(3)
        V(3)=SIN_ANGLE*VECTOR(2)+COS_ANGLE*VECTOR(3)
      ELSE IF(id.EQ.2) THEN !rotate about 2nd axis
        V(1)=COS_ANGLE*VECTOR(1)-SIN_ANGLE*VECTOR(3)
        V(2)=VECTOR(2)
        V(3)=SIN_ANGLE*VECTOR(1)+COS_ANGLE*VECTOR(3)
      ELSE IF(id.EQ.3) THEN !rotate about 3rd axis
        V(1)=COS_ANGLE*VECTOR(1)-SIN_ANGLE*VECTOR(2)
        V(2)=SIN_ANGLE*VECTOR(1)+COS_ANGLE*VECTOR(2)
        V(3)=VECTOR(3)
      ENDIF
      VECTOR(1)=V(1)
      VECTOR(2)=V(2)
      VECTOR(3)=V(3)

      CALL EXITS('ROTATION')
      RETURN
 9999 CALL ERRORS('ROTATION',ERROR)
      CALL EXITS('ROTATION')
      RETURN 1
      END


      REAL*8 FUNCTION TIMER(IHANDLE,T_FLAG)

C**** If T_FLAG=T_INITIALISE initialises timer.
C**** Otherwise returns elapsed time (if T_FLAG=T_ELAPSED) or
C**** CPU time (if T_FLAG=T_CPU) in seconds since initialisation.
C**** NOTE: ETIME uses single precision real numbers.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER IHANDLE,T_FLAG
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR
      REAL*8 RETURN_TIME
      CHARACTER ERROR*255

      ERR=0
      CALL CTIMER(RETURN_TIME,IHANDLE,T_FLAG,ERR,CERROR)
      TIMER = RETURN_TIME
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,OP_STRING(1))
        OP_STRING(1)(CERRLEN:)=' '
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


      SUBROUTINE TRAN(ip,A,B,C)

C#### Subroutine: TRAN
C###  Description:
C###    TRAN transforms 2nd order 3*3 tensor A into B with coordinate
C###    transformation matrix C,
C###              covariantly (C(t)AC) if IP=1,
C###           or contravariantly (CAC(t)) if IP=2,
C###           or mixed variantly (C(-1)AC -similarity trans.) if IP=3.

      IMPLICIT NONE
!     Parameter List
      INTEGER ip
      REAL*8  A(3,3),B(3,3),C(3,3)
!     Local Variables
      INTEGER i,j,k
      REAL*8 CC,D(3,3)

C     CALL ENTERS('TRAN',*9999)
      IF(ip.EQ.1) THEN
        DO i=1,3
          DO j=1,3
            B(i,j)=(A(1,1)*C(1,i)+A(2,1)*C(2,i)+A(3,1)*C(3,i))*C(1,j)
     '            +(A(1,2)*C(1,i)+A(2,2)*C(2,i)+A(3,2)*C(3,i))*C(2,j)
     '            +(A(1,3)*C(1,i)+A(2,3)*C(2,i)+A(3,3)*C(3,i))*C(3,j)
          ENDDO !j
        ENDDO !i

      ELSE IF(IP.EQ.2) THEN
        DO i=1,3
          DO j=1,3
            B(i,j)=(A(1,1)*C(i,1)+A(2,1)*C(i,2)+A(3,1)*C(i,3))*C(j,1)
     '            +(A(1,2)*C(i,1)+A(2,2)*C(i,2)+A(3,2)*C(i,3))*C(j,2)
     '            +(A(1,3)*C(i,1)+A(2,3)*C(i,2)+A(3,3)*C(i,3))*C(j,3)
          ENDDO !j
        ENDDO !i

      ELSE IF(IP.EQ.3) THEN
        CALL INVERT(3,C,D,CC)
        DO i=1,3
          DO j=1,3
            B(i,j)=(A(1,1)*D(i,1)+A(2,1)*D(i,2)+A(3,1)*D(i,3))*C(1,j)
     '            +(A(1,2)*D(i,1)+A(2,2)*D(i,2)+A(3,2)*D(i,3))*C(2,j)
     '            +(A(1,3)*D(i,1)+A(2,3)*D(i,2)+A(3,3)*D(i,3))*C(3,j)
          ENDDO !j
        ENDDO !i
      ENDIF

C     CALL EXITS('TRAN')
      RETURN
      END


      SUBROUTINE TRIM(CDATA,IBEGIN,IEND)

C#### Subroutine: TRIM
C###  Description:
C###   TRIM finds the first and last non-blank characters in a
C###    string CDATA.
C**** On exit IBEGIN is the position of the first non-blank character
C**** and IEND is the position of the last non-(blank/tab/null)
C**** character.
C**** If CDATA consists of only blank characters IEND will be 1 and
C**** IBEGIN will be LEN(CDATA).

      IMPLICIT NONE
!     Parameter List
      INTEGER IBEGIN,IEND
      CHARACTER CDATA*(*)
!     Local Variables
      INTEGER l,LENGTH
      CHARACTER NUL*1,TAB*1

      LENGTH=LEN(CDATA)
      NUL=CHAR(0) !ASCII NULL
      TAB=CHAR(9) !ASCII HT (horizontal TAB)
      DO 1 l=1,LENGTH
C MPN 11Mar96 Also ignore tabs and null chars
        IF(CDATA(l:l).NE.' '.AND.CDATA(l:l).NE.TAB.AND.
     '    CDATA(l:l).NE.NUL) THEN
C old        IF(CDATA(l:l).NE.' ') THEN
          IBEGIN=l
          GOTO 2
        ENDIF
 1    CONTINUE
      IBEGIN=LENGTH
 2    DO 3 l=LENGTH,1,-1
C MPN 11Mar96 Also ignore tabs and null chars
        IF(CDATA(l:l).NE.' '.AND.CDATA(l:l).NE.TAB.AND.
     '    CDATA(l:l).NE.NUL) THEN
C old        IF(CDATA(l:l).NE.' ') THEN
          IEND=l
          GOTO 4
        ENDIF
 3    CONTINUE
      IEND=1
 4    IF(IBEGIN.GT.IEND) THEN
        IBEGIN=1
        IEND=1
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION VTIME()

C#### Function: VTIME
C###  Type: REAL*8
C###  Description:
C###    VTIME returns the current time in seconds.

      IMPLICIT NONE

c     VTIME=DBLE(SECNDS(0.0))
C     VTIME=ITIME(2)*1.0D-3
      VTIME=0.0D0
      RETURN
      END


      INTEGER FUNCTION GETNYR(nc,NPNY,nr,nrc,nrc1,ny,NYNE,NYNP)

C#### Function: GETNYR
C###  Type: INTEGER
C###  Description:
C###    GETNYR returns the corresponding variable nyr with a nrc and
C###    nc for a region nr for a variable ny from nrc1, eg. the
C###    equivalent flux variable.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER nc,NPNY(0:6,NYM,0:NRCM),nr,nrc,nrc1,ny,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
!     Local Variables
      INTEGER na,ne,nh,nk,np,nrr,nv
      CHARACTER ERROR*100

      GETNYR=0
      IF(NPNY(0,ny,nrc1).EQ.1) THEN !variable is nodally based
        nk=NPNY(1,ny,nrc1)
        nv=NPNY(2,ny,nrc1)
        nh=NPNY(3,ny,nrc1)
        np=NPNY(4,ny,nrc1)
        IF(nr.EQ.0) THEN
          nrr=NPNY(6,ny,nrc1)
        ELSE
          nrr=nr
        ENDIF
        GETNYR=NYNP(nk,nv,nh,np,nrc,nc,nrr)
      ELSE IF(NPNY(0,ny,nrc1).EQ.2) THEN !variable is element based
        na=NPNY(1,ny,nrc1)
        nh=NPNY(2,ny,nrc1)
        ne=NPNY(4,ny,nrc1)
        GETNYR=NYNE(na,nh,nrc,nc,ne)
      ELSE
        WRITE(OP_STRING,'('' >>Error in GETNYR: NPNY(0,'',I5,'
     '    //''','',I1,'') is out of range.'')') ny,nrc1
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

      IF(GETNYR.EQ.0) THEN
        WRITE(OP_STRING,'('' >>Error in GETNYR: Could not find nyr'
     '    //' for ny ='',I5)') ny
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


      SUBROUTINE CPCG(ie,nb,NPNE,nr,nx,CE,CG,CP,PG,ERROR,*)

C#### Subroutine: CPCG
C###  Description:
C###    CPCG transfers element parameters CE(il,ne) or interpolates
C###    nodal parameters CP(il,np) with linear basis functions (nb=NMB),
C###    for il=1,ILT(ie,nr,nx), in element ne to return Gauss pt values
C###    CG(il,ng). ILP(il,ie,nr,nx) indicates whether parameter IL is
C###    constant(1),piecewise constant(2) or piecewise linear(3).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER ie,nb,NPNE(NNM,NBFM),nr,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),PG(NSM,NUM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER il,ILMAX,NBL,ng,nk,nm,nn,np,ns
      REAL*8 CPE(64),SUM

      CALL ENTERS('CPCG',*9999)
      DO nm=1,NMM
        DO ng=1,NGM
          CG(nm,ng)=0.0D0
        ENDDO
      ENDDO
      ILMAX=ILT(ie,nr,nx)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' nr='',I2,'' nx='',I2,'' ie='',I4,'
     '    //''' ILMAX='',I4)') nr,nx,ie,ILMAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO il=1,ILMAX
        IF(ABS(ILP(il,ie,nr,nx)).EQ.1.OR.ABS(ILP(il,ie,nr,nx)).EQ.2)
     '    THEN
C       Constant spatially or piecewise constant
C       (Use ABS for KUNLE problem 2-7-91)
          DO ng=1,NGT(nb)
            CG(il,ng)=CE(il)
          ENDDO

        ELSE IF(ABS(ILP(il,ie,nr,nx)).EQ.3) THEN
C       Piecewise linear (nodal parameters)
          NBL=NMB(il,ie,nx)
          ns=0
          DO nn=1,NNT(NBL)
            np=NPNE(nn,NBL)
            DO nk=1,NKT(nn,NBL)
              ns=ns+1
              IF(NKT(nn,NBL).EQ.1) THEN
                CPE(ns)=CP(il,np)
              ELSE IF(NKT(nn,NBL).EQ.2) THEN
!!NOTE: assumes nodes are numbered sequentially - should have nk in CP
                CPE(ns)=CP(il,2*(np-1)+nk)
              ELSE
                ERROR='Only 2 derivs allowed for material param interp'
                GO TO 9999
              ENDIF
            ENDDO
          ENDDO !nn
          DO ng=1,NGT(NBL)
            SUM=0.0D0
            DO ns=1,NST(NBL)
              SUM=SUM+PG(ns,1,ng,NBL)*CPE(ns)
            ENDDO
            CG(il,ng)=SUM
          ENDDO !ng

        ELSE IF(ABS(ILP(il,ie,nr,nx)).EQ.4) THEN
C       Defined by Gauss points
          DO ng=1,NGT(nb)
!            CG(il,ng)=YG()
          ENDDO

        ELSE IF(ABS(ILP(il,ie,nr,nx)).EQ.5) THEN
C       Defined elsewhere
          DO ng=1,NGT(nb)
            CG(il,ng)=CE(il)
          ENDDO
        ENDIF
      ENDDO !il

      IF(DOP) THEN
        DO ng=1,NGT(nb)
          WRITE(OP_STRING,'('' CG(1..,ng='',I3,''): '',6E11.3,'
     '      //'/(17X,6E11.3))') ng,(CG(il,ng),il=1,ILMAX)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF !dop

      CALL EXITS('CPCG')
      RETURN
 9999 CALL ERRORS('CPCG',ERROR)
      CALL EXITS('CPCG')
      RETURN 1
      END


      SUBROUTINE CPXI(ie,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CXI,XI,ERROR,*)

C#### Subroutine: CPXI
C###  Description:
C###    CPXI interpolates element parameters CE(il,ne) or nodal
C###    parameters CP(il,np) with linear basis functions (nb=NMB),
C###    for il=1,ILT(ie,nr,nx), to return Xi position params CXI(il).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER ie,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NPNE(NNM,NBFM),nr,nx
      REAL*8 CE(NMM),CP(NMM,NPM),CXI(NMM),XI(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER il,ILMAX,nb,nk,nn,np,ns
      REAL*8 CPE(64),PSI1

      CALL ENTERS('CPXI',*9999)

      ILMAX=ILT(ie,nr,nx)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' nr='',I2,'' nx='',I2,'' ie='',I4,'
     '    //''' ILMAX='',I4)') nr,nx,ie,ILMAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO il=1,ILMAX
        CXI(il)=0.0d0 !initialise
        IF(ILP(il,1,nr,nx).EQ.1.OR.ILP(il,1,nr,nx).EQ.2) THEN
C         constant spatially or defined by elements
          CXI(il)=CE(il)

        ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
          nb=NMB(il,ie,nx)
C         Put global node material params into element param array CPE
          ns=0
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb)
            DO nk=1,NKT(nn,nb)
              ns=ns+1
              IF(NKT(nn,nb).EQ.1) THEN
                CPE(ns)=CP(il,np)
              ELSE IF(NKT(nn,nb).EQ.2) THEN
C!!!NOTE: assumes nodes are numbered sequentially - should have nk in CP
                CPE(ns)=CP(il,2*(np-1)+nk)
              ELSE
                ERROR='Only 2 derivs allowed for material param interp'
                GO TO 9999
              ENDIF
            ENDDO !nk
          ENDDO !nn
C         Interpolate local element params CPE at Xi
          ns=0
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              ns=ns+1
              CXI(il)=CXI(il)+PSI1(IBT,IDO,INP,nb,1,nk,nn,XI)*CPE(ns)
            ENDDO !nk
          ENDDO !nn

        ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by Gauss points
          CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
     '      //'not implemented',ERROR,*9999)
        ENDIF
      ENDDO !il

      IF(DOP) THEN
        WRITE(OP_STRING,'('' CXI(1..): '',6D12.4,:/(11X,6D12.4))')
     '    (CXI(il),il=1,ILMAX)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF !dop

      CALL EXITS('CPXI')
      RETURN
 9999 CALL ERRORS('CPXI',ERROR)
      CALL EXITS('CPXI')
      RETURN 1
      END


      SUBROUTINE DEFMGRADRC(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DZDX,PG,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: DEFMGRADRC
C###  Description:
C###    DEFMGRADRC calculates components of the deformation
C###    gradient tensor wrt rectangular cartesian coords.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,nr,nx
      REAL*8 DZDX(3,3),PG(NSM,NUM,NGM,NBM),XG(NJM,NUM),XI(3),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mhx,mj,nhx,ni,NITB,nj,NU1(0:3)
      REAL*8 DETERM,DXIXJ(3,3),DXJXI(3,3),dXref_dXrc,DXZ,dZrc_dZref,
     '  DZX,SUM,X(3),Z(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('DEFMGRADRC',*9999)
      NITB=NIT(NBJ(1))

C     Calculate derivatives of Xi wrt Xj (reference) coords, DXIXJ
      DO ni=1,NITB
        DO nj=1,NJ_LOC(NJL_GEOM,0)
          DXJXI(nj,ni)=XG(nj,NU1(ni))
        ENDDO !nj
      ENDDO !ni
      CALL INVERT(NITB,DXJXI,DXIXJ,DETERM)

      IF(ng.EQ.0) THEN
C       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
        CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '    ERROR,*9999)
      ELSE
C       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
        CALL ZEZG(1,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
      ENDIF

C     Put undef/def coords into X/Z for DZX/DXZ function calls below
      DO nhx=1,NJ_LOC(NJL_GEOM,0)
        X(nhx)=XG(nhx,1)
        Z(nhx)=ZG(nhx,1)
      ENDDO !nhx

      DO nhx=1,NJ_LOC(NJL_GEOM,0)
        DO nj=1,NJ_LOC(NJL_GEOM,0)
          SUM=0.0d0
          DO mhx=1,NJ_LOC(NJL_GEOM,0)
            dZrc_dZref=DZX(ITYP11(nr),nhx,mhx,Z)
            DO mj=1,NJ_LOC(NJL_GEOM,0)
              dXref_dXrc=DXZ(ITYP10(nr),mj,nj,X)
              SUM=SUM+dZrc_dZref*ZG(mhx,NU1(mj))*dXref_dXrc
            ENDDO !mj
          ENDDO !mhx
          DZDX(nhx,nj)=SUM
        ENDDO !nj
      ENDDO !nhx

      IF(DOP) THEN
        DO nhx=1,NJ_LOC(NJL_GEOM,0)
          WRITE(OP_STRING,'('' DZDX('',I1,'',nj)   : '',3D12.4)')
     '      nhx,(DZDX(nhx,nj),nj=1,NJ_LOC(NJL_GEOM,0))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nhx
      ENDIF

      CALL EXITS('DEFMGRADRC')
      RETURN
 9999 CALL ERRORS('DEFMGRADRC',ERROR)
      CALL EXITS('DEFMGRADRC')
      RETURN 1
      END


      SUBROUTINE DLZJDX(ICOORD,nb,NJE,DZJDX,XG,ZG,ERROR,*)

C#### Subroutine: DLZJDX
C###  Description:
C###    DLZJDX calculates the covariant derivatives DZJDX of deformed
C###    theta wrt undeformed theta (ie the cpts wrt theta of the
C###    deformation gradient tensor) or wrt undeformed Nu, or wrt Xi
C###    depending on derivatives contained in ZG.

C**** ICOORD=1..5 is coordinate system: Rectangular Cartesian,
C****   Cylindrical Polar, Spherical Polar, Prolate Spheroidal,
C****   Oblate Spheroidal.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER ICOORD,nb,NJE
      REAL*8 DZJDX(3,3),XG(NJM,NUM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,mj,ni,NITB,NU1(0:3)
      REAL*8 CC,CCL,CLX,CLZ,CMX,CMZ,CS,CSL,CT,CX,CZ,DLB,DMB,
     '  DT,DTB,G1,G3,PB,RB,RX,RZ,SC,SCL,SLX,SLZ,SMX,SMZ,SS,SSL,
     '  ST,SX,SZ,TB

      DATA NU1/1,2,4,7/

      CALL ENTERS('DLZJDX',*9999)
      NITB=NIT(nb)
      IF(ICOORD.EQ.1) THEN
        DO mj=1,NJE
          DO ni=1,NITB
            DZJDX(mj,ni)=ZG(mj,NU1(ni))
          ENDDO !ni
        ENDDO !mj

      ELSE IF(ICOORD.EQ.2) THEN
        RX=XG(1,1)
        RZ=ZG(1,1)
        DT=ZG(2,1)-XG(2,1)
        CT=DCOS(DT)
        ST=DSIN(DT)
        DO ni=1,NITB
          DZJDX(1,ni)= ZG(1,NU1(ni))*CT-RZ*ST*ZG(2,NU1(ni))
          DZJDX(2,ni)=(ZG(1,NU1(ni))*ST+RZ*CT*ZG(2,NU1(ni)))/RX
          IF(NJE.EQ.3) DZJDX(3,ni)= ZG(3,NU1(ni))
        ENDDO !ni

      ELSE IF(ICOORD.EQ.3) THEN
        RX=XG(1,1)
        RZ=ZG(1,1)
        DT=ZG(2,1)-XG(2,1)
        CT=DCOS(DT)
        ST=DSIN(DT)
        CX=DCOS(XG(3,1))
        SX=DSIN(XG(3,1))
        CZ=DCOS(ZG(3,1))
        SZ=DSIN(ZG(3,1))
        CC=CX*CZ
        SS=SX*SZ
        CS=CX*SZ
        SC=SX*CZ
        DO ni=1,NITB
          RB=ZG(1,NU1(ni))
          TB=ZG(2,NU1(ni))
          PB=ZG(3,NU1(ni))
          DZJDX(1,ni)=CC*(RB*CT-RZ*ST*TB)-CS*RZ*CT*PB+SC*RZ*PB+SS*RB
          DZJDX(2,ni)=(RZ*CZ*CT*TB+(RB*CZ-RZ*SZ*PB)*ST)/(RX*CX)
          DZJDX(3,ni)=(CC*RZ*PB+CS*RB+SC*(RZ*ST*TB-RB*CT)+SS*RZ*CT*PB)
     '      /RX
        ENDDO !ni

      ELSE IF(ICOORD.EQ.4) THEN
        SLX=DSINH(XG(1,1))
        SLZ=DSINH(ZG(1,1))
        SMX=DSIN(XG(2,1))
        SMZ=DSIN(ZG(2,1))
        CLX=DSQRT(1.0D0+SLX*SLX)
        CLZ=DSQRT(1.0D0+SLZ*SLZ)
        CMX=DSQRT(1.0D0-SMX*SMX)
        CMZ=DSQRT(1.0D0-SMZ*SMZ)
        DT=ZG(3,1)-XG(3,1)
        CT=DCOS(DT)
        ST=DSIN(DT)
        CCL=CLX*CLZ
        CSL=CLX*SLZ
        SCL=SLX*CLZ
        SSL=SLX*SLZ
        CC =CMX*CMZ
        CS =CMX*SMZ
        SC =SMX*CMZ
        SS =SMX*SMZ
        G1=SLX*SLX+SMX*SMX
        G3=SLX*SLX*SMX*SMX
        DO ni=1,NITB
          DLB=ZG(1,NU1(ni))
          DMB=ZG(2,NU1(ni))
          DTB=ZG(3,NU1(ni))
          DZJDX(1,ni)=(( SSL*CC+CCL*SS*CT)*DLB+(-SCL*CS+CSL*SC*CT)*DMB
     '      -CSL*SS*ST*DTB)/G1
          DZJDX(2,ni)=((-CSL*SC+SCL*CS*CT)*DLB+( CCL*SS+SSL*CC*CT)*DMB
     '      -SSL*CS*ST*DTB)/G1
          DZJDX(3,ni)=(SCL*SS*ST*DLB+SSL*SC*ST*DMB+SSL*SS*CT*DTB)/G3
        ENDDO !ni
      ELSE IF(ICOORD.EQ.5) THEN
      ENDIF

      IF(DOP) THEN
        DO mi=1,NITB
          WRITE(OP_STRING,'('' DZJDX('',I1,'',ni)   : '',3D12.4)')
     '      mi,(DZJDX(mi,ni),ni=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
      ENDIF

      CALL EXITS('DLZJDX')
      RETURN
 9999 CALL ERRORS('DLZJDX',ERROR)
      CALL EXITS('DLZJDX')
      RETURN 1
      END


      SUBROUTINE DXIDNU(nb,nr,DXIXN,DXNXI,GL,GU,XG,ERROR,*)

C#### Subroutine: DXIDNU
C###  Description:
C###    DXIDNU evaluates derivatives (DXIXN) of Xi- wrt
C###    undeformed Nu(fibre)-coords, and their inverse DXNXI.

C Rewritten MPN 18-Apr-96
C     MAT_VEC_NG is used to calculate the rectangular cartesian
C     components of the undeformed anatomical material vectors:
C       DXDXN(k,1) is the undef fibre direction vector in rc coords
C       DXDXN(k,2) is the undef sheet direction vector in rc coords
C       DXDXN(k,3) is the undef sheet-normal dirn vector in rc coords

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER nb,nr
      REAL*8 DXIXN(3,3),DXNXI(3,3),GL(3,3),GU(3,3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,mjj,ni,ni2,NITB,njj,NU1(0:3)
      REAL*8 DETERM,DXDXN(3,3),dXrc_dXref,DZX,SUM,X(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('DXIDNU',*9999)
      NITB=NIT(nb)

C     Compute undeformed anatomical fibre vectors wrt rc coordinates
      CALL MAT_VEC_NG(NITB,nr,DXDXN(1,1),DXDXN(1,2),DXDXN(1,3),XG,
     '  ERROR,*9999)

C     Put undeformed coords into X for DZX function call below
      DO njj=1,NJ_LOC(NJL_GEOM,0)
        X(njj)=XG(njj,1)
      ENDDO !njj

C     Calc derivatives of Xi wrt undeformed Nu
      DO ni=1,NITB
        DO mi=1,NITB
          SUM=0.0d0
          DO ni2=1,NITB
            DO njj=1,NJ_LOC(NJL_GEOM,0)
              DO mjj=1,NJ_LOC(NJL_GEOM,0)
                dXrc_dXref=DZX(ITYP10(nr),mjj,njj,X)
                SUM=SUM+GU(ni,ni2)*dXrc_dXref*XG(njj,NU1(ni2))*
     '            DXDXN(mjj,mi)
              ENDDO !mjj
            ENDDO !njj
          ENDDO !ni2
          DXIXN(ni,mi)=SUM
        ENDDO !mi
      ENDDO !ni

      CALL INVERT(NITB,DXIXN,DXNXI,DETERM)

      IF(DOP) THEN
        DO mi=1,NITB
          WRITE(OP_STRING,'('' DXIXN('',I1,'',ni): '',3D12.4)')
     '      mi,(DXIXN(mi,ni),ni=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
      ENDIF

C old way of handling fibre/sheets
C**** 13Sep88: Nu(1) lies in the Xi(1)-Xi(2) plane at angle eta (fibre
C****          angle) to Xi(1). If JTYP9=0 fibre angle is 0 by default.
C****          Nu(2) is orthog to the fibre coord & lies in the Xi(1)-
C****          Xi(2) plane & the remaining Nu coord is orthog to
C****          this plane. The fibre angle eta is eta(1), the angle
C****          between the Nu(1) (fibre) coord & the Xi(1) axis.
C****          The base vectors of the Nu coords have unit length.
C****          Whether the derivs are of deformed or undeformed nu-
C****          coords is determined by the metric tensors with respect
C****          to Xi (GL and GU).
C      REAL*8 C1,C2,DETERM,DXDNU(3,3),DXDT(3,3),ETA1,ETA2,
C     '  RD,RDSQ,RGU33,RK1,RK2
C      IF(NITB.EQ.1) THEN
C        DXINU(1,1)=DSQRT(GL(1,1)) !arclength along 1D element
C
C      ELSE IF(NITB.EQ.2) THEN
C        IF(JTYP9.EQ.0) THEN
C          ETA1=0.0D0
C          ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C          C1=1.0D0
C          C2=DCOS(ETA2)
C        ELSE IF(JTYP9.GE.1) THEN !fibre angle defined
C          IF(JTYP12.EQ.1) THEN !fibre angle defined wrt Xi(1) coord
C            ETA1=XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C          ELSE IF(JTYP12.EQ.2) THEN !fibre angle defined wrt Xi(2) coord
C            ETA2=-XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C          ENDIF
C          C1=DCOS(ETA1)
C          C2=DCOS(ETA2)
C        ENDIF
C        RK1=DSQRT(GL(1,1))*C1
C        RK2=DSQRT(GL(2,2))*C2
C        RDSQ=(GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2))
C        RD=DSQRT(RDSQ)
C        RGU33=DSQRT(GU(3,3))
C        DXINU(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C        DXINU(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C        DXINU(3,1)= 0.0D0
C        DXINU(1,2)=-RK2/RD
C        DXINU(2,2)= RK1/RD
C        DXINU(3,2)= 0.0D0
C        DXINU(1,3)= 0.0D0
C        DXINU(2,3)= 0.0D0
C        DXINU(3,3)= RGU33
C
C      ELSE IF(NITB.EQ.3) THEN
C        IF(JTYP9.EQ.0) THEN
C          ETA1=0.0D0
C          ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C          C1=1.0D0
C          C2=DCOS(ETA2)
C        ELSE IF(JTYP9.GE.1) THEN !fibre angle defined
C          IF(JTYP12.EQ.1) THEN !fibre angle defined wrt Xi(1) coord
C            ETA1=XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C          ELSE IF(JTYP12.EQ.2) THEN !fibre angle defined wrt Xi(2) coord
C            ETA2=-XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C          ENDIF
C          C1=DCOS(ETA1)
C          C2=DCOS(ETA2)
C        ENDIF
C        RK1=DSQRT(GL(1,1))*C1
C        RK2=DSQRT(GL(2,2))*C2
C        RDSQ=(GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2))
C        RD=DSQRT(RDSQ)
C        RGU33=DSQRT(GU(3,3))
C        DXINU(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C        DXINU(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C        DXINU(3,1)= 0.0D0
C        DXINU(1,2)=-RK2/RD
C        DXINU(2,2)= RK1/RD
C        DXINU(3,2)= 0.0D0
C        DXINU(1,3)=RGU33*(GL(1,2)*GL(2,3)-GL(2,2)*GL(1,3))/RDSQ
C        DXINU(2,3)=RGU33*(GL(1,2)*GL(1,3)-GL(1,1)*GL(2,3))/RDSQ
C        DXINU(3,3)=RGU33
C      ENDIF
C      CALL INVERT(NITB,DXINU,DNUXI,DETERM)

      CALL EXITS('DXIDNU')
      RETURN
 9999 CALL ERRORS('DXIDNU',ERROR)
      CALL EXITS('DXIDNU')
      RETURN 1
      END


      SUBROUTINE DXIDZN(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: DXIDZN
C###  Description:
C###    DXIDZN evaluates derivatives (DXIZN) of Xi- wrt
C###    deformed Nu(fibre)-coords, and their inverse DZNXI.

C Written MPN 18-Apr-96
C     MAT_VEC_DEF is used to calculate the rectangular cartesian
C     components of the deformed anatomical material vectors:
C       DZDNU(k,1) is the deformed fibre direction vector in rc coords
C       DZDNU(k,2) is the deformed sheet direction vector in rc coords
C       DZDNU(k,3) is the deformed sheet-normal dirn vector in rc coords

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,nr,nx
      REAL*8 DXIZN(3,3),DZNXI(3,3),
     '  PG(NSM,NUM,NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),XI(3),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,mhx,ni,ni2,nhx,NITB,NU1(0:3)
      REAL*8 DETERM,DXIX(3,3),DZDNU(3,3),dZrc_dZref,DZX,
     '  GZ,GZL(3,3),GZU(3,3),
     '  SUM,Z(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('DXIDZN',*9999)
      NITB=NIT(NBJ(1))

      CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '  PG,XE,XG,XI,ZE,ZG,ERROR,*9999)

      IF(ng.EQ.0) THEN
C ***   Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
        CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG,XI,
     '    ERROR,*9999)
      ELSE
C ***   Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
      ENDIF
C *** Calculate deformed metric tensors wrt Xi (GZL,GZU)
      CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)

C     Put def coords into Z for DZX function call below
      DO nhx=1,NJ_LOC(NJL_GEOM,0)
        Z(nhx)=ZG(nhx,1)
      ENDDO !nhx

C     Calc derivs of Xi wrt deformed Nu
      DO ni=1,NITB
        DO mi=1,NITB
          SUM=0.0d0
          DO ni2=1,NITB
            DO nhx=1,NJ_LOC(NJL_GEOM,0)
              DO mhx=1,NJ_LOC(NJL_GEOM,0)
                dZrc_dZref=DZX(ITYP11(nr),mhx,nhx,Z)
                SUM=SUM+GZU(ni,ni2)*dZrc_dZref*ZG(nhx,NU1(ni2))*
     '            DZDNU(mhx,mi)
              ENDDO !mhx
            ENDDO !nhx
          ENDDO !ni2
          DXIZN(ni,mi)=SUM
        ENDDO !mi
      ENDDO !ni

      CALL INVERT(NITB,DXIZN,DZNXI,DETERM)

      IF(DOP) THEN
        DO mi=1,NITB
          WRITE(OP_STRING,'('' DXIZN('',I1,'',ni): '',3D12.4)')
     '      mi,(DXIZN(mi,ni),ni=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
      ENDIF

      CALL EXITS('DXIDZN')
      RETURN
 9999 CALL ERRORS('DXIDZN',ERROR)
      CALL EXITS('DXIDZN')
      RETURN 1
      END


      SUBROUTINE FIBRE_REF_VECS(NITB,nr,F_VECTOR,G_VECTOR,H_VECTOR,
     '  XG,ERROR,*)

C#### Subroutine: FIBRE_REF_VECS
C###  Description:
C###    FIBRE_REF_VECS calculates direction cosines of undeformed
C###    fibre reference vectors at Gauss point ng.
C###    F_VECTOR coincides with the Xi1 base vector
C###    G_VECTOR lies in the Xi1-Xi2 plane and is normal to Xi1
C###    H_VECTOR is normal to Xi1-Xi2 plane

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 F_VECTOR(3),G_VECTOR(3),H_VECTOR(3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,nj1,nj2,NU1(0:3)
      REAL*8 dXRC_dXI(3,3),dXrc_dXref,DZX,SUM,X(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('FIBRE_REF_VECS',*9999)

C     Put undef coords into X for DZX function call below
C     and initialise all vectors
      DO nj1=1,3
        IF(nj1.LE.NJ_LOC(NJL_GEOM,0)) THEN
          X(nj1)=XG(nj1,1)
        ELSE
          X(nj1)=0.0d0
        ENDIF
        F_VECTOR(nj1)=0.0d0
        G_VECTOR(nj1)=0.0d0
        H_VECTOR(nj1)=0.0d0
        DO ni=1,3
          dXRC_dXI(nj1,ni)=0.0d0
        ENDDO !ni
      ENDDO !nj1

C     Compute derivatives of rc coords wrt Xi coords
      DO ni=1,NITB
        DO nj1=1,NJ_LOC(NJL_GEOM,0)
          SUM=0.0d0
          DO nj2=1,NJ_LOC(NJL_GEOM,0)
            dXrc_dXref=DZX(ITYP10(nr),nj1,nj2,X)
            SUM=SUM+dXrc_dXref*XG(nj2,NU1(ni))
          ENDDO !nj2
          dXRC_dXI(nj1,ni)=SUM
        ENDDO !nj1
      ENDDO !ni

C     F_VECTOR is the normalised undeformed Xi1 base vector
      DO nj1=1,NJ_LOC(NJL_GEOM,0)
        F_VECTOR(nj1)=dXRC_dXI(nj1,1)
      ENDDO !nj1
      CALL NORMALISE(3,F_VECTOR,ERROR,*9999)

      IF(NITB.GE.2) THEN !2D or 3D
C       H_VECTOR is the undeformed Xi1-Xi2 plane normal
        CALL CROSS(dXRC_dXI(1,1),dXRC_dXI(1,2),H_VECTOR)
        CALL NORMALISE(3,H_VECTOR,ERROR,*9999)
C       G_VECTOR lies in the undeformed Xi1-Xi2 plane and is
C       normal to F_VECTOR
        CALL CROSS(H_VECTOR,F_VECTOR,G_VECTOR)
        CALL NORMALISE(3,G_VECTOR,ERROR,*9999)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Undeformed fibre reference vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    f_vector    g_vector    h_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO ni=1,3
          WRITE(OP_STRING,'(X,3D12.3)')
     '      F_VECTOR(ni),G_VECTOR(ni),H_VECTOR(ni)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
      ENDIF

      CALL EXITS('FIBRE_REF_VECS')
      RETURN
 9999 CALL ERRORS('FIBRE_REF_VECS',ERROR)
      CALL EXITS('FIBRE_REF_VECS')
      RETURN 1
      END


      SUBROUTINE FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,ng,
     '  NHE,NITB,nr,nx,
     '  FD_VECTOR,GD_VECTOR,HD_VECTOR,PG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: FIBRE_REF_VECS_DEF
C###  Description:
C###    FIBRE_REF_VECS_DEF calculates direction cosines of deformed
C###    fibre reference vectors at Gauss point ng.
C###    FD_VECTOR coincides with the deformed Xi1 base vector
C###    GD_VECTOR lies in the def Xi1-Xi2 plane and is normal to Xi1
C###    HD_VECTOR is normal to deformed Xi1-Xi2 plane

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),ng,NHE,NITB,nr,nx
      REAL*8 FD_VECTOR(3),GD_VECTOR(3),HD_VECTOR(3),PG(NSM,NUM,NGM,NBM),
     '  XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mhx,ni,nhx,NU1(0:3)
      REAL*8 DXIX(3,3),dZRC_dXI(3,3),dZrc_dZref,DZX,SUM,Z(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('FIBRE_REF_VECS_DEF',*9999)

C     Put def coords into Z for DZX function call below
C     and initialise all vectors
      DO nhx=1,3
        IF(nhx.LE.NJ_LOC(NJL_GEOM,0)) THEN
          Z(nhx)=ZG(nhx,1)
        ELSE
          Z(nhx)=0.0d0
        ENDIF
        FD_VECTOR(nhx)=0.0d0
        GD_VECTOR(nhx)=0.0d0
        HD_VECTOR(nhx)=0.0d0
        DO ni=1,3
          dZRC_dXI(nhx,ni)=0.0d0
        ENDDO !ni
      ENDDO !nhx

C     Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
      IF(ng.EQ.0) THEN
        CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG,XI,
     '    ERROR,*9999)
      ELSE
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
      ENDIF

C     Compute derivatives of deformed rc coords wrt Xi coords
      DO ni=1,NITB
        DO nhx=1,NJ_LOC(NJL_GEOM,0)
          SUM=0.0d0
          DO mhx=1,NJ_LOC(NJL_GEOM,0)
            dZrc_dZref=DZX(ITYP11(nr),nhx,mhx,Z)
            SUM=SUM+dZrc_dZref*ZG(mhx,NU1(ni))
          ENDDO !mhx
          dZRC_dXI(nhx,ni)=SUM
        ENDDO !nhx
      ENDDO !ni

C     FD_VECTOR is the normalised deformed Xi1 base vector
      DO nhx=1,NJ_LOC(NJL_GEOM,0)
        FD_VECTOR(nhx)=dZRC_dXI(nhx,1)
      ENDDO !nhx
      CALL NORMALISE(3,FD_VECTOR,ERROR,*9999)

      IF(NITB.GE.2) THEN !2D or 3D
C       HD_VECTOR is the undeformed Xi1-Xi2 plane normal
        CALL CROSS(dZRC_dXI(1,1),dZRC_dXI(1,2),HD_VECTOR)
        CALL NORMALISE(3,HD_VECTOR,ERROR,*9999)
C       GD_VECTOR lies in the deformed Xi1-Xi2 plane and is
C       normal to FD_VECTOR
        CALL CROSS(HD_VECTOR,FD_VECTOR,GD_VECTOR)
        CALL NORMALISE(3,GD_VECTOR,ERROR,*9999)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Deformed fibre reference vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    fd_vector   gd_vector   hd_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO ni=1,3
          WRITE(OP_STRING,'(X,3D12.3)')
     '      FD_VECTOR(ni),GD_VECTOR(ni),HD_VECTOR(ni)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
      ENDIF

      CALL EXITS('FIBRE_REF_VECS_DEF')
      RETURN
 9999 CALL ERRORS('FIBRE_REF_VECS_DEF',ERROR)
      CALL EXITS('FIBRE_REF_VECS_DEF')
      RETURN 1
      END


      SUBROUTINE MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  AD_VECTOR,BD_VECTOR,CD_VECTOR,
     '  PG,XE,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: MAT_VEC_DEF
C###  Description:
C###    MAT_VEC_DEF calculates direction cosines of deformed
C###    normalised material vectors at Gauss point ng
C###    or at XI (if ng=0).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,nr,nx
      REAL*8 AD_VECTOR(3),BD_VECTOR(3),CD_VECTOR(3),PG(NSM,NUM,NGM,NBM),
     '  XE(NSM,NJM),XG(NJM,NUM),XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mj,ni,NITB,nj
      REAL*8 DXDNU(3,3),DZDNU(3,3),DZDX(3,3),SUM

      CALL ENTERS('MAT_VEC_DEF',*9999)
      NITB=NIT(NBJ(1))

      IF(ng.EQ.0) THEN
C       Compute undeformed anatomical fibre vectors wrt rc coords at XI
        CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '    DXDNU(1,1),DXDNU(1,2),DXDNU(1,3),XE,XG,XI,ERROR,*9999)
      ELSE
C       Compute undeformed anatomical fibre vectors wrt rc coords at ng
        CALL MAT_VEC_NG(NITB,nr,DXDNU(1,1),DXDNU(1,2),DXDNU(1,3),XG,
     '    ERROR,*9999)
      ENDIF
C     Initialise deformed material vectors with undef material vec dirns
      DO nj=1,3
        AD_VECTOR(nj)=DXDNU(nj,1)
        BD_VECTOR(nj)=DXDNU(nj,2)
        CD_VECTOR(nj)=DXDNU(nj,3)
      ENDDO !nj

C     Compute the deformation gradient tensor wrt rc coords
C     ie want derivatives of deformed rc coordinates wrt
C     undeformed rc coordinates.
      CALL DEFMGRADRC(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DZDX,PG,XG,XI,ZE,ZG,ERROR,*9999)

C     Compute deformed material vectors wrt rc coordinates using
C     the deformation gradient tensor and the undeformed material
C     vectors wrt rc coordinates
      DO nj=1,NJ_LOC(NJL_GEOM,0)
        DO ni=1,NITB
          SUM=0.0d0
          DO mj=1,NJ_LOC(NJL_GEOM,0)
            SUM=SUM+DZDX(nj,mj)*DXDNU(mj,ni)
          ENDDO !nj1
          DZDNU(nj,ni)=SUM
        ENDDO !ni
        AD_VECTOR(nj)=DZDNU(nj,1)
        BD_VECTOR(nj)=DZDNU(nj,2)
        IF(NITB.EQ.3) CD_VECTOR(nj)=DZDNU(nj,3)
      ENDDO !nj

C     Normalise the deformed vectors so that deformed material
C     coordinates are measures of physical arc length
      CALL NORMALISE(3,AD_VECTOR,ERROR,*9999)
      CALL NORMALISE(3,BD_VECTOR,ERROR,*9999)
      CALL NORMALISE(3,CD_VECTOR,ERROR,*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Deformed material vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    ad_vector   bd_vector   cd_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj=1,3
          WRITE(OP_STRING,'(X,3D12.3)')
     '      AD_VECTOR(nj),BD_VECTOR(nj),CD_VECTOR(nj)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nj
      ENDIF

      CALL EXITS('MAT_VEC_DEF')
      RETURN
 9999 CALL ERRORS('MAT_VEC_DEF',ERROR)
      CALL EXITS('MAT_VEC_DEF')
      RETURN 1
      END


      SUBROUTINE MAT_VEC_NG(NITB,nr,A_VECTOR,B_VECTOR,C_VECTOR,XG,
     '  ERROR,*)

C#### Subroutine: MAT_VEC_NG
C###  Description:
C###    MAT_VEC_NG calculates direction cosines of undeformed
C###    material vectors at Gauss point ng.

C     MPN rewritten 25-Apr-96

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      REAL*8 FIBRE_ORIENT(3,3)

      CALL ENTERS('MAT_VEC_NG',*9999)

C     Compute components of undeformed fibre reference vectors
C     at Gauss pt NG wrt rc coord system
      CALL FIBRE_REF_VECS(NITB,nr,FIBRE_ORIENT(1,1),
     '  FIBRE_ORIENT(1,2),FIBRE_ORIENT(1,3),XG,ERROR,*9999)

C     Rotate fibre reference vectors into true material fibre vectors
      CALL MAT_VEC_ROTATE(NITB,nr,FIBRE_ORIENT,XG,ERROR,*9999)
      DO nj=1,3
        A_VECTOR(nj)=FIBRE_ORIENT(nj,1)
        B_VECTOR(nj)=FIBRE_ORIENT(nj,2)
        C_VECTOR(nj)=FIBRE_ORIENT(nj,3)
      ENDDO

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Undeformed material vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    a_vector    b_vector    c_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj=1,3
          WRITE(OP_STRING,'(X,3D12.3)')
     '      A_VECTOR(nj),B_VECTOR(nj),C_VECTOR(nj)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
      ENDIF

C old MPN 25-Apr-96
CC**** Delta is angle between h (wall normal) and v (normal to wall in
CC**** yz-plane).
CC**** Phi is angle between v and r.
CC**** Note: If Beta>0 A_VECTOR is calculated with Beta terms included.
C      REAL*8 ALFA,ARG1,ARG2,BETA,GAMA,cosh_LAMDA,
C     '  cos_ALFA,cos_BETA,cos_GAMA,cos_DELTA,cos_PHI,cos_MU,
C     '  cos_THETA,coth_LAMDA,cot_MU,
C     '  D1,D2,D3,DELTA,dLAMDA_dXi1,dLAMDA_dXi2,dMU_dXi1,dMU_dXi2,
C     '  dTHETA_dXi1,dX_dXI1(3),PHI,RLAMDA,RMU,sinh_LAMDA,
C     '  sin_ALFA,sin_BETA,sin_GAMA,sin_DELTA,sin_PHI,sin_MU,sin_THETA,
C     '  SUMSQ,THETA,tanh_LAMDA,tan_MU
C      IF(NITB.EQ.1) THEN    !1d
C        IF(ITYP10(1).EQ.1) THEN  !r.c. coordinates
C          DO nj=1,3
C            dX_dXI1(nj)=XG(nj,2)
C          ENDDO
C        ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar coordinates
C          THETA=XG(2,1)
C          sin_THETA=DSIN(THETA)
C          cos_THETA=DCOS(THETA)
C          dX_dXI1(1)=XG(1,2)*cos_THETA - XG(1,1)*sin_THETA*XG(2,2)
C          dX_dXI1(2)=XG(1,2)*sin_THETA + XG(1,1)*cos_THETA*XG(2,2)
C          dX_dXI1(3)=XG(3,2)
C  ENDIF
C        SUMSQ=0.0D0
C        DO nj=1,3
C          SUMSQ=SUMSQ+dX_dXI1(nj)*dX_dXI1(nj)
C        ENDDO
C        SUMSQ=DSQRT(SUMSQ)
C        DO nj=1,3
C          IF(SUMSQ.GT.1.D-12) A_VECTOR(nj)=dX_dXI1(nj)/SUMSQ
C          B_VECTOR(nj)=0.0d0
C          C_VECTOR(nj)=0.0d0
C        ENDDO
C      ELSE IF(NITB.EQ.2) THEN   !2D
C        IF(ITYP10(1).EQ.1) THEN  !r.c. coordinates only so far
C          nj=NJ_LOC(NJL_FIBR,1)  !Fibre angle position
C         ALFA=XG(nj,1)          !Fibre angle
C         dX_dXI1(1)=XG(1,2)
C         dX_dXI1(2)=XG(2,2)
C         PHI=DATAN2(dX_dXI1(2),dX_dXI1(1)) !Angle betw x and xi
C    A_VECTOR(1)=DCOS(ALFA+PHI)
C    A_VECTOR(2)=DSIN(ALFA+PHI)
C    A_VECTOR(3)=0.0D0
C    B_VECTOR(1)=-DSIN(ALFA+PHI)
C    B_VECTOR(2)=DCOS(ALFA+PHI)
C    B_VECTOR(3)=0.0D0
C    C_VECTOR(1)=0.0D0
C    C_VECTOR(2)=0.0D0
C    C_VECTOR(3)=0.0D0
C        ENDIF
C      ELSE IF(NITB.EQ.3) THEN   !3D
C        IF(ITYP10(1).EQ.1) THEN  !r.c. coordinates
C
C!  Alfa (fibre angle)  : rotation in x,y plane about cross axis
C          nj=NJ_LOC(NJL_FIBR,1)
C          IF(nj.EQ.0) THEN
C            WRITE(OP_STRING,
C     '        '('' >>WARNING !!! Assuming fibre angle of zero'')')
C            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C            ALFA       = 0.0D0
C          ELSE
C            ALFA       = XG(nj,1)
C          ENDIF
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''   ALFA ='',E12.3)') ALFA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_ALFA   = DSIN(ALFA)
C          cos_ALFA   = DCOS(ALFA)
C
C!  Beta (imbrication angle)  : rotation of fibre axis about sheet axis
C          nj=NJ_LOC(NJL_FIBR,2)
C          IF(nj.EQ.0) THEN
C            WRITE(OP_STRING,
C     '        '('' >>WARNING !!! Assuming imbrication angle of zero'')')
C            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C            BETA       = 0.0D0
C          ELSE
C            BETA       = XG(nj,1)
C          ENDIF
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''   BETA ='',E12.3)') BETA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_BETA   = DSIN(BETA)
C          cos_BETA   = DCOS(BETA)
C
C!  Gamma (sheet angle)  : rotation of sheet axis about fibre axis
C          nj=NJ_LOC(NJL_FIBR,3)
C          IF(nj.EQ.0) THEN
C            WRITE(OP_STRING,
C     '        '('' >>WARNING !!! Assuming sheet angle of zero'')')
C            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C            GAMA       = 0.0D0
C          ELSE
C            GAMA       = XG(nj,1)
C          ENDIF
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''  GAMMA ='',E12.3)') GAMA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_GAMA   = DSIN(GAMA)
C          cos_GAMA   = DCOS(GAMA)
CC GBS Removed 22-MAR-1996
Cc          WRITE(OP_STRING,
Cc     '      '('' >>WARNING !!! Glen has hacked this to be orthog'')')
Cc          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
Cc          sin_GAMA   = 1.0D0
Cc          cos_GAMA   = 0.0D0
C          A_VECTOR(1) =  cos_ALFA*cos_BETA + sin_ALFA*sin_BETA*sin_GAMA
C          A_VECTOR(2) =  sin_ALFA*cos_BETA - cos_ALFA*sin_GAMA*sin_BETA
C          A_VECTOR(3) = -cos_GAMA*sin_BETA
C
C          B_VECTOR(1) =  cos_ALFA*sin_BETA - sin_ALFA*sin_GAMA*cos_BETA
C          B_VECTOR(2) =  sin_ALFA*sin_BETA + cos_ALFA*sin_GAMA*cos_BETA
C          B_VECTOR(3) =  cos_GAMA*cos_BETA
C
C     C_VECTOR(1) =  sin_ALFA*cos_GAMA
C      C_VECTOR(2) = -cos_ALFA*cos_GAMA
C      C_VECTOR(3) =  sin_GAMA
CC GMH this is hacked in so that we have orthogonal axes
CC GBS Removed 22-MAR-1996
Cc          CALL CROSS(A_VECTOR,B_VECTOR,C_VECTOR)
C
C        ELSE IF(ITYP10(1).EQ.4) THEN !prolate coordinates
C
C!  Lamda
C          RLAMDA     = XG(1,1)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'('' RLAMDA ='',E12.3)') RLAMDA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          dLAMDA_dXi1= XG(1,2)
C          dLAMDA_dXi2= XG(1,4)
C          sinh_LAMDA = DSINH(RLAMDA)
C          cosh_LAMDA = DCOSH(RLAMDA)
C          IF(DABS(cosh_LAMDA).GT.1.0D-6) THEN
C            tanh_LAMDA = sinh_LAMDA/cosh_LAMDA
C          ELSE
C            tanh_LAMDA = 1.0D8
C          ENDIF
C          IF(DABS(sinh_LAMDA).GT.1.0D-6) THEN
C            coth_LAMDA = cosh_LAMDA/sinh_LAMDA
C          ELSE
C            coth_LAMDA = 1.0D8
C          ENDIF
C
C!  Mu
C          RMU        = XG(2,1)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''    RMU ='',E12.3)') RMU
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          dMU_dXi1   = XG(2,2)
C          dMU_dXi2   = XG(2,4)
C          sin_MU     = DSIN(RMU)
C          cos_MU     = DCOS(RMU)
C          IF(DABS(cos_MU).GT.1.0D-6) THEN
C            tan_MU = sin_MU/cos_MU
C          ELSE
C            tan_MU = 1.0D8
C          ENDIF
C          IF(DABS(sin_MU).GT.1.0D-6) THEN
C            cot_MU     = cos_MU/sin_MU
C          ELSE
C            cot_MU = 1.0D8
C          ENDIF
C
C!  Theta
C          THETA      = XG(3,1)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''  THETA ='',E12.3)') THETA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          dTHETA_dXi1= XG(3,2)
C          sin_THETA  = DSIN(THETA)
C          cos_THETA  = DCOS(THETA)
C
C!  Alfa (fibre angle)
C          nj=NJ_LOC(NJL_FIBR,1)
C          ALFA       = XG(nj,1)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''   ALFA ='',E12.3)') ALFA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_ALFA   = DSIN(ALFA)
C          cos_ALFA   = DCOS(ALFA)
C
C!  Beta (imbrication angle)
C          nj=NJ_LOC(NJL_FIBR,2)
C          BETA       = XG(nj,1)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''   BETA ='',E12.3)') BETA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_BETA   = DSIN(BETA)
C          cos_BETA   = DCOS(BETA)
C
C!  Gamma (sheet angle)
C          nj=NJ_LOC(NJL_FIBR,3)
C          GAMA       = XG(nj,1)
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''   GAMA ='',E12.3)') GAMA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_GAMA   = DSIN(GAMA)
C          cos_GAMA   = DCOS(GAMA)
C
C!  Phi (eqtn 8 in "Laminar structure of the heart II")
C          ARG1=coth_LAMDA*dLAMDA_dXi1+cot_MU*dMU_dXi1
C          ARG2=dTHETA_dXi1
C          PHI =DATAN2(ARG1,ARG2)
C          IF(PHI.GT.PI/2.0D0) THEN
C            PHI=PHI-PI
C          ELSE IF(PHI.LT.-PI/2.0D0) THEN
C            PHI=PHI+PI
C          ENDIF
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''    PHI ='',E12.3)') PHI
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_PHI   = DSIN(PHI)
C          cos_PHI   = DCOS(PHI)
C
C!  Delta (eqtn 12 in "Laminar structure of the heart II")
C          ARG1=(tan_MU*dLAMDA_dXi2+tanh_LAMDA*dMU_dXi2)*cos_PHI
C          ARG2=tanh_LAMDA*dLAMDA_dXi2-tan_MU*dMU_dXi2
C          DELTA=DATAN2(ARG1,ARG2)
C          IF(DELTA.GT.PI/2.0D0) THEN
C            DELTA=DELTA-PI
C          ELSE IF(DELTA.LT.-PI/2.0D0) THEN
C            DELTA=DELTA+PI
C          ENDIF
C          IF(DOP) THEN
C            WRITE(OP_STRING,'(''  DELTA ='',E12.3)') DELTA
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          sin_DELTA  = DSIN(DELTA)
C          cos_DELTA  = DCOS(DELTA)
C
C!  A vector (eqtn 3 in "Laminar structure of the heart II")
C          IF(DABS(BETA).LT.1.0D-5) THEN !no imbrication angle
C            A_VECTOR(1)=-sin_ALFA*cos_DELTA
C            A_VECTOR(2)=-cos_ALFA*sin_PHI-sin_ALFA*sin_DELTA*cos_PHI
C            A_VECTOR(3)=-cos_ALFA*cos_PHI+sin_ALFA*sin_DELTA*sin_PHI
C          ELSE !include imbrication angle
C            D1=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*cos_DELTA
C     '         + sin_BETA*cos_GAMA*sin_DELTA
C            D2=  cos_ALFA*cos_BETA+sin_ALFA*sin_BETA*sin_GAMA
C            D3=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*sin_DELTA
C     '         - sin_BETA*cos_GAMA*cos_DELTA
C            A_VECTOR(1)= D1
C            A_VECTOR(2)=-D2*sin_PHI+D3*cos_PHI
C            A_VECTOR(3)=-D2*cos_PHI-D3*sin_PHI
C          ENDIF
C          CALL ROTATION(1,THETA,A_VECTOR,ERROR,*9999)
C
C!  B vector (eqtn 4 in "Laminar structure of the heart II")
C          ARG1=cos_ALFA*sin_GAMA*sin_DELTA-cos_GAMA*cos_DELTA
C          B_VECTOR(1)=-cos_ALFA*sin_GAMA*cos_DELTA - cos_GAMA*sin_DELTA
C          B_VECTOR(2)= sin_ALFA*sin_GAMA*sin_PHI - ARG1*cos_PHI
C          B_VECTOR(3)= sin_ALFA*sin_GAMA*cos_PHI + ARG1*sin_PHI
C          CALL ROTATION(1,THETA,B_VECTOR,ERROR,*9999)
C
C!  C vector (eqtn 5 in "Laminar structure of the heart II")
C          CALL CROSS(A_VECTOR,B_VECTOR,C_VECTOR)
C        ENDIF
C      ENDIF

      CALL EXITS('MAT_VEC_NG')
      RETURN
 9999 CALL ERRORS('MAT_VEC_NG',ERROR)
      CALL EXITS('MAT_VEC_NG')
      RETURN 1
      END


      SUBROUTINE MAT_VEC_ROTATE(NITB,nr,FIBRE_ORIENT,XG,ERROR,*)

C#### Subroutine: MAT_VEC_ROTATE
C###  Description:
C###    MAT_VEC_ROTATE rotates the fibre reference vector orientations
C###    into the true material fibre directions using the fibre, sheet
C###    and sheet-normal angles stored in XG.

C     MPN 25-Apr-96: written to reflect Laminar Structure
C                    of the Heart paper

C#### Comment: FIBRE, IMBRICATION and SHEET ANGLES
C###  Description:
C###    Eqtn 5 of "Laminar Structure of the Heart" by
C###    Le Grice, Hunter and Smaill defines how the fibre reference
C###    vectors F_VECTOR, G_VECTOR and H_VECTOR are
C###    computed from the microstructural material axis vectors
C###    A_VECTOR, B_VECTOR and C_VECTOR. Three successive coordinate
C###    rotations define the transformation (the order IS important).
C###    Firstly, the material vectors are rotated by GAMA-PI/2 about
C###    the fibre axis until the new sheet axis lies in the
C###    Xi1-Xi2 plane. Secondly, the resulting vectors are rotated
C###    by BETA about the new sheet axis until the new fibre axis also
C###    lies in the Xi1-Xi2 plane. Lastly, the resulting vectors are
C###    rotated by ALFA about the new sheet-normal axis (which now
C###    coincides with the Xi1-Xi2 plane normal) until the new fibre
C###    axis coincides with the Xi1 base vector.

C**** ALFA is fibre angle.
C**** BETA is imbrication angle.
C**** GAMA is sheet angle.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 XG(NJM,NUM),FIBRE_ORIENT(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      REAL*8 ALFA,BETA,GAMA

      CALL ENTERS('MAT_VEC_ROTATE',*9999)

C     To compute undeformed material vectors from the fibre reference
C     apply the inverse of the transformation described above.
      IF(JTYP9.GE.1) THEN !fibre angles defined
        nj=NJ_LOC(NJL_FIBR,1)
        ALFA=XG(nj,1) !fibre angle subtends the Xi1 base vector
C                     !and fibre axis in the Xi1-Xi2 plane.
C                     !ie) angle that rotates the Xi1 base vector
C                     !(about the sheet-normal axis) into the
C                     !direction of the fibre axis given that
C                     !the sheet axis lies in the Xi1-Xi2 plane
C                     !and the sheet-normal coincides with the
C                     !Xi1-Xi2 plane normal.
C       Rotate fibre ref vectors about the sheet-normal axis by ALFA
        CALL ROT_COORDSYS(3,ALFA,FIBRE_ORIENT,ERROR,*9999)
      ELSE
        IF(DOP.AND.NITB.GE.2) THEN
          WRITE(OP_STRING,'('' >>>WARNING: Assuming fibre '
     '      //'angle is zero'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      IF(JTYP9.GE.2) THEN !imbrication angles defined
        nj=NJ_LOC(NJL_FIBR,2)
        BETA=XG(nj,1) !imbrication angle subtends the fibre axis
C                     !and the Xi1-Xi2 plane (about the sheet axis)
C                     !given that the sheet axis lies in
C                     !the Xi1-Xi2 plane
C       Rotate resulting vectors about the new sheet axis by BETA
        CALL ROT_COORDSYS(2,BETA,FIBRE_ORIENT,ERROR,*9999)
      ELSE
        IF(DOP.AND.NITB.EQ.3) THEN
          WRITE(OP_STRING,'('' >>>WARNING: Assuming imbrication '
     '      //'angle is zero'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      IF(JTYP9.GE.3) THEN !sheet angles defined
        nj=NJ_LOC(NJL_FIBR,3)
        GAMA=XG(nj,1) !sheet angle-PI/2 subtends sheet axis and the
C                     !Xi1-Xi2 plane about the true fibre axis
C!!! TEMPORARY
CC         Rotate resulting vectors about the fibre axis by GAMA
C          CALL ROT_COORDSYS(1,GAMA,FIBRE_ORIENT,ERROR,*9999)

C!!!      NOTE: this is temporary until we redefine the way the sheet
C!!!      angle is defined. Better to rotate by GAMA above and not
C!!!      GAMA-PI/2. This means the sheet angles for the full
C!!!      heart mesh will need transforming. See also below.
C         Rotate resulting vectors about the fibre axis by GAMA-PI/2
          CALL ROT_COORDSYS(1,GAMA-PI/2.0d0,FIBRE_ORIENT,ERROR,*9999)
      ELSE
        IF(DOP.AND.NITB.EQ.3) THEN
          WRITE(OP_STRING,'('' >>>WARNING: Assuming sheet '
     '      //'angle is zero'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        ENDIF
C!!! TEMPORARY: see comment 12 lines above
        IF(NITB.EQ.3) THEN !only rotate for 3D elements
          CALL ROT_COORDSYS(1,-PI/2.0d0,FIBRE_ORIENT,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('MAT_VEC_ROTATE')
      RETURN
 9999 CALL ERRORS('MAT_VEC_ROTATE',ERROR)
      CALL EXITS('MAT_VEC_ROTATE')
      RETURN 1
      END


      SUBROUTINE MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '  A_VECTOR,B_VECTOR,C_VECTOR,XE,XG,XI,ERROR,*)

C#### Subroutine: MAT_VEC_XI
C###  Description:
C###    MAT_VEC_XI calculates direction cosines of material vectors at
C###    Xi in element ne.

C     MPN rewritten 25-Apr-96

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBJ(NJM),nr
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  XE(NSM,NJM),XG(NJM,NUM),XI(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,nj,NITB
      REAL*8 FIBRE_ORIENT(3,3)

      CALL ENTERS('MAT_VEC_XI',*9999)

      NITB=NIT(NBJ(1))
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Xi coords: '',3D12.3)') (XI(ni),ni=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C     Interpolate midwall geometric vars XG and derivs wrt Xi
      CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)

C     Compute components of undeformed fibre reference vectors at XI
C     wrt rc coord system
      CALL FIBRE_REF_VECS(NITB,nr,FIBRE_ORIENT(1,1),
     '  FIBRE_ORIENT(1,2),FIBRE_ORIENT(1,3),XG,ERROR,*9999)

C     Rotate fibre reference vectors into true material fibre vectors
      CALL MAT_VEC_ROTATE(NITB,nr,FIBRE_ORIENT,XG,ERROR,*9999)
      DO nj=1,3
        A_VECTOR(nj)=FIBRE_ORIENT(nj,1)
        B_VECTOR(nj)=FIBRE_ORIENT(nj,2)
        C_VECTOR(nj)=FIBRE_ORIENT(nj,3)
      ENDDO

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Undeformed material vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    a_vector    b_vector    c_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj=1,3
          WRITE(OP_STRING,'(X,3D12.3)')
     '      A_VECTOR(nj),B_VECTOR(nj),C_VECTOR(nj)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
      ENDIF

C old MPN 25-Apr-96
CC**** Note: XE must be passed in (calculated in previous call to XPXE).
CC**** A_VECTOR is fibre angle vector (in sheet)
CC**** B_VECTOR is sheet angle vector (in sheet orthog to fibres)
CC**** C_VECTOR   is normal to sheet
CC**** Alfa is fibre angle
CC**** Beta is imbrication angle
CC**** Gama is sheet angle
CC**** Delta is angle between h (wall normal) and v (normal to wall in
CC**** yz-plane)
CC**** Phi is angle between v and r
CC**** GAMA,PHI & DELTA are also returned.
CC**** Note: If Beta>0 A_VECTOR is calculated with Beta terms included.
C      REAL*8 ARG1,ARG2,cosh_LAMDA,
C     '  cos_ALFA,cos_BETA,cos_GAMA,cos_DELTA,cos_PHI,cos_MU,
C     '  coth_LAMDA,cot_MU,
C     '  D1,D2,D3,dLAMDA_dXi1,dLAMDA_dXi2,dMU_dXi1,dMU_dXi2,
C     '  dTHETA_dXi1,
C     '  PXI,RLAMDA,RMU,sinh_LAMDA,sin_ALFA,sin_BETA,sin_GAMA,sin_DELTA,
C     '  sin_PHI,sin_MU,THETA,tanh_LAMDA,tan_MU
C      IF(ITYP10(1).EQ.4) THEN !prolate coordinates
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' Xi coords: '',3E12.3)') (XI(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C!  Lamda
C        nb=NBJ(1)
C        RLAMDA     = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,1))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' RLAMDA ='',E12.3)') RLAMDA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        dLAMDA_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,
C     '    XE(1,1))
C        dLAMDA_dXi2= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,XI,
C     '    XE(1,1))
C        sinh_LAMDA = DSINH(RLAMDA)
C        cosh_LAMDA = DCOSH(RLAMDA)
C        IF(DABS(cosh_LAMDA).GT.1.0D-6) THEN
C          tanh_LAMDA = sinh_LAMDA/cosh_LAMDA
C        ELSE
C          tanh_LAMDA = 1.0D8
C        ENDIF
C        IF(DABS(sinh_LAMDA).GT.1.0D-6) THEN
C          coth_LAMDA = cosh_LAMDA/sinh_LAMDA
C        ELSE
C          coth_LAMDA = 1.0D8
C        ENDIF
C
C!  Mu
C        nb=NBJ(2)
C        RMU        = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,2))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''    RMU ='',E12.3)') RMU
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        dMU_dXi1   = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,
C     '    XE(1,2))
C        dMU_dXi2   = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,XI,
C     '    XE(1,2))
C        sin_MU     = DSIN(RMU)
C        cos_MU     = DCOS(RMU)
C        IF(DABS(cos_MU).GT.1.0D-6) THEN
C          tan_MU = sin_MU/cos_MU
C        ELSE
C          tan_MU = 1.0D8
C        ENDIF
C        IF(DABS(sin_MU).GT.1.0D-6) THEN
C          cot_MU     = cos_MU/sin_MU
C        ELSE
C          cot_MU = 1.0D8
C        ENDIF
C
C!  Theta
C        nb=NBJ(3)
C        THETA = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,3))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''  THETA ='',E12.3)') THETA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        dTHETA_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,
C     '    XE(1,3))
C
C!  Alfa (fibre angle)
C        nj=NJ_LOC(NJL_FIBR,1)
C        nb=NBJ(nj)
C        ALFA = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,nj))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''   ALFA ='',E12.3)') ALFA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_ALFA   = DSIN(ALFA)
C        cos_ALFA   = DCOS(ALFA)
C
C!  Beta (imbrication angle)
C        nj=NJ_LOC(NJL_FIBR,2)
C        nb=NBJ(nj)
C        IF(nb.GT.0) THEN
C          BETA=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '      XE(1,nj))
C        ELSE
C          BETA=0.0D0
C        ENDIF
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''   BETA ='',E12.3)') BETA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_BETA   = DSIN(BETA)
C        cos_BETA   = DCOS(BETA)
C
C!  Gama (sheet angle)
C        nj=NJ_LOC(NJL_FIBR,3)
C        nb=NBJ(nj)
C        GAMA = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,nj))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''   GAMA ='',E12.3)') GAMA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_GAMA   = DSIN(GAMA)
C        cos_GAMA   = DCOS(GAMA)
C
C!  Phi  (eqtn 8 in "Laminar structure of the heart II")
C        ARG1=coth_LAMDA*dLAMDA_dXi1+cot_MU*dMU_dXi1
C        ARG2=dTHETA_dXi1
C        PHI =DATAN2(ARG1,ARG2)
C        IF(PHI.GT.PI/2.0D0) THEN
C          PHI=PHI-PI
C        ELSE IF(PHI.LT.-PI/2.0D0) THEN
C          PHI=PHI+PI
C        ENDIF
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''    PHI ='',E12.3)') PHI
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_PHI   = DSIN(PHI)
C        cos_PHI   = DCOS(PHI)
C
C!  Delta (eqtn 12 in "Laminar structure of the heart II")
C        ARG1=(tan_MU*dLAMDA_dXi2+tanh_LAMDA*dMU_dXi2)*cos_PHI
C        ARG2=tanh_LAMDA*dLAMDA_dXi2-tan_MU*dMU_dXi2
C        DELTA=DATAN2(ARG1,ARG2)
C        IF(DELTA.GT.PI/2.0D0) THEN
C          DELTA=DELTA-PI
C        ELSE IF(DELTA.LT.-PI/2.0D0) THEN
C          DELTA=DELTA+PI
C        ENDIF
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''  DELTA ='',E12.3)') DELTA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_DELTA  = DSIN(DELTA)
C        cos_DELTA  = DCOS(DELTA)
C
C!  A vector (eqtn 3 in "Laminar structure of the heart II")
C        IF(DABS(BETA).LT.1.0D-5) THEN !no imbrication angle
C          A_VECTOR(1)=-sin_ALFA*cos_DELTA
C          A_VECTOR(2)=-cos_ALFA*sin_PHI-sin_ALFA*sin_DELTA*cos_PHI
C          A_VECTOR(3)=-cos_ALFA*cos_PHI+sin_ALFA*sin_DELTA*sin_PHI
C        ELSE !include imbrication angle
C          D1=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*cos_DELTA
C     '       + sin_BETA*cos_GAMA*sin_DELTA
C          D2=  cos_ALFA*cos_BETA+sin_ALFA*sin_BETA*sin_GAMA
C          D3=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*sin_DELTA
C     '       - sin_BETA*cos_GAMA*cos_DELTA
C          A_VECTOR(1)= D1
C          A_VECTOR(2)=-D2*sin_PHI+D3*cos_PHI
C          A_VECTOR(3)=-D2*cos_PHI-D3*sin_PHI
C        ENDIF
C        CALL ROTATION(1,THETA,A_VECTOR,ERROR,*9999)
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' a_vector:'',3E12.3)')
C     '      (A_VECTOR(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C!  B vector (eqtn 4 in "Laminar structure of the heart II")
C        ARG1=cos_ALFA*sin_GAMA*sin_DELTA-cos_GAMA*cos_DELTA
C        B_VECTOR(1)=-cos_ALFA*sin_GAMA*cos_DELTA - cos_GAMA*sin_DELTA
C        B_VECTOR(2)= sin_ALFA*sin_GAMA*sin_PHI - ARG1*cos_PHI
C        B_VECTOR(3)= sin_ALFA*sin_GAMA*cos_PHI + ARG1*sin_PHI
C        CALL ROTATION(1,THETA,B_VECTOR,ERROR,*9999)
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' b_vector:'',3E12.3)')
C     '      (B_VECTOR(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C!  C vector (eqtn 5 in "Laminar structure of the heart II")
C        CALL CROSS(A_VECTOR,B_VECTOR,C_VECTOR)
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' c_vector:'',3E12.3)')
C     '      (C_VECTOR(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C      ENDIF

      CALL EXITS('MAT_VEC_XI')
      RETURN
 9999 CALL ERRORS('MAT_VEC_XI',ERROR)
      CALL EXITS('MAT_VEC_XI')
      RETURN 1
      END


      SUBROUTINE MELGE(LGE,NBH,nc,ne,NHE,NHST,NKH,NPNE,nr,NVHE,nx,
     '  NYNE,NYNP,ERROR,*)

C#### Subroutine: MELGE
C###  Description:
C###    MELGE calculates the row numbers (LGE(*,1)) and column numbers
C###    (LGE(*,2)) in the global matrix nc for element variables nhs
C###    in region nr. It also returns the total number of element
C###    variables NHST(nrc).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER LGE(NHM*NSM,NRCM),NBH(NHM,NCM),nc,ne,NHE,NHST(2),
     '  NKH(NHM,NPM,NCM),NPNE(NNM,NBFM),nr,NVHE(NNM,NBFM,NHM),
     '  nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER na,nb,nh,nhs,nhx,nk,nn,np,nrc,nv

      CALL ENTERS('MELGE',*9999)

      DO nrc=1,2
        NHST(nrc)=0
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx)
C !!!     Use the LHS (nc=1) basis to determine the # of equations
          IF(nrc.EQ.1) THEN
            nb=NBH(nh,1)
          ELSE
            nb=NBH(nh,nc)
          ENDIF
          DO nn=1,NNT(nb)            !nodal variables
            np=NPNE(nn,nb)
            nv=NVHE(nn,nb,nh)
C CPB 12/9/95 Using NKT instead of NKH
C            IF(nrc.EQ.1) THEN
C              NK_TOT=MAX(NKH(nh,np,1)-KTYP93(1,nr),NKH(nh,np,2)-
C     '          KTYP93(2,nr),1)
C            ELSE
C              NK_TOT=MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
C            ENDIF
C            DO nk=1,NK_TOT
            DO nk=1,NKT(nn,nb)
              NHST(nrc)=NHST(nrc)+1
              LGE(NHST(nrc),nrc)=NYNP(nk,nv,nh,np,nrc,nc,nr)
            ENDDO !nk
          ENDDO !nn
          DO na=1,NAT(nb)            !auxillary variables
            NHST(nrc)=NHST(nrc)+1
            LGE(NHST(nrc),nrc)=NYNE(na,nh,nrc,nc,ne)
          ENDDO !na
        ENDDO !nh
      ENDDO !nrc

      IF(DOP) THEN
        DO nrc=1,2
          WRITE(OP_STRING,'('' NHST('',I1,'')='',I4,'' LGE(nhs,'',I1,'
     '      //'''):'',/(20I4))') nrc,NHST(nrc),nrc,(LGE(nhs,nrc),nhs=1,
     '      NHST(nrc))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('MELGE')
      RETURN
 9999 CALL ERRORS('MELGE',ERROR)
      CALL EXITS('MELGE')
      RETURN 1
      END


      SUBROUTINE TOFFEL(nb,NJE,nr,CHTOFF,DBM,GU,XG,X3G,CURVE,ERROR,*)

C#### Subroutine: TOFFEL
C###  Description:
C###    TOFFEL evaluates Christoffel symbol of 2nd kind CHTOFF(ia,ib,ic)
C###    and, if CURVE is .true., the curvature tensor DBM(ia,ib,ic)
C###    in Xi-coordinate system at current Gauss point.
C**** X3G(nj,nd) are 3rd derivs of Xj (see subroutine X3XG)
C**** G(nj,1) is tangent to surface in Xi(1) direction
C**** G(nj,2) is    "     "    "     " Xi(2)     "
C**** G(nj,3) is    "     "    "      "Xi(3)     "
C**** DG(nj,ib,ic) are cpts of base vector G(nj,ib) derivs wrt Xi(ic)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
!     Parameter List
      INTEGER nb,NJE,nr
      REAL*8 CHTOFF(3,3,3),DBM(3,3,3),GU(3,3),XG(NJM,NUM),X3G(4,3)
      CHARACTER ERROR*(*)
      LOGICAL CURVE
!     Local Variables
      INTEGER i,ia,ib,ic,id,ig,ij,il,j,jj,k,kj,lj,m,
     '  ND3(2,2,2),NITB,ni,nj,nk,nl,NU1(0:3),NU2(3,3)
      REAL*8 COSHL,COSM,COSP,COST,DDG(3,3,3,3),DG(3,3,3),dXrc_dXref,DZX,
     '  DZXX,DZXXX(3,3,3,3),G(3,3),X(3),
     '  LAMDA,MU,PHI,RAD,SINHL,SINM,SINP,SINT,SUM,SUM1,SUM2,THETA

      DATA NU1/1,2,4,7/
      DATA NU2/3,6,10,6,5,9,10,9,8/
      DATA ND3/1,2,2,3,2,3,3,4/

      CALL ENTERS('TOFFEL',*9999)

C     Initialise all real arrays to zero.
      DO ni=1,3
        DO nj=1,3
          G(nj,ni)=0.0d0
          DO nk=1,3
            DG(nk,nj,ni)=0.0d0
            DBM(ni,nj,nk)=0.0d0
            DO nl=1,3
              DDG(nl,nk,nj,ni)=0.0d0
              DZXXX(nl,nk,nj,ni)=0.0d0
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      NITB=NIT(nb)
      IF(ITYP10(nr).EQ.1) THEN !rectangular cartesian coords
        DO nj=1,NJE
          DO j=1,NITB
            G(nj,j)=XG(nj,NU1(j))
          ENDDO
          DO j=1,NITB
            DO k=1,NITB
              DG(nj,j,k)=XG(nj,NU2(j,k))
            ENDDO
          ENDDO
        ENDDO
        DO i=1,NITB
          DO j=1,NITB
            DO k=1,NITB
              SUM=0.0D0
              DO nj=1,NJE
                DO m=1,NITB
                  SUM=SUM+DG(nj,j,k)*G(nj,m)*GU(i,m)
                ENDDO
              ENDDO
              IF(DABS(SUM).GT.1.0D-8) THEN
                CHTOFF(i,j,k)=SUM
              ELSE
                CHTOFF(i,j,k)=0.0D0
              ENDIF
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' CHTOFF('',I1,'','',I1,'',k'''
     '          //''') = '',3E13.5)') i,j,(CHTOFF(i,j,k),k=1,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO

      ELSE IF(ITYP10(nr).GT.1) THEN !curvilinear coordinates

C old MPN 24-Apr-96: DZX and DZXX functions in FE01 now used
C        IF(ITYP10(nr).EQ.2) THEN !cylindrical polar coords
C          RAD=XG(1,1)
C          THETA=XG(2,1)
C          SINT=DSIN(THETA)
C          COST=DCOS(THETA)
C          DZX(1,1)=COST
C          DZX(1,2)=-RAD*SINT
C          DZX(2,1)=SINT
C          DZX(2,2)=RAD*COST
C          DZX(3,3)=1.0D0
C          DZXX(1,2,1)=-SINT
C          DZXX(2,2,1)=COST
C          DZXX(1,1,2)=-SINT
C          DZXX(1,2,2)=-RAD*COST
C          DZXX(2,1,2)=COST
C          DZXX(2,2,2)=-RAD*SINT
C
C        ELSE IF(ITYP10(nr).EQ.3) THEN !spherical polar coords
C          RAD  =XG(1,1)
C          THETA=XG(2,1)
C          PHI  =XG(3,1)
C          SINT=DSIN(THETA)
C          COST=DCOS(THETA)
C          SINP=DSIN(PHI)
C          COSP=DCOS(PHI)
C          DZX(1,1)=COST*COSP
C          DZX(1,2)=-RAD*SINT*COSP
C          DZX(1,3)=-RAD*COST*SINP
C          DZX(2,1)=SINT*COSP
C          DZX(2,2)=RAD*COST*COSP
C          DZX(2,3)=-RAD*SINT*SINP
C          DZX(3,1)=SINP
C          DZX(3,3)=RAD*COSP
C          DZXX(1,2,1)=-SINT*COSP
C          DZXX(1,3,1)=-COST*SINP
C          DZXX(2,2,1)=COST*COSP
C          DZXX(2,3,1)=-SINT*SINP
C          DZXX(3,3,1)=COSP
C          DZXX(1,1,2)=-SINT*COSP
C          DZXX(1,2,2)=-RAD*COST*COSP
C          DZXX(1,3,2)=RAD*SINT*SINP
C          DZXX(2,1,2)=COST*COSP
C          DZXX(2,2,2)=-RAD*SINT*COSP
C          DZXX(2,3,2)=-RAD*COST*SINP
C          DZXX(1,1,3)=-COST*SINP
C          DZXX(1,2,3)=RAD*SINT*SINP
C          DZXX(1,3,3)=-RAD*COST*COSP
C          DZXX(2,1,3)=-SINT*SINP
C          DZXX(2,2,3)=-RAD*COST*SINP
C          DZXX(2,3,3)=-RAD*SINT*COSP
C          DZXX(3,1,3)=COSP
C          DZXX(3,3,3)=-RAD*SINP
C
C        ELSE IF(ITYP10(nr).EQ.4) THEN !prolate spheroidal coords
C          LAMDA=XG(1,1)
C          mu=XG(2,1)
C          THETA=XG(3,1)
C          SINHL=DSINH(LAMDA)
C          COSHL=DCOSH(LAMDA)
C          SINM=DSIN(MU)
C          COSM=DCOS(MU)
C          SINT=DSIN(THETA)
C          COST=DCOS(THETA)
C          DZX(1,1)=FOCUS*SINHL*COSM
C          DZX(1,2)=-FOCUS*COSHL*SINM
C          DZX(1,3)=0.0D0
C          DZX(2,1)=FOCUS*COSHL*SINM*COST
C          DZX(2,2)=FOCUS*SINHL*COSM*COST
C          DZX(2,3)=-FOCUS*SINHL*SINM*SINT
C          DZX(3,1)=FOCUS*COSHL*SINM*SINT
C          DZX(3,2)=FOCUS*SINHL*COSM*SINT
C          DZX(3,3)=FOCUS*SINHL*SINM*COST
C          DZXX(1,1,1)=FOCUS*COSHL*COSM
C          DZXX(1,2,1)=-FOCUS*SINHL*SINM
C          DZXX(2,1,1)=FOCUS*SINHL*SINM*COST
C          DZXX(2,2,1)=FOCUS*COSHL*COSM*COST
C          DZXX(2,3,1)=-FOCUS*COSHL*SINM*SINT
C          DZXX(3,1,1)=FOCUS*SINHL*SINM*SINT
C          DZXX(3,2,1)=FOCUS*COSHL*COSM*SINT
C          DZXX(3,3,1)=FOCUS*COSHL*SINM*COST
C          DZXX(1,1,2)=-FOCUS*SINHL*SINM
C          DZXX(1,2,2)=-FOCUS*COSHL*COSM
C          DZXX(2,1,2)=FOCUS*COSHL*COSM*COST
C          DZXX(2,2,2)=-FOCUS*SINHL*SINM*COST
C          DZXX(2,3,2)=-FOCUS*SINHL*COSM*SINT
C          DZXX(3,1,2)=FOCUS*COSHL*COSM*SINT
C          DZXX(3,2,2)=-FOCUS*SINHL*SINM*SINT
C          DZXX(3,3,2)=FOCUS*SINHL*COSM*COST
C          DZXX(2,1,3)=-FOCUS*COSHL*SINM*SINT
C          DZXX(2,2,3)=-FOCUS*SINHL*COSM*SINT
C          DZXX(2,3,3)=-FOCUS*SINHL*SINM*COST
C          DZXX(3,1,3)=FOCUS*COSHL*SINM*COST
C          DZXX(3,2,3)=FOCUS*SINHL*COSM*COST
C          DZXX(3,3,3)=-FOCUS*SINHL*SINM*SINT
C
C        ENDIF !ityp10

        DO nj=1,NJE
          X(nj)=XG(nj,1)
        ENDDO!nj

        DO nj=1,NJE
          DO j=1,NITB
            SUM=0.0D0
            DO ij=1,NJE
              dXrc_dXref=DZX(ITYP10(nr),nj,ij,X)
              SUM=SUM+dXrc_dXref*XG(ij,NU1(j))
            ENDDO
            G(nj,j)=SUM
            DO k=1,NITB
              SUM=0.0D0
              DO ij=1,NJE
                dXrc_dXref=DZX(ITYP10(nr),nj,ij,X)
                SUM=SUM+dXrc_dXref*XG(ij,NU2(j,k))
                DO lj=1,NJE
                  SUM=SUM+DZXX(ITYP10(nr),nj,ij,lj,X)*
     '              XG(ij,NU1(j))*XG(lj,NU1(k))
                ENDDO
              ENDDO
              IF(DABS(SUM).GT.1.0D-8) THEN
                DG(nj,j,k)=SUM
              ELSE
                DG(nj,j,k)=0.0D0
              ENDIF
            ENDDO !k
          ENDDO !j
        ENDDO !nj

        DO i=1,NITB
          DO j=1,NITB
            DO k=1,NITB
              SUM=0.0D0
              DO m=1,NITB
                SUM1=0.0D0
                DO nj=1,NJE
                  SUM1=SUM1+DG(nj,j,k)*G(nj,m)
                ENDDO
                SUM=SUM+SUM1*GU(i,m)
              ENDDO !m
              IF(DABS(SUM).GT.1.0D-8) THEN
                CHTOFF(i,j,k)=SUM
              ELSE
                CHTOFF(i,j,k)=0.0D0
              ENDIF
            ENDDO !k
          ENDDO !j
        ENDDO !i
      ENDIF
C
C  If curvature tensors and curvature change tensors are required,
C  find with CURVE=.TRUE.
C
      IF(CURVE) THEN
        IF(NITB.EQ.2) THEN !surface geometry
          G(1,3)=G(2,1)*G(3,2)-G(3,1)*G(2,2)
          G(2,3)=G(3,1)*G(1,2)-G(1,1)*G(3,2)
          G(3,3)=G(1,1)*G(2,2)-G(2,1)*G(1,2)
          SUM=0.0D0
          DO nj=1,NJE
            SUM=SUM+G(nj,3)**2
          ENDDO
          DO nj=1,NJE
            G(nj,3)=G(nj,3)/DSQRT(SUM) !sqrt(sum)=rg
          ENDDO
          DO j=1,NITB
            DO k=1,NITB
              SUM=0.0D0
              DO nj=1,NJE
                SUM=SUM+DG(nj,j,k)*G(nj,3)
              ENDDO
              IF(DABS(SUM).GT.1.0D-8) THEN
                CHTOFF(3,j,k)=SUM
              ELSE
                CHTOFF(3,j,k)=0.0D0
              ENDIF
            ENDDO !k
          ENDDO !j

          IF(ITYP10(nr).EQ.1) THEN !rectangular cartesian coords
            DO nj=1,NITB
              DO jj=1,NITB
                DO lj=1,NITB
                  DO kj=1,NITB
                    DDG(nj,jj,lj,kj)=X3G(ND3(jj,lj,kj),nj)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ELSE IF(ITYP10(nr).GT.1) THEN !curvilinear coordinates

C!!!        DZXXX should be put into a function like DZXX in FE01

            IF(ITYP10(nr).EQ.2) THEN !cylindrical polar coords
              RAD=XG(1,1)
              THETA=XG(2,1)
              SINT=DSIN(THETA)
              COST=DCOS(THETA)
              DZXXX(1,2,1,2)=-COST
              DZXXX(2,2,1,2)=-SINT
              DZXXX(1,2,2,1)=-COST
              DZXXX(2,2,2,1)=-SINT
              DZXXX(1,1,2,2)=-COST
              DZXXX(1,2,2,2)=RAD*SINT
              DZXXX(2,1,2,2)=-SINT
              DZXXX(2,2,2,2)=-RAD*COST

            ELSE IF(ITYP10(nr).EQ.3) THEN !spherical polar coords
              RAD  =XG(1,1)
              THETA=XG(2,1)
              PHI  =XG(3,1)
              SINT=DSIN(THETA)
              COST=DCOS(THETA)
              SINP=DSIN(PHI)
              COSP=DCOS(PHI)
              DZXXX(1,2,1,2)=-COST*COSP
              DZXXX(1,3,1,2)=SINT*SINP
              DZXXX(2,2,1,2)=-SINT*COSP
              DZXXX(2,3,1,2)=-COST*SINP
              DZXXX(1,2,1,3)=SINT*SINP
              DZXXX(1,3,1,3)=-COST*COSP
              DZXXX(2,2,1,3)=-COST*SINP
              DZXXX(2,3,1,3)=-SINT*COSP
              DZXXX(3,3,1,3)=-SINT
              DZXXX(1,2,2,1)=-COST*COSP
              DZXXX(1,3,2,1)=SINT*SINP
              DZXXX(2,2,2,1)=-SINT*COSP
              DZXXX(2,3,2,1)=-COST*SINP
              DZXXX(1,1,2,2)=-COST*COSP
              DZXXX(1,2,2,2)=RAD*SINT*COSP
              DZXXX(1,3,2,2)=RAD*COST*SINP
              DZXXX(2,1,2,2)=-SINT*COSP
              DZXXX(2,2,2,2)=-RAD*COST*COSP
              DZXXX(2,3,2,2)=RAD*SINT*SINP
              DZXXX(1,1,2,3)=SINT*SINP
              DZXXX(1,2,2,3)=RAD*COST*SINP
              DZXXX(1,3,2,3)=RAD*SINT*COSP
              DZXXX(2,1,2,3)=-COST*SINP
              DZXXX(2,2,2,3)=RAD*SINT*SINP
              DZXXX(2,3,2,3)=-RAD*COST*COSP
              DZXXX(1,2,3,1)=SINT*SINP
              DZXXX(1,3,3,1)=-COST*COSP
              DZXXX(2,2,3,1)=-COST*SINP
              DZXXX(2,3,3,1)=-SINT*COSP
              DZXXX(3,3,3,1)=-SINP
              DZXXX(1,1,3,2)=SINT*SINP
              DZXXX(1,2,3,2)=RAD*COST*SINP
              DZXXX(1,3,3,2)=RAD*SINT*COSP
              DZXXX(2,1,3,2)=-COST*SINP
              DZXXX(2,2,3,2)=RAD*SINT*SINP
              DZXXX(2,3,3,2)=RAD*COST*COSP
              DZXXX(1,1,3,3)=-COST*COSP
              DZXXX(1,2,3,3)=RAD*SINT*COSP
              DZXXX(1,3,3,3)=RAD*COST*SINP
              DZXXX(2,1,3,3)=-SINT*COSP
              DZXXX(2,2,3,3)=-RAD*COST*COSP
              DZXXX(2,3,3,3)=RAD*SINT*SINP
              DZXXX(3,1,3,3)=-SINP
              DZXXX(3,3,3,3)=-RAD*COSP

            ELSE IF(ITYP10(nr).EQ.4) THEN !prolate spheroidal coords
              LAMDA=XG(1,1)
              MU=XG(2,1)
              THETA=XG(3,1)
              SINHL=DSINH(LAMDA)
              COSHL=DCOSH(LAMDA)
              SINM=DSIN(MU)
              COSM=DCOS(MU)
              SINT=DSIN(THETA)
              COST=DCOS(THETA)
              DZXXX(1,1,1,1)=FOCUS*SINHL*COSM
              DZXXX(1,2,1,1)=-FOCUS*COSHL*SINM
              DZXXX(2,1,1,1)=FOCUS*COSHL*SINM*COST
              DZXXX(2,2,1,1)=FOCUS*SINHL*COSM*COST
              DZXXX(2,3,1,1)=-FOCUS*SINHL*SINM*SINT
              DZXXX(3,1,1,1)=FOCUS*COSHL*SINM*SINT
              DZXXX(3,2,1,1)=FOCUS*SINHL*COSM*SINT
              DZXXX(3,3,1,1)=FOCUS*SINHL*SINM*COST
              DZXXX(1,1,1,2)=-FOCUS*COSHL*SINM
              DZXXX(1,2,1,2)=-FOCUS*SINHL*COSM
              DZXXX(2,1,1,2)=FOCUS*SINHL*COSM*COST
              DZXXX(2,2,1,2)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,3,1,2)=-FOCUS*COSHL*COSM*SINT
              DZXXX(3,1,1,2)=FOCUS*SINHL*COSM*SINT
              DZXXX(3,2,1,2)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,3,1,2)=FOCUS*COSHL*COSM*COST
              DZXXX(2,1,1,3)=-FOCUS*SINHL*SINM*SINT
              DZXXX(2,2,1,3)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,3,1,3)=-FOCUS*COSHL*SINM*COST
              DZXXX(3,1,1,3)=FOCUS*SINHL*SINM*COST
              DZXXX(3,2,1,3)=FOCUS*COSHL*COSM*COST
              DZXXX(3,3,1,3)=-FOCUS*COSHL*SINM*SINT
              DZXXX(1,1,2,1)=-FOCUS*COSHL*SINM
              DZXXX(1,2,2,1)=-FOCUS*SINHL*COSHL
              DZXXX(2,1,2,1)=FOCUS*SINHL*COSM*COST
              DZXXX(2,2,2,1)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,3,2,1)=-FOCUS*COSHL*COSM*SINT
              DZXXX(3,1,2,1)=FOCUS*SINHL*COSM*SINT
              DZXXX(3,2,2,1)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,3,2,1)=FOCUS*COSHL*COSM*COST
              DZXXX(1,1,2,2)=-FOCUS*SINHL*COSM
              DZXXX(1,2,2,2)=FOCUS*COSHL*SINM
              DZXXX(2,1,2,2)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,2,2,2)=-FOCUS*SINHL*COSM*COST
              DZXXX(2,3,2,2)=FOCUS*SINHL*SINM*SINT
              DZXXX(3,1,2,2)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,2,2,2)=-FOCUS*SINHL*COSM*SINT
              DZXXX(3,3,2,2)=-FOCUS*SINHL*SINM*COST
              DZXXX(2,1,2,3)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,2,2,3)=FOCUS*SINHL*SINM*SINT
              DZXXX(2,3,2,3)=-FOCUS*SINHL*COSM*COST
              DZXXX(3,1,2,3)=FOCUS*COSHL*COSM*COST
              DZXXX(3,2,2,3)=-FOCUS*SINHL*SINM*COST
              DZXXX(3,3,2,3)=-FOCUS*SINHL*COSM*SINT
              DZXXX(2,1,3,1)=-FOCUS*SINHL*SINM*SINT
              DZXXX(2,2,3,1)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,3,3,1)=-FOCUS*COSHL*SINM*COST
              DZXXX(3,1,3,1)=FOCUS*SINHL*SINM*COST
              DZXXX(3,2,3,1)=FOCUS*COSHL*COSM*COST
              DZXXX(3,3,3,1)=-FOCUS*COSHL*SINM*SINT
              DZXXX(2,1,3,2)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,2,3,2)=FOCUS*SINHL*SINM*SINT
              DZXXX(2,3,3,2)=-FOCUS*SINHL*COSM*COST
              DZXXX(3,1,3,2)=FOCUS*COSHL*COSM*COST
              DZXXX(3,2,3,2)=-FOCUS*SINHL*SINM*COST
              DZXXX(3,3,3,2)=-FOCUS*SINHL*COSM*SINT
              DZXXX(2,1,3,3)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,2,3,3)=-FOCUS*SINHL*COSM*COST
              DZXXX(2,3,3,3)=FOCUS*SINHL*SINM*SINT
              DZXXX(3,1,3,3)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,2,3,3)=-FOCUS*SINHL*COSM*SINT
              DZXXX(3,3,3,3)=-FOCUS*SINHL*SINM*COST
            ENDIF !ityp10(nr)

            DO nj=1,NJE
              DO ia=1,2
                DO ib=1,2
                  DO ic=1,2
                    SUM=0.0D0
                    DO ij=1,NJE
                      dXrc_dXref=DZX(ITYP10(nr),nj,ij,X)
                      SUM=SUM+dXrc_dXref*X3G(ND3(ia,ib,ic),ij)
                      DO kj=1,NJE
                        SUM=SUM+DZXX(ITYP10(nr),nj,ij,kj,X)*
     '                    (XG(ij,NU2(ia,ic))*
     '                    XG(kj,NU1(ib))+XG(ij,NU1(ia))*
     '                    XG(kj,NU2(ib,ic))+XG(kj,NU1(ic))*
     '                    XG(ij,NU2(ia,ib)))
                        DO jj=1,NJE
                          SUM=SUM+DZXXX(nj,ij,jj,kj)*(XG(kj,NU1(ic))
     '                      *XG(ij,NU1(ia))*XG(jj,NU1(ib)))
                        ENDDO !jj
                      ENDDO !kj
                    ENDDO !ij
                    IF(DABS(SUM).LT.1.0D-8) THEN
                      DDG(nj,ia,ib,ic)=0.0D0
                    ELSE
                      DDG(nj,ia,ib,ic)=SUM
                    ENDIF
                  ENDDO !ic
                ENDDO !ib
              ENDDO !ia
            ENDDO !nj
          ENDIF !rc/curvilinear coords

          DO ia=1,NITB
            DO ib=1,NITB
              DO ig=1,NITB
                SUM=0.0d0
                SUM2=0.0d0
                DO il=1,2
                  SUM1=0.0d0
                  DO nj=1,NJE
                    SUM1=SUM1+G(nj,3)*DDG(nj,il,ib,ig)
                  ENDDO
                  SUM2=SUM2+GU(ia,il)*SUM1
                  DO id=1,2
                    SUM=SUM-GU(ia,il)*(CHTOFF(id,il,ig)*CHTOFF(3,id,ib)
     '                +CHTOFF(id,ig,ib)*CHTOFF(3,id,il)
     '                +CHTOFF(id,ib,il)*CHTOFF(3,id,ig))
                  ENDDO
                ENDDO
                IF(DABS(SUM+SUM2).GT.1.0d-8) THEN
                  DBM(ia,ib,ig)=SUM+SUM2
                ELSE
                  DBM(ia,ib,ig)=0.0d0
                ENDIF
              ENDDO !ig
            ENDDO !ib
          ENDDO !ia
        ENDIF !nit=2
      ENDIF !curve

      CALL EXITS('TOFFEL')
      RETURN
 9999 CALL ERRORS('TOFFEL',ERROR)
      CALL EXITS('TOFFEL')
      RETURN 1
      END


      SUBROUTINE XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*)

C#### Subroutine: XEXG
C###  Description:
C###    XEXG evaluates Gauss point array XG from element node array
C###    XE at current Gauss point ng. XG is corrected for JTYP10=2
C###    (isochoric interpolation).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER NBJ(NJM),ng,nr
      REAL*8 PG(NSM,NUM,NGM,NBM),VE(NSM,NKM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,na,nb,ni,nj,njj1,njj2,nk,ns,nu,NU2(3,3)
      REAL*8 COSHX,CSS,D,DDOT,DES,RAD,SINHX,SS,SUM,SUMM,THETA

      DATA NU2/3,6,10,6,5,9,10,9,8/

      CALL ENTERS('XEXG',*9999)
      DO njj1=1,3,2   !Loop over geometry and field
        DO njj2=1,NJ_LOC(njj1,0)
          nj=NJ_LOC(njj1,njj2)
          nb=NBJ(nj)
          DO nu=1,NUT(nb)
            SUM=0.0D0
            IF(nb.GT.0.AND.NNT(nb).GT.0) THEN
c PJH 26Aug95 use BLAS routine
!             DO ns=1,NST(nb)+NAT(nb)
!               SUM=SUM+PG(ns,nu,ng,nb)*XE(ns,nj)
!             ENDDO
              SUM=DDOT(NST(nb)+NAT(nb),PG(1,nu,ng,nb),1,XE(1,nj),1)
            ELSE IF(NAT(nb).GT.0) THEN
              DO nk=1,NKT(0,nb)
                SUMM=0.0D0
                DO na=1,NAT(nb)
                  SUMM=SUMM+VE(na,nk)*XE(na,nj)
                ENDDO
                SUM=SUM+PG(nk,nu,ng,nb)*SUMM
              ENDDO !nk
            ENDIF
            XG(nj,nu)=SUM
          ENDDO !nu
        ENDDO !njj2
      ENDDO !njj1

C     Fibre angle interpolation
      IF(JTYP9.GT.0) THEN !fibre angle defined
        DO njj1=1,NJ_LOC(NJL_FIBR,0)
          nj=NJ_LOC(NJL_FIBR,njj1)
          nb=NBJ(nj)
          IF(nb.GT.0) THEN
            DO nu=1,NUT(nb)
              SUM=0.0D0
              IF(nb.GT.0.AND.NNT(nb).GT.0) THEN
                DO ns=1,NST(nb)+NAT(nb)
                  SUM=SUM+PG(ns,nu,ng,nb)*XE(ns,nj)
                ENDDO
              ENDIF
              XG(nj,nu)=SUM
            ENDDO !nu
          ENDIF !nb
        ENDDO !njj1
      ENDIF !jtyp9

C     Isochoric interpolation
      IF(JTYP10.GE.2) THEN
        nb=NBJ(1)

        IF(ITYP10(nr).EQ.2) THEN      !cyl polar
          RAD=DSQRT(XG(1,1))
          XG(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
          RAD=XG(1,1)**(1.0D0/3.0D0)
          XG(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
          IF(JTYP10.EQ.2) THEN
            SS=XG(1,1)/(FOCUS*FOCUS)
            SINHX=DSQRT(SS)
            COSHX=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=XG(1,1)/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0D0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHX=D+1.0D0/(3.0D0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHX=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
            SINHX=DSQRT(COSHX*COSHX-1.D0)
          ENDIF
          XG(1,1)=DLOG(COSHX+SINHX)
        ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
        ENDIF

        DO ni=1,NIT(nb)
          nu=1+ni*(ni+1)/2
          IF(ITYP10(nr).EQ.2) THEN      !cyl polar
            XG(1,nu)=XG(1,nu)/(2.0D0*RAD)
          ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
            XG(1,nu)=XG(1,nu)/(3.0D0*RAD*RAD)
          ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
            IF(JTYP10.EQ.2) THEN
              XG(1,nu)=XG(1,nu)/(2.0D0*FOCUS*FOCUS*SINHX*COSHX)
            ELSE IF(JTYP10.EQ.3) THEN
              XG(1,nu)=XG(1,nu)/((3.0D0*COSHX*COSHX-1.0D0)*SINHX)
     '                               /FOCUS**3
            ENDIF
          ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
          ENDIF
        ENDDO !ni

        DO ni=1,NIT(nb)
          DO mi=1,ni
            nu=NU2(ni,mi)
            IF(ITYP10(nr).EQ.2) THEN      !cyl polar
              XG(1,nu)=-XG(1,nu)/(4.0D0*RAD**3)
            ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
              XG(1,nu)=-XG(1,nu)*2.0D0/(9.0D0*RAD**5)
            ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
              IF(JTYP10.EQ.2) THEN
                XG(1,nu)=-XG(1,nu)*(COSHX*COSHX+SINHX*SINHX)
     '            /(4.0d0*FOCUS**4*SINHX**3*COSHX**3)
              ELSE IF(JTYP10.EQ.3) THEN
                XG(1,nu)=-XG(1,nu)*(2.0D0+9.0D0*SINHX*SINHX)*COSHX/
     '            (FOCUS**6*(3.0D0*COSHX*COSHX-1)**3*SINHX**3)
              ENDIF
            ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
            ENDIF
          ENDDO !mi
        ENDDO !ni
      ENDIF !jtyp10

      IF(DOP) THEN
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0)
            nj=NJ_LOC(njj1,njj2)
            WRITE(OP_STRING,'('' XG(nj='',I1,'',nu=1..): '',6D12.4,'
     '        //'/(18X,6D12.4))') nj,(XG(nj,nu),nu=1,NUT(NBJ(nj)))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !njj2
        ENDDO !njj1
      ENDIF !dop

      CALL EXITS('XEXG')
      RETURN
 9999 CALL ERRORS('XEXG',ERROR)
      CALL EXITS('XEXG')
      RETURN 1
      END


      SUBROUTINE XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XW,XI,ERROR,*)

C#### Subroutine: XEXW
C###  Description:
C###    XEXW interpolates XE(nj,1) and first derivatives XE(nj,NU1(ni))
C###    at XI for Lagrange/Hermite tensor product basis functions
C###    and returns XW(nj,NU1(ni)). XW is corrected for JTYP10=2
C###    (isochoric interpolation).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBJ(NJM),nr
      REAL*8 XI(3),XE(NSM,NJM),XW(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nj,njj1,njj2,ni,NITB,nu,NU1(0:3)
      REAL*8 COSHZ,CSS,D,DES,PFXI,RAD,SINHZ,SS,THETA

      DATA NU1/1,2,4,7/

      CALL ENTERS('XEXW',*9999)
      NITB=NIT(NBJ(1))
      DO ni=0,NITB
        DO njj1=1,3 !geometry/fibres/field
          DO njj2=1,NJ_LOC(njj1,0)
            nj=NJ_LOC(njj1,njj2)
            nb=NBJ(nj)
!new        MPN 19/7/93 array overrun for NNT in PFXI if nb=0
            IF(nb.GT.0) THEN
              XW(nj,NU1(ni))=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          NAN(1,1,nb),nb,NU1(ni),XE(1,nj),XI)
            ELSE
              XW(nj,NU1(ni))=0.0D0 !should never get into here!!!
            ENDIF
          ENDDO !njj2
        ENDDO !njj1
      ENDDO !ni

      IF(JTYP10.GE.2) THEN
        nb=NBJ(1)
        IF(ITYP10(nr).EQ.2) THEN
          RAD=DSQRT(XW(1,1))
          XW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.3) THEN
          RAD=XW(1,1)**(1.0D0/3.0D0)
          XW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.4) THEN
          IF(JTYP10.EQ.2) THEN
            SS=XW(1,1)/(FOCUS*FOCUS)
            SINHZ=DSQRT(SS)
            COSHZ=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=XW(1,1)/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0D0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHZ=D+1.0D0/(3.0D0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
            SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
          ENDIF
          XW(1,1)=DLOG(COSHZ+SINHZ)
        ENDIF
        DO ni=1,NIT(nb)
          IF(ITYP10(nr).EQ.2) THEN
            XW(1,NU1(ni))=XW(1,NU1(ni))/(2.0D0*RAD)
          ELSE IF(ITYP10(nr).EQ.3) THEN
            XW(1,NU1(ni))=XW(1,NU1(ni))/(3.0D0*RAD*RAD)
          ELSE IF(ITYP10(nr).EQ.4) THEN
            IF(JTYP10.EQ.2) THEN
              XW(1,NU1(ni))=XW(1,NU1(ni))/
     '          (2.0D0*FOCUS*FOCUS*SINHZ*COSHZ)
            ELSE IF(JTYP10.EQ.3) THEN
              XW(1,NU1(ni))=XW(1,NU1(ni))/((3.0D0*COSHZ*COSHZ-1.0D0)*
     '          SINHZ)/FOCUS**3
            ENDIF
          ENDIF
        ENDDO !ni
      ENDIF

      IF(DOP) THEN
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0)
            nj=NJ_LOC(njj1,njj2)
            WRITE(OP_STRING,'('' XW(nj='',I1,'',nu=1..): '',6D12.4,'
     '        //'/(18X,6D12.4))') nj,(XW(nj,nu),nu=1,NUT(NBJ(nj)))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !njj2
        ENDDO !njj1
      ENDIF !dop

      CALL EXITS('XEXW')
      RETURN
 9999 CALL ERRORS('XEXW',ERROR)
      CALL EXITS('XEXW')
      RETURN 1
      END


      SUBROUTINE XGMG(IP,JAC,nb,nr,DXIX,GL,GU,RG,XG,ERROR,*)

C#### Subroutine: XGMG
C###  Description:
C###    XGMG evaluates the covariant (GL) & contravariant (GU) metric
C###    tensors wrt the Xi-coordinate system  and  the derivs  of
C###    the Xi-coords wrt the Xj-coords (DXIX) -if NIT=NJT only - at
C###    current Gauss pt.
C**** If IP=0 DXIX contains derivatives of Xi wrt X(ref)-coords.
C**** If IP=1 DXIX contains derivatives of Xi wrt Nu(fibre)-coords.
C**** The Jacobian RG for a length,area or volume integral is returned
C****   if JAC=1,2 or 3, respec.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IP,JAC,nb,nr
      REAL*8 DXIX(3,3),GL(3,3),GU(3,3),RG,XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,ni,NITB,njj,nu
      REAL*8 AA,D,DXXI(3,3),G,G1,G3,R,RC,RR,RRC,SLX,SMX

      CALL ENTERS('XGMG',*9999)

C     Calculate derivatives of X wrt Xi
      NITB=NIT(nb)
      DO ni=1,NITB
        nu=1+ni*(1+ni)/2
        DO njj=1,NJ_LOC(NJL_GEOM,0)
          DXXI(njj,ni)=XG(njj,nu)
        ENDDO !njj
      ENDDO !ni

C     Initialise metric tensors
      DO mi=1,3
        DO ni=1,3
          GL(mi,ni)=0.0d0
          GU(mi,ni)=0.0d0
        ENDDO !ni
        GL(mi,mi)=1.0d0
        GU(mi,mi)=1.0d0
      ENDDO !mi

C     Calculate covariant metric tensor GL(i,j)
      IF(ITYP10(nr).EQ.2) THEN       !cyl polar
        R=XG(1,1)
        RR=R*R
      ELSE IF(ITYP10(nr).EQ.3) THEN  !sph polar
        R=XG(1,1)
        RR=R*R
        RC=R*DCOS(XG(3,1))
        RRC=RC*RC
      ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate sph
        IF( DABS(XG(2,1)).LT.RDELTA) THEN
          RG=0.0d0
          WRITE(OP_STRING,'('' >>Warning: mu is zero in XGMG'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GO TO 9998
        ENDIF
        AA=FOCUS*FOCUS
        SLX=DSINH(XG(1,1))
        SMX=DSIN(XG(2,1))
        G1=AA*(SLX*SLX+SMX*SMX)
        G3=AA* SLX*SLX*SMX*SMX
      ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate sph
      ENDIF

      DO mi=1,NITB
        DO ni=1,NITB
          IF(ITYP10(nr).NE.4) GL(mi,ni)=DXXI(1,mi)*DXXI(1,ni)
          IF(NJ_LOC(NJL_GEOM,0).GT.1) THEN
            IF(ITYP10(nr).EQ.1) THEN       !rect cart
              DO njj=2,NJ_LOC(NJL_GEOM,0)
                GL(mi,ni)=GL(mi,ni)+DXXI(njj,mi)*DXXI(njj,ni)
              ENDDO !njj
            ELSE IF(ITYP10(nr).EQ.2) THEN  !cyl polar
              GL(mi,ni)=GL(mi,ni)+RR*DXXI(2,mi)*DXXI(2,ni)
              IF(NJ_LOC(NJL_GEOM,0).EQ.3)
     '          GL(mi,ni)=GL(mi,ni)+DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(nr).EQ.3) THEN  !sph polar
              GL(mi,ni)=GL(mi,ni)+RRC*DXXI(2,mi)*DXXI(2,ni)
     '                             +RR *DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate sph
              GL(mi,ni)=G1*(DXXI(1,mi)*DXXI(1,ni)
     '                               +DXXI(2,mi)*DXXI(2,ni))
              IF(NJ_LOC(NJL_GEOM,0).EQ.3) GL(mi,ni)=GL(mi,ni)+G3
     '                               *DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate sph
            ENDIF
          ENDIF
        ENDDO !ni
      ENDDO !mi

C new MPN 17-Apr-96: calc GU and DXIX with calls to INVERT and DXIDNU

C     Calculate contravariant metric tensor GU(i,j)
      CALL INVERT(NITB,GL,GU,G)
      IF(DABS(G).LT.RDELTA) THEN
        RG=0.0d0
        WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
     '    //''' RDELTA='',D12.5)') G,RDELTA
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        GOTO 9998
      ENDIF

C     Calculate derivs DXIX(i,j) of Xi wrt X (IP=0) or Nu (IP=1)
      IF(IP.EQ.0) THEN !DXIX is based on reference coords, X
        IF(NITB.EQ.NJ_LOC(NJL_GEOM,0)) CALL INVERT(NITB,DXXI,DXIX,D)
        IF(DOP) THEN
          DO mi=1,NITB
            WRITE(OP_STRING,'('' DXIX('',I1,'',ni)   : '',3E12.4)')
     '        mi,(DXIX(mi,ni),ni=1,NITB)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !mi
        ENDIF
      ELSE IF(IP.GE.1) THEN !DXIX is based on material fibre coords, Nu
        CALL DXIDNU(nb,nr,DXIX,DXXI,GL,GU,XG,ERROR,*9999)
      ENDIF

C old way of handling fibre/sheets
C old DATA M/1,2,3,1,2/
C
C      IF(NITB.EQ.1) THEN
C        G=GL(1,1)
C        IF(DABS(G).GT.RDELTA) THEN
C          GU(1,1)=1.0d0/G
C        ELSE
C          RG=0.0d0
C          WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
C     '      //''' RDELTA='',D12.5)') G,RDELTA
C          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C          GOTO 9998
C        ENDIF
C        IF(IP.EQ.0) THEN !DXIX is based on reference coords
C          IF(DABS(DXXI(1,1)).GT.RDELTA) THEN
C            DXIX(1,1)=1.D0/DXXI(1,1)
C          ENDIF
C        ELSE IF(IP.GE.1) THEN !DXIX is based on material coords
C                              !.. ie arclength along 1D element
C          DXIX(1,1)=DSQRT(GL(1,1))
C        ENDIF
C      ELSE IF(NITB.EQ.2) THEN
C        G=GL(1,1)*GL(2,2)-GL(1,2)*GL(2,1)
C        IF(DABS(G).LT.RDELTA) THEN
C          RG=0.0d0
C          WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
C     '      //''' RDELTA='',D12.5)') G,RDELTA
C          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C          GOTO 9998
C        ENDIF
C        GU(1,1)= GL(2,2)/G
C        GU(1,2)=-GL(1,2)/G
C        GU(2,1)=-GL(2,1)/G
C        GU(2,2)= GL(1,1)/G
C        GU(3,3)=1.0d0
C        IF(IP.EQ.0) THEN !DXIX is based on reference coords
C          IF(NITB.EQ.NJ_LOC(NJL_GEOM,0)) THEN
C            D=DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)
C            IF(DABS(D).GT.RDELTA) THEN
CC             calc DXIX(ni,nj) from DXXI(nj,ni)
C              DXIX(1,1)= DXXI(2,2)/D
C              DXIX(1,2)=-DXXI(1,2)/D
C              DXIX(2,1)=-DXXI(2,1)/D
C              DXIX(2,2)= DXXI(1,1)/D
C            ENDIF
C          ENDIF
C        ELSE IF(IP.GE.1) THEN !DXIX is based on material coords
C          IF(JTYP9.EQ.0) THEN
C            ETA1=0.0D0
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C            C1=1.0d0
C            C2=DCOS(ETA2)
C          ELSE IF(JTYP9.GE.1) THEN
C            IF(JTYP12.EQ.1) THEN
C              ETA1=XG(NJ_LOC(NJL_FIBR,1),1)
C              ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C            ELSE IF(JTYP12.EQ.2) THEN
C              ETA2=-XG(NJ_LOC(NJL_FIBR,1),1)
C              ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C            ENDIF
C            C1=DCOS(ETA1)
C            C2=DCOS(ETA2)
C          ENDIF
C          RK1=DSQRT(GL(1,1))*C1
C          RK2=DSQRT(GL(2,2))*C2
C          RDSQ=GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2)
C          RD=DSQRT(RDSQ)
C          RGU33=DSQRT(GU(3,3))
C          DXIX(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C          DXIX(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C          DXIX(3,1)= 0.0D0
C          DXIX(1,2)=-RK2/RD
C          DXIX(2,2)= RK1/RD
C          DXIX(3,2)= 0.0D0
C          DXIX(1,3)= 0.0D0
C          DXIX(2,3)= 0.0D0
C          DXIX(3,3)= RGU33
C        ENDIF
C      ELSE IF(NITB.EQ.3) THEN
C        G=DET(GL)         !Note DET should handle dble precision.
C        IF(DABS(G).LT.RDELTA) THEN
C          RG=0.0D0
C          WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
C     '      //''' RDELTA='',D12.5)') G,RDELTA
C          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C          GOTO 9998
C        ENDIF
C        DO mi=1,3
C          DO ni=1,3
C            GU(mi,ni)=(GL(M(ni+1),M(mi+1))*GL(M(ni+2),M(mi+2))
C     '                -GL(M(ni+2),M(mi+1))*GL(M(ni+1),M(mi+2)))/G
C          ENDDO
C        ENDDO
C        IF(IP.EQ.0) THEN
C          IF(NITB.EQ.NJ_LOC(NJL_GEOM,0)) THEN
C            D=DXXI(1,1)*
C     '        (DXXI(2,2)*DXXI(3,3)-DXXI(3,2)*DXXI(2,3))
C     '       +DXXI(1,2)*
C     '        (DXXI(2,3)*DXXI(3,1)-DXXI(3,3)*DXXI(2,1))
C     '       +DXXI(1,3)*
C     '        (DXXI(2,1)*DXXI(3,2)-DXXI(3,1)*DXXI(2,2))
C            DO ni=1,NITB
C              DO njj=1,NJ_LOC(NJL_GEOM,0)
C                nj=NJ_LOC(NJL_GEOM,njj)
C                DXIX(ni,nj)=(DXXI(M(nj+1),M(ni+1))*DXXI(M(nj+2),
C     '              M(ni+2))-DXXI(M(nj+2),M(ni+1))*DXXI(M(nj+1),
C     '              M(ni+2)))/D
C              ENDDO
C            ENDDO
C          ENDIF
C        ELSE IF(IP.GE.1) THEN
C          IF(JTYP9.EQ.0) THEN
C            ETA1=0.0D0
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C            C1=1.0d0
C            C2=DCOS(ETA2)
C          ELSE IF(JTYP9.GE.1) THEN
C            IF(JTYP12.EQ.1) THEN
C              ETA1=XG(NJ_LOC(NJL_FIBR,1),1)
C              ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C            ELSE IF(JTYP12.EQ.2) THEN
C              ETA2=-XG(NJ_LOC(NJL_FIBR,1),1)
C              ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C            ENDIF
C            C1=DCOS(ETA1)
C            C2=DCOS(ETA2)
C          ENDIF
C          RK1=DSQRT(GL(1,1))*C1
C          RK2=DSQRT(GL(2,2))*C2
C          RDSQ=(GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2))
C          RD=DSQRT(RDSQ)
C          RGU33=DSQRT(GU(3,3))
C          DXIX(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C          DXIX(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C          DXIX(3,1)= 0.0D0
C          DXIX(1,2)=-RK2/RD
C          DXIX(2,2)= RK1/RD
C          DXIX(3,2)= 0.0D0
C          DXIX(1,3)=RGU33*(GL(1,2)*GL(2,3)-GL(2,2)*GL(1,3))/RDSQ
C          DXIX(2,3)=RGU33*(GL(1,2)*GL(1,3)-GL(1,1)*GL(2,3))/RDSQ
C          DXIX(3,3)=RGU33
C        ENDIF
C      ENDIF

C     Calculate Jacobian RG
      IF(JAC.GT.0) THEN
        IF(JAC.EQ.1) RG=DSQRT(DABS(GL(1,1)))
        IF(JAC.EQ.2) RG=DSQRT(DABS(G*GU(3,3)))
        IF(JAC.EQ.3) RG=DSQRT(DABS(G))
      ENDIF

 9998 CALL EXITS('XGMG')
      RETURN
 9999 CALL ERRORS('XGMG',ERROR)
      CALL EXITS('XGMG')
      RETURN 1
      END


      SUBROUTINE ZEZG(JP,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*)

C#### Subroutine: ZEZG
C###  Description:
C###    ZEZG evaluates the Gauss point array ZG(nhx,nu) from the
C###    element array ZE(ns,nhx) at the current Gauss point ng.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER JP,NBH(NHM),ng,NHE,nx
      REAL*8 DXIX(3,3),PG(NSM,NUM,NGM,NBM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,nb,nh,nhx,ni,ns,nu,NU1(0:3),NU2(3,3)
      REAL*8 COSHZ,CSS,D,DES,PGG,RAD,SINHZ,SS,SUM,THETA

      DATA NU1/1,2,4,7/,NU2/3,6,10,6,5,9,10,9,8/

      CALL ENTERS('ZEZG',*9999)
      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh)
C MPN 18-Jul-95 Initialising ZG
        DO nu=1,NUT(nb)
          ZG(nhx,nu)=0.0d0
        ENDDO
        DO ni=0,NIT(nb)
          SUM=0.0D0
          DO ns=1,NST(nb)+NAT(nb)
            IF(jp.eq.0.or.ni.EQ.0) THEN
              PGG=PG(ns,NU1(ni),ng,nb)
            ELSE IF(JP.GT.0) THEN
              IF(NIT(nb).EQ.1) THEN
                PGG=PG(ns,2,ng,nb)*DXIX(1,1)
              ELSE IF(NIT(nb).EQ.2) THEN
                PGG=PG(ns,2,ng,nb)*DXIX(1,ni) +
     '              PG(ns,4,ng,nb)*DXIX(2,ni)
              ELSE IF(NIT(nb).EQ.3) THEN
                PGG=PG(ns,2,ng,nb)*DXIX(1,ni) +
     '              PG(ns,4,ng,nb)*DXIX(2,ni) +
     '              PG(ns,7,ng,nb)*DXIX(3,ni)
C              ELSE
C                PGG=PGX(nb,ni,ns,DXIX,PG(1,1,ng,nb))
              ENDIF
            ENDIF
            SUM=SUM+PGG*ZE(ns,nhx)
          ENDDO
          ZG(nhx,NU1(ni))=SUM
        ENDDO
      ENDDO

      IF(JTYP10.GE.2) THEN
        nb=NBH(NH_LOC(1,nx))
        IF(ITYP11(1).EQ.2) THEN
          RAD=DSQRT(ZG(1,1))
          ZG(1,1)=RAD
        ELSE IF(ITYP11(1).EQ.3) THEN
          RAD=ZG(1,1)**(1.0D0/3.0D0)
          ZG(1,1)=RAD
        ELSE IF(ITYP11(1).EQ.4) THEN
          IF(JTYP10.EQ.2) THEN
            SS=ZG(1,1)/(FOCUS*FOCUS)
            SINHZ=DSQRT(SS)
            COSHZ=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=ZG(1,1)/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHZ=D+1.0d0/(3.0d0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
            SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
          ENDIF
          ZG(1,1)=DLOG(COSHZ+SINHZ)
        ENDIF
        DO ni=1,NIT(nb)
          IF(ITYP11(1).EQ.2) THEN
            ZG(1,NU1(ni))=ZG(1,NU1(ni))/(2.0D0*RAD)
          ELSE IF(ITYP11(1).EQ.3) THEN
            ZG(1,NU1(ni))=ZG(1,NU1(ni))/(3.0D0*RAD*RAD)
          ELSE IF(ITYP11(1).EQ.4) THEN
            IF(JTYP10.EQ.2) THEN
              ZG(1,NU1(ni))=ZG(1,NU1(ni))/
     '          (2.0D0*FOCUS*FOCUS*SINHZ*COSHZ)
            ELSE IF(JTYP10.EQ.3) THEN
              ZG(1,NU1(ni))=ZG(1,NU1(ni))/((3.0D0*COSHZ*COSHZ-1.0D0)
     '          *SINHZ)/FOCUS**3
            ENDIF
          ENDIF
        ENDDO
        DO ni=1,NIT(nb)
          DO mi=1,ni
            nu=NU2(ni,mi)
            IF(ITYP11(1).EQ.2) THEN
              ZG(1,nu)=-ZG(1,nu)/(4.0D0*RAD**3)
            ELSE IF(ITYP11(1).EQ.3) THEN
              ZG(1,nu)=-ZG(1,nu)*2.0D0/(9.0D0*RAD**5)
            ELSE IF(ITYP11(1).EQ.4) THEN
              IF(JTYP10.EQ.2) THEN
                ZG(1,nu)=-ZG(1,nu)*(COSHZ*COSHZ+SINHZ*SINHZ)
     '            /(4.0D0*FOCUS**4*SINHZ**3*COSHZ**3)
              ELSE IF(JTYP10.EQ.3) THEN
                ZG(1,nu)=-ZG(1,nu)*(2.0D0+9.0D0*SINHZ*SINHZ)*COSHZ/
     '            (FOCUS**6*(3.0D0*COSHZ*COSHZ-1.0D0)**3*SINHZ**3)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF(DOP) THEN
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,'('' ZG(nhx='',I1,'',nu=1..): '',6D12.4,'
     '      //'/(19X,6D12.4))') nhx,(ZG(nhx,nu),nu=1,NUT(NBH(nh)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('ZEZG')
      RETURN
 9999 CALL ERRORS('ZEZG',ERROR)
      CALL EXITS('ZEZG')
      RETURN 1
      END


      SUBROUTINE ZEZW(IP,IBT,IDO,INP,NAN,NBH,NHTOT,nr,nx,DXIX,ZE,ZW,XI,
     '  ERROR,*)

C#### Subroutine: ZEZW
C###  Description:
C###    ZEZW interpolates ZE(nhx,1) and first derivs ZE(nhx,NU1(ni))
C###    at XI for Lagrange/Hermite tensor product basis functions with
C###    auxiliary element parameters and returns ZW(nhx,NU1(ni)). ZW
C###    is corrected for JTYP10=2. Derivatives are wrt Xi if IP=0 or
C###    wrt nu if IP=1.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IP,
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NHTOT,nr,nx
      REAL*8 DXIX(3,3),XI(3),ZE(NSM,NHM),ZW(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,nb,nh,nhx,ni,NITB,nu,NU1(0:3)
      REAL*8 COSHZ,CSS,D,DES,DZHX(3),PFXI,RAD,SINHZ,SS,THETA

      DATA NU1/1,2,4,7/

      CALL ENTERS('ZEZW',*9999)
      NITB=NIT(NBH(NH_LOC(1,nx)))
      DO ni=0,NITB
        DO nhx=1,NHTOT
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh)
!new      MPN 19/7/93 array overrun for NNT in PFXI if nb=0
          IF(nb.GT.0) THEN
            ZW(nhx,NU1(ni))=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        NAN(1,1,nb),nb,NU1(ni),ZE(1,nhx),XI)
          ELSE
            ZW(nhx,NU1(ni))=0.0d0 !should never get into here!!!
            WRITE(OP_STRING,'('' >>Warning: nb is zero in ZEZW'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nhx
      ENDDO !ni

      IF(IP.EQ.1) THEN
        DO nhx=1,NHTOT
          nh=NH_LOC(nhx,nx)
          DO mi=1,NITB
            DZHX(mi)=0.0D0
            DO ni=1,NITB
              DZHX(mi)=DZHX(mi)+ZW(nhx,NU1(ni))*DXIX(ni,mi)
            ENDDO !ni
          ENDDO !mi
          DO mi=1,NITB
            ZW(nhx,NU1(mi))=DZHX(mi)
          ENDDO !mi
        ENDDO !nhx
      ENDIF

      IF(JTYP10.GE.2) THEN
        nb=NBH(NH_LOC(1,nx))
        IF(ITYP10(nr).EQ.2) THEN
          RAD=DSQRT(ZW(1,1))
          ZW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.3) THEN
          RAD=ZW(1,1)**(1.0D0/3.0D0)
          ZW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.4) THEN
          IF(JTYP10.EQ.2) THEN
            SS=ZW(1,1)/(FOCUS*FOCUS)
            SINHZ=DSQRT(SS)
            COSHZ=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=ZW(1,1)/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0D0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHZ=D+1.0D0/(3.0D0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
            SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
          ENDIF
          ZW(1,1)=DLOG(COSHZ+SINHZ)
        ENDIF
        DO ni=1,NIT(nb)
          IF(ITYP10(nr).EQ.2) THEN
            ZW(1,NU1(ni))=ZW(1,NU1(ni))/(2.0D0*RAD)
          ELSE IF(ITYP10(nr).EQ.3) THEN
            ZW(1,NU1(ni))=ZW(1,NU1(ni))/(3.0D0*RAD*RAD)
          ELSE IF(ITYP10(nr).EQ.4) THEN
            IF(JTYP10.EQ.2) THEN
              ZW(1,NU1(ni))=ZW(1,NU1(ni))/
     '          (2.0D0*FOCUS*FOCUS*SINHZ*COSHZ)
            ELSE IF(JTYP10.EQ.3) THEN
              ZW(1,NU1(ni))=ZW(1,NU1(ni))/((3.0D0*COSHZ*COSHZ-1.0D0)*
     '          SINHZ)/FOCUS**3
            ENDIF
          ENDIF
        ENDDO !ni
      ENDIF

      IF(DOP) THEN
        DO nhx=1,NHTOT
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh)
          WRITE(OP_STRING,'('' ZW(nhx='',I1,'',nu=1..):'',6D12.4,'
     '      //'/(18X,6D12.4))') nhx,(ZW(nhx,nu),nu=1,NUT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nhx
      ENDIF !dop

      CALL EXITS('ZEZW')
      RETURN
 9999 CALL ERRORS('ZEZW',ERROR)
      CALL EXITS('ZEZW')
      RETURN 1
      END


      SUBROUTINE ZGMG(nb,GZ,GZL,GZU,ZG,ERROR,*)

C#### Subroutine: ZGMG
C###  Description:
C###    ZGMG evaluates components of metric tensor in deformed state
C###    at current Gauss point from coordinates ZG.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER nb
      REAL*8 GZ,GZL(3,3),GZU(3,3),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,nhx,ni,NITB,NU1(0:3)
      REAL*8 AA,G1,G3,SLZ,SMZ,SUM
      CHARACTER CHAR1*1

      DATA NU1/1,2,4,7/

      CALL ENTERS('ZGMG',*9999)
      NITB=NIT(nb)
      DO mi=1,NITB
        DO ni=1,NITB
          SUM=ZG(1,NU1(mi))*ZG(1,NU1(ni))
          IF(NJ_LOC(NJL_GEOM,0).GT.1) THEN
            IF(ITYP11(1).EQ.1) THEN
              DO nhx=2,NJ_LOC(NJL_GEOM,0)
                SUM=SUM+ZG(nhx,NU1(mi))*ZG(nhx,NU1(ni))
              ENDDO
            ELSE IF(ITYP11(1).EQ.2) THEN
              SUM=SUM+ZG(1,1)**2*ZG(2,NU1(mi))*ZG(2,NU1(ni))
              IF(NJ_LOC(NJL_GEOM,0).EQ.3)
     '          SUM=SUM+ZG(3,NU1(mi))*ZG(3,NU1(ni))
            ELSE IF(ITYP11(1).EQ.3) THEN
              SUM=SUM+ZG(1,1)**2*(DCOS(ZG(3,1))**2*ZG(2,NU1(mi))
     '               *ZG(2,NU1(ni))+ZG(3,NU1(mi))*ZG(3,NU1(ni)))
            ELSE IF(ITYP11(1).EQ.4) THEN
              AA=FOCUS*FOCUS
              SLZ=DSINH(ZG(1,1))
              SMZ=DSIN(ZG(2,1))
              G1=AA*(SLZ*SLZ+SMZ*SMZ)
              G3=AA*SLZ*SLZ*SMZ*SMZ
              SUM=G1*(SUM+ZG(2,NU1(mi))*ZG(2,NU1(ni)))
              IF(NJ_LOC(NJL_GEOM,0).EQ.3)
     '          SUM=SUM+G3*ZG(3,NU1(mi))*ZG(3,NU1(ni))
            ENDIF
          ENDIF
          GZL(mi,ni)=SUM
        ENDDO
      ENDDO

C new MPN 17-Apr-96: calc GZU with call to INVERT
C     Calculate contravariant metric tensor GZU(i,j)
      CALL INVERT(NITB,GZL,GZU,GZ)
      IF(DABS(GZ).LT.RDELTA) THEN
        WRITE(OP_STRING,'('' >>Warning: zero GZ in ZGMG. GZ='',D12.5,'
     '    //''' RDELTA='',D12.5)') GZ,RDELTA
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF
C old
C      IF(NITB.EQ.1) THEN
C        GZ=GZL(1,1)
C        GZU(1,1)=1.0D0/GZ
C      ELSE IF(NITB.EQ.2) THEN
C        GZ=GZL(1,1)*GZL(2,2)-GZL(1,2)*GZL(2,1)
C        GZU(1,1)= GZL(2,2)/GZ
C        GZU(1,2)=-GZL(1,2)/GZ
C        GZU(2,1)=-GZL(2,1)/GZ
C        GZU(2,2)= GZL(1,1)/GZ
C        GZU(3,3)=1.0d0
C      ELSE IF(NITB.EQ.3) THEN
C        CALL INVERT(3,GZL,GZU,GZ)
C      ENDIF

      IF(DOP) THEN
        WRITE(CHAR1,'(I1)') NITB
        WRITE(OP_STRING,'(''  GZL or AZL:'','//CHAR1(1:1)//'D12.4,'
     '    //'''  GZU or AZU:'','//CHAR1(1:1)//'D12.4,'
     '    //'/(13X,'//CHAR1(1:1)//'D12.4,13X,'//CHAR1(1:1)//'D12.4))')
     '    ((GZL(mi,ni),ni=1,NITB),(GZU(mi,ni),ni=1,NITB),mi=1,NITB)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C old format
C        WRITE(OP_STRING,'('' GZL or AZL by rows: '',9E12.4)')
C     '    ((GZL(mi,ni),ni=1,NITB),mi=1,NITB)
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' GZU or AZU by rows: '',9E12.4)')
C     '    ((GZU(mi,ni),ni=1,NITB),mi=1,NITB)
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ZGMG')
      RETURN
 9999 CALL ERRORS('ZGMG',ERROR)
      CALL EXITS('ZGMG')
      RETURN 1
      END


      REAL*8 FUNCTION PAF(K,na,XI)

C#### Function: PAF
C###  Type: REAL*8
C###  Description:
C###    PAF evaluates 1D Fourier auxiliary basis (sine) functions at Xi.

C**** na is half the wavenumber of the function
C**** K is the derivative reqd. (K=1 -> value reqd)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
!     Parameter List
      INTEGER K,na
      REAL*8 XI
!     Local Variables
      REAL*8 PAL0,RNA

      IF(na.EQ.0) THEN        !Use constant basis fn
        PAF=PAL0(K)
      ELSE
        RNA=DBLE(na)
        GO TO (10,20,30),K
 10       PAF=               DSIN(RNA*PI*XI)
          RETURN
 20       PAF=        RNA*PI*DCOS(RNA*PI*XI)
          RETURN
 30       PAF=-RNA*RNA*PI*PI*DSIN(RNA*PI*XI)
          RETURN
      ENDIF

      RETURN
      END


      REAL*8 FUNCTION PAL0(K)

C#### Function: PAL0
C###  Type: REAL*8
C###  Description:
C###    PAL0 evaluates 1D constant auxiliary basis function.

      IMPLICIT NONE
!     Parameter List
      INTEGER K

      GO TO (10,20,30),K
 10     PAL0=1.0d0
        RETURN
 20     PAL0=0.0d0
        RETURN
 30     PAL0=0.0d0
        RETURN
      END


      REAL*8 FUNCTION PAL1(K,XI)

C#### Function: PAL1
C###  Type: REAL*8
C###    PAL1 evaluates 1D linear Legendre auxiliary basis function
C###    at Xi.

      IMPLICIT NONE
!     Parameter List
      INTEGER K
      REAL*8 XI

      GO TO (10,20,30),K
 10     PAL1=XI
        RETURN
 20     PAL1=1.0d0
        RETURN
 30     PAL1=0.0d0
        RETURN
      END


      REAL*8 FUNCTION PAL2(K,na,XI)

C#### Function: PAL2
C###  Type: REAL*8
C###  Description:
C###    PAL2 evaluates 1D Legendre auxiliary basis functions of
C###    degree na, at Xi, for na>=2

      IMPLICIT NONE
!     Parameter List
      INTEGER K,na
      REAL*8 XI
!     Local Variables
      INTEGER N,NPOWER
      REAL*8 FACT,XIA,XINA
      LOGICAL EVEN

      XIA=2.0d0*XI-1.0d0
      NPOWER=na+1-K
      IF( DABS(XIA).LT.1.0d-08) THEN
        IF(NPOWER.LE.0) THEN
          XINA=1.0d0
        ELSE
          XINA=0.0d0
        ENDIF
      ELSE
        XINA=XIA**NPOWER
      ENDIF
      FACT=1.0d0
      DO N=1,na
        FACT=FACT*N
      ENDDO
      EVEN=.TRUE.
      IF((-1)**na.LT.0) EVEN=.FALSE.
      GO TO (10,20,30),K
 10   IF(EVEN) THEN
        PAL2=(XINA-1.0d0)/FACT
      ELSE
        PAL2=(XINA-XIA)/FACT
      ENDIF
      RETURN
 20   IF(EVEN) THEN
        PAL2=2.0d0*na/FACT*XINA
      ELSE
        PAL2=2.0d0/FACT*(na*XINA-1.0d0)
      ENDIF
      RETURN
 30   PAL2=4.0D0*na*(na-1)/FACT*XINA
      RETURN
      END


      REAL*8 FUNCTION PAP4(K,XI)

C#### Function: PAP4
C###  Type: REAL*8
C###  Description:
C###    PAP4 evaluates 1D Pressure quartic (hat) auxiliary basis
C###    functionn at Xi.  Values and first derivatives are zero at
C###    Xi=0,1.

      IMPLICIT NONE
!     Parameter List
      INTEGER K
      REAL*8 XI

      GO TO (10,20,30),K
 10     PAP4=XI*XI*(XI-1.0D0)*(XI-1.0D0)
        RETURN
 20     PAP4=2.0D0*XI*(XI-1.0D0)*(2.0D0*XI-1.0D0)
        RETURN
 30     PAP4=12.0D0*XI*XI - 12.0D0*XI + 2.0D0
        RETURN
      END


      REAL*8 FUNCTION PF1(na,nu,XI)

C#### Function: PF1
C###  Type: REAL*8
C###  Description:
C###    PF1 evaluates Fourier basis function number na at Xi.  nu
C###    specifies time derivative. OMEGA is angular frequency.
C###    If nu=-1 PF1 returns integral terms (from xi=0 so integral of
C###    sin term has integration constant added).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:four00.cmn'
!     Parameter List
      INTEGER na,nu
      REAL*8 XI

      IF(na.EQ.1) THEN              !constant term
        IF(nu.EQ.1) THEN            !0th-order time derivative
          PF1=1.0D0
        ELSE IF(nu.EQ.2) THEN     !1st-order time derivative
          PF1=0.0D0
        ELSE IF(nu.EQ.3) THEN     !2nd-order time derivative
          PF1=0.0D0
        ELSE IF(nu.EQ.-1) THEN     !time integral
          PF1=XI
        ENDIF

      ELSE IF(MOD(na,2).EQ.0) THEN  !cosine terms
        IF(nu.EQ.1) THEN            !0th-order time derivative
          PF1=DCOS((DBLE(na)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.2) THEN     !1st-order time derivative
          PF1=-(DBLE(na)/2.0D0)*OMEGA*DSIN((DBLE(na)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.3) THEN     !2nd-order time derivative
          PF1=-(DBLE(na)/2.0D0)*OMEGA*(DBLE(na)/2.0D0)*OMEGA*
     '      DCOS((DBLE(na)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.-1) THEN     !time integral
          PF1=DSIN((DBLE(na)/2.0D0)*OMEGA*XI)/(DBLE(na)/2.0D0*OMEGA)
        ENDIF

      ELSE                         !sine terms
        IF(nu.EQ.1) THEN            !0th-order time derivative
          PF1=DSIN((DBLE(na-1)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.2) THEN     !1st-order time derivative
          PF1=(DBLE(na-1)/2.0D0)*OMEGA*DCOS((DBLE(na-1)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.3) THEN     !2nd-order time derivative
          PF1=-(DBLE(na-1)/2.0D0)*OMEGA*(DBLE(na-1)/2.0D0)*OMEGA*
     '      DSIN((DBLE(na-1)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.-1) THEN     !time integral
          PF1=(1.0D0-DCOS(DBLE(na-1)/2.0D0*OMEGA*XI))/
     '      (DBLE(na-1)/2.0D0*OMEGA)
        ENDIF
      ENDIF

      RETURN
      END


      REAL*8 FUNCTION PFXI(IBT,IDO,INP,NAN,nb,nu,XE,XI)

C#### Function: PFXI
C###  Type: REAL*8
C###  Description:
C###    PFXI interpolates nodal array XE at XI.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),NAN(NIM,NAM),
     '  nb,nu
      REAL*8  XE(NSM),XI(3)
!     Local Variables
      INTEGER na,nn,nk,ns
      REAL*8 PSI1,PSI8

      PFXI=0.0D0
      ns=0
      DO nn=1,NNT(nb)
        DO nk=1,NKT(nn,nb)
          ns=ns+1
          PFXI=PFXI+PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)*XE(ns)
        ENDDO
      ENDDO
      DO na=1,NAT(nb)
        ns=ns+1
        PFXI=PFXI+PSI8(NAN(1,na),nb,nu,XI)*XE(ns)
      ENDDO

      RETURN
      END


      REAL*8 FUNCTION PGX(nb,nj,ns,DXIX,PG)

C#### Function: PGX
C###  Type: REAL*8
C###  Description:
C###    PGX calculates the first partial derivative of the ns-th
C###    term basis function PG with respect to the Xj-th coordinate.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER nb,nj,ns
      REAL*8 DXIX(3,3),PG(NSM,NUM)
!     Local Variables
      INTEGER ni,nu

      PGX=0.0D0
      DO ni=1,NIT(nb)
        nu=1+ni*(1+ni)/2
        PGX=PGX+PG(ns,nu)*DXIX(ni,nj)
      ENDDO

      RETURN
      END


      REAL*8 FUNCTION PH2(I,J,K,N,XI)

C#### Function: PH2
C###  Type: REAL*8
C###  Description:
C###    PH2 evaluates 1D quadratic Hermite basis function at XI.

C**** I is node position index
C**** J is nodal derivative index
C**** K is partial derivative index
C**** N is the local node number with no derivative term.
C**** XI is Xi-coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER I,J,K,N
      REAL*8 XI

      GOTO (1000,2000) N
 1000   GOTO (1100,1200,1300),K
 1100     GOTO (1110,1120),I
 1110       GOTO (1111,1112),J
 1111         PH2=(XI-2.0d0)*XI+1.0d0       ! xi^2-2xi+1
              RETURN
 1112         PH2=0.0d0                     ! 0
              RETURN
 1120       GOTO (1121,1122),J
 1121         PH2=(2.0d0-XI)*XI             ! -xi^2+2xi
              RETURN
 1122         PH2=(XI-1.0d0)*XI             ! xi^2-xi
              RETURN
 1200      GOTO (1210,1220),I
 1210        GOTO (1211,1212),J
 1211          PH2=2.0d0*XI-2.0d0           ! 2xi-2
               RETURN
 1212          PH2=0.0d0                    ! 0
               RETURN
 1220        GOTO (1221,1222),J
 1221          PH2=-2.0d0*XI+2.0d0          ! -2xi+2
               RETURN
 1222          PH2=2.0d0*XI-1.0d0           ! 2xi-1
               RETURN
 1300      GOTO (1310,1320),I
 1310        GOTO (1311,1312),J
 1311          PH2=2.0d0                    ! 2
               RETURN
 1312          PH2=0.0d0                    ! 0
               RETURN
 1320        GOTO (1321,1322),J
 1321          PH2=-2.0d0                   ! -2
               RETURN
 1322          PH2=2.0d0                    ! 2
               RETURN
 2000   GOTO (2100,2200,2300),K
 2100     GOTO (2110,2120),I
 2110       GOTO (2111,2112),J
 2111         PH2=1.0d0-XI*XI               ! -xi^2+1
              RETURN
 2112         PH2=XI*(1.0d0-XI)             ! -xi^2+xi
              RETURN
 2120       GOTO (2121,2122),J
 2121         PH2=XI*XI                     ! xi^2
              RETURN
 2122         PH2=0.0d0                     ! 0
              RETURN
 2200      GOTO (2210,2220),I
 2210        GOTO (2211,2212),J
 2211          PH2=-2.0d0*XI                ! -2xi
               RETURN
 2212          PH2=1.0d0-2.0d0*XI           ! -2xi+1
               RETURN
 2220        GOTO (2221,2222),J
 2221          PH2=2.0d0*XI                 ! 2xi
               RETURN
 2222          PH2=0.0d0                    ! 0
               RETURN
 2300      GOTO (2310,2320),I
 2310        GOTO (2311,2312),J
 2311          PH2=-2.0d0                   ! -2
               RETURN
 2312          PH2=-2.0d0                   ! -2
               RETURN
 2320        GOTO (2321,2322),J
 2321          PH2=2.0d0                    ! 2
               RETURN
 2322          PH2=0.0d0                    ! 0
               RETURN
      END


      REAL*8 FUNCTION PH3(I,J,K,XI)

C#### Function: PH3
C###  Type: REAL*8
C###  Description:
C###    PH3 evaluates 1D cubic Hermite basis function at XI.

C**** I is node position index
C**** J is nodal derivative index
C**** K is partial derivative index
C**** XI is Xi-coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER I,J,K
      REAL*8 XI

c cpb 28/6/95 Optimising mulitplication operations

      GO TO (100,200,300),K
 100    GO TO (110,120),I
 110      GO TO (111,112),J
C 111        PH3=1.0d0+XI*XI*(2.0d0*XI-3.0d0)
 111        PH3=(XI+XI-3.0d0)*XI*XI+1.0d0  ! 2xi^3-3xi^2+1
            RETURN
C 112        PH3=XI*(XI-1.0d0)*(XI-1.0d0)
 112        PH3=((XI-2.0d0)*XI+1.0d0)*XI   ! xi^3-2xi^2+xi
            RETURN
 120      GO TO (121,122),J
C 121        PH3=XI*XI*(3.0d0-2.0d0*XI)
 121        PH3=XI*XI*(3.0d0-XI-XI)        ! -2xi^3+3xi^2
            RETURN
 122        PH3=XI*XI*(XI-1.0d0)           ! xi^3-xi^2
            RETURN
 200    GO TO (210,220),I
 210      GO TO (211,212),J
 211        PH3=6.0d0*XI*(XI-1.0d0)        ! 6xi^2-6xi
            RETURN
C 212        PH3=(XI-1.0d0)*(3.0d0*XI-1.0d0)
 212        PH3=(3.0d0*XI-4.0d0)*XI+1.0d0  ! 3xi^2-4xi+1
            RETURN
 220      GO TO (221,222),J
 221        PH3=6.0d0*XI*(1.0d0-XI)        ! -6xi^2+6xi
            RETURN
 222        PH3=XI*(3.0d0*XI-2.0d0)        ! 3xi^2-2xi
            RETURN
 300    GO TO (310,320),I
 310      GO TO (311,312),J
 311        PH3=12.0d0*XI-6.0d0            ! 12xi-6
            RETURN
 312        PH3=6.0d0*XI-4.0d0             ! 6xi-4
            RETURN
 320      GO TO (321,322),J
 321        PH3=6.0d0-12.0d0*XI            ! -12xi+6
            RETURN
 322        PH3=6.0d0*XI-2.0d0             ! 6xi-2
            RETURN
      END


      REAL*8 FUNCTION PL1(I,K,XI)

C#### Function: PL1
C###  Type: REAL*8
C###  Description:
C###    PL1 evaluates 1D linear Lagrange basis function at XI.

      IMPLICIT NONE
!     Parameter List
      INTEGER I,K
      REAL*8 XI

      GO TO (10,20,30),K
 10     GO TO (11,12),I
 11       PL1=1.0d0-XI
          RETURN
 12       PL1=XI
          RETURN
 20     GO TO (21,22),I
 21       PL1=-1.0d0
          RETURN
 22       PL1=1.0d0
          RETURN
 30     PL1=0.0d0
        RETURN
      END


      REAL*8 FUNCTION PL2(I,K,XI)

C#### Function: PL2
C###  Type: REAL*8
C###  Description:
C###    PL2 evaluates 1D quadratic Lagrange basis function at XI.

      IMPLICIT NONE
!     Parameter List
      INTEGER I,K
      REAL*8 XI

      GO TO (10,20,30),K
 10     GO TO (11,12,13),I
 11       PL2=1.0d0-3.0d0*XI+2.0d0*XI*XI
          RETURN
 12       PL2=4.0d0*XI*(1.0d0-XI)
          RETURN
 13       PL2=XI*(XI+XI-1.0d0)
          RETURN
 20     GO TO (21,22,23),I
 21       PL2=4.0d0*XI-3.0d0
          RETURN
 22       PL2=4.0d0-8.0d0*XI
          RETURN
 23       PL2=4.0d0*XI-1.0d0
          RETURN
 30     GO TO (31,32,33),I
 31       PL2=4.0d0
          RETURN
 32       PL2=-8.0d0
          RETURN
 33       PL2=4.0d0
          RETURN
      END


      REAL*8 FUNCTION PL3(I,K,XI)

C#### Function: PL3
C###  Type: REAL*8
C###  Description:
C###    PL3 evaluates 1D cubic Lagrange basis funtion at XI.
C**** 08-Jan-1989: The following section added today (ADM)
C****              also corrected 2nd derivs in fn PL2 above
C****              I haven't checked the derivs very carefully!
C**** I is node position index
C**** K is partial derivative index
C**** XI is Xi-coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER I,K
      REAL*8 XI

      GO TO (10,20,30),K
 10     GO TO (11,12,13,14),I
 11       PL3=0.5D0*(3.0D0*XI-1.0D0)*(3.0D0*XI-2.0D0)*(1.0D0-XI)
          RETURN
 12       PL3=4.5D0*XI*(3.0D0*XI-2.0D0)*(XI-1.0D0)
          RETURN
 13       PL3=4.5D0*XI*(3.0D0*XI-1.0D0)*(1.0D0-XI)
          RETURN
 14       PL3=0.5D0*XI*(3.0D0*XI-1.0D0)*(3.0D0*XI-2.0D0)
          RETURN
 20     GO TO (21,22,23,24),I
 21       PL3=-13.5D0*XI*XI+18.0D0*XI-5.5D0
          RETURN
 22       PL3= 40.5D0*XI*XI-45.0D0*XI+9.0D0
          RETURN
 23       PL3=-40.5D0*XI*XI+36.0D0*XI-4.5D0
          RETURN
 24       PL3= 13.5D0*XI*XI- 9.0D0*XI+1.0D0
          RETURN
 30     GO TO (31,32,33,34),I
 31       PL3=9.0D0*(2.0D0-3.0D0*XI)
          RETURN
 32       PL3=9.0D0*(9.0D0*XI-5.0D0)
          RETURN
 33       PL3=9.0D0*(4.0D0-9.0D0*XI)
          RETURN
 34       PL3=9.0D0*(3.0D0*XI-1.0D0)
          RETURN
C**** 08-Jan-1989: End of addition (ADM)
      END


      REAL*8 FUNCTION PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)

C#### Function: PSI1
C###  Type: REAL*8
C###  Description:
C###    PSI1 evaluates tensor product Lagrange/Hermite/Fourier basis
C###    functions at XI.

C**** IPU(nu,ni),nu=1,NUT(nb) identifies the complete set of partial
C**** derivatives with respect to Xi(ni).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nk,nn,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER IPU(11,3),ni
      REAL*8 PF1,PH2,PH3,PL1,PL2,PL3

      DATA IPU/1,2,3,1,1,2,1,1,1,2,2,
     '         1,1,1,2,3,2,1,1,2,1,2,
     '         1,1,1,1,1,1,2,3,2,2,2/

      PSI1=1.0d0
      DO ni=1,NIT(nb)
        IF(IBT(1,ni).EQ.1) THEN
          IF(IBT(2,ni).EQ.1) PSI1=PSI1*PL1(INP(nn,ni),IPU(nu,ni),XI(ni))
          IF(IBT(2,ni).EQ.2) PSI1=PSI1*PL2(INP(nn,ni),IPU(nu,ni),XI(ni))
          IF(IBT(2,ni).EQ.3) PSI1=PSI1*PL3(INP(nn,ni),IPU(nu,ni),XI(ni))
        ELSE IF(IBT(1,ni).EQ.2) THEN
          IF(IBT(2,ni).EQ.1) THEN
            PSI1=PSI1*PH3(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
          ELSE
            PSI1=PSI1*PH2(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),
     '        IBT(2,ni)-1,XI(ni))
          ENDIF
        ELSE IF(IBT(1,ni).EQ.9) THEN
          PSI1=PSI1*PF1(IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
        ENDIF
      ENDDO

      RETURN
      END


      REAL*8 FUNCTION PSI8(NAN,nb,nu,XI)

C#### Function: PSI8
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    PSI8 evaluates auxiliary basis functions at Xi.
C###    <PRE>
C###    NABTYP(nb)=1 Legendre polynomials
C###               2 Fourier  basis functions
C###               3 Pressure basis functions
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER NAN(NIM),nb,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER IPU(11,3),ni
      REAL*8 PAF,PAL0,PAL1,PAL2,PAP4,PH3

      DATA IPU/1,2,3,1,1,2,1,1,1,2,2,
     '         1,1,1,2,3,2,1,1,2,1,2,
     '         1,1,1,1,1,1,2,3,2,2,2/

      PSI8=1.0D0
      DO ni=1,NIT(nb)
        IF(NABTYP(nb).EQ.1) THEN      !Legendre basis
          IF(NAN(ni).EQ.0) PSI8=PSI8*PAL0(IPU(nu,ni))
          IF(NAN(ni).EQ.1) PSI8=PSI8*PAL1(IPU(nu,ni),XI(ni))
          IF(NAN(ni).GE.2) PSI8=PSI8*PAL2(IPU(nu,ni),NAN(ni),XI(ni))
        ELSE IF(NABTYP(nb).EQ.2) THEN !Fourier  basis
          PSI8=PSI8*PAF(IPU(nu,ni),NAN(ni),XI(ni))
        ELSE IF(NABTYP(nb).EQ.3) THEN !Pressure basis
          IF(NAN(ni).LT.0) PSI8=0.0d0
          IF(NAN(ni).EQ.0) PSI8=PSI8*PAL0(IPU(nu,ni))
          IF(NAN(ni).EQ.1) PSI8=PSI8*PAL1(IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.2) PSI8=PSI8*PAL2(IPU(nu,ni),NAN(ni),XI(ni))
          IF(NAN(ni).EQ.3) PSI8=PSI8*PH3(2,1,IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.4) PSI8=PSI8*PAP4(IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.5) PSI8=PSI8*PH3(1,2,IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.6) PSI8=PSI8*PH3(2,2,IPU(nu,ni),XI(ni))
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE ZEES(IBT,IDO,INP,LGE,NAN,NBH,NBJ,
     '  ne,NFF,NGAP,NHE,NJE,NKF,NMNO,NNF,NPF,NPNE,NPNY,nr,NRE,
     '  NW,nx,NXI,NYNE,NYNP,
     '  CE,CG,CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT,
     '  PG,RE1,RE2,RG,SE,VE,WG,XE,XG,ZE,ZE1,ZG,ZG1,FIX,ERROR,*)

C#### Subroutine: ZEES
C###  Description:
C###    ZEES calculates element tangent stiffness matrix ES from
C###    current dependent variable array ZP.

C**** NW=0 : tangent stiffness matrix calculated algebraically.
C**** " >0 :    "        "        "       "      by finite diffs.
C**** " -1 : element is 'dry' (no contribution to global system).
C**** DELTA is global perturbation. DELTA*SE(ns,nb,ne) is element pert.n
C**** RE(ns,nh) has been corrected with scaling factor SE(ns,nb,ne)
C**** ES(mhs,nhs) has scaling factor correction

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ne,NFF(6),
     '  NGAP(NIM,NBM),NHE,NJE,NKF(0:4,16,6,NBFM),NMNO(1:2,0:NOPM),
     '  NNF(0:17,6,NBFM),NPF(15,NFM),NPNE(NNM,NBFM,NEFM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),NW,nx,NXI(-NIM:NIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),D_RE(NSM,NHM,NOPM),
     '  D_RI3(NHM*NSM),D_TG(3,3,NHM*NSM),D_ZG(NHM,NUM,NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM),PG(NSM,NUM,NGM,NBM),
     '  RE1(NSM,NHM),RE2(NSM,NHM),RG(NGM),
     '  SE(NSM,NBFM,NEFM),VE(NSM,NKM),WG(NGM,NBM),XE(NSM,NJM),
     '  XG(NJM,NUM),ZE(NSM,NHM),ZE1(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ik,GETNYR,mh,mhs,mhx,ms,my,my1,nb,nh,nhs,nhx,nk,nn,
     '  ns,ns1,ny,ny1
      REAL*8 COSHZ,CSS,D,DELTA,DES,R,SINHZ,SS,THETA,Y
      PARAMETER (DELTA=1.0D-4)

      CALL ENTERS('ZEES',*9999)
      IF(NW.EQ.-1) GO TO 9998
      IF(KTYP1D.EQ.1) THEN
C       finite difference calculation of tangent stiffness matrix

        IF(ITYP1(nr,nx).EQ.3) THEN      !fe30 problems
          CALL CPCG(1,NBH(1),NPNE,nr,nx,CE,CG,CP,PG,ERROR,*9999)
          CALL ZERE30(NBH,NBJ,ne,NHE,NJE,
     '      CG,PG,RE1,SE,VE,WG,XE,XG,ZE,ZG,ERROR,*9999)
        ELSE IF(ITYP1(nr,nx).EQ.4) THEN !fe40 problems
          CALL CPCG(NW,NBH(1),NPNE,nr,nx,CE,CG,CP,PG,ERROR,*9999)
          CALL ZERE40(NBH,NBJ,ne,NHE,NJE,NW,nx,CE,CG,PG,RE1,SE,VE,
     '      WG,XE,XG,ZE,ZG,ERROR,*9999)
        ELSE IF(ITYP1(nr,nx).EQ.5) THEN !fe50 problems
          CALL CPCG(1,NBH(1),NPNE,nr,nx,CE,CG,CP,PG,ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
            CALL ZERE55(INP,NBH,ne,NHE,nr,nx,
     '        CG,PG,RE1,WG,ZE,ZE1,ZG,ERROR,*9999)
          ELSE                    !finite elasticity stress analysis
            CALL ZERE50(IBT,IDO,INP,NAN,NBH,NBJ,ne,
     '        NFF,NGAP,NHE,NJE,NKF,NNF,NPF,NPNE,nr,NRE,NW,nx,NXI,
     '        CE,CG,CP,FEXT,PG,RE1,RG,SE,VE,WG,
     '        XE,XG,ZE,ZG,ERROR,*9999)
          ENDIF
        ENDIF

        nhs=0
        DO 40 nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh)
        DO 40 ns=1,NST(nb)+NAT(nb)
          nhs=nhs+1
          ny=IABS(LGE(nhs,2)) !local variable number
          ny1=GETNYR(1,NPNY,nr,0,2,ny,NYNE,NYNP) !global variable #
c PJH 19Jul93 IF(.NOT.FIX(ny1,1)) THEN
          IF(.NOT.FIX(ny1,1)) THEN  !uncommented MPN 14-Sep-94
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' ******* Element '',I5,'
     '          //'''   nh='',I2,'' ns='',I2,'' *******'')') ne,nh,ns
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(nhx.EQ.1.AND.JTYP10.GE.2) THEN
C ***         Special perturbation for isochoric interpolation
C             IF(ITYP10(1).LE.3) THEN
C               R=ZE(ns,1)**(1.0D0/ITYP10(1))
C               ZE(ns,1)=(R+DELTA*SE(ns,nb,ne))**ITYP10(1)
C             ELSE IF(ITYP10(1).EQ.4) THEN
C               IF(JTYP10.EQ.2) THEN
C                 SS=ZE(ns,1)/(FOCUS*FOCUS)
C                 SINHZ=DSQRT(SS)
C                 COSHZ=DSQRT(1.0D0+SS)
C                 R=DLOG(COSHZ+SINHZ)
C                 Y=DSINH(R+DELTA*SE(ns,nb,ne))
C                 ZE(ns,1)=(FOCUS*Y)*(FOCUS*Y)
C               ELSE IF(JTYP10.EQ.3) THEN
C                 CSS=ZE(ns,1)/FOCUS**3
C                 DES=CSS*CSS-4.0D0/27.0D0
C                 IF(DES.GT.0.0D0) THEN
C                   D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
C                   COSHZ=D+1.0D0/(3.0D0*D)
C                 ELSE
C                   THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
C                   COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
C                 ENDIF
C                 SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
C                 R=DLOG(COSHZ+SINHZ)
C                 Y=DCOSH(R+DELTA*SE(ns,nb,ne))
C                 ZE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
C               ENDIF
C             ENDIF
              IF(NKT(0,nb).GT.1) THEN
                DO ik=1,NKT(0,nb)
                  DO nn=1,NNT(nb)
                    IF(ns.EQ.(ik+(nn-1)*NKT(0,nb))) THEN
                      nk=ik
                      ns1=1+(nn-1)*NKT(0,nb)
                      GOTO 10
                    ENDIF
                  ENDDO
                ENDDO
              ELSE
                nk=1
              ENDIF
 10           CONTINUE
              IF(nk.EQ.1) THEN
                IF(ITYP11(1).LE.3) THEN
                  R=ZE(ns,1)**(1.0D0/ITYP11(1))
                  ZE(ns,1)=(R+DELTA*SE(ns,nb,ne))**ITYP11(1)
                ELSE IF(ITYP11(1).EQ.4) THEN
                  IF(JTYP10.EQ.2) THEN
                    SS=ZE(ns,1)/(FOCUS*FOCUS)
                    SINHZ=DSQRT(SS)
                    COSHZ=DSQRT(1.0D0+SS)
                    R=DLOG(COSHZ+SINHZ)
                    Y=DSINH(R+DELTA*SE(ns,nb,ne))
                    ZE(ns,1)=(FOCUS*Y)*(FOCUS*Y)
                  ELSE IF(JTYP10.EQ.3) THEN
                    CSS=ZE(ns,1)/FOCUS**3
                    DES=CSS*CSS-4.0D0/27.0D0
                    IF(DES.GT.0.0D0) THEN
                      D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
                      COSHZ=D+1.0D0/(3.0D0*D)
                    ELSE
                      THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
                      COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
                    ENDIF
                    SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
                    R=DLOG(COSHZ+SINHZ)
                    Y=DCOSH(R+DELTA*SE(ns,nb,ne))
                    ZE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
                  ENDIF
                ELSE IF(ITYP11(1).EQ.5) THEN
                ENDIF
              ELSE IF(nk.GT.1) THEN
C****           18-Feb-89: Further modification needed if nodal params
C****           include second deriv, in which case need IDO(nk,nn,0,nb)
                IF(ITYP11(1).LE.3) THEN
                  ZE(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '              *ITYP11(1)*R**(ITYP11(1)-1)
                ELSE IF(ITYP11(1).EQ.4) THEN
                  IF(JTYP10.EQ.2) THEN
                    ZE(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '                *2.0D0*FOCUS*FOCUS*SINHZ*COSHZ
                  ELSE IF(JTYP10.EQ.3) THEN
                    ZE(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '                *FOCUS**3*SINHZ*(3.0D0*COSHZ*COSHZ-1.0D0)
                  ENDIF
                ELSE IF(ITYP11(1).EQ.5) THEN
                ENDIF
              ENDIF

            ELSE
C ***         Normal perturbation for non-isochoric interpolation
              ZE(ns,nhx)=ZE(ns,nhx)+DELTA*SE(ns,nb,ne)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ******* ZE(ns,nhx)='',D12.5,'
     '            //''' Perturbation='',D12.5)')
     '            ZE(ns,nhx),DELTA*SE(ns,nb,ne)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF

C ***       Evaluate perturbed residual
            IF(ITYP1(nr,nx).EQ.3) THEN
              CALL ZERE30(NBH,NBJ,ne,NHE,NJE,
     '          CG,PG,RE1,SE,VE,WG,XE,XG,ZE,ZG,ERROR,*9999)
            ELSE IF(ITYP1(nr,nx).EQ.4) THEN
              CALL ZERE40(NBH,NBJ,ne,NHE,NJE,NW,nx,CE,CG,PG,RE2,SE,
     '          VE,WG,XE,XG,ZE,ZG,ERROR,*9999)
            ELSE IF(ITYP1(nr,nx).EQ.5) THEN
              IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !cnst vol
                CALL ZERE55(INP,NBH,ne,NHE,nr,nx,
     '            CG,PG,RE2,WG,ZE,ZE1,ZG,ERROR,*9999)
              ELSE
                CALL ZERE50(IBT,IDO,INP,NAN,NBH,NBJ,ne,
     '            NFF,NGAP,NHE,NJE,NKF,NNF,NPF,NPNE,nr,NRE,NW,nx,NXI,
     '            CE,CG,CP,FEXT,PG,RE2,RG,SE,VE,WG,
     '            XE,XG,ZE,ZG,ERROR,*9999)
              ENDIF
            ENDIF

            IF(nh.EQ.1.AND.JTYP10.GE.2) THEN
C ***         Isochoric case
C             IF(ITYP10(1).LE.3) THEN
C               ZE(ns,1)=R**ITYP10(1)
C             ELSE IF(ITYP10(1).EQ.4) THEN
C               Y=DCOSH(R)
C               IF(JTYP10.EQ.2) ZE(ns,1)=FOCUS*FOCUS*(Y*Y-1.0D0)
C               IF(JTYP10.EQ.3) ZE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
C             ENDIF
              IF(nk.EQ.1) THEN
                IF(ITYP11(1).LE.3) THEN
                  ZE(ns,1)=R**ITYP11(1)
                ELSE IF(ITYP11(1).EQ.4) THEN
                  Y=DCOSH(R)
                  IF(JTYP10.EQ.2) ZE(ns,1)=FOCUS*FOCUS*(Y*Y-1.0D0)
                  IF(JTYP10.EQ.3) ZE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
                ENDIF
              ELSE IF(nk.EQ.2) THEN
                IF(ITYP11(1).LE.3) THEN
                  ZE(ns,1)=ZE(ns,1)-DELTA*SE(ns,nb,ne)
     '              *ITYP11(1)*R**(ITYP11(1)-1)
                ELSE IF(ITYP11(1).EQ.4) THEN
                  IF(JTYP10.EQ.2) THEN
                    ZE(ns,1)=ZE(ns,1)-DELTA*SE(ns,nb,ne)
     '                *2.0D0*FOCUS*FOCUS*SINHZ*COSHZ
                  ELSE IF(JTYP10.EQ.3) THEN
                    ZE(ns,1)=ZE(ns,1)-DELTA*SE(ns,nb,ne)
     '                *FOCUS**3*SINHZ*(3.0D0*COSHZ*COSHZ-1.0D0)
                  ENDIF
                ENDIF
              ENDIF

            ELSE
C ***         non-isochoric case
              ZE(ns,nhx)=ZE(ns,nhx)-DELTA*SE(ns,nb,ne)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ******* Diagonal term: '','
     '            //'''RE1(ns,nh)='',D20.10,'' RE2(ns,nh)='',D20.10)')
     '            RE1(ns,nh),RE2(ns,nh)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF

C ***       Assemble element stiffness matrix
            mhs=0
            DO mhx=1,NH_LOC(0,nx)
              mh=NH_LOC(mhx,nx)
              DO ms=1,NST(NBH(mh))+NAT(NBH(mh))
                mhs=mhs+1
                my=IABS(LGE(mhs,1)) !local variable number
                my1=GETNYR(1,NPNY,nr,0,2,my,NYNE,NYNP) !global var #
c PJH 19Jul93   IF(.NOT.FIX(my1,1)) THEN
                IF(.NOT.FIX(my1,1)) THEN  !uncommented MPN 14-Sep-94
C ***             Note: DELTA here is global delta
                  ES(mhs,nhs)=(RE2(ms,mh)-RE1(ms,mh))/DELTA
c PJH 19Jul93   ENDIF
                ENDIF  !uncommented MPN 14-Sep-94
              ENDDO
            ENDDO
c PJH 19Jul93 ENDIF
          ENDIF  !uncommented MPN 14-Sep-94
 40     CONTINUE

      ELSE IF(KTYP1D.EQ.2) THEN
C       algebraic calculation of tangent stiffness matrix

        IF(ITYP1(nr,nx).EQ.3) THEN
          ERROR=' Analytic derivs not implemented for PDE''s'
          GO TO 9999
        ELSE IF(ITYP1(nr,nx).EQ.4) THEN
          ERROR=' Analytic derivs not implemented for lin. elasticity'
          GO TO 9999
        ELSE IF(ITYP1(nr,nx).EQ.5) THEN !fe50 problems
          CALL CPCG(1,NBH(1),NPNE,nr,nx,CE,CG,CP,PG,ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
            ERROR='>>Not implemented'
            GO TO 9999
          ELSE
            CALL D_ZERE50('GEOMETRIC_PARAMETERS',IBT,IDO,INP,
     '        NAN,NBH,NBJ,NGAP,ne,NFF,NHE,NJE,NKF,NMNO,
     '        NNF,NPF,NPNE,nr,NRE,NW,nx,NXI,CE,CG,CP,D_RE,D_RI3,
     '        D_TG,D_ZG,ES,FEXT,PG,RG,SE,VE,WG,XE,XG,
     '        ZE,ZE1,ZG,ZG1,ERROR,*9999)
          ENDIF
        ELSE
          ERROR=' Analytic derivs not implemented'
          GO TO 9999
        ENDIF
      ENDIF

 9998 CALL EXITS('ZEES')
      RETURN
 9999 CALL ERRORS('ZEES',ERROR)
      CALL EXITS('ZEES')
      RETURN 1
      END


      SUBROUTINE ENTERS(NAME,*)

C#### Subroutine: ENTERS
C###  Description:
C###    ENTERS traces entry to a subprogram recording the level of
C###    nesting, the entry time, and writing the subprogram name to a
C###    trace file.  Diagnostic o/p is turned on if DIAGNO=.TRUE. and
C###    ALLSUB=.TRUE. or ALLSUB=.FALSE. and NAME=SUBNAM (a subroutine
C###    name).
C**** TRSB         is subprogram name which turns trace on
C**** IOTR         is trace file number (diagnostic o/p file)
C**** TR01         is true if basic tracing on
C**** TR02         is true if full tracing on
C**** NOLV         is current level number
C**** NTLV         is current total number of levels called
C**** NXLV         is maximum number of levels which can be traced
C**** NOSB         is integer number for current subroutine
C**** NTSB         is current total number of subroutines called
C**** NXSB         is maximum no. of subroutines which can be traced
C**** NOSBLV(nolv) is NOSB number at level NOLV
C**** NOSBSM(nosb) is no. of times subroutine no. NOSB has been called
C**** SB(nosb)     is subroutine name for subroutine no. NOSB
C**** TM           is current time
C**** TMST         is trace start time for entry to level 1
C**** TMEL(nolv)   is elapsed time at level NOLV (ie =0 @ entry to NOLV)
C**** TMELSM(nosb) is sum of elapsed times spent in subroutine no. NOSB
C**** TMEN(nolv)   is elapsed time (from start of trace) to level NOLV
C**** TMTLSM(nosb) is total time spent in NOSB & all subs below it

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:ctrl00.cmn'
      INCLUDE 'cmiss$reference:diag00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:trac00.cmn'
!     Parameter List
      CHARACTER NAME*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,ICODE,IEND,IEND1,IEND2,INPUT_CHAN,
     '  ISTATUS,N1SB,NSUB
      REAL*8 TM,VTIME
      CHARACTER CFROMI*10,COLV*10,C1*(MXCH),ERROR*10
c      EXTERNAL CTRLC_AST

      IF(DIAGNO) THEN
        CALL STRING_TRIM(NAME,IBEG1,IEND1)
        IF(ALLSUB) THEN !turn diagnostics on in all subroutines
          WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE IF(.NOT.ALLSUB) THEN !diags on in selected subrs
          NUM_STACK=NUM_STACK+1
          DOP_STACK(NUM_STACK)=.FALSE.
          DO NSUB=1,NT_SUB
            CALL STRING_TRIM(SUBNAM(NSUB),IBEG2,IEND2)
            CALL CUPPER(NAME(IBEG1:IEND1),C1,ERROR,*9999)
          IF(C1.EQ.SUBNAM(NSUB)(IBEG2:IEND2))
C          IF(CUPPER(NAME(IBEG1:IEND1)).EQ.SUBNAM(NSUB)(IBEG2:IEND2))
     '        DOP_STACK(NUM_STACK)=.TRUE.
          ENDDO
          IF(FROMSUB) THEN
            IF(DOP_STACK(NUM_STACK-1)) DOP_STACK(NUM_STACK)=.TRUE.
          ENDIF
          IF(DOP_STACK(NUM_STACK)) THEN
            WRITE(OP_STRING,'('' *** Enter '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(DOP_STACK(NUM_STACK-1)) THEN
            WRITE(OP_STRING,'('' *** Calls '',A)') NAME(IBEG1:IEND1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          DOP=DOP_STACK(NUM_STACK)
        ENDIF
      ENDIF

      IF(TR01) THEN
        IF(NAME.EQ.TRSB) THEN
          TR02=.TRUE.
        ENDIF
        TR01=.FALSE.
        NOLV=NOLV+1
        IF(NOLV.GT.NTLV) THEN
          NTLV=NOLV
        ENDIF
        IF((NOLV.GT.0).AND.(NOLV.LE.NXLV)) THEN
          TM=VTIME()
          IF(NOLV.EQ.1) THEN
            TMST=TM
            IF(TR02) THEN
              WRITE(OP_STRING,'(/''      Time:    Calls:    Level:'//
     '          ' >Subprogram entered'''//
     '          '/''      Time:    Total:   Actual:'//
     '          ' <Subprogram exited'')')
              CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          TMEN(NOLV)=TM-TMST
          TMEL(NOLV)=0.0D0
          DO N1SB=1,NTSB
            IF(SB(N1SB).EQ.NAME) THEN
              NOSB=N1SB
              GOTO 2
            ENDIF
          ENDDO
          IF(NTSB.LT.NXSB) THEN
            NTSB=NTSB+1
            SB(NTSB)=NAME
            NOSB=NTSB
          ELSE
            NTSB=NXSB+1
            NOSB=NXSB+1
          ENDIF
    2     NOSBLV(NOLV)=NOSB
          NOSBSM(NOSB)=NOSBSM(NOSB)+1
          IF(TR02) THEN
            COLV=CFROMI(NOLV,'(I10)')
            CALL STRING_TRIM(COLV,IBEG,IEND)
            WRITE(OP_STRING,'(1X,F10.3,I10,I10,'//COLV(IBEG:IEND)//
     '        '('' >''),A)')
     '        TMEN(NOLV),NOSBSM(NOSB),NOLV,NAME
            CALL WRITES(IOTR,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        TR01=.TRUE.
      ENDIF

      RETURN
 9999 RETURN 1
      END


      SUBROUTINE GET_COMMAND_LINE(ARGS,NUMARGS)

C#### Subroutine: GET_COMMAND_LINE
C###  Description:
C###    GET_COMMAND_LINE gets the command line arguments.

      IMPLICIT NONE
!     Parameter List
      INTEGER NUMARGS
      CHARACTER ARGS(10)*80
!     Local Variables
      INTEGER CLOCAT,I,IBEG,IEND,LENGTH,OLDPOSN,POSN,IARGC
      CHARACTER COMMAND_LINE*80

      NUMARGS = IARGC() !SGI Fortran routine
      DO I=1,NUMARGS
        CALL GETARG(I,ARGS(I)) !SGI Fortran routine
      ENDDO
!old  CALL LIB$GET_FOREIGN(COMMAND_LINE,,LENGTH,)
!old  CALL STRING_TRIM(COMMAND_LINE,IBEG,IEND)
!old  I=0
!old  OLDPOSN=IBEG
!old  POSN=CLOCAT(' ',COMMAND_LINE(OLDPOSN:))
!old  DO WHILE(POSN.NE.1)
!old    I=I+1
!old    ARGS(I)=COMMAND_LINE(OLDPOSN:POSN-1)
!old    OLDPOSN=POSN+1
!old    POSN=CLOCAT(' ',COMMAND_LINE(OLDPOSN:))
!old    IF(I.EQ.10) POSN=1 !No more than 10 args at present
!old  ENDDO
!old  NUMARGS=I

      RETURN
      END


      REAL*8 FUNCTION AG(nb,nh,nr,ns,PPG,TG,ZG)

C#### Function: AG
C###  Type: REAL*8
C###  Description:
C###    AG evaluates integrand of the domain integral in the virtual
C###    work equation for all 4 coordinate systems and all 8 equation
C###    types.  Note: PPG(ns,1+nj) has derivative of basis function
C###    wrt x(nj), where x coord type is as defined by DXIX.
C**** 19Sep88: Changed IF(JTYP9.EQ.0) to IF(KTYP53(nr).EQ.1) and
C****          changed IF(JTYP9.EQ.1) to IF(KTYP53(nr).GE.2) since
C****          KTYP53(nr)now defines whether stresses in the
C****          constitutive laware referred to theta (KTYP53(nr)=1)
C****          or Nu (KTYP53(nr)>1).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
!     Parameter List
      INTEGER nb,nh,nr,ns
      REAL*8 PPG(64,4),TG(3,3),ZG(NHM,NUM)
!     Local Variables
      INTEGER mi,ni,NITB,NU1(0:3)
      REAL*8 CLZ,CMZ,CSLZ,CSMZ,CZ,
     '  D1,D2,D3,DLA,DMA,DRA,DTA,E1,E2,G1,
     '  PA,PGA,PGG,RA,RZ,SLZ,SMZ,SZ,TA
      DATA NU1/1,2,4,7/

      AG=0.0d0
      NITB=NIT(nb)
      IF(ITYP10(nr).EQ.1) THEN       !rectangular cartesian coords
        DO mi=1,NITB
          DO ni=1,NITB
            AG=AG+TG(mi,ni)*ZG(nh,NU1(ni))*PPG(ns,1+mi)
          ENDDO
        ENDDO

      ELSE IF(ITYP10(nr).EQ.2) THEN  !cylindrical polar coords
        RZ=ZG(1,1)
        PGG=PPG(ns,1)
        DO ni=1,NITB
          D1=ZG(1,NU1(ni))
          D2=ZG(2,NU1(ni))
          D3=ZG(3,NU1(ni))
          DO mi=1,NITB
            DRA=ZG(1,NU1(mi)) !check this shells etc (KTYP51(nr).GT.3)
            DTA=ZG(2,NU1(mi)) !  "     "   "     "    "     "     "
            PGA=PPG(ns,1+mi)
            IF(nh.EQ.1) THEN
              AG=AG+TG(mi,ni)*(D1*PGA+D2*RZ*DTA*PGG)
            ELSE IF(nh.EQ.2) THEN
              AG=AG+TG(mi,ni)*(-D1*DTA/RZ*PGG+D2*(PGA-DRA/RZ*PGG))
            ELSE IF(nh.EQ.3) THEN
              AG=AG+TG(mi,ni)*D3*PGA
            ENDIF
          ENDDO
        ENDDO

      ELSE IF(ITYP10(nr).EQ.3) THEN  !spherical polar coords
        RZ=ZG(1,1)
        CZ=DCOS(ZG(3,1))
        SZ=DSIN(ZG(3,1))
        PGG=PPG(ns,1)
        DO ni=1,NITB
          D1=ZG(1,NU1(ni))
          D2=ZG(2,NU1(ni))
          D3=ZG(3,NU1(ni))
          DO mi=1,NITB
            RA=ZG(1,NU1(mi))
            TA=ZG(2,NU1(mi))
            PA=ZG(3,NU1(mi))
            PGA=PPG(ns,1+mi)
            IF(nh.EQ.1) THEN
              AG=AG+TG(mi,ni)*(D1*PGA+RZ*(D2*CZ*CZ*TA+D3*PA)*PGG)
            ELSE IF(nh.EQ.2) THEN
              AG=AG+TG(mi,ni)*(TA*(D3*SZ/CZ-D1/RZ)*PGG
     '                        +D2*(PGA-(RA/RZ-SZ/CZ*PA)*PGG))
            ELSE IF(nh.EQ.3) THEN
              AG=AG+TG(mi,ni)*(D3*(PGA-RA/RZ*PGG)
     '                       -(D1*PA/RZ+D2*CZ*SZ*TA)*PGG)
            ENDIF
          ENDDO
        ENDDO

      ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate spheroidal coords
        SLZ=DSINH(ZG(1,1))
        SMZ=DSIN (ZG(2,1))
        CLZ=DSQRT(1.0d0+SLZ*SLZ)
        CMZ=DSQRT(1.0d0-SMZ*SMZ)
        CSLZ=CLZ/SLZ
        CSMZ=CMZ/SMZ
        G1=SLZ*SLZ+SMZ*SMZ
        E1=CLZ*SLZ/G1
        E2=CMZ*SMZ/G1
        PGG=PPG(ns,1)
        DO ni=1,NITB
          D1=ZG(1,NU1(ni))
          D2=ZG(2,NU1(ni))
          D3=ZG(3,NU1(ni))
          DO mi=1,NITB
            DLA=ZG(1,NU1(mi))
            DMA=ZG(2,NU1(mi))
            DTA=ZG(3,NU1(mi))
            PGA=PPG(ns,1+mi)
            IF(nh.EQ.1) THEN
              AG=AG+TG(mi,ni)*(D1*(PGA-(E1*DLA+E2*DMA)*PGG)
     '             +D2*(E1*DMA-E2*DLA)*PGG+D3*E1*SMZ*SMZ*DTA*PGG)
            ELSE IF(nh.EQ.2) THEN
              AG=AG+TG(mi,ni)*(D1*(E2*DLA-E1*DMA)*PGG
     '             +D2*(PGA-(E1*DLA+E2*DMA)*PGG)+D3*E2*SLZ*SLZ*DTA*PGG)
            ELSE IF(nh.EQ.3) THEN
              AG=AG+TG(mi,ni)*(-(D1*CSLZ+D2*CSMZ)*DTA*PGG
     '             +D3*(PGA-(CSLZ*DLA+CSMZ*DMA)*PGG))
            ENDIF
          ENDDO
        ENDDO
      ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate spheroidal coords
      ENDIF

      RETURN
      END


      REAL*8 FUNCTION D_AG(PARAMTYPE,nb,nh,nr,ns,
     '  D_TG,D_ZG,PPG,TG,ZG)

C#### Function: D_AG
C###  Type: REAL*8
C###  Description:
C###    D_AG evaluates the derivative wrt PARAMTYPE parameters, of the
C###    integrand of the domain integral in the virtual work
C###    equation for all 4 coordinate systems and all 8 equation types.
C###    Note: PPG(ns,1+nj) has derivative of basis function wrt x(nj),
C###    where x coord type is as defined by DXIX.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
!     Parameter List
      INTEGER nb,nh,nr,ns
      REAL*8 D_TG(3,3),D_ZG(NHM,NUM),PPG(64,4),TG(3,3),ZG(NHM,NUM)
      CHARACTER PARAMTYPE*(*)
!     Local Variables
      INTEGER mi,ni,NITB,NU1(0:3)
      REAL*8 AG,CLZ,CMZ,CSLZ,CSMZ,CZ,D1,D2,D3,D_CLZ,D_CMZ,D_CSLZ,
     '  D_CSMZ,D_CZ,D_D1,D_D2,D_D3,D_DLA,D_DMA,D_DRA,D_DTA,D_E1,
     '  D_E2,D_G1,D_PA,D_RA,D_RZ,D_SLZ,D_SMZ,D_SZ,D_TA,DLA,DMA,DRA,
     '  DTA,E1,E2,G1,PA,PGA,PGG,RA,RZ,SLZ,SMZ,SZ,TA
      DATA NU1/1,2,4,7/

      D_AG=0.0d0
      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        D_AG=AG(nb,nh,nr,ns,PPG,D_TG,ZG)

      ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
        NITB=NIT(nb)

        IF(ITYP10(nr).EQ.1) THEN       !rectangular cartesian coords
          DO mi=1,NITB
            DO ni=1,NITB
              D_AG=D_AG+(D_TG(mi,ni)*ZG(nh,NU1(ni))
     '                  +TG(mi,ni)*D_ZG(nh,NU1(ni)))*PPG(ns,1+mi)
            ENDDO
          ENDDO

        ELSE IF(ITYP10(nr).EQ.2) THEN  !cylindrical polar coords
          RZ  =  ZG(1,1)
          D_RZ=D_ZG(1,1)
          PGG=PPG(ns,1)
          DO ni=1,NITB
            D1  =  ZG(1,NU1(ni))
            D_D1=D_ZG(1,NU1(ni))
            D2  =  ZG(2,NU1(ni))
            D_D2=D_ZG(2,NU1(ni))
            D3  =  ZG(3,NU1(ni))
            D_D3=D_ZG(3,NU1(ni))
            DO mi=1,NITB
              DRA  =  ZG(1,NU1(mi)) !check shells etc(KTYP51(nr).GT.3)
              D_DRA=D_ZG(1,NU1(mi))
              DTA  =  ZG(2,NU1(mi)) !  "    "     "    "     "     "
              D_DTA=D_ZG(2,NU1(mi))
              PGA=PPG(ns,1+mi)
              IF(nh.EQ.1) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*PGA+D2*RZ*DTA*PGG)
     '                   +TG(mi,ni)*(D_D1*PGA
     '                     +(D_D2*RZ*DTA+D2*D_RZ*DTA+D2*RZ*D_DTA)*PGG)
              ELSE IF(nh.EQ.2) THEN
                D_AG=D_AG+D_TG(mi,ni)*(-D1*DTA/RZ*PGG
     '                                 +D2*(PGA-DRA/RZ*PGG))
     '                   +TG(mi,ni)*(D_D2*PGA
     '                            +(-D_D1*DTA/RZ
     '                              -D1*D_DTA/RZ
     '                              +D1*DTA*D_RZ/(RZ*RZ)
     '                              -D_D2*DRA/RZ
     '                              -D2*D_DRA/RZ
     '                              +D2*DRA*D_RZ/(RZ*RZ))*PGG)
              ELSE IF(nh.EQ.3) THEN
                D_AG=D_AG+(D_TG(mi,ni)*D3+TG(mi,ni)*D_D3)*PGA
              ENDIF
            ENDDO
          ENDDO

        ELSE IF(ITYP10(nr).EQ.3) THEN  !spherical polar coords
          RZ  =  ZG(1,1)
          D_RZ=D_ZG(1,1)
          CZ  = DCOS(ZG(3,1))
          D_CZ=-DSIN(ZG(3,1))*D_ZG(3,1)
          SZ  = DSIN(ZG(3,1))
          D_SZ= DCOS(ZG(3,1))*D_ZG(3,1)
          PGG=PPG(ns,1)
          DO ni=1,NITB
            D1  =  ZG(1,NU1(ni))
            D_D1=D_ZG(1,NU1(ni))
            D2  =  ZG(2,NU1(ni))
            D_D2=D_ZG(2,NU1(ni))
            D3  =  ZG(3,NU1(ni))
            D_D3=D_ZG(3,NU1(ni))
            DO mi=1,NITB
              RA  =  ZG(1,NU1(mi))
              D_RA=D_ZG(1,NU1(mi))
              TA  =  ZG(2,NU1(mi))
              D_TA=D_ZG(2,NU1(mi))
              PA  =  ZG(3,NU1(mi))
              D_PA=D_ZG(3,NU1(mi))
              PGA=PPG(ns,1+mi)
              IF(nh.EQ.1) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*PGA+
     '                                 RZ*(D2*CZ*CZ*TA+D3*PA)*PGG)
     '                   +TG(mi,ni)*(D_D1*PGA
     '                             +(D_RZ*D2*CZ*CZ*TA
     '                              +RZ*D_D2*CZ*CZ*TA
     '                              +RZ*D2*2.0d0*CZ*D_CZ*TA
     '                              +RZ*D2*CZ*CZ*D_TA
     '                              +D_RZ*D3*PA
     '                              +RZ*D_D3*PA
     '                              +RZ*D3*D_PA)*PGG)
              ELSE IF(nh.EQ.2) THEN
                D_AG=D_AG+D_TG(mi,ni)*(TA*(D3*SZ/CZ-D1/RZ)*PGG
     '                         +D2*(PGA-(RA/RZ-SZ/CZ*PA)*PGG))
     '                   +TG(mi,ni)*(D_D2*PGA
     '                             +(D_TA*D3*SZ/CZ
     '                              +TA*D_D3*SZ/CZ
     '                              +TA*D3*D_SZ/CZ
     '                              -TA*D3*SZ*D_CZ/(CZ*CZ)
     '                              -D_TA*D1/RZ
     '                              -TA*D_D1/RZ
     '                              +TA*D1*D_RZ/(RZ*RZ)
     '                              -D_D2*RA/RZ
     '                              -D2*D_RA/RZ
     '                              +D2*RA*D_RZ/(RZ*RZ)
     '                              +D_D2*SZ/CZ*PA
     '                              +D2*D_SZ/CZ*PA
     '                              -D2*SZ*D_CZ/(CZ*CZ)*PA
     '                              +D2*SZ/CZ*D_PA)*PGG)
              ELSE IF(nh.EQ.3) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D3*(PGA-RA/RZ*PGG)
     '                         -(D1*PA/RZ+D2*CZ*SZ*TA)*PGG)
     '                   +TG(mi,ni)*(D_D3*PGA
     '                            +(-D_D3*RA/RZ
     '                              -D3*D_RA/RZ
     '                              +D3*RA*D_RZ/(RZ*RZ)
     '                              -D_D1*PA/RZ
     '                              -D1*D_PA/RZ
     '                              +D1*PA*D_RZ/(RZ*RZ)
     '                              -D_D2*CZ*SZ*TA
     '                              -D2*D_CZ*SZ*TA
     '                              -D2*CZ*D_SZ*TA
     '                              -D2*CZ*SZ*D_TA)*PGG)
              ENDIF
            ENDDO
          ENDDO

        ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate spheroidal coords
          SLZ  =DSINH(ZG(1,1))
          D_SLZ=DCOSH(ZG(1,1))*D_ZG(1,1)
          SMZ  =DSIN (ZG(2,1))
          D_SMZ=DCOS (ZG(2,1))*D_ZG(2,1)
          CLZ  = DSQRT(1.0d0+SLZ*SLZ)
          D_CLZ= SLZ*D_SLZ/DSQRT(1.0d0+SLZ*SLZ)
          CMZ  = DSQRT(1.0d0-SMZ*SMZ)
          D_CMZ=-SMZ*D_SMZ/DSQRT(1.0d0-SMZ*SMZ)
          CSLZ  =CLZ/SLZ
          D_CSLZ=D_CLZ/SLZ-CLZ*D_SLZ/(SLZ*SLZ)
          CSMZ  =CMZ/SMZ
          D_CSMZ=D_CMZ/SMZ-CMZ*D_SMZ/(SMZ*SMZ)
          G1  =SLZ*SLZ+SMZ*SMZ
          D_G1=2.0d0*(SLZ*D_SLZ+SMZ*D_SMZ)
          E1  =CLZ*SLZ/G1
          D_E1=D_CLZ*SLZ/G1+CLZ*D_SLZ/G1-CLZ*SLZ*D_G1/(G1*G1)
          E2  =CMZ*SMZ/G1
          D_E2=D_CMZ*SMZ/G1+CMZ*D_SMZ/G1-CMZ*SMZ*D_G1/(G1*G1)
          PGG=PPG(ns,1)
          DO ni=1,NITB
            D1  =  ZG(1,NU1(ni))
            D_D1=D_ZG(1,NU1(ni))
            D2  =  ZG(2,NU1(ni))
            D_D2=D_ZG(2,NU1(ni))
            D3  =  ZG(3,NU1(ni))
            D_D3=D_ZG(3,NU1(ni))
            DO mi=1,NITB
              DLA  =  ZG(1,NU1(mi))
              D_DLA=D_ZG(1,NU1(mi))
              DMA  =  ZG(2,NU1(mi))
              D_DMA=D_ZG(2,NU1(mi))
              DTA  =  ZG(3,NU1(mi))
              D_DTA=D_ZG(3,NU1(mi))
              PGA=PPG(ns,1+mi)
              IF(nh.EQ.1) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*(PGA-(E1*DLA+E2*DMA)*PGG)
     '             +D2*(E1*DMA-E2*DLA)*PGG+D3*E1*SMZ*SMZ*DTA*PGG)
     '                   +TG(mi,ni)*(D_D1*PGA
     '                            +(-D_D1*E1*DLA-D1*D_E1*DLA-D1*E1*D_DLA
     '                              -D_D1*E2*DMA-D1*D_E2*DMA-D1*E2*D_DMA
     '                              +D_D2*E1*DMA+D2*D_E1*DMA+D2*E1*D_DMA
     '                              -D_D2*E2*DLA-D2*D_E2*DLA-D2*E2*D_DLA
     '                              +D_D3*E1*SMZ*SMZ*DTA
     '                              +D3*D_E1*SMZ*SMZ*DTA
     '                              +D3*E1*2.0d0*SMZ*D_SMZ*DTA
     '                              +D3*E1*SMZ*SMZ*D_DTA)*PGG)
              ELSE IF(nh.EQ.2) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*(E2*DLA-E1*DMA)*PGG
     '             +D2*(PGA-(E1*DLA+E2*DMA)*PGG)+D3*E2*SLZ*SLZ*DTA*PGG)
     '                   +TG(mi,ni)*(D_D2*PGA
     '                             +(D_D1*E2*DLA+D1*D_E2*DLA+D1*E2*D_DLA
     '                              -D_D1*E1*DMA-D1*D_E1*DMA-D1*E1*D_DMA
     '                              -D_D2*E1*DLA-D2*D_E1*DLA-D2*E1*D_DLA
     '                              -D_D2*E2*DMA-D2*D_E2*DMA-D2*E2*D_DMA
     '                              +D_D3*E2*SLZ*SLZ*DTA
     '                              +D3*D_E2*SLZ*SLZ*DTA
     '                              +D3*E2*2.0d0*SLZ*D_SLZ*DTA
     '                              +D3*E2*SLZ*SLZ*D_DTA)*PGG)
              ELSE IF(nh.EQ.3) THEN
                D_AG=D_AG+D_TG(mi,ni)*(-(D1*CSLZ+D2*CSMZ)*DTA*PGG
     '             +D3*(PGA-(CSLZ*DLA+CSMZ*DMA)*PGG))
     '                   +TG(mi,ni)*(D_D3*PGA
     '                +(-D_D1*CSLZ*DTA-D1*D_CSLZ*DTA-D1*CSLZ*D_DTA
     '                  -D_D2*CSMZ*DTA-D2*D_CSMZ*DTA-D2*CSMZ*D_DTA
     '                  -D_D3*CSLZ*DLA-D3*D_CSLZ*DLA-D3*CSLZ*D_DLA
     '                  -D_D3*CSMZ*DMA-D3*D_CSMZ*DMA-D3*CSMZ*D_DMA)*PGG)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate spheroidal coords
        ENDIF
      ENDIF

      RETURN
      END


      SUBROUTINE D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,P1,P2,P3,P4,P5,P6,
     '  ERROR,*)

C#### Subroutine: D_ENERGY
C###  Description:
C###    D_ENERGY calculates 2nd derivatives of strain energy function
C###    wrt either the NMNO'th material parameter or the
C###    physical strains (KTYP55(nr)=3) at current Gauss point

C**** P1..P6 are the physical components of Green's strain wrt nu
C**** coordinates
C**** (KTYP53(nr)>1):E(1,1),E(2,2),E(3,3),E(1,2),E(1,3) and E(2,3)
C**** D_DW(1..6) are d2W/dE2(1,1) .. d2W/dE2(2,3)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER NMNO,nr
      REAL*8 CG(NMM),D_DW(6),P1,P2,P3,P4,P5,P6
      CHARACTER ERROR*(*),PARAMTYPE*(*)
!     Local Variables
      INTEGER iterm,iterm1,iterm2,pole1num,pole2num
      REAL*8 alpha,alpha1,alpha2,coeff,coeff1,coeff2,
     '  D_DWslope,D_DWslope1,D_DWslope2,DENOM,DENOM1,DENOM2,
     '  Dpole,DtoAlpha,DtoAlpha1,DtoAlpha2,pole,pole1,pole2,
     '  strain,strain1,strain11,strain2,strain21,TOL,ZERO
      PARAMETER(ZERO=0.0d0)

      CALL ENTERS('D_ENERGY',*9999)

      CALL ASSERT(KTYP55(nr).EQ.3,'Strain energy fn must be given as a'
     '  //' function of fibre and transverse strains',ERROR,*9999)
      CALL ASSERT(KTYP56(nr).EQ.3,'Strain energy fn must be'
     '  //' pole-zero in fibre and transverse strains',ERROR,*9999)

C old
C      L0_fibre=CG(28)           !initial fibre ext ratio
C      L0_sheet=CG(29)           !initial sheet ext ratio
C      L0_sheetnormal=CG(30)    !initial sheet-normal ext ratio
C      E0_fibre=0.5d0*(L0_fibre*L0_fibre-1.0d0) !init fibre strain
C      E0_sheet=0.5d0*(L0_sheet*L0_sheet-1.0d0) !init sheet strain
C      E0_sheetnormal=0.5d0*(L0_sheetnormal*L0_sheetnormal-1.0d0) !init sheet-normal strain
C      PP1=P1+E0_fibre
C      PP2=P2+E0_sheet
C      PP3=P3+E0_sheetnormal

      DO iterm=1,6
        D_DW(iterm)=0.0d0
      ENDDO

C      Calculate the 2nd deriv of the pole-zero law
C      for strains between 0 and 90% of the pole for each term.
C      For strains beyond 90% of their pole the stress/strain is flat
C      with magnitude equal to the stress at a strain of 90% of the pole
C      For compressive strains use linear stress/strain relationship
C      with stiffness equal to the 2nd deriv of W wrt e
C      evaluated at e=0
      TOL=0.90d0

      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN

        CALL ASSERT(NMNO.LE.9,' >>Analytic derivs only defined'
     '    //' wrt material param no.s 1-9',ERROR,*9999)

C       Determine term in constitutive law the depends on material param
C       to be optimised
        iterm=(NMNO-1)/3+1    !NOTE: need integer division here
        IF(iterm.EQ.1) THEN
          strain=P1 !=EG(1,1)
        ELSE IF(iterm.EQ.2) THEN
          strain=P2 !=EG(2,2)
        ELSE IF(iterm.EQ.3) THEN
          strain=P3 !=EG(3,3)
        ENDIF

        IF(NMNO.EQ.1.OR.NMNO.EQ.4.OR.NMNO.EQ.7) THEN
C         Parameter is a coefficient.
C !!!     NOTE: if fibre distribution model is altered to express shear
C !!!     coefficients in terms of axial coefficients then this will
C !!!     need updating
          coeff=CG(NMNO)
          pole =CG(NMNO+1)
          alpha=CG(NMNO+2)
          strain1=strain
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
          ELSE IF(strain.GT.(TOL*pole)) THEN
            strain1=TOL*pole !yielded
          ENDIF
          DENOM=pole-strain1
          DtoAlpha=1.0d0
          IF(alpha.GT.ZERO) DtoAlpha=DENOM**alpha
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole)) THEN
C           Calc slope
            D_DWslope=(2.0d0*DENOM*(DENOM+2.0d0*alpha*strain1)
     '        +alpha*(alpha+1.0d0)*strain1*strain1)
     '        /(DtoAlpha*DENOM*DENOM)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=strain1*(2.0d0+alpha*strain1/DENOM)/DtoAlpha
            IF(strain.GT.(TOL*pole))
     '        D_DW(iterm)=D_DW(iterm)+(strain-strain1)*D_DWslope
          ELSE !Compressive strain
C           Use deriv of DW at strain=0 for compresive slope
            D_DW(iterm)=strain*D_DWslope
          ENDIF

        ELSE IF(NMNO.EQ.2.OR.NMNO.EQ.5.OR.NMNO.EQ.8) THEN
C         Parameter is a pole - set params for main term in
C         constitutive law
C !!!     NOTE: if fibre distribution model is altered for the shear
C !!!     poles then this will need updating
          coeff=CG(NMNO-1)
          pole =CG(NMNO)
          alpha=CG(NMNO+1)

C         Calculate deriv wrt axial pole
          strain1=strain
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
          ELSE IF(strain.GT.(TOL*pole)) THEN
            strain1=TOL*pole !yielded
          ENDIF
          DENOM=pole-strain1
          DtoAlpha=1.0d0
          IF(alpha.GT.ZERO) DtoAlpha=DENOM**alpha
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole)) THEN
C           Calc slope
            D_DWslope=(4.0d0*coeff*DENOM*(DENOM+alpha*strain1)
     '        -coeff*(alpha+2.0d0)*(2.0d0*DENOM*(DENOM
     '        +2.0d0*alpha*strain1)
     '        +alpha*(alpha+1.0d0)*strain1*strain1))
     '        /(DtoAlpha*DENOM*DENOM*DENOM)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=-alpha*coeff*strain1
     '        *(2.0d0+(alpha+1.0d0)*strain1/DENOM)/(DtoAlpha*DENOM)
            IF(strain1.GT.(TOL*pole))
     '        D_DW(iterm)=D_DW(iterm)+(strain-strain1)*D_DWslope
          ELSE !Compressive strain
            D_DW(iterm)=strain*D_DWslope
          ENDIF

C         Set params for shear terms in constitutive law.
          IF(NMNO.EQ.2) THEN !1 dirn
C           params associated with 1-2 deformation mode
            strain1=P4 !EG(1,2)
            iterm1=4
            pole1num=11
C           params associated with 1-3 deformation mode
            strain2=P5 !EG(1,3)
            iterm2=5
            pole2num=14
          ELSE IF(NMNO.EQ.5) THEN  !2 dirn
C           params associated with 2-1 deformation mode
            strain1=P4 !EG(2,1)=EG(1,2)
            iterm1=4
            pole1num=17
C           params associated with 2-3 deformation mode
            strain2=P6 !EG(2,3)
            iterm2=6
            pole2num=20
          ELSE IF(NMNO.EQ.8) THEN  !3 dirn
C           params associated with 3-1 deformation mode
            strain1=P5 !EG(3,1)=EG(1,3)
            iterm1=5
            pole1num=23
C           params associated with 3-2 deformation mode
            strain2=P6 !EG(3,2)=EG(2,3)
            iterm2=6
            pole2num=26
          ENDIF

C         Set parameters appropriate for selected terms of the
C         constitutive law
          coeff1=CG(pole1num-1)
          pole1 =CG(pole1num)
          alpha1=CG(pole1num+1)
          coeff2=CG(pole2num-1)
          pole2 =CG(pole2num)
          alpha2=CG(pole2num+1)

C         Calc deriv of shear poles wrt axial pole to be optimised
          Dpole=2.0d0*(1.0d0+pole)/((1.0d0+2.0d0*pole)**1.5d0)

C         Calc d(DW(iterm1))/d(pole) using chain rule. ie)
C         d(DW(iterm1))/d(pole)=d(DW(iterm1))/d(pole1)*d(pole1)/d(pole)
          strain11=strain1
          IF(strain1.LE.ZERO) THEN !compressive strain
            strain11=ZERO
          ELSE IF(strain1.GT.(TOL*pole1)) THEN
            strain11=TOL*pole1 !yielded
          ENDIF
          DENOM1=pole1-strain11
          DtoAlpha1=1.0d0
          IF(alpha1.GT.ZERO) DtoAlpha1=DENOM1**alpha1
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain1.LE.ZERO.OR.strain1.GT.(TOL*pole1)) THEN
C           Calc slope
            D_DWslope1=Dpole*(4.0d0*coeff1*DENOM1*(DENOM1
     '        +alpha1*strain11)
     '        -coeff1*(alpha1+2.0d0)*(2.0d0*DENOM1*(DENOM1
     '        +2.0d0*alpha1*strain11)
     '        +alpha1*(alpha1+1.0d0)*strain11*strain11))
     '        /(DtoAlpha1*DENOM1*DENOM1*DENOM1)
          ENDIF
          IF(strain1.GT.ZERO) THEN !Tensile strain
            D_DW(iterm1)=0.5d0*Dpole*(-alpha1*coeff1*strain11
     '        *(2.0d0+(alpha1+1.d0)*strain11/DENOM1)/(DtoAlpha1*DENOM1))
            IF(strain1.GT.(TOL*pole1)) D_DW(iterm1)=D_DW(iterm1)
     '        +(strain1-strain11)*D_DWslope1/2.0d0
          ELSE !Compressive strain
            D_DW(iterm1)=strain1*D_DWslope1/2.0d0
          ENDIF

C         Calc d(DW(iterm2))/d(pole) using chain rule. ie)
C         d(DW(iterm2))/d(pole)=d(DW(iterm2))/d(pole2)*d(pole2)/d(pole)
          strain21=strain2
          IF(strain2.LE.ZERO) THEN !compressive strain
            strain21=ZERO
          ELSE IF(strain2.GT.(TOL*pole2)) THEN
            strain21=TOL*pole2 !yielded
          ENDIF
          DENOM2=pole2-strain21
          DtoAlpha2=1.0d0
          IF(alpha2.GT.ZERO) DtoAlpha2=DENOM2**alpha2
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain2.LE.ZERO.OR.strain2.GT.(TOL*pole2)) THEN
C           Calc slope
            D_DWslope2=Dpole*(4.0d0*coeff2*DENOM2*(DENOM2
     '        +alpha2*strain21)
     '        -coeff2*(alpha2+2.0d0)*(2.0d0*DENOM2*(DENOM2
     '        +2.0d0*alpha2*strain21)
     '        +alpha2*(alpha2+1.0d0)*strain21*strain21))
     '        /(DtoAlpha2*DENOM2*DENOM2*DENOM2)
          ENDIF
          IF(strain2.GT.ZERO) THEN !Tensile strain
            D_DW(iterm2)=0.5d0*Dpole*(-alpha2*coeff2*strain21
     '        *(2.0d0+(alpha2+1.d0)*strain21/DENOM2)/(DtoAlpha2*DENOM2))
            IF(strain2.GT.(TOL*pole2)) D_DW(iterm2)=D_DW(iterm2)
     '        +(strain2-strain21)*D_DWslope2/2.0d0
          ELSE !Compressive strain
            D_DW(iterm2)=strain2*D_DWslope2/2.0d0
          ENDIF

        ELSE IF(NMNO.EQ.3.OR.NMNO.EQ.6.OR.NMNO.EQ.9) THEN
C         Parameter is a curvature.
C !!!     NOTE: if fibre distribution model is altered to express shear
C !!!     curvatures in terms of axial curvatures then this will
C !!!     need updating
          coeff=CG(NMNO-2)
          pole =CG(NMNO-1)
          alpha=CG(NMNO)
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
          ELSE IF(strain.GT.(TOL*pole)) THEN
            strain1=TOL*pole !yielded
          ENDIF
          DENOM=pole-strain1
          DtoAlpha=1.0d0
          IF(alpha.GT.ZERO) DtoAlpha=DENOM**alpha
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole)) THEN
C           Calc slope
            D_DWslope=(coeff*strain1*(4.0d0*DENOM
     '        +strain1*(2.0d0*alpha+1.0d0))
     '        +coeff*(2.0d0*DENOM*(DENOM+2.0d0*alpha*strain1)
     '        +alpha*(alpha+1.0d0)*strain1*strain1)*DLOG(DENOM))
     '        /(DtoAlpha*DENOM*DENOM)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=coeff*strain1*(strain1/DENOM
     '        -(2.0d0+alpha*strain1/DENOM)*DLOG(DENOM))/DtoAlpha
            IF(strain.GT.(TOL*pole))
     '        D_DW(iterm)=D_DW(iterm)+(strain-strain1)*D_DWslope
          ELSE !Compressive strain
            D_DW(iterm)=strain*D_DWslope
          ENDIF
        ENDIF

      ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN

        DO iterm=1,6
          coeff1=CG((iterm-1)*3+1)
          pole1 =CG((iterm-1)*3+2)
          alpha1=CG((iterm-1)*3+3)
          IF(iterm.EQ.1.OR.iterm.EQ.2.OR.iterm.EQ.3) THEN  !Axial
            IF(iterm.EQ.1) THEN
              strain=P1
            ELSE IF(iterm.EQ.2) THEN
              strain=P2
            ELSE IF(iterm.EQ.3) THEN
              strain=P3
            ENDIF
            coeff2=coeff1
            pole2 =pole1
            alpha2=alpha1
          ELSE IF(iterm.EQ.4.OR.iterm.EQ.5.OR.iterm.EQ.6) THEN  !Shear
            IF(iterm.EQ.4) THEN
              strain=P4 !=EG(1,2)
            ELSE IF(iterm.EQ.5) THEN
              strain=P5 !=EG(1,3)
            ELSE IF(iterm.EQ.6) THEN
              strain=P6 !=EG(2,3)
            ENDIF
            coeff2=CG((iterm-1)*3+10)
            pole2 =CG((iterm-1)*3+11)
            alpha2=CG((iterm-1)*3+12)
          ENDIF

          strain1=strain
          strain2=strain
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
            strain2=ZERO
          ELSE IF(strain.GT.(TOL*pole1).OR.strain.GT.(TOL*pole2)) THEN
            IF(strain1.GT.(TOL*pole1)) strain1=TOL*pole1 !yielded
            IF(strain1.GT.(TOL*pole2)) strain2=TOL*pole2 !yielded
          ENDIF
          DENOM1=pole1-strain1
          DtoAlpha1=1.0d0
          IF(alpha1.GT.ZERO) DtoAlpha1=DENOM1**alpha1
          DENOM2=pole2-strain2
          DtoAlpha2=1.0d0
          IF(alpha2.GT.ZERO) DtoAlpha2=DENOM2**alpha2
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole1).OR.
     '      strain.GT.(TOL*pole2)) THEN
C           Calc slope for first set of params
            D_DWslope1=(2.0d0*coeff1*DENOM1*(DENOM1
     '        +2.0d0*alpha1*strain1)
     '        +coeff1*alpha1*(alpha1+1.0d0)*strain1*strain1)
     '        /(DtoAlpha1*DENOM1*DENOM1)
C           Calc slope for second set of params
            D_DWslope2=(2.0d0*coeff2*DENOM2*(DENOM2
     '        +2.0d0*alpha2*strain2)
     '        +coeff2*alpha2*(alpha2+1.0d0)*strain2*strain2)
     '        /(DtoAlpha2*DENOM2*DENOM2)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=0.5d0*(
     '        2.0d0*coeff1/DtoAlpha1
     '        +4.0d0*coeff1*alpha1*strain1/(DtoAlpha1*DENOM1)
     '        +coeff1*alpha1*(alpha1+1.0d0)*strain1*strain1
     '        /(DtoAlpha1*DENOM1*DENOM1)
     '        +2.0d0*coeff2/DtoAlpha2
     '        +4.0d0*coeff2*alpha2*strain2/(DtoAlpha2*DENOM2)
     '        +coeff2*alpha2*(alpha2+1.0d0)*strain2*strain2
     '        /(DtoAlpha2*DENOM2*DENOM2))
            IF(strain.GT.(TOL*pole1).OR.strain.GT.(TOL*pole2)) THEN
              D_DW(iterm)=ZERO !term above is a const for these cases
              IF(strain.GT.(TOL*pole1))
     '          D_DW(iterm)=D_DW(iterm)+D_DWslope1/2.0d0
              IF(strain.GT.(TOL*pole2))
     '          D_DW(iterm)=D_DW(iterm)+D_DWslope2/2.0d0
            ENDIF
          ELSE !Compressive strain
C           Use deriv of D_DW wrt strain at strain=0 for compress. slope
            D_DW(iterm)=(D_DWslope1+D_DWslope2)/2.0d0
          ENDIF
        ENDDO !iterm
      ENDIF

      CALL EXITS('D_ENERGY')
      RETURN
 9999 CALL ERRORS('D_ENERGY',ERROR)
      CALL EXITS('D_ENERGY')
      RETURN 1
      END


      SUBROUTINE D_PFRE_NE(PARAMTYPE,IBT,IDO,INP,NAN,NBH,NBJ,
     '  NFF,NGAP,NHE,nhs,NMNO,NPF,NPNE,nr,NRE,NSP,NW,nx,NXI,
     '  CE,CP,D_RE,D_RI3,D_TW,D_ZE,D_ZW,ES,FEXT,XE,XW,ZE,ZW,ZW1,
     '  ERROR,*)

C#### Subroutine: D_PFRE_NE
C###  Description:
C###    D_PFRE_NE evaluates contribution to derivatives of element
C###    residuals wrt material parameters D_RE(na,4) (associated with
C###    the hydrostatic pressure variable) for incompressible materials,
C###    arising from the stress constraintdue to pressure bcs.
C###    OR
C###    Evaluates contribution to derivatives of element residuals wrt
C###    geometric variables ES (associated with the hydrostatic
C###    pressure variable) for incompressible materials, arising from
C###    the stress constraint due to pressure b.c's.

C**** The following code is similar to PFRE_NE and should be kept in
C**** sych with PFRE_NE.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ipma50.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NFF(6),NGAP(NIM,NBM),NHE,
     '  nhs,NMNO(1:2,0:NOPM),NPF(15,NFM),NPNE(NNM,NBFM),nr,NRE(NEM),
     '  NSP(-2:2),NW,nx,NXI(-NIM:NIM)
      REAL*8 CE(NMM),CP(NMM,NPM),D_RE(NSM,NHM,NOPM),D_RI3(NHM*NSM),
     '  D_TW(3,3,NHM*NSM),D_ZE(NSM,NHM),D_ZW(NHM,NUM,NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM),XE(NSM,NJM),XW(NJM,NUM),
     '  ZE(NSM,NHM),ZW(NHM,NUM),ZW1(NHM,NUM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER i,iface,IFE,j,k,mi,mj,
     '  NB,NB1,NBFF,NCW,nf,ngi1,ngi2,ngi3,ng_near,nh1,nhs1,nhx1,ni,
     '  NITB,nj,noopti,ns,ns1,ns2
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AXU(3,3),AZ,AZ1,AZL(3,3),AZL1(3,3),
     '  AZU(3,3),AZU1(3,3),CW(NCW),
     '  D_AZ,D_AZL(3,3),D_AZU(3,3),D_DZNXN(3,3),D_EG(3,3),
     '  DELTA_ZE,D_TC33,DXIX(3,3),
     '  DXIXN(3,3),DXIZN(3,3),DXIZN1(3,3),DXNXI(3,3),DXNZN(3,3),
     '  DZNXI(3,3),DZNXI1(3,3),DZNXN(3,3),DZNXN1(3,3),EG(3,3),
     '  GXL(3,3),GXU(3,3),GZ,GZ1,GZL(3,3),GZL1(3,3),GZU(3,3),GZU1(3,3),
     '  RI1,RI2,RI3,RSIGN,RWX,SUM1,SUM2,TW(3,3),TWA,XI(3)
      DATA DELTA_ZE/1.D-8/

      CALL ENTERS('D_PFRE_NE',*9999)
      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMMX)',ERROR,*9999)
      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)
      DO mi=1,NITB
        XI(MI)=0.5d0
      ENDDO !mi

      DO iface=1,2
        IFE=iface+4
        nf=NFF(IFE)
        NBFF=NPF(10,nf) !basis fn for first geometric variable
        XI(3)=DBLE(iface-1)
C new MPN 15-Apr-96: nodal interpolation of material params
        CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,ERROR,*9999)
C old
C        DO il=1,ILT(1,nr,nx)
C          IF(ILP(il,1,nr,nx).EQ.3) THEN
C            CW(il)=0.0d0
C            DO nnbf=1,NNT(NBFF)
C              CW(il)=CW(il)+CP(IL,NPNE(nnbf,nb))
C            ENDDO
C            CW(il)=CW(il)/DBLE(NNT(nbff))
C          ELSE IF(ILP(il,1,nr,nx).NE.3) THEN
C            CW(il)=CE(il)
C          ENDIF
C        ENDDO !il
C       Interpolate midwall geometric vars XW and derivs wrt Xi
        CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XW,XI,ERROR,*9999)
C       Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C       derivs (DXIXN) of Xi wrt Nu  coords.
        CALL XGMG(1,NITB,nb,nr,DXIXN,GXL,GXU,RWX,XW,ERROR,*9999)
C       Interpolate dependent var.s ZW and derivs wrt Nu
        CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '    DXIXN,ZE,ZW,XI,ERROR,*9999)
C       Calculate deformed metric tensors wrt Nu (AZL,AZU)
        CALL ZGMG(nb,AZ,AZL,AZU,ZW,ERROR,*9999)
        CALL ASSERT(KTYP51(nr).EQ.3.OR.KTYP51(nr).EQ.4,
     '    '>>Analytic derivatives only'
     '    //' implemented for 3D and membrane problems',ERROR,*9999)
        IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
C         Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C         stress tensor (D_TW - undeformed Nu coordinates) wrt
C         each of the material parameters
C         Note: AZ,AZL,AZU,EG are all independent of material params
C         so the calculation of D_AZ,D_AZL,D_AZU is not
C         needed here.
          IF(KTYP51(nr).EQ.3) THEN
            RI3=AZ !this saves calling ZGTG53 for RI3
            DO noopti=1,NTOPTI
              CALL D_ZGTG53(PARAMTYPE,nb,NMNO(1,noopti),nr,nx,
     '          AXU,AZ,AZL,AZU,CW,D_AZ,D_AZL,D_AZU,D_EG,
     '          D_RI3(noopti),D_TW(1,1,noopti),D_ZW(1,1,noopti),EG,ZW,
     '          ERROR,*9999)
            ENDDO !noopti
          ELSE IF(KTYP51(nr).EQ.4) THEN
            RI3=1.0d0 !this saves calling ZGTG54 for RI3
            DO noopti=1,NTOPTI
              CALL D_ZGTG54(PARAMTYPE,NMNO(1,noopti),nr,nx,
     '          AXU,AZ,AZL,AZU,CW,D_AZ,D_AZL,D_AZU,D_EG,
     '          D_RI3(noopti),D_TW(1,1,noopti),EG,ERROR,*9999)
            ENDDO !noopti
          ENDIF

        ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
C         Calc contravariant cpts of 2nd Piola-Kirchhoff stress
C         tensor (TW) wrt undeformed Nu coordinates
          IF(KTYP51(nr).EQ.3) THEN
            CALL ZGTG53(nb,nr,nx,AXU,AZ,AZL,AZU,
     '        CW,EG,RI1,RI2,RI3,TW,XW,ZW,ERROR,*9999)
          ELSE IF(KTYP51(nr).EQ.4) THEN
            CALL ZGTG54(nb,nr,nx,AXU,AZ,AZL,AZU,
     '        CW,EG,RI1,RI2,RI3,TW,ERROR,*9999)
          ENDIF
C         Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C         stress tensor (D_TW - undeformed Nu coordinates) and deriv
C         of third strain invariant (D_RI3) wrt each of the geometric
C         variables
          nhs1=0
          DO nhx1=1,NH_LOC(0,nx)
            nh1=NH_LOC(nhx1,nx)
            NB1=NBH(nh1)
            DO ns1=1,NST(NB1)+NAT(NB1)
              nhs1=nhs1+1
C             Calc D_ZW, deriv of def coords (ZW) wrt ele coords (ZE)
              D_ZE(ns1,nhx1)=1.0d0 !diff'ing wrt ns1,nh1 elem coord
              CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '          DXIXN,D_ZE,D_ZW(1,1,nhs1),XI,ERROR,*9999)
              D_ZE(ns1,nhx1)=0.0d0
C             Calculate derivs of deformed metric tensors, D_AZL and
C             ..D_AZU (Nu coords) and deriv of the det of AZL,
C             ..D_AZ wrt the current geom var using finite diff
C             ..approximations
C             Perturb deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
C             Interpolate perturbed dep vars ZW1 and derivs wrt Nu
              CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '          DXIXN,ZE,ZW1,XI,ERROR,*9999)
C             Reset deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE
C             Calc deformed metric tensors for perturbed coords
              CALL ZGMG(nb,AZ1,AZL1,AZU1,ZW1,ERROR,*9999)
C             Finite diff approx to derivs
              DO i=1,3
                DO j=1,3
                  D_AZL(i,j)=(AZL1(i,j)-AZL(i,j))/DELTA_ZE
                  D_AZU(i,j)=(AZU1(i,j)-AZU(i,j))/DELTA_ZE
                ENDDO !j
              ENDDO !i
              D_AZ=(AZ1-AZ)/DELTA_ZE
C             Calc derivs of TW wrt current geom var analytically
              IF(KTYP51(nr).EQ.3) THEN
                CALL D_ZGTG53(PARAMTYPE,nb,0,nr,nx,
     '            AXU,AZ,AZL,AZU,CW,D_AZ,D_AZL,D_AZU,D_EG,
     '            D_RI3(nhs1),D_TW(1,1,nhs1),D_ZW(1,1,nhs1),EG,ZW,
     '            ERROR,*9999)
              ELSE IF(KTYP51(nr).EQ.4) THEN
                CALL D_ZGTG54(PARAMTYPE,0,nr,nx,
     '            AXU,AZ,AZL,AZU,CW,D_AZ,D_AZL,D_AZU,D_EG,
     '            D_RI3(nhs1),D_TW(1,1,nhs1),EG,ERROR,*9999)
              ENDIF
            ENDDO !ns1
          ENDDO !nhx1
        ENDIF

C       Interpolate dependent var.s ZG and derivs wrt Xj
        CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '    DXIX,ZE,ZW,XI,ERROR,*9999)
        CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZW,ERROR,*9999)
C       Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
        CALL DXIDNU(NBJ(1),nr,DXIXN,DXNXI,GXL,GXU,XW,ERROR,*9999)
C       Get derivs of Xi wrt deformed Nu coords,DXIZN
        CALL DXIDNU(NBJ(1),nr,DXIZN,DZNXI,GZL,GZU,XW,ERROR,*9999)
C       Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
C       and inverse (DXNZN)
        DO ni=1,NITB
          DO mi=1,NITB
            SUM1=0.0d0
            SUM2=0.0d0
            DO k=1,NITB
              SUM1=SUM1+DZNXI(ni,k)*DXIXN(k,mi)
              SUM2=SUM2+DXNXI(ni,k)*DXIZN(k,mi)
            ENDDO
            DZNXN(ni,mi)=SUM1
            DXNZN(ni,mi)=SUM2
          ENDDO !mi
        ENDDO !ni
        IF(KTYP53(nr).EQ.3) THEN !Active stress component included
C!!!      Pick Gauss pt nearest to centre of current face in
C!!!      current element. To be strictly correct need to either
C!!!      store FEXT at face basis Gauss pts or somehow interpolate
C!!!      element Gauss pt values to the central point on the face
          ngi1=NGAP(1,nb)
          ngi2=NGAP(2,nb)
          ngi3=NGAP(3,nb)
          ng_near=(ngi1+1)/2 + ((ngi2-1)/2)*ngi1
     '      + (iface-1)*(ngi3-1)*ngi2*ngi1 !need integer division
          CALL ZGTG5A(nr,FEXT(1,ng_near),DXNZN,DZNXN,
     '      CW(IL_time_delay),TW,TWA,ERROR,*9999)
C!!!      NOTE: The code below assumes that the active fibre stress
C!!!      component is independent of the
C!!!      passive material/geometrical parameters
        ENDIF

        IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
          DO noopti=1,NTOPTI
C           Get deriv of Physical Cauchy stress D_TC33
C           wrt deformed Nu_3 from D_TW
            D_TC33=0.0d0
            DO nj=1,NITB
              DO mj=1,NITB
                D_TC33=D_TC33+DZNXN(3,nj)*DZNXN(3,mj)
     '            *D_TW(nj,mj,noopti)
              ENDDO !mj
            ENDDO !nj
            D_TC33=D_TC33/DSQRT(RI3)
            IF(iface.EQ.1) THEN !Xi3=0 face
              RSIGN= 1.0d0
            ELSE IF(iface.EQ.2) THEN !Xi3=1 face
              RSIGN=-1.0d0
            ENDIF
            ns=NSP(iface)
            IF(iface.EQ.1.AND.(NXI(-3).EQ.0.OR.nr.NE.NRE(NXI(-3))).AND.
     '        (NW.EQ.2.OR.NW.EQ.4).OR. !ext. press bc appl on Xi3=0 face
     '        iface.EQ.2.AND.(NXI(3).EQ.0.OR.nr.NE.NRE(NXI(3))).AND.
     '        (NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc appl on Xi3=1 face
C             match norm stress for press bcs
              D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)=D_TC33
            ELSE !no external pressure bc applied on current face
              IF(KTYP5A(nr).EQ.1) THEN !match hyd. press across elems
                D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)=0.0d0
              ELSE IF(KTYP5A(nr).EQ.2) THEN !match norm stress acr elems
                D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)=D_TC33*RSIGN
              ENDIF
            ENDIF

            ns2=NSP(-iface) !ns for pressure bc params
            IF(KTYP5A(nr).EQ.1) THEN !match hyd. press across elems
              D_RE(ns2,NH_LOC(NH_LOC(0,nx),nx),noopti)=0.0d0
            ELSE IF(KTYP5A(nr).EQ.2) THEN !match norm stress acr elems
              D_RE(ns2,NH_LOC(NH_LOC(0,nx),nx),noopti)=D_TC33*RSIGN
            ENDIF
          ENDDO !noopti

        ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
          IF(NW.EQ.3) THEN
            nhs=nhs+2
          ELSE
            nhs=nhs+1
          ENDIF
          nhs1=0
          DO nhx1=1,NH_LOC(0,nx)
            nh1=NH_LOC(nhx1,nx)
            NB1=NBH(nh1)
            DO ns1=1,NST(NB1)+NAT(NB1)
              nhs1=nhs1+1
C             Perturb deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
C             Interpolate perturbed dep var.s ZW1 and derivs wrt Xi
              CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '          DXIXN,ZE,ZW1,XI,ERROR,*9999)
C             Reset deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE
C             Calculate deformed metric tensors wrt Xi (GZL1,GZU1)
C             ..for perturbed elem coords
              CALL ZGMG(NBJ(1),GZ1,GZL1,GZU1,ZW1,ERROR,*9999)
C             Get derivs of Xi wrt deformed Nu coords, DXIZN1
C             ..for perturbed elem coords
              CALL DXIDNU(NBJ(1),nr,DXIZN1,DZNXI1,GZL1,GZU1,XW,
     '          ERROR,*9999)
C             Calculate derivs of deformed Nu wrt undeformed Nu
C             ..for perturbed elem coords
              ni=3
              DO mi=1,NITB
                SUM1=0.0d0
                SUM2=0.0d0
                DO k=1,NITB
                  SUM1=SUM1+DZNXI1(ni,k)*DXIXN(k,mi)
                ENDDO !k
                DZNXN1(ni,mi)=SUM1
              ENDDO !mi
C             Finite diff approx to derivatives
              DO i=1,3
                DO j=1,3
                  D_DZNXN(i,j)=(DZNXN1(i,j)-DZNXN(i,j))/DELTA_ZE
                ENDDO !j
              ENDDO !i
C             Calc deriv of Cauchy Stress wrt current geometric var.
              D_TC33=0.0d0
              DO nj=1,NITB
                DO mj=1,NITB
                  D_TC33=D_TC33
     '              +D_DZNXN(3,nj)*DZNXN(3,mj)*TW(nj,mj)/DSQRT(RI3)
     '              +DZNXN(3,nj)*D_DZNXN(3,mj)*TW(nj,mj)/DSQRT(RI3)
     '              +DZNXN(3,nj)*DZNXN(3,mj)*D_TW(nj,mj,nhs1)
     '              /DSQRT(RI3)
     '              -DZNXN(3,nj)*DZNXN(3,mj)*TW(nj,mj)*D_RI3(nhs1)
     '              /(2.0d0*DSQRT(RI3*RI3*RI3))
                ENDDO !mj
              ENDDO !nj

              IF(iface.EQ.1) THEN !Xi3=0 face
                RSIGN= 1.0d0
              ELSE IF(iface.EQ.2) THEN !Xi3=1 face
                RSIGN=-1.0d0
              ENDIF

              IF(iface.EQ.1.AND.(NXI(-3).EQ.0.OR.nr.NE.NRE(NXI(-3)))
     '          .AND.
     '          (NW.EQ.2.OR.NW.EQ.4).OR. !ext pres bc appl on Xi3=0 face
     '          iface.EQ.2.AND.(NXI(3).EQ.0.OR.nr.NE.NRE(NXI(3))).AND.
     '          (NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc on Xi3=1 face
                ES(nhs,nhs1)=D_TC33 !match norm stress for press bcs
              ELSE !no external pressure bc applied on current face
                IF(KTYP5A(nr).EQ.1) THEN !match hyd. press across elems
                  ES(nhs,nhs1)=
     '              D_ZW(NH_LOC(NH_LOC(0,nx),nx),1,nhs1)*RSIGN
C!!! I don't think D_ZW has been set up
                ELSE IF(KTYP5A(nr).EQ.2) THEN !match norm strs acr elems
                  ES(nhs,nhs1)=D_TC33*RSIGN
                ENDIF
              ENDIF

              CALL ASSERT(.FALSE.,'>>Needs fixing',ERROR,*9999)
c              ns2=NSP(-iface) !ns for pressure bc params
cC!!! to be fixed: determine nhs2 from NSP(-iface)
c              IF(KTYP5A(nr).EQ.1) THEN !match hyd. press across elems
c                ES(nhs2,nhs1)=D_ZW(NH_LOC(NH_LOC(0,nx),nx),1,nhs1)*RSIGN
cC!!! I don't think D_ZW has been set up
c              ELSE IF(KTYP5A(nr).EQ.2) THEN !match norm stress acr elems
c                ES(nhs2,nhs1)=D_TC33*RSIGN
c              ENDIF

            ENDDO !ns1
          ENDDO !nhx1
        ENDIF
      ENDDO !iface

      CALL EXITS('D_PFRE_NE')
      RETURN
 9999 CALL ERRORS('D_PFRE_NE',ERROR)
      CALL EXITS('D_PFRE_NE')
      RETURN 1
      END


      SUBROUTINE D_PFRF(IBT,IDO,INP,IXF,
     '  NAN,NBH,NBJ,NGAP,NHE,nhx1,NPF,nr,ns1,nx,
     '  D_RF,D_ZE,D_ZG,PF,PG,WG,XE,XG,ZE,ZG,ZG1,ERROR,*)

C#### Subroutine: D_PFRF
C###  Description:
C###    D_PFRF evaluates derivatives wrt geometric variables of
C###    contributions to element residuals RE from the pressure PF(i)
C###    acting on the Xi(3)=IXF face (IXF=iface-1).

C**** Note: GZL & GZU are deformed state metric tensors wrt Xi.
C**** The following code was copied from PFRF and altered and should
C**** be kept in sych with PFRF.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  IXF,NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NGAP(NIM,NBM),
     '  NHE,nhx1,NPF(15),nr,ns1,nx
      REAL*8 D_RF(32,3),D_ZE(NSM,NHM),D_ZG(NHM,NUM),PF,
     '  PG(NSM,NUM,NGM,NBM),WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,NBFF,ng,ng1,ng2,NGI1,NGI2,nh,nhx,NITB,ni,ns,NU1(0:3)
      REAL*8 D(5,5),DELTA_ZE,D_GZ,D_GZU(3,3),D_RGZ,D_RWG,D_SUM,
     '  DXIX(3,3),GZ,GZ1,GZL(3,3),GZL1(3,3),
     '  GZU(3,3),GZU1(3,3),RGZ,RWG,SUM,XI(3)
      DATA NU1/1,2,4,7/
      DATA D/5*0.0d0,-0.288675134594813d0,0.288675134594813d0,3*0.0d0,
     '       -0.387298334620741d0,0.0d0,0.387298334620741d0,2*0.0d0,
     '       -0.430568155797026d0,    -0.169990521792428d0,
     '        0.169990521792428d0,     0.430568155797026d0,0.0d0,
     '       -0.453089922969332d0,    -0.269234655052841d0,0.0d0,
     '        0.269234655052841d0,     0.453089922969332d0/
      DATA DELTA_ZE/1.0d-8/

      CALL ENTERS('D_PFRF',*9999)
      NITB=NIT(NBH(NH_LOC(1,nx)))

      DO nhx=1,NJ_LOC(NJL_GEOM,0)
        nh=NH_LOC(nhx,nx)
        NBFF=NPF(9+nhx)
        DO ns=1,NST(NBFF)+NAT(NBFF)
          D_RF(ns,nh)=0.0d0
        ENDDO
      ENDDO

      XI(3)=DBLE(IXF)
      NGI1=NGAP(1,NPF(10))
      NGI2=NGAP(2,NPF(10))
      ng=0
      DO ng2=1,NGI2
        DO ng1=1,NGI1
          ng=ng+1
          XI(1)=0.5d0+D(ng1,NGI1)
          XI(2)=0.5d0+D(ng2,NGI2)
          CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
          CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG,XI,ERROR,
     '      *9999)
          CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)
          RGZ=DSQRT(GZ)
          RWG=RGZ*WG(ng,NPF(10))
C         Calc D_ZG, deriv of def coords (ZG), wrt elem coords (ZE)
          D_ZE(ns1,nhx1)=1.0d0 !differentiating wrt ns1,nhx1 elem coord
          CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,D_ZE,D_ZG,XI,
     '      ERROR,*9999)
          D_ZE(ns1,nhx1)=0.0d0
C         Perturb deformed element coords
          ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
          CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG1,XI,
     '      ERROR,*9999)
C         Reset deformed element coords
          ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE
          CALL ZGMG(NBH(NH_LOC(1,nx)),GZ1,GZL1,GZU1,ZG1,ERROR,*9999)
C         Finite diff approx to derivatives
          DO I=1,3
            DO J=1,3
              D_GZU(i,j)=(GZU1(i,j)-GZU(i,j))/DELTA_ZE
            ENDDO
          ENDDO
          D_GZ=(GZ1-GZ)/DELTA_ZE
          D_RGZ=D_GZ/(2.d0*DSQRT(GZ))
          D_RWG=D_RGZ*WG(ng,NPF(10))

C old MPN 28-Jun-1995: integrals now wrt deformed coords
C        IF(ITYP10(nr).EQ.1) THEN
C          DO nhx=1,NJ_LOC(NJL_GEOM,0)
C            nh=NH_LOC(nhx,nx)
C            DO ni=1,NITB
C              DZI(nhx,ni)=ZG(nhx,NU1(ni))
C              D_DZI(nhx,ni)=D_ZG(nhx,NU1(ni))
C            ENDDO
C          ENDDO
C        ELSE IF(ITYP10(nr).EQ.2) THEN
C          RX=XG(1,1)
C          RZ=ZG(1,1)
C          D_RZ=D_ZG(1,1)
C          DT=ZG(2,1)-XG(2,1)
C          D_DT=D_ZG(2,1)
C          CT=DCOS(DT)
C          D_CT=-DSIN(DT)*D_DT
C          ST=DSIN(DT)
C          D_ST= DCOS(DT)*D_DT
C          DO ni=1,NITB
C            DZI(1,ni)= ZG(1,NU1(ni))*CT-RZ*ST*ZG(2,NU1(ni))
C            D_DZI(1,ni)= D_ZG(1,NU1(ni))*CT +ZG(1,NU1(ni))*D_CT
C     '        -D_RZ*ST*ZG(2,NU1(ni)) -RZ*D_ST*ZG(2,NU1(ni))
C     '                               -RZ*ST*D_ZG(2,NU1(ni))
C            DZI(2,ni)=(ZG(1,NU1(ni))*ST+RZ*CT*ZG(2,NU1(ni)))/RX
C            D_DZI(2,ni)=(D_ZG(1,NU1(ni))*ST +ZG(1,NU1(ni))*D_ST
C     '        +D_RZ*CT*ZG(2,NU1(ni)) +RZ*D_CT*ZG(2,NU1(ni))
C     '                               +RZ*CT*D_ZG(2,NU1(ni)))/RX
C            DZI(3,ni)= ZG(3,NU1(ni))
C            D_DZI(3,ni)= D_ZG(3,NU1(ni))
C          ENDDO
C        ELSE IF(ITYP10(nr).EQ.3) THEN
C          RX=XG(1,1)
C          RZ=ZG(1,1)
C          D_RZ=D_ZG(1,1)
C          DT=ZG(2,1)-XG(2,1)
C          D_DT=D_ZG(2,1)
C          CT=DCOS(DT)
C          D_CT=-DSIN(DT)*D_DT
C          ST=DSIN(DT)
C          D_ST= DCOS(DT)*D_DT
C          CX=DCOS(XG(3,1))
C          SX=DSIN(XG(3,1))
C          CZ=DCOS(ZG(3,1))
C          D_CZ=-DSIN(ZG(3,1))*D_ZG(3,1)
C          SZ=DSIN(ZG(3,1))
C          D_SZ= DCOS(ZG(3,1))*D_ZG(3,1)
C          CC=CX*CZ
C          D_CC=CX*D_CZ
C          SS=SX*SZ
C          D_SS=SX*D_SZ
C          CS=CX*SZ
C          D_CS=CX*D_SZ
C          SC=SX*CZ
C          D_SC=SX*D_CZ
C          DO ni=1,NITB
C            RB=ZG(1,NU1(ni))
C            D_RB=D_ZG(1,NU1(ni))
C            TB=ZG(2,NU1(ni))
C            D_TB=D_ZG(2,NU1(ni))
C            PB=ZG(3,NU1(ni))
C            D_PB=D_ZG(3,NU1(ni))
C            DZI(1,ni)=CC*(RB*CT-RZ*ST*TB)-CS*RZ*CT*PB+SC*RZ*PB+SS*RB
C            D_DZI(1,ni)=D_CC*RB*CT +CC*D_RB*CT +CC*RB*D_CT
C     '        -D_CC*RZ*ST*TB-CC*D_RZ*ST*TB-CC*RZ*D_ST*TB-CC*RZ*ST*D_TB
C     '        -D_CS*RZ*CT*PB-CS*D_RZ*CT*PB-CS*RZ*D_CT*PB-CS*RZ*CT*D_PB
C     '        +D_SC*RZ*PB+SC*D_RZ*PB+SC*RZ*D_PB  +  D_SS*RB+SS*D_RB
C            DZI(2,ni)=(RZ*CZ*CT*TB+(RB*CZ-RZ*SZ*PB)*ST)/(RX*CX)
C            D_DZI(2,ni)=(D_RZ*CZ*CT*TB +RZ*D_CZ*CT*TB +RZ*CZ*D_CT*TB
C     '                                                +RZ*CZ*CT*D_TB
C     '        +D_RB*CZ*ST +RB*D_CZ*ST +RB*CZ*D_ST
C     '        -D_RZ*SZ*PB*ST -RZ*D_SZ*PB*ST -RZ*SZ*D_PB*ST
C     '                                      -RZ*SZ*PB*D_ST)/(RX*CX)
C            DZI(3,ni)=(CC*RZ*PB+CS*RB+SC*(RZ*ST*TB-RB*CT)+SS*RZ*CT*PB)
C     '                                                 /RX
C            D_DZI(3,ni)=(D_CC*RZ*PB +CC*D_RZ*PB +CC*RZ*D_PB
C     '        +D_CS*RB +CS*D_RB
C     '        +D_SC*RZ*ST*TB+SC*D_RZ*ST*TB+SC*RZ*D_ST*TB+SC*RZ*ST*D_TB
C     '        -D_SC*RB*CT -SC*D_RB*CT -SC*RB*D_CT
C     '        +D_SS*RZ*CT*PB+SS*D_RZ*CT*PB+SS*RZ*D_CT*PB+SS*RZ*CT*D_PB)
C     '                                                 /RX
C          ENDDO
C        ELSE IF(ITYP10(nr).EQ.4) THEN
C          SLX=DSINH(XG(1,1))
C          SLZ=DSINH(ZG(1,1))
C          D_SLZ=DCOSH(ZG(1,1))*D_ZG(1,1)
C          SMX=DSIN (XG(2,1))
C          SMZ=DSIN (ZG(2,1))
C          D_SMZ=DCOS (ZG(2,1))*D_ZG(2,1)
C          CLX=DSQRT(1.0d0+SLX*SLX)
C          CLZ=DSQRT(1.0d0+SLZ*SLZ)
C          D_CLZ= SLZ*D_SLZ/DSQRT(1.0d0+SLZ*SLZ)
C          CMX=DSQRT(1.0d0-SMX*SMX)
C          CMZ=DSQRT(1.0d0-SMZ*SMZ)
C          D_CMZ=-SMZ*D_SMZ/DSQRT(1.0d0-SMZ*SMZ)
C          CSLX=CLX/SLX
C          CSMX=CMX/SMX
C          DT=ZG(3,1)-XG(3,1)
C          D_DT=D_ZG(3,1)
C          CT=DCOS(DT)
C          D_CT=-DSIN(DT)*D_DT
C          ST=DSIN(DT)
C          D_ST= DCOS(DT)*D_DT
C          CCL=CLX*CLZ
C          D_CCL=CLX*D_CLZ
C          CSL=CLX*SLZ
C          D_CSL=CLX*D_SLZ
C          SCL=SLX*CLZ
C          D_SCL=SLX*D_CLZ
C          SSL=SLX*SLZ
C          D_SSL=SLX*D_SLZ
C          CC=CMX*CMZ
C          D_CC=CMX*D_CMZ
C          CS=CMX*SMZ
C          D_CS=CMX*D_SMZ
C          SC=SMX*CMZ
C          D_SC=SMX*D_CMZ
C          SS=SMX*SMZ
C          D_SS=SMX*D_SMZ
C          G1=SLX*SLX+SMX*SMX
C          G3=SLX*SLX*SMX*SMX
C          DO ni=1,NITB
C            DLB=ZG(1,NU1(ni))
C            D_DLB=D_ZG(1,NU1(ni))
C            DMB=ZG(2,NU1(ni))
C            D_DMB=D_ZG(2,NU1(ni))
C            DTB=ZG(3,NU1(ni))
C            D_DTB=D_ZG(3,NU1(ni))
C            DZI(1,ni)=(( SSL*CC+CCL*SS*CT)*DLB+(-SCL*CS+CSL*SC*CT)*DMB
C     '                                              -CSL*SS*ST*DTB)/G1
C            D_DZI(1,ni)=(D_SSL*CC*DLB +SSL*D_CC*DLB +SSL*CC*D_DLB
C     '        +D_CCL*SS*CT*DLB +CCL*D_SS*CT*DLB +CCL*SS*D_CT*DLB
C     '                                          +CCL*SS*CT*D_DLB
C     '        -D_SCL*CS*DMB -SCL*D_CS*DMB -SCL*CS*D_DMB
C     '        +D_CSL*SC*CT*DMB +CSL*D_SC*CT*DMB +CSL*SC*D_CT*DMB
C     '                                          +CSL*SC*CT*D_DMB
C     '        -D_CSL*SS*ST*DTB -CSL*D_SS*ST*DTB -CSL*SS*D_ST*DTB
C     '                                          -CSL*SS*ST*D_DTB)/G1
C            DZI(2,ni)=((-CSL*SC+SCL*CS*CT)*DLB+( CCL*SS+SSL*CC*CT)*DMB
C     '                                              -SSL*CS*ST*DTB)/G1
C            D_DZI(2,ni)=(-D_CSL*SC*DLB -CSL*D_SC*DLB -CSL*SC*D_DLB
C     '        +D_SCL*CS*CT*DLB +SCL*D_CS*CT*DLB +SCL*CS*D_CT*DLB
C     '                                          +SCL*CS*CT*D_DLB
C     '        +D_CCL*SS*DMB +CCL*D_SS*DMB +CCL*SS*D_DMB
C     '        +D_SSL*CC*CT*DMB +SSL*D_CC*CT*DMB +SSL*CC*D_CT*DMB
C     '                                          +SSL*CC*CT*D_DMB
C     '        -D_SSL*CS*ST*DTB -SSL*D_CS*ST*DTB -SSL*CS*D_ST*DTB
C     '                                          -SSL*CS*ST*D_DTB)/G1
C            DZI(3,ni)=(SCL*SS*ST*DLB+SSL*SC*ST*DMB+SSL*SS*CT*DTB)/G3
C            D_DZI(3,ni)=(D_SCL*SS*ST*DLB +SCL*D_SS*ST*DLB
C     '                  +SCL*SS*D_ST*DLB +SCL*SS*ST*D_DLB
C     '        +D_SSL*SC*ST*DMB +SSL*D_SC*ST*DMB +SSL*SC*D_ST*DMB
C     '                                          +SSL*SC*ST*D_DMB
C     '        +D_SSL*SS*CT*DTB +SSL*D_SS*CT*DTB +SSL*SS*D_CT*DTB
C     '                                          +SSL*SS*CT*D_DTB)/G3
C          ENDDO
C        ENDIF
          DO nhx=1,NJ_LOC(NJL_GEOM,0)
            nh=NH_LOC(nhx,nx)
            NBFF=NPF(9+nhx)
            SUM=0.0d0
            D_SUM=0.0d0
            DO ni=1,NITB
C MPN 28-Jun-1995: integrals now wrt deformed coords
              SUM=SUM+GZU(ni,3)*ZG(nhx,NU1(ni))
              D_SUM=D_SUM+D_GZU(ni,3)*ZG(nhx,NU1(ni))
     '          +GZU(ni,3)*D_ZG(nhx,NU1(ni))
C old       SUM=SUM+GZU(ni,3)*DZI(nhx,ni)
C old       D_SUM=D_SUM+D_GZU(ni,3)*DZI(nhx,ni)+GZU(ni,3)*D_DZI(nhx,ni)
            ENDDO
            DO ns=1,NST(NBFF)+NAT(NBFF)
              IF(IXF.EQ.1)THEN !positive face
                D_RF(ns,nh)=D_RF(ns,nh)-PF*PG(ns,1,ng,NBFF)
     '            *(D_SUM*RWG+SUM*D_RWG)
              ELSE !negative face
                D_RF(ns,nh)=D_RF(ns,nh)+PF*PG(ns,1,ng,NBFF)
     '            *(D_SUM*RWG+SUM*D_RWG)
              ENDIF
            ENDDO
          ENDDO
        ENDDO !ngi1
      ENDDO !ngi2

      IF(DOP) THEN
        DO nhx=1,NJ_LOC(NJL_GEOM,0)
          nh=NH_LOC(nhx,nx)
          NBFF=NPF(9+nhx)
          WRITE(OP_STRING,
     '      '('' D_RF(ns,'',I2,''): '',5D12.4)')
     '      nh,(D_RF(ns,nh),NS=1,NST(NBFF)+NAT(NBFF))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('D_PFRF')
      RETURN
 9999 CALL ERRORS('D_PFRF',ERROR)
      CALL EXITS('D_PFRF')
      RETURN 1
      END


      SUBROUTINE D_ZERE50(PARAMTYPE,IBT,IDO,INP,NAN,NBH,NBJ,NGAP,ne,
     '  NFF,NHE,NJE,NKF,NMNO,NNF,NPF,NPNE,nr,NRE,NW,nx,NXI,
     '  CE,CG,CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT,PG,RGX,SE,VE,WG,XE,XG,
     '  ZE,D_ZE,ZG,ZG1,ERROR,*)

C#### Subroutine: D_ZERE50
C###  Description:
C###    D_ZERE50 calculates derivatives of element residual D_RE from
C###    current dependent variable array ZE.

C**** The following code was copied from ZERE50 and altered and should
C**** be kept in synch with ZERE50.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ipma50.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NGAP(NIM,NBM),ne,NFF(6),
     '  NHE,NJE,NKF(0:4,16,6,NBFM),NMNO(1:2,0:NOPM),NNF(0:17,6,NBFM),
     '  NPF(15,NFM),NPNE(NNM,NBFM,NEFM),nr,NRE(NEM),NW,
     '  nx,NXI(-NIM:NIM,0:NEM)
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),D_RE(NSM,NHM,NOPM),
     '  D_RI3(NHM*NSM),
     '  D_TG(3,3,NHM*NSM),D_ZG(NHM,NUM,NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),RGX(NGM),SE(NSM,NBFM,NEFM),VE(NSM,NKM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),ZE(NSM,NHM),D_ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZG1(NHM,NUM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER i,iface,IFE,IXF,j,k,NA,NATB,nb,NB1,NBE,
     '  NBFF,NBP,nf,ng,nh,nh1,nhs,nhs1,NHSGMAX,NHSTART,nhx,nhx1,
     '  NI,NITB,nj,njj,njj1,njj2,nk,nkbf,nn,nnbf,noopti,
     '  ns,ns1,nsa,nse,nsf,NSP(-2:2),NSTB,NU1(0:3)
      REAL*8 AXU(3,3),
     '  AZ,AZ1,AZL(3,3),AZL1(3,3),AZU(3,3),AZU1(3,3),
     '  CHTOFF(3,3,3),CHTOFF1(3,3,3),CLZ,CMZ,CSLZ,CSMZ,
     '  D11,D12,D13,D21,D22,D23,D31,D32,D33,
     '  D_AG,D_AGE,D_AGE1,D_AGE2,D_AGE3,Darcy_Resid,D_AZ,
     '  D_AZL(3,3),D_AZU(3,3),DBM(3,3,3),
     '  D_CHTOFF(3,3,3),D_Darcy_Resid,D_delsqP,
     '  D_EG(3,3),delsqP,DELTA_ZE,
     '  DLA1,DLA2,DLA3,DMA1,DMA2,DMA3,
     '  D_RF(32,3),D_RZWG,D_SUM,D_SUM_MEMB(3,40),
     '  DTA1,DTA2,DTA3,DXIX(3,3),
     '  E1,E2,EG(3,3),G1,GXL(3,3),GXU(3,3),
     '  GZ,GZ1,GZL(3,3),GZL1(3,3),GZU(3,3),GZU1(3,3),
     '  PF(2),PGA1,PGA2,PGA3,PGG,PGX,PPG(64,4),PPGG(4),
     '  RI1,RI2,RI3,RWG,RZWG,RZWG1,
     '  SLZ,SMZ,SUM,TG(3,3),X3G(4,3)
      CHARACTER CHAR1*3,CHAR2*3,TYPE*9
      LOGICAL BCPARAM,ELEMPRESS,NODEPRESS,SAMEDEPBASIS
      DATA DELTA_ZE/1.0D-8/
      DATA NU1/1,2,4,7/

      CALL ENTERS('D_ZERE50',*9999)
      NITB=NIT(NBJ(1))

C *** Test whether to use unrolled loops
      SAMEDEPBASIS=.TRUE.
      nh1=NH_LOC(1,nx)
      DO nhx=2,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        IF(NBH(nh).NE.NBH(nh1)) SAMEDEPBASIS=.FALSE.
      ENDDO !nhx
      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        IF(ITYP10(nr).EQ.1.AND.KTYP51(nr).EQ.3.AND.KTYP53(nr).LE.3.AND.
     '    NJ_LOC(NJL_GEOM,0).EQ.3.AND.NIT(NBH(NH_LOC(1,nx))).EQ.3) THEN
          TYPE='RC3D'
        ELSE IF(ITYP10(nr).EQ.4.AND.KTYP51(nr).EQ.3.AND.
     '      KTYP53(nr).LE.3.AND.SAMEDEPBASIS.AND..NOT.DOP.AND.
     '      NJ_LOC(NJL_GEOM,0).EQ.3.AND.
     '      NIT(NBH(NH_LOC(1,nx))).EQ.3) THEN
          TYPE='PROLATE'
        ELSE
          TYPE=' '
        ENDIF
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
            DO noopti=1,NTOPTI
              D_RE(ns,nh,noopti)=0.0d0
            ENDDO !noopti
          ENDDO !ns
        ENDDO !nhx
        CALL ASSERT(NTOPTI.LE.NHM*NSM,'>>Increase dimension of D_RE',
     '    ERROR,*9999)

      ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
        TYPE=' '
        nhs=0
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
            D_ZE(ns,nhx)=0.0d0
            nhs=nhs+1
            nhs1=0
            DO nhx1=1,NH_LOC(0,nx)
              nh1=NH_LOC(nhx1,nx)
              DO ns1=1,NST(NBH(nh1))+NAT(NBH(nh1))
                nhs1=nhs1+1
                ES(nhs,nhs1)=0.0d0
              ENDDO !ns1
            ENDDO !nhx1
          ENDDO !ns
        ENDDO !nhx
        IF(KTYP51(nr).EQ.4.AND.NW.GT.1) THEN !membrane + press bc(s)
          WRITE(CHAR1,'(I3)') NJ_LOC(NJL_GEOM,0)
          WRITE(CHAR2,'(I3)') nhs1
          CALL ASSERT(NJ_LOC(NJL_GEOM,0).LE.3.AND.nhs1.LE.40,
     '      '>>Increase dimension of D_SUM_MEMB to be ('
     '      //CHAR1(1:3)//','//CHAR2(1:3)//')',ERROR,*9999)
        ENDIF
      ENDIF

C *** Set up PF array for pressure bcs from aux vars in ZE
      IF(KTYP57(nr).GT.1) THEN
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis function for press vars
        NODEPRESS=.FALSE.
        IF(NST(NBP).GT.0) NODEPRESS=.TRUE. !node based hyd press interp
        ELEMPRESS=.FALSE.
C       Put current XI3 face pressures into PF array
        NSP(1)=0
        NSP(2)=0
        DO na=1,NAT(NBP)
C         Check that pressure varies in Xi(3) dirn only for
C         element based pressure interpolation
          IF(NST(NBP).EQ.0.AND.  !elem based press interp
     '      (NAN(1,na,NBP).NE.0.OR.NAN(2,na,NBP).NE.0)) THEN
            ERROR='>>Pressure must vary only with Xi(3) '
     '        //'for elem based pressure interpolation'
            GOTO 9999
          ENDIF
          IF(NAN(3,na,NBP).EQ.0) THEN
C           Pick up param assoc with const pressure term
            NSP(1)=na
          ELSE IF(NAN(3,na,NBP).EQ.1.OR.NAN(3,na,NBP).EQ.3) THEN
C           Pick up param assoc with linear or cubic press term
            NSP(2)=na
          ELSE IF(NAN(3,na,NBP).EQ.-1) THEN
C           Pick up param assoc with Xi3=0 face pressure bc
            PF(1)=ZE(NST(NBP)+na,NH_LOC(0,nx))
            NSP(-1)=na
          ELSE IF(NAN(3,na,NBP).EQ.-2) THEN
C           Pick up param assoc with Xi3=1 face pressure bc
            PF(2)=ZE(NST(NBP)+na,NH_LOC(0,nx))
            NSP(-2)=na
          ENDIF
C         check for element based hyd press interp apart from
C         pressure boundary condition parameters
          BCPARAM=.FALSE.
          DO ni=1,NIT(NBP)
            IF(NAN(ni,na,NBP).LT.0) BCPARAM=.TRUE.
          ENDDO !ni
          IF(.NOT.BCPARAM) ELEMPRESS=.TRUE.
        ENDDO !na
      ENDIF

      CALL ASSERT(KTYP51(nr).EQ.3.OR.KTYP51(nr).EQ.4,
     '  '>>Analytic derivatives only'
     '  //' implemented for 3D and membrane problems',ERROR,*9999)

C *** Main Gauss point loop
      DO 50 ng=1,NGT(NBH(NH_LOC(1,nx)))
C       Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
        CALL XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*9999)
C       Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C       ..derivs (DXIX) of Xi wrt Nu coords.
        CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXIX,GXL,GXU,
     '    RGX(ng),XG,ERROR,*9999)
C       Interpolate dependent var.s ZG and derivs wrt Nu
        CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C       Calculate deformed metric tensors wrt Xj or Nu (AZL,AZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),AZ,AZL,AZU,ZG,ERROR,*9999)

       IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
C         Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C         stress tensor (D_TW - undeformed Nu coordinates) wrt
C         each of the material parameters
C         Note: AZ,AZL,AZU,EG are all independent of material params
C         so the calculation of D_AZ,D_AZL,D_AZU is not
C         needed here.
          IF(KTYP51(nr).EQ.3) THEN
            DO noopti=1,NTOPTI
              CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)),
     '          NMNO(1,noopti),nr,nx,
     '          AXU,AZ,AZL,AZU,CG(1,ng),D_AZ,D_AZL,D_AZU,D_EG,
     '          D_RI3(noopti),D_TG(1,1,noopti),D_ZG(1,1,noopti),EG,ZG,
     '          ERROR,*9999)
            ENDDO
          ELSE IF(KTYP51(nr).EQ.4) THEN
            DO noopti=1,NTOPTI
              CALL D_ZGTG54(PARAMTYPE,NMNO(1,noopti),nr,nx,
     '          AXU,AZ,AZL,AZU,CG(1,ng),D_AZ,D_AZL,D_AZU,D_EG,
     '          D_RI3(noopti),D_TG(1,1,noopti),EG,ERROR,*9999)
            ENDDO
          ENDIF

        ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
C         Calc contravariant cpts of 2nd Piola-Kirchhoff stress
C         tensor (TG) wrt undeformed Nu coordinates
          IF(KTYP51(nr).EQ.3) THEN
            CALL ZGTG53(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '        CG(1,ng),EG,RI1,RI2,RI3,TG,XG,ZG,ERROR,*9999)
          ELSE IF(KTYP51(nr).EQ.4) THEN
            CALL ZGTG54(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '        CG(1,ng),EG,RI1,RI2,RI3,TG,ERROR,*9999)
          ENDIF

C         Calc deriv's of contravariant cpts of 2nd Piola-Kirchhoff
C         stress tensor (D_TG - undeformed Nu coordinates) and deriv
C         of third strain invariant (D_RI3) wrt each of the geometric
C         variables
          nhs1=0
          DO nhx1=1,NH_LOC(0,nx)
            nh1=NH_LOC(nhx1,nx)
            NB1=NBH(nh1)
            DO ns1=1,NST(NB1)+NAT(NB1)
              nhs1=nhs1+1
C             Calc D_ZG, deriv of def coords (ZG + deriv wrt Nu coords)
C             wrt elem coords (ZE)
C             NOTE: code up to the next ! is similar to ZEZG with JP=1
              DO ni=0,NIT(NB1)
                IF(NI.EQ.0) THEN
                  D_ZG(nh1,NU1(ni),nhs1)=PG(ns1,NU1(ni),ng,NB1)
                ELSE  !derivs of basis fns wrt Nu coords
                  D_ZG(nh1,NU1(ni),nhs1)=PG(ns1,2,ng,NB1)*DXIX(1,ni) +
     '                                   PG(ns1,4,ng,NB1)*DXIX(2,ni) +
     '                                   PG(ns1,7,ng,NB1)*DXIX(3,ni)
                ENDIF
              ENDDO
C             end of ZEZG stuff
C             Calculate derivs of deformed metric tensors, D_AZL and
C             ..D_AZU (Nu coords) and deriv of the determinant of AZL,
C             ..D_AZ wrt the current geom. var. using finite
C             ..difference approximations.
C             Perturb deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
C             Interpolate perturbed dep. var.s ZG1 and derivs wrt Nu
              CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG1,ERROR,*9999)
C             Reset deformed element coords
              ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE
C             Calculate deformed metric tensors for perturbed coords
              CALL ZGMG(NBH(NH_LOC(1,nx)),AZ1,AZL1,AZU1,ZG1,ERROR,*9999)
C             Finite diff approx to derivs
              DO I=1,3
                DO J=1,3
                  D_AZL(i,j)=(AZL1(i,j)-AZL(i,j))/DELTA_ZE
                  D_AZU(i,j)=(AZU1(i,j)-AZU(i,j))/DELTA_ZE
                ENDDO
              ENDDO
              D_AZ=(AZ1-AZ)/DELTA_ZE
C             Calc derivs of TG wrt current geom var analytically
              IF(KTYP51(nr).EQ.3) THEN
                CALL D_ZGTG53(PARAMTYPE,NBH(NH_LOC(1,nx)),0,nr,nx,
     '            AXU,AZ,AZL,AZU,CG(1,ng),D_AZ,D_AZL,D_AZU,D_EG,
     '            D_RI3(nhs1),D_TG(1,1,nhs1),D_ZG(1,1,nhs1),EG,ZG,
     '            ERROR,*9999)
              ELSE IF(KTYP51(nr).EQ.4) THEN
                CALL D_ZGTG54(PARAMTYPE,0,nr,nx,
     '            AXU,AZ,AZL,AZU,CG(1,ng),D_AZ,D_AZL,D_AZU,D_EG,
     '            D_RI3(nhs1),D_TG(1,1,nhs1),EG,ERROR,*9999)
              ENDIF

              IF(KTYP51(nr).EQ.4.AND.NW.GT.1) THEN !membrane+press bc(s)
C               calculate this quantity here for use later in Gauss pt
C               loop (saves recalculation of D_AZU later)
                DO nhx=1,NJ_LOC(NJL_GEOM,0)
                  nh=NH_LOC(nhx,nx)
                  nb=NBH(nh)
                  D_SUM_MEMB(nh,nhs1)=0.0d0
                  DO ni=1,NIT(nb)
                    D_SUM_MEMB(nhx,nhs1)=D_SUM_MEMB(nhx,nhs1)
     '                +D_AZU(ni,3)*ZG(nhx,NU1(ni))
     '                +AZU(ni,3)*D_ZG(nhx,NU1(ni),nhs1)
                  ENDDO !ni
                ENDDO !nhx
              ENDIF

            ENDDO !ns1
          ENDDO !nhx1
        ENDIF

        RWG=RGX(ng)*WG(ng,1)
        IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*XG(1,1)    !cyl symmetry about x
        IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*XG(2,1)    !cyl symmetry about y
        IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*XG(1,1)**2 !spherical symmetry

C ***   Derivative of main element residual
        NSTB=NST(NBH(NH_LOC(1,nx)))
        NATB=NAT(NBH(NH_LOC(1,nx)))
        IF(TYPE(1:4).EQ.'RC3D') THEN     !3D rect.cart.
C         Calculate derivs wrt material params using unrolled loops
          IF(DOP) THEN
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DO njj=1,NJ_LOC(NJL_GEOM,0)
            nh=NJ_LOC(NJL_GEOM,njj)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              PPGG(2) = PG(ns,NU1(1),ng,nb)*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,1)
              PPGG(3) = PG(ns,NU1(1),ng,nb)*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,2)
              PPGG(4) = PG(ns,NU1(1),ng,nb)*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,3)
              DO noopti=1,NTOPTI
          D_AGE=(D_TG(1,1,noopti)*ZG(nh,2)+
     '           D_TG(1,2,noopti)*ZG(nh,4)+
     '           D_TG(1,3,noopti)*ZG(nh,7))*PPGG(2)
     '         +(D_TG(2,1,noopti)*ZG(nh,2)+
     '           D_TG(2,2,noopti)*ZG(nh,4)+
     '           D_TG(2,3,noopti)*ZG(nh,7))*PPGG(3)
     '         +(D_TG(3,1,noopti)*ZG(nh,2)+
     '           D_TG(3,2,noopti)*ZG(nh,4)+
     '           D_TG(3,3,noopti)*ZG(nh,7))*PPGG(4)
          D_RE(ns,nh,noopti)=D_RE(ns,nh,noopti)+
     '            D_AGE*RWG*SE(ns,nb,ne)
              ENDDO
            ENDDO
          ENDDO

        ELSE IF(TYPE(1:7).EQ.'PROLATE') THEN  !3D prolate spheroidal
C         Calculate derivs wrt material params using unrolled loops
          IF(DOP) THEN
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          SLZ=SINH(ZG(1,1))
          SMZ=SIN (ZG(2,1))
          CLZ=SQRT(1.0d0+SLZ*SLZ)
          CMZ=SQRT(1.0d0-SMZ*SMZ)
          CSLZ=CLZ/SLZ
          CSMZ=CMZ/SMZ
          G1=SLZ*SLZ+SMZ*SMZ
          E1=CLZ*SLZ/G1
          E2=CMZ*SMZ/G1

          DO ns=1,NSTB+NATB
            PPG(ns,1) = PG(ns,1,ng,NBH(NH_LOC(1,nx)))
            PPG(ns,2) = PG(ns,NU1(1),ng,NBH(NH_LOC(1,nx)))*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,NBH(NH_LOC(1,nx)))*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,NBH(NH_LOC(1,nx)))*DXIX(3,1)
            PPG(ns,3) = PG(ns,NU1(1),ng,NBH(NH_LOC(1,nx)))*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,NBH(NH_LOC(1,nx)))*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,NBH(NH_LOC(1,nx)))*DXIX(3,2)
            PPG(ns,4) = PG(ns,NU1(1),ng,NBH(NH_LOC(1,nx)))*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,NBH(NH_LOC(1,nx)))*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,NBH(NH_LOC(1,nx)))*DXIX(3,3)
          ENDDO

          D11=ZG(1,NU1(1))
          D12=ZG(1,NU1(2))
          D13=ZG(1,NU1(3))
          D21=ZG(2,NU1(1))
          D22=ZG(2,NU1(2))
          D23=ZG(2,NU1(3))
          D31=ZG(3,NU1(1))
          D32=ZG(3,NU1(2))
          D33=ZG(3,NU1(3))

          DLA1=ZG(1,NU1(1))
          DLA2=ZG(1,NU1(2))
          DLA3=ZG(1,NU1(3))
          DMA1=ZG(2,NU1(1))
          DMA2=ZG(2,NU1(2))
          DMA3=ZG(2,NU1(3))
          DTA1=ZG(3,NU1(1))
          DTA2=ZG(3,NU1(2))
          DTA3=ZG(3,NU1(3))

          DO ns=1,NSTB+NATB
            PGG =PPG(ns,1)
            PGA1=PPG(ns,2)
            PGA2=PPG(ns,3)
            PGA3=PPG(ns,4)

        DO noopti=1,NTOPTI
              D_AGE1=D_TG(1,1,noopti)*(D11*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+
     '          D21*(E1*DMA1-E2*DLA1)*PGG+D31*E1*SMZ*SMZ*DTA1*PGG) +
     '          D_TG(2,1,noopti)*(D11*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '          D21*(E1*DMA2-E2*DLA2)*PGG+D31*E1*SMZ*SMZ*DTA2*PGG) +
     '          D_TG(3,1,noopti)*(D11*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '          D21*(E1*DMA3-E2*DLA3)*PGG+D31*E1*SMZ*SMZ*DTA3*PGG) +
     '          D_TG(1,2,noopti)*(D12*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '          D22*(E1*DMA1-E2*DLA1)*PGG+D32*E1*SMZ*SMZ*DTA1*PGG) +
     '          D_TG(2,2,noopti)*(D12*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '          D22*(E1*DMA2-E2*DLA2)*PGG+D32*E1*SMZ*SMZ*DTA2*PGG) +
     '          D_TG(3,2,noopti)*(D12*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '          D22*(E1*DMA3-E2*DLA3)*PGG+D32*E1*SMZ*SMZ*DTA3*PGG) +
     '          D_TG(1,3,noopti)*(D13*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '          D23*(E1*DMA1-E2*DLA1)*PGG+D33*E1*SMZ*SMZ*DTA1*PGG) +
     '          D_TG(2,3,noopti)*(D13*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '          D23*(E1*DMA2-E2*DLA2)*PGG+D33*E1*SMZ*SMZ*DTA2*PGG) +
     '          D_TG(3,3,noopti)*(D13*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '          D23*(E1*DMA3-E2*DLA3)*PGG+D33*E1*SMZ*SMZ*DTA3*PGG)

              D_AGE2 = D_TG(1,1,noopti)*(D11*(E2*DLA1-E1*DMA1)*PGG +
     '          D21*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D31*E2*SLZ*SLZ*DTA1*
     '            PGG)+
     '          D_TG(2,1,noopti)*(D11*(E2*DLA2-E1*DMA2)*PGG +
     '          D21*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D31*E2*SLZ*SLZ*DTA2*
     '            PGG)+
     '          D_TG(3,1,noopti)*(D11*(E2*DLA3-E1*DMA3)*PGG +
     '          D21*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D31*E2*SLZ*SLZ*DTA3*
     '            PGG)+
     '          D_TG(1,2,noopti)*(D12*(E2*DLA1-E1*DMA1)*PGG +
     '          D22*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D32*E2*SLZ*SLZ*DTA1*
     '            PGG)+
     '          D_TG(2,2,noopti)*(D12*(E2*DLA2-E1*DMA2)*PGG +
     '          D22*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D32*E2*SLZ*SLZ*DTA2*
     '            PGG)+
     '          D_TG(3,2,noopti)*(D12*(E2*DLA3-E1*DMA3)*PGG +
     '          D22*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D32*E2*SLZ*SLZ*DTA3*
     '            PGG)+
     '          D_TG(1,3,noopti)*(D13*(E2*DLA1-E1*DMA1)*PGG +
     '          D23*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D33*E2*SLZ*SLZ*DTA1*
     '            PGG)+
     '          D_TG(2,3,noopti)*(D13*(E2*DLA2-E1*DMA2)*PGG +
     '          D23*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D33*E2*SLZ*SLZ*DTA2*
     '            PGG)+
     '          D_TG(3,3,noopti)*(D13*(E2*DLA3-E1*DMA3)*PGG +
     '          D23*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D33*E2*SLZ*SLZ*DTA3*
     '            PGG)

              D_AGE3 = D_TG(1,1,noopti)*(-(D11*CSLZ+D21*CSMZ)*DTA1*PGG +
     '          D31*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '          D_TG(2,1,noopti)*(-(D11*CSLZ+D21*CSMZ)*DTA2*PGG +
     '          D31*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '          D_TG(3,1,noopti)*(-(D11*CSLZ+D21*CSMZ)*DTA3*PGG +
     '          D31*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '          D_TG(1,2,noopti)*(-(D12*CSLZ+D22*CSMZ)*DTA1*PGG +
     '          D32*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '          D_TG(2,2,noopti)*(-(D12*CSLZ+D22*CSMZ)*DTA2*PGG +
     '          D32*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '          D_TG(3,2,noopti)*(-(D12*CSLZ+D22*CSMZ)*DTA3*PGG +
     '          D32*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '          D_TG(1,3,noopti)*(-(D13*CSLZ+D23*CSMZ)*DTA1*PGG +
     '          D33*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '          D_TG(2,3,noopti)*(-(D13*CSLZ+D23*CSMZ)*DTA2*PGG +
     '          D33*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '          D_TG(3,3,noopti)*(-(D13*CSLZ+D23*CSMZ)*DTA3*PGG +
     '          D33*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG))

              D_RE(ns,1,noopti)=D_RE(ns,1,noopti)
     '          +D_AGE1*RWG*SE(ns,NBH(NH_LOC(1,nx)),ne)
              D_RE(ns,2,noopti)=D_RE(ns,2,noopti)
     '          +D_AGE2*RWG*SE(ns,NBH(NH_LOC(2,nx)),ne)
              D_RE(ns,3,noopti)=D_RE(ns,3,noopti)
     '          +D_AGE3*RWG*SE(ns,NBH(NH_LOC(3,nx)),ne)

            ENDDO
          ENDDO

        ELSE
C         Calculate derivs wrt mat or geometric params for general case
          nhs=0
          DO njj1=1,NJ_LOC(NJL_GEOM,0)
            nh=NJ_LOC(NJL_GEOM,njj1)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              PPG(ns,1)=PG(ns,1,ng,nb)
              DO njj2=1,NJ_LOC(NJL_GEOM,0)
                nj=NJ_LOC(NJL_GEOM,njj2)
                PPG(ns,1+nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
              ENDDO
            ENDDO
           IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
              DO ns=1,NST(nb)+NAT(nb)
                DO noopti=1,NTOPTI
            D_AGE=D_AG(PARAMTYPE,nb,nh,nr,ns,D_TG(1,1,noopti),
     '              D_ZG(1,1,noopti),PPG,TG,ZG)
            D_RE(ns,nh,noopti)=D_RE(ns,nh,noopti)+
     '              D_AGE*RWG*SE(ns,nb,ne)
                ENDDO
              ENDDO
            ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
              DO ns=1,NST(nb)+NAT(nb)
                nhs=nhs+1
                nhs1=0
                DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  NB1=NBH(nh1)
                  DO ns1=1,NST(NB1)+NAT(NB1)
                    nhs1=nhs1+1
                    D_AGE=D_AG(PARAMTYPE,nb,nh,nr,ns,D_TG(1,1,nhs1),
     '                D_ZG(1,1,nhs1),PPG,TG,ZG)
                    ES(nhs,nhs1)=ES(nhs,nhs1)+D_AGE*RWG*SE(ns,nb,ne)
             ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          NHSGMAX=nhs   !max row # for geom resids
        ENDIF

C ***   Incompressibilty (+fluid) constraints
        IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).GE.2) THEN
          NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for press vars
          IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
            DO ns=1,NST(NBP)+NAT(NBP)
              DO noopti=1,NTOPTI
                D_RE(ns,NH_LOC(NH_LOC(0,nx),nx),noopti)=0.0d0
              ENDDO !noopti
            ENDDO !ns
          ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
            IF(KTYP52(nr).EQ.2) THEN       !Incompressible only
              DO ns=1,NST(NBP)+NAT(NBP)
                nhs=nhs+1
                nhs1=0
                DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  NB1=NBH(nh1)
                  DO ns1=1,NST(NB1)+NAT(NB1)
                    nhs1=nhs1+1
                    IF(NODEPRESS.AND.ELEMPRESS) THEN
C                     node and element based hydrostatic press interp
C                     NOTE: pressure basis function is not used to
C                           weight the residual here (non-Galerkin)
                      ES(nhs,nhs1)=ES(nhs,nhs1)+D_RI3(nhs1)/
     '                  (2.0d0*DSQRT(RI3))*RWG
                    ELSE
C                     node (only) or element (only) based hyd. pressure
C                     interp (standard Galerkin)
                      ES(nhs,nhs1)=ES(nhs,nhs1)+D_RI3(nhs1)/
     '                  (2.0d0*DSQRT(RI3))*PG(ns,1,ng,NBP)*RWG
                    ENDIF
                  ENDDO !ns1
                ENDDO !nhx1
              ENDDO !ns

            ELSE IF(KTYP52(nr).EQ.3) THEN !Incompressible + fluid
C             Calc deformed Christoffel symbols wrt undef Nu coords
C             NOTE: ZG needs derivs wrt Nu not Xi !
              CALL TOFFEL(NBJ(1),NJE,nr,CHTOFF,DBM,AZU,ZG,X3G,.FALSE.,
     '          ERROR,*9999)
C             Calc the Jacobian for integration wrt def coords:
C             Interpolate dependent var.s ZG and derivs wrt Xi
              CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C             Calculate deformed metric tensors wrt Xi (GZL,GZU)
              CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)
C             Calc deformed Jacobian*Gaussian Quadrature weight
              RZWG=DSQRT(GZ)*WG(ng,NBH(NH_LOC(1,nx)))
              IF(JTYP4.EQ.2) RZWG=RZWG*2.0d0*PI*ZG(1,1)    !x cyl sym
              IF(JTYP4.EQ.3) RZWG=RZWG*2.0d0*PI*ZG(2,1)    !y cyl sym
              IF(JTYP4.EQ.4) RZWG=RZWG*4.0d0*PI*ZG(1,1)**2 !sph sym

              DO ns=1,NST(NBP)+NAT(NBP)
                nhs=nhs+1
                nhs1=0
                DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  NB1=NBH(nh1)
                  DO ns1=1,NST(NB1)+NAT(NB1)
                    nhs1=nhs1+1
C                   Perturb deformed element coords
                    ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
C                   Interp. perturbed dep. var.s ZG1 and derivs wrt Nu
                    CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG1,
     '                ERROR,*9999)
C                   Calc deformed metric tensors at perturbed coords
                    CALL ZGMG(NBH(NH_LOC(1,nx)),AZ1,AZL1,AZU1,ZG1,
     '                ERROR,*9999)
C                   Calculate deformed Christoffel symbols wrt undef
C                   Nu coords at perturbed geom coords
C                   NOTE: ZG1 needs derivs wrt Nu not Xi!
                    CALL TOFFEL(NBJ(1),NJE,nr,CHTOFF1,DBM,AZU,ZG1,X3G,
     '                .FALSE.,ERROR,*9999)
C                   Calc the deriv (wrt current geom var)
C                   of the Jacobian for integration wrt def coords:
C                   Interp. perturbed dep. var.s ZG1 and derivs wrt Xi
                    CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG1,
     '                ERROR,*9999)
C                   Calculate deformed metric tensors wrt Xi at
C                   perturbed coords (GZL1,GZU1)
                    CALL ZGMG(NBH(NH_LOC(1,nx)),GZ1,GZL1,GZU1,ZG1,
     '                ERROR,*9999)
C                   Calc deformed Jacobian*Gaussian Quadrature weight
                    RZWG1=DSQRT(GZ1)*WG(ng,NBH(NH_LOC(1,nx)))
                    IF(JTYP4.EQ.2)RZWG1=RZWG1*2.d0*PI*ZG(1,1) !x cyl sym
                    IF(JTYP4.EQ.3)RZWG1=RZWG1*2.d0*PI*ZG(2,1) !y cyl sym
                    IF(JTYP4.EQ.4)RZWG1=RZWG1*4.d0*PI*ZG(1,1)**2!sph sym
C                   Reset deformed element coords
                    ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE

C                   Finite diff approx to derivs D_AZU,D_CHTOFF,D_RZWG
                    DO I=1,3
                      DO J=1,3
                        D_AZU(i,j)=(AZU1(i,j)-AZU(i,j))/DELTA_ZE
                        DO K=1,3
                          D_CHTOFF(i,j,k)=(CHTOFF1(i,j,k)-CHTOFF(i,j,k))
     '                      /DELTA_ZE
                        ENDDO
                      ENDDO
                    ENDDO
                    D_RZWG=(RZWG1-RZWG)/DELTA_ZE

C                   Calc deriv wrt current geom var of del-squared(p)
C                   (wrt def coords), where p is the hyd press
C                   that varies with Xi3 only.
                    SUM=0.0d0
                    D_SUM=0.0d0
                    DO J=1,NITB
                      DO K=1,NITB
                        SUM=SUM+CHTOFF(3,j,k)*AZU(j,k)
                        D_SUM=D_SUM+D_CHTOFF(3,j,k)*AZU(j,k)
     '                             +CHTOFF(3,j,k)*D_AZU(j,k)
                      ENDDO
                    ENDDO
                    delsqP=0.0d0
                    D_delsqP=0.0d0
                    DO nsa=1,NST(NBP)+NAT(NBP)
                      delsqP=delsqP+ZE(nsa,NH_LOC(0,nx))*
     '                  (PG(nsa,8,ng,NBP)*AZU(3,3)-SUM*PG(nsa,7,ng,NBP))
                      D_delsqP=D_delsqP+ZE(nsa,NH_LOC(0,nx))*
     '                  (PG(nsa,8,ng,NBP)*D_AZU(3,3)-
     '                  D_SUM*PG(nsa,7,ng,NBP))
                    ENDDO
C                   Check if taking deriv wrt the current ZE aux var
                    IF(nh1.EQ.NH_LOC(NH_LOC(0,nx),nx)) THEN
                      D_delsqP=D_delsqP+
     '                  (PG(ns1,8,ng,NBP)*AZU(3,3)-SUM*PG(ns1,7,ng,NBP))
                    ENDIF

C                   Calc derivs of resids associated with Darcy's Law
                    Darcy_Resid=CG(IL_fluid_conductivity,ng)
     '                *DT*delsqP+(DSQRT(RI3)-1.d0)/DSQRT(RI3)
                    D_Darcy_Resid=CG(IL_fluid_conductivity,ng)
     '                *DT*D_delsqP+D_RI3(nhs1)
     '                /(2.0d0*DSQRT(RI3*RI3*RI3))
                    IF(NODEPRESS.AND.ELEMPRESS) THEN
C                     node and element based hydrostatic pressure interp
C                     NOTE: pressure basis function is not used to
C                           weight the residual here (non-Galerkin)
                      ES(nhs,nhs1)=ES(nhs,nhs1)
     '                  +D_Darcy_Resid*RZWG
     '                  +Darcy_Resid*D_RZWG
                    ELSE
C                     node (only) or element (only) based hyd. pressure
C                     interp (standard Galerkin)
                      ES(nhs,nhs1)=ES(nhs,nhs1)
     '                  +(D_Darcy_Resid*RZWG+Darcy_Resid*D_RZWG)
     '                  *PG(ns,1,ng,NBP)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

        ELSE IF(KTYP51(nr).EQ.4.AND.NW.GT.1) THEN !membrane+press bc(s)
          IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
            nhs=0
            DO njj=1,NJ_LOC(NJL_GEOM,0)
              nj=NJ_LOC(NJL_GEOM,njj)
              nb=NBH(nj)
C             NOTE: deriv of SUM, D_SUM_MEMB(nj,nhs1) was calc'ed above
              DO ns=1,NST(nb)+NAT(nb)
                nhs=nhs+1
                nhs1=0
                DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  NB1=NBH(nh1)
                  DO ns1=1,NST(NB1)+NAT(NB1)
                    nhs1=nhs1+1
                    ES(nhs,nhs1)=ES(nhs,nhs1)
     '                +PF(1)*D_SUM_MEMB(nj,nhs1)*PG(ns,1,ng,nb)*RWG
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

        ENDIF

 50   CONTINUE

C *** External pressure loads in 3D case
      IF(KTYP51(nr).EQ.3.AND.NW.GE.2.AND.NW.LE.4) THEN
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for hydr. pressure
C       Boundary pressure constraint equations
        IF(KTYP52(nr).GE.2.AND.KTYP57(nr).GT.1) THEN
C         Boundary pressure constraint equations
C         NOTE: if hydrostatic pressure in nodally (only) based
C               then there are no boundary pressure constraint eqns
          IF(NODEPRESS.AND.ELEMPRESS) THEN
C           node and element based hydrostatic pressure interp
            ERROR='>>Derivatives not implemented'
            GOTO 9999
          ELSE IF(ELEMPRESS) THEN
C           element (only) based hydrostatic pressure interp
            CALL D_PFRE_NE(PARAMTYPE,IBT,IDO,INP,
     '        NAN,NBH,NBJ,NFF,NGAP,NHE,NHSGMAX,NMNO,NPF,
     '        NPNE(1,1,ne),nr,NRE,NSP,NW,nx,NXI(-NIM,ne),
     '        CE,CP,D_RE,D_RI3,D_TG,D_ZE,D_ZG,ES,FEXT,XE,XG,ZE,ZG,ZG1,
     '        ERROR,*9999)
          ENDIF
        ENDIF

C       For pressure loading need to add contributions to the geometric
C       residuals. Note that this needs to be done only for the
C       derivatives wrt geometric parameters as these contributions are
C       independent of the material parameters
        IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
C         External pressure loading
          DO iface=1,2
            IF(iface.EQ.1.AND.(NXI(-3,ne).EQ.0.OR.nr.NE.NRE(NXI(-3,ne)))
     '        .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
     '        iface.EQ.2.AND.(NXI(3,ne).EQ.0.OR.nr.NE.NRE(NXI(3,ne)))
     '        .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc on Xi3=1 face
              IF(DABS(PF(iface)).GT.1.0D-10) THEN
                IXF=iface-1
                IFE=iface+4
                nf=NFF(IFE)
                nhs1=0
           DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  NB1=NBH(nh1)
             DO ns1=1,NST(NB1)+NAT(NB1)
                    nhs1=nhs1+1
C                   Analytic derivs of face press resids wrt current
C                   geom variable
                    CALL D_PFRF(IBT,IDO,INP,IXF,
     '                NAN,NBH,NBJ,NGAP,NHE,nhx1,NPF(1,nf),nr,ns1,nx,
     '                D_RF,D_ZE,D_ZG(1,1,nhs1),PF(iface),PG,
     '                WG,XE,XG,ZE,ZG,ZG1,ERROR,*9999)
                    NHSTART=0
                    DO nhx=1,NJ_LOC(NJL_GEOM,0)
                      nh=NH_LOC(nhx,nx)
                      nsf=0
                      NBFF=NPF(9+nh,nf)
                      NBE=NBH(nh)
                      DO nnbf=1,NNT(NBFF)
                        nn=NNF(1+nnbf,IFE,NBE)
                        DO nkbf=1,NKT(nnbf,NBFF)
                          nk=NKF(nkbf,nnbf,IFE,NBE)
                          nse=nk+(nn-1)*NKT(nn,NBE)
                          nhs=NHSTART+nk+(nn-1)*NKT(nn,NBE)
                          nsf=nsf+1
                          ES(nhs,nhs1)=ES(nhs,nhs1)
     '                      -D_RF(nsf,nh)*SE(nse,NBE,ne)
                        ENDDO !nkbf (nk)
                      ENDDO !nnbf (nn)
                      DO nn=1,NNT(NBE)
                        NHSTART=NHSTART+NKT(nn,NBE)
                      ENDDO !nn
                    ENDDO !nhx
                  ENDDO !ns1
                ENDDO !nhx1
              ENDIF
            ENDIF
          ENDDO !iface
        ENDIF
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' >>D_ZERE50 diagnostic output'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
          WRITE(OP_STRING,'(/('' D_RE('',I2,'','',I1,'
     '      //''',noopti): '',5(1X,D12.4),/(20X,5(1X,D12.4))))')
     '      ((nh,(ns,D_RE(ns,nh,noopti),noopti=1,NTOPTI),
     '      ns=1,NST(NBH(NH_LOC(nhx,nx)))
     '      +NAT(NBH(NH_LOC(nhx,nx)))),nhx=1,NH_LOC(0,nx))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
          nhs=0
          DO nhx=1,NH_LOC(0,nx)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              nhs=nhs+1
              WRITE(OP_STRING,'(/'' ES('',I3,'
     '          //''',nhs1): '',5(1X,D12.4),/(15X,5(1X,D12.4)))')
     '          nhs,((ES(nhs,ns1+(nh1-1)*(NST(NBH(NH_LOC(nhx1,nx)))
     '          +NAT(NBH(NH_LOC(nhx1,nx))))),
     '          ns1=1,NST(NBH(NH_LOC(nhx1,nx)))+
     '          NAT(NBH(NH_LOC(nhx1,nx)))),nhx1=1,NH_LOC(0,nx))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('D_ZERE50')
      RETURN
 9999 CALL ERRORS('D_ZERE50',ERROR)
      CALL EXITS('D_ZERE50')
      RETURN 1
      END


      SUBROUTINE D_ZGTG53(PARAMTYPE,nb,NMNO,nr,nx,
     '  AXU,AZ,AZL,AZU,CG,D_AZ,D_AZL,D_AZU,D_EG,D_RI3,D_TG,D_ZG,
     '  EG,ZG,ERROR,*)

C#### Subroutine: D_ZGTG53
C###  Description:
C###    D_ZGTG53 evaluates derivatives of components of 2nd
C###    Piola-Kirchhoff stress tensor TG wrt material parameters or
C###    geometric variables at current Gauss point for 3-dimensional
C###    problems (ktyp51(nr)=3).

C**** The following code was copied from ZGTG53 and altered and should
C**** be kept in sych with ZGTG53.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER nb,NMNO,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),D_AZ,D_AZL(3,3),
     '  D_AZU(3,3),D_EG(3,3),D_RI3,D_TG(3,3),D_ZG(NHM,NUM),EG(3,3),
     '  ZG(NHM,NUM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER i,j,m,n
      REAL*8 AXL(3,3),AZL_tmp(3,3),D_DW(6),D_W3,Fgrowth(3,3),W3

      CALL ENTERS('D_ZGTG53',*9999)

      CALL ASSERT(KTYP55(nr).EQ.3,'Stresses must be a function of'
     '  //' fibre and transverse strains in the constitutive law',
     '  ERROR,*9999)
      CALL ASSERT(KTYP53(nr).GT.1,'>>Stresses must be wrt body/fibre '
     '  //'coords',ERROR,*9999)

C     Compute undeformed metric tensors
      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
          AXL(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
        AXL(i,i)=1.0d0
      ENDDO

C!!! NOTE can only use resid strains for pole zero law until init extn
C!!!      mat params are set up for other problem types
      IF(KTYP56(nr).EQ.3) THEN !pole zero law
C       Calc growth defm tens for resid strain and copy AZL to AZL_tmp
        DO i=1,3
          DO j=1,3
            Fgrowth(i,j)=0.0d0 !growth defm tensor for resid strains
            AZL_tmp(i,j)=AZL(i,j)
          ENDDO !j
          Fgrowth(i,i)=CG(27+i) !NOTE init extns CG(28),CG(29),CG(30)
        ENDDO !i
C       Apply growth defm to deformed covariant metric tensor AZL
        DO i=1,3
          DO j=1,3
            AZL(i,j)=0.0d0
            DO m=1,3
              DO n=1,3
                AZL(i,j)=AZL(i,j)+
     '            Fgrowth(i,m)*AZL_tmp(m,n)*Fgrowth(n,j)
              ENDDO !n
            ENDDO !m
          ENDDO !j
        ENDDO !i
C       Recompute AZU,AZ from transformed AZL
        CALL INVERT(NIT(nb),AZL,AZU,AZ)
      ENDIF !pole zero law
C     Compute Green strain tensor wrt material fibre coords
      DO i=1,3
        DO j=1,3
          EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
        ENDDO !j
      ENDDO !i
      IF(DOP) THEN
        WRITE(OP_STRING,'('' EG: '',3D12.4,/(5X,3D12.4))')
     '    ((EG(i,j),j=1,3),i=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        D_RI3=0.0d0
C       Calc derivs of strain energy fn wrt material parameters
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)
C       Calculate derivs of TG wrt material parameters
        DO i=1,3
          D_TG(i,i)=D_DW(i)/AXL(i,i)
        ENDDO
        D_TG(1,2)=D_DW(4)/DSQRT(AXL(1,1)*AXL(2,2))
        D_TG(1,3)=D_DW(5)/DSQRT(AXL(1,1)*AXL(3,3))
        D_TG(2,3)=D_DW(6)/DSQRT(AXL(2,2)*AXL(3,3))
        D_TG(2,1)=D_TG(1,2)
        D_TG(3,1)=D_TG(1,3)
        D_TG(3,2)=D_TG(2,3)

      ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN

        IF(KTYP56(nr).EQ.3) THEN !pole zero law
          CALL ASSERT(CG(28).EQ.1.0d0.AND.CG(29).EQ.1.0d0.AND.
     '      CG(30).EQ.1.0d0,'>>Needs updating for growth tensor',
     '      ERROR,*9999)
        ENDIF

C       Set the derivative of the third strain invariant wrt the
C       current geometric variable
        D_RI3=D_AZ
C       Calculate the derivative of the strain matrix, D_EG, wrt the
C       current geometric variable at the current Gauss point
C       Here D_AZL is the derivative of deformed metric tensor
C       (Nu coords) wrt the current geometric variable
        DO i=1,3
          DO j=1,3
            D_EG(i,j)=0.50d0*D_AZL(i,j)/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,'('' D_EG: '',9D12.4)')
     '      ((D_EG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

C       Get deriv's of strain energy fn wrt fibre strains
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)

        IF(KTYP52(nr).EQ.1) THEN !compressible
          W3=0.0d0   !Check this
          D_W3=0.0d0 !Check this
        ELSE                 !incompressible
          W3=ZG(NH_LOC(0,nx),1)
          D_W3=D_ZG(NH_LOC(0,nx),1)
        ENDIF

C       Get deriv's of stress tensor wrt geometric variables
        DO i=1,3
          D_TG(i,i)=D_DW(i)*D_EG(i,i)/AXL(i,i)
     '            + 2.0d0*D_W3*AZU(i,i) + 2.0d0*W3*D_AZU(i,i)
        ENDDO
        D_TG(1,2)=D_DW(4)*D_EG(1,2)/DSQRT(AXL(1,1)*AXL(2,2))
     '            + 2.0d0*D_W3*AZU(1,2) + 2.0d0*W3*D_AZU(1,2)
        D_TG(1,3)=D_DW(5)*D_EG(1,3)/DSQRT(AXL(1,1)*AXL(3,3))
     '            + 2.0d0*D_W3*AZU(1,3) + 2.0d0*W3*D_AZU(1,3)
        D_TG(2,3)=D_DW(6)*D_EG(2,3)/DSQRT(AXL(2,2)*AXL(3,3))
     '            + 2.0d0*D_W3*AZU(2,3) + 2.0d0*W3*D_AZU(2,3)
        D_TG(2,1)=D_TG(1,2)
        D_TG(3,1)=D_TG(1,3)
        D_TG(3,2)=D_TG(2,3)
      ENDIF

      CALL EXITS('D_ZGTG53')
      RETURN
 9999 CALL ERRORS('D_ZGTG53',ERROR)
      CALL EXITS('D_ZGTG53')
      RETURN 1
      END


      SUBROUTINE D_ZGTG54(PARAMTYPE,NMNO,nr,nx,
     '  AXU,AZ,AZL,AZU,CG,D_AZ,D_AZL,D_AZU,D_EG,D_RI3,D_TG,
     '  EG,ERROR,*)

C#### Subroutine: D_ZGTG54
C###  Description:
C###    D_ZGTG54 evaluates derivatives of components of 2nd
C###    Piola-Kirchhoff stress tensor TG wrt material parameters or
C###    geometric variables at current Gauss point for membrane
C###    problems (ktyp51(nr)=4).

C**** The following code was copied from ZGTG54 and altered and should
C**** be kept in sych with ZGTG54.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
!     Parameter List
      INTEGER NMNO,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),D_AZ,D_AZL(3,3),
     '  D_AZU(3,3),D_EG(3,3),D_RI3,D_TG(3,3),EG(3,3)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER i,ie,j,je
      REAL*8 AXL(3,3),D_DW(6),D_H,DW(6),H

      CALL ENTERS('D_ZGTG54',*9999)

      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
          AXL(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
        AXL(i,i)=1.0d0
      ENDDO

      CALL ASSERT(KTYP55(nr).EQ.3,'Stresses must be a function of'
     '  //' fibre  and transverse strains in the constitutive law',
     '  ERROR,*9999)

      IF(KTYP53(nr).EQ.1) THEN      !stresses in theta coords
        ERROR='>>Stresses must be wrt body/fibre coords'
        GO TO 9999
      ELSE IF(KTYP53(nr).GT.1) THEN !stresses in body/fibre Nu coords
        AZL(3,3)=1.0d0/AZ !AZ=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        AZU(3,3)=AZ
        AZL(1,3)=0.0d0
        AZL(2,3)=0.0d0
        AZL(3,1)=0.0d0
        AZL(3,2)=0.0d0
        IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
          D_AZL(3,3)=0.0d0
          D_AZU(3,3)=0.0d0
        ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
          D_AZL(3,3)=-D_AZ/(AZ*AZ)
          D_AZU(3,3)=D_AZ
        ENDIF
        D_AZL(1,3)=0.0d0
        D_AZL(2,3)=0.0d0
        D_AZL(3,1)=0.0d0
        D_AZL(3,2)=0.0d0
        D_RI3=0.0d0  !wrt both material and geometric parameters
      ENDIF

      !Calculate strain matrix, EG, at current Gauss point
      !Here AZL is the deformed metric tensor wrt Nu coords
      DO ie=1,3
        DO je=1,3
          EG(ie,je)=0.50d0*(AZL(ie,je)-AXL(ie,je))/
     '      DSQRT(AXL(ie,ie)*AXL(je,je))
        ENDDO
      ENDDO
      IF(DOP) THEN
        WRITE(OP_STRING,'('' EG: '',9D12.4)')
     '    ((EG(ie,je),je=1,3),ie=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        !Calc derivs of strain energy fn wrt material parameters
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)
        !Calculate derivs of TG wrt material parameters
        IF(KTYP52(nr).EQ.1) THEN !compressible
          !22Jan89: Compressible case should also be included here
        ELSE                 !incompressible
          D_H=-D_DW(3)/AZU(3,3)
        ENDIF
        D_TG(1,1)=D_DW(1)/AXL(1,1)                +D_H*AZU(1,1)
        D_TG(2,2)=D_DW(2)/AXL(2,2)                +D_H*AZU(2,2)
        D_TG(1,2)=D_DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+D_H*AZU(1,2)
        D_TG(2,1)=D_TG(1,2)
        D_TG(1,3)=0.0d0
        D_TG(2,3)=0.0d0
        D_TG(3,1)=0.0d0
        D_TG(3,2)=0.0d0

      ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
        !Calculate the derivative of the strain matrix, D_EG, wrt the
        !current geometric variable at the current Gauss point
        !Here D_AZL is the derivative of deformed metric tensor
        !(Nu coords) wrt the current geometric variable
        DO ie=1,3
          DO je=1,3
            D_EG(ie,je)=0.50d0*D_AZL(ie,je)/DSQRT(AXL(ie,ie)*AXL(je,je))
          ENDDO
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,'('' D_EG: '',9D12.4)')
     '      ((D_EG(ie,je),je=1,3),ie=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '    EG(2,3),ERROR,*9999)

        !Get deriv's of strain energy fn wrt fibre strains
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)

        IF(KTYP52(nr).EQ.1) THEN !compressible
          !22Jan89: Compressible case should also be included here
        ELSE                 !incompressible
          H= -DW(3)/AZU(3,3)
          D_H= -D_DW(3)*D_EG(3,3)/AZU(3,3)
     '         +DW(3)*D_AZU(3,3)/(AZU(3,3)*AZU(3,3))
        ENDIF

        !Get deriv's of stress tensor wrt geometric variables
        D_TG(1,1)=D_DW(1)*D_EG(1,1)/AXL(1,1)
     '    + D_H*AZU(1,1) + H*D_AZU(1,1)
        D_TG(2,2)=D_DW(2)*D_EG(2,2)/AXL(2,2)
     '    + D_H*AZU(2,2) + H*D_AZU(2,2)
        D_TG(1,2)=D_DW(4)*D_EG(1,2)/DSQRT(AXL(1,1)*AXL(2,2))
     '    + D_H*AZU(1,2) + H*D_AZU(1,2)
        D_TG(2,1)=D_TG(1,2)
        D_TG(1,3)=0.0d0
        D_TG(2,3)=0.0d0
        D_TG(3,1)=0.0d0
        D_TG(3,2)=0.0d0
      ENDIF

      CALL EXITS('D_ZGTG54')
      RETURN
 9999 CALL ERRORS('D_ZGTG54',ERROR)
      CALL EXITS('D_ZGTG54')
      RETURN 1
      END


      SUBROUTINE ENERGY(nr,CG,DW,P1,P2,P3,P4,P5,P6,ERROR,*)

C#### Subroutine: ENERGY
C###  Description:
C###    ENERGY calculates derivatives of strain energy function wrt
C###    either principal strain invariants (KTYP55(nr)=1), or
C###    principal extensions (KTYP55(nr)=2), or physical strains
C###    (KTYP55(nr)=3) at current Gauss point.
C**** For KTYP55(nr)=1:
C****      P1 is the First  principal invariant RI1;
C****      P2 is the Second principal invariant RI2;
C****      P3 is the Third  principal invariant RI3;
C****      P4 is the First  transverse isotropic invariant RK1;
C****      P5 is the Second transverse isotropic invariant RK2;
C****      DW(1..5) are dW/dI1,dW/dI2,dW/dI3,dW/dK1,dW/dK2.
C**** For KTYP55(nr)=2:
C****      P1 is the First  principal extension ratio RL1;
C****      P2 is the Second principal extension ratio RL2;
C****      P3 is the Third  principal extension ratio RL3;
C****      DW(1..3) are dW/dL1,dW/dL2,dW/dL3.
C**** For KTYP55(nr)=3:
C****      P1..P6 are the physical components of Green's strain wrt
C****      theta (KTYP53(nr)=1) or nu coords (KTYP53(nr)>1);
C****      E(1,1),E(2,2),E(3,3),E(1,2),E(1,3) and E(2,3);
C****      DW(1..6) are dW/dE(1,1) .. dW/dE(2,3).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
!     Parameter List
      INTEGER nr
      REAL*8 CG(NMM),DW(6),P1,P2,P3,P4,P5,P6
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i1,i2,i3,il,iterm
      REAL*8 alpha1,alpha2,coeff1,coeff2,CEXPQ,DENOM1,DENOM2,
     '  DF1,DF2,DF3,DtoAlpha1,DtoAlpha2,DWslope1,DWslope2,
     '  FI1,FI2,FI3,pole1,pole2,
     '  R1,R2,R3,RI1,RI2,RI3,RK1,RL1,RL2,RL3,
     '  strain,strain1,strain2,TOL,ZERO
      PARAMETER(ZERO=0.0d0)

      CALL ENTERS('ENERGY',*9999)
      IF(KTYP55(nr).EQ.1) THEN      !strain invariants
        RI1=P1
        RI2=P2
        RI3=P3
        RK1=P4
        R1=RI1-3.0d0
        R2=RI2-3.0d0
        R3=RI3-1.0d0
      ELSE IF(KTYP55(nr).EQ.2) THEN !principal extensions
        RL1=P1
        RL2=P2
        RL3=P3
        R1=RL1-1.0d0
        R2=RL2-1.0d0
        R3=RL3-1.0d0
      ELSE IF(KTYP55(nr).EQ.3) THEN !fibre and transverse strains
        R1=P1
        R2=P2
        R3=P3
      ENDIF

      IF(KTYP56(nr).EQ.1) THEN      !W is polynomial
        il=0
        DW(1)=0.0d0
        DW(2)=0.0d0
        DW(3)=0.0d0
        DW(4)=0.0d0
        DW(5)=0.0d0
        DO i3=0,IT(3,nr)
          FI3=1.0d0
          DF3=1.0d0
          IF(i3.GT.0) FI3=R3**i3
          IF(i3.GT.1) DF3=DBLE(i3)*R3**(i3-1)
          DO i2=0,IT(2,nr)
            FI2=1.0d0
            DF2=1.0d0
            IF(i2.GT.0) FI2=R2**i2
            IF(i2.GT.1) DF2=DBLE(i2)*R2**(i2-1)
            DO i1=0,IT(1,nr)
              FI1=1.0d0
              DF1=1.0d0
              IF(i1.GT.0) FI1=R1**i1
              IF(i1.GT.1) DF1=DBLE(i1)*R1**(i1-1)
              il=il+1
              IF(i1.GT.0) DW(1)=DW(1)+CG(il)*DF1*FI2*FI3
              IF(i2.GT.0) DW(2)=DW(2)+CG(il)*FI1*DF2*FI3
              IF(i3.GT.0) DW(3)=DW(3)+CG(il)*FI1*FI2*DF3
            ENDDO
          ENDDO
        ENDDO

      ELSE IF(KTYP56(nr).EQ.2) THEN !W is special function
        IF(KTYP55(nr).EQ.1) THEN      !W is a func of I1,I2,I3
          IF(KTYP52(nr).EQ.1) THEN      !W is Blatz-Ko material
            DW(1)=CG(1)/RI3
            DW(2)=0.0d0
            DW(3)=-CG(1)*(RI1+2.0d0)/(RI3*RI3)
          ELSE                      !W is Exponential in I1,I2 (Demiray)
            DW(1)=CG(1)*CG(2)*DEXP(CG(2)*(RI1-3.0d0))
            DW(2)=0.0d0
          ENDIF
        ELSE IF(KTYP55(nr).EQ.2) THEN !func of princ extn ratios (Ogden)
          IF(KTYP52(nr).EQ.1) THEN
            DW(1)=CG(1)*P1**(CG(2)-1)
            DW(2)=CG(1)*P2**(CG(2)-1)
            DW(3)=CG(1)*P3**(CG(2)-1)
          ELSE
            DW(1)=CG(1)*P1**(CG(2)-1)
            DW(2)=CG(1)*P2**(CG(2)-1)
          ENDIF
        ELSE IF(KTYP55(nr).EQ.3) THEN !W=func of fibre and trans strains
          TOL=1.0D-08         !trans isotrop expon law W=Cexp(Q) (Fung)
          IF((DABS(P1).LT.TOL).AND.(DABS(P2).LT.TOL).AND.
     '      (DABS(P3).LT.TOL).AND.(DABS(P4).LT.TOL).AND.
     '      (DABS(P5).LT.TOL).AND.(DABS(P6).LT.TOL)) THEN
            CEXPQ=CG(1)
          ELSE
      CEXPQ=CG(1)*DEXP(2.0d0*CG(2)*(P1+P2+P3)
     '                    +     CG(3)* P1*P1
     '                    +     CG(4)*(P2*P2+P3*P3+2.0d0*P6*P6)
     '                    + 2.0d0*CG(5)*(P4*P4+P5*P5))
!MPN 4-Feb-1994: replaced old exponential law for D. Bloomgarden
c        CEXPQ1=CG(1)*DEXP(CG(2)*P1*P1+CG(3)*(P4*P4+P5*P5))
c        CEXPQ2=CG(4)*DEXP(CG(5)*(P2+P3)**2+CG(6)*(P2*P3-P6*P6))
          ENDIF
          DW(1)= CEXPQ*(CG(2)+CG(3)*P1)
          DW(2)= CEXPQ*(CG(2)+CG(4)*P2)
          DW(3)= CEXPQ*(CG(2)+CG(4)*P3)
          DW(4)= CEXPQ*CG(5)*P4
          DW(5)= CEXPQ*CG(5)*P5
          DW(6)= CEXPQ*CG(4)*P6
!MPN 4-Feb-1994: replaced old exponential law for D. Bloomgarden
c          DW(1)= CEXPQ1*2.0d0*CG(2)*P1
c          DW(2)= CEXPQ2*(2.0d0*CG(5)*(P2+P3)+CG(6)*P3)
c          DW(3)= CEXPQ2*(2.0d0*CG(5)*(P2+P3)+CG(6)*P2)
c          DW(4)= CEXPQ1*CG(3)*2.0d0*P4
c          DW(5)= CEXPQ1*CG(3)*2.0d0*P5
c          DW(6)=-CEXPQ2*CG(6)*2.0d0*P6
        ENDIF

      ELSE IF(KTYP56(nr).EQ.3) THEN   !W is special function
        IF(KTYP55(nr).EQ.1) THEN        !W is a func of I1,I2,I3
          IF(KTYP52(nr).EQ.1) THEN        !W is
            DW(1)=0.0d0
            DW(2)=0.0d0
            DW(3)=0.0d0
          ELSE                    !W is
            DW(1)=0.0d0
            DW(2)=0.0d0
          ENDIF
        ELSE IF(KTYP55(nr).EQ.2) THEN !func of princ extn ratios (Ogden)
          IF(KTYP52(nr).EQ.1) THEN
            DW(1)=0.0d0
            DW(2)=0.0d0
            DW(3)=0.0d0
          ELSE
            DW(1)=0.0d0
            DW(2)=0.0d0
          ENDIF
        ELSE IF(KTYP55(nr).EQ.3) THEN !func of fibre and transv strains
C         pole-zero law - orthotropic model (incl shear terms)
C old MPN 13-Apr-96: init extns handled by 'growth' defm tensor
C                    see ZGTG53.
C          L0_fibre=CG(28)               !initial fibre ext ratio
C          L0_sheet=CG(29)               !initial sheet ext ratio
C          L0_sheetnormal=CG(30)         !initial sheetnormal ext ratio
C          E0_fibre=0.5d0*(L0_fibre*L0_fibre-1.d0) !initial fibre strain
C          E0_sheet=0.5d0*(L0_sheet*L0_sheet-1.d0) !initial cross strain
C          E0_sheetnormal=0.5d0*(L0_sheetnormal*L0_sheetnormal-1.d0) !initial sheet-normal strain
C          PP1=P1+E0_fibre
C          PP2=P2+E0_sheet
C          PP3=P3+E0_sheetnormal

C         Calculate the 2nd deriv of the pole-zero law
C         for strains between 0 and 90% of the pole for each term.
C         For strains beyond 90% of their pole the stress/strain is
C         flat with magnitude equal to the stress at a strain of 90%
C         of the pole.
C         For compressive strains use linear stress/strain relationship
C         with stiffness equal to the 2nd deriv of W wrt e
C         evaluated at e=0
          TOL=0.90d0

          DO iterm=1,6
            coeff1=CG((iterm-1)*3+1) !coefficient for current term
            pole1 =CG((iterm-1)*3+2) !pole for current term
            alpha1=CG((iterm-1)*3+3) !curvature for current term
            IF(iterm.EQ.1.OR.iterm.EQ.2.OR.iterm.EQ.3) THEN  !Axial
              IF(iterm.EQ.1) THEN
                strain=P1 !=EG(1,1)
              ELSE IF(iterm.EQ.2) THEN
                strain=P2 !=EG(2,2)
              ELSE IF(iterm.EQ.3) THEN
                strain=P3 !=EG(3,3)
              ENDIF
              coeff2=coeff1
              pole2 =pole1
              alpha2=alpha1
            ELSE IF(iterm.EQ.4.OR.iterm.EQ.5.OR.iterm.EQ.6) THEN  !Shear
              IF(iterm.EQ.4) THEN
                strain=P4 !=EG(1,2)
              ELSE IF(iterm.EQ.5) THEN
                strain=P5 !=EG(1,3)
              ELSE IF(iterm.EQ.6) THEN
                strain=P6 !=EG(2,3)
              ENDIF
              coeff2=CG((iterm-1)*3+10)
              pole2 =CG((iterm-1)*3+11)
              alpha2=CG((iterm-1)*3+12)
            ENDIF

            strain1=strain
            strain2=strain
            IF(strain.LE.ZERO) THEN !compressive strain
              strain1=ZERO
              strain2=ZERO
              IF(DOP.AND.strain.LT.ZERO) THEN
                WRITE(OP_STRING,'('' >>Compressive strain in ENERGY'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ELSE IF(strain.GT.(TOL*pole1).OR.strain.GT.(TOL*pole2)) THEN
              IF(strain1.GT.(TOL*pole1)) strain1=TOL*pole1 !yielded
              IF(strain2.GT.(TOL*pole2)) strain2=TOL*pole2 !yielded
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Pole strain exceeded in ENERGY'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            DENOM1=pole1-strain1
            DtoAlpha1=1.0d0
            IF(alpha1.GT.ZERO) DtoAlpha1=DENOM1**alpha1
            DENOM2=pole2-strain2
            DtoAlpha2=1.0d0
            IF(alpha2.GT.ZERO) DtoAlpha2=DENOM2**alpha2
C           If outside pole-zero range - calc slopes for linear relns
            IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole1)) THEN
C             Calc slope for first set of params
              DWslope1=(2.0d0*coeff1*DENOM1*(DENOM1
     '          +2.0d0*alpha1*strain1)
     '          +coeff1*alpha1*(alpha1+1.0d0)*strain1*strain1)
     '          /(DtoAlpha1*DENOM1*DENOM1)
            ENDIF
            IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole2)) THEN
C             Calc slope for second set of params
              DWslope2=(2.0d0*coeff2*DENOM2*(DENOM2
     '          +2.0d0*alpha2*strain2)
     '          +coeff2*alpha2*(alpha2+1.0d0)*strain2*strain2)
     '          /(DtoAlpha2*DENOM2*DENOM2)
            ENDIF
            IF(strain.GT.ZERO) THEN !Tensile strain
              DW(iterm)=0.5d0*
     '          (coeff1*strain1*(2.0d0+alpha1*strain1/DENOM1)/DtoAlpha1
     '          +coeff2*strain2*(2.0d0+alpha2*strain2/DENOM2)/DtoAlpha2)
              IF(strain.GT.(TOL*pole1))
     '          DW(iterm)=DW(iterm)+(strain-strain1)*DWslope1/2.0d0
              IF(strain.GT.(TOL*pole2))
     '          DW(iterm)=DW(iterm)+(strain-strain2)*DWslope2/2.0d0
            ELSE !Compressive strain
C             Use deriv of DW wrt strain at strain=0 for compress. slope
              DW(iterm)=strain*(DWslope1+DWslope2)/2.0d0
            ENDIF
          ENDDO !iterm
        ENDIF

      ELSE IF(KTYP56(nr).EQ.4) THEN   !W is special function
        IF(KTYP55(nr).EQ.1) THEN      !W is a func of I1,I2,I3
          IF(KTYP52(nr).EQ.1) THEN    !W is
            DW(1)=0.0d0
            DW(2)=0.0d0
            DW(3)=0.0d0
          ELSE                    !W is
            DW(1)=0.0d0
            DW(2)=0.0d0
          ENDIF
        ELSE IF(KTYP55(nr).EQ.2) THEN !func of princ extn ratios (Ogden)
          IF(KTYP52(nr).EQ.1) THEN
            DW(1)=0.0d0
            DW(2)=0.0d0
            DW(3)=0.0d0
          ELSE
            DW(1)=0.0d0
            DW(2)=0.0d0
          ENDIF
        ELSE IF(KTYP55(nr).EQ.3) THEN !W=func of fibre and trans strains
          DW(1)=0.0d0
          DW(2)=0.0d0
          DW(3)=0.0d0
          DW(4)=0.0d0
          DW(5)=0.0d0
          DW(6)=0.0d0
        ENDIF

      ELSE IF(KTYP56(nr).EQ.5) THEN   !W is user defined function
        IF(KTYP55(nr).EQ.1) THEN      !in princ strain invariants
          CALL USER51(CG,DW,RK1,ERROR,*9999)
        ELSE IF(KTYP55(nr).EQ.2) THEN !in princ extension ratios
          CALL USER52(ERROR,*9999)
        ELSE IF(KTYP55(nr).EQ.3) THEN !in fibre & transverse strains
          CALL USER53(CG,DW,P1,P2,P3,P4,P5,P6,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('ENERGY')
      RETURN
 9999 CALL ERRORS('ENERGY',ERROR)
      CALL EXITS('ENERGY')
      RETURN 1
      END


      SUBROUTINE PFRE_NE(IBT,IDO,INP,NAN,NBH,NBJ,NFF,NGAP,NHE,
     '  NPF,NPNE,nr,NRE,NSP,NW,nx,NXI,
     '  CE,CP,FEXT,PF,PG,RE,XE,XW,ZE,ZW,ERROR,*)

C#### Subroutine: PFRE_NE
C###  Description:
C###    PFRE_NE evaluates contribution to element residuals RE(na,4)
C###    (associated with the hydrostatic pressure variable) for
C###    incompressible materials, arising from the stress constraint
C###    due to pressure bcs.

C**** Note: XW,ZW,TW,CW etc are element geometric, dependent, stress &
C****       material parameters etc interpolated at the centre of the
C****       pressure-loaded wall.
C**** 30APR89: NSP(i), i=1..2 are the hydrostatic press params coupled
C**** with the pressure contraint on the inside and the outside of the
C**** element, respectively.
C**** MPN 11-Jan-95: This routine is used when the hydrostatic pressure
C****                is interpolated using element based variables.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ipma50.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NFF(6),NGAP(NIM,NBM),NHE,
     '  NPF(15,NFM),NPNE(NNM,NBFM),nr,NRE(NEM),NSP(-2:2),NW,
     '  nx,NXI(-NIM:NIM)
      REAL*8 CE(NMM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),
     '  PF(2),PG(NSM,NUM,NGM,NBM),RE(NSM),
     '  XE(NSM,NJM),XW(NJM,NUM),ZE(NSM,NHM),ZW(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER iface,IFE,k,mi,mj,mz,nb,NBFF,NCW,nf,
     '  ngi1,ngi2,ngi3,ng_near,ni,NITB,nj,ns,ns2,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CW(NCW),DETERM,DET_F_NU,
     '  DNUDZ(3,3),DNUREFDZ(3,3),DXIXN(3,3),DXIZN(3,3),
     '  DXNZN(3,3),DZDNU(3,3),DZDNUREF(3,3),DZNXI(3,3),DZNXN(3,3),
     '  EG(3,3),GXL(3,3),GXU(3,3),RI1,RI2,RI3,RSIGN,RWX,SUM,
     '  TC(3,3),TCNUREF(3,3),TCNUREF33,TCRC(3,3),TW(3,3),TWA,XI(3)
      CHARACTER CHAR1*1

      CALL ENTERS('PFRE_NE',*9999)
      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMMX)',ERROR,*9999)
      CALL ASSERT(KTYP53(nr).GT.1,'stresses must be referred to '
     '  //'Nu coords',ERROR,*9999)

      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)
      DO mi=1,NITB
        XI(mi)=0.50d0
      ENDDO !mi
      IF(DOP) THEN
        WRITE(OP_STRING,'('' >>>PFRE_NE  diagnostic op'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO iface=1,2
        IF(DOP) THEN
          WRITE(OP_STRING,'('' >>> PF('',I1,'')= '',D12.3)')
     '      iface,PF(iface)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IFE=iface+4
        nf=NFF(IFE)
        NBFF=NPF(10,nf) !basis fn for first geometric variable
        XI(3)=DBLE(iface-1)

        CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,ERROR,*9999)
C old
C        DO il=1,ILT(1,nr,nx)
C          IF(ILP(il,1,nr,nx).NE.1.OR.ILP(il,1,nr,nx).NE.2) THEN
CC           constant spatially or defined by elements
C            CW(il)=CE(il)
C          ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
C            CW(il)=0.0d0
C            DO nnbf=1,NNT(NBFF)
C              CW(il)=CW(il)+CP(il,NPNE(nnbf,nb))
C            ENDDO
C            CW(il)=CW(il)/DBLE(NNT(NBFF))
C          ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by Gauss points
C            CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
C     '        //'not implemented',ERROR,*9999)
C          ENDIF
C        ENDDO !il

C       Interpolate midwall geometric vars XW and derivs wrt Xi
        CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XW,XI,ERROR,*9999)
C       Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C       derivs of Xi wrt undef Nu coords, DXIXN, (JP=1).
        CALL XGMG(1,NITB,nb,nr,DXIXN,GXL,GXU,RWX,XW,ERROR,*9999)
C       Interpolate dependent var.s ZW and derivs wrt Nu (JP=1)
        CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '    DXIXN,ZE,ZW,XI,ERROR,*9999)
C       Calculate deformed metric tensors wrt Nu (AZL,AZU)
        CALL ZGMG(nb,AZ,AZL,AZU,ZW,ERROR,*9999)
C       Get contravariant cpts of 2nd Piola-Kirchhoff stress
C       tensor (TW) wrt undeformed Nu coordinates
        IF(KTYP51(nr).EQ.1) THEN
          CALL ZGTG51(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TW,ZW,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.2) THEN
          CALL ZGTG52(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TW,ZW,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.3) THEN
          CALL ZGTG53(nb,nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TW,XW,ZW,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.4) THEN
          CALL ZGTG54(nb,nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TW,ERROR,*9999)
        ENDIF
C new MPN 4-May-96: new way of handling sheets
C       Get derivs of Xi wrt deformed Nu coords, DXIZN, and inv., DZNXI
        CALL DXIDZN(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '    DXIZN,DZNXI,PG,XE,XW,XI,ZE,ZW,ERROR,*9999)
C       Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
C       and inverse (DXNZN)
        DO ni=1,NITB
          DO mi=1,NITB
            SUM=0.0d0
            DO k=1,NITB
              SUM=SUM+DZNXI(ni,k)*DXIXN(k,mi)
            ENDDO
            DZNXN(ni,mi)=SUM
          ENDDO
        ENDDO
        CALL INVERT(NITB,DZNXN,DXNZN,DET_F_NU)
        IF(DOP) THEN
          WRITE(CHAR1,'(I1)') NITB
          WRITE(OP_STRING,'(''  DZNXN:'','//CHAR1(1:1)//'D12.4,'
     '      //'''  DXNZN:'','//CHAR1(1:1)//'D12.4,'
     '      //'/(8X,'//CHAR1(1:1)//'D12.4,8X,'//CHAR1(1:1)//'D12.4))')
     '      ((DZNXN(mi,ni),ni=1,NITB),(DXNZN(mi,ni),ni=1,NITB),
     '      mi=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C old
CC       Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
C        CALL DXIDNU(NBJ(1),nr,DXIXN,DXNXI,GXL,GXU,XW,ERROR,*9999)
CC       Interpolate dependent var.s ZG and derivs wrt Xi
C        CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
C     '    DXIX,ZE,ZW,XI,ERROR,*9999)
C        CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZW,ERROR,*9999)
C        CALL DXIDNU(NBJ(1),nr,DXIZN,DZNXI,GZL,GZU,XW,ERROR,*9999)
CC       Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
CC       and inverse (DXNZN)
C        DO ni=1,NITB
C          DO mi=1,NITB
C            SUM1=0.0d0
C            SUM2=0.0d0
C            DO k=1,NITB
C              SUM1=SUM1+DZNXI(ni,k)*DXIXN(k,mi)
C              SUM2=SUM2+DXNXI(ni,k)*DXIZN(k,mi)
C            ENDDO
C            DZNXN(ni,mi)=SUM1
C            DXNZN(ni,mi)=SUM2
C          ENDDO
C        ENDDO
C end old
        IF(KTYP53(nr).EQ.3) THEN !Active stress component included
C!!!      Pick Gauss pt nearest to centre of current face in
C!!!      current element. To be strictly correct need to either
C!!!      store FEXT at face basis Gauss pts or somehow interpolate
C!!!      element Gauss pt values to the central point on the face
          ngi1=NGAP(1,nb)
          ngi2=NGAP(2,nb)
          ngi3=NGAP(3,nb)
          ng_near=(ngi1+1)/2 + ((ngi2-1)/2)*ngi1
     '      + (iface-1)*(ngi3-1)*ngi2*ngi1 !need integer division
          CALL ZGTG5A(nr,FEXT(1,ng_near),DXNZN,DZNXN,
     '      CW(IL_time_delay),TW,TWA,ERROR,*9999)
        ENDIF
C new MPN 5-May-96: rotate TW from def material coordinates to
C                   deformed fibre reference (wall) coordinates
C       Get Physical Cauchy stress tensor TC wrt deformed nu-material
C       coordinates from TW
        DO mz=1,3
          DO nz=1,3
            SUM=0.0d0
            DO mj=1,NITB
              DO nj=1,NITB
                SUM=SUM+DZNXN(mz,mj)*TW(mj,nj)*DZNXN(nz,nj)
              ENDDO !nj
            ENDDO !mj
            TC(mz,nz)=SUM/DET_F_NU
          ENDDO !nz
        ENDDO !mz
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TC:'',12X,3D12.4,/(16X,3D12.4))')
     '      ((TC(mz,nz),nz=1,3),mz=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(JTYP9.GE.2) THEN !imbric (+ sheet) angles defined
C         Compute def anatomical fibre vects wrt rc coords at XI
          CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '      DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '      PG,XE,XW,XI,ZE,ZW,ERROR,*9999)
          CALL INVERT(NITB,DZDNU,DNUDZ,DETERM)
C         Compute cmpts of Cauchy stress tensor wrt rc coords
C         by rotating def material coord system into rc coords
          DO mi=1,3
            DO ni=1,3
              SUM=0.0d0
              DO mj=1,3
                DO nj=1,3
                  SUM=SUM+DZDNU(mi,mj)*TC(mj,nj)*DNUDZ(nj,ni)
                ENDDO !nj
              ENDDO !mj
              TCRC(mi,ni)=SUM
            ENDDO !ni
          ENDDO !mi
          IF(DOP) THEN
            WRITE(OP_STRING,'('' TCRC:'',12X,3D12.4,'
     '        //'/(18X,3D12.4))') ((TCRC(mi,ni),ni=1,3),mi=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C         Compute deformed fibre ref vectors wrt rc coords at XI
          CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
     '      NHE,NITB,nr,nx,DZDNUREF(1,1),DZDNUREF(1,2),DZDNUREF(1,3),
     '      PG,XI,ZE,ZW,ERROR,*9999)
          CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
C         Compute cmpts of Cauchy stress tensor wrt deformed
C         fibre reference coords by rotating rc coord system
C         into deformed fibre ref coords
C         NOTE: only use TCNUREF(3,3); could delete outer loops
          DO mi=1,3
            DO ni=1,3
              SUM=0.0d0
              DO mj=1,3
                DO nj=1,3
                  SUM=SUM+DNUREFDZ(mi,mj)*TCRC(mj,nj)*DZDNUREF(nj,ni)
                ENDDO !nj
              ENDDO !mj
              TCNUREF(mi,ni)=SUM
            ENDDO !ni
          ENDDO !mi
          IF(DOP) THEN
            WRITE(OP_STRING,'('' TCNUREF:'',12X,3D12.4,'
     '        //'/(21X,3D12.4))')
     '        ((TCNUREF(mi,ni),ni=1,3),mi=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE !at most fibres defined
C         TCNUREF33 is aligned with TC(3,3) for
C         zero imbric/sheet angles
          TCNUREF(3,3)=TC(3,3)
        ENDIF !JTYP9.GE.2
        TCNUREF33=TCNUREF(3,3)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TCNUREF33='',D12.4)') TCNUREF33
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C old
CC       Get Physical Cauchy stress TC33 wrt deformed Nu_3 from TW
C        TC33=0.0d0
C        DO nj=1,NITB
C          DO mj=1,NITB
C            TC33=TC33+DZNXN(3,nj)*DZNXN(3,mj)*TW(nj,mj)
C          ENDDO
C        ENDDO
C        TC33=TC33/DSQRT(RI3)

        IF(iface.EQ.1) THEN !Xi3=0 face
          RSIGN= 1.0d0
        ELSE IF(iface.EQ.2) THEN !Xi3=1 face
          RSIGN=-1.0d0
        ENDIF

        ns=NSP(iface)
        IF(iface.EQ.1.AND.(NXI(-3).EQ.0.OR.nr.NE.NRE(NXI(-3)))
     '    .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
     '     iface.EQ.2.AND.(NXI( 3).EQ.0.OR.nr.NE.NRE(NXI( 3)))
     '    .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext press bc on Xi3=1 face
          RE(ns)=TCNUREF33+PF(iface) !match norm. stress for press bcs
        ELSE !no external pressure bc applied on current face
          IF(KTYP5A(nr).EQ.1) THEN !match hyd. press across elems
            RE(ns)=ZW(NH_LOC(NH_LOC(0,nx),nx),1)*RSIGN
          ELSE IF(KTYP5A(nr).EQ.2) THEN !match norm. stress across elems
            RE(ns)=TCNUREF33*RSIGN
          ENDIF
        ENDIF

        ns2=NSP(-iface) !ns for pressure bc params
        IF(KTYP5A(nr).EQ.1) THEN !match hyd. press across elems
          RE(ns2)=ZW(NH_LOC(NH_LOC(0,nx),nx),1)*RSIGN
        ELSE IF(KTYP5A(nr).EQ.2) THEN !match norm. stress across elems
          RE(ns2)=TCNUREF33*RSIGN
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,'('' RE('',I2,'')= '',D12.3)') ns,RE(ns)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RE('',I2,'')= '',D12.3)') ns2,RE(ns2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !iface

      CALL EXITS('PFRE_NE')
      RETURN
 9999 CALL ERRORS('PFRE_NE',ERROR)
      CALL EXITS('PFRE_NE')
      RETURN 1
      END


      SUBROUTINE PFRE_NP(IBT,IDO,INP,NAN,NBH,NBJ,NFF,NGAP,NHE,NJE,NKF,
     '  NNF,NPF,NPNE,nr,NRE,NW,nx,NXI,
     '  CE,CG,CP,FEXT,PF,PG,RE,VE,WG,XE,XG,ZE,ZG,ERROR,*)

C#### Subroutine: PFRE_NP
C###  Description:
C###    PFRE_NP evaluates contribution to element residuals RE(ns,4)
C###    (associated with the hydrostatic pressure variable) for
C###    incompressible materials, arising from the stress constraint
C###    due to pressure bcs.

C**** MPN 11-Jan-95: This routine is used when the hydrostatic pressure
C****                is interpolated using nodal based variables.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NFF(6),NGAP(NIM,NBM),NHE,
     '  NJE,NKF(0:4,16,6,NBFM),NNF(0:17,6,NBFM),NPF(15,NFM),
     '  NPNE(NNM,NBFM),nr,NRE(NEM),NW,nx,NXI(-NIM:NIM)
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),PF(2),
     '  PG(NSM,NUM,NGM,NBM),RE(NSM),VE(NSM,NKM),WG(NGM,NBM),
     '  XE(NSM,NJM),XG(NJM,NUM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER iface,IFE,mi,mj,NBFF,NBP,nf,ng,ng1,ng2,NGI1,NGI2,
     '  ni,NITB,nj,nk,nkbf,nn,nnbf,nse,nsf
      REAL*8 D(5,5),DETERM,DNUDZ(3,3),DNUREFDZ(3,3),
     '  DZDNU(3,3),DZDNUREF(3,3),PHI(3),PST(3),
     '  RGX2D,RGX_local,RGZ,RGZ2D,RM(3,3),RSIGN,RWG,
     '  SUM,TC(3,3),TCNUREF(3,3),TCNUREF33,TCRC(3,3),
     '  TG(3,3),TN(3,3),XI(3)
      LOGICAL Papplied
      DATA D/5*0.0d0,-0.288675134594813d0,0.288675134594813d0,3*0.0d0,
     '       -0.387298334620741d0,0.0d0,0.387298334620741d0,2*0.0d0,
     '       -0.430568155797026d0,    -0.169990521792428d0,
     '        0.169990521792428d0,     0.430568155797026d0,  0.0d0,
     '       -0.453089922969332d0,    -0.269234655052841d0,  0.0d0,
     '        0.269234655052841d0,     0.453089922969332d0/

      CALL ENTERS('PFRE_NP',*9999)
      NITB=NIT(NBH(NH_LOC(1,nx)))

      IF(DOP) THEN
        WRITE(OP_STRING,'('' >>>PFRE_NP diagnostic op:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis function for pressure vars
      DO iface=1,2
        IF(iface.EQ.1) THEN !Xi3=0 face
          RSIGN= 1.0d0
        ELSE IF(iface.EQ.2) THEN !Xi3=1 face
          RSIGN=-1.0d0
        ENDIF
        IF(iface.EQ.1.AND.(NXI(-3).EQ.0.OR.nr.NE.NRE(NXI(-3))).AND.
     '    (NW.EQ.2.OR.NW.EQ.4).OR. !ext. press bc applied on Xi3=0 face
     '    iface.EQ.2.AND.(NXI(3).EQ.0.OR.nr.NE.NRE(NXI(3))).AND.
     '    (NW.EQ.3.OR.NW.EQ.4)) THEN !ext press bc applied on Xi3=1 face
          Papplied=.TRUE.
        ELSE
          Papplied=.FALSE.
        ENDIF
        IFE=iface+4
        nf=NFF(IFE)
        NBFF=NPF(12,nf) !basis fn for third geometric variable (theta)
C                       !!!WARNING: Need to store pressure face bases
C                       !!!in NPF
        XI(3)=DBLE(iface-1)
        IF(DOP) THEN
          IF(Papplied) THEN
            WRITE(OP_STRING,'(/'' Xi(3)='',I1,'' face (nf='',I3,'').'
     '        //' Applied pressure ='',D12.4)') iface-1,nf,PF(iface)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(/'' Xi(3)='',I1,'' face (nf='',I3,'').'
     '        //' No applied pressure bc'')') iface-1,nf
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' Pressure face basis nbff='',I2)') nbff
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C       Loop over face nodes and derivs for pressure basis
        nsf=0
C       NOTE: nnbf must loop over nn's for face basis fn of first
C       geometric variable since pressure basis faces aren't stored
C       in NPF. Shouldn't cause problems.
        DO nnbf=1,NNT(NBFF)
          nn=NNF(1+nnbf,IFE,NBP)
          DO nkbf=1,NKT(nnbf,NBP)
            nk=NKF(nkbf,nnbf,IFE,NBP)
            nsf=nsf+1
            nse=nk+(nn-1)*NKT(nn,NBP)
            RE(nse)=0.0d0
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' nn='',I2,'' nk='',I2,'
     '          //''' nsf='',I2,'' nse='',I2)') nn,nk,nsf,nse
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
C           Loop over gauss pts on face
            ng=0
            NGI1=NGAP(1,NBP) !Gauss pts in Xi1 for pressure basis
            NGI2=NGAP(2,NBP) !Gauss pts in Xi2 for pressure basis
            DO ng2=1,NGI2
              DO ng1=1,NGI1
                ng=ng+1
                XI(1)=0.5d0+D(ng1,NGI1)
                XI(2)=0.5d0+D(ng2,NGI2)
C               Calculate stress at XI (which is
C               position of current Gauss pt on face)
                CALL ZETX50('Fibre','Cauchy',IBT,IDO,INP,
     '            NAN,NBH,NBJ,0,NHE,NJE,NPNE,nr,nx,
     '            CE,CG,CP,FEXT(1,ng),PG,PHI,PST,
     '            RGX_local,RGX2D,RGZ,RGZ2D,RM,TC,TG,TN,VE,
     '            XE,XG,XI,ZE,ZG,ERROR,*9999)
C               Deformed Jacobian (integrating wrt deformed area)
C               Using the 2D Jacobian for the face integration
                RWG=RGZ2D*WG(ng,NBP)
                IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*ZG(1,1) !cyl symm in x
                IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*ZG(2,1) !cyl symm in y
                IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*ZG(1,1)**2 !sph sym
                IF(JTYP9.GE.2) THEN !imbric (+ sheet) angles defined
C                 Compute def anatomical fibre vects wrt rc coords at XI
                  CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '              DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '              PG,XE,XG,XI,ZE,ZG,ERROR,*9999)
                  CALL INVERT(NITB,DZDNU,DNUDZ,DETERM)
C                 Compute cmpts of Cauchy stress tensor wrt rc coords
C                 by rotating def material coord system into rc coords
                  DO mi=1,3
                    DO ni=1,3
                      SUM=0.0d0
                      DO mj=1,3
                        DO nj=1,3
                          SUM=SUM+DZDNU(mi,mj)*TC(mj,nj)*DNUDZ(nj,ni)
                        ENDDO !nj
                      ENDDO !mj
                      TCRC(mi,ni)=SUM
                    ENDDO !ni
                  ENDDO !mi
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' TCRC:'',12X,3D12.4,'
     '                //'/(18X,3D12.4))') ((TCRC(mi,ni),ni=1,3),mi=1,3)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
C                 Compute deformed fibre ref vectors wrt rc coords at XI
                  CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,0,
     '              NHE,NITB,nr,nx,
     '              DZDNUREF(1,1),DZDNUREF(1,2),DZDNUREF(1,3),
     '              PG,XI,ZE,ZG,ERROR,*9999)
                  CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
C                 Compute cmpts of Cauchy stress tensor wrt deformed
C                 fibre reference coords by rotating rc coord system
C                 into deformed fibre ref coords
C                 NOTE: only use TCNUREF(3,3);could delete outer loops
                  DO mi=1,3
                    DO ni=1,3
                      SUM=0.0d0
                      DO mj=1,3
                        DO nj=1,3
                          SUM=SUM+DNUREFDZ(mi,mj)*TCRC(mj,nj)*
     '                      DZDNUREF(nj,ni)
                        ENDDO !nj
                      ENDDO !mj
                      TCNUREF(mi,ni)=SUM
                    ENDDO !ni
                  ENDDO !mi
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' TCNUREF:'',12X,3D12.4,'
     '                //'/(21X,3D12.4))')
     '                ((TCNUREF(mi,ni),ni=1,3),mi=1,3)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE !at most fibres defined
C                 TCNUREF33 is aligned with TC(3,3) for
C                 zero imbric/sheet angles
                  TCNUREF(3,3)=TC(3,3)
                ENDIF !JTYP9.GE.2
                TCNUREF33=TCNUREF(3,3)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' TCNUREF33='',D12.4)') TCNUREF33
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF

                IF(Papplied) THEN !pressure bc applied on current face
                  RE(nse)=RE(nse)+
     '              (TCNUREF33+PF(iface))*PG(nsf,1,ng,NBP)*RWG
                ELSE !no pressure bc applied on current face
                  RE(nse)=RE(nse)+RSIGN*TCNUREF33*PG(nsf,1,ng,NBP)*RWG
                ENDIF

                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ng='',I2,'' rgz2d='',D12.4,'
     '              //''' wg='',D12.4,'' PG(nsf..)='',D12.4,'
     '              //''' TCNUREF33='',D12.4)')
     '              ng,rgz2d,WG(ng,NBP),PG(nsf,1,ng,NBP),TCNUREF33
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !ng1
            ENDDO !ng2
            IF(DOP) THEN
              WRITE(OP_STRING,'('' RE('',I2,'')='',D12.4)')
     '          nse,RE(nse)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nkbf (nk)
        ENDDO !nnbf (nn)
      ENDDO !iface

      CALL EXITS('PFRE_NP')
      RETURN
 9999 CALL ERRORS('PFRE_NP',ERROR)
      CALL EXITS('PFRE_NP')
      RETURN 1
      END


      SUBROUTINE PFRF(IBT,IDO,INP,IXF,NAN,NBH,NBJ,NGAP,NJE,NPF,nr,nx,
     '  PF,PG,RF,WG,XE,XG,ZE,ZG,ERROR,*)

C#### Subroutine: PFRF
C###  Description:
C###    PFRF evaluates contribution RF(ns,nj) to element residuals RE
C###    from the pressure PF(i) acting on the Xi(3)=IXF face
C###    (IXF=iface-1).

C**** Note: GZL & GZU are deformed state metric tensors wrt Xi.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  IXF,NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NGAP(NIM,NBM),NJE,
     '  NPF(15),nx
      REAL*8 PF,PG(NSM,NUM,NGM,NBM),RF(32,3),WG(NGM,NBM),XE(NSM,NJM),
     '  XG(NJM,NUM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,NBFF,ng,ng1,ng2,NGI1,NGI2,NITB,nh,nhx,ni,
     '  nj,njj,nr,ns,NU1(0:3)
      REAL*8 D(5,5),DXIX(3,3),GZ,GZL(3,3),GZU(3,3),RGZ,RWG,SUM,XI(3)
      DATA NU1/1,2,4,7/
      DATA D/5*0.0d0,-0.288675134594813d0,0.288675134594813d0,3*0.0d0,
     '       -0.387298334620741d0,0.0d0,0.387298334620741d0,2*0.0d0,
     '       -0.430568155797026d0,    -0.169990521792428d0,
     '        0.169990521792428d0,     0.430568155797026d0,  0.0d0,
     '       -0.453089922969332d0,    -0.269234655052841d0,  0.0d0,
     '        0.269234655052841d0,     0.453089922969332d0/

      CALL ENTERS('PFRF',*9999)
      NITB=NIT(NBH(NH_LOC(1,nx)))

      DO njj=1,NJ_LOC(NJL_GEOM,0)
        nj=NJ_LOC(NJL_GEOM,njj)
        NBFF=NPF(9+nj)
        DO ns=1,NST(NBFF)+NAT(NBFF)
          RF(ns,nj)=0.0d0
        ENDDO !ns
      ENDDO !nj

      XI(3)=DBLE(IXF)
      NGI1=NGAP(1,NPF(10))
      NGI2=NGAP(2,NPF(10))
      ng=0
      DO ng2=1,NGI2
        DO ng1=1,NGI1
          ng=ng+1
          XI(1)=0.5d0+D(ng1,NGI1)
          XI(2)=0.5d0+D(ng2,NGI2)
          CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
          CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NJE,nr,nx,DXIX,ZE,ZG,XI,ERROR,
     '      *9999)
          CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' >>>PFRF diagnostic op at Gauss pt '',I2)') NG
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nhx=1,NJ_LOC(NJL_GEOM,0)
              nj=NJ_LOC(NJL_GEOM,nhx)
              WRITE(OP_STRING,'(''  XG('',I1,'',ni): '',4D12.4)')
     '          nj,(XG(nj,NU1(ni)),ni=0,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''  ZG('',I1,'',ni): '',4D12.4)')
     '          nhx,(ZG(nhx,NU1(ni)),ni=0,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nj
            DO mi=1,NITB
              WRITE(OP_STRING,'('' GZU('',I1,'',ni): '',3D12.4)')
     '          MI,(GZU(mi,ni),ni=1,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !mi
          ENDIF
          RGZ=DSQRT(GZ)
          RWG=RGZ*WG(ng,NPF(10))

C old MPN 28-Jun-1995: integrals now wrt deformed coords
C          IF(ITYP10(nr).EQ.1) THEN
C            DO nhx=1,NJ_LOC(NJL_GEOM,0)
C              nx=NH_LOC(nhx,nx)
C              DO ni=1,NITB
C                DZI(nhx,ni)=ZG(nhx,NU1(ni))
C              ENDDO !ni
C            ENDDO !nhx
C          ELSE IF(ITYP10(nr).EQ.2) THEN
C            RX=XG(1,1)
C            RZ=ZG(1,1)
C            DT=ZG(2,1)-XG(2,1)
C            CT=DCOS(DT)
C            ST=DSIN(DT)
C            DO ni=1,NITB
C              DZI(1,ni)= ZG(1,NU1(ni))*CT-RZ*ST*ZG(2,NU1(ni))
C              DZI(2,ni)=(ZG(1,NU1(ni))*ST+RZ*CT*ZG(2,NU1(ni)))/RX
C              DZI(3,ni)= ZG(3,NU1(ni))
C            ENDDO !ni
C          ELSE IF(ITYP10(nr).EQ.3) THEN
C            RX=XG(1,1)
C            RZ=ZG(1,1)
C            DT=ZG(2,1)-XG(2,1)
C            CT=DCOS(DT)
C            ST=DSIN(DT)
C            CX=DCOS(XG(3,1))
C            SX=DSIN(XG(3,1))
C            CZ=DCOS(ZG(3,1))
C            SZ=DSIN(ZG(3,1))
C            CC=CX*CZ
C            SS=SX*SZ
C            CS=CX*SZ
C            SC=SX*CZ
C            DO ni=1,NITB
C              RB=ZG(1,NU1(ni))
C              TB=ZG(2,NU1(ni))
C              PB=ZG(3,NU1(ni))
C              DZI(1,ni)=CC*(RB*CT-RZ*ST*TB)-CS*RZ*CT*PB+SC*RZ*PB+SS*RB
C              DZI(2,ni)=(RZ*CZ*CT*TB+(RB*CZ-RZ*SZ*PB)*ST)/(RX*CX)
C              DZI(3,ni)=(CC*RZ*PB+CS*RB+SC*(RZ*ST*TB-RB*CT)+SS*RZ*CT*PB)
C     '                                                   /RX
C            ENDDO !ni
C          ELSE IF(ITYP10(nr).EQ.4) THEN
C            SLX=DSINH(XG(1,1))
C            SLZ=DSINH(ZG(1,1))
C            SMX=DSIN (XG(2,1))
C            SMZ=DSIN (ZG(2,1))
C            CLX=DSQRT(1.0d0+SLX*SLX)
C            CLZ=DSQRT(1.0d0+SLZ*SLZ)
C            CMX=DSQRT(1.0d0-SMX*SMX)
C            CMZ=DSQRT(1.0d0-SMZ*SMZ)
C            CSLX=CLX/SLX
C            CSMX=CMX/SMX
C            DT=ZG(3,1)-XG(3,1)
C            CT=DCOS(DT)
C            ST=DSIN(DT)
C            CCL=CLX*CLZ
C            CSL=CLX*SLZ
C            SCL=SLX*CLZ
C            SSL=SLX*SLZ
C            CC=CMX*CMZ
C            CS=CMX*SMZ
C            SC=SMX*CMZ
C            SS=SMX*SMZ
C            G1=SLX*SLX+SMX*SMX
C            G3=SLX*SLX*SMX*SMX
C            DO ni=1,NITB
C              DLB=ZG(1,NU1(ni))
C              DMB=ZG(2,NU1(ni))
C              DTB=ZG(3,NU1(ni))
C              DZI(1,ni)=(( SSL*CC+CCL*SS*CT)*DLB+(-SCL*CS+CSL*SC*CT)*DMB
C     '                                                -CSL*SS*ST*DTB)/G1
C              DZI(2,ni)=((-CSL*SC+SCL*CS*CT)*DLB+( CCL*SS+SSL*CC*CT)*DMB
C     '                                                -SSL*CS*ST*DTB)/G1
C              DZI(3,ni)=(SCL*SS*ST*DLB+SSL*SC*ST*DMB+SSL*SS*CT*DTB)/G3
C            ENDDO !ni
C          ENDIF
          DO nhx=1,NJ_LOC(NJL_GEOM,0)
            nh=NH_LOC(nhx,nx)
            NBFF=NPF(9+nhx)
            SUM=0.0d0
            DO ni=1,NITB
C MPN 28-Jun-1995: integrals now wrt deformed coords
              SUM=SUM+GZU(ni,3)*ZG(nhx,NU1(ni))
C old              SUM=SUM+GZU(ni,3)*DZI(nhx,ni)
            ENDDO !ni
            DO ns=1,NST(NBFF)+NAT(NBFF)
              IF(IXF.EQ.1)THEN           !positive face
                RF(ns,nh)=RF(ns,nh)-PF*SUM*PG(ns,1,ng,NBFF)*RWG
              ELSE                       !negative face
                RF(ns,nh)=RF(ns,nh)+PF*SUM*PG(ns,1,ng,NBFF)*RWG
              ENDIF
            ENDDO !ns
          ENDDO !nhx
        ENDDO !ng1
      ENDDO !ng2

      IF(DOP) THEN
        DO nhx=1,NJ_LOC(NJL_GEOM,0)
          nh=NH_LOC(nhx,nx)
          NBFF=NPF(9+nhx)
          WRITE(OP_STRING,
     '      '('' RF(ns,'',I2,''): '',5D12.4,/(12X,5D12.4))')
     '      nh,(RF(ns,nh),ns=1,NST(NBFF)+NAT(NBFF))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nhx
      ENDIF

      CALL EXITS('PFRF')
      RETURN
 9999 CALL ERRORS('PFRF',ERROR)
      CALL EXITS('PFRF')
      RETURN 1
      END


      SUBROUTINE USER51(CG,DW,RK1,ERROR,*)

C#### Subroutine: USER51
C###  Description:
C###    USER51 returns the derivatives of the user-defined strain energy
C###    function of principal strain invariants I1,I2,I3,K1,K2 to
C###    subroutine ENERGY at current Gauss point.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      REAL*8 CG(NMM),DW(6),RK1
      CHARACTER ERROR*(*)

      CALL ENTERS('USER51',*9999)
C     DW(1)=CG(1)*DEXP(CG(2)*(RI1-3.0d0))
C     DW(2)=0.0d0
C     DW(3)=0.0d0
C     DW(4)=CG(3)*(DEXP(CG(4)*RK1)-1.0d0)
C     DW(5)=CG(5)*(DEXP(CG(6)*RK2))
      DW(1)=CG(1)
      DW(2)=CG(2)
      DW(3)=0.0d0
      DW(4)=2.0d0*CG(3)*RK1
      DW(5)=CG(4)
!news 28-MAY-1991 Function proposed by Humphrey/Yin.  JSW.
!     RK12=DSQRT(2.0d0*RK1+1.0d0)
!     DW(1)=CG(3)+CG(4)*(RK12-1.0d0)+2.0d0*CG(5)*(RI1-3.0d0)
!     DW(2)=0.0d0
!     DW(3)=0.0d0
!     DW(4)=2.0d0*CG(1)*(RK12-1.0d0)+3.0d0*CG(2)*(RK12-1.0d0)**2.0d0
!    '  +CG(4)*(RI1-3.0d0)
!     DW(4)=DW(4)/RK12
!     DW(5)=0.0d0
!newe

      CALL EXITS('USER51')
      RETURN
 9999 CALL ERRORS('USER51',ERROR)
      CALL EXITS('USER51')
      RETURN 1
      END


      SUBROUTINE USER52(ERROR,*)

C#### Subroutine: USER52
C###  Description:
C###    USER52 returns the derivatives of the user-defined strain energy
C###    function of principal extension ratios L1, L2, and L3 to
C###    subroutine ENERGY at current Gauss point.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)

      CALL ENTERS('USER52',*9999)

C     DW(1)=
C     DW(2)=
C     DW(3)=

      CALL EXITS('USER52')
      RETURN
 9999 CALL ERRORS('USER52',ERROR)
      CALL EXITS('USER52')
      RETURN 1
      END


      SUBROUTINE USER53(CG,DW,E11,E22,E33,E12,E13,E23,ERROR,*)

C#### Subroutine: USER53
C###  Description:
C###    USER53 returns the derivatives of the user-defined strain energy
C###    function of fibre and transverse physical strains E11...E23
C###    to subroutine ENERGY at current Gauss point.
C###    Double exponential law.
C###    Exponential-dependence on fibre/transverse extension ratios

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      REAL*8 CG(NMM),DW(6),E11,E12,E13,E22,E23,E33
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 CEXPQ1,CEXPQ2,TOL

      CALL ENTERS('USER53',*9999)

      TOL=1.0D-08
      IF((DABS(E11).LT.TOL).AND.(DABS(E22).LT.TOL).AND.
     '  (DABS(E33).LT.TOL).AND.(DABS(E12).LT.TOL).AND.
     '  (DABS(E13).LT.TOL).AND.(DABS(E23).LT.TOL)) THEN
        CEXPQ1=CG(1)
        CEXPQ2=CG(4)
      ELSE
        CEXPQ1=CG(1)*DEXP(CG(2)*E11*E11+CG(3)*(E12*E12+E13*E13))
        CEXPQ2=CG(4)*DEXP(CG(5)*(E22+E33)**2+CG(6)*(E22*E33-E23*E23))
      ENDIF
      DW(1)=CEXPQ1*2.0d0*CG(2)*E11
      DW(2)=CEXPQ2*(2.0d0*CG(5)*(E22+E33)+CG(6)*E33)
      DW(3)=CEXPQ2*(2.0d0*CG(5)*(E22+E33)+CG(6)*E22)
      DW(4)= CEXPQ1*CG(3)*2.0d0*E12
      DW(5)= CEXPQ1*CG(3)*2.0d0*E13
      DW(6)=-CEXPQ2*CG(6)*2.0d0*E23

      CALL EXITS('USER53')
      RETURN
 9999 CALL ERRORS('USER53',ERROR)
      CALL EXITS('USER53')
      RETURN 1
      END


      SUBROUTINE ZERE50(IBT,IDO,INP,NAN,NBH,NBJ,ne,NFF,NGAP,
     '  NHE,NJE,NKF,NNF,NPF,NPNE,nr,NRE,NW,nx,NXI,
     '  CE,CG,CP,FEXT,PG,RE,RGX,SE,VE,WG,XE,XG,ZE,ZG,ERROR,*)

C#### Subroutine: ZERE50
C###  Description:
C###    ZERE50 calculates element residual RE from current dependent
C###    variable array ZE.

C**** X variables refer to orthog curvilinear coords in reference state.
C**** Z     "       "         "        "         "      deformed    "
C**** Material Theta-coordinates (reference for deformation) coincide
C****   with Xj-coords in ref state.
C**** Material Xi-coords are the finite element mesh coordinates.
C**** Material Nu-coordinates (reference for stresses): are orthogonal
C****   and (Nu1,Nu2) lie in the (Xi1-Xi2) plane such that Nu(1)
C****   is aligned with the 'fibres' to which material aeolotropy
C****   is referred; The undeformed base vectors are defined such that
C****   the undeformed metric tensors wrt the Nu are delta(i,j).
C****
C**** ITYP2(nr,nx) is 1..15 for problem type
C**** ITYP4(nr,nx) is 1..4: fem/direct bem/indirect bem/orthog colloc.
C**** ITYP5(nr,nx) is 1..5: static/time integration/modal analysis
C****                       /Fourier analysis/buckling analysis
C**** ITYP6(nr,nx) is 1,2:  linear/nonlinear problem
C**** KTYP5  is 1..3: initial solution zero/read in/restarted
C**** KTYP7  is 1..3: eqn parameters constant wrt time/user defined
C****                 /read from file
C**** KTYP8  is 1..6: geom/fibre/field/potential/Fourier/opt.n fitting
C**** ITYP9(nr,nx) is 1..3: solution by full Newton/modified Newton
C****                 /BFGS inverse/Conjugate gradient
C**** KTYP10 is 1,2 : solution with no search/linear search
C**** KTYP12 is 1,0 : fitting with/without constraints
C**** KTYP13 is 1 if pressure read from file (PRESS.VSAERO)
C**** KTYP16 is 1,2 : lowest/highest eigenvalue
C**** KTYP17 is number of eigenvalue pairs
C**** KTYP18 is number of subspace iteration vectors
C**** KTYP19 is number of starting vectors
C**** KTYP22 is 1..3: time integ.n algorithm linear/quadratic/cubic
C**** KTYP23 is 1,2:  time step fixed/calculated to control error
C**** KTYP25 is 1..3: b.c. in form of impulse/step/sine wave
C**** KTYP26 is 1..2: opt of material params/geometric params
C**** KTYP27 is 1..5: type of minimization objective function
C**** KTYP31 is 1,2:  activation model forwards/backwards
C**** KTYP43 is 1..3: linear elastic material isotropic/trans.isotr.
C****                 /orthotropic
C**** KTYP51(nr) is 1..6: plane stress/plane strain/3D/membrane
C****                 /string/shell
C**** KTYP52(nr) is 1..3: compress/incomp/incomp with fluid perfusate
C**** KTYP53(nr) is 1..3: isotrop/aeleotropic/aeleo + active fibres
C**** KTYP54(nr) is 1..3: hyperelastic/Cauchy-elasticity/creep
C**** KTYP55(nr) is 1..3: strain invariants/ext ratios/fibre strains
C**** KTYP56(nr) is 1..3: polynomial/special function/user defined
C**** KTYP57(nr) is type of pressure loading applied to elements
C**** KTYP58(nr) is 1,2: conventional/isochoric element
C**** KTYP59(nr)  is elastance/Hill-type/fading-memory formulation

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:aero00.cmn'
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ipma50.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ne,NFF(6),NGAP(NIM,NBM),
     '  NHE,NJE,NKF(0:4,16,6,NBFM),NNF(0:17,6,NBFM),NPF(15,NFM),
     '  NPNE(NNM,NBFM,NEFM),nr,NRE(NEM),NW,nx,NXI(-NIM:NIM,0:NEM)
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),RGX(NGM),SE(NSM,NBFM,NEFM),
     '  VE(NSM,NKM),WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER iface,IFE,IHANDLE,IXF,j,JP,k,mi,
     '  na,nb,NBE,NBFF,NBP,nf,ng,nh,nh1,nhx,ni,NITB,nj,njj,njj1,njj2,
     '  nk,nkbf,nn,nnbf,
     '  ns,nse,nsf,NSP(-2:2),nu,NU1(0:3)
      REAL*8 AG,AGE,AGE1,AGE2,AGE3,AXU(3,3),AZ,AZL(3,3),
     '  AZU(3,3),CHTOFF(3,3,3),CLZ,CMZ,
     '  CSLZ,CSMZ,
     '  D11,D12,D13,D21,D22,D23,D31,D32,D33,
     '  Darcy_Resid,DBM(3,3,3),delsqP,
     '  DLA1,DLA2,DLA3,DMA1,DMA2,DMA3,
     '  DTA1,DTA2,DTA3,
     '  DXIX(3,3),DXIXN(3,3),DXIZN(3,3),DXNXI(3,3),DXNZN(3,3),
     '  dZ1_dNu1,dZ1_dNu2,dZ2_dNu1,dZ2_dNu2,DZNXI(3,3),DZNXN(3,3),
     '  E1,E2,EG(3,3),G1,GXL(3,3),GXU(3,3),GZ,GZL(3,3),GZU(3,3),PF(2),
     '  PGA1,PGA2,PGA3,PGG,PGX,PPG(64,4),PPGG(4),RF(32,3),RI1,RI2,RI3,
     '  RWG,RZWG,SLZ,SMZ,SUM,SUM1,SUM2,
     '  TIME1,TIME2,TIMER,TG(3,3),TNA,X3G(4,3),XI(3)
      CHARACTER ERROR_DUMMY*255,TYPE*9
      LOGICAL BCPARAM,ELEMPRESS,NODEPRESS,SAMEDEPBASIS
      DATA NU1/1,2,4,7/

      CALL ENTERS('ZERE50',*9999)
      NITB=NIT(NBJ(1))

c      CALL GETTIMER(IHANDLE,ERROR,*1111)

C *** Test whether to use unrolled loops
      SAMEDEPBASIS=.TRUE.
      nh1=NH_LOC(1,nx)
      DO nhx=2,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        IF(NBH(nh).NE.NBH(nh1)) SAMEDEPBASIS=.FALSE.
      ENDDO
      IF(ITYP10(nr).EQ.1.AND.KTYP51(nr).EQ.3.AND.KTYP53(nr).LE.3.AND.
     '  NJ_LOC(NJL_GEOM,0).EQ.3.AND.NIT(NBH(NH_LOC(1,nx))).EQ.3) THEN
        TYPE='RC3D'
      ELSE IF(ITYP10(nr).EQ.4.AND.KTYP51(nr).EQ.3.AND.
     '    KTYP53(nr).LE.3.AND.SAMEDEPBASIS.AND..NOT.DOP.AND.
     '    NJ_LOC(NJL_GEOM,0).EQ.3.AND.
     '    NIT(NBH(NH_LOC(1,nx))).EQ.3) THEN
        TYPE='PROLATE'
      ELSE IF(NJT.EQ.2.AND.KTYP51(nr).EQ.4) THEN
        TYPE='BIAXIAL'
      ELSE
        TYPE=' '
      ENDIF

      DO nhx=1,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
          RE(ns,nh)=0.0d0
        ENDDO
      ENDDO

C *** Set up PF array for pressure bcs from aux vars in ZE
      IF(KTYP57(nr).GT.1) THEN
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for pressure vars
        NODEPRESS=.FALSE.
        IF(NST(NBP).GT.0) NODEPRESS=.TRUE. !node based hyd press interp
        ELEMPRESS=.FALSE.
C       Put current XI3 face pressures into PF array
        NSP(1)=0
        NSP(2)=0
        DO na=1,NAT(NBP)
C         Check that pressure varies in Xi(3) dirn only for
C         element based pressure interpolation
          IF(NST(NBP).EQ.0.AND.  !elem based press interp
     '      (NAN(1,na,NBP).NE.0.OR.NAN(2,na,NBP).NE.0)) THEN
            ERROR='>>Pressure must vary only with Xi(3) '
     '        //'for elem based pressure interpolation'
            GOTO 9999
          ENDIF
          IF(NAN(3,na,NBP).EQ.0) THEN
C           Pick up param assoc with const pressure term
            NSP(1)=na
          ELSE IF(NAN(3,na,NBP).EQ.1.OR.NAN(3,na,NBP).EQ.3) THEN
C           Pick up param assoc with linear or cubic press term
            NSP(2)=na
          ELSE IF(NAN(3,na,NBP).EQ.-1) THEN
C           Pick up param assoc with Xi3=0 face pressure bc
            PF(1)=ZE(NST(NBP)+na,NH_LOC(0,nx))
            NSP(-1)=na
          ELSE IF(NAN(3,na,NBP).EQ.-2) THEN
C           Pick up param assoc with Xi3=1 face pressure bc
            PF(2)=ZE(NST(NBP)+na,NH_LOC(0,nx))
            NSP(-2)=na
          ENDIF
C         check for element based hyd press interp apart from
C         pressure boundary condition parameters
          BCPARAM=.FALSE.
          DO ni=1,NIT(NBP)
            IF(NAN(ni,na,NBP).LT.0) BCPARAM=.TRUE.
          ENDDO !ni
          IF(.NOT.BCPARAM) ELEMPRESS=.TRUE.
        ENDDO !na
      ENDIF

C *** Main Gauss point loop
      DO 50 ng=1,NGT(NBH(NH_LOC(1,nx)))

c        TIME1=TIMER(IHANDLE,T_INITIALISE)

        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Gauss pt '',I3)') NG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' ------------''/)')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C       Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
        CALL XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*9999)
        !18Sep88: Eventually, replace following 6 stmts with
        ! call xgmg(KTYP53(nr),...) where xgmg is modified so that DXIX
        ! are derivs of Xi wrt Xj if 1st arg=1 or wrt Nu if >1
        IF(KTYP53(nr).EQ.1) THEN
C         stresses are referred to Xj in constitutive law...
          JP=0
        ELSE IF(KTYP53(nr).GE.2) THEN
C         stresses referred to Nu in constitutive law...
          JP=1
        ENDIF
C       Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C       ..derivs (DXIX) of Xi wrt Xj (JP=0) or Nu (JP=1) coords.
        CALL XGMG(JP,NIT(NBJ(1)),NBJ(1),nr,DXIX,GXL,GXU,RGX(ng),XG,
     '    ERROR,*9999)

        IF(KTYP53(nr).EQ.3) THEN  !active stress
C         Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
          CALL DXIDNU(NBJ(1),nr,DXIXN,DXNXI,GXL,GXU,XG,ERROR,*9999)
C new MPN 4-May-96: new way of handling sheets
C         Get derivs of Xi wrt deformed Nu coords, DXIZN
          CALL DXIDZN(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,ERROR,*9999)
C old
CC         Interpolate dependent var.s ZG and derivs wrt Xi
C          CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C          CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)
CC         Get derivs of Xi wrt deformed Nu coords,DXIZN
C          CALL DXIDNU(NBJ(1),nr,DXIZN,DZNXI,GZL,GZU,XG,ERROR,*9999)
C end old
C         Calculate derivs of deformed Nu wrt undeformed Nu (DZNXN)
          DO ni=1,NITB
            DO mi=1,NITB
              SUM1=0.0d0
              SUM2=0.0d0
              DO k=1,NITB
                SUM1=SUM1+DZNXI(ni,k)*DXIXN(k,mi)
                SUM2=SUM2+DXNXI(ni,k)*DXIZN(k,mi)
              ENDDO
              DZNXN(ni,mi)=SUM1
              DXNZN(ni,mi)=SUM2
            ENDDO
          ENDDO
        ENDIF

C       Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
        CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C       Calculate deformed metric tensors wrt Nu (AZL,AZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),AZ,AZL,AZU,ZG,ERROR,*9999)

C       Get contravariant cpts of 2nd Piola-Kirchhoff stress
C       ..tensor (TG) wrt undeformed Nu coordinates
        IF(KTYP51(nr).EQ.1) THEN      !plane stress
          CALL ZGTG51(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '      CG(1,ng),RI1,RI2,RI3,TG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.2) THEN !plane strain
          CALL ZGTG52(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '      CG(1,ng),RI1,RI2,RI3,TG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.3) THEN !3D
          CALL ZGTG53(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '      CG(1,ng),EG,RI1,RI2,RI3,TG,XG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.4) THEN !membrane
          CALL ZGTG54(NBH(NH_LOC(1,nx)),nr,nx,AXU,AZ,AZL,AZU,
     '      CG(1,ng),EG,RI1,RI2,RI3,TG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.5) THEN !string
          CALL ZGTG55(nr,AZL,CG(1,ng),EG,TG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.6) THEN !shell
        ENDIF

C ***   Active stress component
        FEXT(1,ng)=DSQRT(AZL(1,1))
        IF(KTYP53(nr).EQ.3) THEN
          !09-Dec-1989: NOTE: Don't have to define a separate face array
          CALL ZGTG5A(nr,FEXT(1,ng),DXNZN,DZNXN,CG(IL_time_delay,ng),
     '      TG,TNA,ERROR,*9999)
        ENDIF

        RWG=RGX(ng)*WG(ng,NBH(NH_LOC(1,nx)))
        IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*XG(1,1)    !cyl. symm. about x
        IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*XG(2,1)    !cyl. symm. about y
        IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*XG(1,1)**2 !spherical symmetry

C ***   Main element residual
        IF(TYPE(1:7).EQ.'BIAXIAL') THEN !biaxial membrane testing

          IF(DOP) THEN
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          dZ1_dNu1=ZG(1,NU1(1))
          dZ1_dNu2=ZG(1,NU1(2))
          dZ2_dNu1=ZG(2,NU1(1))
          dZ2_dNu2=ZG(2,NU1(2))
          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)+NAT(nb)
            PGA1 = PG(ns,2,ng,nb)*DXIX(1,1) +
     '             PG(ns,4,ng,nb)*DXIX(2,1)
            PGA2 = PG(ns,2,ng,nb)*DXIX(1,2) +
     '             PG(ns,4,ng,nb)*DXIX(2,2)
            AGE1 = (TG(1,1)*dZ1_dNu1+TG(1,2)*dZ1_dNu2)*PGA1 +
     '             (TG(2,1)*dZ1_dNu1+TG(2,2)*dZ1_dNu2)*PGA2
            AGE2 = (TG(1,1)*dZ2_dNu1+TG(1,2)*dZ2_dNu2)*PGA1 +
     '             (TG(2,1)*dZ2_dNu1+TG(2,2)*dZ2_dNu2)*PGA2
            RE(ns,1)=RE(ns,1)+AGE1*RWG*SE(ns,nb,ne)
            RE(ns,2)=RE(ns,2)+AGE2*RWG*SE(ns,nb,ne)
          ENDDO

        ELSE IF(TYPE(1:4).EQ.'RC3D') THEN     !3D rect.cart.

          IF(DOP) THEN
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DO njj=1,NJ_LOC(NJL_GEOM,0)
            nh=NJ_LOC(NJL_GEOM,njj)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)

C             PPGG(1) = PG(ns,1,ng,nb)
              PPGG(2) = PG(ns,NU1(1),ng,nb)*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,1)
              PPGG(3) = PG(ns,NU1(1),ng,nb)*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,2)
              PPGG(4) = PG(ns,NU1(1),ng,nb)*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,3)

              AGE=(TG(1,1)*ZG(nh,2)+TG(1,2)*ZG(nh,4)+TG(1,3)*ZG(nh,7))
     '           *PPGG(2)
     '           +(TG(2,1)*ZG(nh,2)+TG(2,2)*ZG(nh,4)+TG(2,3)*ZG(nh,7))
     '           *PPGG(3)
     '           +(TG(3,1)*ZG(nh,2)+TG(3,2)*ZG(nh,4)+TG(3,3)*ZG(nh,7))
     '           *PPGG(4)

              RE(ns,nh)=RE(ns,nh)+AGE*RWG*SE(ns,nb,ne)
            ENDDO
          ENDDO

        ELSE IF(TYPE(1:7).EQ.'PROLATE') THEN  !3D prolate spheroidal

          IF(DOP) THEN
            WRITE(OP_STRING,'('' >>Using unrolled loops'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          SLZ=SINH(ZG(1,1))
          SMZ=SIN (ZG(2,1))
          CLZ=SQRT(1.0d0+SLZ*SLZ)
          CMZ=SQRT(1.0d0-SMZ*SMZ)
          CSLZ=CLZ/SLZ
          CSMZ=CMZ/SMZ
          G1=SLZ*SLZ+SMZ*SMZ
          E1=CLZ*SLZ/G1
          E2=CMZ*SMZ/G1

          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)+NAT(nb)
            PPG(ns,1) = PG(ns,1,ng,nb)
            PPG(ns,2) = PG(ns,NU1(1),ng,nb)*DXIX(1,1) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,1) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,1)
            PPG(ns,3) = PG(ns,NU1(1),ng,nb)*DXIX(1,2) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,2) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,2)
            PPG(ns,4) = PG(ns,NU1(1),ng,nb)*DXIX(1,3) +
     '                  PG(ns,NU1(2),ng,nb)*DXIX(2,3) +
     '                  PG(ns,NU1(3),ng,nb)*DXIX(3,3)
          ENDDO

          D11=ZG(1,NU1(1))
          D12=ZG(1,NU1(2))
          D13=ZG(1,NU1(3))
          D21=ZG(2,NU1(1))
          D22=ZG(2,NU1(2))
          D23=ZG(2,NU1(3))
          D31=ZG(3,NU1(1))
          D32=ZG(3,NU1(2))
          D33=ZG(3,NU1(3))

          DLA1=ZG(1,NU1(1))
          DLA2=ZG(1,NU1(2))
          DLA3=ZG(1,NU1(3))
          DMA1=ZG(2,NU1(1))
          DMA2=ZG(2,NU1(2))
          DMA3=ZG(2,NU1(3))
          DTA1=ZG(3,NU1(1))
          DTA2=ZG(3,NU1(2))
          DTA3=ZG(3,NU1(3))

          DO ns=1,NST(nb)+NAT(nb)
            PGG =PPG(ns,1)
            PGA1=PPG(ns,2)
            PGA2=PPG(ns,3)
            PGA3=PPG(ns,4)

            AGE1 = TG(1,1)*(D11*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '        D21*(E1*DMA1-E2*DLA1)*PGG+D31*E1*SMZ*SMZ*DTA1*PGG) +
     '        TG(2,1)*(D11*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '        D21*(E1*DMA2-E2*DLA2)*PGG+D31*E1*SMZ*SMZ*DTA2*PGG) +
     '        TG(3,1)*(D11*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '        D21*(E1*DMA3-E2*DLA3)*PGG+D31*E1*SMZ*SMZ*DTA3*PGG) +
     '        TG(1,2)*(D12*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '        D22*(E1*DMA1-E2*DLA1)*PGG+D32*E1*SMZ*SMZ*DTA1*PGG) +
     '        TG(2,2)*(D12*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '        D22*(E1*DMA2-E2*DLA2)*PGG+D32*E1*SMZ*SMZ*DTA2*PGG) +
     '        TG(3,2)*(D12*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '        D22*(E1*DMA3-E2*DLA3)*PGG+D32*E1*SMZ*SMZ*DTA3*PGG) +
     '        TG(1,3)*(D13*(PGA1-(E1*DLA1+E2*DMA1)*PGG) +
     '        D23*(E1*DMA1-E2*DLA1)*PGG+D33*E1*SMZ*SMZ*DTA1*PGG) +
     '        TG(2,3)*(D13*(PGA2-(E1*DLA2+E2*DMA2)*PGG) +
     '        D23*(E1*DMA2-E2*DLA2)*PGG+D33*E1*SMZ*SMZ*DTA2*PGG) +
     '        TG(3,3)*(D13*(PGA3-(E1*DLA3+E2*DMA3)*PGG) +
     '        D23*(E1*DMA3-E2*DLA3)*PGG+D33*E1*SMZ*SMZ*DTA3*PGG)

            AGE2 = TG(1,1)*(D11*(E2*DLA1-E1*DMA1)*PGG +
     '        D21*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D31*E2*SLZ*SLZ*DTA1*PGG)+
     '        TG(2,1)*(D11*(E2*DLA2-E1*DMA2)*PGG +
     '        D21*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D31*E2*SLZ*SLZ*DTA2*PGG)+
     '        TG(3,1)*(D11*(E2*DLA3-E1*DMA3)*PGG +
     '        D21*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D31*E2*SLZ*SLZ*DTA3*PGG)+
     '        TG(1,2)*(D12*(E2*DLA1-E1*DMA1)*PGG +
     '        D22*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D32*E2*SLZ*SLZ*DTA1*PGG)+
     '        TG(2,2)*(D12*(E2*DLA2-E1*DMA2)*PGG +
     '        D22*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D32*E2*SLZ*SLZ*DTA2*PGG)+
     '        TG(3,2)*(D12*(E2*DLA3-E1*DMA3)*PGG +
     '        D22*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D32*E2*SLZ*SLZ*DTA3*PGG)+
     '        TG(1,3)*(D13*(E2*DLA1-E1*DMA1)*PGG +
     '        D23*(PGA1-(E1*DLA1+E2*DMA1)*PGG)+D33*E2*SLZ*SLZ*DTA1*PGG)+
     '        TG(2,3)*(D13*(E2*DLA2-E1*DMA2)*PGG +
     '        D23*(PGA2-(E1*DLA2+E2*DMA2)*PGG)+D33*E2*SLZ*SLZ*DTA2*PGG)+
     '        TG(3,3)*(D13*(E2*DLA3-E1*DMA3)*PGG +
     '        D23*(PGA3-(E1*DLA3+E2*DMA3)*PGG)+D33*E2*SLZ*SLZ*DTA3*PGG)

            AGE3 = TG(1,1)*(-(D11*CSLZ+D21*CSMZ)*DTA1*PGG +
     '        D31*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '        TG(2,1)*(-(D11*CSLZ+D21*CSMZ)*DTA2*PGG +
     '        D31*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '        TG(3,1)*(-(D11*CSLZ+D21*CSMZ)*DTA3*PGG +
     '        D31*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '        TG(1,2)*(-(D12*CSLZ+D22*CSMZ)*DTA1*PGG +
     '        D32*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '        TG(2,2)*(-(D12*CSLZ+D22*CSMZ)*DTA2*PGG +
     '        D32*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '        TG(3,2)*(-(D12*CSLZ+D22*CSMZ)*DTA3*PGG +
     '        D32*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG)) +
     '        TG(1,3)*(-(D13*CSLZ+D23*CSMZ)*DTA1*PGG +
     '        D33*(PGA1-(CSLZ*DLA1+CSMZ*DMA1)*PGG)) +
     '        TG(2,3)*(-(D13*CSLZ+D23*CSMZ)*DTA2*PGG +
     '        D33*(PGA2-(CSLZ*DLA2+CSMZ*DMA2)*PGG)) +
     '        TG(3,3)*(-(D13*CSLZ+D23*CSMZ)*DTA3*PGG +
     '        D33*(PGA3-(CSLZ*DLA3+CSMZ*DMA3)*PGG))

            RE(ns,1)=RE(ns,1)+AGE1*RWG*SE(ns,nb,ne)
            RE(ns,2)=RE(ns,2)+AGE2*RWG*SE(ns,nb,ne)
            RE(ns,3)=RE(ns,3)+AGE3*RWG*SE(ns,nb,ne)

          ENDDO

        ELSE !all other cases

          DO njj1=1,NJ_LOC(NJL_GEOM,0)
            nh=NJ_LOC(NJL_GEOM,njj1)
            nb=NBH(nh)
            DO ns=1,NST(nb)+NAT(nb)
              PPG(ns,1)=PG(ns,1,ng,nb)
              DO njj2=1,NJ_LOC(NJL_GEOM,0)
                nj=NJ_LOC(NJL_GEOM,njj2)
                PPG(ns,1+nj)=PGX(nb,nj,ns,DXIX,PG(1,1,ng,nb))
              ENDDO !njj2
            ENDDO !ns
            DO ns=1,NST(nb)+NAT(nb)
              AGE=AG(nb,nh,nr,ns,PPG,TG,ZG)
              RE(ns,nh)=RE(ns,nh)+AGE*RWG*SE(ns,nb,ne)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' >>Residual calcs:'
     '            //' nh='',I2,'' ns='',I2,'' nb='',I2,'
     '            //''' AG='',D10.3,'' RE='',D10.3/,'
     '            //''' PG(ns,1,ng,nb):'',D10.3,'' RWG='',D10.3)')
     '            nh,ns,nb,AGE,RE(ns,nh),PG(ns,1,ng,nb),RWG
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !ns
          ENDDO !njj1

        ENDIF

C ***   Incompressibilty constraint in 3D case
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for press vars
        IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.2) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,'('' NBP='',I2,''  SQRT(RI3)='',D12.4)')
     '        NBP,DSQRT(RI3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' PG(1..,1,ng,nb='',I2,''):'',5D10.3'
     '        //'/(20X,5D10.3))')
     '        NBP,(PG(ns,1,ng,NBP),NS=1,NST(NBP)+NAT(NBP))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NODEPRESS.AND.ELEMPRESS) THEN
C           node and element based hydrostatic pressure interp
C           NOTE: pressure basis function is not used to
C                 weight the residual here (non-Galerkin)
            DO ns=1,NST(NBP)+NAT(NBP)
              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))+(DSQRT(RI3)-1.0d0)*RWG
            ENDDO
          ELSE
C           node (only) or element (only) based hyd. pressure
C           interp (standard Galerkin)
            DO ns=1,NST(NBP)+NAT(NBP)
              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))
     '          +(DSQRT(RI3)-1.0d0)*PG(ns,1,ng,NBP)*RWG
            ENDDO
          ENDIF

C ***   Incompressibilty + fluid constraint in 3D case
        ELSE IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.3) THEN
C         Calculate deformed Christoffel symbols wrt undef Nu coords
C         NOTE: ZG needs derivs wrt Nu not Xi !
          CALL TOFFEL(NBJ(1),NJE,nr,CHTOFF,DBM,AZU,ZG,X3G,.FALSE.,
     '      ERROR,*9999)
C         Calculate del-squared(p) wrt deformed coords, where p is the
C         hydrostatic pressure that varies with Xi(3) only.
          SUM=0.0d0
          DO j=1,NITB
            DO k=1,NITB
              SUM=SUM+CHTOFF(3,j,k)*AZU(j,k)
            ENDDO
          ENDDO
          delsqP=0.0d0
          DO ns=1,NST(NBP)+NAT(NBP)
            delsqP=delsqP+ZE(ns,NH_LOC(0,nx))*
     '        (PG(ns,8,ng,NBP)*AZU(3,3)-SUM*PG(ns,7,ng,NBP))
          ENDDO

C         Calculate the Jacobian for integration wrt def coords:
C         Interpolate dependent var.s ZG and derivs wrt Xi
          CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C         Calculate deformed metric tensors wrt Xi (GZL,GZU)
          CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)
C         Calc deformed Jacobian*Gaussian Quadrature weight
          RZWG=DSQRT(GZ)*WG(ng,NBH(NH_LOC(1,nx)))
          IF(JTYP4.EQ.2) RZWG=RZWG*2.0d0*PI*ZG(1,1)    !cyl sym about x
          IF(JTYP4.EQ.3) RZWG=RZWG*2.0d0*PI*ZG(2,1)    !cyl sym about y
        IF(JTYP4.EQ.4) RZWG=RZWG*4.0d0*PI*ZG(1,1)**2 !spherical sym

C ***     Calculate residuals associated with Darcy's Law
          Darcy_Resid=CG(IL_fluid_conductivity,ng)*DT*delsqP
     '      +(DSQRT(RI3)-1.d0)/DSQRT(RI3)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' NBP='',I2,''  SQRT(RI3)='',D12.4)')
     '        NBP,DSQRT(RI3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' del-squared(pressure)='',D12.4)')delsqP
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' def Jac * Gauss weight='',D12.4)')RZWG
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' PG(1..,1,ng,nb='',I2,''):'',5D10.3'
     '        //'/(20X,5D10.3))')
     '        NBP,(PG(ns,1,ng,NBP),NS=1,NST(NBP)+NAT(NBP))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Darcy residual='',D12.4)') Darcy_Resid
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NODEPRESS.AND.ELEMPRESS) THEN
C           node and element based hydrostatic pressure interp
C           NOTE: pressure basis function is not used to
C                 weight the residual here (non-Galerkin)
            DO ns=1,NST(NBP)+NAT(NBP)
              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))+Darcy_Resid*RZWG
            ENDDO
          ELSE
C           node (only) or element (only) based hyd. pressure
C           interp (standard Galerkin)
            DO ns=1,NST(NBP)+NAT(NBP)
              RE(ns,NH_LOC(NH_LOC(0,nx),nx))=
     '          RE(ns,NH_LOC(NH_LOC(0,nx),nx))+
     '          Darcy_Resid*PG(ns,1,ng,NBP)*RZWG
            ENDDO
          ENDIF

C ***   External press loads in case of 2D or 3D membrane
        ELSE IF(KTYP51(nr).EQ.4.AND.NW.GT.1) THEN !membrane
          DO nhx=1,NJ_LOC(NJL_GEOM,0)
            nh=NH_LOC(nhm,nx)
            nb=NBH(nh)
            SUM=0.0d0
            DO ni=1,NIT(nb)
              SUM=SUM+AZU(ni,3)*ZG(nhx,NU1(ni))
            ENDDO
            DO ns=1,NST(nb)+NAT(nb)
              RE(ns,nh)=RE(ns,nh)+PF(1)*SUM*PG(ns,1,ng,nb)*RWG
            ENDDO
          ENDDO

C ***   External press loads in case of 2D string
        ELSE IF(KTYP51(nr).EQ.5) THEN !2D string
          IF(NRT.EQ.2) THEN !special case for coupled sail-flow problem
            PF(1)=PRESS_DIFF_AERO(ng,ne-NET(1))
          ENDIF
C         For x residual use ZG(2,2)=dy/dNu
          nb=NBH(NH_LOC(1,nx))
          DO ns=1,NST(nb)+NAT(nb)
            RE(ns,1)=RE(ns,1)+PF(1)*ZG(2,2)*PG(ns,1,ng,nb)*RWG
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'('' ne='',I4,'' ng='',I3,'
     '        //''' PF(1)='',D12.3,'' ZG(2,2)='',D12.3,'
     '        //''' RE(ns,1):'',6D12.3)')
     '        ne,ng,PF(1),ZG(2,2),(RE(ns,1),ns=1,NST(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C         For y residual use ZG(1,2)=dx/dNu. Note -ve.
          nb=NBH(NH_LOC(2,nx))
          DO ns=1,NST(nb)+NAT(nb)
            RE(ns,2)=RE(ns,2)-PF(1)*ZG(1,2)*PG(ns,1,ng,nb)*RWG
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'('' ne='',I4,'' ng='',I3,'
     '        //''' PF(1)='',D12.3,'' ZG(1,2)='',D12.3,'
     '        //''' RE(ns,2):'',6D12.3)')
     '        ne,ng,PF(1),ZG(1,2),(RE(ns,2),ns=1,NST(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

c        TIME2=TIMER(IHANDLE,T_CPU)
        time1=0.0d0
        time2=0.0d0
        IF(DOP) THEN
          WRITE(OP_STRING,'('' All calcs for current gauss pt took'','
     '      //'D11.4,'' secs'')') TIME2-TIME1
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

 50   CONTINUE !end of ng loop

C *** External pressure loads in 3D case
      IF(KTYP51(nr).EQ.3.AND.NW.GE.2.AND.NW.LE.4) THEN
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for hydr. pressure
        IF(KTYP52(nr).GE.2.AND.KTYP57(nr).GT.1) THEN
C         Boundary pressure constraint equations
C         NOTE: if hydrostatic pressure in nodally (only) based
C               then there are no boundary pressure constraint eqns
          IF(NODEPRESS.AND.ELEMPRESS) THEN
C           node and element based hydrostatic pressure interp
            CALL PFRE_NP(IBT,IDO,INP,NAN,NBH,NBJ,NFF,NGAP,NHE,NJE,
     '        NKF,NNF,NPF,NPNE(1,1,ne),nr,NRE,NW,nx,NXI(-NIM,ne),
     '        CE,CG,CP,FEXT,PF,PG,RE(1,NH_LOC(NH_LOC(0,nx),nx)),
     '        VE,WG,XE,XG,ZE,ZG,ERROR,*9999)
          ELSE IF(ELEMPRESS) THEN
C           element (only) based hydrostatic pressure interp
            CALL PFRE_NE(IBT,IDO,INP,NAN,NBH,NBJ,NFF,NGAP,NHE,
     '        NPF,NPNE(1,1,ne),nr,NRE,NSP,NW,nx,NXI(-NIM,ne),
     '        CE,CP,FEXT,PF,PG,RE(1,NH_LOC(NH_LOC(0,nx),nx)),
     '        XE,XG,ZE,ZG,ERROR,*9999)
          ENDIF
        ENDIF

C       Contribution of external pressure loads to stress equilibrium
C       equation residuals at nodes
        DO iface=1,2
          IF(iface.EQ.1.AND.(NXI(-3,ne).EQ.0.OR.nr.NE.NRE(NXI(-3,ne)))
     '      .AND.(NW.EQ.2.OR.NW.EQ.4).OR. !ext press bc on Xi3=0 face
     '      iface.EQ.2.AND.(NXI(3,ne).EQ.0.OR.nr.NE.NRE(NXI(3,ne)))
     '      .AND.(NW.EQ.3.OR.NW.EQ.4)) THEN !ext pres bc on Xi3=1 face
            IF(DABS(PF(iface)).GT.1.0D-10) THEN !non-zero press bc
              IXF=iface-1
              IFE=iface+4
              nf=NFF(IFE)
              CALL PFRF(IBT,IDO,INP,IXF,NAN,NBH,NBJ,NGAP,NJE,NPF(1,nf),
     '          nr,nx,PF(iface),PG,RF,WG,XE,XG,ZE,ZG,ERROR,*9999)
              DO nhx=1,NJ_LOC(NJL_GEOM,0)
                nh=NH_LOC(nhx,nx)
                nsf=0
                NBFF=NPF(9+nh,nf)
                NBE=NBH(nh)
                DO nnbf=1,NNT(NBFF)
                  nn=NNF(1+nnbf,IFE,NBE)
                  DO nkbf=1,NKT(nnbf,NBFF)
                    nk=NKF(nkbf,nnbf,IFE,NBE)
                    nse=nk+(NN-1)*NKT(nn,NBE)
                    nsf=nsf+1
                    RE(nse,nh)=RE(nse,nh)-RF(nsf,nh)*SE(nse,NBE,ne)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' >>ZERE50 diagnostic op:'
     '                    //' RF('',I2,'','',I1,''):'',D10.3,'
     '                  //''' RE('',I2,'','',I1,''):'',D10.3)')
     '                  nse,nh,RF(nsf,nh),
     '                  nsf,nh,RE(nse,nh)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO !nkbf (nk)
                ENDDO !nnbf (nn)
              ENDDO !nhx
            ENDIF
          ENDIF
        ENDDO !iface
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' Residuals for element'',I4,'':'')') ne
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,'('' RE(ns,'',I1,''): '',8D11.3,'
     '      //'/(11X,8D11.3))')
     '      nh,(RE(ns,nh),ns=1,NST(NBH(nh))+NAT(NBH(nh)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

c      CALL FREETIMER(IHANDLE,ERROR,*1111)
      CALL EXITS('ZERE50')
      RETURN
 9999 CALL FREETIMER(IHANDLE,ERROR_DUMMY,*1111)
 1111 CALL ERRORS('ZERE50',ERROR)
      CALL EXITS('ZERE50')
      RETURN 1
      END


      SUBROUTINE ZERE55(INP,NBH,ne,NHE,nr,nx,
     '  CG,PG,RE,WG,ZE,ZEREF,ZG,ERROR,*)

C#### Subroutine: ZERE55
C###  Description:
C###    ZERE55 calculates constant-volume elements for cavity loading.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER INP(NNM,NIM,NBFM),NBH(NHM),ne,NHE,nr,nx
      REAL*8 CG(NMM,NGM),PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),
     '  WG(NGM,NBM),ZE(NSM,NHM),ZEREF(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IXI1,IXI2,IXI3,nb,ng,nh,nhx,nk,nn,ns
      REAL*8 deltaMU,DXIX(3,3),GZ,GZL(3,3),GZU(3,3),
     '  kstif,press_current,press_init,
     '  press_sign,RGZ,RGZREF,RWGZ,RWGZREF,VOLDEF,VOLUND
      LOGICAL FOUNDnn

      CALL ENTERS('ZERE55',*9999)

      VOLDEF=0.0d0  !  Deformed element volume
      VOLUND=0.0d0  !Undeformed element volume

      DO nhx=1,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
          RE(ns,nh)=0.0d0
        ENDDO !ns
      ENDDO !nh

      CALL ASSERT(KTYP52(nr).EQ.2,'Must be incompressible',ERROR,*9999)

C *** Main Gauss point loop
      DO ng=1,NGT(NBH(NH_LOC(1,nx)))
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Gauss pt '',I3)') NG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' ------------''/)')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

C!!! requires a different reference state than the ventricular wall
CC ***   Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
C        CALL XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*9999)
CC ***   Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
CC ***   ..derivs (DXIX) of Xi wrt Xj (JP=0) coords.
C        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GXL,GXU,RGX(ng),XG,
C     '    ERROR,*9999)
C        RWG=RGX(ng)*WG(ng,NBH(NH_LOC(1,nx)))
C        IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*XG(1,1)    !cyl symm about x
C        IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*XG(2,1)    !cyl symm about y
C        IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*XG(1,1)**2 !spherical symmetry

C ***   Interpolate dependent var.s ZG and derivs wrt Xi from ZEREF
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZEREF,ZG,ERROR,*9999)
C ***   Calculate deformed metric tensors wrt Xi (GZL,GZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)
        RGZREF=DSQRT(GZ)
        RWGZREF=RGZREF*WG(ng,NBH(NH_LOC(1,nx)))
        IF(JTYP4.EQ.2) RWGZREF=RWGZREF*2.0d0*PI*ZG(1,1)    !cyl symm (x)
        IF(JTYP4.EQ.3) RWGZREF=RWGZREF*2.0d0*PI*ZG(2,1)    !cyl symm (y)
        IF(JTYP4.EQ.4) RWGZREF=RWGZREF*4.0d0*PI*ZG(1,1)**2 !sph symmetry
        press_init=ZG(NH_LOC(NH_LOC(0,nx),nx),1) !initial cav. press.

C ***   Interpolate dependent var.s ZG and derivs wrt Xi from ZE
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C ***   Calculate deformed metric tensors wrt Xi (GZL,GZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),GZ,GZL,GZU,ZG,ERROR,*9999)
        RGZ=DSQRT(GZ)
        RWGZ=RGZ*WG(ng,NBH(NH_LOC(1,nx)))
        IF(JTYP4.EQ.2) RWGZ=RWGZ*2.0d0*PI*ZG(1,1)    !cyl symm about x
        IF(JTYP4.EQ.3) RWGZ=RWGZ*2.0d0*PI*ZG(2,1)    !cyl symm about y
        IF(JTYP4.EQ.4) RWGZ=RWGZ*4.0d0*PI*ZG(1,1)**2 !spherical symmetry
        press_current=ZG(NH_LOC(NH_LOC(0,nx),nx),1) !current cav. press.

C ***   pressure-displ constraints on lambda,mu,theta
        kstif=CG(1,ng) !GENERALISE!!!
        IF(DOP) THEN
          WRITE(OP_STRING,'('' >>mu stiffness='',D12.4,'
     '      //''' initial cavity pressure='',D12.4,'
     '      //''' current cavity pressure='',D12.4,'
     '      //''' RWGZREF='',D12.4,'' RWGZ='',D12.4)')
     '      kstif,press_init,press_current,RWGZREF,RWGZ
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Residual calcs:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !DOP
        DO nhx=1,NJ_LOC(NJL_GEOM,0) !loop over geom vars
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh) !basis fn for nh
          DO IXI3=1,2
            DO IXI2=1,2
              IF(IXI2.EQ.1) THEN
                press_sign= 1.0d0
              ELSE IF(IXI2.EQ.2) THEN
                press_sign=-1.0d0
              ENDIF
              DO IXI1=1,2
                ns=0
                nn=1
                FOUNDnn=.FALSE.
                DO WHILE(.NOT.FOUNDnn.AND.nn.LE.NNT(nb))
                  IF(INP(nn,1,nb).EQ.IXI1.AND.INP(nn,2,nb).EQ.IXI2
     '              .AND.INP(nn,3,nb).EQ.IXI3) THEN
                    FOUNDnn=.TRUE.
                    DO nk=1,NKT(nn,nb)
                      ns=ns+1
                      deltaMU=ZE(ns,nhx)-ZEREF(ns,nhx) !delta MU;deriv
                      RE(ns,nh)=RE(ns,nh)+(kstif*deltaMU+press_sign*
     '                  (press_current-press_init))*
     '                  PG(ns,1,ng,nb)*RWGZREF
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' nh='',I2,'' ns='',I2,'
     '                    //''' deltaMU='',D12.4,'
     '                    //''' PG(ns,1,ng,nb):'',D12.4,'
     '                    //''' Resid='',D12.4)')
     '                    nh,ns,deltaMU,PG(ns,1,ng,nb),(kstif*deltaMU+
     '                    press_sign*(press_current-press_init))
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF !DOP
                    ENDDO !nk
                  ELSE
                    DO nk=1,NKT(nn,nb)
                      ns=ns+1
                    ENDDO !nk
                    nn=nn+1
                  ENDIF
                ENDDO !while...
              ENDDO !IXI1
            ENDDO !IXI2
          ENDDO !IXI3
        ENDDO !nhx

        VOLUND=VOLUND+RWGZREF
        VOLDEF=VOLDEF+RWGZ

      ENDDO !ng

C *** Const volume constraint
      nb=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for hyd. pressure
      DO ns=1,NST(nb)+NAT(nb)
        RE(ns,NH_LOC(NH_LOC(0,nx),nx))=VOLDEF-VOLUND
      ENDDO !ns

      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' Residuals for element'',I4,'':'')') ne
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,'('' RE(ns,'',I1,''): '',8D11.3,'
     '      //'/(11X,8D11.3))')
     '      nh,(RE(ns,nh),ns=1,NST(NBH(nh))+NAT(NBH(nh)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C MPN 28-Aug-95: old code for pressure/vol relationship
C      RE(1,1)=VOLDEF
C      RE(2,1)=VOLUND
C      DVOL=VOLDEF-VOLUND
C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' Volume change for element '','
C     '    //'I4,'' is '',D12.3)') NE,DVOL
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Total volumes for element '','
C     '    //'I4,'' is '',D12.3,'' (deformed) and '',D12.3,'
C     '    //''' (undeformed)'')') NE,VOLDEF,VOLUND
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF

      CALL EXITS('ZERE55')
      RETURN
 9999 CALL ERRORS('ZERE55',ERROR)
      CALL EXITS('ZERE55')
      RETURN 1
      END


      SUBROUTINE ZETX50(COORDS,CSTYPE,IBT,IDO,INP,NAN,NBH,NBJ,ng,
     '  NHE,NJE,NPNE,nr,nx,CE,CG,CP,FEXT,PG,PHI,PST,
     '  RGX,RGX2D,RGZ,RGZ2D,RM,TC,TG,TN,VE,XE,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: ZETX50
C###  Description:
C###    ZETX50 calculates 2nd Piola-Kirchhoff, Cauchy and  Nominal
C###    stresses with respect to 'Reference' or 'Fibre' coords
C###    (as specified by COORDS) at position XI in current element
C###    if ng=0 else at Gauss point ng.

C**** Since base vectors of theta coords are not orthonormal, cpts
C**** of Cauchy and Nominal stress are converted to 'physical' values.
C**** AZ,AZL,AZU  are deformed metric tensors wrt undeformed coords
C**** RI1,RI2,RI3 are principal invariants of AZL
C**** AXU  are contravariant cpts of undeformed metric tensor
C**** XG   are undeformed theta coords and derivs wrt Xi
C**** ZG   are deformed theta coords and derivs wrt undeformed coords
C**** TG   are tensor cpts of 2nd Piola-Kirchhoff stresses
C**** TN   are physical cpts of Nominal stresses
C**** TC   are physical cpts of Cauchy stress
C**** PST  are principal stresses
C**** RM   is the modal matrix whose cols are the eigenvectors
C****      associated with PST
C**** PHI  are the Euler angles of the principal stresses
C**** AZLZ are deformed Eulerian metrics (wrt deformed theta coords)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbst02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ipma50.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,NJE,NPNE(NNM,NBFM),
     '  nr,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),FEXT(NIFEXTM),
     '  PG(NSM,NUM,NGM,NBM),PHI(3),PST(3),RGX,RGX2D,RGZ,RGZ2D,RM(3,3),
     '  TC(3,3),TG(3,3),TN(3,3),VE(NSM,NKM),XE(NSM,NJM),XG(NJM,NUM),
     '  XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER COORDS*(*),CSTYPE*(*),ERROR*(*)
!    Local Variables
      INTEGER i,IFAIL,il,j,k,mi,mj,mjj,mix,mz,
     '  nb,NCW,ni,NITB,nj,njj,NU1(0:3),nix,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AA,AXU(3,3),AZ,AZL(3,3),AZLZ(3,3),AZU(3,3),CW(NCW),
     '  DET,DET_DZDX,DXDZ(3,3),DXIXJ(3,3),DXIXN(3,3),DXIZN(3,3),
     '  DXJXN(3,3),DXNXI(3,3),DXNXJ(3,3),DZDX(3,3),DZNXI(3,3),EG(3,3),
     '  G1,G3,GX2D,GXL(3,3),GXU(3,3),GZ,GZ2D,GZL(3,3),GZU(3,3),
     '  RC,RI1,RI2,RI3,RR,RWX,SLX,SMX,SUM,TGX(3,3),TNA,TOL,VALTMP,
     '  WK1_LOCAL(3)
      DATA NU1/1,2,4,7/

      CALL ENTERS('ZETX50',*9999)
      CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '  //' must equal dimension of CE array (NMMX)',ERROR,*9999)
      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)

      IF(ng.EQ.0) THEN
C       Interpolate material parameters at XI
        CALL CPXI(1,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CW,XI,ERROR,*9999)
      ELSE
C       Put Gauss pt params into CW array
        DO il=1,ILT(1,nr,nx)
          CW(il)=CG(il,ng)
        ENDDO !il
      ENDIF

      IF(KTYP53(nr).EQ.1) THEN !stress ref'ed to Xj in constitutive law
        IF(ng.EQ.0) THEN
C ***     Interpolate midwall geometric var.s XG and derivs wrt Xi
          CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
        ELSE
C ***     Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
          CALL XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*9999)
        ENDIF
C ***   Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C ***   derivatives of Xi wrt Xj (reference) coords, DXIXJ (IP=0)
        CALL XGMG(0,NITB,NBJ(1),nr,DXIXJ,GXL,GXU,RGX,XG,ERROR,*9999)
C ***   Calculate 2D Jacobian wrt undef coords for face integrals
        GX2D=GXL(1,1)*GXL(2,2)-GXL(1,2)*GXL(2,1)
        RGX2D=DSQRT(GX2D*GXU(3,3))
        IF(DOP) THEN
          WRITE(OP_STRING,'('' RGX='',D12.4,'' RGX2D='',D12.4)')
     '      RGX,RGX2D
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C ***   Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
        CALL DXIDNU(NBJ(1),nr,DXIXN,DXNXI,GXL,GXU,XG,ERROR,*9999)
        IF(ng.EQ.0) THEN
C ***     Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
          CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '      ERROR,*9999)
        ELSE
C ***     Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
          CALL ZEZG(1,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
        ENDIF
C ***   Calculate deformed metric tensors wrt Xj (AZL,AZU)
        CALL ZGMG(nb,AZ,AZL,AZU,ZG,ERROR,*9999)
C ***   Get contravariant cpts of 2nd Piola-Kirchhoff stress
C ***   tensor (TG) wrt (undeformed) Xj coordinates
        IF(KTYP51(nr).EQ.1) THEN
          CALL ZGTG51(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.2) THEN
          CALL ZGTG52(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.3) THEN
          CALL ZGTG53(nb,nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TG,XG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.4) THEN
          CALL ZGTG54(nb,nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TG,ERROR,*9999)
        ENDIF
        IF(COORDS(1:5).EQ.'Refer') THEN
C ***     Put partial derivs of deformed theta wrt undef Xj into DZDX
          CALL DLZJDX(1,nb,NJE,DZDX,XG,ZG,ERROR,*9999)
          CALL INVERT(NITB,DZDX,DXDZ,DET_DZDX)
        ELSE IF(COORDS(1:5).EQ.'Fibre'.OR.COORDS(1:5).EQ.'Princ') THEN
          IF(ng.EQ.0) THEN
C ***       Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
            CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '        ERROR,*9999)
          ELSE
C ***       Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
            CALL ZEZG(0,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
          ENDIF
C ***     Calculate deformed metric tensors wrt Xi (GZL,GZU)
          CALL ZGMG(nb,GZ,GZL,GZU,ZG,ERROR,*9999)
          RGZ=DSQRT(GZ)
C ***     Calculate 2D Jacobian wrt def coords for face integrals
          GZ2D=GZL(1,1)*GZL(2,2)-GZL(1,2)*GZL(2,1)
          RGZ2D=DSQRT(GZ2D*GZU(3,3))
          IF(DOP) THEN
            WRITE(OP_STRING,'('' RGZ='',D12.4,'' RGZ2D='',D12.4)')
     '        RGZ,RGZ2D
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C new MPN 4-May-96: new way of handling sheets
C         Get derivs of Xi wrt deformed Nu coords, DXIZN
          CALL DXIDZN(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,ERROR,*9999)
C old
CC ***     Get derivs of Xi wrt deformed Nu coords, DXIZN
C          CALL DXIDNU(NBJ(1),nr,DXIZN,DZNXI,GZL,GZU,XG,ERROR,*9999)
C ***     Calculate derivs of deformed Nu wrt undeformed Nu (DZDX)
          DO ni=1,NITB
            DO mi=1,NITB
              SUM=0.0d0
              DO k=1,NITB
                SUM=SUM+DZNXI(ni,k)*DXIXN(k,mi)
              ENDDO !k
              DZDX(ni,mi)=SUM
            ENDDO !mi
          ENDDO !ni
          CALL INVERT(NITB,DZDX,DXDZ,DET_DZDX)
          IF(ng.EQ.0) THEN
C ***       Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
            CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXN,ZE,ZG,XI,
     '        ERROR,*9999)
          ELSE
C ***       Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
            CALL ZEZG(1,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
          ENDIF
C ***     Calculate deformed metric tensors wrt Nu (AZL,AZU)
          CALL ZGMG(nb,AZ,AZL,AZU,ZG,ERROR,*9999)
C ***     Calculate derivs of undeformed Nu wrt undeformed Xj (DXNXJ)
          DO ni=1,NITB
            DO mjj=1,NJ_LOC(NJL_GEOM,0)
              mj=NJ_LOC(NJL_GEOM,mjj)
              SUM=0.0d0
              DO k=1,NITB
                SUM=SUM+DXNXI(ni,k)*DXIXJ(k,mj)
              ENDDO
              DXNXJ(ni,mj)=SUM
            ENDDO
          ENDDO
C ***     Transform TG to (undeformed) Nu coords
          DO ni=1,NITB
            DO mi=1,NITB
              SUM=0.0d0
              DO njj=1,NJ_LOC(NJL_GEOM,0)
                nj=NJ_LOC(NJL_GEOM,njj)
                DO mjj=1,NJ_LOC(NJL_GEOM,0)
                  mj=NJ_LOC(NJL_GEOM,mjj)
                  SUM=SUM+DXNXJ(ni,nj)*DXNXJ(mi,mj)*TG(nj,mj)
                ENDDO
              ENDDO
              TGX(ni,mi)=SUM
            ENDDO
          ENDDO
          DO ni=1,NITB
            DO mi=1,NITB
              TG(ni,mi)=TGX(ni,mi)
            ENDDO
          ENDDO
        ENDIF

      ELSE IF(KTYP53(nr).GT.1) THEN !stress ref'ed to Nu in constit law
        IF(ng.EQ.0) THEN
C ***     Interpolate midwall geometric var.s XG and derivs wrt Xi
          CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
        ELSE
C ***     Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
          CALL XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*9999)
        ENDIF
C ***   Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C ***   derivatives of Xi wrt Xj (reference) coords, DXIXJ (IP=0)
        CALL XGMG(0,NITB,NBJ(1),nr,DXIXJ,GXL,GXU,RWX,XG,
     '    ERROR,*9999)
C ***   Calculate 2D Jacobian wrt undef coords for face integrals
        GX2D=GXL(1,1)*GXL(2,2)-GXL(1,2)*GXL(2,1)
        RGX2D=DSQRT(GX2D*GXU(3,3))
        IF(DOP) THEN
          WRITE(OP_STRING,'('' RGX='',D12.4,'' RGX2D='',D12.4)')
     '      RGX,RGX2D
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C ***   Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
        CALL DXIDNU(NBJ(1),nr,DXIXN,DXNXI,GXL,GXU,XG,ERROR,*9999)
        IF(COORDS(1:5).EQ.'Fibre'.OR.COORDS(1:5).EQ.'Princ'
     '     .OR.KTYP53(nr).EQ.3) THEN
          IF(ng.EQ.0) THEN
C ***       Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
            CALL ZEZW(0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '        ERROR,*9999)
          ELSE
C ***       Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
            CALL ZEZG(0,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
          ENDIF
C ***     Calculate deformed metric tensors wrt Xi (GZL,GZU)
          CALL ZGMG(nb,GZ,GZL,GZU,ZG,ERROR,*9999)
          RGZ=DSQRT(GZ)
C ***     Calculate 2D Jacobian wrt def coords for face integrals
          GZ2D=GZL(1,1)*GZL(2,2)-GZL(1,2)*GZL(2,1)
          RGZ2D=DSQRT(GZ2D*GZU(3,3))
          IF(DOP) THEN
            WRITE(OP_STRING,'('' RGZ='',D12.4,'' RGZ2D='',D12.4)')
     '        RGZ,RGZ2D
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C new MPN 4-May-96: new way of handling sheets
C         Get derivs of Xi wrt deformed Nu coords, DXIZN
          CALL DXIDZN(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '      DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,ERROR,*9999)
C old
CC ***     Get derivs of Xi wrt deformed Nu coords, DXIZN
C          CALL DXIDNU(NBJ(1),nr,DXIZN,DZNXI,GZL,GZU,XG,ERROR,*9999)
C end old
C ***     Calculate derivs of deformed Nu wrt undeformed Nu (DZDX)
          DO ni=1,NITB
            DO mi=1,NITB
              SUM=0.0d0
              DO k=1,NITB
                SUM=SUM+DZNXI(ni,k)*DXIXN(k,mi)
              ENDDO !k
              DZDX(ni,mi)=SUM
            ENDDO !mi
          ENDDO !ni
          CALL INVERT(NITB,DZDX,DXDZ,DET_DZDX)
        ENDIF
        IF(ng.EQ.0) THEN
C ***     Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
          CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXN,ZE,ZG,XI,
     '      ERROR,*9999)
        ELSE
C ***     Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
          CALL ZEZG(1,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
        ENDIF
C ***   Calculate deformed metric tensors wrt Nu (AZL,AZU)
        CALL ZGMG(nb,AZ,AZL,AZU,ZG,ERROR,*9999)
C ***   Get contravariant cpts of 2nd Piola-Kirchhoff stress
C ***   tensor (TG) wrt (undeformed) Xj coordinates
        IF(KTYP51(nr).EQ.1) THEN !plane stress
          CALL ZGTG51(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.2) THEN !plane strain
          CALL ZGTG52(nb,nr,nx,AXU,AZ,AZL,AZU,CW,
     '      RI1,RI2,RI3,TG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.3) THEN !3D
          CALL ZGTG53(nb,nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TG,XG,ZG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.4) THEN !membrane
          CALL ZGTG54(nb,nr,nx,AXU,AZ,AZL,AZU,CW,EG,
     '      RI1,RI2,RI3,TG,ERROR,*9999)
        ELSE IF(KTYP51(nr).EQ.5) THEN !string
          CALL ZGTG55(nr,AZL,CW,EG,TG,ERROR,*9999)
        ENDIF
        IF(KTYP53(nr).EQ.3) THEN
C ***     Add active fibre stress component
          FEXT(1)=DSQRT(AZL(1,1))
          CALL ZGTG5A(nr,FEXT,DXDZ,DZDX,CW(IL_time_delay),
     '      TG,TNA,ERROR,*9999)
        ENDIF !KTYP53(nr)=3 (active stresses)

        IF(COORDS(1:5).EQ.'Refer') THEN
          IF(ng.EQ.0) THEN
C ***       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
            CALL ZEZW(1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '        ERROR,*9999)
          ELSE
C ***       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
            CALL ZEZG(1,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
          ENDIF
C ***     Put partial derivs of deformed theta wrt undef Xj into DZDX
          CALL DLZJDX(1,nb,NJE,DZDX,XG,ZG,ERROR,*9999)
          DET_DZDX=DET(DZDX)
C ***     Calculate deformed metric tensors wrt Xj (AZL,AZU)
          CALL ZGMG(nb,AZ,AZL,AZU,ZG,ERROR,*9999)
C ***     Calculate derivs of undeformed Xj wrt undeformed Nu (DXJXN)
          DO njj=1,NJ_LOC(NJL_GEOM,0)
            nj=NJ_LOC(NJL_GEOM,njj)
            DO mi=1,NITB
              SUM=0.0d0
              DO k=1,NITB
                SUM=SUM+DXIXN(k,mi)*XG(nj,NU1(k))
              ENDDO
              DXJXN(nj,mi)=SUM
            ENDDO
          ENDDO
C ***     Transform TG to (undeformed) Xj coords
          DO njj=1,NJ_LOC(NJL_GEOM,0)
            nj=NJ_LOC(NJL_GEOM,njj)
            DO mjj=1,NJ_LOC(NJL_GEOM,0)
              mj=NJ_LOC(NJL_GEOM,mjj)
              SUM=0.0d0
              DO ni=1,NITB
                DO mi=1,NITB
                  SUM=SUM+DXJXN(nj,ni)*DXJXN(mj,mi)*TG(ni,mi)
                ENDDO
              ENDDO
              TGX(nj,mj)=SUM
            ENDDO
          ENDDO
          DO ni=1,NITB
            DO mi=1,NITB
              TG(ni,mi)=TGX(ni,mi)
            ENDDO
          ENDDO
        ENDIF !coords=reference
      ENDIF !KTYP53(nr) (reference axes for stresses)

      IF(KTYP51(nr).NE.5) THEN !plane stress,strain/3D/membrane/shell
        DO mz=1,NITB
          DO nz=1,NITB
            SUM=0.0d0
            DO mix=1,NITB
              DO nix=1,NITB
                SUM=SUM+DZDX(mz,mix)*TG(mix,nix)*DZDX(nz,nix)
              ENDDO !nix
            ENDDO !mix
            TC(mz,nz)=SUM/DET_DZDX
          ENDDO !nz
        ENDDO !mz
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TC:'',12X,3D12.4,/(16X,3D12.4))')
     '      ((TC(mz,nz),nz=1,3),mz=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nix=1,NITB
          DO nz=1,NITB
            SUM=0.0d0
            DO mix=1,NITB
              SUM=SUM+TG(nix,mix)*DZDX(nz,mix)
            ENDDO !mix
            TN(nix,nz)=SUM
          ENDDO !nz
        ENDDO !nix
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TN:'',12X,3D12.4,/(16X,3D12.4))')
     '      ((TN(mz,nz),nz=1,3),mz=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(COORDS(1:5).EQ.'Refer') THEN
C ***     Compute undeformed metric tensor wrt theta coordinates AXU
C ***     and Eulerian metric tensor AZLZ (i.e. wrt deformed theta)
          DO i=1,3
            DO j=1,3
              AXU (i,j)=0.0d0
              AZLZ(i,j)=0.0d0
            ENDDO
            AXU (i,i)=1.0d0
            AZLZ(i,i)=1.0d0
          ENDDO
          IF(NITB.EQ.2) THEN
            IF(ITYP10(nr).EQ.2) THEN
              RR=XG(1,1)**2
              AXU(2,2)=1.0d0/RR
              AZLZ(2,2)=ZG(1,1)**2
            ELSE IF(ITYP10(nr).EQ.3) THEN
              RR=XG(1,1)**2
              AXU(2,2)=1.0d0/RR
              AXU(3,3)=1.0d0
              RR=ZG(1,1)**2
              AZLZ(3,3)=1.0d0
              AZLZ(2,2)=RR
            ELSE IF(ITYP10(nr).EQ.4) THEN
              AA=FOCUS*FOCUS
              SLX=DSINH(XG(1,1))
              SMX=DSIN (XG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AXU(1,1)=1.0d0/G1
              AXU(2,2)=1.0d0/G1
              AXU(3,3)=1.0d0
              SLX=DSINH(ZG(1,1))
              SMX=DSIN (ZG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AZLZ(1,1)=G1
              AZLZ(2,2)=G1
              AZLZ(3,3)=1.0d0
            ELSE IF(ITYP10(nr).EQ.5) THEN
            ENDIF
          ELSE IF(NITB.EQ.3) THEN
            IF(ITYP10(nr).EQ.2) THEN
              RR=XG(1,1)**2
              AXU(2,2)=1.0d0/RR
              AZLZ(2,2)=ZG(1,1)**2
            ELSE IF(ITYP10(nr).EQ.3) THEN
              RR=XG(1,1)**2
              RC=RR*DCOS(XG(3,1))**2
              AXU(2,2)=1.0d0/RC
              AXU(3,3)=1.0d0/RR
              RR=ZG(1,1)**2
              AZLZ(3,3)=RR
              AZLZ(2,2)=RR*DCOS(ZG(3,1))**2
            ELSE IF(ITYP10(nr).EQ.4) THEN
              AA=FOCUS*FOCUS
              SLX=DSINH(XG(1,1))
              SMX=DSIN (XG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AXU(1,1)=1.0d0/G1
              AXU(2,2)=1.0d0/G1
              AXU(3,3)=1.0d0/G3
              SLX=DSINH(ZG(1,1))
              SMX=DSIN (ZG(2,1))
              G1=AA*(SLX*SLX+SMX*SMX)
              G3=AA* SLX*SLX*SMX*SMX
              AZLZ(1,1)=G1
              AZLZ(2,2)=G1
              AZLZ(3,3)=G3
            ELSE IF(ITYP10(nr).EQ.5) THEN
            ENDIF
          ENDIF
C ***     Compute physical cpts of Cauchy & Nominal stresses
          DO njj=1,NJ_LOC(NJL_GEOM,0)
            nj=NJ_LOC(NJL_GEOM,njj)
            DO mjj=1,NJ_LOC(NJL_GEOM,0)
              mj=NJ_LOC(NJL_GEOM,mjj)
              TC(nj,mj)=TC(nj,mj)*DSQRT(AZLZ(nj,nj)*AZLZ(mj,mj))
              TN(nj,mj)=TN(nj,mj)*DSQRT(AZLZ(mj,mj)/ AXU(nj,nj))
            ENDDO
          ENDDO
        ENDIF
C ***   Compute principal Stresses, Euler angles and eigenvectors
        IFAIL=0
        IF(CSTYPE(1:5).EQ.'Piola') THEN
          DO ni=1,NITB
            DO mi=1,NITB
              TGX(ni,mi)=TG(ni,mi)
            ENDDO
          ENDDO
          CALL F02ABF(TGX,3,NITB,PST,RM,3,WK1_LOCAL,IFAIL)
        ELSE IF(CSTYPE(1:5).EQ.'Nomin') THEN
          CALL F02ABF(TN,3,NITB,PST,RM,3,WK1_LOCAL,IFAIL)
        ELSE IF(CSTYPE(1:5).EQ.'Cauch') THEN
          CALL F02ABF(TC,3,NITB,PST,RM,3,WK1_LOCAL,IFAIL)
        ENDIF
!news MPN 30-Jun-94
C ***   Order PST from max->min and change cols of RM accordingly
        IF(PST(2).GT.PST(1)) THEN
          VALTMP=PST(1)
          PST(1)=PST(2)
          PST(2)=VALTMP
          DO ni=1,NITB
            VALTMP=RM(ni,1)
            RM(ni,1)=RM(ni,2)
            RM(ni,2)=VALTMP
          ENDDO
        ENDIF
        IF(NITB.EQ.3) THEN
          DO mi=1,2
            IF(PST(3).GT.PST(MI)) THEN
              VALTMP=PST(MI)
              PST(MI)=PST(3)
              PST(3)=VALTMP
              DO ni=1,NITB
                VALTMP=RM(ni,mi)
                RM(ni,mi)=RM(ni,3)
                RM(ni,3)=VALTMP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
!newe
        IF(NITB.EQ.2) THEN
          IF(DABS(RM(1,1)).LE.1.0d0) THEN
            PHI(1)=DACOS(RM(1,1))
          ELSE
            PHI(1)=0.0d0
          ENDIF
          IF(DABS(RM(1,2)).LE.1.0d0) THEN
            PHI(2)=DACOS(RM(1,2))
          ELSE
            PHI(2)=0.0d0
          ENDIF
        ELSE IF(NITB.EQ.3) THEN
          TOL=1.0d-08
          IF(DABS(RM(1,1)).GT.TOL) THEN
            PHI(1)=DATAN2(RM(2,1),RM(1,1))
          ELSE IF(RM(2,1)*RM(1,1).GT.0.0d0) THEN
            PHI(1)=90.0d0
          ELSE
            PHI(1)=-90.0d0
          ENDIF
          IF(DABS(RM(3,1)).LE.1.0d0) THEN
            PHI(2)=DASIN(RM(3,1))
          ELSE IF(RM(3,1).GT.1.0d0) THEN
            PHI(2)=90.0d0
          ELSE
            PHI(2)=-90.0d0
          ENDIF
           IF(DABS(DCOS(PHI(1))).GT.TOL .AND.
     '       DABS(DCOS(PHI(1))).GE.DABS(RM(3,3))) THEN
             PHI(3)=DACOS(RM(3,3)/DCOS(PHI(1)))
          ELSE
            PHI(3)=0.0d0
          ENDIF
        ENDIF
        DO ni=1,NITB
          PHI(ni)=PHI(ni)*180.0d0/PI
          IF(PHI(ni).GT.90.0d0) PHI(ni)=PHI(ni)-180.0d0
          IF(PHI(ni).LT.-90.0d0) PHI(ni)=PHI(ni)+180.0d0
C         Find Max/Min principal stresses
          IF(PST(ni).GT.PRSTMAX) PRSTMAX=PST(ni)
          IF(PST(ni).LT.PRSTMIN) PRSTMIN=PST(ni)
        ENDDO
      ENDIF !KTYP51(nr)=plane stress,strain/3D/membrane/shell

      CALL EXITS('ZETX50')
      RETURN
 9999 CALL ERRORS('ZETX50',ERROR)
      CALL EXITS('ZETX50')
      RETURN 1
      END


      SUBROUTINE ZGTG51(nb,nr,nx,AXU,AZ,AZL,AZU,CG,RI1,RI2,RI3,
     '  TG,ZG,ERROR,*)

C#### Subroutine: ZGTG51
C###  Description:
C###    ZGTG51 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for plane stress problems
C###    (KTYP51(nr)=1).
C**** Note: Stress components should be referred to (undeformed)
C****       Nu (body/fibre) coords (KTYP53(nr)>1).
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER nb,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),RI1,RI2,RI3,
     '  TG(3,3),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IFAIL,j,mi,mj,ni,NITB,nj
      REAL*8 AXL(3,3),BG(3,3),DW(6),EDG(3,3),EG(3,3),EG12,EG13,
     '  EVAL(3),RK1,RK2,RL1,RL2,RL3,RM(3,3),TGP(3,3),VALTMP,W3,
     '  WK1_LOCAL(3)

      CALL ENTERS('ZGTG51',*9999)

      IF(KTYP53(nr).EQ.1) THEN
        ERROR='Stresses must be referred to body/fibre coords'
        GOTO 9999
      ENDIF

      NITB=NIT(nb)
      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
      ENDDO

      RI3= AZ
      RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
      RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
      IF(KTYP55(nr).EQ.1) THEN
C ***   Invariants K1,K2 are fns of Physical strain that are
C ***   isotropic in the plane transverse to the Nu(1) axis
        EG13=AZL(1,3)/2.0d0
        EG12=AZL(1,2)/2.0d0
        RK1=(AZL(1,1)-1.0d0)/2.0d0
        RK2=EG13*EG13+EG12*EG12
        BG(1,1)=RI1-AZL(1,1)
        BG(2,2)=RI1-AZL(2,2)
        BG(3,3)=RI1-AZL(3,3)
        BG(2,1)=   -AZL(2,1)
        BG(3,1)=   -AZL(3,1)
        BG(3,2)=   -AZL(3,2)
        CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN      !compressible
          W3=RI3*DW(3)
        ELSE IF(KTYP52(nr).GT.1) THEN !incompressible
          W3=ZG(NH_LOC(0,nx),1)
        ENDIF
        DO mj=1,NITB
          DO nj=1,mj
            TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)
     '        +W3*AZU(mj,nj))
          ENDDO
        ENDDO
        TG(3,3)=0.0d0
        TG(3,1)=0.0d0
        TG(3,2)=0.0d0
        TG(2,1)=TG(2,1)+EG12*DW(5)
        TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
        TG(1,2)=TG(2,1)
        TG(1,3)=TG(3,1)
        TG(2,3)=TG(3,2)
      ELSE IF(KTYP55(nr).EQ.2) THEN
C ***   TG is a function of the principal stretches.
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,3
          DO j=1,3
            EDG(i,j)=0.5d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' EDG: '',9D12.4)') ((EDG(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IFAIL=0
        CALL F02ABF(EDG,3,NITB,EVAL,RM,3,WK1_LOCAL,IFAIL)
!news MPN 30-Jun-94
C ***   Order EVAL from max->min and change cols of RM accordingly
        IF(EVAL(2).GT.EVAL(1)) THEN
          VALTMP=EVAL(1)
          EVAL(1)=EVAL(2)
          EVAL(2)=VALTMP
          DO ni=1,NITB
            VALTMP=RM(ni,1)
            RM(ni,1)=RM(ni,2)
            RM(ni,2)=VALTMP
          ENDDO
        ENDIF
        IF(NITB.EQ.3) THEN
          DO mi=1,2
            IF(EVAL(3).GT.EVAL(MI)) THEN
              VALTMP=EVAL(MI)
              EVAL(MI)=EVAL(3)
              EVAL(3)=VALTMP
              DO ni=1,NITB
                VALTMP=RM(ni,mi)
                RM(ni,mi)=RM(ni,3)
                RM(ni,3)=VALTMP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
!newe
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EVAL:'',7X,(3X,3D11.3))')
     '      (EVAL(i),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RM:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((RM(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        RL1=DSQRT(2.0d0*EVAL(1)+1.0d0)
        RL2=DSQRT(2.0d0*EVAL(2)+1.0d0)
        RL3=DSQRT(2.0d0*EVAL(3)+1.0d0)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' RL1,RL2,RL3:'',(3X,3D11.3))')
     '      RL1,RL2,RL3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO i=1,3
          DO j=1,3
            TGP(i,j)=0.0d0
          ENDDO
        ENDDO
        CALL ENERGY(nr,CG,DW,RL1,RL2,RL3,0.0d0,0.0d0,0.0d0,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN      !compressible
          TGP(1,1)=DW(1)/RL1
          TGP(2,2)=DW(2)/RL2
          TGP(3,3)=DW(3)/RL3
        ELSE IF(KTYP52(nr).GT.1) THEN !incompressible
          TGP(1,1)=DW(1)/RL1+2.0d0*ZG(NH_LOC(0,nx),1)/RL1**2
          TGP(2,2)=DW(2)/RL2+2.0d0*ZG(NH_LOC(0,nx),1)/RL2**2
          TGP(3,3)=DW(3)/RL3+2.0d0*ZG(NH_LOC(0,nx),1)/RL3**2
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TGP:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((TGP(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL TRAN(3,TGP,TG,RM)
        TG(3,3)=0.0d0
        TG(1,3)=0.0d0
        TG(3,1)=0.0d0
        TG(2,3)=0.0d0
        TG(3,2)=0.0d0
      ELSE IF(KTYP55(nr).EQ.3) THEN
C ***   TG is a function of the fibre and transverse strains
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,3
          DO j=1,3
            EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EG: '',9D12.4)')
     '      ((EG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '    EG(2,3),ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN      !compressible
          W3=0.0d0 !Check this
        ELSE IF(KTYP52(nr).GT.1) THEN !incompressible
          W3=ZG(NH_LOC(0,nx),1)
        ENDIF
        DO i=1,2
          TG(i,i)=DW(i)/AXL(i,i)+2.0d0*W3*AZU(i,i)
        ENDDO
        TG(2,1)=DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+2.0d0*W3*AZU(2,1)
        TG(3,1)=0.0d0
        TG(3,2)=0.0d0
        TG(3,3)=0.0d0
        TG(1,2)=TG(2,1)
        TG(1,3)=TG(3,1)
        TG(2,3)=TG(3,2)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' TG:'',12X,2D12.4,/(16X,2D12.4))')
     '    ((TG(i,j),j=1,2),i=1,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(4,'('' TG:'',12X,2D12.4,/(16X,2D12.4))')
     '    ((TG(i,j),j=1,2),i=1,2)
      ENDIF

      CALL EXITS('ZGTG51')
      RETURN
 9999 CALL ERRORS('ZGTG51',ERROR)
      CALL EXITS('ZGTG51')
      RETURN 1
      END


      SUBROUTINE ZGTG52(nb,nr,nx,AXU,AZ,AZL,AZU,CG,RI1,RI2,RI3,
     '  TG,ZG,ERROR,*)

C#### Subroutine: ZGTG52
C###  Description:
C###    ZGTG52 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for plane stress problems
C###    (KTYP51(nr)=1).
C**** Note: Stress components should be referred to (undeformed)
C****       Nu (body/fibre) coords (KTYP53(nr)>1).
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER nb,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),RI1,RI2,RI3,
     '  TG(3,3),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IFAIL,j,mi,mj,ni,NITB,nj
      REAL*8 AXL(3,3),BG(3,3),DW(6),EDG(3,3),EG(3,3),EG12,EG13,EVAL(3),
     '  RK1,RK2,RL1,RL2,RL3,RM(3,3),TGP(3,3),VALTMP,W3,WK1_LOCAL(3)

      CALL ENTERS('ZGTG52',*9999)

      IF(KTYP53(nr).EQ.1) THEN
        ERROR='Stresses must be referred to body/fibre coords'
        GOTO 9999
      ENDIF

      NITB=NIT(nb)
      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
      ENDDO

      RI1=AZL(1,1)+AZL(2,2)
      RI2=AZ
C     RI3=(AZU(1,1)+AZU(2,2))*RI3

      IF(KTYP55(nr).EQ.1) THEN
C ***   Invariants K1,K2 are fns of Physical strain that are
C ***   isotropic in the plane transverse to the Nu(1) axis
        EG13=0.0d0
        EG12=AZL(1,2)/2.0d0
        RK1=(AZL(1,1)-1.0d0)/2.0d0
        RK2=EG13*EG13+EG12*EG12
        BG(1,1)=RI1-AZL(1,1)
        BG(2,2)=RI1-AZL(2,2)
        BG(3,3)=RI1-1.0d0
        BG(2,1)=-AZL(2,1)
        BG(3,1)=0.0d0
        BG(3,2)=0.0d0
        CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          W3=RI3*DW(3)
        ELSE IF(KTYP52(nr).GT.1) THEN
          W3=ZG(NH_LOC(0,nx),1)
        ENDIF
        DO mj=1,3
          DO nj=1,mj
            TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)
     '        +W3*AZU(mj,nj))
          ENDDO
        ENDDO
        TG(3,1)=TG(3,1)+EG13*DW(5)
        TG(2,1)=TG(2,1)+EG12*DW(5)
        TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
        TG(1,2)=TG(2,1)
        TG(1,3)=TG(3,1)
        TG(2,3)=TG(3,2)

      ELSE IF(KTYP55(nr).EQ.2) THEN
C ***   TG is a function of the principal stretches.
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,2
          DO j=1,2
            EDG(i,j)=0.5d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' EDG: '',9D12.4)') ((EDG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IFAIL=0
        CALL F02ABF(EDG,3,NITB,EVAL,RM,3,WK1_LOCAL,IFAIL)
!news MPN 30-Jun-94
C ***   Order EVAL from max->min and change cols of RM accordingly
        IF(EVAL(2).GT.EVAL(1)) THEN
          VALTMP=EVAL(1)
          EVAL(1)=EVAL(2)
          EVAL(2)=VALTMP
          DO ni=1,NITB
            VALTMP=RM(ni,1)
            RM(ni,1)=RM(ni,2)
            RM(ni,2)=VALTMP
          ENDDO
        ENDIF
        IF(NITB.EQ.3) THEN
          DO mi=1,2
            IF(EVAL(3).GT.EVAL(MI)) THEN
              VALTMP=EVAL(MI)
              EVAL(MI)=EVAL(3)
              EVAL(3)=VALTMP
              DO ni=1,NITB
                VALTMP=RM(ni,mi)
                RM(ni,mi)=RM(ni,3)
                RM(ni,3)=VALTMP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
!newe
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EVAL:'',7X,(3X,3D11.3))')
     '      (EVAL(i),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RM:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((RM(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        RL1=DSQRT(2.0d0*EVAL(1)+1.0d0)
        RL2=DSQRT(2.0d0*EVAL(2)+1.0d0)
        RL3=1.0d0
        IF(DOP) THEN
          WRITE(OP_STRING,'('' RL1,RL2,RL3:'',(3X,3D11.3))')
     '      RL1,RL2,RL3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO i=1,3
          DO j=1,3
            TGP(i,j)=0.0d0
          ENDDO
        ENDDO
        CALL ENERGY(nr,CG,DW,RL1,RL2,RL3,0.0d0,0.0d0,0.0d0,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          TGP(1,1)=DW(1)/RL1
          TGP(2,2)=DW(2)/RL2
          TGP(3,3)=DW(3)/RL3
        ELSE IF(KTYP52(nr).GT.1) THEN
          TGP(1,1)=DW(1)/RL1+2.0d0*ZG(NH_LOC(0,nx),1)/RL1**2
          TGP(2,2)=DW(2)/RL2+2.0d0*ZG(NH_LOC(0,nx),1)/RL2**2
          TGP(3,3)=DW(3)/RL3+2.0d0*ZG(NH_LOC(0,nx),1)/RL3**2
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TGP:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((TGP(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL TRAN(3,TGP,TG,RM)

      ELSE IF(KTYP55(nr).EQ.3) THEN
C ***   TG is a function of the fibre and transverse strains
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
            EG(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,2
          DO j=1,2
            EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EG: '',9D12.4)')
     '      ((EG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '    EG(2,3),ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          W3=0.0d0 !Check this
        ELSE
          W3=ZG(NH_LOC(0,nx),1)
        ENDIF
        DO i=1,3
          TG(i,i)=DW(i)/AXL(i,i)+2.0d0*W3*AZU(i,i)
        ENDDO
        TG(2,1)=DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+2.0d0*W3*AZU(2,1)
        TG(3,1)=DW(5)/DSQRT(AXL(1,1)*AXL(3,3))+2.0d0*W3*AZU(3,1)
        TG(3,2)=DW(6)/DSQRT(AXL(2,2)*AXL(3,3))+2.0d0*W3*AZU(3,2)
        TG(1,2)=TG(2,1)
        TG(1,3)=TG(3,1)
        TG(2,3)=TG(3,2)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),J=1,3),I=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(4,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),J=1,3),I=1,3)
      ENDIF

      CALL EXITS('ZGTG52')
      RETURN
 9999 CALL ERRORS('ZGTG52',ERROR)
      CALL EXITS('ZGTG52')
      RETURN 1
      END


      SUBROUTINE ZGTG53(nb,nr,nx,AXU,AZ,AZL,AZU,CG,EG,RI1,RI2,RI3,
     '  TG,XG,ZG,ERROR,*)

C#### Subroutine: ZGTG53
C###  Description:
C###    ZGTG53 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for 3-dimensional problems
C###    (KTYP51(nr)=3). Stress components are wrt (undeformed) theta
C###    coordinates (KTYP53(nr)=1) or (undeformed) Nu (body/fibre)
C###    coordinates (KTYP53(nr)>1).

C**** Note: For transversely isotropic strain energy functions of K1,K2
C****       the isotropic plane is transverse to theta(3) (KTYP53(nr)=1)
C****       or to Nu(1) (KTYP53(nr)>1)
C****       For transv. isotropic strain energy functions of strains
C****       EG, the isotropic plane is transverse to theta(1)
C****       (KTYP53(nr)=1) or to Nu(1) (KTYP53(nr)>1)
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER nb,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),EG(3,3),RI1,RI2,RI3,
     '  TG(3,3),XG(NJM,NUM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IFAIL,j,m,mi,mj,n,ni,NITB,nj
      REAL*8 AA,AXL(3,3),AZL_tmp(3,3),BG(3,3),DW(6),EDG(3,3),
     '  EG12,EG13,EG31,EG32,EVAL(3),Fgrowth(3,3),G1,G3,
     '  RC,RK1,RK2,RL1,RL2,RL3,RL33,RM(3,3),RR,SLX,SMX,TGP(3,3),
     '  VALTMP,W3,WK1_LOCAL(3)

      CALL ENTERS('ZGTG53',*9999)
      NITB=NIT(nb)

      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
      ENDDO

      IF(KTYP53(nr).EQ.1) THEN !stresses referred to theta coords
        IF(ITYP10(nr).EQ.1) THEN
          RI3= AZ
          RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
          RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
        ELSE IF(ITYP10(nr).EQ.2) THEN
          RR=XG(1,1)**2
          RI3=AZ/RR
          RI1= AZL(1,1)+AZL(2,2)/RR+AZL(3,3)
          RI2=(AZU(1,1)+AZU(2,2)*RR+AZU(3,3))*RI3
          AXU(2,2)=1.0d0/RR
        ELSE IF(ITYP10(nr).EQ.3) THEN
          RR=XG(1,1)**2
          RC=RR*DCOS(XG(3,1))**2
          RI3=AZ/(RR*RC)
          RI1= AZL(1,1)+AZL(2,2)/RC+AZL(3,3)/RR
          RI2=(AZU(1,1)+AZU(2,2)*RC+AZU(3,3)*RR)*RI3
          AXU(2,2)=1.0d0/RC
          AXU(3,3)=1.0d0/RR
        ELSE IF(ITYP10(nr).EQ.4) THEN
          AA=FOCUS*FOCUS
          SLX=DSINH(XG(1,1))
          SMX=DSIN (XG(2,1))
          G1=AA*(SLX*SLX+SMX*SMX)
          G3=AA* SLX*SLX*SMX*SMX
          RI3=AZ/(G1*G1*G3)
          RI1= (AZL(1,1)+AZL(2,2))/G1+AZL(3,3)/G3
          RI2=((AZU(1,1)+AZU(2,2))*G1+AZU(3,3)*G3)*RI3
          AXU(1,1)=1.0d0/G1
          AXU(2,2)=1.0d0/G1
          AXU(3,3)=1.0d0/G3
        ELSE IF(ITYP10(nr).EQ.5) THEN
        ENDIF

      ELSE IF(KTYP53(nr).GT.1) THEN
        !stress referred to body/fibre Nu coords
        RI3= AZ
        RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
        RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
      ENDIF

      IF(KTYP55(nr).EQ.1) THEN !TG is a fn of the principal invariants
        IF(KTYP53(nr).EQ.1) THEN !isotropic
          !Invariants K1,K2 are fns of Physical strain that are
          !isotropic in the plane transverse to the theta(3) axis
          IF(ITYP10(nr).EQ.1) THEN
            EG31=AZL(3,1)/2.0d0
            EG32=AZL(3,2)/2.0d0
            RK1=(AZL(3,3)-1.0d0)/2.0d0
            RK2=EG31*EG31+EG32*EG32
            BG(1,1)=RI1-AZL(1,1)
            BG(2,2)=RI1-AZL(2,2)
            BG(3,3)=RI1-AZL(3,3)
            BG(2,1)=   -AZL(2,1)
            BG(3,1)=   -AZL(3,1)
            BG(3,2)=   -AZL(3,2)
          ELSE IF(ITYP10(nr).EQ.2) THEN
            EG31=AZL(3,1)/2.0d0
            EG32=AZL(3,2)/(2.0d0*RR)
            RK1=(AZL(3,3)-1.0d0)/2.0d0
            RK2=(AZL(1,3)*EG31+AZL(2,3)*EG32)/2.0d0
            BG(1,1)= RI1-AZL(1,1)
            BG(2,2)=(RI1-AZL(2,2)/RR)/RR
            BG(3,3)= RI1-AZL(3,3)
            BG(2,1)=    -AZL(2,1)/RR
            BG(3,1)=    -AZL(3,1)
            BG(3,2)=    -AZL(3,2)/RR
          ELSE IF(ITYP10(nr).EQ.3) THEN
            EG31=AZL(3,1)/(2.0d0*RR)
            EG32=AZL(3,2)/(2.0d0*RR*RC)
            RK1=(AZL(3,3)-RR)/(2.0d0*RR)
            RK2=(AZL(1,3)*EG31+AZL(2,3)*EG32)/2.0d0
            BG(1,1)= RI1-AZL(1,1)
            BG(2,2)=(RI1-AZL(2,2)/RC)/RC
            BG(3,3)=(RI1-AZL(3,3)/RR)/RR
            BG(2,1)=    -AZL(2,1)/RC
            BG(3,1)=    -AZL(3,1)/RR
            BG(3,2)=    -AZL(3,2)/(RR*RC)
          ELSE IF(ITYP10(nr).EQ.4) THEN
            EG31=AZL(3,1)/(2.0d0*G1*G3)
            EG32=AZL(3,2)/(2.0d0*G1*G3)
            RK1=(AZL(3,3)-G3)/(2.0d0*G3)
            RK2=(AZL(1,3)*EG31+AZL(2,3)*EG32)/2.0d0
            BG(1,1)=(RI1-AZL(1,1)/G1)/G1
            BG(2,2)=(RI1-AZL(2,2)/G1)/G1
            BG(3,3)=(RI1-AZL(3,3)/G3)/G3
            BG(2,1)=    -AZL(2,1)/(G1*G1)
            BG(3,1)=    -AZL(3,1)/(G1*G3)
            BG(3,2)=    -AZL(3,2)/(G1*G3)
          ELSE IF(ITYP10(nr).EQ.5) THEN
          ENDIF
          CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,ERROR,*9999)
          IF(KTYP52(nr).EQ.1) THEN
            W3=RI3*DW(3)
          ELSE IF(KTYP52(nr).GT.1) THEN
            IF(KTYP51(nr).EQ.3) THEN
              W3=ZG(NH_LOC(0,nx),1)
            ELSE IF(KTYP51(nr).EQ.4) THEN
              RL33=1.0d0 !!!!!!!!!!!! Check this PJH 15Jan90
              W3=-(DW(1)+DW(2)*BG(3,3))*RL33
            ENDIF
          ENDIF
          DO mj=1,NITB
            DO nj=1,mj
              TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)
     '          +W3*AZU(mj,nj))
            ENDDO
          ENDDO
          TG(3,1)=TG(3,1)+EG31*DW(5)
          TG(3,2)=TG(3,2)+EG32*DW(5)
          TG(3,3)=TG(3,3)+AXU(3,3)*DW(4)
          TG(1,2)=TG(2,1)
          TG(1,3)=TG(3,1)
          TG(2,3)=TG(3,2)

        ELSE IF(KTYP53(nr).GT.1) THEN !aeolotropic
          !Invariants K1,K2 are fns of Physical strain that are
          !isotropic in the plane transverse to the Nu(1) axis
          EG13=AZL(1,3)/2.0d0
          EG12=AZL(1,2)/2.0d0
          RK1=(AZL(1,1)-1.0d0)/2.0d0
          RK2=EG13*EG13+EG12*EG12
          BG(1,1)=RI1-AZL(1,1)
          BG(2,2)=RI1-AZL(2,2)
          BG(3,3)=RI1-AZL(3,3)
          BG(2,1)=   -AZL(2,1)
          BG(3,1)=   -AZL(3,1)
          BG(3,2)=   -AZL(3,2)
          CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,ERROR,*9999)
          IF(KTYP52(nr).EQ.1) THEN      !compressible
            W3=RI3*DW(3)
          ELSE IF(KTYP52(nr).GT.1) THEN !incompressible
            W3=ZG(NH_LOC(0,nx),1)
          ENDIF
          DO mj=1,NITB
            DO nj=1,mj
              TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)
     '          +W3*AZU(mj,nj))
            ENDDO
          ENDDO
          TG(3,1)=TG(3,1)+EG13*DW(5)
          TG(2,1)=TG(2,1)+EG12*DW(5)
          TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
          TG(1,2)=TG(2,1)
          TG(1,3)=TG(3,1)
          TG(2,3)=TG(3,2)
        ENDIF
!news  21-MAY-1991 Printing strain invariants for transv isotropy. JSW
        IF(DOP) THEN
          WRITE(OP_STRING,'('' RK1= '',D12.4,'' RK2= '',D12.4)')
     '      RK1,RK2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(KTYP53(nr).EQ.1) THEN
            WRITE(OP_STRING,'('' G3= '',D12.4)')G3
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
!newe

      ELSE IF(KTYP55(nr).EQ.2) THEN
C       TG is a function of the princ stretches
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,3
          DO j=1,3
            EDG(i,j)=0.5d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' EDG: '',9D12.4)') ((EDG(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IFAIL=0
        CALL F02ABF(EDG,3,NITB,EVAL,RM,3,WK1_LOCAL,IFAIL)
!news MPN 30-Jun-94
C ***   Order EVAL from max->min and change cols of RM accordingly
        IF(EVAL(2).GT.EVAL(1)) THEN
          VALTMP=EVAL(1)
          EVAL(1)=EVAL(2)
          EVAL(2)=VALTMP
          DO ni=1,NITB
            VALTMP=RM(ni,1)
            RM(ni,1)=RM(ni,2)
            RM(ni,2)=VALTMP
          ENDDO
        ENDIF
        IF(NITB.EQ.3) THEN
          DO mi=1,2
            IF(EVAL(3).GT.EVAL(MI)) THEN
              VALTMP=EVAL(MI)
              EVAL(MI)=EVAL(3)
              EVAL(3)=VALTMP
              DO ni=1,NITB
                VALTMP=RM(ni,mi)
                RM(ni,mi)=RM(ni,3)
                RM(ni,3)=VALTMP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
!newe
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EVAL:'',7X,(3X,3D11.3))')
     '      (EVAL(i),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RM:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((RM(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        RL1=DSQRT(2.0d0*EVAL(1)+1.0d0)
        RL2=DSQRT(2.0d0*EVAL(2)+1.0d0)
        RL3=DSQRT(2.0d0*EVAL(3)+1.0d0)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' RL1,RL2,RL3:'',(3X,3D11.3))')
     '      RL1,RL2,RL3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO i=1,3
          DO j=1,3
            TGP(i,j)=0.0d0
          ENDDO
        ENDDO
        CALL ENERGY(nr,CG,DW,RL1,RL2,RL3,0.0d0,0.0d0,0.0d0,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          TGP(1,1)=DW(1)/RL1
          TGP(2,2)=DW(2)/RL2
          TGP(3,3)=DW(3)/RL3
        ELSE IF(KTYP52(nr).GT.1) THEN
          TGP(1,1)=DW(1)/RL1+2.0d0*ZG(NH_LOC(0,nx),1)/RL1**2
          TGP(2,2)=DW(2)/RL2+2.0d0*ZG(NH_LOC(0,nx),1)/RL2**2
          TGP(3,3)=DW(3)/RL3+2.0d0*ZG(NH_LOC(0,nx),1)/RL3**2
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TGP:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((TGP(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL TRAN(3,TGP,TG,RM)

      ELSE IF(KTYP55(nr).EQ.3) THEN !TG=func of fibre and transv strains
C       Compute undeformed covariant metric tensor
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO !j
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO !i
C!!! NOTE can only use resid strains for pole zero law until init extn
C!!!      mat params are set up for other problem types
        IF(KTYP56(nr).EQ.3) THEN !pole zero law
C         Calc growth defm tens for resid strain and copy AZL to AZL_tmp
          DO i=1,3
            DO j=1,3
              Fgrowth(i,j)=0.0d0 !growth defm tensor for resid strains
              AZL_tmp(i,j)=AZL(i,j)
            ENDDO !j
            Fgrowth(i,i)=CG(27+i) !NOTE init extns CG(28),CG(29),CG(30)
          ENDDO !i
C         Apply growth defm to deformed covariant metric tensor AZL
          DO i=1,3
            DO j=1,3
              AZL(i,j)=0.0d0
              DO m=1,3
                DO n=1,3
                  AZL(i,j)=AZL(i,j)+
     '              Fgrowth(i,m)*AZL_tmp(m,n)*Fgrowth(n,j)
                ENDDO !n
              ENDDO !m
            ENDDO !j
          ENDDO !i
C         Recompute AZU,AZ from transformed AZL
          CALL INVERT(NIT(nb),AZL,AZU,AZ)
        ENDIF !pole zero law
C       Compute Green strain tensor wrt material fibre coords
        DO i=1,3
          DO j=1,3
            EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO !j
        ENDDO !i
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EG: '',3D12.4,/(5X,3D12.4))')
     '      ((EG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '    EG(2,3),ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN !compressible
          W3=0.0d0 !Check this
        ELSE                 !incompressible
          W3=ZG(NH_LOC(0,nx),1)
        ENDIF
        DO i=1,3
          TG(i,i)=DW(i)/AXL(i,i)+2.0d0*W3*AZU(i,i)
        ENDDO
        TG(1,2)=DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+2.0d0*W3*AZU(1,2)
        TG(1,3)=DW(5)/DSQRT(AXL(1,1)*AXL(3,3))+2.0d0*W3*AZU(1,3)
        TG(2,3)=DW(6)/DSQRT(AXL(2,2)*AXL(3,3))+2.0d0*W3*AZU(2,3)
        TG(2,1)=TG(1,2)
        TG(3,1)=TG(1,3)
        TG(3,2)=TG(2,3)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),j=1,3),i=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ZGTG53')
      RETURN
 9999 CALL ERRORS('ZGTG53',ERROR)
      CALL EXITS('ZGTG53')
      RETURN 1
      END


      SUBROUTINE ZGTG54(nb,nr,nx,AXU,AZ,AZL,AZU,CG,EG,RI1,RI2,RI3,
     '  TG,ERROR,*)

C#### Subroutine: ZGTG54
C###  Description:
C###    ZGTG54 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for membrane problems
C###    (KTYP51(nr)=4).
C**** Stress components must always be with respect to
C**** (undeformed) Nu (body/fibre) coords (KTYP53(nr)>1).
C**** Note: For transversely isotropic strain energy functions of K1,K2
C****       the isotropic plane is transverse to Nu(1) (KTYP53(nr)>1)
C****       For transv isotropic strain energy functions of strains EG,
C****       the isotropic plane is transverse to Nu(1) (KTYP53(nr)>1)
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER nb,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),EG(3,3),RI1,RI2,RI3,
     '  TG(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IFAIL,j,m,mi,mj,n,ni,NITB,nj
      REAL*8 AXL(3,3),AZL_tmp(3,3),BG(3,3),DW(6),EG12,EG13,EVAL(3),
     '  Fgrowth(3,3),H,RK1,RK2,RL1,RL2,RL3,RM(3,3),
     '  TGP(3,3),VALTMP,W3,WK1_LOCAL(3)

      CALL ENTERS('ZGTG54',*9999)
      NITB=NIT(nb)

      DO i=1,3
        DO j=1,3
          AXL(i,j)=0.0d0
          AXU(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
        AXL(i,i)=1.0d0
      ENDDO

      IF(KTYP53(nr).EQ.1) THEN      !stresses in theta coords
        ERROR='>>Stresses must be wrt body/fibre coords'
        GO TO 9999
      ELSE IF(KTYP53(nr).GT.1) THEN !stresses in body/fibre Nu coords
        AZL(3,3)=1.0d0/AZ !AZ=azl(1,1)*azl(2,2)-azl(1,2)*azl(2,1)
        AZU(3,3)=AZ
        AZL(1,3)=0.0d0
        AZL(2,3)=0.0d0
        AZL(3,1)=0.0d0
        AZL(3,2)=0.0d0
C       RI3= AZ    !????? I3=det(AZL)=AZ/AZ=1.0d0
        RI3=1.0d0    !13Aug90
        RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
        RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
      ENDIF

      IF(KTYP55(nr).EQ.1) THEN !TG is fn of the principal invariants
        IF(KTYP53(nr).GT.1) THEN
C         Invariants K1,K2 are fns of Physical strain that are
C         isotropic in the plane transverse to the Nu(1) axis
          EG13=AZL(1,3)/2.0d0
          EG12=AZL(1,2)/2.0d0
          RK1=(AZL(1,1)-1.0d0)/2.0d0
          RK2=EG13*EG13+EG12*EG12
          BG(1,1)=RI1-AZL(1,1)
          BG(2,2)=RI1-AZL(2,2)
          BG(3,3)=RI1-AZL(3,3)
          BG(2,1)=   -AZL(2,1)
          BG(3,1)=   -AZL(3,1)
          BG(3,2)=   -AZL(3,2)
          CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,ERROR,*9999)
          IF(KTYP52(nr).EQ.1) THEN !compressible
            W3=RI3*DW(3)
          ELSE                 !incompressible
            W3=-(DW(1)+DW(2)*BG(3,3))*AZL(3,3)  !hydro press=2*W3
          ENDIF
          DO mj=1,NITB
            DO nj=1,mj
              TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)+W3
     '          *AZU(mj,nj))
            ENDDO
          ENDDO
          TG(3,3)=0.0d0
          TG(3,1)=TG(3,1)+EG13*DW(5)
          TG(2,1)=TG(2,1)+EG12*DW(5)
          TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
          TG(1,2)=TG(2,1)
          TG(1,3)=TG(3,1)
          TG(2,3)=TG(3,2)
        ENDIF

      ELSE IF(KTYP55(nr).EQ.2) THEN !TG is fn of the principal stretches
        DO i=1,3
          DO j=1,3
            EG(i,j)=0.5d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO !j
        ENDDO !i
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EG:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((EG(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IFAIL=0
        CALL F02ABF(EG,3,NITB,EVAL,RM,3,WK1_LOCAL,IFAIL)
C news MPN 30-Jun-94
C ***   Order EVAL from max->min and change cols of RM accordingly
        IF(EVAL(2).GT.EVAL(1)) THEN
          VALTMP=EVAL(1)
          EVAL(1)=EVAL(2)
          EVAL(2)=VALTMP
          DO ni=1,NITB
            VALTMP=RM(ni,1)
            RM(ni,1)=RM(ni,2)
            RM(ni,2)=VALTMP
          ENDDO
        ENDIF
        IF(NITB.EQ.3) THEN
          DO mi=1,2
            IF(EVAL(3).GT.EVAL(MI)) THEN
              VALTMP=EVAL(MI)
              EVAL(MI)=EVAL(3)
              EVAL(3)=VALTMP
              DO ni=1,NITB
                VALTMP=RM(ni,mi)
                RM(ni,mi)=RM(ni,3)
                RM(ni,3)=VALTMP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C newe
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EVAL:'',7X,(3X,3D11.3))')
     '      (EVAL(i),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RM:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((RM(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        RL1=DSQRT(2.0d0*EVAL(1)+1.0d0)
        RL2=DSQRT(2.0d0*EVAL(2)+1.0d0)
        IF(KTYP52(nr).EQ.1) THEN      !compressible
          RL3=1.0d0
        ELSE IF(KTYP52(nr).GT.1) THEN !incompressible
          RL3=1.0d0/(RL1*RL2)
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' RL1,RL2,RL3:'',(3X,3D11.3))')
     '      RL1,RL2,RL3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO i=1,3
          DO j=1,3
            TGP(i,j)=0.0d0
          ENDDO
        ENDDO
        CALL ENERGY(nr,CG,DW,RL1,RL2,RL3,0.0d0,0.0d0,0.0d0,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          TGP(1,1)=DW(1)/RL1
          TGP(2,2)=DW(2)/RL2
          TGP(3,3)=DW(3)/RL3
        ELSE IF(KTYP52(nr).GT.1)  THEN
          RL3=1.0d0/(RL1*RL2)
          TGP(1,1)=DW(1)/RL1
          TGP(2,2)=DW(2)/RL2
          TGP(3,3)=DW(3)/RL3
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TGP:'',12X,3D11.3,/(16X,3D11.3))')
     '       ((TGP(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL TRAN(3,TGP,TG,RM)

      ELSE IF(KTYP55(nr).EQ.3) THEN !TG is fn of the fibre trans strains
C!!! NOTE can only use resid strains for pole zero law until init extn
C!!!      mat params are set up for other problem types
        IF(KTYP56(nr).EQ.3) THEN !pole zero law
C         Calc growth defm tens for resid strain and copy AZL to AZL_tmp
          DO i=1,3
            DO j=1,3
              Fgrowth(i,j)=0.0d0 !growth defm tensor for resid strains
              AZL_tmp(i,j)=AZL(i,j)
            ENDDO !j
            Fgrowth(i,i)=CG(27+i) !NOTE init extns CG(28),CG(29),CG(30)
          ENDDO !i
C         Apply growth defm to deformed covariant metric tensor AZL
          DO i=1,3
            DO j=1,3
              AZL(i,j)=0.0d0
              DO m=1,3
                DO n=1,3
                  AZL(i,j)=AZL(i,j)+
     '              Fgrowth(i,m)*AZL_tmp(m,n)*Fgrowth(n,j)
                ENDDO !n
              ENDDO !m
            ENDDO !j
          ENDDO !i
C         Recompute AZU,AZ from transformed AZL
          CALL INVERT(NIT(nb),AZL,AZU,AZ)
        ENDIF !pole zero law
        DO i=1,3
          DO j=1,3
            EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO !j
        ENDDO !i
        IF(DOP) THEN
          WRITE(OP_STRING,'('' EG:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((EG(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '    EG(2,3),ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN !compressible
C ***     22Jan89: Compressible case should also be included here
        ELSE                 !incompressible
          H=-DW(3)/AZU(3,3)
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'('' DW(1..4):'',4D11.3,'' H='',D11.3)')
     '      (DW(i),i=1,4),H
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        TG(1,1)=DW(1)/AXL(1,1)                +H*AZU(1,1)
        TG(2,2)=DW(2)/AXL(2,2)                +H*AZU(2,2)
        TG(1,2)=DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+H*AZU(1,2)
        TG(2,1)=TG(1,2)
        TG(1,3)=0.0d0
        TG(2,3)=0.0d0
        TG(3,1)=0.0d0
        TG(3,2)=0.0d0
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),j=1,3),i=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ZGTG54')
      RETURN
 9999 CALL ERRORS('ZGTG54',ERROR)
      CALL EXITS('ZGTG54')
      RETURN 1
      END


      SUBROUTINE ZGTG55(nr,AZL,CG,EG,TG,ERROR,*)

C#### Subroutine: ZGTG55
C###  Description:
C###    ZGTG55 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for string problems
C###    (KTYP51(nr)=5).
C**** Stress components must always be with respect to
C**** (undeformed) Nu (body/fibre) coords (KTYP53(nr)>1).
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER nr
      REAL*8 AZL(3,3),CG(NMM),EG(3,3),TG(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 STRAIN

      CALL ENTERS('ZGTG55',*9999)

      IF(KTYP53(nr).EQ.1) THEN      !stresses in theta coords
        ERROR='>>Stresses must be wrt body/fibre coords'
        GO TO 9999
      ELSE IF(KTYP53(nr).GT.1) THEN !stresses in body/fibre Nu coords
        STRAIN=0.50d0*(AZL(1,1)-1.0d0)
        IF(DABS(STRAIN).LT.1.E+02) THEN
          EG(1,1)=STRAIN
        ELSE
          EG(1,1)=100.0d0
C          ERROR='>>>Strain cpt out of bounds'
C          GO TO 9999
        ENDIF
      ENDIF

      IF(KTYP55(nr).EQ.1) THEN !TG is fn of the principal invariants
        ERROR='>>Stresses must be fn of fibre strain'
        GO TO 9999
      ELSE IF(KTYP55(nr).EQ.2) THEN
        !TG is a fn of the principal stretches
        ERROR='>>Stresses must be fn of fibre strain'
        GO TO 9999
      ELSE IF(KTYP55(nr).EQ.3) THEN
        !TG is func of the fibre transv strains
        IF(KTYP52(nr).EQ.1) THEN !compressible
          IF(EG(1,1).GT.0.0d0) THEN !stretch
            TG(1,1)=CG(1)*EG(1,1)
          ELSE                     !compression
            TG(1,1)=0.0d0
          ENDIF
        ELSE                 !incompressible
          ERROR='>>Incompressible case not implemented'
          GO TO 9999
        ENDIF
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' AZL(1,1)='',D12.3,'' Strain=EG(1,1)='',D12.3,'
     '    //''' Stress=TG(1,1)='',D12.3)') AZL(1,1),EG(1,1),TG(1,1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('ZGTG55')
      RETURN
 9999 CALL ERRORS('ZGTG55',ERROR)
      CALL EXITS('ZGTG55')
      RETURN 1
      END


      SUBROUTINE ZGTG5A(nr,FEXT,DXNZN,DZNXN,TDELAY,TG,TNA,ERROR,*)

C#### Subroutine: ZGTG5A
C###  Description:
C###    ZGTG5A adds active fibre stress to TG array.
C**** FEXT(1,ng,ne) is current muscle fibre extension ratio
C****   "  2    "   "  previous   "     "       "       "
C****   "  3    "   "  muscle fibre ext. ratio at time of activation
C****   "  4    "   "  zero when gauss point is inactive else 1
C****   "  5    "   " previous hereditary integral for 1st time constant
C****   "  6    "   "     "         "         "     "  2nd  "      "
C****   "  7    "   "     "         "         "     "  3rd  "      "
C**** KTYP59(nr) is elastance/Hill-type/fading-memory formulation
C**** DEL_T  is the time step used (taken from load step loop in FE07)
C**** TV_SLO is slope of tension/vel. in lengthening (before yield)
C**** YIELDR is ratio of yield tension to isometric tension
C**** SNLPA  is static nonlinearity parameter "a"
C**** NTACTV is #dynamic terms in the material response function
C**** ACOEFF(nactv), nactv=1,NTACTV are coeffs   for lin dynamic terms
C**** ALFA(nactv),     "       "     "time constants  "     "      "
C****
C**** KTYP59(nr)=1 SS tension-length-ca relation
C**** ========
C**** FEXT(4,ng,ne) is current intracellular calcium concentration
C**** Tref is the max isometric tension at ext.ratio=1
C**** T0_beta is the non-dimensional slope parameter
C**** Ca_c50 is the c50 for the [Ca]i saturation curve (0<c<1)
C**** Ca_h is the Hill coefficient for the [Ca]i saturation curve
C**** Cai (0<Cai<1) is the Ca variable indicating degree of activation
C****
C**** KTYP59(nr)=3 Fading memory model
C**** ========
C**** PARAMS FOR THE EXTRACELLULAR CA CONC RELATION Cao(t)****
C**** time1  is the time (s)    at junction of 1st & 2nd cubic elements.
C**** Cao1   is the Ca conc. (mM)     "     "   "  "  "    "      "
C**** Caslo1 is the slope of the fn   "     "   "  "  "    "      "
C**** time2  is the time (s)    at junction of 2nd & 3rd cubic elements.
C**** Cao2   is the Ca conc. (mM)     "     "   "  "  "   "  "    "
C**** Caslo2 is the slope of the fn   "     "   "  "  "   "  "    "
C**** t_end  is the time (s) after onset of contraction,
C****                        when Cao drops to zero
C****
C**** PARAMETERS FOR THE ISOMETRIC TENSION RELATION, T0(Cao,FEXT) ****
C**** FEXTo  is the exten ratio at which there is zero active tension
C**** FEXTmx is the max permitted exten ratio (limited by expt'l data)
C**** T2_k   is the constant of the calcium dep nodal value @ XI2=1
C**** T2_a   is the pole     "   "     "     "    "     "   @ XI2=1
C**** dT1_k  is the constant of the calcium dep nodal slope @ XI2=0
C**** dT1_a  is the pole     "   "     "     "    "     "   @ XI2=0
C**** dT2_k  is the constant "   "     "     "    "     "   @ XI2=1
C**** dT2_a  is the pole     "   "     "     "    "     "   @ XI2=1
C****
C**** VARIABLES ****
C**** XI1    is the local [0,1] variable varying with time
C**** XI2    is the local [0,1] variable varying with FEXT
C**** Ca(i), i=1,4 are the time dep nodal params for the Cao function
C**** T(i),  i=1,4 are the Ca and length dep nodal params for the T0 fn

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:acti01.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:nonl00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
!     Parameter List
      INTEGER nr
      REAL*8 DXNZN(3,3),DZNXN(3,3),FEXT(NIFEXTM),
     '  TDELAY,TG(3,3),TN(3,3),TNA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,mix,mz,nactv,nb,nix,nz
      REAL*8 Ca(4),Cai,Cao,Caosqr,DFEXT,PSI10,PSI11,PSI20,PSI21,
     '  SUM,T(4),T0,TIMEAA,Q,VEL,XI,XI1,XI2

      REAL*8 Cao1,Cao2,Caslo1,Caslo2,dT1_a,dT1_k,dT2_k,dT2_a,
     '  FEXTo,FEXTmx,T2_a,T2_k,t_end,time1,time2
      DATA Cao1 /1.347d0/, Cao2 /0.6918d0/, Caslo1 /22.48d0/,
     '  Caslo2 /-1.916d0/,
     '  dT1_a /0.77d0/, dT1_k /449.d0/, dT2_k /31.d0/, dT2_a /0.13d0/,
     '  FEXTo /0.865d0/, FEXTmx /1.189d0/, T2_a /0.13d0/, T2_k /137.d0/,
     '  t_end /1.023d0/, time1 /0.1d0/, time2 /0.4173d0/

      PSI10(XI) = 1.0d0 - 3.0d0*XI*XI + 2.0d0*XI*XI*XI !Cubic H basis fn
      PSI11(XI) = XI*(XI-1.0d0)*(XI-1.0d0)
      PSI20(XI) = XI*XI*(3.0d0-2.0d0*XI)
      PSI21(XI) = XI*XI*(XI-1.0d0)

      CALL ENTERS('ZGTG5A',*9999)

      TNA=0.0d0
      IF(KTYP59(nr).EQ.1) THEN      !SS tension-length-Ca relation

        Cai = FEXT(4)
        IF(Cai.LT.0.0d0) THEN
          WRITE(OP_STRING,'('' Cannot handle negative Cai='','
     '      //'D12.3,''; continuing using zero active stress'')') Cai
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GOTO 9998
        ENDIF
        TNA = Ca_max*Cai**Ca_h/(Cai**Ca_h + Ca_c50**Ca_h)
     '        *Tref*(1.d0+T0_beta*(FEXT(1)-1.d0))
        IF(DOP) THEN
          WRITE(OP_STRING,'('' Cai='',D12.3,'' TNA='',D12.3)') Cai,TNA
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(KTYP59(nr).EQ.3) THEN !Fading memory model

C****   Determine the time after activation.
        TIMEAA=DEL_T*DBLE(NOSTEP)-TDELAY
        IF(DOP) THEN
          WRITE(OP_STRING,'('' TIMEAA='',D13.6)') TIMEAA
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(TIMEAA.GE.0.0d0)THEN
          IF(FEXT(4).EQ.0.0d0) THEN
            FEXT(3)=FEXT(1)
            FEXT(4)=1.0d0
          ENDIF
C****     Determine change in sarcomere length and current estimate of
C****     the velocity. Note: positive velocity for lengthening.
          DFEXT = (FEXT(1)-FEXT(2))
          VEL = DFEXT/DEL_T
          IF(DOP) THEN
            WRITE(OP_STRING,'('' FEXT(1..7)='',7D10.3/'' DFEXT='','
     '        //'D11.4,'' DFEXT='',D11.4)') (FEXT(i),I=1,7),DFEXT,VEL
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C****     Calculate Cao as a function of time after activation.
C****     Uses three cubic hermite elements.
          XI1 = 0.0d0
          DO nb=1,4
            Ca(nb) = 0.0d0
          ENDDO
          IF((TIMEAA.GE.0.0d0).AND.(TIMEAA.LT.time1))THEN
            XI1 = TIMEAA/time1
            Ca(2) = Cao1
            Ca(4) = Caslo1*time1
          ELSE IF((TIMEAA.GE.time1).AND.(TIMEAA.LT.time2))THEN
            XI1 = (TIMEAA-time1)/(time2-time1)
            Ca(1) = Cao1
            Ca(2) = Cao2
            Ca(3) = Caslo1*(time2-time1)
            Ca(4) = Caslo2*(time2-time1)
          ELSE IF((TIMEAA.GE.time2).AND.(TIMEAA.LT.t_end))THEN
            XI1 = (TIMEAA-time2)/(t_end-time2)
            Ca(1) = Cao2
            Ca(3) = Caslo2*(t_end-time2)
          ENDIF
          Cao = PSI10(XI1)*Ca(1) + PSI20(XI1)*Ca(2) + PSI11(XI1)*Ca(3) +
     '          PSI21(XI1)*Ca(4)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Cao='',D13.6)') Cao
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C****     T0 calculated as a fn of FEXT(1) and Cao.
C****     Uses 1 cubic hermite element.
          XI2 = 0.0d0
          DO nb=1,4
            T(nb) = 0.0d0
          ENDDO
          IF(FEXT(1).GE.FEXTo)THEN
            CALL ASSERT(FEXT(1).LE.FEXTmx,'>>Extension Ratio greater'
     '        //' than maximum FEXT ',ERROR,*9999)
            XI2 = (FEXT(1)-FEXTo)/(FEXTmx-FEXTo)  !Local Variable XI2.
            Caosqr=Cao*Cao
            T(1) = 0.0d0                          !Ca dep parameters
            T(2) = T2_k  * Caosqr/(Caosqr+T2_a )  !(Hill type relations).
            T(3) = dT1_k * Caosqr/(Caosqr+dT1_a)
            T(4) = dT2_k * Caosqr/(Caosqr+dT2_a)
          ENDIF
          T0 = PSI10(XI2)*T(1) + PSI20(XI2)*T(2) +PSI11(XI2)*T(3)
     '       +PSI21(XI2)*T(4)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' T0='',D13.6)') T0
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C****     Determine the magnitude of the active muscle fibre stress.
          IF(VEL.GE.0.0d0) THEN        !lengthing
            TNA = TV_SLO*VEL + T0
            IF(TNA.GT.(YIELDR*T0)) TNA = YIELDR*T0

          ELSE IF(VEL.LT.0.0d0) THEN   !shortening
C****       Calculate hereditary integral, summing over each rate constant
            Q=0.0d0
            DO nactv=1,NTACTV
              CALL ASSERT(ALFA(nactv)*DEL_T.LT.80.0d0, !to avoid underflow
     '          '>>Reduce time step using the DEFINE ACTIVE command',
     '          ERROR,*9999)
              Q = Q + ACOEFF(nactv)*(DEXP(-ALFA(nactv)*DEL_T)
     '          *FEXT(nactv+4)+DEXP(-ALFA(nactv)*DEL_T/2.0d0)*DFEXT)
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Q='',D13.6)') Q
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            TNA=T0*(1.0d0+SNLPA*Q)/(1.0d0-Q)
          ENDIF
        ENDIF !timeaa

      ENDIF !KTYP59(nr)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' TNA='',D13.6)') TNA
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C**** Rotate TG into TN whose 1 dirn lines up with the def fibre axis
C new MPN 5-May-96: rotating the coordinate system so must use
C                   tensor transformation (ie rotate both indices)
      DO mz=1,3
        DO nz=1,3
          SUM=0.0d0
          DO mix=1,3
            DO nix=1,3
              SUM=SUM+TG(mix,nix)*DZNXN(mz,mix)*DZNXN(nz,nix)
            ENDDO !nix
          ENDDO !mix
          TN(mz,nz)=SUM
        ENDDO !nz
      ENDDO !mz
C old
C      DO nix=1,3
C        DO nz=1,3
C          SUM=0.0d0
C          DO mix=1,3
C            SUM=SUM+TG(nix,mix)*DZNXN(nz,mix)
C          ENDDO
C          TN(nix,nz)=SUM
C        ENDDO
C      ENDDO

C**** Add active stress (in the 1,1 dirn as fibres generate the force)
      TN(1,1)=TN(1,1)+TNA

C**** Rotate TN back into TG with the original orientation
C new MPN 5-May-96: rotating the coordinate system so must use
C                   tensor transformation (ie rotate both indices)
      DO mix=1,3
        DO nix=1,3
          SUM=0.0d0
          DO mz=1,3
            DO nz=1,3
              SUM=SUM+TN(mz,nz)*DXNZN(mix,mz)*DXNZN(nix,nz)
            ENDDO !nix
          ENDDO !mix
          TG(mix,nix)=SUM
        ENDDO !nz
      ENDDO !mz
C old
C      DO nix=1,3
C        DO mix=1,3
C          SUM=0.0d0
C          DO mz=1,3
C            SUM=SUM+TN(nix,mz)*DXNZN(mix,mz)
C          ENDDO
C          TG(nix,mix)=SUM
C        ENDDO
C      ENDDO

      IF(DOP) THEN
        WRITE(OP_STRING,'('' TN(1,1)='',D13.6)') TN(1,1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),j=1,3),i=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

 9998 CALL EXITS('ZGTG5A')
      RETURN
 9999 CALL ERRORS('ZGTG5A',ERROR)
      CALL EXITS('ZGTG5A')
      RETURN 1
      END


      SUBROUTINE ZERE30(A,*)
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      DIMENSION A(*)
      CHARACTER ERROR*10
      WRITE(OP_STRING(1),*) '>>Link with module FE30'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
9999  RETURN 1
      END


      SUBROUTINE ZERE40(A,*)
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      DIMENSION A(*)
      CHARACTER ERROR*10
      WRITE(OP_STRING(1),*) '>>Link with module FE40'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
9999  RETURN 1
      END
