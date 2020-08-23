C $LARGE
C $NOFLOATCALLS
C                            APPENDIX 1                          1/15/81
C
C        SUBROUTINES FOR SOLVING SPARSE NONSYMMETRIC SYSTEMS
C        OF LINEAR EQUATIONS  (UNCOMPRESSED POINTER STORAGE)
C
C        REAL*8 VERSION. NOTE: THE ORIGINAL SUBROUTINES
C
C            NDRV, NSF, NNF, NNS AND NNT
C
C        HAVE BEEN RENAMED TO
C
C            YALE8, NSF8, NNF8, NNS8 AND NNT8, RESPECTIVELY.
C
C*** SUBROUTINE YALE8 (OLD NAME: NDRV)
C*** DRIVER FOR SUBROUTINES FOR SOLVING SPARSE NONSYMMETRIC SYSTEMS OF
C       LINEAR EQUATIONS (UNCOMPRESSED POINTER STORAGE)
C
C       SUBROUTINE  NDRV  (= OLD NAME)
C       SUBROUTINE  YALE8 (= NEW NAME)
        SUBROUTINE  NDRV
     *     (N, R,C,IC, IA,JA,A, B, Z, NSP,ISP,RSP,ESP, PATH, FLAG)
C
C    PARAMETERS
C    CLASS ABBREVIATIONS ARE --
C       N - INTEGER VARIABLE
C       F - REAL VARIABLE
C       V - SUPPLIES A VALUE TO THE DRIVER
C       R - RETURNS A RESULT FROM THE DRIVER
C       I - USED INTERNALLY BY THE DRIVER
C       A - ARRAY
C
C CLASS   PARAMETER
C ------+----------
C
C         THE NONZERO ENTRIES OF THE COEFFICIENT MATRIX M ARE STORED
C    ROW-BY-ROW IN THE ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO
C    ENTRIES IN EACH ROW, WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY
C    LIES.  THE COLUMN INDICES WHICH CORRESPOND TO THE NONZERO ENTRIES
C    OF M ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
C    JA(K) = J.  IN ADDITION, WE NEED TO KNOW WHERE EACH ROW STARTS AND
C    HOW LONG IT IS.  THE INDEX POSITIONS IN JA AND A WHERE THE ROWS OF
C    M BEGIN ARE STORED IN THE ARRAY IA;  I.E., IF M(I,J) IS THE FIRST
C    NONZERO ENTRY (STORED) IN THE I-TH ROW AND A(K) = M(I,J),  THEN
C    IA(I) = K.  MOREOVER, THE INDEX IN JA AND A OF THE FIRST LOCATION
C    FOLLOWING THE LAST ELEMENT IN THE LAST ROW IS STORED IN IA(N+1).
C    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS GIVEN BY
C    IA(I+1) - IA(I),  THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED
C    CONSECUTIVELY IN
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),
C    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).
C    FOR EXAMPLE, THE 5 BY 5 MATRIX
C                ( 1. 0. 2. 0. 0.)
C                ( 0. 3. 0. 0. 0.)
C            M = ( 0. 4. 5. 6. 0.)
C                ( 0. 0. 0. 7. 0.)
C                ( 0. 0. 0. 8. 9.)
C    WOULD BE STORED AS
C                 1  2  3  4  5  6  7  8  9
C            ---+--------------------------
C            IA   1  3  4  7  8 10
C            JA   1  3  2  2  3  4  4  4  5
C             A   1. 2. 3. 4. 5. 6. 7. 8. 9.         .
C
C NV      N     - NUMBER OF VARIABLES/EQUATIONS.
C FVA     A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED
C                   BY ROWS.
C                   SIZE = NUMBER OF NONZERO ENTRIES IN M.
C NVA     IA    - POINTERS TO DELIMIT THE ROWS IN A.
C                   SIZE = N+1.
C NVA     JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A.
C                   SIZE = SIZE OF A.
C FVA     B     - RIGHT-HAND SIDE B;  B AND Z CAN THE SAME ARRAY.
C                   SIZE = N.
C FRA     Z     - SOLUTION X;  B AND Z CAN BE THE SAME ARRAY.
C                   SIZE = N.
C
C         THE ROWS AND COLUMNS OF THE ORIGINAL MATRIX M CAN BE
C    REORDERED (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY)
C    BEFORE CALLING THE DRIVER.  IF NO REORDERING IS DONE, THEN SET
C    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE SOLUTION Z IS RETURNED
C    IN THE ORIGINAL ORDER.
C
C NVA     R     - ORDERING OF THE ROWS OF M.
C                   SIZE = N.
C NVA     C     - ORDERING OF THE COLUMNS OF M.
C                   SIZE = N.
C NVA     IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E.,
C                   IC(C(I)) = I  FOR I=1,...,N.
C                   SIZE = N.
C
C         THE SOLUTION OF THE SYSTEM OF LINEAR EQUATIONS IS DIVIDED INTO
C    THREE STAGES --
C      NSF -- THE MATRIX M IS PROCESSED SYMBOLICALLY TO DETERMINE WHERE
C              FILLIN WILL OCCUR DURING THE NUMERIC FACTORIZATION.
C      NNF -- THE MATRIX M IS FACTORED NUMERICALLY INTO THE PRODUCT LDU
C              OF A UNIT LOWER TRIANGULAR MATRIX L, A DIAGONAL MATRIX D,
C              AND A UNIT UPPER TRIANGULAR MATRIX U, AND THE SYSTEM
C              MX = B  IS SOLVED.
C      NNS -- THE LINEAR SYSTEM  MX = B  IS SOLVED USING THE LDU
C  OR          FACTORIZATION FROM NNF.
C      NNT -- THE TRANSPOSED LINEAR SYSTEM  MT X = B  IS SOLVED USING
C              THE LDU FACTORIZATION FROM NNF.
C    FOR SEVERAL SYSTEMS WHOSE COEFFICIENT MATRICES HAVE THE SAME
C    NONZERO STRUCTURE, NSF NEED BE DONE ONLY ONCE (FOR THE FIRST
C    SYSTEM);  THEN NNF IS DONE ONCE FOR EACH ADDITIONAL SYSTEM.  FOR
C    SEVERAL SYSTEMS WITH THE SAME COEFFICIENT MATRIX, NSF AND NNF NEED
C    BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN NNS OR NNT IS DONE
C    ONCE FOR EACH ADDITIONAL RIGHT-HAND SIDE.
C
C NV      PATH  - PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE --
C                   1  PERFORM NSF AND NNF.
C                   2  PERFORM NNF ONLY  (NSF IS ASSUMED TO HAVE BEEN
C                       DONE IN A MANNER COMPATIBLE WITH THE STORAGE
C                       ALLOCATION USED IN THE DRIVER).
C                   3  PERFORM NNS ONLY  (NSF AND NNF ARE ASSUMED TO
C                       HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE
C                       STORAGE ALLOCATION USED IN THE DRIVER).
C                   4  PERFORM NNT ONLY  (NSF AND NNF ARE ASSUMED TO
C                       HAVE BEEN DONE IN A MANNER COMPATIBLE WITH THE
C                       STORAGE ALLOCATION USED IN THE DRIVER).
C                   5  PERFORM NSF ONLY.
C
C         VARIOUS ERRORS ARE DETECTED BY THE DRIVER AND THE INDIVIDUAL
C    SUBROUTINES.
C
C NR      FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE --
C                     0     NO ERRORS DETECTED
C                     N+K   NULL ROW IN A  --  ROW = K
C                    2N+K   DUPLICATE ENTRY IN A  --  ROW = K
C                    3N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K
C                    4N+1   INSUFFICIENT STORAGE IN NNF
C                    5N+K   NULL PIVOT  --  ROW = K
C                    6N+K   INSUFFICIENT STORAGE IN NSF  --  ROW = K
C                    7N+1   INSUFFICIENT STORAGE IN NNF
C                    8N+K   ZERO PIVOT  --  ROW = K
C                   10N+1   INSUFFICIENT STORAGE IN NDRV
C                   11N+1   ILLEGAL PATH SPECIFICATION
C
C         WORKING STORAGE IS NEEDED FOR THE FACTORED FORM OF THE MATRIX
C    M PLUS VARIOUS TEMPORARY VECTORS.  THE ARRAYS ISP AND RSP SHOULD BE
C    EQUIVALENCED;  INTEGER STORAGE IS ALLOCATED FROM THE BEGINNING OF
C    ISP AND REAL STORAGE FROM THE END OF RSP.
C
C NV      NSP   - DECLARED DIMENSION OF RSP;  NSP GENERALLY MUST
C                   BE LARGER THAN  5N+3 + 2K  (WHERE  K = (NUMBER OF
C                   NONZERO ENTRIES IN M)).
C NVIRA   ISP   - INTEGER WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS
C                   NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE
C                   EQUIVALENCED.
C                   SIZE = LRATIO*NSP
C FVIRA   RSP   - REAL WORKING STORAGE DIVIDED UP INTO VARIOUS ARRAYS
C                   NEEDED BY THE SUBROUTINES;  ISP AND RSP SHOULD BE
C                   EQUIVALENCED.
C                   SIZE = NSP.
C NR      ESP   - IF SUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE
C                   SYMBOLIC FACTORIZATION (NSF), THEN ESP IS SET TO THE
C                   AMOUNT OF EXCESS STORAGE PROVIDED (NEGATIVE IF
C                   INSUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE
C                   NUMERIC FACTORIZATION (NNF)).
C
C
C  CONVERSION TO DOUBLE PRECISION
C
C    TO CONVERT THESE ROUTINES FOR DOUBLE PRECISION ARRAYS, SIMPLY USE
C    THE DOUBLE PRECISION DECLARATIONS IN PLACE OF THE REAL DECLARATIONS
C    IN EACH SUBPROGRAM;  IN ADDITION, THE DATA VALUE OF THE INTEGER
C    VARIABLE LRATIO MUST BE SET AS INDICATED IN SUBROUTINE NDRV
C
        INTEGER  R(1), C(1), IC(1),  IA(1), JA(1),  ISP(1), ESP,
     *     PATH, FLAG,  Q, IM, D, U, ROW, TMP,  UMAX
C       REAL  A(1),  B(1),  Z(1),  RSP(1)
        DOUBLE PRECISION  A(1),  B(1),  Z(1),  RSP(1)
C
C  SET LRATIO EQUAL TO THE RATIO BETWEEN THE LENGTH OF FLOATING POINT
C  AND INTEGER ARRAY DATA;  E. G., LRATIO = 1 FOR (REAL, INTEGER),
C  LRATIO = 2 FOR (DOUBLE PRECISION, INTEGER)
C
        DATA LRATIO/2/
C
        IF (PATH.LT.1 .OR. 5.LT.PATH)  GO TO 111
C  ******  INITIALIZE AND DIVIDE UP TEMPORARY STORAGE  *****************
        IL = 1
        IU = IL + N+1
        JL = IU + N+1
C
C  ******  CALL NSF IF FLAG IS SET  ************************************
        IF ((PATH-1) * (PATH-5) .NE. 0)  GO TO 2
          MAX = (LRATIO*NSP + 1 - JL) - (N+1) - N
          JLMAX = MAX/2
          Q     = JL  + JLMAX
          IM    = Q   + (N+1)
          JUTMP = IM  +   N
          JUMAX = LRATIO*NSP + 1 - JUTMP
          ESP = MAX/LRATIO
          IF (JLMAX.LE.0 .OR. JUMAX.LE.0)  GO TO 110
          CALL  NSF8
     *       (N,  R, IC,  IA, JA,
     *        ISP(IL), ISP(JL), JLMAX,  ISP(IU), ISP(JUTMP), JUMAX,
     *        ISP(Q),  ISP(IM),  FLAG)
          IF (FLAG.NE.0)  GO TO 100
C  ******  MOVE JU NEXT TO JL  *****************************************
          JLMAX = ISP(IL+N)-1
          JU    = JL + JLMAX
          JUMAX = ISP(IU+N)-1
          IF (JUMAX.LE.0)  GO TO 2
          DO 1 J=1,JUMAX
   1        ISP(JU+J-1) = ISP(JUTMP+J-1)
C
C  ******  CALL REMAINING SUBROUTINES  *********************************
   2    JLMAX = ISP(IL+N)-1
        JU    = JL  + JLMAX
        JUMAX = ISP(IU+N)-1
        L     = (JU + JUMAX - 2 + LRATIO)  /  LRATIO    +    1
        LMAX  = JLMAX
        D     = L   + LMAX
        U     = D   + N
        ROW   = NSP + 1 - N
        TMP   = ROW - N
        UMAX  = TMP - U
        ESP = UMAX - JUMAX
C
        IF ((PATH-1) * (PATH-2) .NE. 0)  GO TO 3
          IF (UMAX.LE.0)  GO TO 110
          CALL  NNF8
     *       (N,  R, C, IC,  IA, JA, A,  Z,  B,
     *        ISP(IL), ISP(JL), RSP(L), LMAX,   RSP(D),
     *           ISP(IU), ISP(JU), RSP(U), UMAX,
     *        RSP(ROW),  RSP(TMP),  FLAG)
          IF (FLAG.NE.0)  GO TO 100
          RETURN
C
   3    IF ((PATH-3) .NE. 0)  GO TO 4
          CALL  NNS8
     *       (N,  R, C,
     *        ISP(IL), ISP(JL), RSP(L),  RSP(D),
     *          ISP(IU), ISP(JU), RSP(U),
     *        Z,  B,  RSP(TMP))
C
   4    IF ((PATH-4) .NE. 0)  GO TO 5
          CALL  NNT8
     *       (N,  R, C,
     *        ISP(IL), ISP(JL), RSP(L),  RSP(D),
     *          ISP(IU), ISP(JU), RSP(U),
     *        Z,  B,  RSP(TMP))
   5    RETURN
C
C ** ERROR:  ERROR DETECTED IN NSF, NNF, NNS, OR NNT
 100    RETURN
C ** ERROR:  INSUFFICIENT STORAGE
 110    FLAG = 10*N + 1
        RETURN
C ** ERROR:  ILLEGAL PATH SPECIFICATION
 111    FLAG = 11*N + 1
        RETURN
        END
C
C       ----------------------------------------------------------------
C
C               YALE SPARSE MATRIX PACKAGE - NONSYMMETRIC CODES
C                    SOLVING THE SYSTEM OF EQUATIONS MX = B
C                        (UNCOMPRESSED POINTER STORAGE)
C
C    I.   CALLING SEQUENCES
C         THE COEFFICIENT MATRIX CAN BE PROCESSED BY AN ORDERING ROUTINE
C    (E.G., TO REDUCE FILLIN OR ENSURE NUMERICAL STABILITY) BEFORE USING
C    THE REMAINING SUBROUTINES.  IF NO REORDERING IS DONE, THEN SET
C    R(I) = C(I) = IC(I) = I  FOR I=1,...,N.  THE CALLING SEQUENCE IS --
C        (      (MATRIX ORDERING))
C         NSF   (SYMBOLIC FACTORIZATION TO DETERMINE WHERE FILLIN WILL
C                 OCCUR DURING NUMERIC FACTORIZATION)
C         NNF   (NUMERIC FACTORIZATION INTO PRODUCT LDU OF UNIT LOWER
C                 TRIANGULAR MATRIX L, DIAGONAL MATRIX D, AND UNIT UPPER
C                 TRIANGULAR MATRIX U, AND SOLUTION OF LINEAR SYSTEM)
C         NNS   (SOLUTION OF LINEAR SYSTEM FOR ADDITIONAL RIGHT-HAND
C     OR          SIDE USING LDU FACTORIZATION FROM NNF)
C         NNT   (SOLUTION OF TRANSPOSED LINEAR SYSTEM FOR ADDITIONAL
C                 RIGHT-HAND SIDE USING LDU FACTORIZATION FROM NNF)
C
C    II.  STORAGE OF SPARSE MATRICES
C         THE NONZERO ENTRIES OF THE COEFFICIENT MATRIX M ARE STORED
C    ROW-BY-ROW IN THE ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO
C    ENTRIES IN EACH ROW, WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY
C    LIES.  THE COLUMN INDICES WHICH CORRESPOND TO THE NONZERO ENTRIES
C    OF M ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
C    JA(K) = J.  IN ADDITION, WE NEED TO KNOW WHERE EACH ROW STARTS AND
C    HOW LONG IT IS.  THE INDEX POSITIONS IN JA AND A WHERE THE ROWS OF
C    M BEGIN ARE STORED IN THE ARRAY IA;  I.E., IF M(I,J) IS THE FIRST
C    NONZERO ENTRY (STORED) IN THE I-TH ROW AND A(K) = M(I,J),  THEN
C    IA(I) = K.  MOREOVER, THE INDEX IN JA AND A OF THE FIRST LOCATION
C    FOLLOWING THE LAST ELEMENT IN THE LAST ROW IS STORED IN IA(N+1).
C    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS GIVEN BY
C    IA(I+1) - IA(I),  THE NONZERO ENTRIES OF THE I-TH ROW ARE STORED
C    CONSECUTIVELY IN
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),
C    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).
C    FOR EXAMPLE, THE 5 BY 5 MATRIX
C                ( 1. 0. 2. 0. 0.)
C                ( 0. 3. 0. 0. 0.)
C            M = ( 0. 4. 5. 6. 0.)
C                ( 0. 0. 0. 7. 0.)
C                ( 0. 0. 0. 8. 9.)
C    WOULD BE STORED AS
C                 1  2  3  4  5  6  7  8  9
C            ---+--------------------------
C            IA   1  3  4  7  8 10
C            JA   1  3  2  2  3  4  4  4  5
C             A   1. 2. 3. 4. 5. 6. 7. 8. 9.         .
C
C         THE STRICT TRIANGULAR PORTIONS OF THE MATRICES L AND U ARE
C    STORED IN THE SAME FASHION USING THE ARRAYS  IL, JL, L  AND
C    IU, JU, U  RESPECTIVELY.  THE DIAGONAL ENTRIES OF L AND U ARE
C    ASSUMED TO BE EQUAL TO ONE AND ARE NOT STORED.  THE ARRAY D
C    CONTAINS THE RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D.
C
C    III. ADDITIONAL STORAGE SAVINGS
C         IN NSF, R AND IC CAN BE THE SAME ARRAY IN THE CALLING
C    SEQUENCE IF NO REORDERING OF THE COEFFICIENT MATRIX HAS BEEN DONE.
C         IN NNF, R, C AND IC CAN ALL BE THE SAME ARRAY IF NO REORDERING
C    HAS BEEN DONE.  IF ONLY THE ROWS HAVE BEEN REORDERED, THEN C AND IC
C    CAN BE THE SAME ARRAY.  IF THE ROW AND COLUMN ORDERINGS ARE THE
C    SAME, THEN R AND C CAN BE THE SAME ARRAY.  Z AND ROW CAN BE THE
C    SAME ARRAY.
C         IN NNS OR NNT, R AND C CAN BE THE SAME ARRAY IF NO REORDERING
C    HAS BEEN DONE OR IF THE ROW AND COLUMN ORDERINGS ARE THE SAME.  Z
C    AND B CAN BE THE SAME ARRAY;  HOWEVER, THEN B WILL BE DESTROYED.
C
C    IV.  PARAMETERS
C         FOLLOWING IS A LIST OF PARAMETERS TO THE PROGRAMS.  NAMES ARE
C    UNIFORM AMONG THE VARIOUS SUBROUTINES.  CLASS ABBREVIATIONS ARE --
C       N - INTEGER VARIABLE
C       F - REAL VARIABLE
C       V - SUPPLIES A VALUE TO A SUBROUTINE
C       R - RETURNS A RESULT FROM A SUBROUTINE
C       I - USED INTERNALLY BY A SUBROUTINE
C       A - ARRAY
C
C CLASS   PARAMETER
C ------+----------
C FVA     A     - NONZERO ENTRIES OF THE COEFFICIENT MATRIX M, STORED
C                   BY ROWS.
C                   SIZE = NUMBER OF NONZERO ENTRIES IN M.
C FVA     B     - RIGHT-HAND SIDE B.
C                   SIZE = N.
C NVA     C     - ORDERING OF THE COLUMNS OF M.
C                   SIZE = N.
C FVRA    D     - RECIPROCALS OF THE DIAGONAL ENTRIES OF THE MATRIX D.
C                   SIZE = N.
C NR      FLAG  - ERROR FLAG;  VALUES AND THEIR MEANINGS ARE --
C                    0     NO ERRORS DETECTED
C                    N+K   NULL ROW IN A  --  ROW = K
C                   2N+K   DUPLICATE ENTRY IN A  --  ROW = K
C                   3N+K   INSUFFICIENT STORAGE FOR JL  --  ROW = K
C                   4N+1   INSUFFICIENT STORAGE FOR L
C                   5N+K   NULL PIVOT  --  ROW = K
C                   6N+K   INSUFFICIENT STORAGE FOR JU  --  ROW = K
C                   7N+1   INSUFFICIENT STORAGE FOR U
C                   8N+K   ZERO PIVOT  --  ROW = K
C NVA     IA    - POINTERS TO DELIMIT THE ROWS IN A.
C                   SIZE = N+1.
C NVA     IC    - INVERSE OF THE ORDERING OF THE COLUMNS OF M;  I.E.,
C                   IC(C(I) = I  FOR I=1,...N.
C                   SIZE = N.
C NVRA    IL    - POINTERS TO DELIMIT THE ROWS IN L.
C                   SIZE = N+1.
C NVRA    IU    - POINTERS TO DELIMIT THE ROWS IN U.
C                   SIZE = N+1.
C NVA     JA    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF A.
C                   SIZE = SIZE OF A.
C NVRA    JL    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF L.
C                   SIZE = JLMAX.
C NV      JLMAX - DECLARED DIMENSION OF JL;  JLMAX MUST BE LARGER THAN
C                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER
C                   TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF).
C NVRA    JU    - COLUMN NUMBERS CORRESPONDING TO THE ELEMENTS OF U.
C                   SIZE = JUMAX.
C NV      JUMAX - DECLARED DIMENSION OF JU;  JUMAX MUST BE LARGER THAN
C                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER
C                   TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF).
C FVRA    L     - NONZERO ENTRIES IN THE STRICT LOWER TRIANGULAR PORTION
C                   OF THE MATRIX L, STORED BY ROWS.
C                   SIZE = LMAX
C NV      LMAX  - DECLARED DIMENSION OF L;  LMAX MUST BE LARGER THAN
C                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT LOWER
C                   TRIANGLE OF M PLUS FILLIN  (IL(N+1)-1 AFTER NSF).
C NV      N     - NUMBER OF VARIABLES/EQUATIONS.
C NVA     R     - ORDERING OF THE ROWS OF M.
C                   SIZE = N.
C FVRA    U     - NONZERO ENTRIES IN THE STRICT UPPER TRIANGULAR PORTION
C                   OF THE MATRIX U, STORED BY ROWS.
C                   SIZE = UMAX.
C NV      UMAX  - DECLARED DIMENSION OF U;  UMAX MUST BE LARGER THAN
C                   THE NUMBER OF NONZERO ENTRIES IN THE STRICT UPPER
C                   TRIANGLE OF M PLUS FILLIN  (IU(N+1)-1 AFTER NSF).
C FRA     Z     - SOLUTION X.
C                   SIZE = N.
C
C
C       ----------------------------------------------------------------
C
C*** SUBROUTINE NSF
C*** SYMBOLIC LDU-FACTORIZATION OF A NONSYMMETRIC SPARSE MATRIX
C      (UNCOMPRESSED POINTER STORAGE)
C
        SUBROUTINE  NSF8
     *     (N, R,IC, IA,JA, IL,JL,JLMAX, IU,JU,JUMAX, Q, IM, FLAG)
C
C       INPUT VARIABLES:   N, R,IC, IA,JA, JLMAX, JUMAX.
C       OUTPUT VARIABLES:  IL,JL, IU,JU, FLAG.
C
C       PARAMETERS USED INTERNALLY:
C NIA     Q     - SUPPOSE M' IS THE RESULT OF REORDERING M;  IF
C                   PROCESSING OF THE KTH ROW OF M' (HENCE THE KTH ROWS
C                   OF L AND U) IS BEING DONE, THEN Q(J) IS INITIALLY
C                   NONZERO IF M'(K,J) IS NONZERO;  SINCE VALUES NEED
C                   NOT BE STORED, EACH ENTRY POINTS TO THE NEXT
C                   NONZERO;  FOR EXAMPLE, IF  N=9  AND THE 5TH ROW OF
C                   M' IS
C                           0 X X 0 X 0 0 X 0,
C                   THEN Q WILL INITIALLY BE
C                           A 3 5 A 8 A A 10 A 2        (A - ARBITRARY);
C                   Q(N+1) POINTS TO THE FIRST NONZERO IN THE ROW AND
C                   THE LAST NONZERO POINTS TO  N+1;  AS THE ALGORITHM
C                   PROCEEDS, OTHER ELEMENTS OF Q ARE INSERTED IN THE
C                   LIST BECAUSE OF FILLIN.
C                   SIZE = N+1.
C NIA     IM    - AT EACH STEP IN THE FACTORIZATION, IM(I) IS THE LAST
C                   ELEMENT IN THE ITH ROW OF U WHICH NEEDS TO BE
C                   CONSIDERED IN COMPUTING FILLIN.
C                   SIZE = N.
C
C  INTERNAL VARIABLES--
C    JLPTR - POINTS TO THE LAST POSITION USED IN  JL.
C    JUPTR - POINTS TO THE LAST POSITION USED IN  JU.
C
        INTEGER  R(1), IC(1),  IA(1), JA(1),  IL(1), JL(1),
     *     IU(1), JU(1),  Q(1),  IM(1),  FLAG,  QM, VJ
C
C  ******  INITIALIZE POINTERS  ****************************************
        JLPTR = 0
        IL(1) = 1
        JUPTR = 0
        IU(1) = 1
C
C  ******  FOR EACH ROW OF L AND U  ************************************
        DO 10 K=1,N
C  ******  SET Q TO THE REORDERED ROW OF A  ****************************
          Q(N+1) = N+1
          JMIN = IA(R(K))
          JMAX = IA(R(K)+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 101
          DO 2 J=JMIN,JMAX
            VJ = IC(JA(J))
            QM = N+1
   1        M = QM
            QM = Q(M)
            IF (QM.LT.VJ)  GO TO 1
            IF (QM.EQ.VJ)  GO TO 102
              Q(M) = VJ
              Q(VJ) = QM
   2        CONTINUE
C
C  ******  FOR EACH ENTRY IN THE LOWER TRIANGLE  ***********************
          I = N+1
   3      I = Q(I)
          IF (I.GE.K)  GO TO 7
C  ******  L(K,I) WILL BE NONZERO, SO ADD IT TO JL  ********************
            JLPTR = JLPTR+1
            IF (JLPTR.GT.JLMAX)  GO TO 103
            JL(JLPTR) = I
            QM = I
C  ******  INSPECT ITH ROW FOR FILLIN, ADJUST IM IF POSSIBLE  **********
            JMIN = IU(I)
            JMAX = IM(I)
            IF (JMIN.GT.JMAX)  GO TO 6
            DO 5 J=JMIN,JMAX
              VJ = JU(J)
              IF (VJ.EQ.K)  IM(I) = J
   4          M = QM
              QM = Q(M)
              IF (QM.LT.VJ)  GO TO 4
              IF (QM.EQ.VJ)  GO TO 5
                Q(M) = VJ
                Q(VJ) = QM
                QM = VJ
   5          CONTINUE
   6        GO TO 3
C
C  ******  CHECK FOR NULL PIVOT  ***************************************
   7      IF (I.NE.K)  GO TO 105
C  ******  REMAINING ELEMENTS OF Q DEFINE STRUCTURE OF U(K, )  *********
   8      I = Q(I)
          IF (I.GT.N)  GO TO 9
            JUPTR = JUPTR+1
            IF (JUPTR.GT.JUMAX)  GO TO 106
            JU(JUPTR) = I
            GO TO 8
C  ******  GET READY FOR NEXT ROW  *************************************
   9      IM(K) = JUPTR
          IL(K+1) = JLPTR+1
  10      IU(K+1) = JUPTR+1
C
        FLAG = 0
        RETURN
C
C ** ERROR:  NULL ROW IN A
 101    FLAG = N + R(K)
        RETURN
C ** ERROR:  DUPLICATE ENTRY IN A
 102    FLAG = 2*N + R(K)
        RETURN
C ** ERROR:  INSUFFICIENT STORAGE FOR JL
 103    FLAG = 3*N + K
        RETURN
C ** ERROR:  NULL PIVOT
 105    FLAG = 5*N + K
        RETURN
C ** ERROR:  INSUFFICIENT STORAGE FOR JU
 106    FLAG = 6*N + K
        RETURN
        END
C
C       ----------------------------------------------------------------
C
C*** SUBROUTINE NNF
C*** NUMERIC LDU-FACTORIZATION OF SPARSE NONSYMMETRIC MATRIX AND
C      SOLUTION OF SYSTEM OF LINEAR EQUATIONS (UNCOMPRESSED POINTER
C      STORAGE)
C
        SUBROUTINE  NNF8
     *     (N, R,C,IC, IA,JA,A, Z, B, IL,JL,L,LMAX, D, IU,JU,U,UMAX,
     *      ROW, TMP, FLAG)
C
C       INPUT VARIABLES:   N, R,C,IC, IA,JA,A, B, IL,JL,LMAX, IU,JU,UMAX
C       OUTPUT VARIABLES:  Z, L,D,U, FLAG
C
C       PARAMETERS USED INTERNALLY:
C FIA     ROW   - HOLDS INTERMEDIATE VALUES IN CALCULATION OF L, D, U.
C                   SIZE = N.
C FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE
C                   EQUATION  UX = B'.
C                   SIZE = N.
C
        INTEGER  R(1), C(1), IC(1),  IA(1), JA(1),
     *     IL(1), JL(1), LMAX,  IU(1), JU(1), UMAX,  FLAG
C       REAL  A(1), Z(1), B(1),  L(1), D(1), U(1),
C    *     ROW(1), TMP(1),  LI, SUM, DK
        DOUBLE PRECISION  A(1), Z(1), B(1),  L(1), D(1), U(1),
     *     ROW(1), TMP(1),  LI, SUM, DK
C
C  ******  CHECK STORAGE  **********************************************
        IF (IL(N+1)-1 .GT. LMAX)  GO TO 104
        IF (IU(N+1)-1 .GT. UMAX)  GO TO 107
C
C  ******  FOR EACH ROW  ***********************************************
        DO 10 K=1,N
C  ******  SET THE INITIAL STRUCTURE OF ROW  ***************************
          JMIN = IL(K)
          JMAX = IL(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 2
C  ******  IF L(K,M) .NE. 0, ROW(M)=0  *********************************
          DO 1 J=JMIN,JMAX
   1        ROW(JL(J)) = 0
   2      ROW(K) = 0
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 4
C  ******  IF U(K,M) .NE. 0, ROW(M)=0  *********************************
          DO 3 J=JMIN,JMAX
   3        ROW(JU(J)) = 0
   4      JMIN = IA(R(K))
          JMAX = IA(R(K)+1) - 1
C  ******  SET ROW TO KTH ROW OF REORDERED A  **************************
          DO 5 J=JMIN,JMAX
   5        ROW(IC(JA(J))) = A(J)
C  ******  INITIALIZE SUM  *********************************************
          SUM = B(R(K))
C
C  ******  ASSIGN THE KTH ROW OF L AND ADJUST ROW, SUM  ****************
          IMIN = IL(K)
          IMAX = IL(K+1) - 1
          IF (IMIN.GT.IMAX)  GO TO 8
          DO 7 I=IMIN,IMAX
            LI = - ROW(JL(I))
C  ******  IF L IS NOT REQUIRED, THEN COMMENT OUT THE FOLLOWING LINE  **
            L(I) = - LI
            SUM = SUM + LI * TMP(JL(I))
            JMIN = IU(JL(I))
            JMAX = IU(JL(I)+1) - 1
            IF (JMIN.GT.JMAX)  GO TO 7
            DO 6 J=JMIN,JMAX
   6          ROW(JU(J)) = ROW(JU(J)) + LI * U(J)
   7        CONTINUE
C
C  ******  ASSIGN DIAGONAL D AND KTH ROW OF U, SET TMP(K)  *************
   8      IF (ROW(K).EQ.0)  GO TO 108
          DK = 1 / ROW(K)
          D(K) = DK
          TMP(K) = SUM * DK
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 10
          DO 9 J=JMIN,JMAX
   9        U(J) = ROW(JU(J)) * DK
  10      CONTINUE
C
C  ******  SOLVE  UX = TMP  BY BACK SUBSTITUTION  **********************
        K = N
        DO 13 I=1,N
          SUM = TMP(K)
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 12
          DO 11 J=JMIN,JMAX
  11        SUM = SUM - U(J) * TMP(JU(J))
  12      TMP(K) = SUM
          Z(C(K)) = SUM
  13      K = K-1
C
        FLAG = 0
        RETURN
C
C ** ERROR:  INSUFFICIENT STORAGE FOR L
 104    FLAG = 4*N + 1
        RETURN
C ** ERROR:  INSUFFICIENT STORAGE FOR U
 107    FLAG = 7*N + 1
        RETURN
C ** ERROR:  ZERO PIVOT
 108    FLAG = 8*N + K
        RETURN
        END
C
C       ----------------------------------------------------------------
C
C*** SUBROUTINE NNS
C*** NUMERIC SOLUTION OF A SPARSE NONSYMMETRIC SYSTEM OF LINEAR
C      EQUATIONS GIVEN LDU-FACTORIZATION (UNCOMPRESSED POINTER STORAGE)
C
        SUBROUTINE  NNS8
     *     (N, R,C, IL,JL,L, D, IU,JU,U, Z, B, TMP)
C
C       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B
C       OUTPUT VARIABLES:  Z
C
C       PARAMETERS USED INTERNALLY:
C FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE
C                   EQUATION UX = B'.
C                   SIZE = N.
C
        INTEGER  R(1), C(1),  IL(1), JL(1),  IU(1), JU(1)
C       REAL  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), SUM
        DOUBLE PRECISION  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), SUM
C
C  ******  SOLVE LDY = B  BY FORWARD SUBSTITUTION  *********************
        DO 2 K=1,N
          SUM = B(R(K))
          JMIN = IL(K)
          JMAX = IL(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 2
          DO 1 J=JMIN,JMAX
   1        SUM = SUM - L(J) * TMP(JL(J))
   2      TMP(K) = SUM * D(K)
C
C  ******  SOLVE  UX = Y  BY BACK SUBSTITUTION  ************************
        K = N
        DO 5 I=1,N
          SUM = TMP(K)
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 4
          DO 3 J=JMIN,JMAX
   3        SUM = SUM - U(J) * TMP(JU(J))
   4      TMP(K) = SUM
          Z(C(K)) = SUM
   5      K = K-1
        RETURN
        END
C       ----------------------------------------------------------------
C
C*** SUBROUTINE NNT
C*** NUMERIC SOLUTION OF THE TRANSPOSE OF A SPARSE NONSYMMETRIC SYSTEM
C      OF LINEAR EQUATIONS GIVEN LDU-FACTORIZATION (UNCOMPRESSED POINTER
C      STORAGE)
C
        SUBROUTINE  NNT8
     *     (N, R,C, IL,JL,L, D, IU,JU,U, Z, B, TMP)
C
C       INPUT VARIABLES:   N, R,C, IL,JL,L, D, IU,JU,U, B
C       OUTPUT VARIABLES:  Z
C
C       PARAMETERS USED INTERNALLY:
C FIA     TMP   - HOLDS NEW RIGHT-HAND SIDE B' FOR SOLUTION OF THE
C                   EQUATION LX = B'.
C                   SIZE = N.
C
        INTEGER  R(1), C(1),  IL(1), JL(1),  IU(1), JU(1)
C       REAL  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), TMPK
        DOUBLE PRECISION  L(1), D(1), U(1),  Z(1), B(1),  TMP(1), TMPK
C
C  ******  SOLVE  UT Y = B  BY FORWARD SUBSTITUTION  *******************
        DO 1 K=1,N
   1      TMP(K) = B(C(K))
        DO 3 K=1,N
          TMPK = - TMP(K)
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 3
          DO 2 J=JMIN,JMAX
   2        TMP(JU(J)) = TMP(JU(J)) + U(J) * TMPK
   3      CONTINUE
C
C  ******  SOLVE  D LT X = Y  BY BACK SUBSTITUTION  ********************
        K = N
        DO 6 I=1,N
          TMPK = - (TMP(K) * D(K))
          JMIN = IL(K)
          JMAX = IL(K+1) - 1
          IF (JMIN.GT.JMAX)  GO TO 5
          DO 4 J=JMIN,JMAX
   4        TMP(JL(J)) = TMP(JL(J)) + L(J) * TMPK
   5      Z(R(K)) = - TMPK
   6      K = K-1
        RETURN
        END
