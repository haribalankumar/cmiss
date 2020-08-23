      SUBROUTINE CALCULATE_BASIS_STRING(IBT,nb,DATAFILE,BASIS_STRING,
     '  SPECIAL_BASIS_FLAG,ERROR,*)

C#### Subroutine: CALCULATE_BASIS_STRING
C###  Description:
C###    CALCULATE_BASIS_STRING calculates the BASIS_STRING and sets
C###    the SPECIAL_BASIS_FLAG for the specified basis (nb) when
C###    DATAFILE is true. When DATAFILE is false, the description is
C###    written to the socket and the SPECIAL_BASIS_FLAG is set for
C###    the specifed basis (nb)

      IMPLICIT NONE
      INCLUDE 'cmis00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
!     Parameter List
      CHARACTER BASIS_STRING*(*),ERROR*(*)
      INTEGER IBT(3,NIM,NBFM),nb,SPECIAL_BASIS_FLAG
      LOGICAL DATAFILE
!     SPECIAL_BASIS_FLAG is zero except where the basis is not
!     a tensor product basis
C       1 is Hermite simplex with apex at node 1
C       2 is Hermite simplex with apex at node 3
C       3 is sector element
C       4 other bases with cross derivs set to zero
!     Local Variables
      CHARACTER BASE(9)*15,CHAR1*1,SIMPLEX_BASIS_LIST*6
      INTEGER i,IBASE(9),IBEG,IEND,IBEG2,IEND2,INDEX_ni,ISOCKET(2),ni,
     '  SIMPLEX_BASIS_DIRNS(0:3)

      DATA BASE/
     '  'l.Lagrange     ',
     '  'q.Lagrange     ',
     '  'c.Lagrange     ',
     '  'c.Hermite      ',
     '  'LagrangeHermite',
     '  'HermiteLagrange',
     '  'l.simplex      ',
     '  'q.simplex      ',
     '  'c.simplex      '/
      DATA IBASE/1,2,3,4,5,6,7,8,9/

      CALL ENTERS('CALCULATE_BASIS_STRING',*9999)

      SPECIAL_BASIS_FLAG=0
C New Simplex elements CS 12/11/97
      SIMPLEX_BASIS_DIRNS(0)=0
      SIMPLEX_BASIS_DIRNS(1)=0
      DO ni=1,NIT(nb)
        IF(IBT(1,ni,nb).EQ.3) THEN
          SIMPLEX_BASIS_DIRNS(0)=SIMPLEX_BASIS_DIRNS(0)+1
          SIMPLEX_BASIS_DIRNS(SIMPLEX_BASIS_DIRNS(0))=ni
        ENDIF
      ENDDO
      IF(SIMPLEX_BASIS_DIRNS(0).NE.0) THEN
        SIMPLEX_BASIS_LIST='('
        DO i=1,SIMPLEX_BASIS_DIRNS(0)
          WRITE(CHAR1,'(I1)') SIMPLEX_BASIS_DIRNS(i+1)
          CALL STRING_TRIM(SIMPLEX_BASIS_LIST,IBEG,IEND)
          SIMPLEX_BASIS_LIST=SIMPLEX_BASIS_LIST(IBEG:IEND)//CHAR1//';'
        ENDDO
        CALL STRING_TRIM(SIMPLEX_BASIS_LIST,IBEG2,IEND2)
        SIMPLEX_BASIS_LIST=SIMPLEX_BASIS_LIST(IBEG:IEND-1)//')'
      ENDIF
      ni=1
      DO WHILE((ni.LE.NIT(nb)).AND.((SPECIAL_BASIS_FLAG.EQ.0).OR.
     '  (SPECIAL_BASIS_FLAG.EQ.3)))
        IF(IBT(1,ni,nb).EQ.1) THEN !Lagrange
          INDEX_ni=IBT(2,ni,nb) !linear,quad or cubic
        ELSE IF(IBT(1,ni,nb).EQ.2) THEN !cubic Hermite
c cpb 3/8/95 Adding quadratic Hermites
          IF(IBT(2,ni,nb).EQ.1) THEN
            INDEX_ni=4
          ELSE IF(IBT(2,ni,nb).EQ.2) THEN
            INDEX_ni=5
          ELSE
            INDEX_ni=6
          ENDIF
        ELSE IF(IBT(1,ni,nb).EQ.3) THEN !simplex
          IF(IBT(2,ni,nb).LT.4) THEN ! l.q.or c.
            INDEX_ni=IBT(2,ni,nb)+6 !linear,quad or cubic
          ELSE IF(IBT(2,ni,nb).EQ.4) THEN !Hermite
C Note these are hermite simplexes and will be dusted in the future.
C These are not what is currently refered to as simplex elements,
C they belong to the superset called sectors but haven't been changed
C yet. CS 12/11/97
            IF(NKT(1,nb).EQ.1) THEN !apex node 1
              IF(NIT(nb).EQ.2) THEN
                SPECIAL_BASIS_FLAG=1
                IF(DATAFILE) THEN
                  BASIS_STRING='c.Hermite*LagrangeHermite'
                ELSE
                  ISOCKET(1)=4
                  ISOCKET(2)=5
                  IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
              ELSE
                SPECIAL_BASIS_FLAG=-1
              ENDIF
            ELSE IF(NKT(3,nb).EQ.1) THEN !apex node 3
              IF(NIT(nb).EQ.2) THEN
                SPECIAL_BASIS_FLAG=2
                IF(DATAFILE) THEN
                  BASIS_STRING='c.Hermite*HermiteLagrange'
                ELSE
                  ISOCKET(1)=4
                  ISOCKET(2)=6
                  IF(FSKWRITE(ISOCKET,SK_LONG_INT,2,CONNID2).EQ.-1)
     '              GOTO 9999
                ENDIF
              ELSE
                SPECIAL_BASIS_FLAG=-1
              ENDIF
            ENDIF
          ENDIF
C CPB 3/8/95 extending sector bases
        ELSE IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN !sector
C        ELSE IF(IBT(1,ni,nb).EQ.5) THEN    !sector
C          IF((2.EQ.NIT(nb)).OR.(3.EQ.NIT(nb))) THEN
C            SPECIAL_BASIS_FLAG=3
C            INDEX_ni=IBT(2,ni,nb)
C          ELSE
C            SPECIAL_BASIS_FLAG=-1
C          ENDIF
          SPECIAL_BASIS_FLAG=3
          INDEX_ni=IBT(2,ni,nb)
        ELSE
          SPECIAL_BASIS_FLAG=-1
        ENDIF
        IF((SPECIAL_BASIS_FLAG.EQ.0).OR.(SPECIAL_BASIS_FLAG.EQ.3)) THEN
          IF(DATAFILE) THEN
            IF(ni.EQ.1) THEN
              BASIS_STRING=BASE(INDEX_ni)
            ELSE IF(ni.GT.1) THEN
              CALL STRING_TRIM(BASIS_STRING,IBEG,IEND)
              BASIS_STRING=BASIS_STRING(IBEG:IEND)//'*'//BASE(INDEX_ni)
            ENDIF
            IF(ni.EQ.SIMPLEX_BASIS_DIRNS(1)) THEN
              CALL STRING_TRIM(BASIS_STRING,IBEG,IEND)
              BASIS_STRING=BASIS_STRING(IBEG:IEND)//SIMPLEX_BASIS_LIST(
     '          IBEG2:IEND2)
            ENDIF
          ELSE
            IF(FSKWRITE(IBASE(INDEX_ni),SK_LONG_INT,1,CONNID2).EQ.-1)
     '        GOTO 9999
          ENDIF
          ni=ni+1
        ENDIF
      ENDDO
C KAT 26Jan00: Fergusson bases
      IF(SPECIAL_BASIS_FLAG.EQ.0.AND.NBCD(nb).EQ.1) SPECIAL_BASIS_FLAG=4

      CALL EXITS('CALCULATE_BASIS_STRING')
      RETURN
 9999 CALL ERRORS('CALCULATE_BASIS_STRING',ERROR)
      CALL EXITS('CALCULATE_BASIS_STRING')
      RETURN 1

      END


