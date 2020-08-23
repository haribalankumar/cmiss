      INTEGER FUNCTION SEGORDER(PARENT_ORDER,END,CONECT,ISEED)

C#### FUNCTION: ORDER
C#### Type: INTEGER
C###  Description:
C###    ORDER calcultes the order of a daughter vessel branching
C###    from a parent vessel

      IMPLICIT NONE
!     Parameter List
      INTEGER PARENT_ORDER,iseed
      REAL*8 CONECT(0:13),RANDOM,CM_RANDOM_NUMBER,SUM
      LOGICAL END

      SEGORDER=0
      SUM=0

      RANDOM=CM_RANDOM_NUMBER(iseed)

      IF (.NOT.END) THEN
        SEGORDER=0
        SUM=CONECT(0)
        DO WHILE (RANDOM.GT.SUM)
          SEGORDER=SEGORDER+1
          SUM=SUM+CONECT(SEGORDER)

        ENDDO
      ELSE !vessel end
        IF (RANDOM.LT.CONECT(12)) THEN
          SEGORDER=PARENT_ORDER-1
        ELSE
          SEGORDER=PARENT_ORDER
        ENDIF
      ENDIF

      RETURN
      END


