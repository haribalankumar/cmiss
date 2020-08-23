      SUBROUTINE ITERATE_DYNAM3(ISEG,CSEG,ECHO,END,STRING,INTWORK,
     '  REALWORK,ERROR,*)

C#### Subroutine: ITERATE_DYNAM3
C###  Description:
C###    ITERATE_DYNAM3 handles all coupled bidomain iteration
C###    calls to FEM.
C*** Written by Martin Buist, 28th February 2000

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'iter00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER ISEG(*),INTWORK(*)
      REAL*8 REALWORK(*)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL ECHO,END
!     Local Variables
      INTEGER I,IBEG(14),IBEG3,IEND(14),IEND3,it_count
      REAL*8 H1,H2
      CHARACTER CHAR2*2,CHAR3*12,USERDEFVARS(20)*(9)
      LOGICAL FINISHED,FINISHED1,FINISHED2

      CALL ENTERS('ITERATE_DYNAM3',*9999)

      IF(KTYP20.EQ.3) THEN !Circular Bidomain Iteration

        !Find equivalent numbers to user defined expressions
        USERDEFVARS(1)='PLAPLACE'
        USERDEFVARS(2)='TLAPLACE'
        USERDEFVARS(3)='BLAPLACE'
        USERDEFVARS(4)='TMEMBRNE'
        USERDEFVARS(5)='EXTCELLR'
        USERDEFVARS(6)=' '
        USERDEFVARS(7)='CAVITY'
        USERDEFVARS(8)='BLOOD'
        USERDEFVARS(9)='BEMHT'
        USERDEFVARS(10)='MUSCLE'
        USERDEFVARS(11)='GRIDHT'

        CALL USER(USERDEFVARS(1),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(2),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(3),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(4),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(5),' ;,.:=+-*/()[]^',ERROR,*9999)
C        CALL USER(USERDEFVARS(6),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(7),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(8),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(9),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(10),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(11),' ;,.:=+-*/()[]^',ERROR,*9999)

        IF(DOP) THEN
          WRITE(OP_STRING,'('' PLAPLACE = '',A)') USERDEFVARS(1)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' TLAPLACE = '',A)') USERDEFVARS(2)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' BLAPLACE = '',A)') USERDEFVARS(3)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' TMEMBRNE = '',A)') USERDEFVARS(4)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' EXTCELLR = '',A)') USERDEFVARS(5)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          WRITE(OP_STRING,'('' EXTUPDTE = '',A)') USERDEFVARS(6)(1:1)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' CAVITY   = '',A)') USERDEFVARS(7)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' BLOOD    = '',A)') USERDEFVARS(8)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' BEMHT    = '',A)') USERDEFVARS(9)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' MUSCLE   = '',A)') USERDEFVARS(10)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' GRIDHT   = '',A)') USERDEFVARS(11)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IT_COUNT=1
        FINISHED=.FALSE.
        FINISHED1=.FALSE.
        FINISHED2=.FALSE.
        DO WHILE(.NOT.FINISHED)

          !=============================================!
          !fem solve activ coupled extracellular restart!
          !  reg GRIDHT,CAVITY,BLOOD                    !
          !  class TMEMBRNE,EXTCELLR,TLAPLACE,BLAPLACE  !
          !  noupdate static alpha ??                   !
          !=============================================!

          WRITE(CHAR3,'(D11.4)') ITER_ALPHA2
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='ACTIVATION'
          CO(4)='COUPLED'
          CO(5)='EXTRACELLULAR'
          CO(6)='RESTART'
          CO(7)='REGION'
          CO(8)=USERDEFVARS(11)(1:1)//','//USERDEFVARS(7)(1:1)
     '      //','//USERDEFVARS(8)(1:1)
          CO(9)='CLASS'
          CO(10)=USERDEFVARS(4)(1:1)//','//USERDEFVARS(5)(1:1)
     '      //','//USERDEFVARS(2)(1:1)
     '      //','//USERDEFVARS(3)(1:1)
          CO(11)='NOUPDATE_VM'
          CO(12)='STATIC'
          CO(13)='ALPHA'
          CO(14)=CHAR3(IBEG3:IEND3)
          NTCO=14
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !===========================================!
          !fem solve activ coupled both restart       !
          !  reg GRIDHT,CAVITY,BLOOD                  !
          !  class TMEMBRNE,EXTCELLR,TLAPLACE,BLAPLACE!
          !  noupdate static alpha ??                 !
          !===========================================!

          WRITE(CHAR3,'(D11.4)') ITER_ALPHA2
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='ACTIVATION'
          CO(4)='COUPLED'
          CO(5)='BOTH'
          CO(6)='RESTART'
          CO(7)='REGION'
          CO(8)=USERDEFVARS(11)(1:1)//','//USERDEFVARS(7)(1:1)
     '      //','//USERDEFVARS(8)(1:1)
          CO(9)='CLASS'
          CO(10)=USERDEFVARS(4)(1:1)//','//USERDEFVARS(5)(1:1)
     '      //','//USERDEFVARS(2)(1:1)
     '      //','//USERDEFVARS(3)(1:1)
          CO(11)='NOUPDATE_VM'
          CO(12)='STATIC'
          CO(13)='ALPHA'
          CO(14)=CHAR3(IBEG3:IEND3)
          NTCO=14
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !====================================================!
          !fem update init from grid region CAVITY,BEMHT,GRIDHT!
          !  class EXTCELLR,TLAPLACE,EXTCELLR alpha ??         !
          !====================================================!

          WRITE(CHAR3,'(D11.4)') ITER_ALPHA1
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='INITIAL'
          CO(4)='FROM'
          CO(5)='GRID'
          CO(6)='REGION'
          CO(7)=USERDEFVARS(7)(1:1)//','//USERDEFVARS(9)(1:1)//','
     '      //USERDEFVARS(11)(1:1)
          CO(8)='CLASS'
          CO(9)=USERDEFVARS(5)(1:1)//','//USERDEFVARS(2)(1:1)//','
     '      //USERDEFVARS(5)(1:1)
          CO(10)='ALPHA'
          CO(11)=CHAR3(IBEG3:IEND3)
          NTCO=11
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !===================================================!
          !fem update init from grid region BLOOD,BEMHT,GRIDHT!
          !  class EXTCELLR,BLAPLACE,EXTCELLR alpha ??        !
          !===================================================!

          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='INITIAL'
          CO(4)='FROM'
          CO(5)='GRID'
          CO(6)='REGION'
          CO(7)=USERDEFVARS(8)(1:1)//','//USERDEFVARS(9)(1:1)//','
     '      //USERDEFVARS(11)(1:1)
          CO(8)='CLASS'
          CO(9)=USERDEFVARS(5)(1:1)//','//USERDEFVARS(3)(1:1)//','
     '      //USERDEFVARS(5)(1:1)
          CO(10)='ALPHA'
          CO(11)=CHAR3(IBEG3:IEND3)
          NTCO=11
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !==================================================!
          !fem solve coupled reg CAVITY,MUSCLE class TLAPLACE!
          !==================================================!

          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='COUPLED'
          CO(4)='REGION'
          CO(5)=USERDEFVARS(7)(1:1)//','//USERDEFVARS(10)(1:1)
          CO(6)='CLASS'
          CO(7)=USERDEFVARS(2)(1:1)
          NTCO=7
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !==================================!
          !fem solve reg BLOOD class BLAPLACE!
          !==================================!

          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='REGION'
          CO(4)=USERDEFVARS(8)(1:1)
          CO(5)='CLASS'
          CO(6)=USERDEFVARS(3)(1:1)
          NTCO=6
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !====================================================!
          !fem check conv grid_region BEMHT grid_class EXTCELLR!
          !  bem_region BLOOD bem_class BLAPLACE               !
          !  grid_update_class EXTCELLR                        !
          !====================================================!

          CO(1)='FEM'
          CO(2)='CHECK'
          CO(3)='CONVERGENCE'
          CO(4)='GRID_REGION'
          CO(5)=USERDEFVARS(9)(1:1)
          CO(6)='GRID_CLASS'
          CO(7)=USERDEFVARS(5)(1:1)
          CO(8)='BEM_REGION'
          CO(9)=USERDEFVARS(8)(1:1)
          CO(10)='BEM_CLASS'
          CO(11)=USERDEFVARS(3)(1:1)
          CO(12)='GRID_UPDATE_CLASS'
          CO(13)=USERDEFVARS(5)(1:1)
          NTCO=13
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          IF((YPYQPOTE.LE.ITER_TOL(1)).AND.
     '      (YPYQFLUX.LE.ITER_TOL(2))) THEN
            FINISHED1=.TRUE.
          ELSE
            FINISHED1=.FALSE.
          ENDIF
          IF(ECHO) THEN
            WRITE(OP_STRING,'(''  Potential,Flux '',2F12.6)')
     '        YPYQPOTE,YPYQFLUX
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          H1=YPYQPOTE
          H2=YPYQFLUX

          !====================================================!
          !fem check conv grid_region BEMHT grid_class EXTCELLR!
          !  bem_region CAVITY bem_class TLAPLACE              !
          !  grid_update_class EXTCELLR                        !
          !====================================================!

          CO(1)='FEM'
          CO(2)='CHECK'
          CO(3)='CONVERGENCE'
          CO(4)='GRID_REGION'
          CO(5)=USERDEFVARS(9)(1:1)
          CO(6)='GRID_CLASS'
          CO(7)=USERDEFVARS(5)(1:1)
          CO(8)='BEM_REGION'
          CO(9)=USERDEFVARS(7)(1:1)
          CO(10)='BEM_CLASS'
          CO(11)=USERDEFVARS(2)(1:1)
          CO(12)='GRID_UPDATE_CLASS'
          CO(13)=USERDEFVARS(5)(1:1)
          NTCO=13
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          IF((YPYQPOTE.LE.ITER_TOL(1)).AND.
     '      (YPYQFLUX.LE.ITER_TOL(2))) THEN
            FINISHED2=.TRUE.
          ELSE
            FINISHED2=.FALSE.
          ENDIF
          IF(ECHO) THEN
            WRITE(OP_STRING,'(''  Potential,Flux '',2F12.6)')
     '        YPYQPOTE,YPYQFLUX
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          IF(FINISHED1.AND.FINISHED2) THEN
            FINISHED=.TRUE.
            WRITE(OP_STRING,'('' Solution converged in '
     '        //' in '',I4,'' iterations'')') IT_COUNT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSEIF(IT_COUNT.GE.NUM_ITS) THEN
            FINISHED=.TRUE.
            WRITE(OP_STRING,'('' Solution did not converge '
     '        //' in '',I4,'' iterations'')') NUM_ITS
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Inside  '',2F12.6)') H1,H2
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Outside '',2F12.6)') YPYQPOTE,YPYQFLUX
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          IT_COUNT=IT_COUNT+1

        ENDDO

      ELSEIF(KTYP20.EQ.4) THEN !Simple torso slice iteration
        CALL ASSERT(.FALSE.,'>>Not implemented',ERROR,*9999)
      ELSEIF(KTYP20.EQ.5) THEN !Full torso slice iteration
        !Find equivalent numbers to user defined expressions
        USERDEFVARS(1)='PLAPLACE'
        USERDEFVARS(2)='TLAPLACE'
        USERDEFVARS(3)='BLAPLACE'
        USERDEFVARS(4)='TMEMBRNE'
        USERDEFVARS(5)='EXTCELLR'
        USERDEFVARS(6)='CAVITY'
        USERDEFVARS(7)='VENTCL'
        USERDEFVARS(8)='BEMHRT'
        USERDEFVARS(9)='GRIDHT'

        CALL USER(USERDEFVARS(1),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(2),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(3),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(4),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(5),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(6),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(7),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(8),' ;,.:=+-*/()[]^',ERROR,*9999)
        CALL USER(USERDEFVARS(9),' ;,.:=+-*/()[]^',ERROR,*9999)

        IF(DOP) THEN
          WRITE(OP_STRING,'('' PLAPLACE = '',A)') USERDEFVARS(1)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' TLAPLACE = '',A)') USERDEFVARS(2)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' BLAPLACE = '',A)') USERDEFVARS(3)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' TMEMBRNE = '',A)') USERDEFVARS(4)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' EXTCELLR = '',A)') USERDEFVARS(5)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' CAVITY   = '',A)') USERDEFVARS(6)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' VENTCK   = '',A)') USERDEFVARS(7)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' BEMHRT   = '',A)') USERDEFVARS(8)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' GRIDHT   = '',A)') USERDEFVARS(9)(1:1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IT_COUNT=1
        FINISHED=.FALSE.
        FINISHED1=.FALSE.
        FINISHED2=.FALSE.
        DO WHILE(.NOT.FINISHED)

          !=============================================!
          !fem solve activ coupled extracellular restart!
          !  reg GRIDHT,CAVITY,VENTCL                   !
          !  class TMEMBRNE,EXTCELLR,TLAPLACE,BLAPLACE  !
          !  noupdate static alpha ??                   !
          !=============================================!

          WRITE(CHAR3,'(D11.4)') ITER_ALPHA2
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='ACTIVATION'
          CO(4)='COUPLED'
          CO(5)='EXTRACELLULAR'
          CO(6)='RESTART'
          CO(7)='REGION'
          CO(8)=USERDEFVARS(9)(1:1)//','//USERDEFVARS(6)(1:1)
     '      //','//USERDEFVARS(7)(1:1)
          CO(9)='CLASS'
          CO(10)=USERDEFVARS(4)(1:1)//','//USERDEFVARS(5)(1:1)
     '      //','//USERDEFVARS(2)(1:1)
     '      //','//USERDEFVARS(3)(1:1)
          CO(11)='NOUPDATE_VM'
          CO(12)='STATIC'
          CO(13)='ALPHA'
          CO(14)=CHAR3(IBEG3:IEND3)
          NTCO=14
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !===========================================!
          !fem solve activ coupled both restart       !
          !  reg GRIDHT,CAVITY,VENTCL                 !
          !  class TMEMBRNE,EXTCELLR,TLAPLACE,BLAPLACE!
          !  noupdate static alpha ??                 !
          !===========================================!

          WRITE(CHAR3,'(D11.4)') ITER_ALPHA2
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='ACTIVATION'
          CO(4)='COUPLED'
          CO(5)='BOTH'
          CO(6)='RESTART'
          CO(7)='REGION'
          CO(8)=USERDEFVARS(9)(1:1)//','//USERDEFVARS(6)(1:1)
     '      //','//USERDEFVARS(7)(1:1)
          CO(9)='CLASS'
          CO(10)=USERDEFVARS(4)(1:1)//','//USERDEFVARS(5)(1:1)
     '      //','//USERDEFVARS(2)(1:1)
     '      //','//USERDEFVARS(3)(1:1)
          CO(11)='NOUPDATE_VM'
          CO(12)='STATIC'
          CO(13)='ALPHA'
          CO(14)=CHAR3(IBEG3:IEND3)
          NTCO=14
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !=====================================================!
          !fem update init from grid region CAVITY,BEMHRT,GRIDHT!
          !  class EXTCELLR,TLAPLACE,EXTCELLR alpha ??          !
          !=====================================================!

          WRITE(CHAR3,'(D11.4)') ITER_ALPHA1
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='INITIAL'
          CO(4)='FROM'
          CO(5)='GRID'
          CO(6)='REGION'
          CO(7)=USERDEFVARS(6)(1:1)//','//USERDEFVARS(8)(1:1)//','
     '      //USERDEFVARS(9)(1:1)
          CO(8)='CLASS'
          CO(9)=USERDEFVARS(5)(1:1)//','//USERDEFVARS(2)(1:1)//','
     '      //USERDEFVARS(5)(1:1)
          CO(10)='ALPHA'
          CO(11)=CHAR3(IBEG3:IEND3)
          NTCO=11
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !=====================================================!
          !fem update init from grid region VENTCL,BEMHRT,GRIDHT!
          !  class EXTCELLR,BLAPLACE,EXTCELLR alpha ??          !
          !=====================================================!

          CO(1)='FEM'
          CO(2)='UPDATE'
          CO(3)='INITIAL'
          CO(4)='FROM'
          CO(5)='GRID'
          CO(6)='REGION'
          CO(7)=USERDEFVARS(7)(1:1)//','//USERDEFVARS(8)(1:1)//','
     '      //USERDEFVARS(9)(1:1)
          CO(8)='CLASS'
          CO(9)=USERDEFVARS(5)(1:1)//','//USERDEFVARS(3)(1:1)//','
     '      //USERDEFVARS(5)(1:1)
          CO(10)='ALPHA'
          CO(11)=CHAR3(IBEG3:IEND3)
          NTCO=11
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !===================================!
          !fem solve reg CAVITY class TLAPLACE!
          !===================================!

          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='REGION'
          CO(4)=USERDEFVARS(6)(1:1)
          CO(5)='CLASS'
          CO(6)=USERDEFVARS(2)(1:1)
          NTCO=6
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !===================================!
          !fem solve reg VENTCL class BLAPLACE!
          !===================================!

          CO(1)='FEM'
          CO(2)='SOLVE'
          CO(3)='REGION'
          CO(4)=USERDEFVARS(7)(1:1)
          CO(5)='CLASS'
          CO(6)=USERDEFVARS(3)(1:1)
          NTCO=6
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          !=====================================================!
          !fem check conv grid_region BEMHRT grid_class EXTCELLR!
          !  bem_region VENCTL bem_class BLAPLACE               !
          !  grid_update_class EXTCELLR                         !
          !=====================================================!

          CO(1)='FEM'
          CO(2)='CHECK'
          CO(3)='CONVERGENCE'
          CO(4)='GRID_REGION'
          CO(5)=USERDEFVARS(8)(1:1)
          CO(6)='GRID_CLASS'
          CO(7)=USERDEFVARS(5)(1:1)
          CO(8)='BEM_REGION'
          CO(9)=USERDEFVARS(7)(1:1)
          CO(10)='BEM_CLASS'
          CO(11)=USERDEFVARS(3)(1:1)
          CO(12)='GRID_UPDATE_CLASS'
          CO(13)=USERDEFVARS(5)(1:1)
          NTCO=13
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          IF((YPYQPOTE.LE.ITER_TOL(1)).AND.
     '      (YPYQFLUX.LE.ITER_TOL(2))) THEN
            FINISHED1=.TRUE.
          ELSE
            FINISHED1=.FALSE.
          ENDIF
          IF(ECHO) THEN
            WRITE(OP_STRING,'(''  Potential,Flux '',2F12.6)')
     '        YPYQPOTE,YPYQFLUX
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          H1=YPYQPOTE
          H2=YPYQFLUX

          !=====================================================!
          !fem check conv grid_region BEMHRT grid_class EXTCELLR!
          !  bem_region CAVITY bem_class TLAPLACE               !
          !  grid_update_class EXTCELLR                         !
          !=====================================================!

          CO(1)='FEM'
          CO(2)='CHECK'
          CO(3)='CONVERGENCE'
          CO(4)='GRID_REGION'
          CO(5)=USERDEFVARS(8)(1:1)
          CO(6)='GRID_CLASS'
          CO(7)=USERDEFVARS(5)(1:1)
          CO(8)='BEM_REGION'
          CO(9)=USERDEFVARS(6)(1:1)
          CO(10)='BEM_CLASS'
          CO(11)=USERDEFVARS(2)(1:1)
          CO(12)='GRID_UPDATE_CLASS'
          CO(13)=USERDEFVARS(5)(1:1)
          NTCO=13
          NTCOQU(3)=0
          DO i=1,NTCO
            CALL STRING_TRIM(CO(i),IBEG(i),IEND(i))
          ENDDO
          WRITE(CHAR2,'(I2)') NTCO
          IF(ECHO) THEN
            WRITE(OP_STRING,'('' *'','//CHAR2(1:2)//'A)')
     '        (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(STRING,'('//CHAR2(1:2)//'A)')
     '      (CO(i)(IBEG(i):IEND(i)+1),i=1,NTCO)
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          IF((YPYQPOTE.LE.ITER_TOL(1)).AND.
     '      (YPYQFLUX.LE.ITER_TOL(2))) THEN
            FINISHED2=.TRUE.
          ELSE
            FINISHED2=.FALSE.
          ENDIF
          IF(ECHO) THEN
            WRITE(OP_STRING,'(''  Potential,Flux '',2F12.6)')
     '        YPYQPOTE,YPYQFLUX
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          IF(FINISHED1.AND.FINISHED2) THEN
            FINISHED=.TRUE.
            WRITE(OP_STRING,'('' Solution converged in '
     '        //' in '',I4,'' iterations'')') IT_COUNT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSEIF(IT_COUNT.GE.NUM_ITS) THEN
            FINISHED=.TRUE.
            WRITE(OP_STRING,'('' Solution did not converge '
     '        //' in '',I4,'' iterations'')') NUM_ITS
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Inside  '',2F12.6)') H1,H2
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Outside '',2F12.6)') YPYQPOTE,YPYQFLUX
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          IT_COUNT=IT_COUNT+1

        ENDDO
      ENDIF

      CALL EXITS('ITERATE_DYNAM3')
      RETURN
 9999 CALL ERRORS('ITERATE_DYNAM3',ERROR)
      CALL EXITS('ITERATE_DYNAM3')
      RETURN 1
      END


