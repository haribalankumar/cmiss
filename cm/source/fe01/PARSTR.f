      SUBROUTINE PARSTR(STRING,NXLV,NTLV,NTITLV,NXRA,RA,ERROR,*)

C#### Subroutine: PARSTR
C###  Description:
C###    PARSTR parses STRING to evaluate various operators (+-*/) via
C###    call to CALCUL, and functions (sin,cos,tan,asin,acos,atan,
C###    sinh,cosh,tanh,log,exp,abs,sqrt). Then the nested list of
C###    items in STRING is broken into separate REAL*8 number stored
C###    in the REAL*8 array RA. The maximum level of nesting is
C###    returned as NTLV. The size of each dimension is carried in the
C###    first NTLV components of the integer array NTITLV. Each
C###    nesting is delimited by left and right brackets while each
C###    item is separated by commas.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NXLV,NTLV,NXRA,NTITLV(NXLV)
      REAL*8 RA(NXRA)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER MXLV
      PARAMETER (MXLV=16)
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,n1ch,n2ch,NOCHAR,
     '  NOITLV(MXLV),nolv,nora,NTRA
      REAL*8 R,RFROMC
      CHARACTER C3*10,C4*10,CHARTEMP*20,CHAR,STR*255,SUBST*500
      LOGICAL SCALAR

      CALL ENTERS('PARSTR',*9999)
      CALL STRING_TRIM(STRING,IBEG,IEND)
      STRING(1:)='('//STRING(IBEG:IEND)//')'
      CALL STRING_TRIM(STRING,IBEG1,IEND1)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Initial string in PARSTR = '',A)')
     '    STRING(IBEG1:IEND1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      nolv=0
      n2ch=IBEG1-1
      DO WHILE(nolv.GE.0)
        n2ch=n2ch+1
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' *** nolv='',I2,'' *** n2ch= '',I4)')
     '      nolv,n2ch
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        CHAR=STRING(n2ch:n2ch)
        IF(CHAR.EQ.'(') THEN
          nolv=nolv+1
          IF((nolv.GT.NXLV).OR.(nolv.GT.MXLV)) THEN
            ERROR=' Maximum array nesting exceeded'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=n2ch
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(''  left bracket at '',I4,'
     '        //''' nolv='',I2)') n2ch,nolv
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ELSE IF(CHAR.EQ.')') THEN
          n1ch=NOITLV(nolv)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' right bracket at '',I4,'
     '        //''' left bracket at '',I4,'' nolv= '',I2)')
     '        n2ch,n1ch,nolv
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) 'string to enter CALCUL= ',
     '        STRING(n1ch:n2ch)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          STR=STRING(n1ch:n2ch)
          CALL CALCUL(STR(1:n2ch-n1ch+1)//' ',SUBST,ERROR,*9999)
          CALL STRING_TRIM(SUBST,IBEG,IEND)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*) 'full string after CALCUL= ',STRING
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) 'SUBST=',SUBST(IBEG:IEND)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(SUBST(IBEG:IBEG).EQ.'[') THEN
            SCALAR=.FALSE.
            CALL PARSRA(SUBST,NXLV,NTLV,NTITLV,NXRA,RA,ERROR,*9999)
            SUBST(1:1)='['
          ELSE
            SCALAR=.TRUE.
            NTITLV(1)=1
            RA(1)=RFROMC(SUBST(IBEG:IBEG+14))
            SUBST(1:1)=' '
          ENDIF
          NTRA=NTITLV(1)
          IEND=1
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*) 'RA:',(RA(nora),nora=1,NTRA)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) 'full string to enter function loop=',
     '        STRING(1:60)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NTRA='',I2,'' n1ch='',I2,'' n2ch='','
     '        //'I2)') NTRA,n1ch,n2ch
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DO nora=1,NTRA
            R=RA(nora)
            NOCHAR=0
            IF(n1ch.GT.3)THEN
              CALL CUPPER(STRING(n1ch-3:n1ch-1),C3)
C              IF(CUPPER(STRING(n1ch-3:n1ch-1)).EQ.'EXP') THEN
              IF(C3.EQ.'EXP') THEN
                R=DEXP(R)
                NOCHAR=3
C              ELSE IF(CUPPER(STRING(n1ch-3:n1ch-1)).EQ.'ABS') THEN
              ELSE IF(C3.EQ.'ABS') THEN
                R=DABS(R)
                NOCHAR=3
C              ELSE IF(CUPPER(STRING(n1ch-3:n1ch-1)).EQ.'SIN') THEN
              ELSE IF(C3.EQ.'SIN') THEN
                R=DSIN(R)
                NOCHAR=3
C              ELSE IF(CUPPER(STRING(n1ch-3:n1ch-1)).EQ.'COS') THEN
              ELSE IF(C3.EQ.'COS') THEN
                R=DCOS(R)
                NOCHAR=3
C              ELSE IF(CUPPER(STRING(n1ch-3:n1ch-1)).EQ.'TAN') THEN
              ELSE IF(C3.EQ.'TAN') THEN
                R=TAN(R)
                NOCHAR=3
C              ELSE IF(CUPPER(STRING(n1ch-3:n1ch-1)).EQ.'LOG') THEN
              ELSE IF(C3.EQ.'LOG') THEN
                R=DLOG(R)
                NOCHAR=3
              ENDIF
              IF(NOCHAR.EQ.0.AND.n1ch.GT.4)THEN
                CALL CUPPER(STRING(n1ch-4:n1ch-1),C4)
C                IF(CUPPER(STRING(n1ch-4:n1ch-1)).EQ.'SQRT') THEN
                IF(C4.EQ.'SQRT') THEN
                  R=DSQRT(R)
                  NOCHAR=4
C                ELSE IF(CUPPER(STRING(n1ch-4:n1ch-1)).EQ.'ASIN') THEN
                ELSE IF(C4.EQ.'ASIN') THEN
                  R=DASIN(R)
                  NOCHAR=4
C                ELSE IF(CUPPER(STRING(n1ch-4:n1ch-1)).EQ.'ACOS') THEN
                ELSE IF(C4.EQ.'ACOS') THEN
                  R=DACOS(R)
                  NOCHAR=4
C                ELSE IF(CUPPER(STRING(n1ch-4:n1ch-1)).EQ.'ATAN') THEN
                ELSE IF(C4.EQ.'ATAN') THEN
                  R=DATAN(R)
                  NOCHAR=4
C                ELSE IF(CUPPER(STRING(n1ch-4:n1ch-1)).EQ.'SINH') THEN
                ELSE IF(C4.EQ.'SINH') THEN
                  R=DSINH(R)
                  NOCHAR=4
C                ELSE IF(CUPPER(STRING(n1ch-4:n1ch-1)).EQ.'COSH') THEN
                ELSE IF(C4.EQ.'COSH') THEN
                  R=DCOSH(R)
                  NOCHAR=4
C                ELSE IF(CUPPER(STRING(n1ch-4:n1ch-1)).EQ.'TANH') THEN
                ELSE IF(C4.EQ.'TANH') THEN
                  R=DTANH(R)
                  NOCHAR=4
                ENDIF
              ENDIF
            ENDIF
            WRITE(CHARTEMP,'(E15.8)') R
            SUBST=SUBST(1:IEND)//CHARTEMP//','
C news MPN 28Jul97: CHARTEMP is now 20 chars long (not 15)
            IEND=IEND+21
C old            IEND=IEND+16
          ENDDO
          IF(SCALAR) THEN
            SUBST(IEND:IEND)=' '
          ELSE
            SUBST(IEND:IEND)=']'
          ENDIF
          CALL STRING_TRIM(SUBST,IBEG,IEND)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*) 'SUBST=',SUBST(IBEG:IEND)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(NOCHAR.EQ.3) THEN
            STRING(n1ch-3:n1ch)='    '
          ELSE IF(NOCHAR.EQ.4) THEN
            STRING(n1ch-4:n1ch)='     '
          ENDIF
          CALL STRING_TRIM(STRING,IBEG1,IEND1)
          IF(n1ch.GT.IBEG1.AND.n2ch.LT.IEND1)THEN
            STRING(1:)=STRING(IBEG1:n1ch-1)//SUBST(IBEG:IEND)
     '        //STRING(n2ch+1:IEND1)
          ELSE IF(n2ch.LT.IEND1)THEN
            STRING(1:)=SUBST(IBEG:IEND)//STRING(n2ch+1:IEND1)
          ELSE IF(n1ch.GT.IBEG1)THEN
            STRING(1:)=STRING(IBEG1:n1ch-1)//SUBST(IBEG:IEND)
          ELSE
            STRING(1:)=SUBST(IBEG:IEND)
          ENDIF
          CALL STRING_TRIM(STRING,IBEG2,IEND2)
          n2ch=n2ch+IEND2-IEND1
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*)' full string at end of loop= ',
     '        STRING(1:IEND2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          nolv=nolv-1
          IF(nolv.EQ.0) nolv=-1
        ENDIF
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        OP_STRING(1)=STRING(1:IEND2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(STRING(1:1).EQ.'[') THEN
        CALL PARSRA(STRING,NXLV,NTLV,NTITLV,NXRA,RA,ERROR,*9999)
      ENDIF

      CALL EXITS('PARSTR')
      RETURN
 9999 CALL ERRORS('PARSTR',ERROR)
      CALL EXITS('PARSTR')
      RETURN 1
      END


