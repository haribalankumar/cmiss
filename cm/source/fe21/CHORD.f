      SUBROUTINE CHORD(NCOY,NKJ,NLCHOR,NPL,
     '  DL,XP,ZPOS,STRING,X1STR,X2STR,ZCHAR,ERROR,*)

C#### Subroutine: CHORD
C###  Description:
C###    CHORD defines individual chord at ZCHAR=ZPOS.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'user00.cmn'
!     Parameter List
      INTEGER NCOY,NKJ(NJM,NPM),NLCHOR(0:10,NRM),
     '  NPL(5,0:3,NLM)
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM),ZPOS
      CHARACTER ERROR*(*),STRING*(MXCH),X1STR*(*),X2STR*(*),ZCHAR*(*)
!     Local Variables
      INTEGER IBEG,IEND,N,nb,nj,np(10),NTITLV(1),NTLV,NTYSTR
      REAL*8 DY(10),RA(1),X1,X2,XDIFF,Z(3,10)
      CHARACTER YSTR*100

      CALL ENTERS('CHORD',*9999)
      CO(noco+1)=ZCHAR
      WRITE(CO(noco+2),'(E11.4)') ZPOS
      CALL ASSIGN(STRING,noco,CO,COQU,ERROR,*9999)

      CALL USER(X1STR,' ;,.:=+-*/()[]^',ERROR,*9999)
      CALL PARSTR(X1STR,1,NTLV,NTITLV,1,RA,ERROR,*9999)
      X1=RA(1)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) ' ZCHAR= ',ZCHAR
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) ' CO(noco+2)= ',CO(noco+2)(1:11)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) ' X1STR= ',X1STR
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) ' X1= ',X1
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL USER(X2STR,' ;,.:=+-*/()[]^',ERROR,*9999)
      CALL PARSTR(X2STR,1,NTLV,NTITLV,1,RA,ERROR,*9999)
      X2=RA(1)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) ' X2STR= ',X2STR
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) ' X2= ',X2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(NCOY.GT.0) THEN
        NTYSTR=1+NTCOQU(NCOY+1)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,*) ' no of y_strings= ',NTYSTR
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DO N=1,NTYSTR
          IF(N.EQ.1) THEN
            CALL STRING_TRIM(CO(NCOY+1),IBEG,IEND)
            YSTR=CO(NCOY+1)(IBEG:IEND)
          ELSE
            CALL STRING_TRIM(COQU(NCOY+1,N-1),IBEG,IEND)
            YSTR=COQU(NCOY+1,N-1)(IBEG:IEND)
          ENDIF
          CALL USER(YSTR,' ;,.:=+-*/()[]^',ERROR,*9999)
          CALL PARSTR(YSTR,1,NTLV,NTITLV,1,RA,ERROR,*9999)
          Z(2,N)=RA(1)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*) ' YSTR= ',YSTR
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) ' Y= ',Z(2,N)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDDO
      ELSE
        NTYSTR=2
        Z(2,1)=0.d0
        Z(2,2)=0.d0
      ENDIF
      DY(1)=Z(2,2)-Z(2,1)
      DY(NTYSTR)=Z(2,NTYSTR)-Z(2,NTYSTR-1)
      DO N=2,NTYSTR-1
        DY(N)=0.d0
      ENDDO

      NTUS=NTUS-1

      nb=1
      XDIFF=(X2-X1)/DBLE(NTYSTR-1)
      DO N=1,NTYSTR
        NPT(1)=NPT(1)+1
        np(N)=NPT(1)
        DO nj=1,NJT
          NKJ(nj,NPT(1))=NKT(0,nb)
        ENDDO
        Z(1,N)=X1+DBLE(N-1)*XDIFF
        Z(3,N)=ZPOS
        XP(1,1,1,NPT(1))=Z(1,N)
        XP(2,1,1,NPT(1))=XDIFF
        XP(1,1,2,NPT(1))=Z(2,N)
        XP(2,1,2,NPT(1))=DY(N)
        XP(1,1,3,NPT(1))=Z(3,N)
        XP(2,1,3,NPT(1))=0.d0
      ENDDO
      NTCHOR=NTCHOR+1
      NLCHOR(0,NTCHOR)=NTYSTR-1
      DO N=1,NTYSTR-1
        NLT=NLT+1
        NLCHOR(N,NTCHOR)=NLT
        NPL(1,0,NLT)=1
        NPL(1,1,NLT)=4
        NPL(1,2,NLT)=4
        NPL(1,3,NLT)=4
        NPL(2,1,NLT)=np(N)
        NPL(3,1,NLT)=np(N+1)
        NPL(4,1,NLT)=2
        NPL(5,1,NLT)=2
        DL(1,NLT)=1.d0
        DL(2,NLT)=1.d0
        DL(3,NLT)=1.d0
      ENDDO
      CALL POLYLINE(1,1,NTYSTR,Z,ERROR,*9999)

      CALL EXITS('CHORD')
      RETURN
 9999 CALL ERRORS('CHORD',ERROR)
      CALL EXITS('CHORD')
      RETURN 1
      END


