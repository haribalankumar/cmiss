      SUBROUTINE NORMAL(ne,nr,NW,XG,XN,INTERFACE,ERROR,*)

C#### Subroutine: NORMAL
C###  Description:
C###    NORMAL finds the unit normal vector XN(nj) in the cartesian
C###    reference frame. Note: Normal is outward if nodes defined such
C###    that domain is kept to the right. The normal can be reversed
C###    in specified elements using the define normal command.

C**** If INTERFACE is true then the domain has been kept to the left
C**** (as in a coupled problem for regions sharing a node) so the
C**** normal must be corrected for this.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Paramter List
      INTEGER ne,nr,NW(NEM,3)
      REAL*8 XG(NJM,NUM),XN(NJM)
      CHARACTER ERROR*(*)
      LOGICAL INTERFACE
!     Local Variables
      INTEGER nj
      REAL*8 XGT1(3),XGT2(3),VLENGTH

      CALL ENTERS('NORMAL',*9999)

      GOTO (10,20,30)ITYP10(nr)
10    CONTINUE
        XGT1(1)=XG(1,2)
        XGT1(2)=XG(2,2)
        IF(NJ_LOC(NJL_GEOM,0,nr).GT.2) THEN
          XGT1(3)=XG(3,2)
          XGT2(1)=XG(1,4)
          XGT2(2)=XG(2,4)
          XGT2(3)=XG(3,4)
        ENDIF
        GOTO 100
20    CONTINUE
C  stuff from here down to the next comment needs to be checked.
C  16/5/89 R.H.H.
        XGT1(1)=XG(1,2)*DCOS(XG(1,2))-XG(1,1)*DSIN(XG(1,2))*XG(2,2)
        XGT1(2)=XG(1,2)*DSIN(XG(1,2))+XG(1,1)*DCOS(XG(1,2))*XG(2,2)
        XGT1(3)=XG(3,2)
        XGT2(1)=XG(1,4)*DCOS(XG(1,2))-XG(1,1)*DSIN(XG(1,2))*XG(2,4)
        XGT2(2)=XG(1,4)*DSIN(XG(1,2))+XG(1,1)*DCOS(XG(1,2))*XG(2,4)
        XGT2(3)=XG(3,4)
        GOTO 100
30    CONTINUE
        XGT1(1)=XG(1,2)*DCOS(XG(1,3))*DCOS(XG(1,2))
     '         -XG(1,1)*DSIN(XG(1,3))*DCOS(XG(1,2))*XG(3,2)
     '         -XG(1,1)*DCOS(XG(1,3))*DSIN(XG(1,2))*XG(2,2)
        XGT1(2)=XG(1,2)*DCOS(XG(1,3))*DSIN(XG(1,2))
     '         -XG(1,1)*DSIN(XG(1,3))*DSIN(XG(1,2))*XG(3,2)
     '         +XG(1,1)*DCOS(XG(1,3))*DCOS(XG(1,2))*XG(2,2)
        XGT1(3)=XG(1,2)*DSIN(XG(1,3))+XG(1,1)*DCOS(XG(1,3))*XG(3,2)
        XGT2(1)=XG(1,4)*DCOS(XG(1,3))*DCOS(XG(1,2))
     '         -XG(1,1)*DSIN(XG(1,3))*DCOS(XG(1,2))*XG(3,4)
     '         -XG(1,1)*DCOS(XG(1,3))*DSIN(XG(1,2))*XG(2,4)
        XGT2(2)=XG(1,4)*DCOS(XG(1,3))*DSIN(XG(1,2))
     '         -XG(1,1)*DSIN(XG(1,3))*DSIN(XG(1,2))*XG(3,4)
     '         +XG(1,1)*DCOS(XG(1,3))*DCOS(XG(1,2))*XG(2,4)
        XGT2(3)=XG(1,4)*DSIN(XG(1,3))+XG(1,1)*DCOS(XG(1,3))*XG(3,4)
        GOTO 100
100   CONTINUE
C
C  Finds normal components XN(nj)
C
      GOTO (110,120,130)NJ_LOC(NJL_GEOM,0,nr)
110   CONTINUE
        GOTO 200
120   CONTINUE
        XN(1)=-XGT1(2)
        XN(2)=XGT1(1)
        VLENGTH=DSQRT(XN(1)*XN(1)+XN(2)*XN(2))
        IF(VLENGTH.LE.ZERO_TOL) THEN
          ERROR='>>Zero normal vector length'
          GOTO 9999
        ENDIF
        XN(1)=XN(1)/VLENGTH
        XN(2)=XN(2)/VLENGTH
        IF(INTERFACE) THEN
          XN(1)=-XN(1)
          XN(2)=-XN(2)
        ENDIF
        GOTO 200
130   CONTINUE
        XN(1)=XGT1(2)*XGT2(3)-XGT1(3)*XGT2(2)
        XN(2)=XGT1(3)*XGT2(1)-XGT1(1)*XGT2(3)
        XN(3)=XGT1(1)*XGT2(2)-XGT1(2)*XGT2(1)
        VLENGTH=DSQRT(XN(1)*XN(1)+XN(2)*XN(2)+XN(3)*XN(3))
C New AJP 15-3-94
        IF(VLENGTH.LE.ZERO_TOL) THEN
          ERROR='>>Zero normal vector length'
          GOTO 9999
        ENDIF
        XN(1)=XN(1)/VLENGTH
        XN(2)=XN(2)/VLENGTH
        XN(3)=XN(3)/VLENGTH
        IF(INTERFACE) THEN
          XN(1)=-XN(1)
          XN(2)=-XN(2)
          XN(3)=-XN(3)
        ENDIF
        GOTO 200
200   CONTINUE

      IF(NW(ne,3).EQ.1) THEN !normal reversal specified
        XN(1)=-XN(1)
        XN(2)=-XN(2)
        IF(NJT.GE.3) XN(3)=-XN(3)
      ENDIF

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(NORMAL_1)
         WRITE(OP_STRING,*)' ne=',ne,' normal=',(XN(nj),nj=1,NJT)
         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(NORMAL_1)

      ENDIF

      CALL EXITS('NORMAL')
      RETURN
9999  CALL ERRORS('NORMAL',ERROR)
      CALL EXITS('NORMAL')
      RETURN 1
      END


