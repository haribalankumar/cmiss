      SUBROUTINE DEXI_POINT(IBT,IDO,INP,LD,NBJ,
     '  ne,NITB,nr,FACTOR_TOL,XE,XI,XID,ZD,USE_LOOSE_TOLERANCE,ERROR,*)

C#### Subroutine: DEXI_POINT
C###  Description:
C###    DEXI_POINT defines data Xi points from a global position
C###    in cartesian coordinates within a host mesh.
C###  USE_LOOSE_TOLERANCE allows a lower convergence tolerance, where
C###  the default is 10*converg_tol, or FACTOR_TOL if specified by user.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),LD,NBJ(NJM,NEM),
     '  ne,NITB,nr
      REAL*8 FACTOR_TOL,XE(NSM,NJM),XID(NIM),ZD(NJM)
      LOGICAL USE_LOOSE_TOLERANCE
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER count,i,j,D_INDEX(3),INFO,IPIV(3),MAX_ITER,nb,
     '  NB_XI(3),ni,nj,NJ_XI(3),nn,ns
      REAL*8  A(3,3),B(3),CONVERG_TEST,
     '  DATA(3),OLD_SUM,PXI,SUM,XI(3)
      LOGICAL FITTED,OVER,UNDER
      CHARACTER TRANS
      PARAMETER (TRANS='N')
      EXTERNAL DGETRS,DGETRF

      CALL ENTERS('DEXI_POINT',*9999)

      DO ni=1,NITB
        NJ_XI(ni)=ni
        NB_XI(ni)=NBJ(NJ_XI(ni),ne)
      ENDDO

      D_INDEX(1)=2 !index when calling PXI for derivative wrt Xi1
      D_INDEX(2)=4 !index when calling PXI for derivative wrt Xi2
      D_INDEX(3)=7 !index when calling PXI for derivative wrt Xi3

      CALL COORD(1,ITYP10(nr),ZD,DATA,ERROR,*9999)

      IF(ITYP10(nr).EQ.4) THEN !Prolate Coordinates
        nj=NJ_LOC(NJL_GEOM,3,1) !theta
        nb=NBJ(nj,ne)
        OVER=.FALSE.
        UNDER=.FALSE.
        ns=1
        DO nn=1,NNT(nb)
          IF (XE(ns,nj).GT.(2.0d0*PI))THEN
            OVER=.TRUE.
          ENDIF
          IF (XE(ns,nj).LT.DATA(nj))THEN
            UNDER=.TRUE.
          ENDIF
          ns=ns+NKT(nn,nb)
        ENDDO
        IF((.NOT.UNDER).AND.OVER) THEN
          DATA(nj)= DATA(nj)+(2.0d0*PI)
        ENDIF
      ENDIF

      SUM=1.0d0
      MAX_ITER=40
! PJH 17Mar97 add count check for nonconvergence
      COUNT=0
      CONVERG_TEST=2.0d0*CONVERG_TOL

C initialize CONVERG_TEST so loop will be entered

      DO WHILE ((CONVERG_TEST.GT.CONVERG_TOL)
     '  .AND.COUNT.LT.MAX_ITER)
        COUNT=COUNT+1
        DO i=1,NITB
          DO j=1,NITB
            A(i,j)=PXI(IBT(1,1,NB_XI(i)),IDO(1,1,0,NB_XI(i)),
     '        INP(1,1,NB_XI(i)),NB_XI(i),D_INDEX(j),XI,XE(1,NJ_XI(i)))
          ENDDO !j
          B(i)=-(DATA(i)-PXI(IBT(1,1,NB_XI(i)),
     '      IDO(1,1,0,NB_XI(i)),
     '      INP(1,1,NB_XI(i)),NB_XI(i),1,XI,XE(1,NJ_XI(i))))
        ENDDO !i
C news MPN 7Nov97
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(/'' count='',I2)') COUNT
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' A(i,j):'',8X,3D12.4,/(16X,3D12.4))')
     '      ((A(i,j),j=1,NITB),i=1,NITB)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' B:'',13X,3D12.4)') (B(i),i=1,NITB)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C newe
        ENDIF !DOP
        CALL DGETRF(NITB,NITB,A,NIM,IPIV,INFO)
        IF (INFO.NE.0) THEN
          COUNT=MAX_ITER+1
        ELSE
          CALL DGETRS(TRANS,NITB,1,A,NIM,IPIV,B,NIM,INFO)
          OLD_SUM=SUM
          SUM=0.0d0
          DO ni=1,NITB
            XI(ni)=XI(ni)-B(ni)
             SUM=SUM+(B(ni)**2.0d0)
          ENDDO !ni
          SUM=SUM**0.5d0
          CONVERG_TEST=DABS(SUM-OLD_SUM)/(DABS(OLD_SUM)+1.0d0)
C news MPN 7Nov97
          IF (NITB.EQ.3) THEN
            IF ((DABS(XI(1)).GT.50).OR.(DABS(XI(2)).GT.50)
     '        .OR.(DABS(XI(3)).GT.50)) THEN
              COUNT=MAX_ITER+1
            ENDIF
          ELSE
            IF ((DABS(XI(1)).GT.50).OR.(DABS(XI(2)).GT.50)) THEN
              COUNT=MAX_ITER+1
            ENDIF
          ENDIF
        ENDIF

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' XI increms:'',4X,3D12.4,'',  abs sum:'','
     '      //'D12.4)') (B(i),i=1,NITB),SUM
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' updated XI:'',4X,3D12.4)')
     '      (XI(i),i=1,NITB)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP
C newe
      ENDDO !tol & count

      IF(COUNT.LT.MAX_ITER) THEN !point found inside tolerance
        FITTED=.TRUE.
        nb=NBJ(1,ne)
        IF((IBT(1,1,nb).EQ. 3).AND.(IBT(3,1,nb).EQ.2)) THEN
! CS 18/9/97 If simplex with xi coordinates
          DO i=1,NITB
            IF(.NOT.((XI(i).GE.0.0d0).AND.(XI(i).LE.1.0d0)))THEN
              FITTED=.FALSE.
            ENDIF
          ENDDO
          IF(NITB.EQ.2) XI(3)=0.0d0
C**** CS 14/5/98 added tolerance
          IF(FITTED.AND.(XI(1)+XI(2)+XI(3).GT.(1.0d0+ZERO_TOL))) THEN
            FITTED=.FALSE.
          ENDIF
        ELSE
          DO i=1,NITB
C**** CS 14/5/98 added tolerance
C            IF(.NOT.((XI(i).GE.0.0d0).AND.(XI(i).LE.1.0d0)))THEN
            IF(.NOT.((XI(i).GE.-LOOSE_TOL).AND.(XI(i).LE.
     '        (1.0d0+LOOSE_TOL))))THEN
              FITTED=.FALSE.
            ENDIF
          ENDDO
        ENDIF
        IF(FITTED) THEN
          DO ni=1,NITB
            XID(ni)=XI(ni)
          ENDDO
          LD=ne
        ELSE
          LD=0
C Cannot have this warning as point may be in different element
C          WRITE(OP_STRING,'('' WARNING: No element found'')')
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !fitted

C LKC 9-AUG-2002 Adding USE_LOOSE_TOLERANCE for convergence problems
C with dexi nodes
      ELSEIF(USE_LOOSE_TOLERANCE.AND.
     '    CONVERG_TEST.LT.CONVERG_TOL*FACTOR_TOL)
     '    THEN
        FITTED=.TRUE.
        nb=NBJ(1,ne)

        FITTED=.TRUE.
        nb=NBJ(1,ne)
        IF((IBT(1,1,nb).EQ. 3).AND.(IBT(3,1,nb).EQ.2)) THEN
! CS 18/9/97 If simplex with xi coordinates
          DO i=1,NITB
            IF(.NOT.((XI(i).GE.0.0d0).AND.(XI(i).LE.1.0d0)))THEN
              FITTED=.FALSE.
            ENDIF
          ENDDO
          IF(NITB.EQ.2) XI(3)=0.0d0
C**** CS 14/5/98 added tolerance
          IF(FITTED.AND.(XI(1)+XI(2)+XI(3).GT.(1.0d0+ZERO_TOL))) THEN
            FITTED=.FALSE.
          ENDIF
        ELSE
          DO i=1,NITB
C**** CS 14/5/98 added tolerance
C            IF(.NOT.((XI(i).GE.0.0d0).AND.(XI(i).LE.1.0d0)))THEN
            IF(.NOT.((XI(i).GE.-LOOSE_TOL).AND.(XI(i).LE.
     '        (1.0d0+LOOSE_TOL))))THEN
              FITTED=.FALSE.
            ENDIF
          ENDDO
        ENDIF
        IF(FITTED) THEN
          DO ni=1,NITB
            XID(ni)=XI(ni)
          ENDDO
          LD=ne
          WRITE(OP_STRING,'('' WARNING: Point only converged'',
     '      '' to within '',F6.2,''x of CONVERG_TOL'')')FACTOR_TOL
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        ELSE
          LD=0
C Cannot have this warning as point may be in different element
C          WRITE(OP_STRING,'('' WARNING: No element found'')')
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !lose fitted
      ELSE
          LD=0
      ENDIF !count

      CALL EXITS('DEXI_POINT')
      RETURN
 9999 CALL ERRORS('DEXI_POINT',ERROR)
      CALL EXITS('DEXI_POINT')
      RETURN 1
      END


