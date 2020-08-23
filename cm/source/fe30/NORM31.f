      SUBROUTINE NORM31(NDIM,nq,NXQ,DXDXIQ,DXDXIQ2,XNLOCAL,ERROR,*)

C#### Subroutine: NORM31
C###  Description:
C###    NORM31 finds the normal vector XNLOCAL(nj) in the cartesian
C###    reference frame. Given a boundary grid point this
C###    routine will return a unit outward normal.
C**** Written by Martin Buist, February 1999.

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NDIM,nq,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),XNLOCAL(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER METHOD,nj
      REAL*8 VECTOR1(3),VECTOR2(3),VLENGTH

      CALL ENTERS('NORM31',*9999)

      XNLOCAL(1)=0.0d0
      XNLOCAL(2)=0.0d0
      XNLOCAL(3)=0.0d0

      IF(NDIM.EQ.1) THEN
        IF(NXQ(-1,1,nq).EQ.0) THEN
          XNLOCAL(1)=-1.0d0
        ELSEIF(NXQ(1,1,nq).EQ.0) THEN
          XNLOCAL(1)=1.0d0
        ELSE
          ERROR='Grid point is not a boundary point'
          GOTO 9999
        ENDIF
      ELSEIF(NDIM.EQ.2) THEN
        METHOD=3

        IF(METHOD.EQ.1) THEN
          IF(NXQ(-1,1,nq).EQ.0) THEN
            XNLOCAL(1)=-DXDXIQ(2,2,nq)
            XNLOCAL(2)=DXDXIQ(1,2,nq)
          ELSEIF(NXQ(1,1,nq).EQ.0) THEN
            XNLOCAL(1)=DXDXIQ(2,2,nq)
            XNLOCAL(2)=-DXDXIQ(1,2,nq)
          ELSEIF(NXQ(-2,1,nq).EQ.0) THEN
            XNLOCAL(1)=-DXDXIQ(2,1,nq)
            XNLOCAL(2)=DXDXIQ(1,1,nq)
          ELSEIF(NXQ(2,1,nq).EQ.0) THEN
            XNLOCAL(1)=DXDXIQ(2,1,nq)
            XNLOCAL(2)=-DXDXIQ(1,1,nq)
          ELSE
            !Not quite this
            XNLOCAL(1)=-1.0d0
            XNLOCAL(2)=0.0d0
          ENDIF
          VLENGTH=DSQRT((XNLOCAL(1)**2)+(XNLOCAL(2)**2))
          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(1)=XNLOCAL(1)/VLENGTH
          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(2)=XNLOCAL(2)/VLENGTH

        ELSE IF(METHOD.EQ.2) THEN

C Please leave
C          IF(NXQ(-1,1,nq).EQ.0) THEN
C            nq1=NXQ(-2,1,nq)
C            nq2=NXQ(2,1,nq)
C          ELSEIF(NXQ(1,1,nq).EQ.0) THEN
C            nq1=NXQ(2,1,nq)
C            nq2=NXQ(-2,1,nq)
C          ELSEIF(NXQ(-2,1,nq).EQ.0) THEN
C            nq1=NXQ(1,1,nq)
C            nq2=NXQ(-1,1,nq)
C          ELSEIF(NXQ(2,1,nq).EQ.0) THEN
C            nq1=NXQ(-1,1,nq)
C            nq2=NXQ(1,1,nq)
C          ELSE
C            ERROR='>>Not a boundary point'
C            GOTO 9999
C          ENDIF
C
C          DO nj=1,NJT
C            TN1(nj)=XQ(nj,nq)-XQ(nj,nq1)
C            TN2(nj)=XQ(nj,nq2)-XQ(nj,nq)
C          ENDDO !nj
C
C          TN1W=0.0d0
C          TN2W=0.0d0
C          DO nj=1,NJT
C            TN1W=TN1W+TN1(nj)**2.0d0
C            TN2W=TN2W+TN2(nj)**2.0d0
C          ENDDO !nj
C          TN1W=DSQRT(TN1W)
C          TN2W=DSQRT(TN2W)
C
C          DO nj=1,NJT
C            TNLOCAL(nj)=(TN1(nj)/TN2W)+(TN2(nj)/TN1W)
C          ENDDO
C
C          XNLOCAL(1)=TNLOCAL(2)
C          XNLOCAL(2)=-TNLOCAL(1)
C
C          VLENGTH=DSQRT((XNLOCAL(1)**2)+(XNLOCAL(2)**2))
C          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(1)=XNLOCAL(1)/VLENGTH
C          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(2)=XNLOCAL(2)/VLENGTH

        ELSE IF(METHOD.EQ.3) THEN

          IF(NXQ(-1,1,nq).EQ.0) THEN
            XNLOCAL(1)=-DXDXIQ2(2,2,nq)
            XNLOCAL(2)=DXDXIQ2(1,2,nq)
          ELSEIF(NXQ(1,1,nq).EQ.0) THEN
            XNLOCAL(1)=DXDXIQ2(2,2,nq)
            XNLOCAL(2)=-DXDXIQ2(1,2,nq)
          ELSEIF(NXQ(-2,1,nq).EQ.0) THEN
            XNLOCAL(1)=-DXDXIQ2(2,1,nq)
            XNLOCAL(2)=DXDXIQ2(1,1,nq)
          ELSEIF(NXQ(2,1,nq).EQ.0) THEN
            XNLOCAL(1)=DXDXIQ2(2,1,nq)
            XNLOCAL(2)=-DXDXIQ2(1,1,nq)
          ELSE
            !Not quite this
            XNLOCAL(1)=-1.0d0
            XNLOCAL(2)=0.0d0
          ENDIF
          VLENGTH=DSQRT((XNLOCAL(1)**2)+(XNLOCAL(2)**2))
          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(1)=XNLOCAL(1)/VLENGTH
          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(2)=XNLOCAL(2)/VLENGTH

C Values now calculated in upgrid and stored in DXDXIQ2
C          ne=NENQ(1,nq)
C          DO ni2=1,NQXI(2,NQS(ne))
C            DO ni1=1,NQXI(1,NQS(ne))
C              neq=ni1+((ni2-1)*NQXI(1,NQS(ne)))
C              nq1=NQNE(ne,neq)
C              IF(nq1.EQ.nq) THEN
C                XI(1)=DBLE(ni1-1)/DBLE(NQXI(1,NQS(ne))-1)
C                XI(2)=DBLE(ni2-1)/DBLE(NQXI(2,NQS(ne))-1)
C                XI(3)=0.0d0
C              ENDIF
C            ENDDO
C          ENDDO
C
C          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
C         '      NVJE(1,1,1,ne),SE(1,1,ne),XE,XP,ERROR,*9999)
C          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,XG,XI,
C     '      ERROR,*9999)
C
C          XNLOCAL(1)=-XG(1,4)
C          XNLOCAL(2)=XG(1,2)
C
C          VLENGTH=DSQRT((XNLOCAL(1)**2)+(XNLOCAL(2)**2))
C          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(1)=XNLOCAL(1)/VLENGTH
C          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(2)=XNLOCAL(2)/VLENGTH

        ENDIF

      ELSEIF(NDIM.EQ.3) THEN
        IF(NXQ(-1,1,nq).EQ.0) THEN
          DO nj=1,NJT
            VECTOR1(nj)=-DXDXIQ(nj,2,nq)
            VECTOR2(nj)=DXDXIQ(nj,3,nq)
          ENDDO
        ELSEIF(NXQ(1,1,nq).EQ.0) THEN
          DO nj=1,NJT
            VECTOR1(nj)=DXDXIQ(nj,2,nq)
            VECTOR2(nj)=DXDXIQ(nj,3,nq)
          ENDDO
        ELSEIF(NXQ(-2,1,nq).EQ.0) THEN
          DO nj=1,NJT
            VECTOR1(nj)=-DXDXIQ(nj,3,nq)
            VECTOR2(nj)=DXDXIQ(nj,1,nq)
          ENDDO
        ELSEIF(NXQ(2,1,nq).EQ.0) THEN
          DO nj=1,NJT
            VECTOR1(nj)=DXDXIQ(nj,3,nq)
            VECTOR2(nj)=DXDXIQ(nj,1,nq)
          ENDDO
        ELSEIF(NXQ(-3,1,nq).EQ.0) THEN
          DO nj=1,NJT
            VECTOR1(nj)=-DXDXIQ(nj,1,nq)
            VECTOR2(nj)=DXDXIQ(nj,2,nq)
          ENDDO
        ELSEIF(NXQ(3,1,nq).EQ.0) THEN
          DO nj=1,NJT
            VECTOR1(nj)=DXDXIQ(nj,1,nq)
            VECTOR2(nj)=DXDXIQ(nj,2,nq)
          ENDDO
        ELSE
          WRITE(OP_STRING,'('' nq='',I9)') nq
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ERROR='Grid point is not a boundary point'
          GOTO 9999
        ENDIF
        XNLOCAL(1)=VECTOR1(2)*VECTOR2(3)-VECTOR1(3)*VECTOR2(2)
        XNLOCAL(2)=VECTOR1(3)*VECTOR2(1)-VECTOR1(1)*VECTOR2(3)
        XNLOCAL(3)=VECTOR1(1)*VECTOR2(2)-VECTOR1(2)*VECTOR2(1)
        VLENGTH=DSQRT((XNLOCAL(1)**2)+(XNLOCAL(2)**2)+(XNLOCAL(3)**2))
        IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(1)=XNLOCAL(1)/VLENGTH
        IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(2)=XNLOCAL(2)/VLENGTH
        IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(3)=XNLOCAL(3)/VLENGTH
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Normal at '',I6, ''is '',3F8.4)') nq,
     '    XNLOCAL(1),XNLOCAL(2),XNLOCAL(3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Vector length '',F12.8)') VLENGTH
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('NORM31')
      RETURN
 9999 CALL ERRORS('NORM31',ERROR)
      CALL EXITS('NORM31')
      RETURN 1
      END



