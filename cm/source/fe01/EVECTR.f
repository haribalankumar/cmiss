      SUBROUTINE EVECTR(NDIM,A,EVAL,EVEC,ERROR,*)

C#### Subroutine: EVECTR
C###  Description:
C###    EVECTR returns normalized eigenvector EVEC(i),i=1,NDIM
C###    corresponding to eigenvalue EVAL of symmetric matrix A.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NDIM
      REAL*8 A(3,3),EVAL,EVEC(NDIM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,i1,i2,ICYCL(3,2),icol,irow
      REAL*8 AL,B1,B2,S,SUM,TOL,U(3,3),X(3)
      LOGICAL SINGLR

      DATA ICYCL/2,3,1,3,1,2/
      DATA TOL/1.0D-6/

      CALL ENTERS('EVECTR',*9999)
      IF(NDIM.EQ.2) THEN
        IF(DABS(A(1,2)).GT.TOL) THEN
          IF(DABS(A(1,1)-EVAL).GT.DABS(A(2,2)-EVAL)) THEN
            AL=DSQRT(A(1,2)**2+(A(1,1)-EVAL)**2)
            EVEC(1)=A(1,2)/AL
            EVEC(2)=(EVAL-A(1,1))/AL
          ELSE
            AL=DSQRT(A(1,2)**2+(A(2,2)-EVAL)**2)
            EVEC(1)=(EVAL-A(2,2))/AL
            EVEC(2)=A(1,2)/AL
          ENDIF
        ELSE IF(EVAL.EQ.A(1,1)) THEN
          EVEC(1)=1.0D0
          EVEC(2)=0.0D0
        ELSE IF(EVAL.EQ.A(2,2)) THEN
          EVEC(1)=0.0D0
          EVEC(2)=1.0D0
        ENDIF
C mh 8/01/96  why is this here?
C       EVEC(3)=0.0D0
      ELSE IF(NDIM.EQ.3) THEN
        IF(DABS(A(1,2)).LT.TOL.AND.DABS(A(1,3)).LT.TOL.AND.
     '    DABS(A(2,3)).LT.TOL) THEN
          GO TO 99
        ELSE
          DO irow=1,3
            DO icol=1,3
              U(irow,icol)=A(irow,icol)
            ENDDO
            U(irow,irow)=U(irow,irow)-EVAL
          ENDDO
          DO i=1,3
            X(i)=1.0D0
            i1=ICYCL(i,1)
            i2=ICYCL(i,2)
            B1=-1.0D0*U(i1,i)
            B2=-1.0D0*U(i2,i)
            CALL SOLV2(U(i1,i1),U(i1,i2),U(i2,i1),U(i2,i2),
     '         B1,B2,X(i1),X(i2),SINGLR,ERROR,*9999)
            SUM=0.0D0
            DO icol=1,3
              SUM=SUM+U(i,icol)*X(icol)
            ENDDO
            IF(DABS(SUM).LT.TOL.AND.(.NOT.SINGLR)) THEN
              S=0.0D0
              DO icol=1,3
                S=S+X(icol)**2
              ENDDO
              S=DSQRT(S)
              DO icol=1,3
                EVEC(icol)=X(icol)/S
              ENDDO
              GO TO 99
            ENDIF
          ENDDO
          WRITE(OP_STRING,'('' EVECTR: Eigenvectors undetermined'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          EVEC(1)=0.0D0
          EVEC(2)=0.0D0
          EVEC(3)=0.0D0
          GO TO 99
        ENDIF
      ENDIF
99    IF(DOP) THEN
        WRITE(OP_STRING,'('' EVEC: '',3E11.4)') (EVEC(i),i=1,NDIM)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('EVECTR')
      RETURN
 9999 CALL ERRORS('EVECTR',ERROR)
      CALL EXITS('EVECTR')
      RETURN 1
      END




