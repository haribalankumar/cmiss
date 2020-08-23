      SUBROUTINE GAUSSLOG(IBT,IDO,INP,nb,ISIGN,N,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSSLOG
C###  Description:
C###    <HTML>
C###    GAUSSLOG returns the N Gauss points (XIG) and weights (WG) and
C###    evaluates the basis functions at those Gauss points for a
C###    Gauss quadrature scheme with a logarithmic weight function.
C###    The limits are 0 to 1. The function form that is to be indicated
C###    depends on ISIGN. The values of ISIGN are:
C###    <UL>
C###    <LI> 0 : ln(x)f(x)
C###    <LI> 1 : -ln(x)f(x)
C###    <LI> 2 : ln(1-x)f(1-x)
C###    <LI> 3 : -ln(1-x)f(1-x)
C###    </UL>
C###    The algorithms used in this subroutine were derrived from
C###    Numerical Recipies in Fortran (2nd edition), section 4.5
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
C
      INTEGER NMAX
      PARAMETER (NMAX=64)
C
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,ISIGN,N
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k,ng,nk,nn,ns,nu,INFO
      REAL*8 A(NMAX),ALPHA(2*NMAX-1),AMU0,ANU(2*NMAX),B(NMAX),
     '  BETA(2*NMAX-1),P,
     '  SIGMA(2*NMAX+1,2*NMAX+1),SUM,Z(NMAX,NMAX),PSI1,WORK(2*NMAX-2)

      CALL ENTERS('GAUSSLOG',*9999)

      CALL ASSERT(N.GT.0,'>>Number of Gauss points must be positive',
     '  ERROR,*9999)
      CALL ASSERT(N.LE.NMAX,'>>Increase NMAX in GAUSSLOG',ERROR,*9999)
      CALL ASSERT(ISIGN.GE.0.AND.ISIGN.LE.3,'>>Invalid ISIGN',
     '  ERROR,*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Number of Gauss points, N = '',I5)') N
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISIGN = '',I2)') ISIGN
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C Find the alpha, beta an nu coefficients for the recurrance relation

      ALPHA(1)=0.5d0
      BETA(1)=1.0d0
      ANU(1)=1.0d0
      ANU(2)=-0.25d0
      DO i=2,2*N-1
        ALPHA(i)=0.5d0
        BETA(i)=0.25d0/(4.0d0-1.0d0/DBLE(i-1)**2)
        ANU(i+1)=-ANU(i)*DBLE(i)*DBLE(i-1)/(2.0d0*DBLE(i+1)*DBLE(2*i-1))
      ENDDO !i

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Reccurrance coefficients: alpha, beta '
     '    //'and nu'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,2*N-1
          WRITE(OP_STRING,'('' i='',I3,3D12.4)') i,ALPHA(i),BETA(i),
     '      ANU(i)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'('' i='',I3,24X,D12.4)') 2*N,ANU(2*N)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C Compute the a and b coeffients of the recurrence relation for the
C monic orthogonal polynomials by Wheeler's algorthim

      DO i=3,2*N
        SIGMA(1,i)=0.0d0
      ENDDO !i
      DO i=2,2*N+1
        SIGMA(2,i)=ANU(i-1)
      ENDDO !i
      A(1)=ALPHA(1)+ANU(2)/ANU(1)
      B(1)=0.0d0
      DO i=3,N+1
        DO j=i,2*N-i+3
          SIGMA(i,j)=SIGMA(i-1,j+1)+(ALPHA(j-1)-A(i-2))*SIGMA(i-1,j)-
     '      B(i-2)*SIGMA(i-2,j)+BETA(j-1)*SIGMA(i-1,j-1)
        ENDDO !j
        A(i-1)=ALPHA(i-1)+SIGMA(i,i+1)/SIGMA(i,i)-
     '    SIGMA(i-1,i)/SIGMA(i-1,i-1)
        B(i-1)=SIGMA(i,i)/SIGMA(i-1,i-1)
      ENDDO !i

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Polynomial coefficients: a, b'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,N
          WRITE(OP_STRING,'('' i='',I3,2D12.4)') i,A(i),B(i)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

C cpb 29/4/97 Using LAPACK routine to find eigenvalues and eigenvectors
CC Compute the tridiagonal Jacobi matrix
C
C      DO i=1,N
C        IF(i.NE.1) B(i)=DSQRT(B(i))
C        DO j=1,N
C          IF(i.EQ.j) THEN
C            Z(i,j)=1.0d0
C          ELSE
C            Z(i,j)=0.0d0
C          ENDIF
C        ENDDO !j
C      ENDDO !i
C
CC Find the eigenvalues and eigenvectors of the tridiagonal matrix with
CC A(i) on the diagonal and B(i) on the sub-diagonal by means of a QL
CC algorithm
C
C      DO i=2,N
C        B(i-1)=B(i)
C      ENDDO !i
C      B(N)=0.0d0
C 1    DO l=1,N
C        ITER=0
C        DO m=L,N-1
C          DD=DABS(A(m))+DABS(A(m+1))
C          IF(DABS(B(m))+DD.EQ.DD) GOTO 2
C        ENDDO !m
C        m=N
C 2      IF(m.NE.l) THEN
C          IF(ITER.EQ.MAXITS) THEN
C            ERROR='>>Too many QL iterations'
C            GOTO 9999
C          ENDIF
C          ITER=ITER+1
C          G=(A(l+1)-A(l))/(2.0d0*B(l))
C          R=PYTHAG(G,1.0d0)
C          G=A(m)-A(l)+B(l)/(G+SIGN(R,G))
C          S=1.0d0
C          C=1.0d0
C          P=0.0d0
C          DO i=M-1,l,-1
C            F=S*B(i)
C            E=C*B(i)
C            R=PYTHAG(F,G)
C            B(i+1)=R
C            IF(R.EQ.0.0d0) THEN
C              A(i+1)=A(i+1)-P
C              B(m)=0.0d0
C              GOTO 1
C            ENDIF
C            S=F/R
C            C=G/R
C            G=A(i+1)-P
C            R=(A(i)-G)*S+2.0d0*C*E
C            P=S*R
C            A(i+1)=G+P
C            G=C*R-E
C            DO k=1,N
C              F=Z(k,i+1)
C              Z(k,i+1)=S*Z(k,i)+C*F
C              Z(k,i)=C*Z(k,i)-S*F
C            ENDDO !k
C          ENDDO !i
C          A(l)=A(l)-P
C          B(l)=G
C          B(m)=0.0d0
C          GOTO 1
C        ENDIF
C      ENDDO !l
C

C     Find the eigenvalues and eigenvectors
      DO i=2,N
        B(i)=DSQRT(B(i))
      ENDDO !i
      CALL DSTEV('V',N,A,B(2),Z,NMAX,WORK,INFO)

      IF(INFO.NE.0) THEN
        WRITE(OP_STRING,'('' >>INFO='',I5,'' in DSTEV'')') INFO
        GOTO 9999
      ENDIF

C Sort the eigenvalues and eigenvectors into descending order

      DO i=1,N-1
        k=i
        P=A(i)
        DO j=i+1,N
          IF(A(j).GE.P) THEN
            k=j
            P=A(j)
          ENDIF
        ENDDO !j
        IF(k.NE.i) THEN
          A(k)=A(i)
          A(i)=P
          DO j=1,N
            P=Z(j,i)
            Z(j,i)=Z(j,k)
            Z(j,k)=P
          ENDDO !j
        ENDIF
      ENDDO !i

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Eigenvalues:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,N
          WRITE(OP_STRING,'('' i='',I3,D12.4)') i,A(i)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'('' First component of Eigenvectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,N
          WRITE(OP_STRING,'('' i='',I3,D12.4)') i,Z(1,i)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

C Find the positions and weights from the eigenvalues and eigenvectors

      IF(ISIGN.EQ.0.OR.ISIGN.EQ.3) THEN
        AMU0=-1.0d0 !the integral of the +'VE weight func
      ELSE IF(ISIGN.EQ.1.OR.ISIGN.EQ.2) THEN
        AMU0=1.0d0 !the integral of the -'VE weight func
      ENDIF

      DO i=1,N
        IF(ISIGN.EQ.0.OR.ISIGN.EQ.1) THEN
          XIG(1,i)=A(i)
        ELSE IF(ISIGN.EQ.2.OR.ISIGN.EQ.3) THEN
          XIG(1,i)=1.0d0-A(i)
        ENDIF
        WG(i)=AMU0*Z(1,i)**2
      ENDDO !i

      DO ng=1,N
        ns=0
        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            ns=ns+1
            DO nu=1,NUT(nb)
              PG(ns,nu,ng)=PSI1(IBT,IDO,INP,nb,nu,nk,nn,XIG(1,ng))
            ENDDO !nu
          ENDDO !nk
        ENDDO !nn
      ENDDO !ng


      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' The '',I2,'' Gauss point locations '
     '    //'are:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,N
          WRITE(OP_STRING,'('' XIG('',I2,'')='',D17.10)') i,XIG(1,i)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !i
        WRITE(OP_STRING,'('' The '',I2,'' Gauss point weights '
     '    //'are:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,N
          WRITE(OP_STRING,'('' WG('',I2,'')='',D17.10)') i,WG(i)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !i
        SUM=0.0d0
        DO i=1,N
          SUM=SUM+WG(i)
        ENDDO !i
        WRITE(OP_STRING,'('' The sum of the weights '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' is '',D17.10,'' and should be '',D17.10)')
     '    SUM,AMU0
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        SUM=0.0d0
        IF(ISIGN.EQ.0.OR.ISIGN.EQ.1) THEN
          DO i=1,N
            SUM=SUM+WG(i)/(1.0d0+XIG(1,i))**2
          ENDDO !i
        ELSE IF(ISIGN.EQ.2.OR.ISIGN.EQ.3) THEN
          DO i=1,N
            SUM=SUM+WG(i)/(2.0d0-XIG(1,i))**2
          ENDDO !i
        ENDIF
        IF(ISIGN.EQ.0) THEN
          WRITE(OP_STRING,'('' Integral of ln(x)/(1+x)**2 from '
     '      //'0 to 1'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' is '',D17.10,'' and should be '','
     '      //'D17.10)') SUM,-LOG(2.0d0)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(ISIGN.EQ.1) THEN
          WRITE(OP_STRING,'('' Integral of -ln(x)/(1+x)**2 from '
     '      //'0 to 1'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' is '',D17.10,'' and should be '','
     '      //'D17.10)') SUM,LOG(2.0d0)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(ISIGN.EQ.2) THEN
          WRITE(OP_STRING,'('' Integral of ln(1-x)/(2-x)**2 from '
     '      //'0 to 1'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' is '',D17.10,'' and should be '','
     '      //'D17.10)') SUM,LOG(2.0d0)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(ISIGN.EQ.3) THEN
          WRITE(OP_STRING,'('' Integral of -ln(1-x)/(2-x)**2 from '
     '      //'0 to 1'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' is '',D17.10,'' and should be '','
     '      //'D17.10)') SUM,-LOG(2.0d0)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('GAUSSLOG')
      RETURN
 9999 CALL ERRORS('GAUSSLOG',ERROR)
      CALL EXITS('GAUSSLOG')
      RETURN 1
      END

C MLB 18/3/97 *** ARCHIVED ***
C
C      SUBROUTINE GAUSSPWB(ALIM,BLIM,N,NNSING,WEIGHT,ABSCIS,ICALL,
C     '  ERROR,*)
C
