      SUBROUTINE OSLO(j,NTKN,TAU,T,K,ALPHA)

C#### Subroutine: OSLO
C###  Description:
C###    OSLO implements Oslo algorithm for calculating discrete
C###    B-splines for refinement and rendering. This is algorithm 1 in
C###    Lyche and Morken, Siam J. Numer. Anal. 23:3:663-675.

      IMPLICIT NONE
!     Parameter List
      INTEGER J,K,NTKN
      REAL*8 ALPHA(*),TAU(*),T(*)
!     Local Variables
      INTEGER i,i1,i2,ih,il,ip,iu,iv,mu,MUDASH,n1,n2,RIGHT_MULT,
     '  LEFT_MULT
      REAL*8 AH(4,4),B,B1,D1,D2,TJ,XI(3)

C     initialise to zero
      n1=1
      n2=NTKN-k
      DO i=n1,n2
        ALPHA(i)=0.0D0
      ENDDO
C     check if any can be non-zero
      IF(T(j).LT.TAU(n1)) RETURN
      IF(T(j).EQ.TAU(n1))THEN
        IF(RIGHT_MULT(j,k,T).GT.RIGHT_MULT(n1,k,TAU))RETURN
      ENDIF
      IF(T(j+k).GT.TAU(n2+k)) RETURN
      IF(T(j+k).EQ.TAU(n2+k))THEN
        IF(LEFT_MULT(j+k,k,T).GT.LEFT_MULT(n2+k,k,TAU))RETURN
      ENDIF
C     calculate mu s.t. tau(mu)<=t(j)<tau(mu+1)
      DO mu=1,NTKN
        IF(TAU(mu).LE.T(j).AND.T(j).LT.TAU(mu+1)) GOTO 10
      ENDDO
 10   i=j+1
      MUDASH=mu
      DO WHILE (T(i).EQ.TAU(MUDASH).AND.i.LT.j+k)
        i=i+1
        MUDASH=MUDASH-1
      ENDDO
      ih=MUDASH+1
      iv=0
      DO ip=1,k-1
        IF(T(j+ip).EQ.TAU(ih))THEN
          ih=ih+1
        ELSE
          iv=iv+1
          XI(iv)=T(j+ip)
        ENDIF
      ENDDO
      AH(k,1)=1.0D0
      DO ip=1,iv
        B1=0.0D0
        TJ=XI(ip)
        IF(ip.GE.MUDASH)THEN
          B1=(TJ-TAU(n1))*AH(1+k-MUDASH,ip)/(TAU(ip+k-iv)-TAU(n1))
        ENDIF
        il=MAX(n1+1,MUDASH-ip+1)
        iu=MIN(MUDASH,n2+iv-ip)
        DO i=il,IU
          D1=TJ-TAU(i)
          D2=TAU(i+ip+k-iv-1)-TJ
          B=AH(i+k-MUDASH,ip)/(D1+D2)
          AH(i+k-MUDASH-1,ip+1)=D2*B+B1
          B1=D1*B
        ENDDO
        AH(iu+k-MUDASH,ip+1)=B1
        IF(iu.LT.MUDASH)THEN
          AH(iu+k-MUDASH,ip+1)=B1+(TAU(n2+k)-TJ)*AH(iu+k-MUDASH+1,ip)
     '      /(TAU(n2+k)-TAU(IU+1))
        ENDIF
      ENDDO
      i1=MAX(MUDASH-iv,n1)
      i2=MIN(MUDASH,n2)
      DO i=i1,i2
        ALPHA(i)=AH(i+k-MUDASH,1+iv)
      ENDDO

      RETURN
      END


