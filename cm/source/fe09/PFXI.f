      REAL*8 FUNCTION PFXI(IBT,IDO,INP,NAN,nb,nu,XE,XI)

C#### Function: PFXI
C###  Type: REAL*8
C###  Description:
C###    PFXI interpolates nodal array XE at XI.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),NAN(NIM,NAM),
     '  nb,nu
      REAL*8  XE(NSM),XI(3)
!     Local Variables
      INTEGER na,nn,nk,ns
      REAL*8 PSI1,PSI8

      PFXI=0.0D0
      ns=0
      DO nn=1,NNT(nb)
        DO nk=1,NKT(nn,nb)
          ns=ns+1
          PFXI=PFXI+PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)*XE(ns)
        ENDDO
      ENDDO
      DO na=1,NAT(nb)
        ns=ns+1
        PFXI=PFXI+PSI8(NAN(1,na),nb,nu,XI)*XE(ns)
      ENDDO

      RETURN
      END


