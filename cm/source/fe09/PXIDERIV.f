
      REAL*8 FUNCTION PXIDERIV(IBT,IDO,INP,nb,nu,XI,XE)

C#### Function: PXIDERIV
C###  Type: REAL*8
C###  Description:
C###    PXIDERIV interpolates nodal array XE at XI with basis functions as
C###    originally defined (jtyp5=1) or transformed to monomial form
C###    (jtyp5=2).

C#### Variable: nu
C###  Type: INTEGER
C###  Set_up: PXIDERIV
C###  Description:
C###    <HTML>
C###    nu is the calculated derivative loop variable.
C###    <PRE>
C###    nu=1 function
C###    nu=2 derivative wrt Xi1
C###    nu=3 derivative wrt Xi1^2
C###    nu=4 derivative wrt Xi2
C###    nu=5 derivative wrt Xi2^2
C###    nu=6 derivative wrt Xi1 Xi2
C###    nu=7 derivative wrt Xi3
C###    nu=8 derivative wrt Xi3^2
C###    nu=9 derivative wrt Xi1 Xi3
C###    nu=10 derivative wrt Xi2 Xi3
C###    nu=11 derivative wrt Xi1 XI2 Xi3
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nu
      REAL*8 XE(NSM),XI(3)
!     Local Variables
      INTEGER i1,i2,i3,IP1,IP2,IP3,ni,nk,nn,ns
      REAL*8 COSHX,CSS,D,DES,PSI1,PSI2,PSI2_HERMITE,PSI2_XI,
     '  PSI5,RAD,SINHX,SS,SUM1,SUM2,THETA,XL(4)
      LOGICAL SECTOR
      SAVE RAD,COSHX,SINHX

        PXIDERIV = 0.0d0
            ns=0
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
               ns=ns+1
               PXIDERIV=PXIDERIV+PSI1(IDO,INP,nb,nu,nk,nn,XI)*XE(ns)
!              write(*,*) nu,nn,nk,XE(ns),PSI1(IDO,INP,nb,nu,nk,nn,XI)
              ENDDO
            ENDDO

           write(*,*) '     ------------'

      RETURN
      END

