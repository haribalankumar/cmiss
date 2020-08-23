      SUBROUTINE MATMOD(NCYCLE,NPROB,FINISH,m,n,nb,ne,NKA,ns,NSCL,
     '  NNAME,A,HA,KA,BL,BU,ASCALE,HS,NAME1,NAME2,X,PI,RC,Z,NWCORE,
     '  IUSER,USER)

C#### Subroutine: MATMOD
C###  Description:
C###    MATMOD routine (called by MINOS).

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER NKA,NNAME,IUSER(*),KA(NKA),m,n,NAME1(NNAME),
     '  NAME2(NNAME),nb,NCYCLE,ne,NPROB,ns,NSCL,NWCORE
      INTEGER*4 HA(ne),HS(nb)
      REAL*8 A(ne),ASCALE(NSCL),BL(nb),BU(nb),PI(M),RC(nb),X(nb),
     '  USER(*),Z(NWCORE)
      LOGICAL FINISH

      CALL ENTERS('MATMOD',*9999)

C**** This subroutine does nothing

      CALL EXITS('MATMOD')
      RETURN
 9999 CALL ERRORIN('MATMOD')
      CALL ERRORIN(' ') !flush call stack
      CALL EXITS('MATMOD')
      RETURN
      END


