      SUBROUTINE XIXL(IBT,IDO,INP,ISTEP,iw,NBJ,XI,XE,XL,ERROR,*)

C#### Subroutine: XIXL
C###  Description:
C###    XIXL returns mapped world coords XL(istep,nj) at given Xi(ni)
C###    coordinates.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISTEP,iw,NBJ(NJM)
      REAL*8 XI(4),XE(NSM,NJM),XL(3,500)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nj
      REAL*8 PXI,X(3)

      CALL ENTERS('XIXL',*9999)
      DO nj=1,NJT
        nb=NBJ(nj)
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '    nb,1,XI,XE(1,nj))
      ENDDO
      IF(iw.NE.4) THEN
        CALL XZ(ITYP10(1),X,XL(1,ISTEP))
      ELSE
        DO nj=1,NJT
          XL(nj,ISTEP)=X(nj)
        ENDDO
      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' XL(nj,'',I3,''): '',3E10.3)')
     '    ISTEP,XL(1,ISTEP),XL(2,ISTEP),XL(3,ISTEP)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('XIXL')
      RETURN
 9999 CALL ERRORS('XIXL',ERROR)
      CALL EXITS('XIXL')
      RETURN 1
      END


