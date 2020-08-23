      SUBROUTINE ZGTG53GAUSFROMGRID(ne,ng,NGLIST,NQLIST,NQNE,YG,YQS,
     &  RET_ERROR,*)

C#### Subroutine: ZGTG53GAUSFROMGRID
C###  Description:
C###    <html><pre>ZGTG53GAUSFROMGRID updates YG(nj,ng,ne) array from
C###         cellular grid point
C###    </pre></html>

      IMPLICIT NONE
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
!     Parameter List
      INTEGER ne,ng,NGLIST(0:NGM),NQLIST(0:NQM),NQNE(NEQM,NQEM)
      REAL*8 YG(NIYGM,NGM,NEM),YQS(NIQSM,NQM)
      CHARACTER RET_ERROR*(*)
!     Local Variables
      INTEGER nig,niqq,niqs

      CALL ENTERS('ZGTG53GAUSFROMGRID',*9999)
C currently assuming Tdevij are always in YQS 2,3,4,5,6,7
      NQLIST(0)=6
      NQLIST(1)=2
      NQLIST(2)=3
      NQLIST(3)=4
      NQLIST(4)=5
      NQLIST(5)=6
      NQLIST(6)=7
C Also assuming that stresses are filled into YG array 1,2,3,4,5,6
      NGLIST(0)=6
      NGLIST(1)=1
      NGLIST(2)=2
      NGLIST(3)=3
      NGLIST(4)=4
      NGLIST(5)=5
      NGLIST(6)=6
      DO niqq=1,NQLIST(0)
        niqs=NQLIST(niqq)
        nig=NGLIST(niqq)
        YG(nig,ng,ne)=0.0d0
        YG(nig,ng,ne)=YQS(niqs,NQNE(ne,ng))
      ENDDO !niqq (niqs)


      CALL EXITS('ZGTG53GAUSFROMGRID')
      RETURN
 9999 CALL ERRORS('ZGTG53GAUSFROMGRID',RET_ERROR)
      CALL EXITS('ZGTG53GAUSFROMGRID')
      RETURN 1
      END


