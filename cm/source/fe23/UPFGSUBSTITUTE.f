      SUBROUTINE UPFGSUBSTITUTE(FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &  ERROR,*)

C#### Subroutine: UPFGSUBSTITUTE
C###  Description:
C###    FG=G for values and derivatives


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NUMVALUES,FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES)
      REAL*8 FG(NKM,NVM,NUMVALUES),G(NKM,NVM,NUMVALUES)
      CHARACTER ERROR*(*)
      LOGICAL MISMATCH
!     Local Variables
      INTEGER nk,np,nv
!      REAL*8
!      LOGICAL
!      CHARACTER

      CALL ENTERS('UPFGSUBSTITUTE',*9999)

      DO np=1,FGNK(1,0)
        DO nv=1,min(FGNV(1,np),FGNV(2,np))
          DO nk=1,min(FGNK(1,np),FGNK(2,np))
            FG(nk,nv,np)=G(nk,nv,np)
          ENDDO
        ENDDO
        ! Check inconsistency between nv for F and G
        DO nv=min(FGNV(1,np),FGNV(2,np))+1,FGNV(1,np)
          DO nk=1,min(FGNK(1,np),FGNK(2,np))
            FG(nk,nv,np)=0.0D0
            MISMATCH=.TRUE.
          ENDDO
        ENDDO
        ! Check inconsistency between nk for F and G
        DO nk=min(FGNK(1,np),FGNK(2,np))+1,FGNK(1,np)
          DO nv=1,min(FGNV(1,np),FGNV(2,np))
            FG(nk,nv,np)=0.0D0
            MISMATCH=.TRUE.
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('UPFGSUBSTITUTE')
      RETURN
 9999 CALL ERRORS('UPFGSUBSTITUTE',ERROR)
      CALL EXITS('UPFGSUBSTITUTE')
      RETURN 1
      END


