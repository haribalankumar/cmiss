      SUBROUTINE UPFGMULTIPLY(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &  ERROR,*)

C#### Subroutine: UPFGMULTIPLY
C###  Description:
C###    FG=F*G at values and derivatives


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NUMVALUES,FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES)
      REAL*8 F(NKM,NVM,NUMVALUES),FG(NKM,NVM,NUMVALUES),
     '  G(NKM,NVM,NUMVALUES)
      CHARACTER ERROR*(*)
      LOGICAL MISMATCH
!     Local Variables
      INTEGER I(1:8,0:16),nk,np,nv,t
!      REAL*8
!      LOGICAL
!      CHARACTER

      DATA I/1,2,2,4,2,4,4,8,
     '       1,2,3,4,5,6,7,8,
     '       1,1,1,1,1,1,1,1,
     '       0,1,1,2,1,2,3,4,
     '       0,2,3,3,5,5,5,5,
     '       0,0,0,3,0,5,5,6,
     '       0,0,0,2,0,2,3,3,
     '       0,0,0,1,0,1,1,2,
     '       0,0,0,4,0,6,7,7,
     '       0,0,0,0,0,0,0,7,
     '       0,0,0,0,0,0,0,2,
     '       0,0,0,0,0,0,0,3,
     '       0,0,0,0,0,0,0,6,
     '       0,0,0,0,0,0,0,5,
     '       0,0,0,0,0,0,0,4,
     '       0,0,0,0,0,0,0,1,
     '       0,0,0,0,0,0,0,8/

      CALL ENTERS('UPFGMULTIPLY',*9999)

      DO np=1,FGNK(1,0)
        DO nv=1,min(FGNV(1,np),FGNV(2,np))
          DO nk=1,min(FGNK(1,np),FGNK(2,np))
            FG(nk,nv,np)=0.0D0
            DO t=0,I(nk,0)-1
              FG(nk,nv,np)=FG(nk,nv,np)+F(I(nk,t*2+1),nv,np)*
     '         G(I(nk,t*2+2),nv,np)
            ENDDO
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

      CALL EXITS('UPFGMULTIPLY')
      RETURN
 9999 CALL ERRORS('UPFGMULTIPLY',ERROR)
      CALL EXITS('UPFGMULTIPLY')
      RETURN 1
      END

