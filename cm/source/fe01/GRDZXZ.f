      SUBROUTINE GRDZXZ(INDEX,ERROR,*)

C#### Subroutine: GRDZXZ
C###  Description:
C###    GRDZXZ draws grid lines of constant z in x,z coordinates.

      IMPLICIT NONE
      INCLUDE 'graf00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER INDEX
      REAL*8 DIFF,P(3,3),RFROMC,ZINCR,ZINCRT,ZZ
      CHARACTER CHAR*9

      CALL ENTERS('GRDZXZ',*9999)
      ZINCR=DBLE(ZMAX-ZMIN)/10.0D0
C      CHAR=CFROMR(ZINCR,'(E8.1)')
      WRITE(CHAR,'(E8.1)') ZINCR
      ZINCRT=RFROMC(CHAR)
      DIFF=DBLE(XMAX-XMIN)
      P(1,2)=DBLE(XMAX)-0.05D0*DIFF
      DO ZZ=0.0D0,DBLE(ZMAX)-ZINCRT,ZINCRT
        P(1,1)=DBLE(XMIN)+0.05D0*DIFF
        P(3,1)=ZZ
        P(3,2)=ZZ
        CALL POLYLINE(INDEX,1,2,P,ERROR,*9999)
C        CHAR=CFROMR(ZZ,'(E9.2)')
        WRITE(CHAR,'(E9.2)') ZZ
        P(1,1)=P(1,1)-0.025D0*DIFF
        CALL TEXT(1,1,CHAR(3:5),P,ERROR,*9999)
      ENDDO
C      CHAR=CFROMR(10.0D0*ZINCR,'(E8.1)')
      WRITE(CHAR,'(E8.1)') 10.0D0*ZINCR
      P(1,1)=DBLE(XMIN)+0.015D0*DIFF
      P(2,1)=P(3,1)+0.5D0*ZINCRT
      CALL TEXT(1,1,'*'//CHAR(5:8),P,ERROR,*9999)
      DO ZZ=-ZINCRT,DBLE(ZMIN)+ZINCRT,-ZINCRT
        P(1,1)=DBLE(XMIN)+0.05D0*DIFF
        P(3,1)=ZZ
        P(3,2)=ZZ
        CALL POLYLINE(INDEX,1,2,P,ERROR,*9999)
C        CHAR=CFROMR(ZZ,'(E9.2)')
        WRITE(CHAR,'(E9.2)') ZZ
        P(1,1)=P(1,1)-0.025D0*DIFF
        CALL TEXT(1,1,'-'//CHAR(3:5),P,ERROR,*9999)
      ENDDO

      CALL EXITS('GRDZXZ')
      RETURN
 9999 CALL ERRORS('GRDZXZ',ERROR)
      CALL EXITS('GRDZXZ')
      RETURN 1
      END


