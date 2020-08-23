      SUBROUTINE MOMENTS(NRLIST,SMOM,XP,ZD,ERROR,*)

C#### Subroutine: MOMENTS
C###  Description:
C###    MOMENTS calculates the second moments of area for the
C###    current node and or data set.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      !Parameter List
      INTEGER NRLIST(0:NRM)
      REAL*8 SMOM(2,6),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
      !Local Variables
      INTEGER nd,nj,nomom,np,nr,notyp

      CALL ENTERS('MOMENTS',*9999)

      nr=NRLIST(1)
      DO notyp=1,2
        DO nomom=1,NJT+2*NJT-3
          SMOM(notyp,nomom)=0.0d0
        ENDDO
      ENDDO

      DO np=1,NPT(nr)
        DO nj=1,NJT
          SMOM(1,nj)=SMOM(1,nj)+XP(1,1,nj,np)**2
        ENDDO
        SMOM(1,NJT+1)=SMOM(1,NJT+1)+XP(1,1,1,np)*XP(1,1,2,np)
        IF(NJT.GT.2) THEN
          SMOM(1,NJT+2)=SMOM(1,NJT+2)+XP(1,1,2,np)*XP(1,1,3,np)
          SMOM(1,NJT+3)=SMOM(1,NJT+3)+XP(1,1,3,np)*XP(1,1,1,np)
        ENDIF
      ENDDO

      DO nd=1,NDT
        DO nj=1,NJT
          SMOM(2,nj)=SMOM(2,nj)+ZD(nj,nd)**2
        ENDDO
        SMOM(2,NJT+1)=SMOM(2,NJT+1)+ZD(1,nd)*ZD(2,nd)
        IF(NJT.GT.2) THEN
          SMOM(2,NJT+2)=SMOM(2,NJT+2)+ZD(2,nd)*ZD(3,nd)
          SMOM(2,NJT+3)=SMOM(2,NJT+3)+ZD(3,nd)*ZD(1,nd)
        ENDIF
      ENDDO

      DO nomom=1,NJT+2*NJT-3
        SMOM(1,nomom)= SMOM(1,nomom)/NPT(nr)
        SMOM(2,nomom)= SMOM(2,nomom)/NDT
      ENDDO

      CALL EXITS('MOMENTS')
      RETURN
 9999 CALL ERRORS('MOMENTS',ERROR)
      CALL EXITS('MOMENTS')
      RETURN 1
      END

