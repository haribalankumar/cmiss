      SUBROUTINE UPFGDIVIDE(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &  ERROR,*)

C#### Subroutine: UPFGDIVIDE
C###  Description:
C###    FG=F/G at values and derivatives


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NUMVALUES,FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES)
      REAL*8 F(NKM,NVM,NUMVALUES),FG(NKM,NVM,NUMVALUES),
     '  G(NKM,NVM,NUMVALUES)
      CHARACTER ERROR*(*)
      LOGICAL MISMATCH
!     Local Variables
      INTEGER nk,np,nv
!      INTEGER I(1:8,0:16),t
      REAL*8 DF,DG,F1,F2,G1,P1,P2,Q1,Q2,Q3
!      LOGICAL
!      CHARACTER

      CALL ENTERS('UPFGDIVIDE',*9999)
      nv=1
      DO np=1,FGNK(1,0)
        DO nv=1,min(FGNV(1,np),FGNV(2,np))
          DO nk=1,min(FGNK(1,np),FGNK(2,np))
            IF(nk.EQ.1)THEN
              FG(1,nv,np)=F(1,nv,np)/G(1,nv,np)
            ELSEIF(nk.EQ.2)THEN
              FG(2,nv,np)=(F(2,nv,np)*G(1,nv,np)-F(1,nv,np)
     '                   *G(2,nv,np))/(G(1,nv,np)**2)
            ELSEIF(nk.EQ.3)THEN
              FG(3,nv,np)=(F(3,nv,np)*G(1,nv,np)-F(1,nv,np)
     '                   *G(3,nv,np))/(G(1,nv,np)**2)
            ELSEIF(nk.EQ.4)THEN
              FG(4,nv,np)=((F(4,nv,np)*G(1,nv,np)+F(2,nv,np)
     '                   *G(3,nv,np)-F(3,nv,np)*G(2,nv,np)
     '                  -F(1,nv,np)*G(4,nv,np))*G(1,nv,np)**2
     '                  -(F(2,nv,np)*G(1,nv,np)
     '                  -F(1,nv,np)*G(2,nv,np))*2*G(1,nv,np)
     '                  *G(3,nv,np))/G(1,nv,np)**4
            ELSEIF(nk.EQ.5)THEN
              FG(5,nv,np)=(F(5,nv,np)*G(1,nv,np)-F(1,nv,np)
     '                  *G(5,nv,np))/(G(1,nv,np)**2)
            ELSEIF(nk.EQ.6)THEN
              FG(6,nv,np)=((F(6,nv,np)*G(1,nv,np)+F(2,nv,np)*G(5,nv,np)
     '                  -F(5,nv,np)*G(2,nv,np)-F(1,nv,np)*G(6,nv,np))
     '                  *G(1,nv,np)**2-(F(2,nv,np)*G(1,nv,np)
     '                  -F(1,nv,np)*G(2,nv,np))*2*G(1,nv,np)
     '                  *G(5,nv,np))/G(1,nv,np)**4
            ELSEIF(nk.EQ.7)THEN
              FG(7,nv,np)=((F(7,nv,np)*G(1,nv,np)+F(3,nv,np)*G(5,nv,np)
     '                  -F(5,nv,np)*G(3,nv,np)-F(1,nv,np)*G(7,nv,np))
     '                  *G(1,nv,np)**2-(F(3,nv,np)*G(1,nv,np)
     '                  -F(1,nv,np)*G(3,nv,np))*2*G(1,nv,np)
     '                  *G(5,nv,np))/G(1,nv,np)**4
            ELSEIF(nk.EQ.8)THEN
              F1 =  (F(4,nv,np)*G(1,nv,np)+F(2,nv,np)*G(3,nv,np)
     '              -F(3,nv,np)*G(2,nv,np)-F(1,nv,np)*G(4,nv,np))
     '              *G(1,nv,np)**2
              F2 = -(F(2,nv,np)*G(1,nv,np)-F(1,nv,np)*G(2,nv,np))
     '              *2*G(1,nv,np)*G(3,nv,np)
              P1 = (F(8,nv,np)*G(1,nv,np)+F(4,nv,np)*G(5,nv,np)
     '              +F(6,nv,np)*G(3,nv,np)+F(2,nv,np)*G(7,nv,np)
     '              -F(7,nv,np)*G(2,nv,np)-F(3,nv,np)*G(6,nv,np)
     '              -F(5,nv,np)*G(4,nv,np)-F(1,nv,np)*G(8,nv,np))
     '              *G(1,nv,np)**2
              P2 = (F(4,nv,np)*G(1,nv,np)+F(2,nv,np)*G(3,nv,np)
     '             -F(3,nv,np)*G(2,nv,np)-F(1,nv,np)*G(4,nv,np))
     '             *2*G(1,nv,np)*G(5,nv,np)
              Q1 = -(F(6,nv,np)*G(1,nv,np)+F(2,nv,np)*G(5,nv,np)
     '              -F(5,nv,np)*G(2,nv,np)-F(1,nv,np)*G(6,nv,np))
     '              *2*G(1,nv,np)*G(3,nv,np)
              Q2 = -(F(2,nv,np)*G(1,nv,np)-F(1,nv,np)*G(2,nv,np))
     '              *2*G(5,nv,np)*G(3,nv,np)
              Q3 = -(F(2,nv,np)*G(1,nv,np)-F(1,nv,np)*G(2,nv,np))
     '              *2*G(1,nv,np)*G(7,nv,np)
              DF = P1+P2+Q1+Q2+Q3
              G1 = G(1,nv,np)**4
              DG = 4*(G(1,nv,np)**3)*G(5,nv,np)
            
              FG(8,nv,np)= (DF*G1-(F1+F2)*DG)/G(1,nv,np)**8
c              FG(8,nv,np)=((F(8,nv,np)*G(1,nv,np)+F(4,nv,np)*G(5,nv,np)
c     '                  +F(6,nv,np)*G(3,nv,np)+F(2,nv,np)*G(7,nv,np)
c     '                  -F(7,nv,np)*G(2,nv,np)-F(3,nv,np)*G(6,nv,np)
c     '                  -F(5,nv,np)*G(4,nv,np)-F(1,nv,np)*G(8,nv,np))
c     '                  *G(1,nv,np)**2+(F(4,nv,np)*G(1,nv,np)+F(2,nv,np)
c     '                  *G(3,nv,np)-F(3,nv,np)*G(2,nv,np)-F(1,nv,np)
c     '                  *G(4,nv,np))*2*G(1,nv,np)*G(5,nv,np)
c     '                  -((F(2,nv,np)*G(1,nv,np)-F(1,nv,np)*G(2,nv,np))
c     '                  *(2*G(5,nv,np)*G(3,nv,np)+2*G(1,nv,np)
c     '                  *G(7,nv,np)))*G(1,nv,np)**4
c     '                  -((F(4,nv,np)*G(1,nv,np)+F(2,nv,np)*G(3,nv,np)
c     '                  -F(3,nv,np)*G(2,nv,np)-F(1,nv,np)*G(4,nv,np))
c     '                  *G(1,nv,np)**2-(F(2,nv,np)*G(1,nv,np)-F(1,nv,np)
c     '                  *G(2,nv,np))*2*G(1,nv,np)*G(3,nv,np))*4
c     '                  *G(1,nv,np)**3*G(5,nv,np))/(G(1,nv,np)**8)
            ELSE
              ERROR='>> Derivative for NK not evaluated. Update code'
              GOTO 9999
            ENDIF
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

      CALL EXITS('UPFGDIVIDE')
      RETURN
 9999 CALL ERRORS('UPFGDIVIDE',ERROR)
      CALL EXITS('UPFGDIVIDE')
      RETURN 1
      END

