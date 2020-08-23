      SUBROUTINE FLOW_AVERAGE(np1,np2,steps,vol,vratio,CE,
     &  dPl,dt,alpha,Qinit,XAB,XP,ERROR,*)

C#### Subroutine: FLOW_AVERAGE
C###  Description:
C###    

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER np1,np2,steps
      REAL*8 alpha,CE(NMM),dPl,dt,
     &  vratio,vol,XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 beta,C,dPadt,dPe,Q,Qinit,R,Ra

      CALL ENTERS('FLOW_AVERAGE',*9999)

      Ra=0.d0 !acinar resistance
      R = CE(nm_R) + Ra ! Pa.s.mm-3
      C = CE(nm_C) ! mm^3/Pa
      beta = dPl*CE(nm_dPl)/dt ! == dPpl/dt (-ve for inspiration)
      Q = C*(alpha-beta)+(Qinit-C*(alpha-beta))*DEXP(-dt/(C*R))

c      XP(2,1,nj_flow,np2)=Q

      XAB(7,np2)=XAB(6,np2)
      XAB(6,np2)=XP(2,1,nj_flow,np2)

! include previous two iterations
      XP(2,1,nj_flow,np2)=0.75d0*XAB(7,np2)+0.25d0*(Q+XAB(6,np2))*0.5d0

      CALL EXITS('FLOW_AVERAGE')
      RETURN
 9999 CALL ERRORS('FLOW_AVERAGE',ERROR)
      CALL EXITS('FLOW_AVERAGE')
      RETURN 1
      END
