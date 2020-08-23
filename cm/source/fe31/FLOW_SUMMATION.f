      SUBROUTINE FLOW_SUMMATION(NBJ,NEELEM,NPNE,NVJE,NVJP,NXI,XP,
     &  ERROR,*)

C#### Subroutine: FLOW_SUMMATION
C###  Description:
C###    

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ne2,noelem,noelem2,np,np1,np2,nv,nv1,nv2
      REAL*8 flow

      CALL ENTERS('FLOW_SUMMATION',*9999)

      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        np=NPNE(2,1,ne)
        IF(NXI(1,0,ne).EQ.0)THEN !terminal, flow already set
!          XP(1,1,nj_flow,np)=XP(2,1,nj_flow,np)
C ARC 20-1-12 THIS IS NEEDED FOR THE ITERATIVE STEP - ABOVE JUST FORCES IT TO CONVERGE WHEN ITS NOT 
          XP(2,1,nj_flow,np)=(XP(2,1,nj_flow,np)
     &      +XP(1,1,nj_flow,np))/2.d0
          XP(1,1,nj_flow,np)=(XP(2,1,nj_flow,np)
     &      +XP(1,1,nj_flow,np))/2.d0
        ELSE
          flow=0.d0
          DO noelem2=1,NXI(1,0,ne) !for each daughter branch
            ne2=NXI(1,noelem2,ne) !the daughter element number
            nb=NBJ(1,ne2)
            np2=NPNE(2,nb,ne2)
            flow=flow+XP(2,1,nj_flow,np2) !sum daughter flows
            XP(2,noelem2+1,nj_flow,np)=XP(2,1,nj_flow,np2)
            XP(1,noelem2+1,nj_flow,np)=XP(1,1,nj_flow,np2)
          ENDDO
          XP(2,1,nj_flow,np)=flow
          XP(1,1,nj_flow,np)=flow
        ENDIF
      ENDDO
      ne=NEELEM(1)
      np=NPNE(1,nb,ne)
      np2=NPNE(2,nb,ne)
      XP(2,1,nj_flow,np)=XP(2,1,nj_flow,np2) 
      XP(1,1,nj_flow,np)=XP(1,1,nj_flow,np2) 

C.....Store the normalised flows in nj_sV field      
c      XP(1,1,nj_sV,1)=1.d0
c      DO noelem=1,NEELEM(0)
c        ne=NEELEM(noelem)
c        nb=NBJ(1,ne)
c        np1=NPNE(1,1,ne)
c        np2=NPNE(2,1,ne)
c        nv1=NVJE(1,nb,nj_flow,ne)
c        nv2=NVJE(2,nb,nj_flow,ne)
c        XP(1,nv1,nj_sV,np1)=XP(2,nv1,nj_flow,np1)/XP(2,1,nj_flow,1)
c        XP(1,nv2,nj_sV,np2)=XP(2,nv2,nj_flow,np2)/XP(2,1,nj_flow,1)
c      ENDDO

      CALL EXITS('FLOW_SUMMATION')
      RETURN
 9999 CALL ERRORS('FLOW_SUMMATION',ERROR)
      CALL EXITS('FLOW_SUMMATION')
      RETURN 1
      END
