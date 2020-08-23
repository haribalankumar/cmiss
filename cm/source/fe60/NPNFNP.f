      SUBROUTINE NPNFNP(IBT,IDO,INP,in1,in2,NBJ,nf,NFNP,np,XE,XI,
     '  XPNP,ADD_NF,ERROR,*)

C#### Subroutine: NPNFNP
C###  Description:
C###    NPNFNP sets up NFNP and XPNP.  That is, the faces surrounding
C###    node np are stored in NFNP, and a normal projection from the
C###    node is stored in XPNP.


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  in1,in2,NBJ(NJM),nf,NFNP(0:10,NP_R_M),np
      REAL*8 XE(NSM,NJM),XI(3),XPNP(3,NP_R_M)
      LOGICAL ADD_NF
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ncount,nj
      REAL*8 X(3),XNORM(3)

      CALL ENTERS('NPNFNP',*9999)

      IF(ADD_NF)THEN
        ncount=NFNP(0,np)+1
        NFNP(0,np)=ncount
        NFNP(ncount,np)=nf
        CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ,X,XE,XI,XNORM,ERROR,*9999)
        DO nj=1,3
          XPNP(nj,np)=XPNP(nj,np)+X(nj)-XNORM(nj)
        ENDDO !nj
      ENDIF !ADD_NF
      ADD_NF=.TRUE.

      CALL EXITS('NPNFNP')
      RETURN
 9999 CALL ERRORS('NPNFNP',ERROR)
      CALL EXITS('NPNFNP')
      RETURN 1
      END


