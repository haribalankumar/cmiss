      SUBROUTINE RUNGE_KUTTA_CURRENT_LINES(NEELEM,NLL,nr,DL,FIRST,ERROR,
     '  *)

C#### Subroutine: RUNGE_KUTTA_CURRENT_LINES
C###  Description:
C###    RUNGE_KUTTA_CURRENT_LINES uses a 4th order Runge-Kutta scheme to
C###    track  current lines. Currently has a constant step size.

C *** DPN 22 February 2001 - Appears to be unused, renaming from
C     RUNGE_KUTTA to avoid routine name conflict

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NLL(12,NEM),nr
      REAL*8 DL(3,NLM)
      LOGICAL FIRST
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nae,ne,nl
      REAL*8 AV_DIM,STEPSIZE
      LOGICAL STOP

      CALL ENTERS('RUNGE_KUTTA_CURRENT_LINES',*9999)

      STEPSIZE=0.5d0
      !Scale stepsize by average element dimension.
      !Note: A constant step size is probably preferable since we
      !want a good distribution of points to plot.
      ne=NEELEM(1,nr) !use the first element for size
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN !1d element
        nl=NLL(1,ne)
        IF(nl.GT.0) THEN
          STEPSIZE=STEPSIZE*DL(3,nl) !Scale by length of element side
        ELSE
          STEPSIZE=0.0d0
        ENDIF
      ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN !2d element
        AV_DIM=0.0d0
        DO nae=1,4
          nl=NLL(nae,ne)
          IF(nl.GT.0) AV_DIM=AV_DIM+DL(3,nl)
        ENDDO
        AV_DIM=AV_DIM/4.0d0
        CALL ASSERT(AV_DIM.GE.RDELTA,'  >>Zero size element ?',
     '              ERROR,*9999)
        STEPSIZE=STEPSIZE*AV_DIM
      ENDIF

      !useful for saying how far along current line we have gone.
      STOP=.FALSE.
      DO WHILE(.NOT.STOP)
        IF(.NOT.FIRST) THEN
C???
C??? LKC 7-OCT-97 ???  Since when has GRADPHI been archived
C???

C GRADPHI is archived
          ERROR='>>GRADPHI_I is archived, see RUNGE_KUTTA_CURRENT_LINES'
          GOTO 9999
C         CALL GRADPHI_I(IBT,IDO,INP,NBH,NBJ,
C    '      NEELEM,NGAP,NHE,NHP,NKE,NKH,NLL,NP_INTERFACE,
C    '      NPF,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,
C    '      NYNE,NYNP,CE,DET,DL,DRDN,
C    '      GRADPHI,PG,RAD,RG,SE,STEPSIZE,WG,XA,XE,XG1,
C    '      XIG,XN,XP,XPFP,XR,YP,ZA,ZE,ZF,ZP,STOP,ERROR,*9999)
        ELSE
          !If first is .true. then the starting point is the node
          !NPSTART.
          !Find initial gradient at node NPSTART
C GRADPHI is archived
          ERROR='>>GRADPHI_N is archived, see RUNGE_KUTTA_CURRENT_LINES'
          GOTO 9999
C         CALL GRADPHI_N(NKH,NP_INTERFACE,NPSTART,nr,NYNP,
C    '     GRADPHI,XG,XN,XP,YP,ERROR,*9999)
C         FIRST=.FALSE.
C         nv=1 ! temporary
C         DO nj=1,NJT
C           XPFP(nj)=XP(1,nv,nj,NPSTART)
C         ENDDO
        ENDIF
        !write out info at first point
C       WRITE(IO4,*)(XPFP(nj),nj=1,NJT),(GRADPHI(nj),nj=1,NJT)
        !Find k1
C       CALL FINDK(K1,GRADPHI,STEPSIZE,ERROR,*9999)
C       DO nj=1,NJT
C         XPFP(nj)=XPFP(nj)+K1(nj)/2.0d0
C       ENDDO
C       IF(DOP) THEN
CC$        call mp_setlock()
C         WRITE(OP_STRING,*)' S=',S
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' Finding K1'
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' GRADPHI(nj)=',(GRADPHI(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' K1(nj)=',(K1(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' XPFP(nj)=',(XPFP(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C       ENDIF
        !Find k2
C GRADPHI is archived
C       CALL GRADPHI_I(IBT,IDO,INP,NBH,NBJ,
C    '    NEELEM,NGAP,NHE,NHP,NKE,NKH,NLL,NP_INTERFACE,
C    '    NPF,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,
C    '    NYNE,NYNP,CE,DET,DL,DRDN,
C    '    GRADPHI,PG,RAD,RG,SE,STEPSIZE,WG,XA,XE,XG1,
C    '    XIG,XN,XP,XPFP,XR,YP,ZA,ZE,ZF,ZP,STOP,ERROR,*9999)
C       CALL FINDK(K2,GRADPHI,STEPSIZE,ERROR,*9999)
C       DO nj=1,NJT
C         XPFP(nj)=XPFP(nj)-K1(nj)/2.0d0+K2(nj)/2.0d0
C       ENDDO
C       IF(DOP) THEN
CC$        call mp_setlock()
C         WRITE(OP_STRING,*)' Finding K2'
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' GRADPHI(nj)=',(GRADPHI(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' K2(nj)=',(K2(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' XPFP(nj)=',(XPFP(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C       ENDIF
        !Find k3
C GRADPHI is archived
C       CALL GRADPHI_I(IBT,IDO,INP,NBH,NBJ,
C    '    NEELEM,NGAP,NHE,NHP,NKE,NKH,NLL,NP_INTERFACE,
C    '    NPF,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,
C    '    NYNE,NYNP,CE,DET,DL,DRDN,
C    '    GRADPHI,PG,RAD,RG,SE,STEPSIZE,WG,XA,XE,XG1,
C    '    XIG,XN,XP,XPFP,XR,YP,ZA,ZE,ZF,ZP,STOP,ERROR,*9999)
C       CALL FINDK(K3,GRADPHI,STEPSIZE,ERROR,*9999)
C       DO nj=1,NJT
C         XPFP(nj)=XPFP(nj)-K2(nj)/2.0d0+K3(nj)
C       ENDDO
C       IF(DOP) THEN
CC$        call mp_setlock()
C         WRITE(OP_STRING,*)' Finding K3'
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' GRADPHI(nj)=',(GRADPHI(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' K3(nj)=',(K3(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' XPFP(nj)=',(XPFP(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C       ENDIF
        !Find k4
C GRADPHI is archived
C       CALL GRADPHI_I(IBT,IDO,INP,NBH,NBJ,
C    '    NEELEM,NGAP,NHE,NHP,NKE,NKH,NLL,NP_INTERFACE,
C    '    NPF,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,
C    '    NYNE,NYNP,CE,DET,DL,DRDN,
C    '    GRADPHI,PG,RAD,RG,SE,STEPSIZE,WG,XA,XE,XG1,
C    '    XIG,XN,XP,XPFP,XR,YP,ZA,ZE,ZF,ZP,STOP,ERROR,*9999)
C       CALL FINDK(K4,GRADPHI,STEPSIZE,ERROR,*9999)
C       IF(DOP) THEN
CC$        call mp_setlock()
C         WRITE(OP_STRING,*)' Finding K4'
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' GRADPHI(nj)=',(GRADPHI(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C         WRITE(OP_STRING,*)' K4(nj)=',(K4(nj),nj=1,NJT)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C       ENDIF
C       DO nj=1,NJT
C         XPFP(nj)=XPFP(nj)-K3(nj) !Back to original point
C       ENDDO
        !Find new point
C       DO nj=1,NJT
C         XPFP(nj)=XPFP(nj)+(K1(nj)+2.0D0*K2(nj)+2.0D0*K3(nj)+
C    '             K4(nj))/6.0D0
C       ENDDO
C       S=S+STEPSIZE !Distance along current line.
C       IF(S.LE.4*STEPSIZE) THEN
C         STOP=.FALSE.
          !If we start from a nodal point then we need to ensure that
          !we travel away from the element before getting to another
          !point where we may want to stop (i.e. we should travel a
          !minimum distance before stopping).
C       ENDIF
C       IF(S.GT.MAXSTEP*STEPSIZE) THEN
C         STOP=.TRUE.
          !Don't want to travel forever
C       ENDIF
      ENDDO

      CALL EXITS('RUNGE_KUTTA_CURRENT_LINES')
      RETURN
9999  CALL ERRORS('RUNGE_KUTTA_CURRENT_LINES',ERROR)
      CALL EXITS('RUNGE_KUTTA_CURRENT_LINES')
      RETURN 1
      END


