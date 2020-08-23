      SUBROUTINE XINORM(IBT,IDO,INP,in1,in2,NBJ,XI1,RADIUS,X,XE,XI,
     '  XIS,XNORM,XP_IB,KEEP,ERROR,*)

C#### Subroutine: XINORM
C###  Description:
C###    XINORM locates the XI coordinates of a point on a surface
C###    (in1-in2 face) that is normal to an internal point
C###    (coordinates in X).


      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  in1,in2,NBJ(NJM),XI1
      REAL*8 RADIUS,X(3),XE(NSM,NJM),XI(3),XIS,XNORM(3),XP_IB(3)
      LOGICAL KEEP
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER IT,ITMAX,nj
      REAL*8 ERR,K,KIN,KOUT,USER_TOL,XNORM_OLD(3),Z_LEN
      LOGICAL FOUND,FOUND_E,INELEM

      CALL ENTERS('XINORM',*9999)

      K=RADIUS
      KOUT=K
      KIN=0.d0
      USER_TOL=LOOSE_TOL
      FOUND=.FALSE.
      FOUND_E=.TRUE.
      INELEM=.TRUE.
      ITMAX=200
      CALL CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ,Z_LEN,USER_TOL,XE,
     '  XI,XP_IB,FOUND_E,INELEM,ERROR,*9999) !get XI coords for IB
      XI(XI1)=XIS !set to boundary face
      CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ,X,XE,XI,XNORM,ERROR,*9999)
      !gets normal to surface point
      DO WHILE(.NOT.FOUND) !until normal point located
        DO nj=1,3
          X(nj)=XP_IB(nj)+K*XNORM(nj) !moving away from IB point
          XNORM_OLD(nj)=XNORM(nj)
        ENDDO !nj
        FOUND_E=.TRUE.
        CALL CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ,Z_LEN,USER_TOL*1.d2,XE,
     '    XI,X,FOUND_E,INELEM,ERROR,*9999) !XI coords for new X
        IF(FOUND_E.AND.DABS(XI(XI1)-XIS).LE.LOOSE_TOL)THEN
          CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ,X,XE,
     '      XI,XNORM,ERROR,*9999) !XNORM=normal to new X
          ERR=0.d0 !check error for normals
          DO nj=1,3
            ERR=ERR+(XNORM(nj)-XNORM_OLD(nj))**2.d0
            XNORM_OLD(nj)=XNORM(nj)
          ENDDO !nj
          ERR=DSQRT(ERR)
          IF(ERR.GT.LOOSE_TOL)THEN
            K=RADIUS
            KOUT=K
            KIN=0.d0
            CALL CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ,Z_LEN,USER_TOL,XE,
     '        XI,XP_IB,FOUND_E,INELEM,ERROR,*9999) !get XI coords for IB
            XI(XI1)=XIS !set to boundary face
            ! resetting initial XI position
          ELSE
            FOUND=.TRUE.
          ENDIF
        ELSE !not surface point yet, increment K
          IF(FOUND_E)THEN !K to point within element
            KIN=K
          ELSE IF(.NOT.FOUND_E)THEN !K to point outside element
            KOUT=K
          ENDIF
          K=0.5d0*(KIN+KOUT)
          CALL CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ,Z_LEN,USER_TOL,XE,
     '      XI,XP_IB,FOUND_E,INELEM,ERROR,*9999) !get XI coords for IB
          XI(XI1)=XIS !set to boundary face
 ! resetting initial XI position
        ENDIF
        IF(K.LT.LOOSE_TOL)THEN !can't find a surface point normal to IB
          FOUND=.TRUE.
          KEEP=.FALSE.
        ELSE IF(DABS(KIN-KOUT).LT.ZERO_TOL)THEN !write warning
          FOUND=.TRUE.
          KEEP=.TRUE.
          WRITE(OP_STRING,'('' CONVERGENCE WARNING '')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !WHILE

      CALL EXITS('XINORM')
      RETURN
 9999 CALL ERRORS('XINORM',ERROR)
      CALL EXITS('XINORM')
      RETURN 1
      END



