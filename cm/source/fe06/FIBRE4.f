      SUBROUTINE FIBRE4(INDEX,IBT,IDO,INP,IW,NAN,NBJ,nr,XE,XG,XI,
     '  ERROR,*)

C#### Subroutine: FIBRE4
C###  Description:
C###    FIBRE4 draws fitted fibre helix.
C**** Note: DIRECTION (=+1 or -1) is set initially in DRFIBR

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fibr00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER INDEX,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IW,NAN(NIM,NAM,NBFM),NBJ(NJM),nr
      REAL*8 XE(NSM,NJM),XG(NJM,NUM),XI(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ni,nj,njj,ns,NT_POINTS
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  DXI,dXi_1,dXi_2,dXi_3,
     '  G(3,3),POINTS(3,500),PXI,X(3),XI_1,XI_2,XI_3
      LOGICAL CONTINUE

      CALL ENTERS('FIBRE4',*9999)

      IF(DOP) THEN
        DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,njj,nr)
          nb=NBJ(nj)
          WRITE(OP_STRING,'('' XE: '',8E11.3)')
     '      (XE(ns,nj),ns=1,NST(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

! Compute initial plotting point
      DO nj=1,3
        nb=NBJ(nj)
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
      ENDDO
      CALL XZ(ITYP10(1),X,POINTS(1,1))
      NT_POINTS=1

! Loop until
      CONTINUE=.TRUE.
      DO WHILE (CONTINUE)
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Starting  Xi: '',3F8.5)')
     '      (XI(ni),ni=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
!   Compute normalized Xi coord base vectors
        CALL BASE_XI(IBT,IDO,INP,NBJ,G,XE,XI,ERROR,*9999)
        IF(DOP) THEN
          write(OP_STRING,'('' g1: '',3E12.3)') G(1,1),G(1,2),G(1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          write(OP_STRING,'('' g2: '',3E12.3)') G(2,1),G(2,2),G(2,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          write(OP_STRING,'('' g3: '',3E12.3)') G(3,1),G(3,2),G(3,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
!   Compute fibre vector
        CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '    A_VECTOR,B_VECTOR,C_VECTOR,XE,XG,XI,.TRUE.,ERROR,*9999)
!   Compute projection of fibre vector onto Xi coords
        dXi_1=G(1,1)*A_VECTOR(1)+G(1,2)*A_VECTOR(2)
     '    +G(1,3)*A_VECTOR(3)
        dXi_2=G(2,1)*A_VECTOR(1)+G(2,2)*A_VECTOR(2)
     '    +G(2,3)*A_VECTOR(3)
CMPN 2/10/92 - temporary to keep fibres at correct Xi(3) value
        dXi_3=0.d0
C will need to replace these statements for nonzero imbrecation angle
C        dXi_3=G(3,1)*A_VECTOR(1)+G(3,2)*A_VECTOR(2)
C     '    +G(3,3)*A_VECTOR(3)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' dXi_1='',E12.3,'' dXi_2='',E12.3,'
     '      //''' dXi_3='',E12.3)') dXi_1,dXi_2,dXi_3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
!   Update position of current Xi point
        DXI=0.02D0/DSQRT(dXi_1*dXi_1+dXi_2*dXi_2+dXi_3*dXi_3)
        XI_1=XI(1)+DIRECTION*dXi_1*DXI
        XI_2=XI(2)+DIRECTION*dXi_2*DXI
        XI_3=XI(3)+DIRECTION*dXi_3*DXI
        IF(DOP) THEN
          WRITE(OP_STRING,'('' Direction='',F4.0,'' xi_2='',F9.7)')
     '      DIRECTION,Xi_2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(XI_2.GT.1.d0) THEN !change direction
          DIRECTION=-DIRECTION
          IF(DOP) THEN
            write(OP_STRING,'('' Direction changed'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          XI_1=XI(1)+DIRECTION*dXi_1*DXI
          XI_2=XI(2)+DIRECTION*dXi_2*DXI
          XI_3=XI(3)+DIRECTION*dXi_3*DXI
          IF(DOP) THEN
            write(OP_STRING,'('' Direction='',F4.0)') DIRECTION
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        XI(1)=XI_1
        XI(2)=XI_2
        XI(3)=XI_3
        IF(DOP) THEN
          WRITE(OP_STRING,'('' Finishing Xi: '',3F8.5)')
     '      (XI(ni),ni=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
!   Update position of current plotting point
        DO nj=1,NJT
          nb=NBJ(nj)
          X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '      XE(1,nj))
        ENDDO
        NT_POINTS=NT_POINTS+1
        IF(DOP) THEN
          WRITE(OP_STRING,'('' no points='',I4)') NT_POINTS
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL XZ(ITYP10(1),X,POINTS(1,NT_POINTS))
        IF(XI(1).LT.0.d0.OR.XI(1).GT.1.d0.OR
     '    .XI(2).LT.0.d0.OR
     '    .XI(3).LT.0.d0.OR.XI(3).GT.1.d0) CONTINUE=.FALSE.
        IF(NT_POINTS.GE.200) THEN
          WRITE(OP_STRING,*) ' Too many points'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CONTINUE=.FALSE.
        ENDIF
      ENDDO

! Plot trajectory
      CALL POLYLINE(INDEX,IW,NT_POINTS,POINTS,ERROR,*9999)

      CALL EXITS('FIBRE4')
      RETURN
 9999 CALL ERRORS('FIBRE4',ERROR)
      CALL EXITS('FIBRE4')
      RETURN 1
      END


