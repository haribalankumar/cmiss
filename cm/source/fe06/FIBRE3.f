      SUBROUTINE FIBRE3(INDEX,IBT,IDO,INP,IW,NAN,NBJ,nr,XE,XG,ERROR,*)

C#### Subroutine: FIBRE3
C###  Description:
C###    FIBRE3 draws fitted fibre sheet.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fibr00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER INDEX,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IW,NAN(NIM,NAM,NBFM),NBJ(NJM),nr
      REAL*8 XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ni,nj,njj,ns
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  AXES(3,2),dXi_2,dXi_3,G(3,3),PXI,SCALE,X(3),XI(3)
      LOGICAL CONTINUE

      CALL ENTERS('FIBRE3',*9999)

      SCALE=0.05D0*DBLE(DIAG)
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
      XI(1)=DXI1   !is initial Xi_1 postion
      XI(2)=DXI2   !is initial Xi_2 postion
      XI(3)=0.d0  !is initial Xi_3 postion

! Loop until
      CONTINUE=.TRUE.
      DO WHILE (CONTINUE)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' Xi: '',3F8.5)') (XI(ni),ni=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
!   Compute start of fibre vector
        DO nj=1,3
          nb=NBJ(nj)
          X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '      XE(1,nj))
        ENDDO
        CALL XZ(ITYP10(nr),X,AXES(1,1))
!   Compute fibre vector and sheet vector
        CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '    A_VECTOR,B_VECTOR,C_VECTOR,XE,XG,XI,.TRUE.,ERROR,*9999)
!   Plot fibre vector
        DO nj=1,NJT
          AXES(nj,2)=AXES(nj,1)+SCALE*A_VECTOR(nj)
        ENDDO
        CALL POLYLINE(INDEX,IW,2,AXES,ERROR,*9999)
!   Compute normalized Xi coord base vectors
        CALL BASE_XI(IBT,IDO,INP,NBJ,G,XE,XI,ERROR,*9999)
        write(OP_STRING,'('' g2: '',3E12.3)') G(2,1),G(2,2),G(2,3)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        write(OP_STRING,'('' g3: '',3E12.3)') G(3,1),G(3,2),G(3,3)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!   Compute projection of sheet vector onto Xi coords
        dXi_2=G(2,1)*B_VECTOR(1)+G(2,2)*B_VECTOR(2)
     '    +G(2,3)*B_VECTOR(3)
        dXi_3=G(3,1)*B_VECTOR(1)+G(3,2)*B_VECTOR(2)
     '    +G(3,3)*B_VECTOR(3)
        write(OP_STRING,'('' dXi_2='',E12.3,'' dXi_3='',E12.3)')
     '    dXi_2,dXi_3
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!   Update position of current Xi point
        IF(DABS(dXi_3).GT.0.d0) THEN
          XI(2)=XI(2)+dXi_2/dXi_3*0.1D0
        ELSE
          XI(2)=XI(2)
        ENDIF
        XI(3)=XI(3)+0.1D0
        IF(DOP) THEN
          WRITE(OP_STRING,'('' Xi: '',3F8.5)') (XI(ni),ni=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(XI(3).GT.1.d0) CONTINUE=.FALSE.
      ENDDO

      CALL EXITS('FIBRE3')
      RETURN
 9999 CALL ERRORS('FIBRE3',ERROR)
      CALL EXITS('FIBRE3')
      RETURN 1
      END


