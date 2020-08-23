
      SUBROUTINE SET_FD3D_EQUATIONS(A,LDA,N,B,NX,NY,NZ,MTYPE,ISC_A,
     '  ISR_A,NZA,TRANSIENT,ERROR,*)

C#### Subroutine: SET_FD3D_EQUATIONS
C###  Description:
C###    Set a system of 3D Finite Difference equations to test the solvers.
C###  Written by Stuart Norris 02/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NX,NY,NZ,MTYPE,ISC_A(*),ISR_A(*),NZA
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
      LOGICAL TRANSIENT
!     Local Variables
      INTEGER I,J,K,L,M
      INTEGER NNX,NNY,NNZ
      REAL*8 X(NX),Y(NY),Z(NZ),DX(NX),DY(NY),DZ(NZ)
      REAL*8 AE,AW,AN,AS,AT,AB,DP,DE,DW,DN,DS,DT,DB
      REAL*8 VE,VW,VN,VS,VT,VB,PI,AREAX,AREAY,AREAZ,VOLUME,RHO,DELT

      DATA VE,VN,VT,VW,VS,VB,RHO / 3*1.0D0,3*-1.0D0,1.0D0 /


      NNX=NX-2
      NNY=NY-2
      NNZ=NZ-2

C     Set the mesh
      ! Regular mesh
      IF(MTYPE.EQ.1) THEN
        DO I=1,NX
          X(I)=DBLE(I-1)/DBLE(NX-1)
        ENDDO
        DO J=1,NY
          Y(J)=DBLE(J-1)/DBLE(NY-1)
        ENDDO
        DO K=1,NZ
          Z(K)=DBLE(K-1)/DBLE(NZ-1)
        ENDDO

      ! Sin mesh
      ELSE IF(MTYPE.EQ.2) THEN
        PI=4.0D0*ATAN(1.0D0)
        DO I=1,NX
          X(I)=0.5D0*( 1.0D0 - COS(PI*DBLE(I-1)/DBLE(NX-1)) )
        ENDDO
        DO J=1,NY
          Y(J)=0.5D0*( 1.0D0 - COS(PI*DBLE(J-1)/DBLE(NY-1)) )
        ENDDO
        DO K=1,NZ
          Z(K)=0.5D0*( 1.0D0 - COS(PI*DBLE(K-1)/DBLE(NZ-1)) )
        ENDDO

      ! Error
      ELSE
        ERROR='>>Unknown mesh type'
        GOTO 9999
      ENDIF

      DO I=2,NX
        DX(I)=X(I)-X(I-1)
      ENDDO
      DO J=2,NY
        DY(J)=Y(J)-Y(J-1)
      ENDDO
      DO K=2,NZ
        DZ(K)=Z(K)-Z(K-1)
      ENDDO
      DELT=2.0D0*MIN(X(NX)/DBLE(NX),Y(NY)/DBLE(NY),Z(NZ)/DBLE(NZ))**2

C     Set the system
      M=0
      L=1
      DO K=2,NZ-1
        DO J=2,NY-1
          DO I=2,NX-1
            M=M+1
            ISR_A(M)=L

            AREAX=0.25D0*(DY(J+1)+DY(J))*(DZ(K+1)+DZ(K))
            AREAY=0.25D0*(DX(I+1)+DX(I))*(DZ(K+1)+DZ(K))
            AREAZ=0.25D0*(DX(I+1)+DX(I))*(DY(J+1)+DY(J))
            VOLUME=SQRT(AREAX*AREAY*AREAZ)

            AE=AREAX/DX(I+1)
            AW=AREAX/DX(I)
            AN=AREAY/DY(J+1)
            AS=AREAY/DY(J)
            AT=AREAZ/DZ(K+1)
            AB=AREAZ/DZ(K)

            DE=-AE
            DW=-AW
            DN=-AN
            DS=-AS
            DT=-AT
            DB=-AB
            DP=(AE+AW+AN+AS+AT+AB)

            IF(TRANSIENT) DP=DP+VOLUME*RHO/DELT
            B(M)=0.0D0

            IF(K.GT.2) THEN
              A(L)=DB
              ISC_A(L)=M-NNX*NNY
              L=L+1
            ELSE
              B(M)=B(M)-DB*VB
            ENDIF

            IF(J.GT.2) THEN
              A(L)=DS
              ISC_A(L)=M-NNX
              L=L+1
            ELSE
              B(M)=B(M)-DS*VS
            ENDIF

            IF(I.GT.2) THEN
              A(L)=DW
              ISC_A(L)=M-1
              L=L+1
            ELSE
              B(M)=B(M)-DW*VW
            ENDIF

            A(L)=DP
            ISC_A(L)=M
            L=L+1

            IF(I.LT.NX-1) THEN
              A(L)=DE
              ISC_A(L)=M+1
              L=L+1
            ELSE
              B(M)=B(M)-DE*VE
            ENDIF

            IF(J.LT.NY-1) THEN
              A(L)=DN
              ISC_A(L)=M+NNX
              L=L+1
            ELSE
              B(M)=B(M)-DN*VN
            ENDIF

            IF(K.LT.NZ-1) THEN
              A(L)=DT
              ISC_A(L)=M+NNX*NNY
              L=L+1
            ELSE
              B(M)=B(M)-DT*VT
            ENDIF

          ENDDO
        ENDDO
      ENDDO
      M=M+1
      ISR_A(M)=L

C     NNX=NX-4
C     NNY=NY-4
C     NNZ=NZ-4
C     NZA=8*4 + 4*(NNX+NNY+NNZ)*5 + 2*(NNX*NNY+NNY*NNZ+NNZ*NNX)*6
C    '  + NNX*NNY*NNZ*7

      RETURN

 9999 CALL ERRORS('SET_FD3D_EQUATIONS',ERROR)
      RETURN 1
      END


      SUBROUTINE SET_FD2D_EQUATIONS(A,LDA,N,B,NX,NY,MTYPE,ISC_A,ISR_A,
     '  NZA,TRANSIENT,ERROR,*)

C#### Subroutine: SET_FD2D_EQUATIONS
C###  Description:
C###    Set a system of 2D Finite Difference equations to test the solvers.
C###  Written by Stuart Norris 12/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NX,NY,MTYPE,ISC_A(*),ISR_A(*),NZA
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
      LOGICAL TRANSIENT
!     Local Variables
      INTEGER I,J,L,M
      INTEGER NNX,NNY
      REAL*8 X(NX),Y(NY),DX(NX),DY(NY)
      REAL*8 AE,AW,AN,AS,DP,DE,DW,DN,DS
      REAL*8 VE,VW,VN,VS,PI,AREAX,AREAY,VOLUME,RHO,DELT

      DATA VE,VN,VW,VS,RHO / 2*1.0D0,2*-1.0D0,1.0D0 /


      NNX=NX-2
      NNY=NY-2

C     Set the mesh
      ! Regular mesh
      IF(MTYPE.EQ.1) THEN
        DO I=1,NX
          X(I)=DBLE(I-1)/DBLE(NX-1)
        ENDDO
        DO J=1,NY
          Y(J)=DBLE(J-1)/DBLE(NY-1)
        ENDDO

      ! Sin mesh
      ELSE IF(MTYPE.EQ.2) THEN
        PI=4.0D0*ATAN(1.0D0)
        DO I=1,NX
          X(I)=0.5D0*( 1.0D0 - COS(PI*DBLE(I-1)/DBLE(NX-1)) )
        ENDDO
        DO J=1,NY
          Y(J)=0.5D0*( 1.0D0 - COS(PI*DBLE(J-1)/DBLE(NY-1)) )
        ENDDO

      ! Error
      ELSE
        ERROR='>>Unknown mesh type'
        GOTO 9999
      ENDIF

      DO I=2,NX
        DX(I)=X(I)-X(I-1)
      ENDDO
      DO J=2,NY
        DY(J)=Y(J)-Y(J-1)
      ENDDO
      DELT=2.0D0*MIN(X(NX)/DBLE(NX-1),Y(NY)/DBLE(NY-1))**2

C     Set the system
      M=0
      L=1
      DO J=2,NY-1
        DO I=2,NX-1
          M=M+1
          ISR_A(M)=L

          AREAX=0.5D0*(DY(J+1)+DY(J))
          AREAY=0.5D0*(DX(I+1)+DX(I))
          VOLUME=AREAX*AREAY

          AE=AREAX/DX(I+1)
          AW=AREAX/DX(I)
          AN=AREAY/DY(J+1)
          AS=AREAY/DY(J)

          DE=-AE
          DW=-AW
          DN=-AN
          DS=-AS
          DP=(AE+AW+AN+AS)

          IF(TRANSIENT) DP=DP+VOLUME*RHO/DELT
          B(M)=0.0D0

          IF(J.GT.2) THEN
            A(L)=DS
            ISC_A(L)=M-NNX
            L=L+1
          ELSE
            B(M)=B(M)-DS*VS
          ENDIF

          IF(I.GT.2) THEN
            A(L)=DW
            ISC_A(L)=M-1
            L=L+1
          ELSE
            B(M)=B(M)-DW*VW
          ENDIF

          A(L)=DP
          ISC_A(L)=M
          L=L+1

          IF(I.LT.NX-1) THEN
            A(L)=DE
            ISC_A(L)=M+1
            L=L+1
          ELSE
            B(M)=B(M)-DE*VE
          ENDIF

          IF(J.LT.NY-1) THEN
            A(L)=DN
            ISC_A(L)=M+NNX
            L=L+1
          ELSE
            B(M)=B(M)-DN*VN
          ENDIF

        ENDDO
      ENDDO
      M=M+1
      ISR_A(M)=L

C     NNX=NX-4
C     NNY=NY-4
C     NZA=4*3 + 2*(NNX+NNY)*4 + NNX*NNY*5

      RETURN

 9999 CALL ERRORS('SET_FD2D_EQUATIONS',ERROR)
      RETURN 1
      END


      SUBROUTINE SET_FV2D_EQUATIONS(A,LDA,N,B,NX,NY,MTYPE,ISC_A,ISR_A,
     '  NZA,TRANSIENT,ERROR,*)

C#### Subroutine: SET_FV2D_EQUATIONS
C###  Description:
C###    Set a system of 2D Finite Volume equations to test the solvers.
C###  Written by Stuart Norris 12/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NX,NY,MTYPE,ISC_A(*),ISR_A(*),NZA
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
      LOGICAL TRANSIENT
!     Local Variables
      INTEGER I,J,K,L,II,JJ
      INTEGER NNX,NNY,INDEX(-1:1,-1:1)
      REAL*8 X(NX),Y(NY),DX(NX),DY(NY),BC(NX,NY),D(-1:1,-1:1)
      REAL*8 AE,AW,AN,AS,VE,VW,VN,VS,PI,AREAX,AREAY,VOLUME,RHO,DELT
      LOGICAL LB(NX,NY)

      DATA VE,VN,VW,VS,RHO / 2*1.0D0,2*-1.0D0,1.0D0 /


      NNX=NX-2
      NNY=NY-2

C     Set the mesh
      ! Regular mesh
      IF(MTYPE.EQ.1) THEN
        DO I=1,NX
          X(I)=DBLE(I-1)/DBLE(NX-1)
        ENDDO
        DO J=1,NY
          Y(J)=DBLE(J-1)/DBLE(NY-1)
        ENDDO

      ! Sin mesh
      ELSE IF(MTYPE.EQ.2) THEN
        PI=4.0D0*ATAN(1.0D0)
        DO I=1,NX
          X(I)=0.5D0*( 1.0D0 - COS(PI*DBLE(I-1)/DBLE(NX-1)) )
        ENDDO
        DO J=1,NY
          Y(J)=0.5D0*( 1.0D0 - COS(PI*DBLE(J-1)/DBLE(NY-1)) )
        ENDDO

      ! Error
      ELSE
        ERROR='>>Unknown mesh type'
        GOTO 9999
      ENDIF

      DO I=2,NX
        DX(I)=X(I)-X(I-1)
      ENDDO
      DO J=2,NY
        DY(J)=Y(J)-Y(J-1)
      ENDDO
      DELT=2.0D0*MIN(X(NX)/DBLE(NX),Y(NY)/DBLE(NY))**2

C     Set the boundary condition flags
      DO J=1,NY
        DO I=1,NX
          BC(I,J)=0.0D0
          LB(I,J)=.FALSE.
        ENDDO
      ENDDO

      DO I=1,NX
        BC(I, 1)=VS
        BC(I,NY)=VN
        LB(I, 1)=.TRUE.
        LB(I,NY)=.TRUE.
      ENDDO
      DO J=1,NY
        BC( 1,J)=VW
        BC(NX,J)=VE
        LB( 1,J)=.TRUE.
        LB(NX,J)=.TRUE.
      ENDDO
      BC( 1, 1)=0.5D0*(VS+VW)
      BC(NX, 1)=0.5D0*(VS+VE)
      BC( 1,NY)=0.5D0*(VN+VW)
      BC(NX,NY)=0.5D0*(VN+VE)

      DO J=-1,1
        DO I=-1,1
          INDEX(I,J)=J*NNX+I
        ENDDO
      ENDDO

C     Set the system
      K=0
      L=1
      DO J=2,NY-1
        DO I=2,NX-1
          K=K+1
          ISR_A(K)=L

          AREAX=0.5D0*(DY(J+1)+DY(J))
          AREAY=0.5D0*(DX(I+1)+DX(I))
          VOLUME=AREAX*AREAY

          AE=AREAX/DX(I+1)
          AW=AREAX/DX(I)
          AN=AREAY/DY(J+1)
          AS=AREAY/DY(J)

          D( 0, 0)=(AE+AW+AN+AS)*6.0D0/8.0D0 ! DP

          D( 1, 0)=((AN+AS)-6.0D0*AE)/8.0D0 ! DE
          D(-1, 0)=((AN+AS)-6.0D0*AW)/8.0D0 ! DW
          D( 0, 1)=((AE+AW)-6.0D0*AN)/8.0D0 ! DN
          D( 0,-1)=((AE+AW)-6.0D0*AS)/8.0D0 ! DS

          D( 1, 1)=-((AE+AN)/8.0D0) ! DNE
          D(-1, 1)=-((AW+AN)/8.0D0) ! DNW
          D( 1,-1)=-((AE+AS)/8.0D0) ! DSE
          D(-1,-1)=-((AW+AS)/8.0D0) ! DSW

          IF(TRANSIENT) D(0,0)=D(0,0)+VOLUME*RHO/DELT
          B(K)=0.0D0

          DO JJ=-1,1
            DO II=-1,1
              IF(LB(I+II,J+JJ)) THEN
                B(K)=B(K)-D(II,JJ)*BC(I+II,J+JJ)
              ELSE
                A(L)=D(II,JJ)
                ISC_A(L)=K+INDEX(II,JJ)
                L=L+1
              ENDIF
            ENDDO
          ENDDO

        ENDDO
      ENDDO
      K=K+1
      ISR_A(K)=L

C     NNX=NX-4
C     NNY=NY-4
C     NZA=4*4 + 2*(NNX+NNY)*6 + NNX*NNY*9

      RETURN

 9999 CALL ERRORS('SET_FV2D_EQUATIONS',ERROR)
      RETURN 1
      END


      SUBROUTINE SET_FV3D_EQUATIONS(A,LDA,N,B,NX,NY,NZ,MTYPE,ISC_A,
     '  ISR_A,NZA,TRANSIENT,ERROR,*)

C#### Subroutine: SET_FV3D_EQUATIONS
C###  Description:
C###    Set a system of 3D Finite Volume equations to test the solvers.
C###  Written by Stuart Norris 03/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NX,NY,NZ,MTYPE,ISC_A(*),ISR_A(*),NZA
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
      LOGICAL TRANSIENT
!     Local Variables
      INTEGER I,J,K,L,M,II,JJ,KK
      INTEGER NNX,NNY,NNZ,INDEX(-1:1,-1:1,-1:1)
      REAL*8 X(NX),Y(NY),Z(NZ),DX(NX),DY(NY),DZ(NZ),BC(NX,NY,NZ)
      REAL*8 AE,AW,AN,AS,AT,AB,VE,VW,VN,VS,VT,VB,D(-1:1,-1:1,-1:1)
      REAL*8 RHO,DELT,PI,AREAX,AREAY,AREAZ,VOLUME
      LOGICAL LB(NX,NY,NZ)

      DATA VE,VN,VT,VW,VS,VB,RHO / 3*1.0D0,3*-1.0D0,1.0D0 /


      NNX=NX-2
      NNY=NY-2
      NNZ=NZ-2

C     Set the mesh
      ! Regular mesh
      IF(MTYPE.EQ.1) THEN
        DO I=1,NX
          X(I)=DBLE(I-1)/DBLE(NX-1)
        ENDDO
        DO J=1,NY
          Y(J)=DBLE(J-1)/DBLE(NY-1)
        ENDDO
        DO K=1,NZ
          Z(K)=DBLE(K-1)/DBLE(NZ-1)
        ENDDO

      ! Sin mesh
      ELSE IF(MTYPE.EQ.2) THEN
        PI=4.0D0*ATAN(1.0D0)
        DO I=1,NX
          X(I)=0.5D0*( 1.0D0 - COS(PI*DBLE(I-1)/DBLE(NX-1)) )
        ENDDO
        DO J=1,NY
          Y(J)=0.5D0*( 1.0D0 - COS(PI*DBLE(J-1)/DBLE(NY-1)) )
        ENDDO
        DO K=1,NZ
          Z(K)=0.5D0*( 1.0D0 - COS(PI*DBLE(K-1)/DBLE(NZ-1)) )
        ENDDO

      ! Error
      ELSE
        ERROR='>>Unknown mesh type'
        GOTO 9999
      ENDIF

      DO I=2,NX
        DX(I)=X(I)-X(I-1)
      ENDDO
      DO J=2,NY
        DY(J)=Y(J)-Y(J-1)
      ENDDO
      DO K=2,NZ
        DZ(K)=Z(K)-Z(K-1)
      ENDDO
      DELT=2.0D0*MIN(X(NX)/DBLE(NX),Y(NY)/DBLE(NY),Z(NZ)/DBLE(NZ))**2

C     Set the boundary condition flags
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
            BC(I,J,K)=0.0D0
            LB(I,J,K)=.FALSE.
          ENDDO
        ENDDO
      ENDDO

      DO J=1,NY
        DO I=1,NX
          BC(I,J, 1)=VB
          BC(I,J,NZ)=VT
          LB(I,J, 1)=.TRUE.
          LB(I,J,NZ)=.TRUE.
        ENDDO
      ENDDO
      DO K=1,NZ
        DO I=1,NX
          BC(I, 1,K)=VS
          BC(I,NY,K)=VN
          LB(I, 1,K)=.TRUE.
          LB(I,NY,K)=.TRUE.
        ENDDO
      ENDDO
      DO K=1,NZ
        DO J=1,NY
          BC( 1,J,K)=VW
          BC(NX,J,K)=VE
          LB( 1,J,K)=.TRUE.
          LB(NX,J,K)=.TRUE.
        ENDDO
      ENDDO

      DO I=1,NX
        BC(I, 1, 1)=0.5D0*(VS+VB)
        BC(I,NY, 1)=0.5D0*(VN+VB)
        BC(I, 1,NZ)=0.5D0*(VS+VT)
        BC(I,NY,NZ)=0.5D0*(VN+VT)
      ENDDO
      DO J=1,NY
        BC( 1,J, 1)=0.5D0*(VW+VB)
        BC(NX,J, 1)=0.5D0*(VE+VB)
        BC( 1,J,NZ)=0.5D0*(VW+VT)
        BC(NX,J,NZ)=0.5D0*(VE+VT)
      ENDDO
      DO K=1,NZ
        BC( 1, 1,K)=0.5D0*(VW+VS)
        BC(NX, 1,K)=0.5D0*(VE+VS)
        BC( 1,NY,K)=0.5D0*(VW+VN)
        BC(NX,NY,K)=0.5D0*(VE+VN)
      ENDDO

      BC( 1, 1, 1)=(VW+VS+VB)/3.0D0
      BC(NX, 1, 1)=(VE+VS+VB)/3.0D0
      BC( 1,NY, 1)=(VW+VN+VB)/3.0D0
      BC(NX,NY, 1)=(VE+VN+VB)/3.0D0
      BC( 1, 1,NZ)=(VW+VS+VT)/3.0D0
      BC(NX, 1,NZ)=(VE+VS+VT)/3.0D0
      BC( 1,NY,NZ)=(VW+VN+VT)/3.0D0
      BC(NX,NY,NZ)=(VE+VN+VT)/3.0D0

      DO K=-1,1
        DO J=-1,1
          DO I=-1,1
            INDEX(I,J,K)=K*NNX*NNY+J*NNX+I
          ENDDO
        ENDDO
      ENDDO

C     Set the system
      M=0
      L=1
      DO K=2,NZ-1
        DO J=2,NY-1
          DO I=2,NX-1
            M=M+1
            ISR_A(M)=L

            AREAX=0.25D0*(DY(J+1)+DY(J))*(DZ(K+1)+DZ(K))
            AREAY=0.25D0*(DX(I+1)+DX(I))*(DZ(K+1)+DZ(K))
            AREAZ=0.25D0*(DX(I+1)+DX(I))*(DY(J+1)+DY(J))
            VOLUME=SQRT(AREAX*AREAY*AREAZ)

            AE=AREAX/DX(I+1)
            AW=AREAX/DX(I)
            AN=AREAY/DY(J+1)
            AS=AREAY/DY(J)
            AT=AREAZ/DZ(K+1)
            AB=AREAZ/DZ(K)

            D( 0, 0, 0)=(AE+AW+AN+AS+AT+AB)*36.0D0/64.0D0 ! DP

            D( 1, 0, 0)=(6.0D0*(AN+AS+AT+AB) - 36.0D0*AE)/64.0D0 ! DE
            D(-1, 0, 0)=(6.0D0*(AN+AS+AT+AB) - 36.0D0*AW)/64.0D0 ! DW
            D( 0, 1, 0)=(6.0D0*(AE+AW+AT+AB) - 36.0D0*AN)/64.0D0 ! DN
            D( 0,-1, 0)=(6.0D0*(AE+AW+AT+AB) - 36.0D0*AS)/64.0D0 ! DS
            D( 0, 0, 1)=(6.0D0*(AE+AW+AN+AS) - 36.0D0*AT)/64.0D0 ! DT
            D( 0, 0,-1)=(6.0D0*(AE+AW+AN+AS) - 36.0D0*AB)/64.0D0 ! DB

            D( 1, 1, 0)=((AT+AB)-6.0D0*(AE+AN))/64.0D0 ! DNE
            D(-1, 1, 0)=((AT+AB)-6.0D0*(AW+AN))/64.0D0 ! DNW
            D( 1,-1, 0)=((AT+AB)-6.0D0*(AE+AS))/64.0D0 ! DSE
            D(-1,-1, 0)=((AT+AB)-6.0D0*(AW+AS))/64.0D0 ! DSW

            D( 1, 0, 1)=((AN+AS)-6.0D0*(AE+AT))/64.0D0 ! DET
            D( 1, 0,-1)=((AN+AS)-6.0D0*(AE+AB))/64.0D0 ! DEB
            D(-1, 0, 1)=((AN+AS)-6.0D0*(AW+AT))/64.0D0 ! DWT
            D(-1, 0,-1)=((AN+AS)-6.0D0*(AW+AB))/64.0D0 ! DWB

            D( 0, 1, 1)=((AE+AW)-6.0D0*(AT+AN))/64.0D0 ! DNT
            D( 0, 1,-1)=((AE+AW)-6.0D0*(AB+AN))/64.0D0 ! DNB
            D( 0,-1, 1)=((AE+AW)-6.0D0*(AT+AS))/64.0D0 ! DST
            D( 0,-1,-1)=((AE+AW)-6.0D0*(AB+AS))/64.0D0 ! DSB

            D( 1, 1, 1)=-((AN+AE+AT)/64.0D0) ! DNET
            D( 1, 1,-1)=-((AN+AE+AB)/64.0D0) ! DNEB
            D(-1, 1, 1)=-((AN+AW+AT)/64.0D0) ! DNWT
            D(-1, 1,-1)=-((AN+AW+AB)/64.0D0) ! DNWB
            D( 1,-1, 1)=-((AS+AE+AT)/64.0D0) ! DSET
            D( 1,-1,-1)=-((AS+AE+AB)/64.0D0) ! DSEB
            D(-1,-1, 1)=-((AS+AW+AT)/64.0D0) ! DSWT
            D(-1,-1,-1)=-((AS+AW+AB)/64.0D0) ! DSWB

            IF(TRANSIENT) D(0,0,0)=D(0,0,0)+VOLUME*RHO/DELT
            B(M)=0.0D0

            DO KK=-1,1
              DO JJ=-1,1
                DO II=-1,1
                  IF(LB(I+II,J+JJ,K+KK)) THEN
                    B(M)=B(M)-D(II,JJ,KK)*BC(I+II,J+JJ,K+KK)
                  ELSE
                    A(L)=D(II,JJ,KK)
                    ISC_A(L)=M+INDEX(II,JJ,KK)
                    L=L+1
                  ENDIF
                ENDDO
              ENDDO
            ENDDO

          ENDDO
        ENDDO
      ENDDO
      M=M+1
      ISR_A(M)=L

C     NNX=NX-4
C     NNY=NY-4
C     NNZ=NZ-4
C     NZA = 8*8 + 4*(NNX+NNY+NNZ)*12 + 2*(NNX*NNY+NNY*NNZ+NNZ*NNX)*18
C    '  + NNX*NNY*NNZ*27

      RETURN

 9999 CALL ERRORS('SET_FV3D_EQUATIONS',ERROR)
      RETURN 1
      END


      SUBROUTINE PRINT_ONE_LINE(A,ISR_A,ISC_A,N)

C#### Subroutine: PRINT_ONE_LINE
C###  Description:
C###    Print one line of a system of equations -- for checking a system is
C###    diagonally dominant.
C###  Written by Stuart Norris 03/07/02

      IMPLICIT NONE
!     Subroutine parameters
      INTEGER N,ISC_A(*),ISR_A(*)
      REAL*8 A(*)
!     Local variables
      INTEGER J,M
      REAL*8 SUM

      M=N/2
      SUM=0.0D0
      DO J=ISR_A(M),ISR_A(M+1)-1
        PRINT '(2(I5,1X),F8.5)',J,ISC_A(J),A(J)
        SUM=SUM+A(J)
      ENDDO
      PRINT '(A,F8.5)', 'Sum = ',SUM

      RETURN
      END


      SUBROUTINE SET_FV2D_EQUATIONS_BOUNDS(A,LDA,N,B,NX,NY,MTYPE,ISC_A,
     '  ISR_A,NZA,TRANSIENT,ERROR,*)

C#### Subroutine: SET_FV2D_EQUATIONS_BOUNDS
C###  Description:
C###    Set a system of 2D Finite Volume equations to test the solvers.
C###  Written by Stuart Norris 12/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NX,NY,MTYPE,ISC_A(*),ISR_A(*),NZA
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
      LOGICAL TRANSIENT
!     Local Variables
      INTEGER I,J,K,L,II,JJ
      INTEGER NNX,NNY,INDEX(-1:1,-1:1)
      REAL*8 X(NX),Y(NY),DX(NX+1),DY(NY+1),BC(NX,NY),D(-1:1,-1:1)
      REAL*8 AE,AW,AN,AS,VE,VW,VN,VS,PI,AREAX,AREAY,VOLUME,RHO,DELT
      LOGICAL LB(NX,NY)

      DATA VE,VN,VW,VS,RHO / 2*1.0D0,2*-1.0D0,1.0D0 /


      NNX=NX
      NNY=NY

C     Set the mesh
      ! Regular mesh
      IF(MTYPE.EQ.1) THEN
        DO I=1,NX
          X(I)=DBLE(I-1)/DBLE(NX-1)
        ENDDO
        DO J=1,NY
          Y(J)=DBLE(J-1)/DBLE(NY-1)
        ENDDO

      ! Sin mesh
      ELSE IF(MTYPE.EQ.2) THEN
        PI=4.0D0*ATAN(1.0D0)
        DO I=1,NX
          X(I)=0.5D0*( 1.0D0 - COS(PI*DBLE(I-1)/DBLE(NX-1)) )
        ENDDO
        DO J=1,NY
          Y(J)=0.5D0*( 1.0D0 - COS(PI*DBLE(J-1)/DBLE(NY-1)) )
        ENDDO

      ! Error
      ELSE
        ERROR='>>Unknown mesh type'
        GOTO 9999
      ENDIF

      DO I=2,NX
        DX(I)=X(I)-X(I-1)
      ENDDO
      DX(1)=DX(2)
      DX(NX+1)=DX(NX)
      DO J=2,NY
        DY(J)=Y(J)-Y(J-1)
      ENDDO
      DY(1)=DY(2)
      DY(NY+1)=DX(NY)
      DELT=2.0D0*MIN(X(NX)/DBLE(NX),Y(NY)/DBLE(NY))**2

C     Set the boundary condition flags
      DO J=1,NY
        DO I=1,NX
          BC(I,J)=0.0D0
          LB(I,J)=.FALSE.
        ENDDO
      ENDDO

      DO I=1,NX
        BC(I, 1)=VS
        BC(I,NY)=VN
        LB(I, 1)=.TRUE.
        LB(I,NY)=.TRUE.
      ENDDO
      DO J=1,NY
        BC( 1,J)=VW
        BC(NX,J)=VE
        LB( 1,J)=.TRUE.
        LB(NX,J)=.TRUE.
      ENDDO
      BC( 1, 1)=0.5D0*(VS+VW)
      BC(NX, 1)=0.5D0*(VS+VE)
      BC( 1,NY)=0.5D0*(VN+VW)
      BC(NX,NY)=0.5D0*(VN+VE)

      DO J=-1,1
        DO I=-1,1
          INDEX(I,J)=J*NNX+I
        ENDDO
      ENDDO

C     Set the system
      K=0
      L=1
      DO J=1,NY
        DO I=1,NX
          K=K+1
          ISR_A(K)=L

          AREAX=0.5D0*(DY(J+1)+DY(J))
          AREAY=0.5D0*(DX(I+1)+DX(I))
          VOLUME=AREAX*AREAY

          AE=AREAX/DX(I+1)
          AW=AREAX/DX(I)
          AN=AREAY/DY(J+1)
          AS=AREAY/DY(J)

          D( 0, 0)=(AE+AW+AN+AS)*6.0D0/8.0D0 ! DP

          D( 1, 0)=((AN+AS)-6.0D0*AE)/8.0D0 ! DE
          D(-1, 0)=((AN+AS)-6.0D0*AW)/8.0D0 ! DW
          D( 0, 1)=((AE+AW)-6.0D0*AN)/8.0D0 ! DN
          D( 0,-1)=((AE+AW)-6.0D0*AS)/8.0D0 ! DS

          D( 1, 1)=-((AE+AN)/8.0D0) ! DNE
          D(-1, 1)=-((AW+AN)/8.0D0) ! DNW
          D( 1,-1)=-((AE+AS)/8.0D0) ! DSE
          D(-1,-1)=-((AW+AS)/8.0D0) ! DSW

          IF(TRANSIENT) D(0,0)=D(0,0)+VOLUME*RHO/DELT
          B(K)=0.0D0

          IF(LB(I,J)) THEN
            DO JJ=-1,1
              DO II=-1,1
                IF(II.EQ.0.AND.JJ.EQ.0) THEN            
                  A(L)=1.0D0
                  B(K)=BC(I,J)
                  ISC_A(L)=K
                  L=L+1
                ELSE IF( ((I+II).GE.1).AND.((I+II).LE.NX)
     '              .AND.((J+JJ).GE.1).AND.((J+JJ).LE.NY)) THEN
                  A(L)=0.0D0
                  ISC_A(L)=K+INDEX(II,JJ)
                  L=L+1
                ENDIF
              ENDDO
            ENDDO
          ELSE
            DO JJ=-1,1
              DO II=-1,1
                A(L)=D(II,JJ)
                ISC_A(L)=K+INDEX(II,JJ)
                L=L+1
              ENDDO
            ENDDO
          ENDIF

        ENDDO
      ENDDO
      K=K+1
      ISR_A(K)=L
      NZA=L-1

C     NNX=NX-4
C     NNY=NY-4
C     NZA=4*4 + 2*(NNX+NNY)*6 + NNX*NNY*9

      RETURN

 9999 CALL ERRORS('SET_FV2D_EQUATIONS_BOUNDS',ERROR)
      RETURN 1
      END
