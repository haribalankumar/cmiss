      SUBROUTINE XECURV(IBT,IDO,INP,NBJ,EG,PCU,XE,XI,ERROR,*)

C#### Subroutine: XECURV
C###  Description:
C###    XECURV calculates curvatures for element with geometric nodal
C###    parameters XE and local coordinate values XI.
C###    Calculates principal curvatures in XI_1 and XI_2 plane and
C###    curvatures in XI_1 and XI_2 directions  at XI.
C###    Z refers to the global rectangular Cartesian coordinate system,
C###    X refers to the global curvilinear coordinate system, and
C###    XI refers to the local curvilinear coordinate system.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM)
      REAL*8 EG(3,3),PCU(3),XE(NSM,NJM),XI(NIM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,NITB,NU1(3),NU2(3,3),lj,mi,mj,ni,nj
      REAL*8 dXrc_dXref,DXXI(3,3),D2XXI(3,3,3),D2ZX,D2ZXI(3,3,3),
     '  DZX,DZXI(3,3),LEN,NORMAL(3),PXI,SUM,X(3),Z(3),
     '  A,B,C,E,F,G,L,M,N
      DATA NU1/2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('XECURV',*9999)

C***  first calculate the first and second derivs of rect cart coords
C     in the xi_1 and xi_2 directions

C**** I call this for curvature calculation from XECURV AAY 2 June 1994
C**** Calculates X and Z coordinates for element with geometric nodal
C**** parameters XE and local coordinate values XI.
C**** Calculates DX(I)/DXI(J), D2X(I)/DXI(J)DXI(K), DZ(I)/DXI(J), and
C**** D2Z(I)/DXI(J)DXI(K) at XI.
C**** Z refers to the global rectangular Cartesian coordinate system,
C**** X refers to the global curvilinear coordinate system, and
C**** XI refers to the local curvilinear coordinate system.

      !new AAY 4 Oct 94 changed calls to PXI
      NITB=NIT(NBJ(1)) !new AAY take NITB from nj=1 2 June 1994
      DO nj=1,NJT
        nb=NBJ(nj)
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
      ENDDO
      CALL XZ(JTYP3,X,Z)
      DO ni=1,NITB
        DO nj=1,NJT
          nb=NBJ(nj)
          DXXI(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      NU1(ni),XI,XE(1,nj))
        ENDDO
      ENDDO
      DO ni=1,NITB
        DO mi=1,NITB
          DO nj=1,NJT
            nb=NBJ(nj)
            D2XXI(nj,ni,mi)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,NU2(ni,mi),XI,XE(1,nj))
          ENDDO
        ENDDO
      ENDDO
      DO ni=1,NITB
        DO nj=1,NJT
          SUM=0.d0
          DO mj=1,NJT
            dXrc_dXref=DZX(JTYP3,nj,mj,X)
            SUM=SUM+dXrc_dXref*DXXI(mj,ni)
          ENDDO
          DZXI(nj,ni)=SUM
        ENDDO !nj
      ENDDO !ni
      DO ni=1,NITB
        DO mi=1,NITB
          DO nj=1,NJT
            SUM=0.d0
            DO mj=1,NJT
              dXrc_dXref=DZX(JTYP3,nj,mj,X)
              SUM=SUM+dXrc_dXref*D2XXI(mj,ni,mi)
              DO lj=1,NJT
                SUM=SUM+D2ZX(JTYP3,nj,mj,lj,X)*DXXI(mj,ni)*DXXI(lj,mi)
              ENDDO !lj
            ENDDO !mj
            D2ZXI(nj,ni,mi)=SUM
          ENDDO !nj
        ENDDO !mi
      ENDDO !ni

      E=0.d0
      DO nj=1,NJT
        E = E+DZXI(nj,1)*DZXI(nj,1)
      ENDDO
      F=0.d0
      DO nj=1,NJT
        F = F+DZXI(nj,1)*DZXI(nj,2)
      ENDDO
      G=0
      DO nj=1,NJT
        G = G+DZXI(nj,2)*DZXI(nj,2)
      ENDDO
      NORMAL(1)=DZXI(2,1)*DZXI(3,2)-DZXI(3,1)*DZXI(2,2)
      NORMAL(2)=DZXI(3,1)*DZXI(1,2)-DZXI(1,1)*DZXI(3,2)
      NORMAL(3)=DZXI(1,1)*DZXI(2,2)-DZXI(2,1)*DZXI(1,2)
      LEN = SQRT(NORMAL(1)**2+NORMAL(2)**2+NORMAL(3)**2)
      NORMAL(1)=NORMAL(1)/LEN
      NORMAL(2)=NORMAL(2)/LEN
      NORMAL(3)=NORMAL(3)/LEN
      L=0
      DO nj=1,NJT
        L = L + D2ZXI(nj,1,1)*NORMAL(nj)
      ENDDO
      M=0
      DO nj=1,NJT
        M = M + D2ZXI(nj,1,2)*NORMAL(nj)
      ENDDO
      N=0
      DO nj=1,NJT
        N = N + D2ZXI(nj,2,2)*NORMAL(nj)
      ENDDO

      IF(DABS(E).LT.1.d-5.OR.DABS(G).LT.1.d-5)THEN
        EG(1,1)=0.d0
        EG(2,2)=0.d0
        PCU(1) =0.d0
        PCU(2) =0.d0
        RETURN
      ENDIF
      EG(1,1)=L/E
      EG(2,2)=N/G

      A = E*G-F*F
      B = 2.d0*F*M-E*N-G*L
      C = L*N-M*M
      PCU(1)= (-B+DSQRT(B*B-4.d0*A*C))/(2.d0*A)
      PCU(2)= (-B-DSQRT(B*B-4.d0*A*C))/(2.d0*A)

      CALL EXITS('XECURV')
      RETURN
 9999 CALL ERRORS('XECURV',ERROR)
      CALL EXITS('XECURV')
      RETURN 1
      END


