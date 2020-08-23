      SUBROUTINE GAUSS10(IBT,IDO,INP,nb,NGAP,ALIM,BLIM,PG,WG,XIG,
     '  ERROR,*)
C#### Subroutine: GAUSS10
C###  Description:
C###    GAUSS10 defines the Gaussian quadrature coordinates XIG and
C###    weights WG and evaluates the basis function Gauss point array
C###    PG for Lagrange/Hermite  tensor product type basis function when
C###    using (high order) boundary element quadrature.

C**** Currently only 64 quadrature points can be used in each direction.
C**** ALIM and BLIM are defined in IPBASE and refer to the lower and
C**** upper limits of the integral for which quadrature schemes are
C**** required.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,NGAP(NIM)
      REAL*8 ALIM(2),BLIM(2),PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,ng,ng1,ng2,NGMAX,ni,nk,nn,ns,nu
      PARAMETER (NGMAX = 64)
      REAL*8 D1(NGMAX),D2(NGMAX),PSI1,SUM,SUML,SUMC,SUMPG1,SUMPG2,
     ' SUMPG3,W1(NGMAX),W2(NGMAX),XIGG(NGMAX,NGMAX,2),XI(3),TEST(16)
C ***  The array XIG(i,ng)contains the positions of the gauss points
C ***  and the respective derivatives.
C ***              i = 1 XI1 coordinate
C ***                  2 XI2 coordinate
C ***                  3 derivative in XI1 direction
C ***                  4 derivative in XI2 direction

      CALL ENTERS('GAUSS10',*9999)
      ng1=NGAP(1)
C old MLB 18/3/97 CALL GAUSSPWB(ALIM(1),BLIM(1),ng1,-1,W1,D1,10,ERROR,*9999)
      CALL GAUSSLEG(ALIM(1),BLIM(1),D1,W1,ng1,ERROR,*9999)

      ng2=1
      W2(1)=1
      IF(NIT(nb).GT.1)THEN
        ng2=NGAP(2)
C old MLB 18/3/97 CALL GAUSSPWB(ALIM(2),BLIM(2),ng2,-1,W2,D2,10,ERROR,*9999)
         CALL GAUSSLEG(ALIM(2),BLIM(2),D2,W2,ng2,ERROR,*9999)
      ENDIF

      IF(DOP) THEN
C GMH 18/12/96 Initialise variables
        SUM=0.0d0
        SUML=0.0d0
        SUMC=0.0d0
        SUMPG1=0.0d0
        SUMPG2=0.0d0
        SUMPG3=0.0d0
        DO I=1,16
          TEST(I)=0.0d0
        ENDDO !I
      ENDIF

C LKC 4-SEP-97 add assert line on NGM
      CALL ASSERT(NGM.GE.ng1+(ng2-1)*ng1,'>>Increase NGM',ERROR,*9999)
      DO J=1,ng2
        DO I=1,ng1
          XIGG(I,J,1)=D1(I)
          IF(NIT(nb).GT.1)THEN
            XIGG(I,J,2)=D2(J)
          ENDIF
          ng=I+(J-1)*ng1
          WG(ng)=W1(I)*W2(J)
          DO ni=1,NIT(nb)
            XI(ni)=XIGG(I,J,ni)
            XIG(ni,ng)=XI(ni)
          ENDDO
          ns=0

C LKC 24-APR-1998 added assert

          DO nn=1,NNT(nb)
          CALL ASSERT(NSM.GE.NKT(NNT(nb),nb),
     '        '>>Increase NSM',ERROR,*9999)
            DO nk=1,NKT(nn,nb)
              ns=ns+1
              DO nu=1,NUT(nb)
                PG(ns,nu,ng)=
     '            PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)
              ENDDO
            ENDDO
          ENDDO
          IF(DOP)THEN
C (for 1d)   SUML=SUML+(DLOG(XIG(1,ng)))*WG(ng)
C (for 1d)   SUMC=SUMC+(3.0D0+(XIG(1,ng)/2.0D0)-XIG(2,ng)-
C                 XIG(1,ng)*XIG(2,ng))*WG(ng)

            SUML=SUML+1.0D0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '           (1.0D0-XIG(1,ng))*(1.0D0-XIG(2,ng))*WG(ng)
            SUM=SUM+WG(ng)
            SUMPG1=SUMPG1+PG(1,1,ng)*WG(ng)
            SUMPG2=SUMPG2+PG(2,1,ng)*WG(ng)
            SUMPG3=SUMPG3+PG(3,1,ng)*WG(ng)
            TEST(1)=TEST(1)+1.0D0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        (1.0D0-XIG(1,ng))*(1.0D0-XIG(2,ng))*WG(ng)
            TEST(2)=TEST(2)+1.0D0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        (1.0D0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(3)=TEST(3)+1.0D0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        XIG(1,ng)*(1.0D0-XIG(2,ng))*WG(ng)
            TEST(4)=TEST(4)+1.0D0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
            TEST(5)=TEST(5)+1.0D0/(DSQRT(XIG(1,ng)**2+
     '        (1.0D0-XIG(2,ng))**2))*(1.0D0-XIG(1,ng))*
     '        (1.0D0-XIG(2,ng))*WG(ng)
            TEST(6)=TEST(6)+1.0D0/(DSQRT(XIG(1,ng)**2+
     '        (1.0D0-XIG(2,ng))**2))*
     '        (1.0D0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(7)=TEST(7)+1.0D0/(DSQRT(XIG(1,ng)**2+
     '        (1.0D0-XIG(2,ng))**2))*
     '        XIG(1,ng)*(1.0D0-XIG(2,ng))*WG(ng)
            TEST(8)=TEST(8)+1.0D0/(DSQRT(XIG(1,ng)**2+
     '        (1.0D0-XIG(2,ng))**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
            TEST(9)=TEST(9)+1.0D0/(DSQRT(XIG(2,ng)**2+
     '        (1.0D0-XIG(1,ng))**2))*
     '        (1.0D0-XIG(1,ng))*(1.0D0-XIG(2,ng))*WG(ng)
            TEST(10)=TEST(10)+1.0D0/(DSQRT(XIG(2,ng)**2+
     '        (1.0D0-XIG(1,ng))**2))*
     '        (1.0D0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(11)=TEST(11)+1.0D0/(DSQRT(XIG(2,ng)**2+
     '        (1.0D0-XIG(1,ng))**2))*
     '        XIG(1,ng)*(1.0D0-XIG(2,ng))*WG(ng)
            TEST(12)=TEST(12)+1.0D0/(DSQRT(XIG(2,ng)**2+
     '        (1.0D0-XIG(1,ng))**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
            TEST(13)=TEST(13)+1.0D0/(DSQRT((1.0D0-XIG(1,ng))**2+
     '        (1.0D0-XIG(2,ng))**2))*
     '        (1.0D0-XIG(1,ng))*(1.0D0-XIG(2,ng))*WG(ng)
            TEST(14)=TEST(14)+1.0D0/(DSQRT((1.0D0-XIG(1,ng))**2+
     '        (1.0D0-XIG(2,ng))**2))*
     '        (1.0D0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(15)=TEST(15)+1.0D0/(DSQRT((1.0d0-XIG(1,ng))**2+
     '        (1.0D0-XIG(2,ng))**2))*
     '        XIG(1,ng)*(1.0D0-XIG(2,ng))*WG(ng)
            TEST(16)=TEST(16)+1.0D0/(DSQRT((1.0D0-XIG(1,ng))**2+
     '        (1.0D0-XIG(2,ng))**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
           ENDIF
        ENDDO
      ENDDO

C ***  These two lines update the actual order of the quadrature scheme
C ***  that is being used.
      NGAP(1)=ng1
      IF(NIT(nb).GT.1) NGAP(2)=ng2

      IF(DOP)THEN
        WRITE(OP_STRING,*)' SUM OF WEIGHTS=',SUM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,*)' LOGE 0-1 =',SUML
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' 1/(4*pi*R)*(1-X)*(1-Y) =',
     '    SUML/(4.0D0*4.0D0*DATAN(1.0D0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' SUM OF FUNCTION = ',SUMC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' SUM BASIS FUNCTION 1=',SUMPG1
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' SUM BASIS FUNCTION 2=',SUMPG2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' SUM BASIS FUNCTION 3=',SUMPG3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO I=1,16
          WRITE(OP_STRING,*)' TEST(',I,')=',TEST(I)/
     '      (4.0D0*4.0D0*DATAN(1.0D0))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('GAUSS10')
      RETURN
9999  CALL ERRORS('GAUSS10',ERROR)
      CALL EXITS('GAUSS10')
      RETURN 1
      END


