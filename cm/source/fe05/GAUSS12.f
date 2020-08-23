      SUBROUTINE GAUSS12(IBT,IDO,INP,nb,NGAP,ALIM,BLIM,PG,WG,XIG,
     '  HERMSECTOR,ERROR,*)

C#### Subroutine: GAUSS12
C###  Description:
C###    GAUSS12 defines the Gaussian quadrature coords XIG and weights
C###    WG and evaluates the basis function Gauss point array PG for
C###    simplex/sector boundary elements.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,NGAP(NIM)
      REAL*8 ALIM(2),BLIM(2),PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
      LOGICAL HERMSECTOR
!     Local Variables
      INTEGER NGMAX
      PARAMETER (NGMAX=64)
      INTEGER i,j,ng,ng1,ng2,ni,nk,nn,ns,nu
      REAL*8 D1(NGMAX),D2(NGMAX),PSI2_HERMITE,PSI5,SUM,SUML,SUMC,
     '  SUMPG1,SUMPG2,SUMPG3,W1(NGMAX),W2(NGMAX),XIGG(NGMAX,NGMAX,2),
     '  XI(2),TEST(16)

c cpb 1/11/96 Generalling for sector elements. This routine used to be
C GAUSS6_HERMITE

      CALL ENTERS('GAUSS12',*9999)

      ng1=NGAP(1)
C old MLB 18/3/97 CALL GAUSSPWB(ALIM(1),BLIM(1),ng1,-1,W1,D1,10,ERROR,*9999)
      CALL GAUSSLEG(ALIM(1),BLIM(1),D1,W1,ng1,ERROR,*9999)
      ng2=NGAP(2)
C old MLB 18/3/97 CALL GAUSSPWB(ALIM(2),BLIM(2),ng2,-1,W2,D2,10,ERROR,*9999)
      CALL GAUSSLEG(ALIM(2),BLIM(2),D2,W2,ng2,ERROR,*9999)
      IF(DOP)THEN
        DO i=1,16
          TEST(i)=0.0d0
        ENDDO
        SUM=0.0d0
        SUML=0.0d0
        SUMC=0.0d0
        SUMPG1=0.0d0
        SUMPG2=0.0d0
        SUMPG3=0.0d0
      ENDIF
      DO j=1,ng2
        DO i=1,ng1
          XIGG(i,j,1)=D1(i)
          XIGG(i,j,2)=D2(j)
          ng=i+(j-1)*ng1
          WG(ng)=W1(i)*W2(j)
          DO ni=1,NIT(nb)
            XI(ni)=XIGG(i,j,ni)
            XIG(ni,ng)=XI(ni)
          ENDDO !ni
          ns=0
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              ns=ns+1
              DO nu=1,NUT(nb)
                IF(HERMSECTOR) THEN
                  PG(ns,nu,ng)=PSI2_HERMITE(IDO,INP,nb,nu,nk,nn,XI)
                ELSE
                  PG(ns,nu,ng)=PSI5(IBT,IDO,INP,nb,nu,nk,nn,XI)
                ENDIF
              ENDDO !nu
            ENDDO !nk
          ENDDO !nn
          IF(DOP)THEN
            SUML=SUML+1.0d0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        (1.0d0-XIG(1,ng))*(1.0d0-XIG(2,ng))*WG(ng)
            SUM=SUM+WG(ng)
            SUMPG1=SUMPG1+PG(1,1,ng)*WG(ng)
            SUMPG2=SUMPG2+PG(2,1,ng)*WG(ng)
            SUMPG3=SUMPG3+PG(3,1,ng)*WG(ng)
            TEST(1)=TEST(1)+1.0d0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        (1.0d0-XIG(1,ng))*(1.0d0-XIG(2,ng))*WG(ng)
            TEST(2)=TEST(2)+1.0d0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        (1.0d0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(3)=TEST(3)+1.0d0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        XIG(1,ng)*(1.0d0-XIG(2,ng))*WG(ng)
            TEST(4)=TEST(4)+1.0d0/(DSQRT(XIG(1,ng)**2+XIG(2,ng)**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
            TEST(5)=TEST(5)+1.0d0/(DSQRT(XIG(1,ng)**2+
     '        (1.0d0-XIG(2,ng))**2))*(1.0d0-XIG(1,ng))*
     '        (1.0d0-XIG(2,ng))*WG(ng)
            TEST(6)=TEST(6)+1.0d0/(DSQRT(XIG(1,ng)**2+
     '        (1.0d0-XIG(2,ng))**2))*
     '        (1.0d0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(7)=TEST(7)+1.0d0/(DSQRT(XIG(1,ng)**2+
     '        (1.0d0-XIG(2,ng))**2))*
     '        XIG(1,ng)*(1.0d0-XIG(2,ng))*WG(ng)
            TEST(8)=TEST(8)+1.0d0/(DSQRT(XIG(1,ng)**2+
     '        (1.0d0-XIG(2,ng))**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
            TEST(9)=TEST(9)+1.0d0/(DSQRT(XIG(2,ng)**2+
     '        (1.0d0-XIG(1,ng))**2))*
     '        (1.0d0-XIG(1,ng))*(1.0d0-XIG(2,ng))*WG(ng)
            TEST(10)=TEST(10)+1.0d0/(DSQRT(XIG(2,ng)**2+
     '        (1.0d0-XIG(1,ng))**2))*
     '        (1.0d0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(11)=TEST(11)+1.0d0/(DSQRT(XIG(2,ng)**2+
     '        (1.0d0-XIG(1,ng))**2))*
     '        XIG(1,ng)*(1.0d0-XIG(2,ng))*WG(ng)
            TEST(12)=TEST(12)+1.0d0/(DSQRT(XIG(2,ng)**2+
     '        (1.0d0-XIG(1,ng))**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
            TEST(13)=TEST(13)+1.0d0/(DSQRT((1.0d0-XIG(1,ng))**2+
     '        (1.0d0-XIG(2,ng))**2))*
     '        (1.0d0-XIG(1,ng))*(1.0d0-XIG(2,ng))*WG(ng)
            TEST(14)=TEST(14)+1.0d0/(DSQRT((1.0d0-XIG(1,ng))**2+
     '        (1.0d0-XIG(2,ng))**2))*
     '        (1.0d0-XIG(1,ng))*XIG(2,ng)*WG(ng)
            TEST(15)=TEST(15)+1.0d0/(DSQRT((1.0d0-XIG(1,ng))**2+
     '        (1.0d0-XIG(2,ng))**2))*
     '        XIG(1,ng)*(1.0d0-XIG(2,ng))*WG(ng)
            TEST(16)=TEST(16)+1.0d0/(DSQRT((1.0d0-XIG(1,ng))**2+
     '        (1.0d0-XIG(2,ng))**2))*
     '        XIG(1,ng)*XIG(2,ng)*WG(ng)
           ENDIF
        ENDDO !i
      ENDDO !j

C Update NGAP to indicate the acutal number of Gauss points used.

      NGAP(1)=ng1
      NGAP(2)=ng2

      IF(DOP)THEN
        WRITE(OP_STRING,'('' Sum of weights = '',D12.4)') SUM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' 1/(4*pi*R)*(1-X)*(1-Y) = '',D12.4)')
     '    SUML/(4.0d0*4.0d0*DATAN(1.0d0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Sum of function = '',D12.4)') SUMC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Sum of basis function 1 = '',D12.4)')
     '    SUMPG1
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Sum of basis function 2 = '',D12.4)')
     '    SUMPG2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Sum of basis function 3 = '',D12.4)')
     '    SUMPG3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,16
          WRITE(OP_STRING,'('' TEST('',I2,'') = '',D12.4)') i,
     '      TEST(i)/(4.0d0*4.0d0*DATAN(1.0d0))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !i
      ENDIF

      CALL EXITS('GAUSS12')
      RETURN
 9999 CALL ERRORS('GAUSS12',ERROR)
      CALL EXITS('GAUSS12')
      RETURN 1
      END


C MLB 18/3/97 Replacement quadrature routine

