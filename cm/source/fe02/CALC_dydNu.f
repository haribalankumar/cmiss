      SUBROUTINE CALC_dydNu(IBT,IDO,INP,NAN,NBJ,nr,
     '  dydNu,dydx,XE,XG,XI,ZE,ERROR,*)

C#### Subroutine: CALC_dydNu
C###  Description:
C###    CALC_dydNu evaluates deformation tensor dydNu, deformed
C###    cartesian vessel coordinates- wrt undeformed Nu(fibre)-coords
C###    at XI

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBJ(NJM),nr
      REAL*8 dxdNu(3,3),dydx(3,3),dydNu(3,3),
     '  XE(NSM,NJM),XG(NJM,NUM),
     '  XI(3),ZE(NSM,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,mjj,nb,ni,ni2,nj,njj,NITB,nu
      REAL*8 DCART_DREF(3,3),DCART_DXI(3,3),
     '  DXDXN(3,3),DXIXN(3,3),dXrc_dXref,
     '  DXXI(3,0:3),DZXI(3,3),DZX,
     '  GL(3,3),GU(3,3),PXI,SUM,X(3)

      CALL ENTERS('CALC_dydNu',*9999)

      NITB=NIT(NBJ(1))

      CALL XMG(IBT,IDO,INP,NBJ,nr,GL,GU,XE,XI,ERROR,*9999)

      DO ni=0,NITB
        nu=1+ni*(1+ni)/2
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          nb=NBJ(nj)
          DXXI(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,nu,
     '      XI,XE(1,nj))
        ENDDO !nj
      ENDDO !ni

      IF(NITB.GT.1) THEN
! Compute undeformed anatomical fibre vectors wrt rc coordinates

        CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '    DXDXN(1,1),DXDXN(1,2),DXDXN(1,3),XE,XG,XI,.TRUE.,ERROR,*9999)

! Put undeformed coords into X for DZX function call below
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          nb=NBJ(nj)
          X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '      INP(1,1,nb),nb,1,XI,XE(1,nj))
        ENDDO !njj


! Calc derivatives of Xi wrt undeformed Nu
        DO ni=1,NITB
          DO mi=1,NITB
            SUM=0.0d0
            DO ni2=1,NITB
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                DO mjj=1,NJ_LOC(NJL_GEOM,0,nr)
                  dXrc_dXref=DZX(ITYP10(nr),mjj,njj,X)
                  SUM=SUM+GU(ni,ni2)*dXrc_dXref*DXXI(njj,ni2)*
     '              DXDXN(mjj,mi)
                ENDDO !mjj
              ENDDO !njj
            ENDDO !ni2
            DXIXN(ni,mi)=SUM
          ENDDO !mi
        ENDDO !ni
      ENDIF !NITB>1

! Calculate DXNXI from inverse

      DO ni=1,NITB
        nu=1+ni*(1+ni)/2
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          nb=NBJ(nj)
          DZXI(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,nu,
     '      XI,ZE(1,nj))
        ENDDO !nj
      ENDDO !ni

C derivative of reference z wrt xi

      DO ni=1,NITB
        DO nj=1,NITB
          DCART_DREF(ni,nj)=DZX(ITYP10(nr),ni,nj,X)
        ENDDO
      ENDDO

      DO ni=1,NITB
        DO nj=1,NITB
          SUM=0.0d0
          DO ni2=1,NITB
            SUM=SUM+DCART_DREF(ni,ni2)*DZXI(ni2,nj)
          ENDDO
          DCART_DXI(ni,nj)=SUM
        ENDDO
      ENDDO

C derivative of cartesian x wrt xi

      DO ni=1,NITB
        DO mi=1,NITB
          SUM=0.0d0
          DO ni2=1,NITB
            SUM=SUM+DCART_DXI(ni,ni2)*DXIXN(ni2,mi)
          ENDDO
          dxdNu(ni,mi)=SUM
        ENDDO
      ENDDO

      DO ni=1,NITB
        DO mi=1,NITB
          SUM=0.0d0
          DO ni2=1,NITB
            SUM=SUM+dydx(ni,ni2)*dxdNu(ni2,mi)
          ENDDO
          dydNu(ni,mi)=SUM
        ENDDO
      ENDDO

      CALL EXITS('CALC_dydNu')
      RETURN
 9999 CALL ERRORS('CALC_dydNu',ERROR)
      CALL EXITS('CALC_dydNu')
      RETURN 1
      END


