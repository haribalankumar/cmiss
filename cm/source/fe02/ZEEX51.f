      SUBROUTINE ZEEX51(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DXIXN,DXNXI,EG,PG,PHI,PST,RGX,RM,XE,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: ZEEX51
C###  Description:
C###    ZEEX51 is a cut-down version of ZEEX50 for calculating Greens
C###    strains EG, principal strains PST and spectral matrix RM
C###    wrt 'Fibre' coords at position XI in current element.
C###    AZ,AZL,AZU  are deformed metric tensors wrt (undeformed) COORDS

C**** XG   are undeformed theta coords and derivs wrt Xi
C**** ZG   are deformed theta coords and derivs wrt undeformed COORDS
C**** EG   are physical cpts of Green's strain
C**** PST  are principal strains
C**** RM   is the spectral matrix whose cols are the eigenvectors

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,nr,nx
      REAL*8 DXIXN(3,3),DXNXI(3,3),EG(3,3),PG(NSM,NUM,NGM,NBM),PHI(3),
     '  PST(3),RM(3,3),RGX,XE(NSM,NJM),XG(NJM,NUM),
     '  XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IFAIL,mi,nb,ni,NITB
      REAL*8 AZ,AZL(3,3),AZU(3,3),DETERM,DXIXJ(3,3),GXL(3,3),GXU(3,3),
     '  TOL,VALTMP,WK1_LOCAL(10)

      DATA TOL /1.0D-12/

      CALL ENTERS('ZEEX51',*9999)

      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)

      IF(ng.EQ.0) THEN
C ***   Interpolate Xi-coord geometric var.s XG and derivs wrt Xi
        CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
      ELSE
C ***   Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
      ENDIF
C *** Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C *** ..derivatives of Xi wrt Xj (reference) coords, DXIXJ (IP=0)
      CALL XGMG(0,NITB,NBJ(1),nr,DXIXJ,GXL,GXU,RGX,XG,ERROR,*9999)
C *** Get derivs of Xi wrt undeformed Nu (body/fibre) coords,DXIXN
      CALL DXIDXM(NITB,nr,DXIXN,DETERM,XG,'Fibre',ERROR,*9999)
C *** Calculate DXNXI from inverse of DXIDN although this doesn't seem
C *** ..to get used anywhere
      CALL INVERT(NITB,DXIXN,DXNXI,DETERM)
      IF(ng.EQ.0) THEN
C ***   Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
        CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,
     '    DXIXN,ZE,ZG,XI,ERROR,*9999)
      ELSE
C ***   Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
        CALL ZEZG(1,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
      ENDIF
C *** Calculate deformed metric tensors wrt Nu (AZL,AZU)
      CALL ZGMG(nb,nr,AZ,AZL,AZU,ZG,ERROR,*9999)

C!!! If using initial extns for residual strains then need to pass mat
C!!! param arrays and uncomment next section of code (see ZEEX50).

CC!!! NOTE can only use resid strains for pole zero law until init extn
CC!!!      mat params are set up for other problem types
C      IF(KTYP56(nr).EQ.3) THEN !pole zero law
CC       Calc growth defm tens for resid strain and copy AZL to AZL_tmp
C        DO i=1,3
C          DO j=1,3
C            Fgrowth(i,j)=0.0d0 !growth defm tensor for resid strains
C            AZL_tmp(i,j)=AZL(i,j)
C            DZDX_tmp(i,j)=DZDX(i,j)
C          ENDDO !j
C          Fgrowth(i,i)=CG(27+i) !NOTE init extns CG(28),CG(29),CG(30)
C        ENDDO !i
CC       Apply growth defm to deformed covariant metric tensor AZL
C        DO i=1,3
C          DO j=1,3
C            AZL(i,j)=0.0d0
C            DZDX(i,j)=0.0d0
C            DO m=1,3
C              DO n=1,3
C                AZL(i,j)=AZL(i,j)+
C     '            Fgrowth(i,m)*AZL_tmp(m,n)*Fgrowth(n,j)
C              ENDDO !n
C              DZDX(i,j)=DZDX(i,j)+DZDX_tmp(i,m)*Fgrowth(m,j)
C            ENDDO !m
C          ENDDO !j
C        ENDDO !i
CC       Recompute AZU,AZ from transformed AZL
C        CALL INVERT(NIT(nb),AZL,AZU,AZ)
C      ENDIF !pole zero law

      DO mi=1,NITB
        DO ni=1,NITB
C ***     Calculate physical shear cpts of Green's strain
          EG(mi,ni)=0.5d0*AZL(mi,ni)
        ENDDO
C ***   Calculate diagonal physical cpts of Green's strain
        EG(mi,mi)=0.5d0*(AZL(mi,mi)-1.0D0)
      ENDDO
      IFAIL=0
C      CALL F02ABF(EG,3,NITB,PST,RM,3,WK1_LOCAL,IFAIL)
      DO ni=1,3
        DO mi=1,3
          RM(ni,mi)=EG(ni,mi)
        ENDDO
      ENDDO
C MLB 19/3/97
C This may not give evectors as accurately as NAG
      CALL DSYEV('V','L',NITB,RM,3,PST,WK1_LOCAL,10,IFAIL)
!news MPN 30-Jun-94
C *** Order PST from max->min and change cols of RM accordingly
      IF(PST(2).GT.PST(1)) THEN
        VALTMP=PST(1)
        PST(1)=PST(2)
        PST(2)=VALTMP
        DO ni=1,NITB
          VALTMP=RM(ni,1)
          RM(ni,1)=RM(ni,2)
          RM(ni,2)=VALTMP
        ENDDO
      ENDIF
      IF(NITB.EQ.3) THEN
        DO mi=1,2
          IF(PST(3).GT.PST(mi)) THEN
            VALTMP=PST(mi)
            PST(mi)=PST(3)
            PST(3)=VALTMP
            DO ni=1,NITB
              VALTMP=RM(ni,mi)
              RM(ni,mi)=RM(ni,3)
              RM(ni,3)=VALTMP
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!newe
      IF(NITB.EQ.2) THEN
        IF(DABS(RM(1,1)).GT.TOL.OR.DABS(RM(2,1)).GT.TOL) THEN
          PHI(1)=DATAN2(RM(2,1),RM(1,1))
        ELSE
          PHI(1)=0.0D0
        ENDIF
        IF(DABS(RM(1,2)).GT.TOL.OR.DABS(RM(2,2)).GT.TOL) THEN
          PHI(2)=DATAN2(RM(2,2),RM(1,2))
        ELSE
          PHI(2)=0.0D0
        ENDIF
      ELSE IF(NITB.EQ.3) THEN
        IF(DABS(RM(1,1)).GT.TOL.OR.DABS(RM(2,1)).GT.TOL) THEN
          PHI(1)=DATAN2(RM(2,1),RM(1,1))
        ELSE IF(RM(2,1)*RM(1,1).GT.0.0D0) THEN
          PHI(1)=90.0D0
        ELSE
          PHI(1)=-90.0D0
        ENDIF
        IF(DABS(RM(3,1)).LE.1.0D0) THEN
          PHI(2)=DASIN(RM(3,1))
        ELSE IF(RM(3,1).GT.1.0D0) THEN
          PHI(2)=90.0D0
        ELSE
          PHI(2)=-90.0D0
        ENDIF
         IF(DABS(DCOS(PHI(1))).GT.TOL .AND.
     '     DABS(DCOS(PHI(1))).GE.DABS(RM(3,3))) THEN
           PHI(3)=DACOS(RM(3,3)/DCOS(PHI(1)))
        ELSE
          PHI(3)=0.0D0
        ENDIF
      ENDIF

      CALL EXITS('ZEEX51')
      RETURN
 9999 CALL ERRORS('ZEEX51',ERROR)
      CALL EXITS('ZEEX51')
      RETURN 1
      END


