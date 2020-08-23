      SUBROUTINE OPEG55(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  EG,PG,RGX,XE,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: OPEG55
C###  Description:
C###    OPEG55 does string theory calculation of Green's strains EG
C###    wrt 'Fibre' coords at Gauss point ng or Xi position XI in
C###    current element.

C**** AZ,AZL,AZU  are deformed metric tensors wrt (undeformed) COORDS
C**** XG   are undeformed theta coords and derivs wrt Xi
C**** ZG   are deformed theta coords and derivs wrt undeformed COORDS
C**** EG   are physical cpts of Green's strain

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,nr,nx
      REAL*8 EG(3,3),PG(NSM,NUM,NGM,NBM),RGX,XE(NSM,NJM),
     '  XG(NJM,NUM),XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NITB
      REAL*8 AZ,AZL(3,3),AZU(3,3),DXIXN(3,3),EXR,GXL(3,3),GXU(3,3)

      CALL ENTERS('OPEG55',*9999)
      NITB=NIT(NBH(NH_LOC(1,nx)))

      IF(ng.EQ.0) THEN !Xi point
C       Interpolate Xi-coord geometric var.s XG and derivs wrt Xi
        CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
      ELSE             !Gauss points
C       Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
      ENDIF

C     Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C     derivatives of Xi wrt undeformed Nu coords, DXIXN (JP=1)
      CALL XGMG(1,NITB,NBJ(1),nr,DXIXN,GXL,GXU,RGX,XG,ERROR,*9999)

      IF(ng.EQ.0) THEN !Xi point
C       Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
        CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXN,
     '    ZE,ZG,XI,ERROR,*9999)
      ELSE             !Gauss points
C       Interpolate dependent var.s ZG and derivs wrt Nu (JP=1)
        CALL ZEZG(1,NBH,ng,NHE,nx,DXIXN,PG,ZE,ZG,ERROR,*9999)
      ENDIF

C     Calculate deformed metric tensors wrt Nu (AZL,AZU)
      CALL ZGMG(NBH(NH_LOC(1,nx)),nr,AZ,AZL,AZU,ZG,ERROR,*9999)

C     Calculate Green strain and extension ratio
      EG(1,1)=0.5d0*(AZL(1,1)-1.0d0)
      EXR=DSQRT(AZL(1,1))
      WRITE(OP_STRING,'('' AZL(1,1) ='',D12.4,'' EG(1,1) ='',D12.4,'
     '  //''' Extension ratio ='',D12.4)') AZL(1,1),EG(1,1),EXR
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      CALL EXITS('OPEG55')
      RETURN
 9999 CALL ERRORS('OPEG55',ERROR)
      CALL EXITS('OPEG55')
      RETURN 1
      END


