      SUBROUTINE MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  AD_VECTOR,BD_VECTOR,CD_VECTOR,PG,XE,XG,XI,ZE,ZG,CALC_XG,
     '  ERROR,*)

C#### Subroutine: MAT_VEC_DEF
C###  Description:
C###    MAT_VEC_DEF calculates direction cosines of deformed
C###    normalised orthogonal material vectors at Gauss point ng or at
C###    XI if ng=0.
C###    If ng=0 this routine assumes that XE/ZE contain element vertex
C###    coordinates and microstructural orientations for the
C###    undeformed/deformed state resp.
C###    If ng>0 this routine assumes XG contains Gauss pt coordinates,
C###    microstructural orientations, and derivatives wrt Xi.

C     written MPN 25-Apr-96

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,nr,nx
      REAL*8 AD_VECTOR(3),BD_VECTOR(3),CD_VECTOR(3),PG(NSM,NUM,NGM,NBM),
     '  XE(NSM,NJM),XG(NJM,NUM),XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
      LOGICAL CALC_XG
!     Local Variables
      INTEGER mj,ni,NITB,nj
      REAL*8 DXDNU(3,3),DZDNU(3,3),DZDX(3,3),SUM

      CALL ENTERS('MAT_VEC_DEF',*9999)
      NITB=NIT(NBJ(1))

      IF(ng.EQ.0) THEN
C       Compute undeformed anatomical fibre vectors wrt rc coords at XI
        CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '    DXDNU(1,1),DXDNU(1,2),DXDNU(1,3),XE,XG,XI,CALC_XG,ERROR,*9999)
      ELSE
C       Compute undeformed anatomical fibre vectors wrt rc coords at ng
        CALL MAT_VEC_NG(NITB,nr,DXDNU(1,1),DXDNU(1,2),DXDNU(1,3),XG,
     '    ERROR,*9999)
      ENDIF
C     Initialise deformed material vectors with undef material vec dirns
      DO nj=1,3
        AD_VECTOR(nj)=DXDNU(nj,1)
        BD_VECTOR(nj)=DXDNU(nj,2)
        CD_VECTOR(nj)=DXDNU(nj,3)
      ENDDO !nj

C     Compute the deformation gradient tensor wrt rc coords
C     ie want derivatives of deformed rc coordinates wrt
C     undeformed rc coordinates.
      CALL DEFMGRADRC(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DZDX,PG,XG,XI,ZE,ZG,ERROR,*9999)

C     Compute deformed material vectors wrt rc coordinates using
C     the deformation gradient tensor and the undeformed material
C     vectors wrt rc coordinates
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        DO ni=1,NITB ! we don't actually need the third component
          SUM=0.0d0
          DO mj=1,NJ_LOC(NJL_GEOM,0,nr)
            SUM=SUM+DZDX(nj,mj)*DXDNU(mj,ni)
          ENDDO !nj1
          DZDNU(nj,ni)=SUM
        ENDDO !ni
        AD_VECTOR(nj)=DZDNU(nj,1)
      ENDDO !nj

C     Normalise the deformed vectors so that deformed material
C     coordinates are measures of physical arc length
      CALL NORMALISE(3,AD_VECTOR,ERROR,*9999)

C     Make the second vector orthogonal to the first (still within the sheet)
      SUM=0.0d0
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        SUM=SUM+AD_VECTOR(nj)*DZDNU(nj,2)
      ENDDO !nj
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        BD_VECTOR(nj)=DZDNU(nj,2)-SUM*AD_VECTOR(nj)
      ENDDO !nj
      
      CALL NORMALISE(3,BD_VECTOR,ERROR,*9999)

C     The third vector is orthogonal to the others (normal to the sheet)
      IF(NITB.EQ.3) CALL CROSS(AD_VECTOR,BD_VECTOR,CD_VECTOR)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Deformed material vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    ad_vector   bd_vector   cd_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nj=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      AD_VECTOR(nj),BD_VECTOR(nj),CD_VECTOR(nj)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nj
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('MAT_VEC_DEF')
      RETURN
 9999 CALL ERRORS('MAT_VEC_DEF',ERROR)
      CALL EXITS('MAT_VEC_DEF')
      RETURN 1
      END


