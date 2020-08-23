      SUBROUTINE DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,COORDS,ERROR,*)

C#### Subroutine: DXIDZM
C###  Description:
C###    DXIDZM evaluates derivatives (DXIZN) of Xi- wrt
C###    deformed Nu(fibre)-coords     if COORDS='Fibre' or
C###    deformed Wall(cardiac)-coords if COORDS='Wall'
C###    (and their inverse DZNXI), at Gauss point ng or at XI if ng=0.
C###    This routine assumes that XE/ZE contain element vertex
C###    coordinates and microstructural orientations for the
C###    undeformed/deformed state resp.

C Written MPN 18-Apr-96
C     If COORDS='Fibre'
C       MAT_VEC_DEF is used to calculate the rectangular cartesian
C       components of the deformed anatomical material vectors:
C         DZDNU(k,1) is the deformed fibre        vector in rc coords
C         DZDNU(k,2) is the deformed sheet        vector in rc coords
C         DZDNU(k,3) is the deformed sheet-normal vector in rc coords
C     If COORDS='Wall'
C       WALL_VEC_DEF is used to calculate the rectangular cartesian
C       components of the deformed (cardiac) wall vectors:
C         DZDNU(k,1) is the deformed circumferential vector in rc coords
C         DZDNU(k,2) is the deformed longitudunal    vector in rc coords
C         DZDNU(k,3) is the deformed radial          vector in rc coords

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),ng,NHE,nr,nx
      REAL*8 DXIZN(3,3),DZNXI(3,3),PG(NSM,NUM,NGM,NBM),
     '  XE(NSM,NJM),XG(NJM,NUM),XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER COORDS*(*),ERROR*(*)
!     Local Variables
      INTEGER i,j,mi,mhx,ni,ni2,nhx,NITB,NU1(0:3)
      REAL*8 DETERM,DXIX(3,3),DZDNU(3,3),dZrc_dZref,DZX,
     '  GZ,GZL(3,3),GZU(3,3),SUM,Z(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('DXIDZM',*9999)
      NITB=NIT(NBJ(1))

      DO i=1,3
        DO j=1,3
          DXIX(i,j)=0.0d0
        ENDDO
      ENDDO

      IF(NITB.LT.3) THEN
C       Initialise DXIZN (not needed for NITB=3 since all elements
C       of DXIZN are set).
        DO ni=1,3
          DO mi=1,3
            DXIZN(ni,mi)=0.0d0
            DZNXI(ni,mi)=0.0d0
            IF(mi.EQ.ni) THEN
              DXIZN(ni,mi)=1.0d0
              DZNXI(ni,mi)=1.0d0
            ENDIF
            IF(mi.EQ.ni) DXIZN(ni,mi)=1.0d0
          ENDDO !mi
        ENDDO !ni
      ENDIF

      IF(COORDS(1:5).EQ.'Fibre') THEN
C       Compute deformed anatomical fibre vectors wrt rc coordinates
        CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '    DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '    PG,XE,XG,XI,ZE,ZG,.TRUE.,ERROR,*9999)
      ELSE IF(COORDS(1:4).EQ.'Wall') THEN
C       Compute deformed (cardiac) wall vectors wrt rc coordinates
        CALL WALL_VEC_DEF(IBT,IDO,INP,NAN,NBH,ng,NHE,NITB,nr,nx,
     '    DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),PG,XI,ZE,ZG,ERROR,*9999)
      ENDIF

      IF(ng.EQ.0) THEN
C ***   Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
        CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG,XI,
     '    ERROR,*9999)
      ELSE
C ***   Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
      ENDIF
C *** Calculate deformed metric tensors wrt Xi (GZL,GZU)
      CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZG,ERROR,*9999)

C     Put def coords into Z for DZX function call below
      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
        Z(nhx)=ZG(nhx,1)
      ENDDO !nhx

C     Calc derivs of Xi wrt deformed Nu/Wall coords
      DO ni=1,NITB
        DO mi=1,NITB
          SUM=0.0d0
          DO ni2=1,NITB
            DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
              DO mhx=1,NJ_LOC(NJL_GEOM,0,nr)
                dZrc_dZref=DZX(ITYP11(nr),mhx,nhx,Z)
                SUM=SUM+GZU(ni,ni2)*dZrc_dZref*ZG(nhx,NU1(ni2))*
     '            DZDNU(mhx,mi)
              ENDDO !mhx
            ENDDO !nhx
          ENDDO !ni2
          DXIZN(ni,mi)=SUM
        ENDDO !mi
      ENDDO !ni

      CALL INVERT(NITB,DXIZN,DZNXI,DETERM)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO mi=1,NITB
          WRITE(OP_STRING,'('' DXIZN('',I1,'',ni): '',3D12.4)')
     '      mi,(DXIZN(mi,ni),ni=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('DXIDZM')
      RETURN
 9999 CALL ERRORS('DXIDZM',ERROR)
      CALL EXITS('DXIDZM')
      RETURN 1
      END


