      SUBROUTINE DEFMGRADRC(IBT,IDO,INP,NAN,NBH,NBJ,ng,NHE,nr,nx,
     '  DZDX,PG,XG,XI,ZE,ZG,ERROR,*)

C#### Subroutine: DEFMGRADRC
C###  Description:
C###    DEFMGRADRC calculates components of the deformation
C###    gradient tensor wrt rectangular cartesian coords.
C###    This routines assumes that XG contains the Gauss pt position,
C###    interpolated material axis orientations, and derivatives
C###    with respect to Xi coordinates, and that ZE contains deformed
C###    element vertex coordinates.

C     written MPN 25-Apr-96:

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
      REAL*8 DZDX(3,3),PG(NSM,NUM,NGM,NBM),XG(NJM,NUM),XI(3),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mhx,mj,nhx,ni,NITB,nj,NU1(0:3)
      REAL*8 DETERM,DXIXJ(3,3),DXJXI(3,3),dXref_dXrc,DXZ,dZrc_dZref,
     '  DZX,SUM,X(3),Z(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('DEFMGRADRC',*9999)
      NITB=NIT(NBJ(1))

C     Calculate derivatives of Xi wrt Xj (reference) coords, DXIXJ
      DO ni=1,NITB
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          DXJXI(nj,ni)=XG(nj,NU1(ni))
        ENDDO !nj
      ENDDO !ni
      CALL INVERT(NITB,DXJXI,DXIXJ,DETERM)

      IF(ng.EQ.0) THEN
C       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
        CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIXJ,ZE,ZG,XI,
     '    ERROR,*9999)
      ELSE
C       Interpolate dependent var.s ZG and derivs wrt Xj (JP=1)
        CALL ZEZG(1,NBH,ng,NHE,nx,DXIXJ,PG,ZE,ZG,ERROR,*9999)
      ENDIF

C     Put undef/def coords into X/Z for DZX/DXZ function calls below
      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
        X(nhx)=XG(nhx,1)
        Z(nhx)=ZG(nhx,1)
      ENDDO !nhx

      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          SUM=0.0d0
          DO mhx=1,NJ_LOC(NJL_GEOM,0,nr)
            dZrc_dZref=DZX(ITYP11(nr),nhx,mhx,Z)
            DO mj=1,NJ_LOC(NJL_GEOM,0,nr)
              dXref_dXrc=DXZ(ITYP10(nr),mj,nj,X)
              SUM=SUM+dZrc_dZref*ZG(mhx,NU1(mj))*dXref_dXrc
            ENDDO !mj
          ENDDO !mhx
          DZDX(nhx,nj)=SUM
        ENDDO !nj
      ENDDO !nhx

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          WRITE(OP_STRING,'('' DZDX('',I1,'',nj)   : '',3D12.4)')
     '      nhx,(DZDX(nhx,nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nhx
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('DEFMGRADRC')
      RETURN
 9999 CALL ERRORS('DEFMGRADRC',ERROR)
      CALL EXITS('DEFMGRADRC')
      RETURN 1
      END

