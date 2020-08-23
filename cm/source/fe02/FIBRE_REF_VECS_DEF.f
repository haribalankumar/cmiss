      SUBROUTINE FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH,ng,
     '  NHE,NITB,nr,nx,FD_VECTOR,GD_VECTOR,HD_VECTOR,PG,XI,
     '  ZE,ZG,ERROR,*)

C#### Subroutine: FIBRE_REF_VECS_DEF
C###  Description:
C###    FIBRE_REF_VECS_DEF calculates direction cosines of deformed
C###    fibre reference vectors at Gauss point ng or at XI if ng=0.
C###    FD_VECTOR coincides with the deformed Xi1 base vector
C###    GD_VECTOR lies in the def Xi1-Xi2 plane and is normal to Xi1
C###    HD_VECTOR is normal to deformed Xi1-Xi2 plane
C###    This routine assumes that ZE contains element vertex
C###    coordinates and material axis angles.

C     MPN 25-Apr-96: written for fibre/imbric/sheet angle reference

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),ng,NHE,NITB,nr,nx
      REAL*8 FD_VECTOR(3),GD_VECTOR(3),HD_VECTOR(3),PG(NSM,NUM,NGM,NBM),
     '  XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,mhx,ni,nhx,NU1(0:3)
      REAL*8 DXIX(3,3),dZRC_dXI(3,3),dZrc_dZref,DZX,SUM,Z(3)
      LOGICAL ZERO_FD_VECTOR

      DATA NU1/1,2,4,7/

      CALL ENTERS('FIBRE_REF_VECS_DEF',*9999)

      DO i=1,3
        DO j=1,3
          DXIX(i,j)=0.0d0
        ENDDO
      ENDDO

C     Interpolate dependent var.s ZG and derivs wrt Xi (JP=0)
      IF(ng.EQ.0) THEN
        CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG,XI,
     '    ERROR,*9999)
      ELSE
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
      ENDIF

C     Put def coords into Z for DZX function call below
C     and initialise all vectors
      DO nhx=1,3
        IF(nhx.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
          Z(nhx)=ZG(nhx,1)
        ELSE
          Z(nhx)=0.0d0
        ENDIF
        FD_VECTOR(nhx)=0.0d0
        GD_VECTOR(nhx)=0.0d0
        HD_VECTOR(nhx)=0.0d0
        DO ni=1,3
          dZRC_dXI(nhx,ni)=0.0d0
        ENDDO !ni
      ENDDO !nhx

C     Compute derivatives of deformed rc coords wrt Xi coords
      DO ni=1,NITB
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          SUM=0.0d0
          DO mhx=1,NJ_LOC(NJL_GEOM,0,nr)
            dZrc_dZref=DZX(ITYP11(nr),nhx,mhx,Z)
            SUM=SUM+dZrc_dZref*ZG(mhx,NU1(ni))
          ENDDO !mhx
          dZRC_dXI(nhx,ni)=SUM
        ENDDO !nhx
      ENDDO !ni

C     FD_VECTOR is the normalised deformed Xi1 base vector
      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
        FD_VECTOR(nhx)=dZRC_dXI(nhx,1)
      ENDDO !nhx
      CALL NORMALISE(3,FD_VECTOR,ERROR,*9999)

      IF(NITB.GE.2) THEN !2D or 3D
C news MPN 29June2000: check for zero length F_VECTOR
C                      and correct if necesary
        ZERO_FD_VECTOR=.TRUE.
        DO nhx=1,3
          IF(DABS(FD_VECTOR(nhx)).GT.ZERO_TOL) ZERO_FD_VECTOR=.FALSE.
        ENDDO !nhx
        IF(ZERO_FD_VECTOR) THEN
C         ...so set GD_VECTOR to be the normalised def Xi2 base vector
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
            GD_VECTOR(nhx)=dZRC_dXI(nhx,2)
          ENDDO !nhx
          CALL NORMALISE(3,GD_VECTOR,ERROR,*9999)
C         ...then FD_VECTOR is the deformed Xi2-Xi3 plane normal
          CALL CROSS(dZRC_dXI(1,2),dZRC_dXI(1,3),FD_VECTOR)
          CALL NORMALISE(3,FD_VECTOR,ERROR,*9999)
C         ...and HD_VECTOR lies in the deformed Xi2-Xi3 plane and is
C         normal to GD_VECTOR
          CALL CROSS(FD_VECTOR,GD_VECTOR,HD_VECTOR)
        ELSE !FD_VECTOR is not all zero
C         HD_VECTOR is the deformed Xi1-Xi2 plane normal
          CALL CROSS(dZRC_dXI(1,1),dZRC_dXI(1,2),HD_VECTOR)
          CALL NORMALISE(3,HD_VECTOR,ERROR,*9999)
C         GD_VECTOR lies in the deformed Xi1-Xi2 plane and is
C         normal to FD_VECTOR
          CALL CROSS(HD_VECTOR,FD_VECTOR,GD_VECTOR)
        ENDIF
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Deformed fibre reference vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    fd_vector   gd_vector   hd_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO ni=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      FD_VECTOR(ni),GD_VECTOR(ni),HD_VECTOR(ni)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('FIBRE_REF_VECS_DEF')
      RETURN
 9999 CALL ERRORS('FIBRE_REF_VECS_DEF',ERROR)
      CALL EXITS('FIBRE_REF_VECS_DEF')
      RETURN 1
      END


