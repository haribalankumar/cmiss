      SUBROUTINE WALL_VEC_DEF(IBT,IDO,INP,NAN,NBH,ng,
     '  NHE,NITB,nr,nx,WD1_VECTOR,WD2_VECTOR,WD3_VECTOR,PG,XI,
     '  ZE,ZG,ERROR,*)

C#### Subroutine: WALL_VEC_DEF
C###  Description:
C###    WALL_VEC_DEF calculates direction cosines of deformed
C###    wall vectors at Gauss point ng or at XI if ng=0.
C###    WD1_VECTOR coincides with the circumferential base vector
C###    WD2_VECTOR lies in the def Xi1-Xi2 plane, normal to WD1_VECTOR
C###    WD3_VECTOR is normal to deformed Xi1-Xi2 plane
C###    This routine assumes that ZE contains element vertex coords.

C     MPN 16-Jun-97: written for LIST STRAIN/STRESS output
C     THIS ROUTINE HAS NOT BEEN FULLY DEBUGGED.
C     MPN Feb2004 now works for rc coords

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM),ng,NHE,NITB,nr,nx
      REAL*8 WD1_VECTOR(3),WD2_VECTOR(3),WD3_VECTOR(3),
     '  PG(NSM,NUM,NGM,NBM),XI(3),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,mhx,mhx_theta,ni,nhx,NU1(0:3)
      REAL*8 DXIX(3,3),dZrc_dXi(3,3),dZrc_dZ(3,3),dZrc_dZref,DZX,
     '  SIGN,SUM,Z(3)

      DATA NU1/1,2,4,7/

      CALL ENTERS('WALL_VEC_DEF',*9999)

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
        WD1_VECTOR(nhx)=0.0d0
        WD2_VECTOR(nhx)=0.0d0
        WD3_VECTOR(nhx)=0.0d0
        DO ni=1,3
          dZrc_dXi(nhx,ni)=0.0d0
          dZrc_dZ(nhx,ni)=0.0d0
        ENDDO !ni
      ENDDO !nhx

C     Compute derivatives of deformed rc coords wrt deformed ref coords
C         and derivatives of deformed rc coords wrt Xi coords.
      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
        DO ni=1,NITB
          SUM=0.0d0
          DO mhx=1,NJ_LOC(NJL_GEOM,0,nr)
            IF(ni.EQ.1) THEN
              dZrc_dZref=DZX(ITYP11(nr),nhx,mhx,Z)
              dZrc_dZ(nhx,mhx)=dZrc_dZref
            ENDIF
            SUM=SUM+dZrc_dZ(nhx,mhx)*ZG(mhx,NU1(ni))
          ENDDO !mhx
          dZrc_dXi(nhx,ni)=SUM
        ENDDO !ni
      ENDDO !nhx

C     WD1_VECTOR is the normalised deformed circumferential base vector
C     MPN Feb2004 now works for rc coords
      IF(ITYP11(nr).EQ.1) THEN !rc
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          WD1_VECTOR(nhx)=dZrc_dXi(nhx,1)
        ENDDO !nhx
C!!!    MPN March2010, following comments by me stuff up, does not work
C!!!    for rc with long-axis = X   
C!!!    MPN Feb2004: long axis of LV is assumed to be the z-axis and 
C!!!    circumferential dirn is in X-Y plane (thus set Z component to zero)
C!!!    This could eventually be generalised for a specified long axis
C        WD1_VECTOR(3)=0.0d0
      ELSE IF(ITYP11(nr).EQ.2.OR.ITYP11(nr).EQ.3) THEN !cylindrical/spherical polar
        mhx_theta=2
        SIGN=1.0d0
      ELSE IF(ITYP11(nr).GE.4) THEN !prolate/oblate sph.
        mhx_theta=3
        SIGN=-1.0d0 !compute circumf dirn in same quadrant as Xi1
      ENDIF
      IF(ITYP11(nr).NE.1) THEN !not rc
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          WD1_VECTOR(nhx)=dZrc_dZ(nhx,mhx_theta)*SIGN
        ENDDO !nhx
      ENDIF
      CALL NORMALISE(3,WD1_VECTOR,ERROR,*9999)

      IF(NITB.GE.2) THEN !2D or 3D
C       WD3_VECTOR is the deformed Xi1-Xi2 plane normal
        CALL CROSS(dZrc_dXi(1,1),dZrc_dXi(1,2),WD3_VECTOR)
        CALL NORMALISE(3,WD3_VECTOR,ERROR,*9999)
C       WD2_VECTOR lies in the deformed Xi1-Xi2 plane and is
C       normal to WD1_VECTOR
        CALL CROSS(WD3_VECTOR,WD1_VECTOR,WD2_VECTOR)
C        not necessary to normalise WD2_VECTOR, since WD3_VECTOR and
C        WD1_VECTOR have already been normalised
C        CALL NORMALISE(3,WD2_VECTOR,ERROR,*9999)
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Deformed wall vectors:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''   wd1_vector  wd2_vector  wd3_vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO ni=1,3
          WRITE(OP_STRING,'(1X,3D12.3)')
     '      WD1_VECTOR(ni),WD2_VECTOR(ni),WD3_VECTOR(ni)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ni
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('WALL_VEC_DEF')
      RETURN
 9999 CALL ERRORS('WALL_VEC_DEF',ERROR)
      CALL EXITS('WALL_VEC_DEF')
      RETURN 1
      END


