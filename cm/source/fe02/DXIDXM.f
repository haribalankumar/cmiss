      SUBROUTINE DXIDXM(NITB,nr,DXIXN,RG,XG,COORDS,ERROR,*)

C#### Subroutine: DXIDXM
C###  Description:
C###    DXIDXM evaluates derivatives (DXIXN) of Xi- wrt
C###    undeformed Nu(fibre)-coords     if COORDS='Fibre' or
C###    undeformed Wall(cardiac)-coords if COORDS='Wall'
C###    (for NIT=NJT only).
C###    RG returns the domain Jacobian.
C###    This routine assumes that XG contains
C###    the Gauss pt position, interpolated material axis
C###    orientations, and derivatives with respect to Xi coordinates.

C Rewritten MPN 18-Apr-96
C     If COORDS='Fibre'
C       MAT_VEC is used to calculate the rectangular cartesian
C       components of the undeformed anatomical material vectors:
C         DXDXN(k,1) is the undef fibre direction vector in rc coords
C         DXDXN(k,2) is the undef sheet direction vector in rc coords
C         DXDXN(k,3) is the undef sheet-normal dirn vect in rc coords
C     If COORDS='Wall'
C       WALL_VEC is used to calculate the rectangular cartesian
C       components of the undeformed (cardiac) wall vectors:
C         DXDXN(k,1) is the undef circumferential dirn vect in rc coords
C         DXDXN(k,2) is the undef longitudunal    dirn vect in rc coords
C         DXDXN(k,3) is the undef radial          dirn vect in rc coords

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NITB,nr
      REAL*8 DXIXN(3,3),RG,XG(NJM,NUM)
      CHARACTER COORDS*(*),ERROR*(*)
!     Local Variables
      INTEGER mjj,ni,nj,njj
      REAL*8 DXRCXN(3,3),DXRCXI(3,3),DXIXRC(3,3)

      CALL ENTERS('DXIDXM',*9999)

      IF(NITB.EQ.1.AND.ITYP10(nr).EQ.1) THEN
! Handle 1D as special case
C KAT 20Nov97:  This seems to be reciprocal of what we want.
C        DXIXN(1,1)=DSQRT(GL(1,1)) !arclength along 1D element
        IF(NJT.EQ.1) THEN
          nj=NJ_LOC(NJL_GEOM,1,nr)
          DXRCXI(1,1)=XG(nj,2) !arclength along 1D element
        ELSE
C CS 11/4/2003 Generalising for the 1d element in 2d or 3d space case
C Note: The variable DXRCXI is used here for consistancy with 
C       other cases
C       but we are infact calculating ds/dxi here.
          DXRCXI(1,1)=0.0d0
          DO njj=1,NJT
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            DXRCXI(1,1)=DXRCXI(1,1)+XG(nj,2)**2
          ENDDO
          DXRCXI(1,1)=DSQRT(DXRCXI(1,1))
        ENDIF
        CALL INVERT(1,DXRCXI,DXIXN,RG)
      ELSE
C       Calculate dXrc/dXI
        CALL DXRCDXI(NITB,nr,DXRCXI,XG,ERROR,*9999)

        IF(COORDS(1:5).EQ.'Fibre') THEN
C         Compute undeformed anatomical fibre vectors wrt rc coordinates
          CALL MAT_VEC(NITB,nr,DXRCXN(1,1),DXRCXN(1,2),DXRCXN(1,3),
     '      DXRCXI,XG,ERROR,*9999)
        ELSE IF(COORDS(1:4).EQ.'Wall') THEN
C         Compute undeformed wall (cardiac) vectors wrt rc coordinates
          CALL WALL_VEC(NITB,nr,DXRCXN(1,1),DXRCXN(1,2),DXRCXN(1,3),
     '      DXRCXI,XG,ERROR,*9999)
        ENDIF

C       Calculate dXI/dXrc
        CALL INVERT(NITB,DXRCXI,DXIXRC,RG)

C       Calc derivatives of Xi wrt undeformed Nu/Wall
        DO njj=1,3
          DO ni=1,3
            DXIXN(ni,njj)=0.0d0
          ENDDO !mi
        ENDDO !ni
        DO ni=NITB+1,3
          DXIXN(ni,ni)=1.0d0
        ENDDO !ni
        DO njj=1,NITB
          DO mjj=1,NITB
            DO ni=1,NITB
              DXIXN(ni,njj)=DXIXN(ni,njj)+DXIXRC(ni,mjj)*DXRCXN(mjj,njj)
            ENDDO !ni
          ENDDO !mjj
        ENDDO !njj

C KAT 19Nov97:  Old and slow:
C! Put undeformed coords into X for DZX function call below
C        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
C          X(njj)=XG(njj,1)
C        ENDDO !njj
C
C! Calc derivatives of Xi wrt undeformed Nu
C        DO ni=1,NITB
C          DO mi=1,NITB
C            SUM=0.0d0
C            DO ni2=1,NITB
C              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
C                DO mjj=1,NJ_LOC(NJL_GEOM,0,nr)
C                  dXrc_dXref=DZX(ITYP10(nr),mjj,njj,X)
C                  SUM=SUM+GU(ni,ni2)*dXrc_dXref*XG(njj,NU1(ni2))*
C     '              DXXN(mjj,mi)
C                ENDDO !mjj
C              ENDDO !njj
C            ENDDO !ni2
C            DXIXN(ni,mi)=SUM
C          ENDDO !mi
C        ENDDO !ni
      ENDIF !NITB>1

! Calculate DXNXI from inverse
C      CALL INVERT(NITB,DXIXN,DXNXI,DETERM)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO ni=1,NITB
          WRITE(OP_STRING,'('' DXIXN('',I1,'',ni): '',3D12.4)')
     '      ni,(DXIXN(ni,njj),njj=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('DXIDXM')
      RETURN
 9999 CALL ERRORS('DXIDXM',ERROR)
      CALL EXITS('DXIDXM')
      RETURN 1
      END


