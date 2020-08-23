      SUBROUTINE DLZJDX(ICOORD,nb,nr,DZJDX,XG,ZG,ERROR,*)

C#### Subroutine: DLZJDX
C###  Description:
C###    DLZJDX calculates the covariant derivatives DZJDX of deformed
C###    theta wrt undeformed theta (ie the cpts wrt theta of the
C###    deformation gradient tensor) or wrt undeformed Nu, or wrt Xi
C###    depending on derivatives contained in ZG.

C**** ICOORD=1..5 is coordinate system: Rectangular Cartesian,
C****   Cylindrical Polar, Spherical Polar, Prolate Spheroidal,
C****   Oblate Spheroidal.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ICOORD,nb,nr
      REAL*8 DZJDX(3,3),XG(NJM,NUM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,mj,ni,NITB,NU1(0:3)
      REAL*8 CC,CCL,CLX,CLZ,CMX,CMZ,CS,CSL,CT,CX,CZ,DLB,DMB,
     '  DT,DTB,G1,G3,PB,RB,RX,RZ,SC,SCL,SLX,SLZ,SMX,SMZ,SS,SSL,
     '  ST,SX,SZ,TB

      DATA NU1/1,2,4,7/

      CALL ENTERS('DLZJDX',*9999)
      NITB=NIT(nb)
      IF(ICOORD.EQ.1) THEN
        DO mj=1,NJ_LOC(NJL_GEOM,0,nr)
          DO ni=1,NITB
            DZJDX(mj,ni)=ZG(mj,NU1(ni))
          ENDDO !ni
        ENDDO !mj

      ELSE IF(ICOORD.EQ.2) THEN
        RX=XG(1,1)
        RZ=ZG(1,1)
        DT=ZG(2,1)-XG(2,1)
        CT=DCOS(DT)
        ST=DSIN(DT)
        DO ni=1,NITB
          DZJDX(1,ni)= ZG(1,NU1(ni))*CT-RZ*ST*ZG(2,NU1(ni))
          DZJDX(2,ni)=(ZG(1,NU1(ni))*ST+RZ*CT*ZG(2,NU1(ni)))/RX
          IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) DZJDX(3,ni)= ZG(3,NU1(ni))
        ENDDO !ni

      ELSE IF(ICOORD.EQ.3) THEN
        RX=XG(1,1)
        RZ=ZG(1,1)
        DT=ZG(2,1)-XG(2,1)
        CT=DCOS(DT)
        ST=DSIN(DT)
        CX=DCOS(XG(3,1))
        SX=DSIN(XG(3,1))
        CZ=DCOS(ZG(3,1))
        SZ=DSIN(ZG(3,1))
        CC=CX*CZ
        SS=SX*SZ
        CS=CX*SZ
        SC=SX*CZ
        DO ni=1,NITB
          RB=ZG(1,NU1(ni))
          TB=ZG(2,NU1(ni))
          PB=ZG(3,NU1(ni))
          DZJDX(1,ni)=CC*(RB*CT-RZ*ST*TB)-CS*RZ*CT*PB+SC*RZ*PB+SS*RB
          DZJDX(2,ni)=(RZ*CZ*CT*TB+(RB*CZ-RZ*SZ*PB)*ST)/(RX*CX)
          DZJDX(3,ni)=(CC*RZ*PB+CS*RB+SC*(RZ*ST*TB-RB*CT)+SS*RZ*CT*PB)
     '      /RX
        ENDDO !ni

      ELSE IF(ICOORD.EQ.4) THEN
        SLX=DSINH(XG(1,1))
        SLZ=DSINH(ZG(1,1))
        SMX=DSIN(XG(2,1))
        SMZ=DSIN(ZG(2,1))
        CLX=DSQRT(1.0D0+SLX*SLX)
        CLZ=DSQRT(1.0D0+SLZ*SLZ)
        CMX=DSQRT(1.0D0-SMX*SMX)
        CMZ=DSQRT(1.0D0-SMZ*SMZ)
        DT=ZG(3,1)-XG(3,1)
        CT=DCOS(DT)
        ST=DSIN(DT)
        CCL=CLX*CLZ
        CSL=CLX*SLZ
        SCL=SLX*CLZ
        SSL=SLX*SLZ
        CC =CMX*CMZ
        CS =CMX*SMZ
        SC =SMX*CMZ
        SS =SMX*SMZ
        G1=SLX*SLX+SMX*SMX
        G3=SLX*SLX*SMX*SMX
        DO ni=1,NITB
          DLB=ZG(1,NU1(ni))
          DMB=ZG(2,NU1(ni))
          DTB=ZG(3,NU1(ni))
          DZJDX(1,ni)=(( SSL*CC+CCL*SS*CT)*DLB+(-SCL*CS+CSL*SC*CT)*DMB
     '      -CSL*SS*ST*DTB)/G1
          DZJDX(2,ni)=((-CSL*SC+SCL*CS*CT)*DLB+( CCL*SS+SSL*CC*CT)*DMB
     '      -SSL*CS*ST*DTB)/G1
          DZJDX(3,ni)=(SCL*SS*ST*DLB+SSL*SC*ST*DMB+SSL*SS*CT*DTB)/G3
        ENDDO !ni
      ELSE IF(ICOORD.EQ.5) THEN
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO mi=1,NITB
          WRITE(OP_STRING,'('' DZJDX('',I1,'',ni)   : '',3D12.4)')
     '      mi,(DZJDX(mi,ni),ni=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('DLZJDX')
      RETURN
 9999 CALL ERRORS('DLZJDX',ERROR)
      CALL EXITS('DLZJDX')
      RETURN 1
      END


