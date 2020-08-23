      SUBROUTINE ZGMG(nb,nr,GZ,GZL,GZU,ZG,ERROR,*)

C#### Subroutine: ZGMG
C###  Description:
C###    ZGMG evaluates components of metric tensor in deformed state
C###    at current Gauss point from coordinates ZG.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER nb,nr
      REAL*8 GZ,GZL(3,3),GZU(3,3),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,nhx,ni,NITB,NU1(0:3)
      REAL*8 AA,G1,G3,SLZ,SMZ,SUM
      CHARACTER CHAR1*1

      DATA NU1/1,2,4,7/

      CALL ENTERS('ZGMG',*9999)
      NITB=NIT(nb)
      DO mi=1,NITB
        DO ni=1,NITB
          SUM=ZG(1,NU1(mi))*ZG(1,NU1(ni))
          IF(NJ_LOC(NJL_GEOM,0,nr).GT.1) THEN
            IF(ITYP11(nr).EQ.1) THEN
              DO nhx=2,NJ_LOC(NJL_GEOM,0,nr)
                SUM=SUM+ZG(nhx,NU1(mi))*ZG(nhx,NU1(ni))
              ENDDO
            ELSE IF(ITYP11(nr).EQ.2) THEN
              SUM=SUM+ZG(1,1)**2*ZG(2,NU1(mi))*ZG(2,NU1(ni))
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3)
     '          SUM=SUM+ZG(3,NU1(mi))*ZG(3,NU1(ni))
            ELSE IF(ITYP11(nr).EQ.3) THEN
              SUM=SUM+ZG(1,1)**2*(DCOS(ZG(3,1))**2*ZG(2,NU1(mi))
     '               *ZG(2,NU1(ni))+ZG(3,NU1(mi))*ZG(3,NU1(ni)))
            ELSE IF(ITYP11(nr).EQ.4) THEN
              AA=FOCUS*FOCUS
              SLZ=DSINH(ZG(1,1))
              SMZ=DSIN(ZG(2,1))
              G1=AA*(SLZ*SLZ+SMZ*SMZ)
              G3=AA*SLZ*SLZ*SMZ*SMZ
              SUM=G1*(SUM+ZG(2,NU1(mi))*ZG(2,NU1(ni)))
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3)
     '          SUM=SUM+G3*ZG(3,NU1(mi))*ZG(3,NU1(ni))
            ENDIF
          ENDIF
          GZL(mi,ni)=SUM
        ENDDO
      ENDDO

C new MPN 17-Apr-96: calc GZU with call to INVERT
C     Calculate contravariant metric tensor GZU(i,j)
      CALL INVERT(NITB,GZL,GZU,GZ)
      IF(DABS(GZ).LT.RDELTA) THEN
        WRITE(OP_STRING,'('' >>Warning: zero GZ in ZGMG. GZ='',D12.5,'
     '    //''' RDELTA='',D12.5)') GZ,RDELTA
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF
C old
C      IF(NITB.EQ.1) THEN
C        GZ=GZL(1,1)
C        GZU(1,1)=1.0D0/GZ
C      ELSE IF(NITB.EQ.2) THEN
C        GZ=GZL(1,1)*GZL(2,2)-GZL(1,2)*GZL(2,1)
C        GZU(1,1)= GZL(2,2)/GZ
C        GZU(1,2)=-GZL(1,2)/GZ
C        GZU(2,1)=-GZL(2,1)/GZ
C        GZU(2,2)= GZL(1,1)/GZ
C        GZU(3,3)=1.0d0
C      ELSE IF(NITB.EQ.3) THEN
C        CALL INVERT(3,GZL,GZU,GZ)
C      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(CHAR1,'(I1)') NITB
        WRITE(OP_STRING,'(''  GZL or AZL:'','//CHAR1(1:1)//'D12.4,'
     '    //'''  GZU or AZU:'','//CHAR1(1:1)//'D12.4,'
     '    //'/(13X,'//CHAR1(1:1)//'D12.4,13X,'//CHAR1(1:1)//'D12.4))')
     '    ((GZL(mi,ni),ni=1,NITB),(GZU(mi,ni),ni=1,NITB),mi=1,NITB)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C old format
C        WRITE(OP_STRING,'('' GZL or AZL by rows: '',9E12.4)')
C     '    ((GZL(mi,ni),ni=1,NITB),mi=1,NITB)
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' GZU or AZU by rows: '',9E12.4)')
C     '    ((GZU(mi,ni),ni=1,NITB),mi=1,NITB)
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZGMG')
      RETURN
 9999 CALL ERRORS('ZGMG',ERROR)
      CALL EXITS('ZGMG')
      RETURN 1
      END


