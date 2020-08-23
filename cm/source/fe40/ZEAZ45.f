      SUBROUTINE ZEAZ45(NBJ,ng,nr,AXL,AZ,AZL,PG,PPG,
     '  XE,XG,ZE,ZG,ERROR,*)

C#### Subroutine: ZEAZ45
C###  Description:
C###    ZEAZ45 is for membrane elements only.
C###    XE and ZE are used to calculate XG,ZG,GL,AXL,AZL,AZ.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM),ng,nr
      REAL*8 AXL(3,*),AZ,AZL(3,*),PG(NSM,NUM,NGM,NBM),PPG(9,2,*),
     '  XE(NSM,NJM),XG(NJM,NUM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nhx,ni,nj,ns,NU1(0:3)
      REAL*8 SUM1,SUM2
      DATA NU1/1,2,4,7/

      CALL ENTERS('ZEAZ45',*9999)
C *** Calculate undeformed (XG) & deformed (ZG) derivs wrt Nu coords
      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
        nj=NJ_LOC(NJL_GEOM,nhx,nr)
        nb=NBJ(nj)
        SUM1=0.d0
        SUM2=0.d0
        DO ns=1,NST(nb)
          SUM1=SUM1+PG(ns,1,ng,nb)*XE(ns,nj)
          SUM2=SUM2+PG(ns,1,ng,nb)*ZE(ns,nhx)
        ENDDO
        XG(nj,1)=SUM1
        ZG(nhx,1)=SUM2
        DO ni=1,2
          SUM1=0.d0
          SUM2=0.d0
          DO ns=1,NST(nb)
            SUM1=SUM1+PPG(nb,ni,ns)*XE(ns,nj)
            SUM2=SUM2+PPG(nb,ni,ns)*ZE(ns,nhx)
          ENDDO
          XG(nj,NU1(ni))=SUM1
          ZG(nhx,NU1(ni))=SUM2
        ENDDO
      ENDDO

C *** Calculate metric tensors wrt Nu coords
      AXL(1,1)=XG(1,2)*XG(1,2)+XG(2,2)*XG(2,2)+XG(3,2)*XG(3,2)
      AXL(2,2)=XG(1,4)*XG(1,4)+XG(2,4)*XG(2,4)+XG(3,4)*XG(3,4)
      AXL(2,1)=XG(1,4)*XG(1,2)+XG(2,4)*XG(2,2)+XG(3,4)*XG(3,2)
      AXL(1,2)=AXL(2,1)
      AZL(1,1)=ZG(1,2)*ZG(1,2)+ZG(2,2)*ZG(2,2)+ZG(3,2)*ZG(3,2)
      AZL(2,2)=ZG(1,4)*ZG(1,4)+ZG(2,4)*ZG(2,4)+ZG(3,4)*ZG(3,4)
      AZL(2,1)=ZG(1,4)*ZG(1,2)+ZG(2,4)*ZG(2,2)+ZG(3,4)*ZG(3,2)
      AZL(1,2)=AZL(2,1)
      AZ=AZL(1,1)*AZL(2,2)-AZL(2,1)*AZL(1,2)

      CALL EXITS('ZEAZ45')
      RETURN
 9999 CALL ERRORS('ZEAZ45',ERROR)
      CALL EXITS('ZEAZ45')
      RETURN 1
      END


