      SUBROUTINE DBOX(iw,XDIM,YDIM,XWC,YWC,ZWC,ERROR,*)

C#### Subroutine: DBOX
C###  Description:
C###    DBOX draws box of dimensions 2*XDIM,2*YDIM around the point
C###    XWC,YWC,ZWC.

      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER iw
      REAL*8 XDIM,XWC,YDIM,YWC,ZWC
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 Z(3,5)

      CALL ENTERS('DBOX',*9999)
      IF(iw.EQ.1) THEN
        IF(NJT.EQ.2) THEN !plot x against y
          Z(1,1)=XWC-XDIM
          Z(2,1)=YWC-YDIM
          Z(1,2)=XWC+XDIM
          Z(2,2)=Z(2,1)
          Z(1,3)=Z(1,2)
          Z(2,3)=YWC+YDIM
          Z(1,4)=Z(1,1)
          Z(2,4)=Z(2,3)
          Z(1,5)=Z(1,1)
          Z(2,5)=Z(2,1)
          CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
        ELSE IF(NJT.EQ.3) THEN !plot x against z
          Z(1,1)=XWC-XDIM
          Z(3,1)=ZWC-YDIM
          Z(1,2)=XWC+XDIM
          Z(3,2)=Z(3,1)
          Z(1,3)=Z(1,2)
          Z(3,3)=ZWC+YDIM
          Z(1,4)=Z(1,1)
          Z(3,4)=Z(3,3)
          Z(1,5)=Z(1,1)
          Z(3,5)=Z(3,1)
          CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
        ENDIF
      ELSE IF(iw.EQ.2) THEN !plot Y against Z
        Z(2,1)=YWC-XDIM
        Z(3,1)=ZWC-YDIM
        Z(2,2)=YWC+XDIM
        Z(3,2)=Z(3,1)
        Z(2,3)=Z(2,2)
        Z(3,3)=ZWC+YDIM
        Z(2,4)=Z(2,1)
        Z(3,4)=Z(3,3)
        Z(2,5)=Z(2,1)
        Z(3,5)=Z(3,1)
        CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
      ELSE IF(iw.EQ.3) THEN !plot X against Y
        Z(1,1)=XWC-XDIM
        Z(2,1)=YWC-YDIM
        Z(1,2)=XWC+XDIM
        Z(2,2)=Z(2,1)
        Z(1,3)=Z(1,2)
        Z(2,3)=YWC+YDIM
        Z(1,4)=Z(1,1)
        Z(2,4)=Z(2,3)
        Z(1,5)=Z(1,1)
        Z(2,5)=Z(2,1)
        CALL POLYLINE(1,iw,5,Z,ERROR,*9999)
      ELSE IF(iw.EQ.4) THEN
      ENDIF

      CALL EXITS('DBOX')
      RETURN
 9999 CALL ERRORS('DBOX',ERROR)
      CALL EXITS('DBOX')
      RETURN 1
      END


