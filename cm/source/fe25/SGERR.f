      SUBROUTINE SGERR(INDEX,ISEG,ISERR,iw,MXI,ne,
     '  nr,CSEG,NEERR,Z,ERROR,*)

C#### Subroutine: SGERR
C###  Description:
C###    SGERR creates new element segment ISERR(iw,ne). The
C###    coordinates of the element error label are stored in
C###    rectangular cartesian coordinates in Z(nj). These are
C###    transformed to curvilinear coordinates for iw=4.

C**** Created by Carey Stevens 22 August 1997

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISERR,iw,MXI(2),ne,nr
      REAL*8 NEERR(NEM,3),Z(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD
      CHARACTER CHAR*20

      CALL ENTERS('SGERR',*9999)
      CALL OPEN_SEGMENT(ISERR,ISEG,iw,'ELEM',INDEX,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

C *** Draw element error
      IF(iw.EQ.4) THEN
        PROJEC=MAP_PROJEC
        IF(PROJEC(1:2).EQ.'XI') THEN
          MXI1=MXI(1) !bottom left coords
          MXI2=MXI(2) !..for Xi map projection
          Z(1)=0.5D0    !Xi coords of element
          Z(2)=0.5D0    !..centres
          Z(3)=0.5D0
        ELSE
          CALL ZX(ITYP10(nr),Z,Z) !transform to curvilinear coords
        ENDIF
      ENDIF
      WRITE(CHAR,'(F8.2)') NEERR(ne,1)
      CALL STRING_TRIM(CHAR,IBEG,IEND)
      CHAR=CHAR(IBEG:IEND)//'%'
      CALL STRING_TRIM(CHAR,IBEG,IEND)
      CALL TEXT(INDEX,iw,CHAR(IBEG:IEND),Z,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISERR,iw,ERROR,*9999)
      CALL EXITS('SGERR')
      RETURN
 9999 CALL ERRORS('SGERR',ERROR)
      CALL EXITS('SGERR')
      RETURN 1
      END


