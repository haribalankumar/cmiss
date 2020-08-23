      SUBROUTINE TEXT(IBUNDLE,iw,STRING,D_PT,ERROR,*)

C#### Subroutine: TEXT
C###  Description:
C###    TEXT draws text on iw with non-geometric attributes bundled in
C###    IBUNDLE, and geometric attributes according to IGEOM.

C**** D_PT(1..3) contains the REAL*8 3D coords of each point.
C**** IF iw=4 coords are curvilinear coords, else rect. cart.
C**** STRING contains the text to be shown at D_PT(nj)
C**** If IBUNDLE is 0 the primitive will use the previously
C**** defined text index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER IBUNDLE,iw
      REAL*8 D_PT(3)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER INDEX,nj,NJ1,NJ2
      REAL HEIGHT,R_PT(2)

      CALL ENTERS('TEXT',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' D_PT(nj): '',3E12.3)')
     '    (D_PT(nj),nj=1,NJT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      INDEX=IBUNDLE

      HEIGHT=0.01*DIAG
      IF(iw.EQ.1) THEN
        NJ1=1
        NJ2=2
      ELSE IF(iw.EQ.2) THEN
        NJ1=2
        NJ2=3
      ELSE IF(iw.EQ.3) THEN
        NJ1=1
        NJ2=3
      ELSE IF(iw.EQ.7) THEN
        NJ1=1
        NJ2=2
      ELSE IF(iw.EQ.9) THEN
        NJ1=1
        NJ2=2
      ELSE IF(iw.EQ.10) THEN
        NJ1=1
        NJ2=2
      ELSE IF(iw.EQ.11) THEN
        NJ1=1
        NJ2=2
      ELSE IF(iw.EQ.12) THEN
        NJ1=1
        NJ2=2
      ELSE IF(iw.EQ.13) THEN
        NJ1=1
        NJ2=2
      ENDIF

      R_PT(1)=REAL(D_PT(NJ1))
      R_PT(2)=REAL(D_PT(NJ2))

      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        CALL TEXT_GX(INDEX,STRING,HEIGHT,R_PT,ERROR,*9999)
      ENDIF

      CALL EXITS('TEXT')
      RETURN
 9999 CALL ERRORS('TEXT',ERROR)
      CALL EXITS('TEXT')
      RETURN 1
      END


