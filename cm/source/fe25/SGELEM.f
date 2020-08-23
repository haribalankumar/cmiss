      SUBROUTINE SGELEM(INDEX,ISEG,ISELNO,iw,MXI,NBJ,ne,
     '  NLL,NPL,nr,BOUNDARY,BOX,CSEG,DL,XP,Z,ERROR,*)

C#### Subroutine: SGELEM
C###  Description:
C###    SGELEM creates new element segment ISELNO(iw,ne). The
C###    coordinates of the element label are stored in rectangular
C###    cartesian coordinates in Z(nj). These are transformed to
C###    curvilinear coordinates for iw=4.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISELNO,iw,MXI(2),NBJ(NJM),ne,NLL(12),
     '  NPL(5,0:3,NLM),nr
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM),Z(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL BOUNDARY,BOX
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD,nae,nb,nl,NTDX
      REAL*8 SCALE,XL(3,20)
      CHARACTER CHAR4*4

      CALL ENTERS('SGELEM',*9999)
      CALL OPEN_SEGMENT(ISELNO,ISEG,iw,'ELEM',INDEX,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

C *** Draw element number
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
      WRITE(CHAR4,'(I4)') ne
      CALL STRING_TRIM(CHAR4,IBEG,IEND)
      SCALE=DBLE(IEND-IBEG+1)
      CALL TEXT(1,iw,CHAR4(IBEG:IEND),Z,ERROR,*9999)
      IF(BOX) THEN !Draw box around element
        IF(iw.NE.4) THEN
          CALL DBOX(iw,0.015D0*SCALE*DBLE(DIAG),0.015D0*DBLE(DIAG),
     '      Z(1),Z(2),Z(3),ERROR,*9999)
        ELSE IF(iw.EQ.4) THEN
          CALL DBOX(iw,0.03D0*SCALE,0.03D0,Z(1),Z(2),Z(3),ERROR,*9999)
        ENDIF
      ENDIF

C *** Draw element boundaries
      IF(BOUNDARY) THEN
C cpb 12/3/97 NTDX was not initialised before use so just set it to
C 20 for now. In the future this should be set according to the line
C type.
        NTDX=20
        nb=NBJ(1)
        DO nae=1,NLE(nb)
          nl=NLL(nae)
          IF(nl.NE.0) THEN
            IF(iw.EQ.1) THEN
              CALL XPXL('UNDEFORMED',1,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,
     '          XP,ERROR,*9999)
            ELSE IF(iw.EQ.4) THEN
              CALL XPXL('UNDEFORMED',0,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,
     '          XP,ERROR,*9999)
            ENDIF
          ENDIF
          CALL POLYLINE(INDEX,iw,NTDX,XL,ERROR,*9999)
        ENDDO
      ENDIF

      CALL CLOSE_SEGMENT(ISELNO,iw,ERROR,*9999)
      CALL EXITS('SGELEM')
      RETURN
 9999 CALL ERRORS('SGELEM',ERROR)
      CALL EXITS('SGELEM')
      RETURN 1
      END


