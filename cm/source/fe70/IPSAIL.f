      SUBROUTINE IPSAIL(NBJ,NEELEM,NKJE,NKJ,NPLIST,NPNE,NPNODE,SE,XP,
     '  ERROR,*)

C#### Subroutine: IPSAIL
C###  Description:
C###    IPSAIl performs sail parameter input.
C**** Luff & leech (chords) are assumed to be in x,z plane.
C**** Luff & leech curves are two element cubic Hermites
C**** with independent x,y and x,z descriptions.
C**** ZTACK(nj),ZCLEW(nj),ZHEADL(nj),ZHEADT(nj) are in global coords
C**** CLUFF is length of luff  chord
C**** CLECH is length of leech chord
C**** ELUFF is distance along luff  to central control pt (set by draft)
C**** ELECH is distance along leech to central control pt (set by draft)
C**** BLUFF is height above luff  chord of central control pt (set by camber)
C**** BLECH is height above leech chord of central control pt (set by camber)
C**** XYLUF1(nk,j,nn) is 1st element of luff  curve in local x,y plane
C**** XYLUF2(nk,j,nn)  " 2nd    "     " luff  curve in local x,y plane
C**** XZLUF1(nk,j,nn) is 1st element of luff  curve in local x,z plane
C**** XZLUF2(nk,j,nn)  " 2nd    "     " luff  curve in local x,z plane
C**** XYLEC1(nk,j,nn) is 1st element of leech curve in local x,y plane
C**** XYLEC2(nk,j,nn)  " 2nd    "     " leech curve in local x,y plane
C**** XZLEC1(nk,j,nn) is 1st element of leech curve in local x,z plane
C**** XZLEC2(nk,j,nn)  " 2nd    "     " leech curve in local x,z plane
C**** XYCHO1(nk,j,nn) is 1st element of section     in local x,y plane
C**** XYCHO2(nk,j,nn)  " 2nd    "     " section     in local x,y plane

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'inout00.cmn'
c     INCLUDE 'SYS$LIBRARY:GKSDEFS.BND'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPLIST(0:NPM),
     '  NKJ(NJM,NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,INFO,IPFILE,j,
     '  nb,NBX,NBY,NBZ,ne,nj,nk,nn,noelem,
     '  nonode,NOQUES,nosect,np,np1,np2,np3,nr,ns,NTEL,NTSECT,nv
      REAL ZCLEW(3),ZHEADL(3),ZHEADT(3),ZTACK(3)
      REAL*8 BCHORD,BLECH,BLUFF,BOOM,CCHORD,CLECH,CLUFF,COS,DH3,DRAFT,
     '  DXLECH,DXLUFF,ECHORD,ELECH,ELUFF,SIN,S1,S2,S3LUFF,
     '  S3LECH,XDIFF1,XDIFF2,XH3,XI,XLECH,XLENG,XLUFF,
     '  XYCHO1(2,2,2),XYCHO2(2,2,2),
     '  XYLEC1(2,2,2),XYLEC2(2,2,2),XYLUF1(2,2,2),XYLUF2(2,2,2),
     '  XZLEC1(2,2,2),XZLEC2(2,2,2),XZLUF1(2,2,2),XZLUF2(2,2,2)
      CHARACTER CHAR*11
      LOGICAL FILEIP

      CALL ENTERS('IPSAIL',*9999)
      nv=1 ! Temporary MPN 12-Nov-94
      IPFILE=1 !is input file version number on 24-Jan-1990
      FILEIP=.FALSE.
      NOQUES=0

      IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.3) THEN
        WRITE(UNIT=IFILE,REC=1,FMT='(A,I2)') 'CMISS Version '//CMISS
     '    //' IPCOOR File Version ',IPFILE
        WRITE(UNIT=IFILE,REC=2,FMT='(A)') 'Heading: '//HEADING
        WRITE(UNIT=IFILE,REC=3,FMT='(1X)')
      ELSE IF(IOTYPE.EQ.2.OR.IOTYPE.EQ.4) THEN
        READ(UNIT=IFILE,REC=1,FMT='(39X,I2)') IPFILE
        IF(DOP) THEN
          WRITE(OP_STRING,'('' File version number is '',I2)') IPFILE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        READ(UNIT=IFILE,REC=2,FMT='(9X,A)') HEADING
        WRITE(OP_STRING,'('' File heading: '',A)') HEADING
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        READ(UNIT=IFILE,REC=3,FMT='(1X)')
      ENDIF

c     CALL GSLN(GLSOLI)
      FORMAT='($,'' Enter x,z-coords of tack [0,0]: '',2E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=ZTACK(1)
        RDATA(2)=ZTACK(3)
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ZTACK(1)=RDATA(1)
        ZTACK(2)=0.D0
        ZTACK(3)=RDATA(2)
      ENDIF

      FORMAT='($,'' Enter x,z-coords of clew [0,0]: '',2E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=ZCLEW(1)
        RDATA(2)=ZCLEW(3)
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ZCLEW(1)=RDATA(1)
        ZCLEW(2)=0.D0
        ZCLEW(3)=RDATA(2)
      ENDIF

      FORMAT='($,'' Enter x,z-coords of head leading  edge [0,0]:'
     '  //' '',2E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=ZHEADL(1)
        RDATA(2)=ZHEADL(3)
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ZHEADL(1)=RDATA(1)
        ZHEADL(2)=0.D0
        ZHEADL(3)=RDATA(2)
      ENDIF
      CLUFF=SQRT((ZHEADL(1)-ZTACK(1))**2+(ZHEADL(3)-ZTACK(3))**2)

      FORMAT='($,'' Enter x,z-coords of head trailing edge [0,0]:'
     '  //' '',2E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=ZHEADT(1)
        RDATA(2)=ZHEADT(3)
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ZHEADT(1)=RDATA(1)
        ZHEADT(2)=0.D0
        ZHEADT(3)=RDATA(2)
      ENDIF
      CLECH=SQRT((ZHEADT(1)-ZCLEW(1))**2+(ZHEADT(3)-ZCLEW(3))**2)

c     CALL ACWK(1,0,ERROR,*9999)
c     CALL GPM(1,ZTACK(1),ZTACK(3))
c     CALL GPM(1,ZCLEW(1),ZCLEW(3))
c     CALL GPM(1,ZHEADL(1),ZHEADL(3))
c     CALL GPM(1,ZHEADT(1),ZHEADT(3))
c     CALL DAWK(1,0,ERROR,*9999)
c     CALL ACWK(2,0,ERROR,*9999)
c     CALL GPM(1,ZTACK(2),ZTACK(3))
c     CALL GPM(1,ZCLEW(2),ZCLEW(3))
c     CALL GPM(1,ZHEADL(2),ZHEADL(3))
c     CALL GPM(1,ZHEADT(2),ZHEADT(3))
c     CALL DAWK(2,0,ERROR,*9999)
C     CALL ACWK(3,0,ERROR,*9999)
C     CALL PHIGS$POLYMARKER3(1,ZTACK)
C     CALL PHIGS$POLYMARKER3(1,ZCLEW)
C     CALL PHIGS$POLYMARKER3(1,ZHEADL)
C     CALL PHIGS$POLYMARKER3(1,ZHEADT)
C     CALL DAWK(3,0,ERROR,*9999)

      DO j=1,2
        XYLUF1(1,j,1)=0.D0
        XZLUF1(1,j,1)=0.D0
        XYLUF2(1,j,2)=0.D0
        XZLUF2(1,j,2)=0.D0
        XYLEC1(1,j,1)=0.D0
        XZLEC1(1,j,1)=0.D0
        XYLEC2(1,j,2)=0.D0
        XZLEC2(1,j,2)=0.D0
      ENDDO
      XYLUF2(1,1,2)=CLUFF
      XZLUF2(1,1,2)=CLUFF
      XYLEC2(1,1,2)=CLECH
      XZLEC2(1,1,2)=CLECH

      FORMAT='($,'' Luff  (fa): Enter draft,camber,bot% and top%'
     '  //' [0,0,0,0]: '',4E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=0.D0
        RDATA(2)=0.D0
        RDATA(3)=0.D0
        RDATA(4)=0.D0
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,4,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ELUFF=RDATA(1)*CLUFF
        BLUFF=RDATA(2)*CLUFF
        XZLUF1(1,1,2)=ELUFF
        XZLUF2(1,1,1)=ELUFF
        XZLUF1(1,2,2)=BLUFF
        XZLUF2(1,2,1)=BLUFF
        XZLUF1(2,2,1)= 8.D0*RDATA(3)*BLUFF-4.D0*BLUFF
        XZLUF1(2,2,2)= 0.D0
        XZLUF2(2,2,1)= 0.D0
        XZLUF2(2,2,2)=-8.D0*RDATA(4)*BLUFF+4.D0*BLUFF
        XDIFF1=XZLUF1(1,1,2)-XZLUF1(1,1,1)
        XZLUF1(2,1,1)=XDIFF1
        XZLUF1(2,1,2)=XDIFF1
        XDIFF2=XZLUF2(1,1,2)-XZLUF2(1,1,1)
        XZLUF2(2,1,1)=XDIFF2
        XZLUF2(2,1,2)=XDIFF2
      ENDIF

      FORMAT='($,'' Luff  (sw): Enter draft,camber,bot% and top%'
     '  //' [0,0,0,0]: '',4E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=0.D0
        RDATA(2)=0.D0
        RDATA(3)=0.D0
        RDATA(4)=0.D0
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,4,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ELUFF=RDATA(1)*CLUFF
        BLUFF=RDATA(2)*CLUFF
        XYLUF1(1,1,2)=ELUFF
        XYLUF2(1,1,1)=ELUFF
        XYLUF1(1,2,2)=BLUFF
        XYLUF2(1,2,1)=BLUFF
        XYLUF1(2,2,1)= 8.D0*RDATA(3)*BLUFF-4.D0*BLUFF
        XYLUF1(2,2,2)= 0.D0
        XYLUF2(2,2,1)= 0.D0
        XYLUF2(2,2,2)=-8.D0*RDATA(4)*BLUFF+4.D0*BLUFF
        XDIFF1=XYLUF1(1,1,2)-XYLUF1(1,1,1)
        XYLUF1(2,1,1)=XDIFF1
        XYLUF1(2,1,2)=XDIFF1
        XDIFF2=XYLUF2(1,1,2)-XYLUF2(1,1,1)
        XYLUF2(2,1,1)=XDIFF2
        XYLUF2(2,1,2)=XDIFF2
      ENDIF

      COS=(ZHEADL(1)-ZTACK(1))/CLUFF
      SIN=(ZHEADL(3)-ZTACK(3))/CLUFF
c Note if reuse this code need REAL*8 XX(40,2)
c     CALL ACWK(1,0,ERROR,*9999)
c     CALL XHERM(COS,-SIN,SIN,COS,XZLUF1,XZLUF2,ZTACK(1),ZTACK(3),XX,
c    '  ERROR,*9999)
c     CALL GPL(40,XX(1,1),XX(1,2))
c     CALL DAWK(1,0,ERROR,*9999)
c     CALL ACWK(2,0,ERROR,*9999)
c     CALL XHERM(0.,1.,SIN,0.,XYLUF1,XYLUF2,ZTACK(2),ZTACK(3),XX,ERROR,
c    '  *9999)
c     CALL GPL(40,XX(1,1),XX(1,2))
c     CALL DAWK(2,0,ERROR,*9999)
C     CALL ACWK(3,0,ERROR,*9999)
C     CALL XHERM(COS,0.,0.,1.,XYLUF1,XYLUF2,ZTACK(1),ZTACK(2),XX,ERROR,
C    '  *9999)
C     CALL GPL(40,XX(1,1),XX(1,2))
C     CALL DAWK(3,0,ERROR,*9999)

      FORMAT='($,'' Leech (fa): Enter draft,camber,bot% and top%'
     '  //' [0,0,0,0]: '',4E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=0.D0
        RDATA(2)=0.D0
        RDATA(3)=0.D0
        RDATA(4)=0.D0
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,4,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ELECH=RDATA(1)*CLECH
        BLECH=RDATA(2)*CLECH
        XZLEC1(1,1,2)= ELECH
        XZLEC2(1,1,1)= ELECH
        XZLEC1(1,2,2)=-BLECH
        XZLEC2(1,2,1)=-BLECH
        XZLEC1(2,2,1)=-( 8.D0*RDATA(3)*BLECH-4.D0*BLECH)
        XZLEC1(2,2,2)= 0.D0
        XZLEC2(2,2,1)= 0.D0
        XZLEC2(2,2,2)=-(-8.D0*RDATA(4)*BLECH+4.D0*BLECH)
        XDIFF1=XZLEC1(1,1,2)-XZLEC1(1,1,1)
        XZLEC1(2,1,1)=XDIFF1
        XZLEC1(2,1,2)=XDIFF1
        XDIFF2=XZLEC2(1,1,2)-XZLEC2(1,1,1)
        XZLEC2(2,1,1)=XDIFF2
        XZLEC2(2,1,2)=XDIFF2
      ENDIF

      FORMAT='($,'' Leech (sw): Enter draft,camber,bot% and top%'
     '  //' [0,0,0,0]: '',4E11.4)'
      IF(IOTYPE.EQ.3) THEN
        RDATA(1)=0.D0
        RDATA(2)=0.D0
        RDATA(3)=0.D0
        RDATA(4)=0.D0
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,4,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) THEN
        ELECH=RDATA(1)*CLECH
        BLECH=RDATA(2)*CLECH
        XYLEC1(1,1,2)=ELECH
        XYLEC2(1,1,1)=ELECH
        XYLEC1(1,2,2)=BLECH
        XYLEC2(1,2,1)=BLECH
        XYLEC1(2,2,1)= 8.D0*RDATA(3)*BLECH-4.D0*BLECH
        XYLEC1(2,2,2)= 0.D0
        XYLEC2(2,2,1)= 0.D0
        XYLEC2(2,2,2)=-8.D0*RDATA(4)*BLECH+4.D0*BLECH
        XDIFF1=XYLEC1(1,1,2)-XYLEC1(1,1,1)
        XYLEC1(2,1,1)=XDIFF1
        XYLEC1(2,1,2)=XDIFF1
        XDIFF2=XYLEC2(1,1,2)-XYLEC2(1,1,1)
        XYLEC2(2,1,1)=XDIFF2
        XYLEC2(2,1,2)=XDIFF2
      ENDIF

      COS=(ZHEADT(1)-ZCLEW(1))/CLECH
      SIN=(ZHEADT(3)-ZCLEW(3))/CLECH
c     CALL ACWK(1,0,ERROR,*9999)
c     CALL XHERM(COS,-SIN,SIN,COS,XZLEC1,XZLEC2,ZCLEW(1),ZCLEW(3),XX,
c    '  ERROR,*9999)
c     CALL GPL(40,XX(1,1),XX(1,2))
c     CALL DAWK(1,0,ERROR,*9999)
c     CALL ACWK(2,0,ERROR,*9999)
c     CALL XHERM(0.,1.,SIN,0.,XYLEC1,XYLEC2,ZCLEW(2),ZCLEW(3),XX,ERROR,
c    '  *9999)
c     CALL GPL(40,XX(1,1),XX(1,2))
c     CALL DAWK(2,0,ERROR,*9999)
C     CALL ACWK(3,0,ERROR,*9999)
C     CALL XHERM(COS,0.,0.,1.,XYLEC1,XYLEC2,ZCLEW(1),ZCLEW(2),XX,
C    '  ERROR,*9999)
C     CALL GPL(40,XX(1,1),XX(1,2))
C     CALL DAWK(3,0,ERROR,*9999)

      WRITE(CHAR,'(I3)') NPT(1)+1
      CALL STRING_TRIM(CHAR,IBEG,IEND)
      IDEFLT(1)=NPT(1)+1
      FORMAT='($,'' Enter node number at tack ['//CHAR(IBEG:IEND)//
     '  ']: '',I3)'
      IF(IOTYPE.EQ.3) IDATA(1)=np
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) np=IDATA(1)

      WRITE(CHAR,'(I3)') NET(1)+1
      CALL STRING_TRIM(CHAR,IBEG,IEND)
      IDEFLT(1)=NET(1)+1
      FORMAT='($,'' Enter element number adjacent to tack ['//
     '  CHAR(IBEG:IEND)//']: '',I3)'
      IF(IOTYPE.EQ.3) IDATA(1)=ne
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NEM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) ne=IDATA(1)

      FORMAT='($,'' Enter basis number for x-coordinate [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NBX
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NBFT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) NBX=IDATA(1)
      FORMAT='($,'' Enter basis number for y-coordinate [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NBY
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBFT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) NBY=IDATA(1)
      FORMAT='($,'' Enter basis number for z-coordinate [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NBZ
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NBFT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) NBZ=IDATA(1)

      FORMAT='($,'' Enter number of sections in sail [1]: '',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=NTSECT
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.NE.3) NTSECT=IDATA(1)

      DXLUFF=CLUFF/DBLE(NTSECT-1)
      DXLECH=CLECH/DBLE(NTSECT-1)

      DO 800 nosect=1,NTSECT
        np1=np
C GMH 8/1/97 Update cmgui link
        CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
        np2=np1+1
        np3=np1+2
        XLUFF=DXLUFF*DBLE(nosect-1)
        XLECH=DXLECH*DBLE(nosect-1)
        COS=(ZHEADL(1)-ZTACK(1))/CLUFF
        SIN=(ZHEADL(3)-ZTACK(3))/CLUFF
        IF(XLUFF.LE.XZLUF1(1,1,2)) THEN
          XI=XLUFF/XZLUF1(1,1,2)
          XP(1,nv,1,np1)=ZTACK(1)+XH3(1,XZLUF1,XI)*COS-
     '      XH3(2,XZLUF1,XI)*SIN
          XP(1,nv,2,np1)=ZTACK(2)+XH3(2,XYLUF1,XI)
          XP(1,nv,3,np1)=ZTACK(3)+XH3(1,XZLUF1,XI)*SIN+
     '      XH3(2,XZLUF1,XI)*COS
          S3LUFF=XLENG(XYLUF1(1,1,1),XYLUF1(1,1,2))
          XP(3,nv,2,np1)=DH3(2,XYLUF1,XI)/S3LUFF
        ELSE IF(XLUFF.GT.XZLUF1(1,1,2)) THEN
          XI=(XLUFF-XZLUF2(1,1,1))/(XZLUF2(1,1,2)-XZLUF2(1,1,1))
          XP(1,nv,1,np1)=ZTACK(1)+XH3(1,XZLUF2,XI)*COS-
     '      XH3(2,XZLUF2,XI)*SIN
          XP(1,nv,2,np1)=ZTACK(2)+XH3(2,XYLUF2,XI)
          XP(1,nv,3,np1)=ZTACK(3)+XH3(1,XZLUF2,XI)*SIN+
     '      XH3(2,XZLUF2,XI)*COS
          S3LUFF=XLENG(XYLUF2(1,1,1),XYLUF2(1,1,2))
          XP(3,nv,2,np1)=DH3(2,XYLUF2,XI)/S3LUFF
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Luff : s3luff='',E11.3,'' xp(3,2,np1)='',E11.3)')
     '      S3LUFF,XP(3,nv,2,np1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
c       CALL ACWK(1,0,ERROR,*9999)
c       CALL GPM(1,REAL(XP(1,nv,1,np1)),REAL(XP(1,nv,3,np1)))
c       CALL DAWK(1,0,ERROR,*9999)
c       CALL ACWK(2,0,ERROR,*9999)
c       CALL GPM(1,REAL(XP(1,nv,2,np1)),REAL(XP(1,nv,3,np1)))
c       CALL DAWK(2,0,ERROR,*9999)

        COS=(ZHEADT(1)-ZCLEW(1))/CLECH
        SIN=(ZHEADT(3)-ZCLEW(3))/CLECH
        IF(XLECH.LE.XZLEC1(1,1,2)) THEN
          XI=XLECH/XZLEC1(1,1,2)
          XP(1,nv,1,np3)=ZCLEW(1)+XH3(1,XZLEC1,XI)*COS-
     '      XH3(2,XZLEC1,XI)*SIN
          XP(1,nv,2,np3)=ZCLEW(2)+XH3(2,XYLEC1,XI)
          XP(1,nv,3,np3)=ZCLEW(3)+XH3(1,XZLEC1,XI)*SIN+
     '      XH3(2,XZLEC1,XI)*COS
          S3LECH=XLENG(XYLEC1(1,1,1),XYLEC1(1,1,2))
          XP(3,nv,2,np3)=DH3(2,XYLEC1,XI)/S3LECH
        ELSE IF(XLECH.GT.XZLEC1(1,1,2)) THEN
          XI=(XLECH-XZLEC2(1,1,1))/(XZLEC2(1,1,2)-XZLEC2(1,1,1))
          XP(1,nv,1,np3)=ZCLEW(1)+XH3(1,XZLEC2,XI)*COS-
     '      XH3(2,XZLEC2,XI)*SIN
          XP(1,nv,2,np3)=ZCLEW(2)+XH3(2,XYLEC2,XI)
          XP(1,nv,3,np3)=ZCLEW(3)+XH3(1,XZLEC2,XI)*SIN+
     '      XH3(2,XZLEC2,XI)*COS
          S3LECH=XLENG(XYLEC2(1,1,1),XYLEC2(1,1,2))
          XP(3,nv,2,np3)=DH3(2,XYLEC2,XI)/S3LECH
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Leech: s3lech='',E11.3,'' xp(3,2,np3)='',E11.3)')
     '      S3LECH,XP(3,nv,2,np3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        XP(3,nv,2,np2)=0.5D0*(XP(3,nv,2,np1)+XP(3,nv,2,np3))
c       CALL ACWK(1,0,ERROR,*9999)
c       CALL GPM(1,REAL(XP(1,nv,1,np3)),REAL(XP(1,nv,3,np3)))
c       CALL DAWK(1,0,ERROR,*9999)
c       CALL ACWK(2,0,ERROR,*9999)
c       CALL GPM(1,REAL(XP(1,nv,2,np3)),REAL(XP(1,nv,3,np3)))
c       CALL DAWK(2,0,ERROR,*9999)

        CCHORD=SQRT((XP(1,nv,1,np3)-XP(1,nv,1,np1))**2+
     '    (XP(1,nv,2,np3)-XP(1,nv,2,np1))**2+
     '    (XP(1,nv,3,np3)-XP(1,nv,3,np1))**2)
        IF(nosect.EQ.1) BOOM=CCHORD
        XYCHO1(1,1,1)=0.D0
        XYCHO1(1,2,1)=0.D0
        XYCHO2(1,1,2)=CCHORD
        XYCHO2(1,2,2)=0.D0

        FORMAT='($,'' Enter draft,camber,front% and back%'
     '    //' [0,0,0,0]: '',4E11.4)'
        IF(IOTYPE.EQ.3) THEN
          RDATA(1)=0.D0
          RDATA(2)=0.D0
          RDATA(3)=0.D0
          RDATA(4)=0.D0
        ENDIF
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,4,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(iotype.NE.3) THEN
          DRAFT =RDATA(1)
          ECHORD=DRAFT*CCHORD
          BCHORD=RDATA(2)*CCHORD
          XYCHO1(1,1,2)=ECHORD
          XYCHO2(1,1,1)=ECHORD
          XYCHO1(1,2,2)=BCHORD
          XYCHO2(1,2,1)=BCHORD
          XYCHO1(2,2,1)= 8.D0*RDATA(3)*BCHORD-4.D0*BCHORD
          XYCHO1(2,2,2)= 0.D0
          XYCHO2(2,2,1)= 0.D0
          XYCHO2(2,2,2)=-8.D0*RDATA(4)*BCHORD+4.D0*BCHORD
          XDIFF1=XYCHO1(1,1,2)-XYCHO1(1,1,1)
          XYCHO1(2,1,1)=XDIFF1
          XYCHO1(2,1,2)=XDIFF1
          XDIFF2=XYCHO2(1,1,2)-XYCHO2(1,1,1)
          XYCHO2(2,1,1)=XDIFF2
          XYCHO2(2,1,2)=XDIFF2
        ENDIF

c       CALL ACWK(1,0,ERROR,*9999)
c       XX(1,1)=XP(1,nv,1,np1)
c       XX(1,2)=XP(1,nv,3,np1)
c       XX(2,1)=XP(1,nv,1,np3)
c       XX(2,2)=XP(1,nv,3,np3)
c       CALL GPL(2,XX(1,1),XX(1,2))
c       CALL DAWK(1,0,ERROR,*9999)
c       CALL SCHORD(XYCHO1,XYCHO2,1./CCHORD,ERROR,*9999)

        COS=(XP(1,nv,1,np3)-XP(1,nv,1,np1))/CCHORD
        SIN=(XP(1,nv,2,np3)-XP(1,nv,2,np1))/CCHORD
C       CALL ACWK(3,0,ERROR,*9999)
C       CALL XHERM(COS,-SIN,SIN,COS,XYCHO1,XYCHO2,
C    '    XP(1,nv,1,np1),XP(1,nv,2,np1),XX,ERROR,*9999)
C       CALL GPL(40,XX(1,1),XX(1,2))
C       CALL DAWK(3,0,ERROR,*9999)

        S1=XLENG(XYCHO1(1,1,1),XYCHO1(1,1,2))
        S2=XLENG(XYCHO2(1,1,1),XYCHO2(1,1,2))
        XP(1,nv,1,np2)=XP(1,nv,1,np1)+XYCHO1(1,1,2)*COS-
     '    XYCHO1(1,2,2)*SIN
        XP(1,nv,2,np2)=XP(1,nv,2,np1)+XYCHO1(1,1,2)*SIN+
     '    XYCHO1(1,2,2)*COS
        XP(1,nv,3,np2)=(1.-DRAFT)*XP(1,nv,3,np1)+DRAFT*XP(1,nv,3,np3)
        XP(2,nv,1,np1)=(XYCHO1(2,1,1)*COS-XYCHO1(2,2,1)*SIN)/S1
        XP(2,nv,2,np1)=(XYCHO1(2,1,1)*SIN+XYCHO1(2,2,1)*COS)/S1
        XP(2,nv,3,np1)=(XP(1,nv,3,np3)-XP(1,nv,3,np1))/CCHORD
        XP(2,nv,1,np2)=(XYCHO1(2,1,2)*COS-XYCHO1(2,2,2)*SIN)/S1
        XP(2,nv,2,np2)=(XYCHO1(2,1,2)*SIN+XYCHO1(2,2,2)*COS)/S1
        XP(2,nv,3,np2)=(XP(1,nv,3,np3)-XP(1,nv,3,np1))/CCHORD
        XP(2,nv,1,np3)=(XYCHO2(2,1,2)*COS-XYCHO2(2,2,2)*SIN)/S2
        XP(2,nv,2,np3)=(XYCHO2(2,1,2)*SIN+XYCHO2(2,2,2)*COS)/S2
        XP(2,nv,3,np3)=(XP(1,nv,3,np3)-XP(1,nv,3,np1))/CCHORD

        IF(nosect.GT.1) THEN

          NPNE(3,NBX,ne-2)=np1
          NPNE(4,NBX,ne-2)=np2
          NPNE(3,NBX,ne-1)=np2
          NPNE(4,NBX,ne-1)=np3
          NPNE(3,NBY,ne-2)=np1
          NPNE(4,NBY,ne-2)=np2
          NPNE(3,NBY,ne-1)=np2
          NPNE(4,NBY,ne-1)=np3
          NPNE(3,NBZ,ne-2)=np1
          NPNE(4,NBZ,ne-2)=np2
          NPNE(3,NBZ,ne-1)=np2
          NPNE(4,NBZ,ne-1)=np3


          DO nn=3,4
            DO nk=1,NKT(nn,NBX)
              ns=nk+(nn-1)*NKT(0,NBX) !Needs updating AJP 25-5-93
              IF(nk.EQ.1) THEN
                SE(ns,NBX,ne-2)=1.D0
                SE(ns,NBX,ne-1)=1.D0
              ELSE IF(nk.EQ.2) THEN
                SE(ns,NBX,ne-2)=S1
                SE(ns,NBX,ne-1)=S2
              ENDIF
            ENDDO
            DO nk=1,NKT(nn,NBY)
              ns=nk+(nn-1)*NKT(0,NBY)
              IF(nk.EQ.1) THEN
                SE(ns,NBY,ne-2)=1.D0
                SE(ns,NBY,ne-1)=1.D0
              ELSE IF(nk.EQ.2) THEN
                SE(ns,NBY,ne-2)=S1
                SE(ns,NBY,ne-1)=S2
              ENDIF
            ENDDO
            DO nk=1,NKT(nn,NBZ)
              ns=nk+(nn-1)*NKT(0,NBZ)
              IF(nk.EQ.1) THEN
                SE(ns,NBZ,ne-2)=1.D0
                SE(ns,NBZ,ne-1)=1.D0
              ELSE IF(nk.EQ.2) THEN
                SE(ns,NBZ,ne-2)=S1
                SE(ns,NBZ,ne-1)=S2
              ENDIF
            ENDDO
          ENDDO

          NTEL=(NTSECT-1)/2
          SE(3+(1-1)*NKT(0,NBY),NBY,ne-2)=S3LUFF/DBLE(NTEL)
          SE(3+(3-1)*NKT(0,NBY),NBY,ne-2)=S3LUFF/DBLE(NTEL)
          SE(3+(2-1)*NKT(0,NBY),NBY,ne-2)=0.5D0*(S3LUFF+S3LECH)
     '                                         /DBLE(NTEL)
          SE(3+(4-1)*NKT(0,NBY),NBY,ne-2)=0.5D0*(S3LUFF+S3LECH)
     '                                         /DBLE(NTEL)
          SE(3+(1-1)*NKT(0,NBY),NBY,ne-1)=0.5D0*(S3LUFF+S3LECH)
     '                                         /DBLE(NTEL)
          SE(3+(3-1)*NKT(0,NBY),NBY,ne-1)=0.5D0*(S3LUFF+S3LECH)
     '                                         /DBLE(NTEL)
          SE(3+(2-1)*NKT(0,NBY),NBY,ne-1)=S3LECH/DBLE(NTEL)
          SE(3+(4-1)*NKT(0,NBY),NBY,ne-1)=S3LECH/DBLE(NTEL)

        ENDIF

        IF(nosect.LT.NTSECT) THEN

          NPNE(1,NBX,ne  )=np1
          NPNE(2,NBX,ne  )=np2
          NPNE(1,NBX,ne+1)=np2
          NPNE(2,NBX,ne+1)=np3
          NPNE(1,NBY,ne  )=np1
          NPNE(2,NBY,ne  )=np2
          NPNE(1,NBY,ne+1)=np2
          NPNE(2,NBY,ne+1)=np3
          NPNE(1,NBZ,ne  )=np1
          NPNE(2,NBZ,ne  )=np2
          NPNE(1,NBZ,ne+1)=np2
          NPNE(2,NBZ,ne+1)=np3


          NBJ(1,ne  )=NBX
          NBJ(1,ne+1)=NBX
          NBJ(2,ne  )=NBY
          NBJ(2,ne+1)=NBY
          NBJ(3,ne  )=NBZ
          NBJ(3,ne+1)=NBZ

          DO nj=1,3
            nb=NBJ(nj)
            DO nn=1,NKT(0,nb)
              DO nk=1,NKT(0,nb)
                NKJE(nk,nn,nj,ne  )=nk
                NKJE(nk,nn,nj,ne+1)=nk
              ENDDO
            ENDDO
          ENDDO

          DO nn=1,2
            DO nk=1,NKT(nn,NBX)
              ns=nk+(nn-1)*NKT(0,NBX)
              IF(nk.EQ.1) THEN
                SE(ns,NBX,ne  )=1.D0
                SE(ns,NBX,ne+1)=1.D0
              ELSE IF(nk.EQ.2) THEN
                SE(ns,NBX,ne  )=S1
                SE(ns,NBX,ne+1)=S2
              ENDIF
            ENDDO
            DO nk=1,NKT(nn,NBY)
              ns=nk+(nn-1)*NKT(0,NBY)
              IF(nk.EQ.1) THEN
                SE(ns,NBY,ne  )=1.D0
                SE(ns,NBY,ne+1)=1.D0
              ELSE IF(nk.EQ.2) THEN
                SE(ns,NBY,ne  )=S1
                SE(ns,NBY,ne+1)=S2
              ENDIF
            ENDDO
            DO nk=1,NKT(nn,NBZ)
              ns=nk+(nn-1)*NKT(0,NBZ)
              IF(nk.EQ.1) THEN
                SE(ns,NBZ,ne  )=1.D0
                SE(ns,NBZ,ne+1)=1.D0
              ELSE IF(nk.EQ.2) THEN
                SE(ns,NBZ,ne  )=S1
                SE(ns,NBZ,ne+1)=S2
              ENDIF
            ENDDO
          ENDDO

          ne=ne+2
          np=np+3
        ENDIF

 800  CONTINUE

      nr=1 !may need generalizing
      NPT(nr)=np+2
      NPNODE(0,nr)=NPT(nr)
      NPNODE(0,0) =NPT(nr)
      DO nonode=1,NPNODE(0,nr)
        NPNODE(nonode,nr)=nonode
      ENDDO
      NET(nr)=ne-1
      NEELEM(0,nr)=NET(nr)
      NEELEM(0,0) =NET(nr)
      DO noelem=1,NEELEM(0,nr)
        NEELEM(noelem,nr)=noelem
      ENDDO
      CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST,NPNE,NPNODE,ERROR,*9999)

c     CALL ACWK(18,0,ERROR,*9999)
c     CALL GCLRWK(18,1)
c     CALL DAWK(18,0,ERROR,*9999)

      CALL EXITS('IPSAIL')
      RETURN
 9999 CALL ERRORS('IPSAIL',ERROR)
      CALL EXITS('IPSAIL')
      RETURN 1
      END


