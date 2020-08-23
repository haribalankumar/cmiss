      SUBROUTINE IPSOLT(nr,nx,FILE_NAME,LINE_NAME,ERROR,*)

C#### Subroutine: IPSOLT
C###  Description:
C###    IPSOLT inputs time integration parameters.

C**** T0 & T1 are the initial & final times,respectively.
C**** DT is the time step (initial only if automatic stepping).
C**** LRESID is .true. if half-step residual is calculated.

C#### Variable: FILE02
C###  Type: CHARACTER*(MXCH)
C###  Set_up: IPSOLT, IPQUASI
C###  Description:
C###    The history file name (sans extension) for
C###    time-integration and quasi-static problems.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'time01.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nr,nx
      LOGICAL LINE_NAME
      CHARACTER ERROR*(*),FILE_NAME*(MXCH)
!     Local Variables
      INTEGER I,ICHAR,INFO,NOQUES !SMAR009 22/12/98 ,N
      CHARACTER TIMESTEP*4 !SMAR009 22/12/98 CHAR1*1,
      LOGICAL FILEIP

      CALL ENTERS('IPSOLT',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(ITYP2(nr,nx).EQ.9) THEN
       FORMAT='('' Specify whether time integration algorithm'//
     '    ' is [1]:'''//
     '    '/''   (1) Linear'''//
     '    '/''  *(2) Quadratic'''//
     '    '/''  *(3) Cubic'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP22
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,1,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP22=IDATA(1)

C PM 26-JUL-01
        ELSE IF (ITYP4(nr,nx).NE.3) THEN  ! not finite difference methods
        FORMAT='('' Specify whether time integration algorithm'//
     '    ' is [1]:'''//
     '    '/''   (1) Linear'''//
     '    '/''   (2) Quadratic'''//
     '    '/''   (3) Cubic'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP22
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP22=IDATA(1)
      ELSE
        KTYP22=0    ! for fluid in elastic tube
      ENDIF

C SGM 30 Oct 2000 ignore grid-based Finite element also
C MLT 2Dec2002 ignore grid finite volume also
C IF Not Collocation or grid-based FE problems
C PM 26-JUL-01 : If not finite difference
C MHT 30-10-01 : If not pulmonary
      IF(ITYP4(nr,nx).NE.4.AND.ITYP4(nr,nx).NE.6.AND.
     '   ITYP4(nr,nx).NE.7
     '  .AND.ITYP4(nr,nx).NE.3.AND.ITYP2(nr,nx).NE.11) THEN
        IF(KTYP22.EQ.1) THEN !linear
          CALL ASSERT(NIYM.GE.10,'>>Increase NIYM sb >= 10',
     '      ERROR,*9999)
        ELSE IF(KTYP22.EQ.2) THEN !quadratic
          CALL ASSERT(NIYM.GE.13,'>>Increase NIYM sb >= 13',
     '      ERROR,*9999)
        ELSE IF(KTYP22.EQ.3) THEN !cubic
          CALL ASSERT(NIYM.GE.16,'>>Increase NIYM sb >= 16',
     '      ERROR,*9999)
        ENDIF
      ENDIF

      FORMAT='('' Specify whether [1]: '''//
     '  '/''   (1) Fixed time step    '''//
     '  '/''   (2) Automatic stepping '''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP23
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP23=IDATA(1)

      IF(KTYP22.EQ.1) THEN
        RDEFLT(1)=2.0d0/3.0d0
        FORMAT='($,'' Specify the time integration parameter'//
     '    ' [2/3]: '',G11.4)'
      ELSE IF(KTYP22.EQ.2) THEN
        RDEFLT(1)=5.0d0/6.0d0
        RDEFLT(2)=8.0d0/9.0d0
        FORMAT='($,'' Specify the 2 time integration parameters'//
     '    ' [5/6,8/9]: '',2G11.4)'
      ELSE IF(KTYP22.EQ.3) THEN
        RDEFLT(1)=0.80d0
        RDEFLT(2)=1.03d0
        RDEFLT(3)=1.29d0
        FORMAT='($,'' Specify the 3 time integration parameters'//
     '    ' [0.80,1.03,1.29]: '',3G11.4)'
      ENDIF

      IF(IOTYPE.EQ.3) THEN
        DO I=1,KTYP22
          RDATA(I)=THETA(I)
        ENDDO
      ENDIF

C PM 26-JUL-01
      IF (KTYP22.NE.0) THEN
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    KTYP22,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO I=1,KTYP22
            THETA(I)=RDATA(I)
          ENDDO
        ENDIF
      ENDIF

      IF((KTYP38.EQ.2).AND.(THETA(1).LT.(1.0d0-ZERO_TOL))) THEN
        WRITE(OP_STRING,
     '    '('' >>WARNING: zero int. flux only works for implicit '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(ITYP2(nr,nx).EQ.9) THEN
        TIMESTEP='msec'
        CALL ASSERT(THETA(1).GT.-ZERO_TOL,'>>Theta must be [0,1]',
     '    ERROR,*9999)
        CALL ASSERT(THETA(1).LT.(1.0d0+ZERO_TOL),
     '    '>>Theta must be [0,1]',ERROR,*9999)
      ELSE
        TIMESTEP='secs'
      ENDIF

      FORMAT='($,'' Specify the initial time ('//TIMESTEP//') [0]:'
     '  //' '',D11.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=T0
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) T0=RDATA(1)

      FORMAT='($,'' Specify the  final time ('//TIMESTEP//') [10]:'
     '  //' '',D11.4)'
      RDEFLT(1)=10.0D0
      IF(IOTYPE.EQ.3) RDATA(1)=T1
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) T1=RDATA(1)

      FORMAT='($,'' Specify the time increment (initial if'//
     '  ' automatic stepping) [1.0]: '',D11.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=TINCR
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RONE,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TINCR=RDATA(1)

      FORMAT='($,'' Enter #time step intervals for history file'
     '  //' o/p (0 for no o/p)[1]: '',I3)'
      IF(IOTYPE.EQ.3) IDATA(1)=HIST_file_intervals
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) HIST_file_intervals=IDATA(1)

      IF(HIST_file_intervals.GT.0) THEN
        FORMAT='($,'' Enter filename for the output history file '
     '    //'[current]: '',A30)'
        IF(LINE_NAME)THEN!defined on command line
          CDEFLT(1)=FILE_NAME
          FILE02(1:)=FILE_NAME(1:)
          IF(IOTYPE.EQ.3) CDATA(1)=FILE02
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        ELSE
          CDEFLT(1)=FILE00
          IF(IOTYPE.EQ.3) CDATA(1)=FILE02
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) FILE02(1:100)=CDATA(1)(1:100)
        ENDIF

        FORMAT='($,'' Is the history file to be stored as a binary '
     '    //'file [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(BINTIMEFILE.GT.0) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            BINTIMEFILE=1
          ELSE IF(ADATA(1).EQ.'N') THEN
            BINTIMEFILE=0
          ENDIF
        ENDIF


CMLB 15-Oct-1998 no longer needed as we now have proper history files
C        IF(ITYP2(nr,nx).EQ.9) THEN !Cardiac activation problems only
C          FORMAT='($,'' Do you want time-history o/p file for'
C     '      //' specified nodes [N]? '',A)'
C          IF(IOTYPE.EQ.3) THEN
C            IF(NODE_HISTORY_OP) THEN
C              ADATA(1)='Y'
C            ELSE
C              ADATA(1)='N'
C            ENDIF
C          ENDIF
C          CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            IF(ADATA(1).EQ.'Y') THEN
C              NODE_HISTORY_OP=.TRUE.
C            ELSE IF(ADATA(1).EQ.'N') THEN
C              NODE_HISTORY_OP=.FALSE.
C            ENDIF
C          ENDIF
C
C          IF(NODE_HISTORY_OP) THEN
C            FORMAT='($,'' Enter the number of nodes (up to 8)[1]: '''
C     '        //',I1)'
C            IF(IOTYPE.EQ.3) IDATA(1)=NODE_HISTORY(0)
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,8,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(IOTYPE.NE.3) NODE_HISTORY(0)=IDATA(1)
C
C            WRITE(CHAR1,'(I1)') NODE_HISTORY(0)
C            FORMAT='($,'' Enter the '//CHAR1//' node numbers [1..'
C     '        //CHAR1//']:'',8I8)'
C            DO N=1,NODE_HISTORY(0)
C              IDEFLT(N)=N
C            ENDDO
C            IF(IOTYPE.EQ.3) THEN
C              DO N=1,NODE_HISTORY(0)
C                IDATA(N)=NODE_HISTORY(N)
C              ENDDO
C            ENDIF
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '        NODE_HISTORY(0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '        IDEFLT,1,MAX(NPM,NQM),
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(IOTYPE.NE.3) THEN
C              DO N=1,NODE_HISTORY(0)
C                NODE_HISTORY(N)=IDATA(N)
C              ENDDO
C            ENDIF
C          ENDIF !node history output
C        ENDIF !ITYP2 is cardiac activation
      ENDIF !history file output

      CALL EXITS('IPSOLT')
      RETURN
 9999 CALL ERRORS('IPSOLT',ERROR)
      CALL EXITS('IPSOLT')
      RETURN 1
      END


