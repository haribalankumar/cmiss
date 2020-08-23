      SUBROUTINE IPQUAS(nr,nx,ERROR,*)

C#### Subroutine: IPQUAS
C###  Description:
C###    IPQUAS inputs quasi static analysis parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'quas00.cmn'
!     Parameter List
      INTEGER nr,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      CHARACTER TIMESTEP*4
      LOGICAL FILEIP

      CALL ENTERS('IPQUAS',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(1,nx)) THEN
        FORMAT='('' Specify the type of quasi-static problem [1]: '''//
     '    '/''   (1) Changing sources                  '''//
     '    '/''   (2) Changing boundary condition values'''//
     '    '/''   (3) Unused                            '''//
     '    '/''   (4) Unused                            '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=QUASIPROB
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) QUASIPROB=IDATA(1)

        IF(QUASIPROB.EQ.2) THEN
          FORMAT='($,'' Enter filename for the input history file '
     '      //'[current]: '',A30)'
          CDEFLT(1)=FILE00
          IF(IOTYPE.EQ.3) CDATA(1)=QUASI_INHISTFILE
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) QUASI_INHISTFILE=CDATA(1)(1:100)
        ENDIF

        FORMAT='('' Specify whether [1]: '''//
     '    '/''   (1) Fixed time step    '''//
     '    '/''   (2) Automatic stepping '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=QUASITIMESTEP
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) QUASITIMESTEP=IDATA(1)

        TIMESTEP='secs'

        FORMAT='($,'' Specify the initial time ('//TIMESTEP//') [0.0]:'
     '    //' '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=T0
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) T0=RDATA(1)

        FORMAT='($,'' Specify the  final time ('//TIMESTEP//') [10.0]:'
     '    //' '',D11.4)'
        RDEFLT(1)=10.0d0
        IF(IOTYPE.EQ.3) RDATA(1)=T1
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) T1=RDATA(1)

        FORMAT='($,'' Specify the time increment (initial if'//
     '    ' automatic stepping) [1.0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=TINCR
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RONE,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TINCR=RDATA(1)

        FORMAT='($,'' Enter filename for the output history file '
     '    //'[current]: '',A30)'
        CDEFLT(1)=FILE00
        IF(IOTYPE.EQ.3) CDATA(1)=QUASI_OUTHISTFILE
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) QUASI_OUTHISTFILE=CDATA(1)(1:100)
        FILE02=QUASI_OUTHISTFILE

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

        FORMAT='($,'' Do you want timing information for each time '
     '    //'step [Y]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(OUTPUT_SOLTIMES) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            OUTPUT_SOLTIMES=.TRUE.
          ELSE IF(ADATA(1).EQ.'N') THEN
            OUTPUT_SOLTIMES=.FALSE.
          ENDIF
        ENDIF

C        FORMAT='($,'' Do you want time-history o/p file for specified'
C     '    //' nodes [N]? '',A)'
C        IF(IOTYPE.EQ.3) THEN
C          IF(NODE_HISTORY_OP) THEN
C            ADATA(1)='Y'
C          ELSE
C            ADATA(1)='N'
C          ENDIF
C        ENDIF
C        CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) THEN
C          IF(ADATA(1).EQ.'Y') THEN
C            NODE_HISTORY_OP=.TRUE.
C          ELSE IF(ADATA(1).EQ.'N') THEN
C            NODE_HISTORY_OP=.FALSE.
C          ENDIF
C        ENDIF
C
C        IF(NODE_HISTORY_OP) THEN
C          FORMAT='($,'' Enter the number of nodes (up to 8) [1]: '','
C     '      I1)'
C          IF(IOTYPE.EQ.3) IDATA(1)=NODE_HISTORY(0)
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,8,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) NODE_HISTORY(0)=IDATA(1)
C
C          CHAR1=CFROMI(NODE_HISTORY(0),'(I1)')
C          FORMAT='($,'' Enter the '//CHAR1//' node numbers [1..'
C     '      //CHAR1//']:'',8I8)'
C          DO N=1,NODE_HISTORY(0)
C            IDEFLT(N)=N
C          ENDDO
C          IF(IOTYPE.EQ.3) THEN
C            DO N=1,NODE_HISTORY(0)
C              IDATA(N)=NODE_HISTORY(N)
C            ENDDO
C          ENDIF
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
C     '      NODE_HISTORY(0),
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
C     '      MAX(NPM,NQM),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
C     '      INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            DO N=1,NODE_HISTORY(0)
C              NODE_HISTORY(N)=IDATA(N)
C            ENDDO
C          ENDIF
C        ENDIF
      ENDIF

      CALL EXITS('IPQUAS')
      RETURN
 9999 CALL ERRORS('IPQUAS',ERROR)
      CALL EXITS('IPQUAS')
      RETURN 1
      END


