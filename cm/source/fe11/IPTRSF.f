      SUBROUTINE IPTRSF(NELIST3,NHP,NKH,NPLIST3,NPLIST4,NPLIST5,NPNY,
     '  NVHP,NYNP,NYNR,YP,FIX,ERROR,*)

C#### Subroutine: IPTRSF
C###  Description:
C###    IPTRSF inputs transfer matrix (from epicardium to body surface)
C###    parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER NELIST3(0:NEM),NHP(NPM,0:NRM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPLIST3(0:NPM),NPLIST4(0:NPM),NPLIST5(0:NPM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,INFO,nc,nh,nhx,nk,n1list,nolist,
     '  no_nynr,NOQUES,np,nr,nv,nx,ny
      CHARACTER CHAR1*7,CHAR2*1
      LOGICAL FILEIP,INLIST

      CALL ENTERS('IPTRSF',*9999)
      nx=1 !temporary
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='(/$,'' Enter the number of first surface nodes '
     '  //'[1]: '',I7)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=NPLIST3(0)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NPLIST3(0)=IDATA(1)
      CALL ASSERT(NPLIST3(0).LE.NPM,'>>Increase NPM',ERROR,*9999)
      WRITE(CHAR1,'(I7)') NPLIST3(0)
      CALL STRING_TRIM(CHAR1,IBEG,IEND)
      FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//' first surface '
     '  //'nodes: '',2I7,/:(3X,9I7))'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        DO nolist=1,NPLIST3(0)
          IDATA(nolist)=NPLIST3(nolist)
        ENDDO
      ENDIF
      CDATA(1)='NODES' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  NPLIST3(0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        DO nolist=1,NPLIST3(0)
          NPLIST3(nolist)=IDATA(nolist)
        ENDDO
      ENDIF

      FORMAT='($,'' Enter the region number of these nodes [1]: '',I5)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=TRSF_NR_FIRST
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NRT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TRSF_NR_FIRST=IDATA(1)

      FORMAT='($,'' Enter the number of second surface nodes '
     '  //'[1]: '',I7)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=NPLIST4(0)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NPLIST4(0)=IDATA(1)
      CALL ASSERT(NPLIST4(0).LE.NPM,'>>Increase NPM',ERROR,*9999)
      WRITE(CHAR1,'(I7)') NPLIST4(0)
      CALL STRING_TRIM(CHAR1,IBEG,IEND)
      FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//' second surface '
     '  //'nodes: '',2I7,/:(3X,9I7))'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        DO nolist=1,NPLIST4(0)
          IDATA(nolist)=NPLIST4(nolist)
        ENDDO
      ENDIF
      CDATA(1)='NODES' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  NPLIST4(0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        DO nolist=1,NPLIST4(0)
          NPLIST4(nolist)=IDATA(nolist)
        ENDDO
      ENDIF

      FORMAT='($,'' Enter the region number of these nodes [1]: '',I5)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=TRSF_NR_SECOND
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NRT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TRSF_NR_SECOND=IDATA(1)

      FORMAT='($,'' Enter the total # of regions involved in the '
     '  //'transfer matrix [1]: '',I5)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=TRSF_NRLIST(0)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TRSF_NRLIST(0)=IDATA(1)
      WRITE(CHAR2,'(I1)') TRSF_NRLIST(0)
      FORMAT='($,'' Enter the '//CHAR2(1:1)//' region numbers'
     '  //': '',9(1X,I1))'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        DO nolist=1,TRSF_NRLIST(0)
          IDATA(nolist)=TRSF_NRLIST(nolist)
        ENDDO
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  TRSF_NRLIST(0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '  1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        DO nolist=1,TRSF_NRLIST(0)
          TRSF_NRLIST(nolist)=IDATA(nolist)
        ENDDO
      ENDIF

      FORMAT='($,'' Enter the sampling frequency (kHz) [1.0]: '',D12.4)'
      RDEFLT(1)=1.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=TRSF_FREQUENCY
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,ZERO_TOL,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TRSF_FREQUENCY=RDATA(1)

      FORMAT='('' Specify how transfer matrix is to be calc [1]: '''
     '  //'/''   (1) Directly'''
     '  //'/''   (2) Algebraicly'''
     '  //'/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP94
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP94=IDATA(1)

      FORMAT='('' Select type of transfer matrix [1]: '''
     '  //'/''   (1) Single layer (epi. potls to body surface)'''
     '  //'/''   (2) Double layer (transmembrane to body surface)'''
     '  //'/''   (3) Double layer (transmembrane to epi surface)'''
     '  //'/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP95
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP95=IDATA(1)

      IF(KTYP95.GT.1) THEN
        IF(KTYP94.EQ.2) THEN
          ERROR='>>Double layer matrix can only be calc directly'
          GOTO 9999
        ENDIF
      ENDIF

      IF(KTYP95.EQ.2.OR.KTYP95.EQ.3) THEN

        FORMAT='(/'' Select type of wavefront function [2]: '''
     '    //'/''   (1) Heaviside step function'''
     '    //'/''   (2) Sigmoidal function'''
     '    //'/''   (3) Arctan function'''
     '    //'/$,''    '',I1)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=TRSF_ACTN_WAVE_TYPE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_TYPE=IDATA(1)

        IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN

          FORMAT='($,'' Enter the resting potential (mV) [-80.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=-80.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_REST
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_REST=RDATA(1)

          FORMAT='($,'' Enter the transmembrane jump (mV) [100.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=100.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_JUMP
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_JUMP=RDATA(1)

        ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN

          FORMAT='($,'' Enter the resting potential (mV) [-80.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=-80.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_REST
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_REST=RDATA(1)

          FORMAT='($,'' Enter the transmembrane jump (mV) [100.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=100.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_JUMP
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_JUMP=RDATA(1)

          FORMAT='($,'' Enter the wavefront width (ms) [3.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=3.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_WIDTH
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_WIDTH=RDATA(1)

        ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN
          FORMAT='($,'' Enter the resting potential (mV) [-80.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=-80.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_REST
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_REST=RDATA(1)

          FORMAT='($,'' Enter the transmembrane jump (mV) [100.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=100.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_JUMP
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_JUMP=RDATA(1)

          FORMAT='($,'' Enter the wavefront width (ms) [3.0]: '','
     '      //'D12.4)'
          RDEFLT(1)=3.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=TRSF_ACTN_WAVE_WIDTH
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_WIDTH=RDATA(1)

        ENDIF

        FORMAT='(/'' Is the activation field [1]: '''
     '    //'/''   (1) Discrete'''
     '    //'/''   (2) Interpolated'''
     '    //'/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) THEN
          IF(TRSF_ACTN_WAVE_INTERPOLATE) THEN
            IDATA(1)=2
          ELSE
            IDATA(1)=1
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TRSF_ACTN_WAVE_INTERPOLATE=IDATA(1).EQ.2

        IF(TRSF_ACTN_WAVE_INTERPOLATE) THEN
          FORMAT='(/$,'' Enter the number of first surface elements '
     '      //'[1]: '',I7)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=NELIST3(0)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NEM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NELIST3(0)=IDATA(1)
          CALL ASSERT(NELIST3(0).LE.NEM,'>>Increase NEM',ERROR,*9999)
          WRITE(CHAR1,'(I7)') NELIST3(0)
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//' first surface '
     '      //'elements: '',2I7,/:(3X,9I7))'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) THEN
            DO nolist=1,NELIST3(0)
              IDATA(nolist)=NELIST3(nolist)
            ENDDO
          ENDIF
          CDATA(1)='ELEMENTS' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      NELIST3(0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '      1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nolist=1,NELIST3(0)
              NELIST3(nolist)=IDATA(nolist)
            ENDDO
          ENDIF

        ENDIF

      ENDIF

      IF(KTYP95.EQ.3) THEN

        FILEIP=.FALSE.
        NOQUES=0
        FORMAT='(/$,'' Enter the number of outer surface nodes '
     '    //'[1]: '',I7)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=NPLIST5(0)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,10000,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NPLIST5(0)=IDATA(1)
        CALL ASSERT(NPLIST5(0).LE.NP_R_M,'>>Increase NP_R_M',
     '    ERROR,*9999)
        WRITE(CHAR1,'(I7)') NPLIST5(0)
        CALL STRING_TRIM(CHAR1,IBEG,IEND)
        FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//' outer surface '
     '    //'nodes: '',2I7,/:(3X,9I7))'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) THEN
          DO nolist=1,NPLIST5(0)
            IDATA(nolist)=NPLIST5(nolist)
          ENDDO
        ENDIF
        CDATA(1)='NODES' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    NPLIST5(0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '    IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO nolist=1,NPLIST5(0)
            NPLIST5(nolist)=IDATA(nolist)
          ENDDO
        ENDIF
        FORMAT='($,'' Enter the region number of these nodes '
     '    //'[1]: '',I5)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=TRSF_NR_OUTER
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NRT,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TRSF_NR_OUTER=IDATA(1)
      ELSE
        TRSF_NR_OUTER=TRSF_NR_SECOND
        NPLIST5(0)=NPLIST4(0)
        DO nolist=1,NPLIST4(0)
          NPLIST5(nolist)=NPLIST4(nolist)
        ENDDO
      ENDIF

      IF(KTYP94.EQ.1) THEN
C       Need to set FIX and YP for direct method, since this method
C       involves solving one forward problem for each
C       inner surface/heart dof.
C       Initialise FIX,YP and temporarily store current information
C       in iy=6 position.
C
CC AJPs 10/11/97
C       If one wants to have the reference value specified at the
C       outer surface (and e.g. use Salu and LU in the solution) then
C       one needs to leave the outer surface bcs alone.  In the direct
C       evaluation of the transfer matrix, GRR is constructed explicitly
C       and this only allows zero values to be fixed on the outer
C       surface (could be generalised later if needed).



C LKC 15-JUL-2000 Do not want to free (set FIX to FALSE) if we
C   are doing an analytic problem and want to solve the whole
C   thing in one hit. If creating a transfer matrix then we
C   want to free all nodes except one node at a time
C

CC JMBs 04-AUG-2000 Added Check on NIYM
        CALL ASSERT(NIYM.GE.6,'>>Increase NIYM >= 6',ERROR,*9999)
CC JMBe
        IF(.NOT.CALL_ANAL_TRSF) THEN
          DO nr=1,NRT
            DO nc=1,NCT(nr,nx)
              DO no_nynr=1,NYNR(0,0,nc,nr)
                ny=NYNR(no_nynr,0,nc,nr)
                IF(NPNY(0,ny,0,nx).EQ.1) THEN
                  np=NPNY(4,ny,0,nx)
C                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF
                YP(ny,6)=YP(ny,1)
                YP(ny,1)=0.0d0
                FIX(ny,5)=FIX(ny,1) !What is fix(ny,5) used for?
                FIX(ny,1)=.FALSE.
              ENDDO !no_nynr
            ENDDO !nc
          ENDDO !nr
        ELSE
          DO nr=1,NRT
            DO nc=1,NCT(nr,nx)
              DO no_nynr=1,NYNR(0,0,nc,nr)
                ny=NYNR(no_nynr,0,nc,nr)
                IF(NPNY(0,ny,0,nx).EQ.1) THEN
                  np=NPNY(4,ny,0,nx)
C                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF
                YP(ny,6)=YP(ny,1)
C Do not set
C                YP(ny,1)=0.0d0
                FIX(ny,5)=FIX(ny,1) !What is fix(ny,5) used for?
                FIX(ny,1)=.FALSE.
              ENDDO !no_nynr
            ENDDO !nc
          ENDDO !nr
        ENDIF

C MLB 11-feb-98 fix for AJP
C This ensureds the outer surface is set to a flux of 0
C if a heart node from outer region copies solution back
C
        DO no_nynr=1,NYNR(0,0,1,TRSF_NR_OUTER)
          ny=NYNR(no_nynr,0,1,TRSF_NR_OUTER)
          IF(NPNY(0,ny,0,nx).EQ.1) THEN
            np=NPNY(4,ny,0,nx)
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            IF(INLIST(np,NPLIST5(1),NPLIST5(0),N1LIST)) THEN
              !np is on the outer surface
              IF(FIX(ny,1))THEN
                CALL ASSERT(DABS(YP(ny,1)).LE.RDELTA,
     '            '>>Outer surface bcs must be zero',ERROR,*9999)
              ELSE
                YP(ny,1)=YP(ny,6)
                FIX(ny,1)=FIX(ny,5)
              ENDIF !fix
            ENDIF !inlist
          ENDIF !NPNY
        ENDDO !no_nynr

C       Set up FIX array for first surface (heart) dof
C       Only needs to be done for the single layer case.  For the
C       double layer case a special manipulation of the governing
C       equations is performed in evtrsf to account for the fact that
C       one is solving a Poisson equation inside the heart, but one has
C       not explicitly constructed the source term in GKK.
        IF(KTYP95.EQ.1) THEN !Single layer transfer matrix
          nr=TRSF_NR_FIRST
          DO nolist=1,NPLIST3(0) !list of first surface nodes
            np=NPLIST3(nolist)
            DO nhx=1,NHP(np,nr)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,1,nr)
                DO nk=1,MAX(NKH(nh,np,1,nr)-KTYP93(1,nr),1)
                  ny=NYNP(nk,nv,nh,np,0,1,nr)
C                 global ny # for the first surface (heart) potential
                  FIX(ny,1)=.TRUE.
                ENDDO !nk
              ENDDO !nv
            ENDDO !nh
          ENDDO !nolist (np)
        ENDIF !ktyp95


C This doing the FIX for all second surface nodes. Should
C   covered by ipinit
C       Set up FIX and YP for outer surface (body) dof


        nr=TRSF_NR_OUTER
        DO nolist=1,NPLIST5(0) !list of body nodes
          np=NPLIST5(nolist)
          DO nhx=1,NHP(np,nr)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,2,nr)
              DO nk=1,MAX(NKH(nh,np,2,nr)-KTYP93(2,nr),1)
                ny=NYNP(nk,nv,nh,np,0,2,nr)
C               global ny # for the second surface (torso) flux
                FIX(ny,1)=.TRUE.
              ENDDO !nk
            ENDDO !nv
          ENDDO !nh
        ENDDO !nolist (np)


      ENDIF !KTYP94

      CALL EXITS('IPTRSF')
      RETURN
 9999 CALL ERRORS('IPTRSF',ERROR)
      CALL EXITS('IPTRSF')
      RETURN 1
      END



