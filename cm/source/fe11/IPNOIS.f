      SUBROUTINE IPNOIS(ERROR,*)

C#### Subroutine: IPNOIS
C###  Description:
C###    IPNOISE inputs noise parameters for applying to signals

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'nois00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPNOIS',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='(/'' Specify where the noise is to be applied [1]: '''
     '  //'/''   (1) Signal Values'''
     '  //'/''   (2) Electrode Locations'''
     '  //'/''   (3) Both Signals and Electrodes'''
     '  //'/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=DEFNOIS_LOC
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) DEFNOIS_LOC=IDATA(1)


      FORMAT='(/'' Specify the noise type [1]: '''
     '  //'/''   (1) Absolute'''
     '  //'/''  *(2) Relative'''
     '  //'/$,''    '',I1)'

      IF(IOTYPE.EQ.3) IDATA(1)=DEFNOIS_TYPE
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,1,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) DEFNOIS_TYPE=IDATA(1)

C*** Absolute noise
      IF(DEFNOIS_TYPE.EQ.1) THEN !absolute noise

        FORMAT='(/'' Specify distribution type [1]: '''
     '    //'/''   (1) Gaussian '''
     '    //'/''   (2) Uniform - ran0'''
     '    //'/''   (3) Uniform - ran1'''
     '    //'/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NOIS_DIST
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NOIS_DIST=IDATA(1)

C*** Define noise values of signals
        IF(DEFNOIS_LOC.EQ.1.OR.DEFNOIS_LOC.EQ.3) THEN

CC AJPs
C          FORMAT='(/'' Specify absolute noise level (microV) [10]: '''
C     '      //'/$,''    '',I2)'
          FORMAT='(/'' Specify absolute noise level (mV) [1]: '''
     '    //'/$,''    '',D12.4)'
          RDEFLT(1)=1.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=NOIS_LEV
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,10,-99,99,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
CC AJPe
          IF(IOTYPE.NE.3) NOIS_LEV=RDATA(1)
        ENDIF
C*** Define noise values for electrode locations
        IF(DEFNOIS_LOC.EQ.2.OR.DEFNOIS_LOC.EQ.3) THEN
          FORMAT='(/'' Specify absolute displacement [20.0]: '''
     '      //'/$,''    '',F8.2)'
          RDEFLT(1)=20
          IF(IOTYPE.EQ.3) RDATA(1)=NOIS_DISP
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMAX,IMIN,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NOIS_DISP=RDATA(1)
        ENDIF

C*** Relative Noise
      ELSE !relative noise
        CALL ASSERT(DEFNOIS_TYPE.NE.3,
     '      '>>Relative noise not implemented',ERROR,*9999)
C***    Insert Relative noise input here
      ENDIF

      CALL EXITS('IPNOIS')
      RETURN
 9999 CALL ERRORS('IPNOIS',ERROR)
      CALL EXITS('IPNOIS')
      END


