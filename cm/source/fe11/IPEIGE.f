      SUBROUTINE IPEIGE(NYNR,FIX,ERROR,*)

C#### Subroutine: IPEIGE
C###  Description:
C###    IPEIGE inputs eigenvalue analysis parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'eige00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      INTEGER NYNR(0:NY_R_M,0:NRCM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER ICHAR,INFO,IB4,IE4,no_nynr,NOQUES,ny,NYTOT
      CHARACTER CHAR4*4
      LOGICAL FILEIP

      CALL ENTERS('IPEIGE',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      NYTOT=0
      DO no_nynr=1,NYNR(0,0,1) !loop over global variables
        ny=NYNR(no_nynr,0,1) !is global variable number
        IF(.NOT.FIX(ny,1)) NYTOT=NYTOT+1
      ENDDO
C      FORMAT='($,'' Specify whether lowest(1) or highest(2)'//
C     '  ' eigenvalues are required [1]: '',I1)'
C      IF(IOTYPE.EQ.3) IDATA(1)=KTYP16
C      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C      IF(IOTYPE.NE.3) KTYP16=IDATA(1)

      WRITE(CHAR4,'(I4)') NYTOT
      CALL STRING_TRIM(CHAR4,IB4,IE4)
      IDEFLT(1)=NYTOT
      FORMAT='($,'' Specify the number of eigenpairs required (le ['
     '  //CHAR4(IB4:IE4)//']) : '',I3)'
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP17
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NYTOT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP17=IDATA(1)

      CALL ASSERT(KTYP17.LE.NTM,'>>Increase NTM',ERROR,*9999)

      RDEFLT(1)=0.0d0
      FORMAT='($,'' Enter the minimum modal frequency required (kHz) '
     '  //'[0.0]: '',D11.4))'
      IF(IOTYPE.EQ.3) RDATA(1)=DSQRT(EIGEN_LIMITS(1))/(2.0d0*PI)
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) EIGEN_LIMITS(1)=(RDATA(1)*2.0d0*PI)**2

      RDEFLT(1)=1.0d0
      FORMAT='($,'' Enter the maximum modal frequency required (kHz) '
     '  //'[1.0]: '',D11.4))'
      IF(IOTYPE.EQ.3) RDATA(1)=DSQRT(EIGEN_LIMITS(2))/(2.0d0*PI)
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) EIGEN_LIMITS(2)=(RDATA(1)*2.0d0*PI)**2

      CALL EXITS('IPEIGE')
      RETURN
 9999 CALL ERRORS('IPEIGE',ERROR)
      CALL EXITS('IPEIGE')
      RETURN 1
      END


