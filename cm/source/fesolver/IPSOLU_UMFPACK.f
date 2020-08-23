      SUBROUTINE IPSOLU_UMFPACK(UMF_PARAM,NOQUES,FILEIP,ERROR,*)

C#### Subroutine: IPSOLU_UMFPACK
C###  Description:
C###    IPSOLU_UMFPACK gets solution parameters for the Umfpack solver

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER NOQUES
      REAL*8 UMF_PARAM(40)
      CHARACTER ERROR*(*)
      LOGICAL FILEIP
!     Local Variables
      INTEGER INFO


      CALL ENTERS('IPSOLU_UMFPACK',*9999)

C     Optionally reset the settings to their default values.
      IF(IOTYPE.NE.3) CALL UMFPACK4_RESET(UMF_PARAM,ERROR,*9999)

C     Start reading in the values ....
      RDEFLT(1)=0.1D0
      FORMAT='($,'' Specify the pivot threshold [0.1D0]: '',D11.4)'
C DAH 07-01-2003 Added in a value for writing files.
      IF(IOTYPE.EQ.3) RDATA(1)=UMF_PARAM(1)
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,99999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,1.0D0,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) UMF_PARAM(1)=RDATA(1)

      CALL EXITS('IPSOLU_UMFPACK')
      RETURN
 9999 CALL ERRORS('IPSOLU_UMFPACK',ERROR)
      CALL EXITS('IPSOLU_UMFPACK')
      RETURN 1
      END


