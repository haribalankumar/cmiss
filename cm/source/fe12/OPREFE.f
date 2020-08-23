      SUBROUTINE OPREFE(ERROR,*)

C#### Subroutine: OPREFE
C###  Description:
C###    OPREGI outputs reference node(s)/location(s)

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'ref000.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 DUMMY(1)

      CALL ENTERS('OPREFE',*9999)

      DUMMY(1)=0.D0
      WRITE(OP_STRING,'(''  The total number of reference locations'
     '  //' is '',I2)') NP_REF_LOC(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' The node numbers of the reference'
     '  //' locations are:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      CALL WRITE_LONG(INTTYPE,1,1,IOFI,NP_REF_LOC(0),10,10,
     '    NP_REF_LOC(1),DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
      CALL EXITS('OPREFE')
      RETURN
 9999 CALL ERRORS('OPREFE',ERROR)
      CALL EXITS('OPREFE')
      RETURN 1
      END


