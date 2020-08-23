      SUBROUTINE OPREGI(NEELEM,NPNODE,ERROR,*)

C#### Subroutine: OPREGI
C###  Description:
C###    OPREGI outputs region data.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nr  !SMAR009 22/12/98 noelem,nonode,
      REAL*8 DUMMY(1)

      CALL ENTERS('OPREGI',*9999)

C LKC 2-MAY-1999 Initialise Dummy
      DUMMY(1)=0.D0

      WRITE(OP_STRING,'(''  The total number of regions is '',I1)') NRT
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nr=1,NRT
        WRITE(OP_STRING,'(/'' #Elements in region '',I1,'' is '','
     '    //'I5,'', Elements are:'')') nr,NEELEM(0,nr)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL WRITE_LONG(INTTYPE,1,1,IOFI,NEELEM(0,nr),10,10,
     '    NEELEM(1,nr),DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
        WRITE(OP_STRING,'(/''  #Nodes   in region '',I1,'' is '','
     '    //'I5,'', Nodes are:'')') nr,NPNODE(0,nr)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL WRITE_LONG(INTTYPE,1,1,IOFI,NPNODE(0,nr),10,10,
     '    NPNODE(1,nr),DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
      ENDDO !nr

      CALL EXITS('OPREGI')
      RETURN
 9999 CALL ERRORS('OPREGI',ERROR)
      CALL EXITS('OPREGI')
      RETURN 1
      END


