      SUBROUTINE OPELEM_INTERFACE(NBJ,NEELEM,NELIST,NP_INTERFACE,
     '  NPNE,nr_interface,NRLIST,ERROR,*)

C#### Subroutine: OPELEM_INTERFACE
C###  Description:
C###    LIELEM_INTEFACE lists elements which have been defined
C###    by nodes which are interface nodes in # regions specified
C###    by the variable nr_interface.

C*** Created by Leo Cheng 1-JUL-1999

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),nr_interface,NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),NRLIST(0:NRM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER IEND,nb,ne,nonode,noelem,noregion,np,nr
      REAL*8 DUMMY(1)
      LOGICAL FOUND

! Function
      LOGICAL INLIST


      CALL ENTERS('OPELEM_INTERFACE',*9999)

      DUMMY(1)=0.D0 !initialise dummy

      NELIST(0)=0
      DO noregion=1,NRLIST(0)
        nr=NRLIST(noregion)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          FOUND=.FALSE.
          DO nonode=1,NNT(nb)
            np=NPNE(nonode,nb,ne)

C 27-SEP-1999 Change so that a specific region can be specified.
C            IF((NP_INTERFACE(np,0).EQ.nr_interface)
C     '        .AND.(FOUND.EQ..FALSE.)) THEN

            IF(INLIST(nr_interface,NP_INTERFACE(np,1),
     '        NP_INTERFACE(np,0),IEND )
     '        .AND..NOT.FOUND) THEN

              NELIST(0)=NELIST(0)+1
              NELIST(NELIST(0))=ne  ! creating a list of elements
              FOUND=.TRUE.
            ENDIF !NP_INTERFACE
          ENDDO !nonode (np)
        ENDDO !noelem (ne)
      ENDDO !nolist (nr)

      WRITE(OP_STRING(1),'('' '')')
      IEND=0
      IF(NELIST(0).GT.0) THEN
        CALL APPENDC(IEND,' #Elements interfacing with region '
     '    ,OP_STRING(2))
        CALL APPENDI(IEND,nr_interface,OP_STRING(2))
        CALL APPENDC(IEND,' is ',OP_STRING(2))
        CALL APPENDI(IEND,NELIST(0),OP_STRING(2))
        CALL APPENDC(IEND,', Elements are: ',OP_STRING(2))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL WRITE_LONG(INTTYPE,1,1,IOFI,NELIST(0),10,10,
     '    NELIST(1),DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
      ELSE
        CALL APPENDC(IEND,' No elements interfacing with region'
     '    ,OP_STRING(2))
        CALL APPENDI(IEND,nr_interface,OP_STRING(2))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPELEM_INTERFACE')
      RETURN
 9999 CALL ERRORS('OPELEM_INTERFACE',ERROR)
      CALL EXITS('OPELEM_INTERFACE')
      RETURN 1
      END


