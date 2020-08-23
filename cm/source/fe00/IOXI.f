      SUBROUTINE IOXI(FD,IUNIT,NDP,LD,LD_NP,NITB,NRE,NRLIST,XID,
     '  NODE_ASSOC,REFERENCE,ERROR,*)

C#### Subroutine: IOXI
C###  Description:
C###    IOXI reads or writes xi positions of data points.
C**** Note: NDTOLD is 0 if only one data set read in
C****          "    " previous value of NDT if another data set has been
C****                 read in by IODATA

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ioda00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'ktyp00.cmn'

!     Parameter List
      INTEGER FD(NDM),IUNIT,LD(NDM),LD_NP(NDM),NDP(NDM),NITB,NRE(NEM),
     &  NRLIST(0:NRM)
      REAL*8 XID(NIM,NDM)
      LOGICAL NODE_ASSOC,REFERENCE
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER DUMMY,IEND,nd,ndd,ne,ni,nr
      CHARACTER FMT*100,FMT1*100
      LOGICAL FAIL    

!     Functions
      LOGICAL INLIST

      CALL ENTERS('IOXI',*9999)

      DUMMY=0
      FAIL=.FALSE.

      IF(IOTYPE.EQ.2) THEN
        DO nd=NDTOLD+1,NDM
C PM 14-AUG-02 : Added for face projection
          IF(KTYP11.EQ.2) THEN !face projection
            READ(IUNIT,*,END=200) ndd,LD(nd),FD(nd),
     '        (XID(ni,nd),ni=1,NITB)
          ELSE IF(NODE_ASSOC)THEN !includes coupled nodes
cc            READ(IUNIT,*,END=200) ndd,LD(nd),LD_NP(nd),
            FMT='(2I7,4E25.16)'
            READ(IUNIT,FMT,END=200) NDP(nd),LD(nd),
     &        (XID(ni,nd),ni=1,NITB)
          ELSE
            READ(IUNIT,*,END=200) ndd,LD(nd),(XID(ni,nd),ni=1,NITB)
          ENDIF
C LKC 22-JUL-1999 Adding check that the element is associated with
C the correct region
C JWF 9/05/03 Extended to check multiple regions
          CALL ASSERT(LD(nd).LE.NEM,
     '      '>>LD(nd) is greater than NEM',ERROR,*9999)
          ne=LD(nd)  
          IF(ne.NE.0) THEN
            nr=NRE(ne)
            IF(.NOT.INLIST(nr,NRLIST(1),NRLIST(0),DUMMY)) THEN
              FAIL=.TRUE.
            ENDIF
            IF(FAIL) THEN
C             output a warning  
              IEND=0
              CALL APPENDC(IEND,'WARNING: Element # ',OP_STRING(1))
              CALL APPENDI(IEND,LD(nd),OP_STRING(1))
              CALL APPENDC(IEND,' associated with data pt ',
     '          OP_STRING(1))
              CALL APPENDI(IEND,ndd,OP_STRING(1))
              CALL APPENDC(IEND,' does not exist',OP_STRING(1))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              FAIL=.FALSE.    
            ENDIF !FAIL
          ENDIF !ne 
        ENDDO
 200    CONTINUE
        IF(NDT.NE.nd-1) THEN
          OP_STRING(1)='WARNING: #xi pts NE #data pts, Setting NDT=#xi'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          NDT=nd-1
        ENDIF
      ELSE IF(IOTYPE.EQ.3) THEN
C new MPN 9Mar08: more digits required
        FMT='(2I7,4E25.16)'
        FMT1='(3I7,3E25.16)'
C old        FMT='(2I5,4E25.16)'
C old        FMT1='(3I5,3E25.16)'
        DO nd=1,NDT
C PM 14-AUG-02 : Added for face projection
          IF(KTYP11.EQ.2) THEN !face projection
            WRITE(IUNIT,FMT1) nd,LD(nd),FD(nd),(XID(ni,nd),ni=1,NITB)
          ELSE IF(NODE_ASSOC)THEN !includes coupled nodes
            FMT='(2I7,4E25.16)'
            WRITE(IUNIT,FMT) NDP(nd),LD(nd),(XID(ni,nd),ni=1,NITB)
          ELSE
            IF(REFERENCE) THEN
              WRITE(IUNIT,FMT) NDP(nd),LD(nd),(XID(ni,nd),ni=1,NITB)
            ELSE
              WRITE(IUNIT,FMT) nd,LD(nd),(XID(ni,nd),ni=1,NITB)
            ENDIF
          ENDIF
        ENDDO !nd
        
        IF(NDT.GT.9999999) THEN
          WRITE(OP_STRING,'('' >>Warning!! fe00/IOXI.f format cannot '
     &      //'handle  >9999999 data points.'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
C        CALL ASSERT(NDT.LT.10000000,
C     &    '>>Format cannot handle >9999999 points',ERROR,*9999)
      ENDIF !iotype
      CALL EXITS('IOXI')
      RETURN
 9999 CALL ERRORS('IOXI',ERROR)
      CALL EXITS('IOXI')
      RETURN 1
      END


