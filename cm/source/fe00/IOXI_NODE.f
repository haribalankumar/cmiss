      SUBROUTINE IOXI_NODE(IUNIT,NEP,NITB,NPNODE,XIP,nr,ERROR,*)

C#### Subroutine: IOXI_NODE
C###  Description:
C###    IOXI reads or writes xi positions of data points in region nr
C###    within a host mesh.
C**** Note: NDTOLD is 0 if only one data set read in
C****          "    " previous value of NDT if another data set has been
C****                 read in by IODATA

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IUNIT,NEP(NPM),NITB,NPNODE(0:NP_R_M,0:NRM)
      REAL*8 XIP(NIM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IEND,ni,np,npp,nr
      CHARACTER FMT*100

      CALL ENTERS('IOXI_NODE',*9999)
C new MPN 7Nov97: more accuracy required
      FMT='(2I5,4E25.16)'
C old      FMT='(2I5,4E18.10)'
      IF(IOTYPE.EQ.2) THEN
        DO np=1,NPNODE(0,nr)
          npp=NPNODE(np,nr)
          READ(IUNIT,*,END=200) npp,NEP(npp),(XIP(ni,npp),ni=1,NITB)
C LKC 7-AUG-2002 Putting in a warning, hopefully it doesn't slow down
C  the code too much
          IF(NEP(npp).EQ.0) THEN

            IEND=0
            CALL APPENDC(IEND,' WARNING: no element for node ',
     '        OP_STRING(1))
            CALL APPENDI(IEND,npp,OP_STRING(1))
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)

          ENDIF

        ENDDO
 200    CONTINUE

      ELSE IF(IOTYPE.EQ.3) THEN
        DO np=1,NPNODE(0,nr)
          npp=NPNODE(np,nr)
          WRITE(IUNIT,FMT) npp,NEP(npp),(XIP(ni,npp),ni=1,NITB)
        ENDDO
      ENDIF

      CALL EXITS('IOXI_NODE')
      RETURN
 9999 CALL ERRORS('IOXI_NODE',ERROR)
      CALL EXITS('IOXI_NODE')
      RETURN 1
      END


