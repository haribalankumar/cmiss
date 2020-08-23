      SUBROUTINE IPSING(ERROR,*)

C#### Subroutine: IPSING
C###  Description:
C###    IPSING defines the node and element number at which a physical
C###    singularity occurs.  Examples are at the corner of surface
C###    cavities (governed by the modified Helmholtz equation) and
C###    stress singularities in some corners.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'sing00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES,NTOTAL
      LOGICAL CONTINUE,FILEIP

      CALL ENTERS('IPSING',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      CONTINUE=.TRUE.
      IDEFLT(1)=-1
      NTOTAL=0
      DO WHILE(CONTINUE)
        NTOTAL=NTOTAL+1
        FORMAT='(/$,'' Enter a singularity node number [Exit]: '',I2)'
        IF(IOTYPE.EQ.3) THEN
          IF(NTOTAL.EQ.(NPSING(0)+1)) THEN
            IDATA(1)=IDEFLT(1)
          ELSE
            IDATA(1)=NPSING(NTOTAL)
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPT(1),
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IDATA(1).GT.IDEFLT(1)) THEN
          IF(IOTYPE.NE.3) THEN
            NPSING(NTOTAL)=IDATA(1)
            NPSING(0)=NPSING(0)+1
          ENDIF
          FORMAT=
     '    '(/$,'' Enter the corresponding element number [Exit]: '',I2)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=NESING(NTOTAL)
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPT(1),
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NESING(NTOTAL)=IDATA(1)
            NESING(0)=NESING(0)+1
          ENDIF
        ELSE
          CONTINUE=.FALSE.
        ENDIF
      ENDDO
      CALL EXITS('IPSING')
      RETURN
 9999 CALL ERRORS('IPSING',ERROR)
      CALL EXITS('IPSING')
      RETURN 1
      END


