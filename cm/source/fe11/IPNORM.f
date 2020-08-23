      SUBROUTINE IPNORM(NELIST,NW,ERROR,*)

C#### Subroutine: IPNORM
C###  Description:
C###    IPNONL inputs normal reversal parameters for boundary
C###    elements.
C***  Created by Martin Buist April 1998

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'inout00.cmn'

!     Parameter List
      INTEGER NELIST(0:NEM),NW(NEM,3,NXM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ICHAR,INFO,n,ne,NOQUES,nx,nxx
      LOGICAL FILEIP

      CALL ENTERS('IPNORM',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(IOTYPE.NE.3) THEN
 802    FORMAT='($,'' Enter element #s [EXIT]: '',I5)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

        IF(IDATA(1).NE.0) THEN !not default exit
          NELIST(0)=IDATA(0)
          DO n=1,IDATA(0)
            NELIST(n)=IDATA(n)
          ENDDO
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO n=1,NELIST(0)
              ne=NELIST(n)
              NW(ne,3,nx)=1
            ENDDO
          ENDDO
          GO TO 802
        ENDIF

      ELSE IF(IOTYPE.EQ.3) THEN
        DO ne=1,NET(0)
          IF(NW(ne,3,1).EQ.1) THEN
            FORMAT='($,'' Enter element #s [EXIT]: '',I6)'
            IDATA(1)=ne
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          ENDIF
        ENDDO

        FORMAT='($,'' Enter element #s [EXIT]: '',I6)'
        IDATA(1)=0
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      ENDIF

      CALL EXITS('IPNORM')
      RETURN
 9999 CALL ERRORS('IPNORM',ERROR)
      CALL EXITS('IPNORM')
      RETURN 1
      END


