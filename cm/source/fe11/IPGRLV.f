      SUBROUTINE IPGRLV(NLQ,NQLIST,ERROR,*)

C#### Subroutine: IPGRLV
C###  Description:
C###    IPGRLV inputs parameters NLQ(nq) to identify grid pt adaptive level
C###     for residual calculation.
C**** Written by PJH, June 1999

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NLQ(NQM),NQLIST(0:NQM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ICHAR,INFO,NOQUES,n,na,nq
      CHARACTER CFROMI*2,CHAR2*2
      LOGICAL FILEIP

      CALL ENTERS('IPGRLV',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      DO nq=1,NQT !initialise arrays
        NQLIST(nq)=0
        NLQ(nq)=0
      ENDDO

      DO na=1,NMGT !loop over grid levels
        IF(IOTYPE.NE.3) THEN
          CHAR2=CFROMI(na,'(I2)')
 100      FORMAT='($,'' Enter adaptive level '//CHAR2
     '      //' grid#s/name [EXIT]: '',I5)'
          CDATA(1)='GRIDS' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      0,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,1,NQT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(NQLIST(1).NE.0) THEN !not default exit
            DO n=1,NQLIST(0)
              NLQ(NQLIST(n))=na
            ENDDO !n
            GO TO 100
          ENDIF !NQLIST(1)
        ELSE IF(IOTYPE.EQ.3) THEN
        ENDIF
      ENDDO !na

      CALL EXITS('IPGRLV')
      RETURN
 9999 CALL ERRORS('IPGRLV',ERROR)
      CALL EXITS('IPGRLV')
      RETURN 1
      END


