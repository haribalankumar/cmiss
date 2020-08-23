      SUBROUTINE IPREFE(ERROR,*)

C#### Subroutine: IPREFE
C###  Description:
C###    IPREFE inputs reference electrode/node location
C###  See-Also APREFE

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ref000.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,ICHAR,INFO,no_ref_loc,NOQUES
      CHARACTER CHAR1*5
      LOGICAL FILEIP

      CALL ENTERS('IPREFE',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='($,'' The total number of reference locations is'
     '    //' [1]: '',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=NP_REF_LOC(0)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NPM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NP_REF_LOC(0)=IDATA(1)
      CALL ASSERT(NP_REF_LOC(0).LE.10,'>> Increase NP_REF_LOC size',
     '  ERROR,*9999)
      WRITE(CHAR1,'(I2)') NP_REF_LOC(0)
      CALL STRING_TRIM(CHAR1,IBEG,IEND)
      FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//
     '  ' node numbers of the reference locations : '',10(1X,I5))'
      DO no_ref_loc=1,NP_REF_LOC(0)
        IDEFLT(no_ref_loc)=NP_REF_LOC(no_ref_loc)
      ENDDO
      IF(IOTYPE.EQ.3) THEN
        DO no_ref_loc=1,NP_REF_LOC(0)
          IDATA(no_ref_loc)=NP_REF_LOC(no_ref_loc)
        ENDDO
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  NP_REF_LOC(0),
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        DO no_ref_loc=1,NP_REF_LOC(0)
          NP_REF_LOC(no_ref_loc)=IDATA(no_ref_loc)
        ENDDO
      ENDIF

      CALL EXITS('IPREFE')
      RETURN
 9999 CALL ERRORS('IPREFE',ERROR)
      CALL EXITS('IPREFE')
      RETURN 1
      END


