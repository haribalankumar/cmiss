      SUBROUTINE SHSECT(ISEG,ISSECT,STRING,ERROR,*)

C#### Subroutine: SHSECT
C###  Description:
C###    SHSECT shows section segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSECT(NGRSEGM),NSLIST(0:20)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,N3CO,nolist,nosect
      LOGICAL CBBREV

      CALL ENTERS('SHSECT',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show section
C###  Parameter:    <at (all/SECTION#s)[all]>
C###    Specify the section(s) to show. The all keyword shows all
C###    defined sections.
C###  Description:
C###    Make the specified section segments visible.

        OP_STRING(1)=STRING(1:IEND) //' <at (all/SECTION#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHSECT',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'AT',1,noco+1,noco+1,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),20,NSLIST(0),NSLIST(1),ERROR,*9999)
        ELSE
          NSLIST(0)=NTSECT
          DO nolist=1,NSLIST(0)
            NSLIST(nolist)=nolist
          ENDDO
        ENDIF
        IF(IWKS(11).GT.0) THEN
C LKC use iw 2-MAY-1998
C          CALL ACWK(11,1,ERROR,*9999)
          iw=11
          CALL ACWK(iw,1,ERROR,*9999)
          DO nolist=1,NSLIST(0)
            nosect=NSLIST(nolist)
            IF(ISEG(ISSECT(nosect)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISSECT(nosect),'VISIBLE',ERROR,*9999)
            ELSE IF(ISEG(ISSECT(nosect)).EQ.0) THEN
              WRITE(OP_STRING,
     '          '('' >>Section is not defined at '',I3)') nosect
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
C LKC use iw 2-MAY-1998
C          CALL DAWK(11,1,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SHSECT')
      RETURN
 9999 CALL ERRORS('SHSECT',ERROR)
      CALL EXITS('SHSECT')
      RETURN 1
      END


