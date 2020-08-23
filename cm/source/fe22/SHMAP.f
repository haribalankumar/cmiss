      SUBROUTINE SHMAP(ISEG,ISMAP,STRING,ERROR,*)

C#### Subroutine: SHMAP
C###  Description:
C###    SHMAP shows increment segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISMAP(NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,iw,N3CO,nomap
      LOGICAL CBBREV

      CALL ENTERS('SHMAP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show map
C###  Parameter:    <at (last/MAP#)[last]>
C##     Specify map segment to show
C###  Description:
C###    Make the specified map segment visible.

        OP_STRING(1)=STRING(1:IEND) //' <at (last/MAP#)[last]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHMAP',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
          nomap=IFROMC(CO(N3CO+1))
        ELSE
          nomap=NTMAP
        ENDIF
C LKC 24-APR-1998 init variables iw
        iw=1

        IF(ISEG(ISMAP(nomap)).EQ.1) THEN
          CALL ACWK(4,1,ERROR,*9999)
          CALL VISIB(iw,ISEG,ISMAP(nomap),'VISIBLE',ERROR,*9999)
          CALL DAWK(4,1,ERROR,*9999)
        ELSE IF(ISEG(ISMAP(nomap)).EQ.0) THEN
          WRITE(OP_STRING,'('' >>Map is not defined at '',I3)') nomap
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SHMAP')
      RETURN
 9999 CALL ERRORS('SHMAP',ERROR)
      CALL EXITS('SHMAP')
      RETURN 1
      END


