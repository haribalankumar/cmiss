      SUBROUTINE LDLDR(LD,NDDATA,NDDL,NDLT,ERROR,*)

C#### Subroutine: LDLDR
C###  Description:
C###    LDLDR calulates the element data point arrays (NDDL and NDLT)
C###    for a data region from the array LD.

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER LD(NDM),NDDATA(0:NDM,0:NRM),NDDL(NEM,NDEM),NDLT(NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nd,ne,nodata

      CALL ENTERS('LDLDR',*9999)

      DO ne=1,NEM
        NDLT(ne)=0
      ENDDO !ne

C LKC 24-MAY-1998
      CALL ASSERT(NDDATA(0,DATA_REGION).GT.0,
     '  '>> No data defined',ERROR,*9999)

      DO nodata=1,NDDATA(0,DATA_REGION)
        nd=NDDATA(nodata,DATA_REGION)
        ne=LD(nd)
        IF(ne.NE.0) THEN
          NDLT(ne)=NDLT(ne)+1
          IF(NDLT(ne).LE.NDEM) THEN
            NDDL(ne,NDLT(ne))=nd
          ENDIF
        ENDIF
      ENDDO !nodata (nd)

      DO ne=1,NEM
        CALL ASSERT(NDLT(ne).LE.NDEM,'>>Increase NDEM',ERROR,*9999)
      ENDDO !ne

      CALL EXITS('LDLDR')
      RETURN
 9999 CALL ERRORS('LDLDR',ERROR)
      CALL EXITS('LDLDR')
      RETURN 1
      END




