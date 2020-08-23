      SUBROUTINE CAISOC(ISEG,ISISOC,STRING,ERROR,*)

C#### Subroutine: CAISOC
C###  Description:
C###    CAISOC cancels isochrone segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISISOC(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('CAISOC',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel isochrones;s
C###  Description:
C###    Cancel isochrone segments.
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAISOC',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
C ISISOC wrongly dimensioned AJP 16/1/96
C            DO nr=1,NRT
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                nb=NBJ(1,ne)
C                DO ng=1,NGT(nb)
C                  IF(ISISOC(iw,ng,ne).GT.0) THEN
c                    CALL DELETE_SEGMENT(ISISOC(iw,ng,ne),ISEG,iw,
c     '                ERROR,*9999)
c                  ENDIF
c                ENDDO
c              ENDDO
c            ENDDO
            IF(ISISOC(iw).GT.0) THEN
              CALL DELETE_SEGMENT(ISISOC(iw),ISEG,iw,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CAISOC')
      RETURN
 9999 CALL ERRORS('CAISOC',ERROR)
      CALL EXITS('CAISOC')
      RETURN 1
      END


