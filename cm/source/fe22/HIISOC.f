      SUBROUTINE HIISOC(ISEG,ISISOC,STRING,ERROR,*)

C#### Subroutine: HIISOC
C###  Description:
C###    HIISOC hides isochrone segments.

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

      CALL ENTERS('HIISOC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide isochrones
C###  Description:
C###    Hide isochrones.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIISOC',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
C ISISOC wrongly dimensioned AJP 16/1/96
c            DO nr=1,NRT
c              DO noelem=1,NEELEM(0,nr)
c                ne=NEELEM(noelem,nr)
c                nb=NBJ(1,ne)
c                DO ng=1,NGT(nb)
c                  IF(ISEG(ISISOC(iw,ng,ne)).EQ.2) THEN
c                    CALL VISIB(iw,ISEG,ISISOC(iw,ng,ne),'INVISIBLE',
c     '                ERROR,*9999)
c                  ENDIF
c                ENDDO
c              ENDDO
c            ENDDO
            IF(ISEG(ISISOC(iw)).EQ.2) THEN
              CALL VISIB(iw,ISEG,ISISOC(iw),'INVISIBLE',
     '          ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIISOC')
      RETURN
 9999 CALL ERRORS('HIISOC',ERROR)
      CALL EXITS('HIISOC')
      RETURN 1
      END


