      SUBROUTINE SHISOC(ISEG,ISISOC,STRING,ERROR,*)

C#### Subroutine: SHISOC
C###  Description:
C###    SHISOC shows isochrone segments.

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

      CALL ENTERS('SHISOC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show isochrones
C###  Description:
C###    Make the isochrone segments visible.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHISOC',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
C ISISOC wrongly dimensioned
c            DO nr=1,NRT
c              DO noelem=1,NEELEM(0,nr)
c                ne=NEELEM(noelem,nr)
c                nb=NBJ(1,ne)
c                DO ng=1,NGT(nb)
c                  IF(ISEG(ISISOC(iw,ng,ne)).EQ.1) THEN
c                    CALL VISIB(iw,ISEG,ISISOC(iw,ng,ne),'VISIBLE',
c     '                ERROR,*9999)
c                  ELSE IF(ISEG(ISISOC(iw,ng,ne)).EQ.0) THEN
c                    WRITE(OP_STRING,
c     '                 '('' >>Isochrone at ng='',I3,'' ne='',I3,'
c     '                //''' is not defined on '',I1)') ng,ne,iw
c                     CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c                  ENDIF
c                ENDDO
c              ENDDO
c            ENDDO
            IF(ISEG(ISISOC(iw)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISISOC(iw),'VISIBLE',ERROR,*9999)
            ELSE IF(ISEG(ISISOC(iw)).EQ.0) THEN
              WRITE(OP_STRING,
     '          '('' >>Isochrone is not defined on '',I1)') iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHISOC')
      RETURN
 9999 CALL ERRORS('SHISOC',ERROR)
      CALL EXITS('SHISOC')
      RETURN 1
      END


