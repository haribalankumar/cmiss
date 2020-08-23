      SUBROUTINE HISTRA(ISEG,ISSTRA,NEELEM,STRING,ERROR,*)

C#### Subroutine: HISTRA
C###  Description:
C###    HISTRA hides element strain segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSTRA(NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,ne,noelem,nostra,nr

      CALL ENTERS('HISTRA',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide strain
C###  Description:
C###    Hide principal strain vectors.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HISTRA',ERROR,*9999)
      ELSE
        iw=2*NJT-3
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nostra=1,NTSTRA
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISEG(ISSTRA(ne,nostra)).EQ.2) THEN
                  CALL VISIB(iw,ISEG,ISSTRA(ne,nostra),'INVISIBLE',
     '              ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('HISTRA')
      RETURN
 9999 CALL ERRORS('HISTRA',ERROR)
      CALL EXITS('HISTRA')
      RETURN 1
      END


