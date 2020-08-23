      SUBROUTINE CAGAUS(ISEG,ISGAUS,NBJ,NEELEM,STRING,ERROR,*)

C#### Subroutine: CAGAUS
C###  Description:
C###    CAGAUS cancels Gauss point segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER ISEG(*),ISGAUS(NWM,NGM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,ne,ng,noelem,nr

      CALL ENTERS('CAGAUS',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel gauss;s
C###  Description:
C###    Cancel gauss segment.
C###    If the s qualified is present the field segments
C###    are canceled on specified workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAGAUS',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        DO iw=1,2*NJT-3+IMAP
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO ng=1,NGT(NBJ(1,ne))
                  IF(ISGAUS(iw,ng,ne).GT.0) THEN
                    CALL DELETE_SEGMENT(ISGAUS(iw,ng,ne),ISEG,iw,
     '                ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CAGAUS')
      RETURN
 9999 CALL ERRORS('CAGAUS',ERROR)
      CALL EXITS('CAGAUS')
      RETURN 1
      END


