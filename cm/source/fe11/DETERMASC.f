      SUBROUTINE  DETERMASC(nr,CELL_ICQS_NAMES,CELL_RCQS_NAMES
     &     ,CELL_YQS_NAMES,ASCOMP,ASC_ARRAYNAME,ASC_CELLVARINDEX,ERROR,
     &     *)

C#### Subroutine: DETERMASC
C###  Description:
C###    <html><pre> DETERMASC determines for the in IPACTI specified 
C###       CellML variable the type and its position (CVELLVARINDEX) 
C###       within the respective CellML value vector
C###    </pre></html>

C
C     WRITTEN BY OR 15-08-06
C
      IMPLICIT NONE
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(CELL_NAME_LENGTH),
     &     CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     &     CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),ERROR*200
      CHARACTER ASCOMP(NRM)*500
      INTEGER nr
!     Local Variables
      INTEGER nqv,niq,IBEG1,IEND1,IBEG2,IEND2,ASC_CELLVARINDEX(NRM)
      LOGICAL FOUNDVAR
      CHARACTER ASCOMP_UC*500,CELL_NAME_UC*(CELL_NAME_LENGTH),
     &     ASC_ARRAYNAME(NRM)*8

      CALL ENTERS('DETERMASC',*9999)

      CALL CUPPER(ASCOMP(nr),ASCOMP_UC)
      CALL STRING_TRIM(ASCOMP_UC,IBEG1,IEND1)
      FOUNDVAR=.FALSE.

      nqv=0
      DO WHILE(nqv.LT.NQVT.AND..NOT.FOUNDVAR)
        nqv=nqv+1
C     Search CELL_YQS_NAMES
        niq=0
        DO WHILE(niq.LT.NIQST.AND..NOT.FOUNDVAR)
          niq=niq+1
          CALL CUPPER(CELL_YQS_NAMES(niq,nqv),CELL_NAME_UC)
          CALL STRING_TRIM(CELL_NAME_UC,IBEG2,IEND2)
          IF(ASCOMP_UC(IBEG1:IEND1).EQ.CELL_NAME_UC(IBEG2:IEND2)) THEN
            ASC_ARRAYNAME(nr)='YQS'
            ASC_CELLVARINDEX(nr)=niq 
            FOUNDVAR=.TRUE. 
          ENDIF 
        ENDDO                   !while niq.LE.NIQST
C     Search CELL_ICQS_NAMES
        niq=0
        DO WHILE(niq.LT.NQIT.AND..NOT.FOUNDVAR)
          niq=niq+1
          CALL CUPPER(CELL_ICQS_NAMES(niq,nqv),CELL_NAME_UC)
          CALL STRING_TRIM(CELL_NAME_UC,IBEG2,IEND2)
          IF(ASCOMP_UC(IBEG1:IEND1).EQ.CELL_NAME_UC(IBEG2:IEND2)) THEN
            ASC_ARRAYNAME(nr)='ICQS'
            ASC_CELLVARINDEX(nr)=niq
            FOUNDVAR=.TRUE.
          ENDIF
        ENDDO                   !while niq.LE.NQIT
C     Search CELL_ICQS_NAMES
        niq=0
        DO WHILE(niq.LT.NQRT.AND..NOT.FOUNDVAR)
          niq=niq+1
          CALL CUPPER(CELL_RCQS_NAMES(niq,nqv),CELL_NAME_UC)
          CALL STRING_TRIM(CELL_NAME_UC,IBEG2,IEND2)
          IF(ASCOMP_UC(IBEG1:IEND1).EQ.CELL_NAME_UC(IBEG2:IEND2)) THEN
            ASC_ARRAYNAME(nr)='RCQS'
            ASC_CELLVARINDEX(nr)=niq
            FOUNDVAR=.TRUE.
          ENDIF
        ENDDO                   !while niq.LE.NQRT
      ENDDO                     !while nqv.LE.NQVT
      CALL ASSERT(FOUNDVAR,'Cell variable was not'
     &     //' found in YQS, ICQS or RCQS. Make sure you'
     &     //' have already defined your CellML variables.'
     &     //' You should move the cmd: FEM define acti'
     &     //' after the Cellml initialisation.',ERROR,*9999)
      

      CALL EXITS('DETERMASC')
      RETURN
 9999 CALL ERRORS('DETERMASC',ERROR)
      CALL EXITS('DETERMASC')
      RETURN 1
      END


