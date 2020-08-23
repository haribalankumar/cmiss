      SUBROUTINE EVALASC(ne,ng,NQNE,YQS,RCQS_SPATIAL,ICQS_SPATIAL,
     &  ASC_NAME,ASC_CVI,ASC,ERROR,*)

C#### Subroutine: EVALASC
C###  Description:
C###    <html><pre> EVALASC
C###       Evaluates the ACTIVE_STRESS value at a particular gauss point
C###       within element ne depending on the inputs in IPACTI, e.g.
C###       CellML function name.
C###    </pre></html>

      IMPLICIT NONE
      
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER ne,ng,NQNE(NEQM,NQEM),ASC_CVI,ICQS_SPATIAL(NQISVM,NQM)
      REAL*8 ASC,YQS(NIQSM,NQM),RCQS_SPATIAL(NQRSVM,NQM)
      CHARACTER ASC_NAME*8,ERROR*(*)
!     Local Variables
      
      CALL ENTERS('EVALASC',*9999)

     
      IF (ASC_NAME(1:3).EQ.'YQS') THEN
        ASC=YQS(ASC_CVI,NQNE(ne,ng))
      ELSEIF (ASC_NAME(1:3).EQ.'RCQ') THEN
        ASC=RCQS_SPATIAL(ASC_CVI,NQNE(ne,ng))
        WRITE(OP_STRING,'('' >>WARNING: Did you really mean to add'
     &       //' a value of a RCQS_SPATIAL field ???'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSEIF (ASC_NAME(1:3).EQ.'ICQ') THEN
        ASC=ICQS_SPATIAL(ASC_CVI,NQNE(ne,ng))
        WRITE(OP_STRING,'('' >>WARNING: Should not happen. Nothing'
     &       //' added !!!! '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      
      
      CALL EXITS('EVALASC')
      RETURN
 9999 CALL ERRORS('EVALASC',ERROR)
      CALL EXITS('EVALASC')
      RETURN 1
      END


