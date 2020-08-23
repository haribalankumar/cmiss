      SUBROUTINE CELL_ASSIGN_SPATIAL(CELL_ICQS_VALUE,ICQS,
     '  ICQS_SPATIAL,IICQS_SPATIAL,
     '  IRCQS_SPATIAL,nq,CELL_RCQS_VALUE,RCQS,RCQS_SPATIAL,ERROR,*)

C#### Subroutine: CELL_ASSIGN_SPATIAL
C###  See-Also: IPCELL*
C###  See-Also: IPMAT3_CELL
C###  See-Also: IPMAT3_CELL1
C###  See-Also: IICQS_SPATIAL
C###  See-Also: IRCQS_SPATIAL
C###  See-Also: ICQS_SPATIAL
C###  See-Also: RCQS_SPATIAL
C###  Description:
C###  Populates the ICQS and RCQS arrays (which will generally be local
C###  arrays for model solution) with the correct values for grid point
C###  nq. The values set are determined by the grid point's cellular
C###  variant, which sets the default values for all parameters from the
C###  CELL_ICQS_VALUE and CELL_RCQS_VALUE arrays, and then any spatial
C###  variation defined for this variant.

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'

      !Parameter list
      INTEGER CELL_ICQS_VALUE(NQIM,NQVM),ICQS(NQIM),
     '  ICQS_SPATIAL(NQISVM,NQM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM),nq
      REAL*8 CELL_RCQS_VALUE(NQRM,NQVM),RCQS(NQRM),
     '  RCQS_SPATIAL(NQRSVM,NQM)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER nqsv,VARIANT

      CALL ENTERS('CELL_ASSIGN_SPATIAL',*9999)

      VARIANT = ICQS_SPATIAL(1,nq)
C *** First, copy over all parameters from the appropriate variant.
      DO nqsv=1,NQIT
        ICQS(nqsv)=CELL_ICQS_VALUE(nqsv,VARIANT)
      ENDDO
      DO nqsv=1,NQRT
        RCQS(nqsv)=CELL_RCQS_VALUE(nqsv,VARIANT)
      ENDDO
C *** Then overwrite the spatially varying parameters, this also ensures
C     that the appropriate variant number is placed into ICQS.
      DO nqsv=1,IICQS_SPATIAL(0,VARIANT)
        ICQS(IICQS_SPATIAL(nqsv,VARIANT))=ICQS_SPATIAL(nqsv,nq)
      ENDDO
      DO nqsv=1,IRCQS_SPATIAL(0,VARIANT)
        RCQS(IRCQS_SPATIAL(nqsv,VARIANT))=RCQS_SPATIAL(nqsv,nq)
      ENDDO

      CALL EXITS('CELL_ASSIGN_SPATIAL')
      RETURN
 9999 CALL ERRORS('CELL_ASSIGN_SPATIAL',ERROR)
      CALL EXITS('CELL_ASSIGN_SPATIAL')
      RETURN 1
      END


