      SUBROUTINE EVGRID(NAQ,NLQ,NWQ,NXQ,
     '  GCHQ,GUQ,PROPQ,YQ,STRING,ERROR,*)

C#### Subroutine: EVGRID
C###  Description:
C###    EVGRID evaluates grid residuals.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NAQ(NQM,NAM),NLQ(NQM),NWQ(8,0:NQM,NAM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 GCHQ(3,NQM),GUQ(3,3,NQM),PROPQ(3,3,4,2,NQM),
     '  YQ(NYQM,NIQM,NAM)
      CHARACTER STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,na,NITB,nq
      REAL*8 DT,RESID
      LOGICAL CONTINUE

      CALL ENTERS('EVGRID',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate grid<;FILENAME>
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <output LEVEL#[0]>
C###  Description:
C###
        OP_STRING(1)=STRING(1:IEND)
!        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
!        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
!        OP_STRING(3)=BLANK(1:15)//'<output LEVEL#[0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVAERO',ERROR,*9999)
      ELSE
        NITB=2  !temporary
        DT=0.d0 !temporary
        CONTINUE=.TRUE.

        DO WHILE(CONTINUE) !until all residuals are below threshold
          DO na=NMGT-1,1,-1
            CALL MG_INTERPOL(na,1,NAQ,NLQ,NXQ,YQ,.TRUE.,ERROR,*9999)
            CALL MG_RELAX( 1,na,1,NITB,NAQ(1,na),NLQ,NWQ(1,0,na),
     '        NXQ(-NIM,0,0,na),.TRUE.,.TRUE.,'STATIC',
     '        DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*9999)
            CALL MG_RESIDUAL(na,2,NITB,NAQ(1,na),NLQ,NWQ(1,0,na),
     '        NXQ(-NIM,0,0,na),.TRUE.,.TRUE.,'STATIC',
     '        DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*9999)
          ENDDO !na

          CONTINUE=.FALSE.
          DO nq=1,NQT !Loop over grid pts to check for convergence
            IF(NLQ(nq).EQ.1) THEN !nq active in grid 1
              RESID=DABS(YQ(nq,2,1))
              IF(RESID.GT.1.d-3) THEN !exceeds threshold
                WRITE(OP_STRING,'('' Residual at nq='',I8,'
     '            //''' is '',E12.3)') nq,RESID
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                CONTINUE=.TRUE.
                GOTO 100
              ENDIF !residual>tolerance
            ENDIF !active pt
          ENDDO !nq

 100      CONTINUE
        ENDDO !while CONTINUE
      ENDIF !?

      CALL EXITS('EVGRID')
      RETURN
 9999 CALL ERRORS('EVGRID',ERROR)
      CALL EXITS('EVGRID')
      RETURN 1
      END

