      SUBROUTINE RESFUN_DIPOLE(ISIZE_MFI,ISIZE_PHI,LIST_RESID,
     '  TSTART,TEND,MFI,PHI,PHI_H,RESID,ERROR,*)

C#### Subroutine: RESFUN_DIPOLE
C###  Description:
C###    RESFUN_DIPOLE returns residuals of magnetic fields
C###    for optimising the inverse problem.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER ISIZE_MFI(3,NSSM),ISIZE_PHI(2),LIST_RESID
      REAL*8 MFI(NDM,NTSM,3,NSSM),PHI(NY_TRANSFER_M,NTSM),
     '  PHI_H(NY_TRANSFER_M,NTSM),RESID(*),TSTART,TEND
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ERR,nd,nj,nres_start,nss_mast,nt,SSTART,SEND,SAMPLES
      REAL*8 SUMSQR
      PARAMETER (nss_mast=1)

!     Functions
      INTEGER CALC_SAMPLE_FROM_TIME


      CALL ENTERS('RESFUN_DIPOLE',*9999)

      CALL ASSERT(KTYP27B.GT.0,
     '  '>> Residual field not defined -- define optimisation',
     '  ERROR,*9999)


C*** Calculate residuals which are the electrodes
      ERR=0
      SSTART=CALC_SAMPLE_FROM_TIME(TSTART,ERR,ERROR)
      IF(ERR.NE.0) GOTO 9999
      SEND=CALC_SAMPLE_FROM_TIME(TEND,ERR,ERROR)
      IF(ERR.NE.0) GOTO 9999

      IF(KTYP27B.EQ.1.OR.KTYP27B.EQ.3) THEN !difference of MFI

C        WRITE(*,*) '***'
C        WRITE(*,*) '*** Magnetic field differences'
C        WRITE(*,*) '***'

C LKC 11-DEC-2002 New check that the residual matrices are the same
        CALL ASSERT(ISIZE_MFI(1,1).EQ.ISIZE_MFI(1,2),
     '    '>> The number of electrodes for magnetic residuals differ',
     '    ERROR,*9999)

        SUMSQR=0.d0
        DO nd=1,ISIZE_MFI(1,nss_mast)
          DO nt=SSTART,SEND
            RESID(nd)=0.d0
            DO nj=1,3
              RESID(nd)=RESID(nd)+(MFI(nd,nt,nj,1)-MFI(nd,nt,nj,2))**2
            ENDDO
            SUMSQR=SUMSQR+RESID(nd)
          ENDDO
        ENDDO
      ENDIF



      IF(KTYP27B.EQ.2.OR.KTYP27B.EQ.3) THEN !difference of PHI

C        WRITE(*,*) '***'
C        WRITE(*,*) '*** Potential field differences'
C        WRITE(*,*) '***'

        IF(KTYP27B.EQ.2) THEN
          nres_start=0
        ELSEIF(KTYP27B.EQ.3) THEN
          nres_start=ISIZE_MFI(1,nss_mast)
        ENDIF

        SUMSQR=0.d0
        DO nd=1,ISIZE_PHI(1)
          DO nt=SSTART,SEND
            RESID(nd+nres_start)=(PHI(nd,nt)-PHI_H(nd,nt))**2
            SUMSQR=SUMSQR+RESID(nd)
          ENDDO
        ENDDO

      ENDIF

      IF(LIST_RESID.GT.0) THEN
        CALL ASSERT(ISIZE_MFI(1,nss_mast).GT.0,
     '    '>> No electrodes/residuals for nss 1',ERROR,*9999)

        SAMPLES=SEND-SSTART+1
        WRITE(OP_STRING(1),'('' Residual Ouput Summary'')')
        WRITE(OP_STRING(2),'(''   Start Time  : '',F12.2,
     '    ''     End Time   : '',F12.2)') TSTART,TEND
        WRITE(OP_STRING(3),'(''   Samples     : '',I12,
     '    ''     Frequency  : '',F12.2)') SAMPLES,TRSF_FREQUENCY

        WRITE(OP_STRING(4),'(''   Electrodes  : '',I12,
     '    ''     Components : '',I12)')
     '    ISIZE_MFI(1,nss_mast),
     '    ISIZE_MFI(3,nss_mast)
        WRITE(OP_STRING(5),'(''   RMS Error   : '',E12.5,
     '    ''     Sum Sqrs   : '',E12.5)')
     '    DSQRT(SUMSQR/ISIZE_MFI(1,nss_mast)/
     '    SAMPLES/ISIZE_MFI(3,nss_mast)),SUMSQR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('RESFUN_DIPOLE')
      RETURN
 9999 CALL ERRORS('RESFUN_DIPOLE',ERROR)
      CALL EXITS('RESFUN_DIPOLE')
      RETURN 1
      END


