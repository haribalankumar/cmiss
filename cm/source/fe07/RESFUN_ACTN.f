      SUBROUTINE RESFUN_ACTN(ISIZE_TBH,LIST_RESID,
     '  MODE,NPLIST3,nr,NYNO,NYNR,
     '  LAPL,LAPLSQR,PHI,PHI_H,RESID,RESJAC,T_BH,YP,
     '  ERROR,*)

C#### Subroutine: RESFUN_ACTN
C###  Description:
C###    RESFUN_ACTN returns residuals for optimising the inverse
C###    problem.

C This has been modified to speed up calculations for an activation
C inverse solution

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER ISIZE_TBH(2),LIST_RESID,MODE,NPLIST3(0:NPM),nr,
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      REAL*8 LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  RESID(*),RESJAC(NREM,*),T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  YP(NYM,NIYM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER CALC_SAMPLE_FROM_TIME,
     '  ERR,icol,noopti,no_nynr,no_nynr1,
     '  nores,nts,NPHI_RES,ny1,S_START,S_END
      REAL*8 CALC_TIME_FROM_SAMPLE,
     '  SS_RESID,SS_SUM,SUM2,SUM3,TOT_REG,TOT_RES,u,WIN2

      CALL ENTERS('RESFUN_ACTN',*9999)

      IF(ACTN_IREGULARISE.LE.1) THEN !no regul
        NPHI_RES=NT_RES
      ELSEIF(ACTN_IREGULARISE.EQ.2) THEN !surface lapl
        NPHI_RES=NT_RES-1
C LKC 20-MAR-2002 Invalid option
C      ELSEIF(ACTN_IREGULARISE.EQ.3) THEN !unused
C        NPHI_RES=NT_RES-1
      ELSE
        ERROR='>>Update Code: Unknown Regularisation'
        GOTO 9999
      ENDIF
      CALL ASSERT(NPHI_RES.EQ.ISIZE_TBH(1),
     '  '>> Number of resids <> # Torso nodes',ERROR,*9999)


C The following will need changing when we have derivatives.
      IF(NTOPTI .NE. (NPLIST3(0)-NPJOIN(0,1)) ) THEN
        WRITE(OP_STRING,'(I5,'' Optim variables'')') NTOPTI
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(I5,'' Transfer nodes'')') NPLIST3(0)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ERROR='#Optim Variables <> First surface nodes-NPJOIN'
        GOTO 9999
      ENDIF

      TOT_RES=0.0d0
      WIN2=TRSF_ACTN_WAVE_WIDTH/2.0d0

      DO nores=1,NPHI_RES
!       these should correspond to the surface nodes
!       PHI(nytr,nts) is the measured signal at torso dof
!                     nytr (=nores if linear elements) and time nts
!       PHI_H(nytr,nts) is the estimated signal at torso dof
!                     nytr (=nores if linear elements) and time nts

!       Need to calculate squared difference and integrate over time
!       (or just use sum of squares of differences over time for now)
!       i.e. resid function = sum_over_time{phi-phi_h}^2

        IF(MODE.EQ.0.OR.MODE.EQ.2) THEN !Calculate Function
          SS_RESID=0.0d0
          DO nts=SOPTI_START,SOPTI_END
            SS_RESID=SS_RESID+(PHI(nores,nts)-PHI_H(nores,nts))**2
          ENDDO
          RESID(nores)=SS_RESID
          TOT_RES=TOT_RES+RESID(nores)
        ENDIF !MODE .EQ. 0 or 2

        IF(MODE.EQ.1.OR.MODE.EQ.2) THEN !Calculate Jacobian
          DO noopti=1,NPLIST3(0)
            ny1=NYNO(1,noopti,2,nr)
            icol=noopti

C*** Looping from the lower window to the end of sequence
C*** so that d(phi_h[nores,nts])/d(ui) can also
C*** be calculated
            SS_SUM=0.0d0
C              DO nts=MAX(1,INT(YP(HEART_ny,1)-WIN2)),
C     '            INT(YP(HEART_ny,1)+WIN2)
            S_START=MAX(SOPTI_START,
     '        CALC_SAMPLE_FROM_TIME(YP(ny1,1)-WIN2,ERR,ERROR)+1)
            S_END=MIN(
     '        CALC_SAMPLE_FROM_TIME(YP(ny1,1)+WIN2,ERR,ERROR),
     '        SOPTI_END)
            DO nts=S_START,S_END
              u=CALC_TIME_FROM_SAMPLE(nts)-YP(ny1,1)
              IF(u.GE.-WIN2.AND.u.LE.0.0d0) THEN
                SS_SUM=SS_SUM+(PHI(nores,nts)-PHI_H(nores,nts))*
     '            T_BH(nores,icol)*(u/WIN2+1.0d0)
              ELSE IF(u.GE.0.0d0.AND.u.LE.WIN2) THEN
                SS_SUM=SS_SUM-(PHI(nores,nts)-PHI_H(nores,nts))*
     '            T_BH(nores,icol)*(u/WIN2-1.0d0)
              ENDIF !End of nts choices.
            ENDDO !nts
            RESJAC(nores,noopti)=SS_SUM*2.0d0*TRSF_ACTN_WAVE_JUMP/WIN2
          ENDDO !nolist1

        ENDIF !mode=1 or 2 (calculate Jacobian)

      ENDDO !nores

C*** Add an additional constraint as nores=NT_RES
      TOT_REG=0.0d0
      IF(ACTN_IREGULARISE.EQ.1) THEN
        IF(LIST_RESID.GE.2) THEN
          WRITE(OP_STRING,'(/'' *** No additional constraints'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSEIF(ACTN_IREGULARISE.EQ.2) THEN

C LKC 15-JUL-2002 New check
        CALL ASSERT(CALL_CALC_LAPL,' >> Evaluate Laplacian first',
     '    ERROR,*9999)

        IF(LIST_RESID.GE.2) THEN
          WRITE(OP_STRING,'(/'' *** Doing surface laplacian'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

C        CALL SURFACELAPLACIAN(IBT,LIST_RESID,NBH,NENP,NPLIST3,NPNE,
C     '    nx,NXI,NYNP,RESID,RESJAC,TOT_REG,XP,YP,ERROR,*9999)

          nr=TRSF_NR_FIRST !heart surface
          TOT_REG=0.D0
          DO no_nynr=1,NYNR(0,0,1,nr) !same as heart nodes
            SUM2=0.D0
            SUM3=0.D0
            DO no_nynr1=1,NYNR(0,0,1,nr) !same as heart nodes
              ny1=NYNR(no_nynr1,0,1,nr)
              SUM2=SUM2+LAPL(no_nynr,no_nynr1)*YP(ny1,1)
              SUM3=SUM3+LAPLSQR(no_nynr,no_nynr1)*YP(ny1,1)
            ENDDO
            TOT_REG=TOT_REG+SUM2*SUM2 !sum of squares
            RESJAC(NT_RES,no_nynr)=SUM3*ACTN_REG_PARAM_LAPLACE*2
          ENDDO
          RESID(NT_RES)=TOT_REG*ACTN_REG_PARAM_LAPLACE
C          WRITE(*,*)
C          WRITE(*,*) '***DER', NT_RES, no_nynr-2,
C     '      RESJAC(NT_RES,no_nynr-2),RESJAC(NT_RES-1,no_nynr-2)
C          WRITE(*,*) '***DER', NT_RES, no_nynr-3,
C     '      RESJAC(NT_RES,no_nynr-3),RESJAC(NT_RES-1,no_nynr-3)
C          WRITE(*,*) '***REG', NT_RES, RESID(NT_RES), RESID(NT_RES-1)

      ENDIF ! ACTN_IREGULARISE

C*** Output the residual values (have already been regularised)
      IF(LIST_RESID.GT.2) THEN
        WRITE(OP_STRING,'(/'' YP(noopti)'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(5D14.6)') (YP(NYNO(1,noopti,2,nr),1),
     '    noopti=1,NPLIST3(0))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING,'(/'' RESID(nts)'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(5D14.6)') (RESID(nts),nts=1,NT_RES)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        DO nts=1,NT_RES
          WRITE(OP_STRING,'(/'' RESJAC('',I4,'',opt)'')') nts
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(5D14.6)') (RESJAC(nts,noopti),
     '      noopti=1,NPLIST3(0))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C*** Output some statistics
      IF(LIST_RESID.GE.1) THEN
        WRITE(OP_STRING(1),'('' '')')
        WRITE(OP_STRING(2),'('' Activation Residual Summary'')')
        WRITE(OP_STRING(3),'('' '')')

C*** Ouput Phi Residuals
        WRITE(OP_STRING(4),
     '    '('' Total Phi Residual              : '',F11.4)') TOT_RES
        WRITE(OP_STRING(5),
     '    '('' Avg. Phi  Residual              : '',F11.4)')
     '    TOT_RES/DBLE(NT_RES)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ACTN_IREGULARISE.GE.2) THEN
C*** Output Laplacian Residuals
          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),
     '      '('' Total Laplacian Residual        : '',F11.4)')
     '      TOT_REG
          WRITE(OP_STRING(3),
     '      '('' Avg. Laplacian Residual         : '',F11.4)')
     '      TOT_REG/DBLE(NTOPTI)

          WRITE(OP_STRING(4),
     '      '('' Reg. Avg. Laplacian Residual    : '',F11.4)')
     '      TOT_REG/DBLE(NTOPTI)*ACTN_REG_PARAM_LAPLACE
          WRITE(OP_STRING(5),
     '      '('' Reg. Tot. Laplacian Residual    : '',F11.4)')
     '      TOT_REG*ACTN_REG_PARAM_LAPLACE
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Output ratios
          WRITE(OP_STRING(1),'('' '')')
          IF(DABS(TOT_REG).GT.ZERO_TOL) THEN
            WRITE(OP_STRING(2),
     '        '('' Avg. Phi-Laplacian Ratio        : '',E11.4)')
     '        TOT_RES/DBLE(NT_RES)/TOT_REG
          ELSE
            WRITE(OP_STRING(2),
     '        '('' Avg. Phi-Laplacian Ratio        : No Laplacian'')')
          ENDIF

          IF(DABS(ACTN_REG_PARAM_LAPLACE*TOT_REG).GT.ZERO_TOL) THEN
            WRITE(OP_STRING(3),
     '        '('' Reg. Avg. Phi-Laplacian Ratio   : '',E11.4)')
     '        TOT_RES/DBLE(NT_RES)/
     '        TOT_REG*ACTN_REG_PARAM_LAPLACE
          ELSE
            WRITE(OP_STRING(3),
     '        '('' Regularised Phi-Laplacian Ratio : No Laplacian'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF !LIST_RESID

      CALL EXITS('RESFUN_ACTN')
      RETURN
 9999 CALL ERRORS('RESFUN_ACTN',ERROR)
      CALL EXITS('RESFUN_ACTN')
      RETURN 1
      END


