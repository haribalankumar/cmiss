      SUBROUTINE RESFUN_ACTN_CC(PARAMTYPE,IBT,ISIZE_TBH,LIST_RESID,
     '  MODE,NBH,NENP,NPLIST3,NPNE,nr,nx,NXI,NYNO,NYNP,PHI,PHI_H,
     '  RESID,RESJAC,T_BH,WK1_INV,XP,YP,ERROR,*)

C#### Subroutine: RESFUN_ACTN_CC
C###  Description:
C###    RESFUN_ACTN_CC returns residuals for optimising the activation times
C###    when there is a Corelation Coefficent component to the objective
C###    funciton.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ISIZE_TBH(2),LIST_RESID,MODE,
     '  NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),NPLIST3(0:NPM),
     '  NPNE(NNM,NBFM,NEM),nr,nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  RESID(*),RESJAC(NREM,*),T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)

!     Local Variables
      INTEGER COUNT,icol,noopti,nores,nts,ny1
      REAL*8 TOT_LAPL,TOT_RES,u,WIN2
      REAL*8 PHI_MEAN,PHI_H_MEAN,PHI_NORM,PHI_H_NORM,SS_RESID,CC_RESID,
     '  PHI_PHI_H_SUM,SS_SUM,CC_SUM,D_PHI_H_MEAN

C not referenced MHT 09-11-00
C      INTEGER nb,nb1,ne,NE1,nh1,ni,nk1,nolist1,np,NP1,NP2,NP3,nv1,
C      REAL*8 B_CC_SUM,B_SS_SUM,CENT_DIFF2,delta1,delta2,SUM,SUM2,SUM3,
C             yp0,yp1,yp2,B_D_PHI_H_MEAN
C      LOGICAL FOUND
C      LOGICAL INLIST
C!     Functions
C      REAL*8 NODE_DIST,SQ_DIST


      CALL ENTERS('RESFUN_ACTN_CC',*9999)

      IF(PARAMTYPE(1:9).EQ.'POTENTIAL') THEN

        IF(MODE.EQ.1.OR.MODE.EQ.2) THEN !Calculate Jacobian

          CALL ASSERT(TOPTI_END-TOPTI_START.LE.NY_TRANSFER_M,
     '      '>>WK1_INV is to small for temporal info.'
     '      //' Increase NY_TRANSFER_M',
     '      ERROR,*9999)
          CALL ASSERT(NT_RES.LE.NY_TRANSFER_M,
     '      '>>WK1_INV is too small for resids. Increase NY_TRANSFER_M',
     '      ERROR,*9999)

          IF(ACTN_IREGULARISE.LE.1) THEN
            CALL ASSERT(NT_RES.EQ.ISIZE_TBH(1),
     '        '>> Number of resids <> # Torso nodes',ERROR,*9999)
          ELSEIF(ACTN_IREGULARISE.EQ.2) THEN
            CALL ASSERT(NT_RES-1.EQ.ISIZE_TBH(1),
     '        '>> Number of resids-1 <> # Torso nodes',ERROR,*9999)
          ELSE
            ERROR='Update Code: Unknown Regularisation'
            GOTO 9999
          ENDIF

C The following will need changing when we have derivatives.
          IF(NTOPTI .NE. (NPLIST3(0)-NPJOIN(0,1)) ) THEN
            WRITE(OP_STRING,'(I5,'' Optim variables'')') NTOPTI
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(I5,'' Transfer nodes'')') NPLIST3(0)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ERROR='#Optim Variables <> First surface nodes-NPJOIN'
            GOTO 9999
          ENDIF
        ENDIF !mode=1 or 2

        TOT_RES=0.0d0

        DO nores=1,NT_RES-ACTN_IREGULARISE+1
!         these should correspond to the surface nodes
!         PHI(nytr,nts) is the measured signal at torso dof
!                       nytr (=nores if linear elements) and time nts
!         PHI_H(nytr,nts) is the estimated signal at torso dof
!                       nytr (=nores if linear elements) and time nts

!         Need to calculate squared difference and integrate over time
!         (or just use sum of squares of differences over time for now)
!         i.e. resid function = sum_over_time{phi-phi_h}^2

!         Need to calculate the correlation coefficient between PHI and
!         PHI_H over time.

          PHI_MEAN=0.0d0
          PHI_H_MEAN=0.0d0
          count=0
          DO nts=1,INT(TOPTI_END-TOPTI_START)
            count=count+1
            PHI_MEAN=PHI_MEAN+PHI(nores,nts)
            PHI_H_MEAN=PHI_H_MEAN+PHI_H(nores,nts)
          ENDDO
          IF(COUNT.GT.0) THEN
            PHI_MEAN=PHI_MEAN/DBLE(count)
            PHI_H_MEAN=PHI_H_MEAN/DBLE(count)
          ENDIF

          PHI_NORM=0.0d0
          PHI_H_NORM=0.0d0
          DO nts=1,INT(TOPTI_END-TOPTI_START)
            PHI_NORM=PHI_NORM+(PHI(nores,nts)-PHI_MEAN)**2
            PHI_H_NORM=PHI_H_NORM+(PHI_H(nores,nts)-PHI_H_MEAN)**2
          ENDDO
          PHI_NORM=DSQRT(PHI_NORM)
          PHI_H_NORM=DSQRT(PHI_H_NORM)

          SS_RESID=0.0d0
          CC_RESID=0.0d0
          DO nts=1,INT(TOPTI_END-TOPTI_START)
            SS_RESID=SS_RESID+(PHI(nores,nts)-PHI_H(nores,nts))**2
            CC_RESID=CC_RESID+(PHI(nores,nts)-PHI_MEAN)*
     '        (PHI_H(nores,nts)-PHI_H_MEAN)
          ENDDO
          IF(DABS(PHI_NORM*PHI_H_NORM).GT.ZERO_TOL)
     '      CC_RESID=CC_RESID/(PHI_NORM*PHI_H_NORM)
          RESID(nores)=SS_OBJ_WEIGHT*SS_RESID+
     '      CC_OBJ_WEIGHT*DSQRT(1.0d0-CC_RESID)
          TOT_RES=TOT_RES+RESID(nores)

          IF(MODE.EQ.1.OR.MODE.EQ.2) THEN !Calculate Jacobian

            PHI_PHI_H_SUM=0.0d0
            DO nts=1,INT(TOPTI_END-TOPTI_START)
              PHI_PHI_H_SUM=PHI_PHI_H_SUM+(PHI(nores,nts)-PHI_MEAN)*
     '          (PHI_H(nores,nts)-PHI_H_MEAN)
            ENDDO

            DO noopti=1,NPLIST3(0)
              ny1=NYNO(1,noopti,2,nr)
              icol=noopti

              IF(YP(ny1,1).LT.0.0d0) THEN
                WRITE(OP_STRING,
     '            '('' Node  '',I5,'' has -ve activation time'
     '            //' '',F8.2)')
     '            NPLIST3(noopti),YP(ny1,1)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                GOTO 9999
              ENDIF
              IF(YP(ny1,1).GT.DBLE(NTSM)) THEN
                WRITE(OP_STRING,
     '            '('' Need to increase NTSM to at least'
     '            //' '',F8.2,'' for node  '',I5)')
     '            YP(ny1,1),noopti
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                GOTO 9999
              ENDIF

C*** Looping from the lower window to the end of sequence
C*** so that d(phi_h[nores,nts])/d(ui) can also
C*** be calculated
              WIN2=TRSF_ACTN_WAVE_WIDTH/2.0d0
              DO nts=MAX(1,INT(YP(ny1,1)-WIN2)),INT(YP(ny1,1)+WIN2)
                u=DBLE(nts)-YP(ny1,1)
                IF(u.GE.-WIN2.AND.u.LE.0.0d0) THEN
                  WK1_INV(1,nts)=WK1_INV(1,nts)-
     '              T_BH(nores,icol)*(u/WIN2+1.0d0)/WIN2
     '              *TRSF_ACTN_WAVE_JUMP
                ELSE IF(u.GE.0.0d0.AND.u.LE.WIN2) THEN
                  WK1_INV(1,nts)=WK1_INV(1,nts)+
     '              T_BH(nores,icol)*(u/WIN2-1.0d0)/WIN2
     '              *TRSF_ACTN_WAVE_JUMP
                ENDIF !End of nts choices.
              ENDDO !nts
 !                     WK1_INV(i,nts) now contains
 !                        i=1  d(phi_h[nores,nts])/d(taui)
              D_PHI_H_MEAN=0.0d0
C              B_D_PHI_H_MEAN=0.0d0
              DO nts=1,INT(TOPTI_END-TOPTI_START)
                D_PHI_H_MEAN=D_PHI_H_MEAN+WK1_INV(1,nts)
              ENDDO
              IF(count.GT.0) THEN
                D_PHI_H_MEAN=D_PHI_H_MEAN/DBLE(count)
              ENDIF
              SS_SUM=0.0d0
              CC_SUM=0.0d0
              DO nts=1,INT(TOPTI_END-TOPTI_START)
                SS_SUM=SS_SUM-2.0d0*(PHI(nores,nts)-
     '            PHI_H(nores,nts))*WK1_INV(1,nts)
                CC_SUM=CC_SUM+(PHI_H_NORM**2*PHI(nores,nts)-
     '            PHI_PHI_H_SUM*PHI_H(nores,nts))*
     '            (WK1_INV(1,nts)-D_PHI_H_MEAN)
              ENDDO
              CC_SUM=-CC_SUM/(2.0d0*DSQRT(1.0d0-CC_RESID))
              IF(DABS(PHI_NORM*PHI_H_NORM).GT.ZERO_TOL) THEN
                CC_SUM=CC_SUM/(PHI_NORM*PHI_H_NORM**3)
              ENDIF
              RESJAC(nores,noopti)=SS_OBJ_WEIGHT*SS_SUM+
     '          CC_OBJ_WEIGHT*CC_SUM
            ENDDO !nolist1

          ENDIF !mode=1 or 2 (calculate Jacobian)

        ENDDO !nores

C*** Add an additional constraint as nores=NT_RES
        IF(ACTN_IREGULARISE.EQ.1) THEN
          IF(LIST_RESID.GE.2) THEN
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),'('' *** No additional constraints'')')
            WRITE(OP_STRING(3),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE IF(ACTN_IREGULARISE.EQ.2) THEN
          IF(LIST_RESID.GE.2) THEN
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),'('' *** Doing surface laplacian'')')
            WRITE(OP_STRING(3),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          TOT_LAPL=0.0d0
          CALL SURFACELAPLACIAN(IBT,LIST_RESID,
     '      NBH,NENP,NPLIST3,NPNE,nx,NXI,NYNP,
     '      RESID,RESJAC,TOT_LAPL,XP,YP,
     '      ERROR,*9999)




C*** Output the residual values (have already been regularised)
          IF(LIST_RESID.GT.2) THEN
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),'('' RESID(nts)'')')
            WRITE(OP_STRING(3),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(5D14.6)') (RESID(nts),nts=1,NT_RES)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            DO nts=1,NT_RES
              WRITE(OP_STRING(1),'('' '')')
              WRITE(OP_STRING(2),'('' RESJAC('',I4,'',opt)'')') nts
              WRITE(OP_STRING(3),'('' '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(5D14.6)') (RESJAC(nts,noopti),
     '          noopti=1,NPLIST3(0))
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
     '        '('' Total Phi Residual              : '',F11.4)') TOT_RES
            WRITE(OP_STRING(5),
     '        '('' Avg.  Phi Residual              : '',F11.4)')
     '        TOT_RES/DBLE(NT_RES)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Output Laplacian Residuals
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),
     '        '('' Total Laplacian Residual        : '',F11.4)')
     '        TOT_LAPL
            WRITE(OP_STRING(3),
     '        '('' Avg. Laplacian Residual         : '',F11.4)')
     '        TOT_LAPL/DBLE(NT_RES)
            WRITE(OP_STRING(4),
     '        '('' Regularised Laplacian Residual  : '',F11.4)')
     '        TOT_LAPL/DBLE(NT_RES)*ACTN_REG_PARAM_LAPLACE
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Output ratios
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),
     '        '('' Phi-Laplacian Ratio             : '',F11.4)')
     '        TOT_RES/DBLE(NT_RES) / (TOT_LAPL/DBLE(NT_RES))

            IF(DABS(ACTN_REG_PARAM_LAPLACE*TOT_LAPL).GT.ZERO_TOL) THEN
              WRITE(OP_STRING(3),
     '          '('' Regularised Phi-Laplacian Ratio : '',F11.4)')
     '          TOT_RES/DBLE(NT_RES)/
     '          (TOT_LAPL*ACTN_REG_PARAM_LAPLACE/DBLE(NT_RES))
            ELSE
              WRITE(OP_STRING(3),
     '          '('' Regularised Phi-Laplacian Ratio : No Laplacian'')')
            ENDIF

            WRITE(OP_STRING(4),
     '        '('' Alpha                           : '',F11.4)')
     '        DSQRT( (TOT_LAPL/(DBLE(NT_RES)))**2 )
            WRITE(OP_STRING(5),
     '        '('' Regularised Alpha               : '',F11.4)')
     '        DSQRT( (TOT_LAPL*ACTN_REG_PARAM_LAPLACE/
     '        (DBLE(NT_RES)))**2 )
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF !LIST_RESID
        ENDIF ! ACTN_IREGULARISE
      ENDIF


      CALL EXITS('RESFUN_ACTN_CC')
      RETURN
 9999 CALL ERRORS('RESFUN_ACTN_CC',ERROR)
      CALL EXITS('RESFUN_ACTN_CC')
      RETURN 1
      END


