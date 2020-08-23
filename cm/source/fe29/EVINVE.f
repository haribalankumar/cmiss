      SUBROUTINE EVINVE(IPIV,ISEG,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,
     '  ISPLOTXY,NXLIST,PHI,PHI_H,PHI_H_EXACT,REG_PARAMETER,SIGMA_PHI,
     '  SIGMA_T_BH,T_BH,T_BH_INV,U_PHI,U_T_BH,VT_PHI,VT_T_BH,
     '  WK1_INV,WK2_INV,WK3_INV,WK4_INV,CSEG,STRING,ERROR,*)

C#### Subroutine: EVINVE
C###  Description:
C###    EVINVE evaluates inverse of transfer matrix (with
C###    regularisation).
C**** If a traditional inverse procedure is used then the inverse
C**** matrix can be either calculated explicitly or the inverse
C**** solution can be found by solving the system of equations in
C**** the least squares sense.  If an explicit inverse transfer
C**** matrix is requested, this routine calculates this and
C**** stores it in T_BH_INV.  This matrix can then be written out
C**** and also used with the "apply inverse" command.
C**** If the inverse solution is found indirectly from the matrix
C**** equations then this routine does nothing (T_BH_INV will be all
C**** zero) and it is not until the "apply inverse" stage that an
C**** inverse solution is found (NB An explicit inverse matrix is
C**** never constructed in this case and the user will not be able to
C**** write out T_BH_INV).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER IPIV(NY_TRANSFER_M),ISEG(*),ISIZE_PHI(2),ISIZE_PHIH(2),
     '  ISIZE_TBH(2),ISPLOTXY(2),NXLIST(0:NXM)
      REAL*8 PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  PHI_H_EXACT(NY_TRANSFER_M,NTSM),REG_PARAMETER(0:NTSM),
     '  SIGMA_PHI(NY_TRANSFER_M),SIGMA_T_BH(NY_TRANSFER_M),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  T_BH_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  U_PHI(NY_TRANSFER_M,NY_TRANSFER_M),
     '  U_T_BH(NY_TRANSFER_M,NY_TRANSFER_M),VT_PHI(NTSM,NTSM),
     '  VT_T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK2_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK3_INV(NY_TRANSFER_M,NY_TRANSFER_M),WK4_INV(NY_TRANSFER_M)
      CHARACTER CSEG(*)*(*),STRING*(MXCH),ERROR*(*)

!     Local Variables
      INTEGER IBEG,icol,ICOL_T_BH_INV,ICUTOFF,IEND,IFAIL,ip,irow,
     '  IROW_T_BH_INV,isize,LWORK,N3CO,nx,nxc,nytr,nytr2
      INTEGER*4 WORK_PHI_PTR
      REAL*8 CONDITION_NUMBER,P(1),Q(1),SUM,SUM2
      CHARACTER FILE*(MXCH)
      LOGICAL CBBREV,CUTOFF,OPFILE,SVD

      CALL ENTERS('EVINVE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate inverse<;FILENAME>
C###  Parameter:        <svd>
C###    Calculate singluar values and condition number of inverse
C###    matrix
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Evaluates the inverse solution. i.e. relates body surface
C###    potentials to epicardial potentials. If FILENAME is included
C###    the inverse matrix (T_BH_INV) is written to FILENAME.ipmatr
C###    (ASCII) or FILENAME.binmat (binary) and the condition number
C###    and the singular values of T_BH_INV are written to
C###    FILENAME.opinve, if no FILENAME given only the condition number
C###    and the singular values are written to the screen.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<svd>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVINVE',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opinve','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '    ERROR,*9999)
          IF(ICALC_TRANSFER.EQ.2) THEN
            CALL ASSERT(EVALUATE_PHI,'>>Evaluate PHI first',ERROR,*9999)
            CALL ASSERT(ISIZE_TBH(1).EQ.ISIZE_PHI(1),
     '        '>>Transfer matrix and signal matrix different sizes ?',
     '        ERROR,*9999)
            IF(IREGULARISE.EQ.3) THEN ! Picard Criterion
              CALL ASSERT(USE_GRAPHICS.EQ.1,'>>Set USE_GRAPHICS=1',
     '          ERROR,*9999)
            ELSEIF(IREGULARISE.EQ.5) THEN ! Optimal
              CALL ASSERT(EVALUATE_PHI_H_EXACT,
     '          '>>Evaluate PHI_H_EXACT first',ERROR,*9999)
            ENDIF
          ENDIF
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'SVD',2,noco+1,NTCO,N3CO)) THEN
          SVD=.TRUE.
        ELSE
          SVD=.FALSE.
        ENDIF

        IF(ICALC_TRANSFER.EQ.1) THEN !not Greensite's inverse method

C***      Initialise work arrays and T_BH_INV
          DO irow=1,NY_TRANSFER_M
            WK4_INV(irow)=0.0d0
            DO icol=1,NY_TRANSFER_M
              WK1_INV(irow,icol)=0.0d0
              WK2_INV(irow,icol)=0.0d0
              WK3_INV(irow,icol)=0.0d0
              T_BH_INV(irow,icol)=0.0d0
            ENDDO !icol
          ENDDO !irow

          IF(ICALC_TRANSFER.EQ.1) THEN !T_bh^t*T_bh needed in all cases
C***        Work out T_bh^t
            DO irow=1,ISIZE_TBH(1)
              DO icol=1,ISIZE_TBH(2)
                WK2_INV(icol,irow)=T_BH(irow,icol) !transpose of T_bh
              ENDDO
            ENDDO
C***        Work out T_bh^t*T_bh
            DO irow=1,ISIZE_TBH(2)
              DO icol=1,ISIZE_TBH(2)
                SUM=0.0d0
                DO isize=1,ISIZE_TBH(1)
                  SUM=SUM+WK2_INV(irow,isize)*T_BH(isize,icol)
                ENDDO
                WK3_INV(irow,icol)=SUM !T_BH^t*T_BH
              ENDDO
            ENDDO
          ENDIF !ICALC_TRANSFER=1

          IF(IREGULARISE.EQ.1) THEN
            IF(ICALC_TRANSFER.EQ.1) THEN !Explicit inv of T_bh^t*T_bh
              IFAIL=1

              CALL DGETRF(ISIZE_TBH(2),ISIZE_TBH(2),WK3_INV,
     '          NY_TRANSFER_M,IPIV,IFAIL)

              IF(IFAIL.NE.0)THEN
                WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in inverse calc'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ERROR='>>IFAIL<>0 in DGETRF'
                GOTO 9999
              ENDIF

              CALL DGETRI(ISIZE_TBH(2),WK3_INV,NY_TRANSFER_M,IPIV,
     '          WK4_INV,NY_TRANSFER_M,IFAIL)

              IF(IFAIL.NE.0)THEN
                WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in inverse calc'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ERROR='>>IFAIL<>0 in DGETRI'
                GOTO 9999
              ENDIF
              DO irow=1,ISIZE_TBH(2)
                DO icol=1,ISIZE_TBH(2)
                  WK1_INV(irow,icol)=WK3_INV(irow,icol)
                ENDDO
              ENDDO

!             Form (T_bh^t*T_bh)^(-1)*(T_bh^t)
              DO irow=1,ISIZE_TBH(2)
                DO icol=1,ISIZE_TBH(1)
                  SUM=0.0d0
                  DO isize=1,ISIZE_TBH(2)
                    SUM=SUM+WK1_INV(irow,isize)*WK2_INV(isize,icol)
                  ENDDO
                  T_BH_INV(irow,icol)=SUM !(T_bh^t*T_bh)^(-1)*(T_bh^t)
                ENDDO
              ENDDO
              IROW_T_BH_INV=ISIZE_TBH(2)
              ICOL_T_BH_INV=ISIZE_TBH(1)
            ENDIF !ICALC_TRANSFER

          ELSE IF(IREGULARISE.EQ.2) THEN
            IF(ICALC_TRANSFER.EQ.1) THEN
C***        Find SVD of T_bh^t*T_bh and remove singular values smaller
C***        than RMIN_SVD. Then calculate inverse of this new matrix.

C***         Check array sizes for SVD of T_BH
              LWORK=MAX(3*MIN(ISIZE_TBH(1),ISIZE_TBH(2))+
     '          MAX(ISIZE_TBH(1),ISIZE_TBH(2)),
     '          5*MIN(ISIZE_TBH(1),ISIZE_TBH(2))-4)
              CALL ASSERT(LWORK.LE.NY_TRANSFER_M*NY_TRANSFER_M,
     '          '>>Work array too small for SVD of T_BH',ERROR,*9999)
              DO irow=1,ISIZE_TBH(1)
                DO icol=1,ISIZE_TBH(2)
                  WK1_INV(irow,icol)=T_BH(irow,icol)
                ENDDO
              ENDDO
              IFAIL=1
              CALL DGESVD('A','A',ISIZE_TBH(1),ISIZE_TBH(2),WK1_INV,
     '          NY_TRANSFER_M,SIGMA_T_BH,U_T_BH,NY_TRANSFER_M,VT_T_BH,
     '          NY_TRANSFER_M,WK3_INV,NY_TRANSFER_M*NY_TRANSFER_M,
     '          IFAIL)
              IF(IFAIL.NE.0)THEN
                WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in SVD for PHI'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ERROR='>>IFAIL<>0 in DGESVD'
                GOTO 9999
              ENDIF
C***          Find cutoff of each body surface equation.
              nytr=0
              CUTOFF=.FALSE.
              CALL ASSERT(SIGMA_T_BH(1).GT.0,'>>All svs are zero?',
     '          ERROR,*9999)
              DO WHILE(.NOT.CUTOFF.AND.nytr.LT.ISIZE_TBH(1))
                IF(SIGMA_T_BH(nytr+1)/SIGMA_T_BH(1).GT.RMIN_SVD) THEN
                  nytr=nytr+1
                ELSE
                  CUTOFF=.TRUE.
                ENDIF
              ENDDO
              ICUTOFF=nytr
              DO nytr=1,ISIZE_TBH(2) !number of heart points (x)
                DO nytr2=1,ISIZE_TBH(1) !number of body points (y)
                  SUM2=0.0d0
                  DO ip=1,ICUTOFF !j
!                   Form TSVD inverse of T_BH -(x,y)th entry
                    SUM2=SUM2+VT_T_BH(ip,nytr)*U_T_BH(nytr2,ip)/
     '                SIGMA_T_BH(ip)
                  ENDDO !ip (j)
                  T_BH_INV(nytr,nytr2)=SUM2
                ENDDO !nytr2 (y)
              ENDDO !nytr (x)
            ENDIF
          ELSE IF(IREGULARISE.GE.4.AND.IREGULARISE.LE.6) THEN !Tikhonov
            IF(IREGULARISE.EQ.4) THEN !zero-order Tikhonov
              IF(ICALC_TRANSFER.EQ.1) THEN !Explicit inverse
                DO irow=1,ISIZE_TBH(2)
                  WK3_INV(irow,irow)=WK3_INV(irow,irow)+TIKH_VALUE
                ENDDO
                IFAIL=1

                CALL DGETRF(ISIZE_TBH(2),ISIZE_TBH(2),WK3_INV,
     '            NY_TRANSFER_M,IPIV,IFAIL)

                IF(IFAIL.NE.0)THEN
                  WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in inverse calc'
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR='>>IFAIL<>0 in DGETRF'
                  GOTO 9999
                ENDIF

                CALL DGETRI(ISIZE_TBH(2),WK3_INV,NY_TRANSFER_M,IPIV,
     '            WK4_INV,NY_TRANSFER_M,IFAIL)

                IF(IFAIL.NE.0)THEN
                  WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in inverse calc'
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR='>>IFAIL<>0 in DGETRI'
                  GOTO 9999
                ENDIF
                DO irow=1,ISIZE_TBH(2)
                  DO icol=1,ISIZE_TBH(2)
                    WK1_INV(irow,icol)=WK3_INV(irow,icol)
                  ENDDO
                ENDDO

!               Form (T_bh^t*T_bh+TIKH_VALUE)^(-1)*(T_bh^t)
                DO irow=1,ISIZE_TBH(2)
                  DO icol=1,ISIZE_TBH(1)
                    SUM=0.0d0
                    DO isize=1,ISIZE_TBH(2)
                      SUM=SUM+WK1_INV(irow,isize)*WK2_INV(isize,icol)
                    ENDDO
                    T_BH_INV(irow,icol)=SUM !(T_bh^t*T_bh)^(-1)*(T_bh^t)
                  ENDDO
                ENDDO
                IROW_T_BH_INV=ISIZE_TBH(2)
                ICOL_T_BH_INV=ISIZE_TBH(1)
              ENDIF
            ENDIF
          ELSE
            ERROR='>>That option not implemented yet'
            GOTO 9999
          ENDIF

          IF(SVD.AND.ICALC_TRANSFER.EQ.1) THEN
C           Find singular values and cond # of T_BH_INV (can only be
C           done when T_BH_INV has been explicitly evaluated.
            IFAIL=1
            DO irow=1,IROW_T_BH_INV
              DO icol=1,ICOL_T_BH_INV
                WK1_INV(irow,icol)=T_BH_INV(irow,icol)
              ENDDO
            ENDDO
            ISIZE=MIN(IROW_T_BH_INV,ICOL_T_BH_INV)

            CALL DGESVD('N','N',IROW_T_BH_INV,ICOL_T_BH_INV,WK1_INV,
     '        NY_TRANSFER_M,WK4_INV,Q,1,P,1,WK3_INV,
     '        NY_TRANSFER_M*NY_TRANSFER_M,IFAIL)

            IF(IFAIL.NE.0)THEN
              WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in SVD of T_BH_INV'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ERROR='>>IFAIL<>0 in DGESVD'
              GOTO 9999
            ENDIF
            IF(WK4_INV(ISIZE).GT.RDELTA)THEN
              CONDITION_NUMBER=WK4_INV(1)/WK4_INV(ISIZE)
            ELSE
              CONDITION_NUMBER=RMAX !default value
            ENDIF
            IF(OPFILE)THEN
              WRITE(OP_STRING,'(/'' Condition num of T_bh_inv :'',
     '          D12.4)') CONDITION_NUMBER
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Singular values of T_bh_inv'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(13X,5D12.4,'
     '          //'/:(13X,5D12.4))')(WK4_INV(icol),icol=1,ISIZE)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(OP_STRING,'(/'' Condition # of T_bh_inv :'',D12.4)')
     '        CONDITION_NUMBER
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Singular values of T_bh_inv'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5D12.4,'
     '        //'/:(13X,5D12.4))')(WK4_INV(icol),icol=1,ISIZE)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !svd

        ELSE IF(ICALC_TRANSFER.EQ.2) THEN ! Use this for new work

!         Dynamic memory allocation
          WORK_PHI_PTR=0
          CALL ALLOCATE_MEMORY(NY_TRANSFER_M*NTSM,0,DPTYPE,
     '      WORK_PHI_PTR,MEM_INIT,ERROR,*9999)

          CALL EVINVE_DYNAM(ISEG,ISIZE_TBH,ISPLOTXY,PHI,PHI_H,
     '      PHI_H_EXACT,REG_PARAMETER,SIGMA_PHI,SIGMA_T_BH,T_BH,U_PHI,
     '      U_T_BH,VT_PHI,VT_T_BH,WK1_INV,%VAL(WORK_PHI_PTR),WK2_INV,
     '      CSEG,ERROR,*9999)

          CALL FREE_MEMORY(WORK_PHI_PTR,ERROR,*9999)
C GBS  9-Aug-2000 - PHI_H is set up for heart nodes
C GBS 25-FEB-2002 moved here
          PHI_H_IS_HEART=.TRUE.
          ISIZE_PHIH(1)=ISIZE_TBH(2)
          ISIZE_PHIH(2)=ISIZE_PHI(2)
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
        EVALUATE_INVERSE=.TRUE.
      ENDIF

      CALL EXITS('EVINVE')
      RETURN
 9999 CALL ERRORS('EVINVE',ERROR)
      CALL EXITS('EVINVE')
      RETURN 1
      END

