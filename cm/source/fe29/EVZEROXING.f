      SUBROUTINE EVZEROXING(ISIZE_PHI,ISIZE_TBH,NPLIST4,
     '  NXLIST,PHI,SIGMA_PHI,
     '  T_BH,U_PHI,VT_PHI,WK1_INV,WK2_INV,WK4_INV,ZCROSSING,STRING,
     '  ERROR,*)

C#### Subroutine: EVZEROXING
C###  Description:
C###    EVZEROXING evaluates the ZCROSSING array using the
C###    zero crossing algorithm of Fred Greensite and Geertjan Huiskamp

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ISIZE_PHI(2),ISIZE_TBH(2),NPLIST4(0:NPM),NXLIST(0:NXM)
      REAL*8 PHI(NY_TRANSFER_M,NTSM),SIGMA_PHI(NY_TRANSFER_M),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  U_PHI(NY_TRANSFER_M,NY_TRANSFER_M),VT_PHI(NTSM,NTSM),
     '  WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK2_INV(NY_TRANSFER_M,NY_TRANSFER_M),WK4_INV(NY_TRANSFER_M),
     '  ZCROSSING(NY_TRANSFER_M,NTSM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER FULLRANK,IBEG,IBEG1,IEND,IEND1,IFAIL,
     '  k,LWORK,music,no_svt,nts,nx,
     '  nxc,nytr,rank
      REAL*8 ONE,SCALE,SIGMA_PHI_MAX,SSQ,TMP,ZERO
      CHARACTER FILE*(MXCH)
      LOGICAL OPFILE
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)

      CALL ENTERS('EVZEROXING',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM evaluate zeroxing<;FILENAME>
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Creates the ZCROSSING(nytr,nts) matrix from the double layer
C###    transfer matrix stored in T_BH and signals stored in PHI.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVZEROXING',ERROR,*9999)
      ELSE
        CALL ASSERT(EVALUATE_TRANSFER,
     '    '>>Evaluate the transfer first',ERROR,*9999)
        CALL ASSERT(EVALUATE_PHI,
     '    '>>Evaluate PHI first',ERROR,*9999)
        CALL ASSERT(EVALUATE_PHI_SVD,
     '   '>>Evaulate PHI SVD first',ERROR,*9999)

        CALL ASSERT(NPLIST4(0).EQ.ISIZE_PHI(1),
     '    '>>Transfer matrix and signal matrix different sizes ?',
     '    ERROR,*9999)

        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opzxng','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        SIGMA_PHI_MAX=SIGMA_PHI(1)
        CALL ASSERT(SIGMA_PHI_MAX.GT.ZERO_TOL,
     '   '>> ZCROSSING matrix is zero ?',ERROR,*9999)
        CALL ASSERT(SVD_CUTOFF_RATIO.GE.ZERO,
     '    '>>Evalute PHI svd_cutoff # first',ERROR,*9999)

! Check array sizes for SVD
        LWORK=MAX(3*MIN(ISIZE_PHI(1),NTST)+MAX(ISIZE_PHI(1),NTST),
     '    5*MIN(ISIZE_PHI(1),NTST)-4)
        CALL ASSERT(LWORK.LE.NY_TRANSFER_M*NY_TRANSFER_M,
     '    '>>Work array too small for SVD of PHI',ERROR,*9999)

        DO nytr=1,ISIZE_TBH(2)
          DO nts=1,NTST
            ZCROSSING(nytr,nts)=0.0d0
          ENDDO !nts
        ENDDO !nytr

        ! Calculate the squared Euclidean norm for each column of T_BH
        DO nytr=1,ISIZE_TBH(2)
          SCALE=ZERO
          SSQ=ONE
          CALL DLASSQ(ISIZE_TBH(1),T_BH(1,nytr),1,SCALE,SSQ)
          WK4_INV(nytr)=(SCALE**2)*SSQ
          IF(WK4_INV(nytr).LT.ZERO_TOL) WK4_INV(nytr)=ONE
        ENDDO


CC For each time find the SVD of PHI(:,1..nts) and construct M+ function
CC and also find the SVD of PHI(:,nts..NTST) and construct M- function
CC ZCROSSING(nytr,nts) = M+(nytr,nts) - M-(nytr,nts)

        DO nts=1,NTST
          DO music=0,1
            k=(1-music)*nts+music*(NTST-nts+1)

            ! Construct M+/M- functions for the intervals [0,t] and
            !  [t,T] respectively.
            CALL DLACPY('F',ISIZE_PHI(1),k,PHI(1,MAX(1,music*nts)),
     '        NY_TRANSFER_M,WK1_INV,NY_TRANSFER_M)

            CALL DGESVD('A','N',ISIZE_PHI(1),k,WK1_INV,
     '        NY_TRANSFER_M,SIGMA_PHI,U_PHI,NY_TRANSFER_M,VT_PHI,NTSM,
     '        WK2_INV,NY_TRANSFER_M*NY_TRANSFER_M,IFAIL)
            IF(IFAIL.NE.0)THEN
              WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in SVD for M+/M-'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ERROR='>>IFAIL<>0 in DGESVD'
              GOTO 9999
            ENDIF

            FULLRANK=MIN(ISIZE_PHI(1),k)
            no_svt=0
            DO rank=1,FULLRANK
              IF(SIGMA_PHI(rank)/SIGMA_PHI_MAX.GT.SVD_CUTOFF_RATIO)
     '          THEN
                no_svt=no_svt+1
              ENDIF
            ENDDO

            CALL DGEMM('T','N',no_svt,ISIZE_TBH(2),ISIZE_TBH(1),
     '        ONE,U_PHI,NY_TRANSFER_M,T_BH,NY_TRANSFER_M,ZERO,
     '        WK1_INV,NY_TRANSFER_M)
            DO nytr=1,ISIZE_TBH(2)
              SCALE=ZERO
              SSQ=ONE
              CALL DLASSQ(no_svt,WK1_INV(1,nytr),1,SCALE,SSQ)
              TMP=(SCALE**2)*SSQ/WK4_INV(nytr)
              CALL ASSERT(DABS(TMP-ONE).GT.RDELTA,
     '            '>>MUSIC distance is 0 ?',ERROR,*9999)
              ZCROSSING(nytr,nts)=DBLE(music)*ZCROSSING(nytr,nts)+
     '          DBLE((-1)**music)*(ONE/DABS(ONE-TMP))
            ENDDO
          ENDDO
        ENDDO

CC Note: No critical times/points have been calculated yet - just the
CC  whole zero-crossing matrix.

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVZEROXING')
      RETURN
 9999 CALL ERRORS('EVZEROXING',ERROR)
      CALL EXITS('EVZEROXING')
      RETURN 1
      END

C---------------------------------------------------------------------
