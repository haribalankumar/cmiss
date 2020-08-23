      SUBROUTINE LIMATR(ISC_GD,ISC_GKK,ISC_GM,ISC_GMM,ISIZE_MFI,
     '  ISIZE_PHI,ISIZE_TBH,ISR_GD,ISR_GKK,ISR_GM,ISR_GMM,
     '  LD_NP,NRLIST,
     '  NXLIST,GD,GKK,GM,GMM,GR,GRR,MFI,PHI,PHI_H,PHI_H_EXACT,T_BH,
     '  T_BH_INV,YP,ZCROSSING,STRING,ERROR,*)

C#### Subroutine: LIMATR
C###  Description:
C###    LIMATR lists matrices.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'matr00.cmn'
      INCLUDE 'matr00.inc'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER ISC_GD(NISC_GDM),ISC_GKK(NISC_GKKM,NXM),
     '  ISC_GM(NISC_GMM),ISC_GMM(NISC_GMMM),ISIZE_MFI(3,NSSM),
     '  ISIZE_PHI(2),ISIZE_TBH(2),ISR_GD(NISR_GDM),
     '  ISR_GKK(NISR_GKKM,NXM),ISR_GM(NISR_GMM),ISR_GMM(NISR_GMMM),
     '  LD_NP(NDM),NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 GD(NZ_GD_M),GKK(NZ_GKK_M,NXM),GM(NZ_GM_M),
     '  GMM(NZ_GMM_M),GR(NYROWM),GRR(NOM),MFI(NDM,NTSM,3,NSSM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  PHI_H_EXACT(NY_TRANSFER_M,NTSM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  T_BH_INV(NY_TRANSFER_M,NY_TRANSFER_M),YP(NYM,NIYM,NXM),
     '  ZCROSSING(NY_TRANSFER_M,NTSM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,matr,N3CO,nc,nomatr,
     '  no_nrlist,no_row,novect,nr,nss,nx,nxc,NYTT,row,ROWLIST(0:20)
      CHARACTER FILE*100,INOUTTYPE*9,NXTYPE*10
      LOGICAL ALL_REGIONS,ALLROWS,ABBREV,CBBREV,FOUND,OPFILE,SPARSITY,
     '  VALUES

!     Functions
      INTEGER IFROMC

      CALL ENTERS('LIMATR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list matrix<;FILENAME>
C###  Description:
C###    Lists matrices to the screen or to the file FILENAME.opmatr if
C###    the qualifier is present.
C###  Parameter:    <matrices NAMELIST[GK]>
C###    Specifies the list of matrices to output. Matrices available to
C###    be output are: GK, GQ, GD, GM, GR, GKK, GMM, GRR, YP, T_BH,
C###    T_BH_INV, PHI, ZCROSSING, PHI_H, YQ, YQS, PHI_H_EXACT, MFI,
C###    LD_NP.
C###  Parameter:    <nc (#s/all)[1]>
C###    Specify the nc (dependent variable derivative type) values of
C###    the dependent variable matrices. If the option 'all' is
C###    specified then nc values for that dependent variable matrix are
C###    used.
C###  Parameter:    <sparsity_pattern>
C###    Specify that the sparsity pattern of the matrix (if it is a
C###    sparse matrix), as opposed to the actual matrix values, is
C###    to be listed.
C###  Parameter:    <values>
C###    Specify that the actual values of the matrix, as opposed to
C###    the sparsity pattern, is to be listed.
C###  Parameter:    <rows (#s/all)[all]>
C###    Specify the rows of the matrix to list. If the 'all' option is
C###    specified then all the rows of the matrix are output.
C###  Parameter:    <nss (#)[1]>
C###    Specify the signal set for the MFI array
C###  Parameter:    <rows (#s/all)[all]>
C###    Specify the rows of the matrix to list. If the 'all' option is
C###    specified then all the rows of the matrix are output.

C###  Parameter:    <region (#s/all)[1]>
C###    Specify the list of region numbers which contain the variables
C###    to be output.
C###  Parameter:    <using (fit/optimise/solve)[solve]>
C###    Specify the type of solution problem to use.
C###  Parameter:    <class #[1]>
C###    Specify the class of solution problem to use.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<matrices NAMELIST[GK]>'
        OP_STRING(3)=BLANK(1:15)//'<nc (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<sparsity_pattern>'
        OP_STRING(5)=BLANK(1:15)//'<values>'
        OP_STRING(3)=BLANK(1:15)//'<nss (#)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<rows (#s/all)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<using (fit/optimise/solve)[solve]>'
        OP_STRING(9)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIMATR',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opmatr','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '    ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'MATRICES',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSSL(CO(N3CO+1),NUMMATRMX,NUMMATRNAMES,MATRNAMELIST,
     '      ERROR,*9999)
          DO nomatr=1,NUMMATRNAMES
C            MATRNAMELIST(nomatr)=CUPPER(MATRNAMELIST(nomatr))
            CALL CUPPER(MATRNAMELIST(nomatr),MATRNAMELIST(nomatr))
          ENDDO !nomatr
        ELSE
          NUMMATRNAMES=1
          MATRNAMELIST(1)='GK '
        ENDIF

        MATRLIST(0)=0
        VECTLIST(0)=0
        DO nomatr=1,NUMMATRNAMES
          CALL STRING_TRIM(MATRNAMELIST(nomatr),IBEG1,IEND1)
          matr=1
          FOUND=.FALSE.
          DO WHILE(.NOT.FOUND.AND.matr.LE.NUMMATRMX)
            CALL STRING_TRIM(MATRNAME(matr),IBEG2,IEND2)
            IF(MATRNAMELIST(nomatr)(IBEG1:IEND1).EQ.
     '        MATRNAME(matr)(IBEG2:IEND2)) THEN
              FOUND=.TRUE.
            ELSE
              matr=matr+1
            ENDIF
          ENDDO
          IF(FOUND) THEN
            IF(MATRTYPE(matr).EQ.0) THEN !is a vector
              VECTLIST(0)=VECTLIST(0)+1
              VECTLIST(VECTLIST(0))=matr
            ELSE IF(MATRTYPE(matr).EQ.1) THEN !is a matrix
              MATRLIST(0)=MATRLIST(0)+1
              MATRLIST(MATRLIST(0))=matr
            ENDIF
          ELSE
            WRITE(ERROR,'('' >>Matrix '',A,'' is unknown'')')
     '        MATRNAMELIST(nomatr)(IBEG1:IEND1)
            GOTO 9999
          ENDIF
        ENDDO !matr

        IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'SOLVE',2)) THEN
            CALL ASSERT(CALL_SOLV,'>>No solve defined',ERROR,*9999)
            NXTYPE='SOLVE'
          ELSE IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
            CALL ASSERT(CALL_FIT,'>>No fit defined',ERROR,*9999)
            NXTYPE='FIT'
          ELSE IF(ABBREV(CO(N3CO+1),'OPTIMISE',2)) THEN
            IF(KTYP8.EQ.6) THEN !data fitting by optimisation
              CALL ASSERT(CALL_FIT,'>>No fit defined',ERROR,*9999)
            ELSE
              CALL ASSERT(CALL_OPTI,'>>No optimisation defined',
     '          ERROR,*9999)
            ENDIF
            NXTYPE='OPTIMISE'
          ELSE
            NXTYPE='SOLVE'
            CALL ASSERT(CALL_SOLV,'>>No solve defined',ERROR,*9999)
          ENDIF
        ELSE
C LKC 17-APR-2002 do not have to have a solve defined for all matrices..
C          CALL ASSERT(CALL_SOLV,'>>No solve defined',ERROR,*9999)
          NXTYPE='SOLVE'
        ENDIF

        IF(NXTYPE(1:3).EQ.'FIT') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
        ELSE IF(NXTYPE(1:8).EQ.'OPTIMISE') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx.NE.0,
     '      '>>No nx defined for this optimisation class',
     '      ERROR,*9999)
        ELSE IF(NXTYPE(1:5).EQ.'SOLVE') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)

C LKC 17-APR-2002 do not need nx's for all matrices
C          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
C     '      ERROR,*9999)
          IF(nx.LE.0) THEN
            OP_STRING(1)='>>WARNING: No nx defined for this solve class'
            OP_STRING(2)='>>WARNING: Forcing nx to 1'
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            nx=1
          ENDIF
        ENDIF

        IF(CBBREV(CO,'SPARSITY_PATTERN',2,noco+1,NTCO,N3CO)) THEN
          SPARSITY=.TRUE.
          CALL ASSERT(USE_SPARSE.NE.0,
     '      '>>Set USE_SPARSE=1 to use sparse arrays',ERROR,*9999)
        ELSE
          SPARSITY=.FALSE.
        ENDIF

        IF(CBBREV(CO,'VALUES',2,noco+1,NTCO,N3CO)) THEN
          VALUES=.TRUE.
        ELSE
          VALUES=.FALSE.
        ENDIF

        IF(CBBREV(CO,'NSS',3,noco+1,NTCO,N3CO)) THEN
          nss=IFROMC(CO(N3CO+1))
        ELSE
          nss=1
        ENDIF
        CALL ASSERT(nss.LE.NSSM,
     '    '>> nss must be smaller than NSSM',ERROR,*9999)

        IF(SPARSITY.AND.VALUES) THEN
          INOUTTYPE='ALL_INOUT'
        ELSE IF(SPARSITY) THEN
          INOUTTYPE='SPARSITY'
        ELSE IF(VALUES) THEN
          INOUTTYPE='ARRAYS'
        ENDIF

        IF(CBBREV(CO,'ROWS',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ALL',2)) THEN
            ALLROWS=.TRUE.
          ELSE
            ALLROWS=.FALSE.
            CALL PARSIL(CO(N3CO+1),20,ROWLIST(0),ROWLIST(1),ERROR,
     '        *9999)
            DO nomatr=1,MATRLIST(0)
              matr=MATRLIST(nomatr)
              IF(matr.EQ.MATR_GK) THEN
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.NYT(1,1,nx),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_GQ) THEN
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.NYT(1,2,nx),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_GD) THEN
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.NYT(1,3,nx),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_GM) THEN
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.NYT(1,4,nx),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_GKK) THEN
                DO no_nrlist=1,NRLIST(0)
                  nr=NRLIST(no_nrlist)
                  DO no_row=1,ROWLIST(0)
                    row=ROWLIST(no_row)
                    CALL ASSERT(row.GT.0.AND.row.LE.NOT(1,1,nr,nx),
     '                '>>Invalid row number',ERROR,*9999)
                  ENDDO !no_row (row)
                ENDDO !nr
              ELSE IF(matr.EQ.MATR_GMM) THEN
                DO no_nrlist=1,NRLIST(0)
                  nr=NRLIST(no_nrlist)
                  DO no_row=1,ROWLIST(0)
                    row=ROWLIST(no_row)
                    CALL ASSERT(row.GT.0.AND.row.LE.NOT(1,4,nr,nx),
     '                '>>Invalid row number',ERROR,*9999)
                  ENDDO !no_row (row)
                ENDDO !nr
              ELSE IF(matr.EQ.MATR_YP) THEN
                NYTT=0
                DO nc=1,NCT(1,nx)
                  NYTT=NYTT+NYT(1,nc,nx)
                ENDDO !nc
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.NYTT,
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_T_BH) THEN
                CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '            ERROR,*9999)
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.ISIZE_TBH(1),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_T_BH_INV) THEN
                CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '            ERROR,*9999)
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.ISIZE_TBH(2),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_PHI) THEN
                CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '            ERROR,*9999)
                CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '            ERROR,*9999)
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.ISIZE_PHI(1),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_ZCROSSING) THEN
                CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '            ERROR,*9999)
                CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '            ERROR,*9999)
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.ISIZE_TBH(1),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_PHI_H) THEN
                CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '            ERROR,*9999)
                CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '            ERROR,*9999)
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.ISIZE_TBH(2),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_PHI_H_EXACT) THEN
                CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '            ERROR,*9999)
                CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '            ERROR,*9999)
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.ISIZE_TBH(2),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ENDIF
            ENDDO !nomatr (matr)
            DO novect=1,VECTLIST(0)
C              vect=VECTLIST(novect)
              IF(matr.EQ.MATR_GR) THEN
                DO no_row=1,ROWLIST(0)
                  row=ROWLIST(no_row)
                  CALL ASSERT(row.GT.0.AND.row.LE.NYT(1,1,nx),
     '              '>>Invalid row number',ERROR,*9999)
                ENDDO !no_row (row)
              ELSE IF(matr.EQ.MATR_GRR) THEN
                DO no_nrlist=1,NRLIST(0)
                  nr=NRLIST(no_nrlist)
                  DO no_row=1,ROWLIST(0)
                    row=ROWLIST(no_row)
                    CALL ASSERT(row.GT.0.AND.row.LE.NOT(1,1,nr,nx),
     '                '>>Invalid row number',ERROR,*9999)
                  ENDDO !no_row (row)
                ENDDO !nr
              ENDIF
            ENDDO !novect (vect)
          ENDIF
        ELSE
          ALLROWS=.TRUE.
        ENDIF
        CALL OPMATR(ISC_GD,%VAL(ISC_GK_PTR(nx)),ISC_GKK(1,nx),
     '    ISC_GM,ISC_GMM,%VAL(ISC_GQ_PTR(nx)),ISIZE_MFI(1,nss),
     '    ISIZE_PHI,ISIZE_TBH,ISR_GD,%VAL(ISR_GK_PTR(nx)),
     '    ISR_GKK(1,nx),ISR_GM,ISR_GMM,%VAL(ISR_GQ_PTR(nx)),LD_NP,
     &    NRLIST,nx,ROWLIST,GD,%VAL(GK_PTR(nx)),GKK(1,nx),GM,GMM,
     '    %VAL(GQ_PTR(nx)),GR,GRR,MFI(1,1,1,nss),
     '    PHI,PHI_H,PHI_H_EXACT,T_BH,T_BH_INV,ZCROSSING,
     '    YP(1,1,nx),INOUTTYPE,ALLROWS,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIMATR')
      RETURN
 9999 CALL ERRORS('LIMATR',ERROR)
      CALL EXITS('LIMATR')
      RETURN 1
      END


