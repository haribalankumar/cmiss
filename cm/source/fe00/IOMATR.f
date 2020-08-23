      SUBROUTINE IOMATR(COMAND,ISC_GK,ISC_GKK,ISC_GQ,
     '  ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,
     '  ISR_GK,ISR_GKK,ISR_GQ,IUNIT,
     '  LD_NP,MAP_ART_VEIN,nr,nx,GK,GKK,GQ,GR,GRR,MFI,
     '  PHI,PHI_H,PHI_H_EXACT,T_BH,T_BH_INV,YP,YQ,YQS,ZCROSSING,
     '  FILEFORMAT,FILE_NAME,INOUTTYPE,ERROR,*)
      
C#### Subroutine: IOMATR
C###  Description:
C###    IOMATR reads and writes .ipmatr files. COMAND determines
C###    whether data is to be read from ('READ') or written to
C###    ('WRITE') IUNIT. FILEFORMAT indicates whether or not the file
C###    is an 'ASCII' or a 'BINARY' file. FILE_NAME specifies the
C###    file name of the matrix file e.g. FILE_NAME.ipmatr or
C###    FILE_NAME.binmat

C###  Comment: BINARY MATRIX FILE FORMAT
C###  Description:
C###    <HTML>
C###    For matrix tags each tag identifies a separate matrix in the
C###    file. Each tag identifer is an integer to identify the cmiss
C###    matrix. Current matrix/vector indentification codes are defined
C###    as parameters of the form MATR_x (where x is the matrix name)
C###    in the file matr00.cmn. The are currently defined as :
C###    <BR>
C###    1 (GK), 2 (GQ), 3 (GD), 4 (GM), 5 (GR), 6 (GKK), 7 (GMM),
C###    8 (GRR), 9 (YP), 10  (T_BH), 11  (T_BH_INV),
C###    12 (PHI), 13 (ZCROSSING),
C###    14 (PHI_H), 15 (YQ), 16 (YQS), 17 (PHI_H_EXACT), 18 (MFI),
C####   19 (LD_NP) 
C###    <BR>
C###    In addition the tag header is the name of the matrix.
C###    The format of each tag is as follows:
C###    <UL>
C###    <LI>An integer for the number of indicies
C###    <LI>The number of indices x integers to specify the index sizes.
C###    <LI>An integer to specify the size of the matrix.
C###    <LI>An integer to specify the number of matrix values stored
C###        in the file.
C###    <LI>An integer to specify whether or not this matrix has (=1) or
C###        has not (=0) a 'time' field associated with it.
C###    <LI>A double precision number indication the 'time' if the
C###        matrix does have a time field.
C###    <LI>An integer to specify the sparsity pattern of the matrix;
C###        0 if the matrix is not sparse, 1 if the matrix has
C###        compressed row sparsity, 2 or 4 if the matrix has row-column
C###        sparsity and 3 if the matrix has compressed column sparsity.
C###    <LI>If the matrix is a sparse matrix there will be a sequence
C###        of integers indicating the sparsity. If the sparsity is
C###        compressed row the sequence will be:
C###        <UL>
C###        <LI>(number of rows in the matrix + 1) integers indicating
C###            the row offsets for the matrix.
C###        <LI>For each row (row_offset(row+1)-row_offset(row) - 1)
C###            integers indicating the column numbers.
C###        </UL>
C###        If the sparsity is compressed column the sequence will be:
C###        <UL>
C###        <LI>(number of columns in the matrix + 1) integers
C###            indicating the column offsets for the matrix.
C###        <LI>For each column (column_offset(column+1)-
C###            column_offset(column) - 1) integers indicating the
C###            row numbers.
C###        </UL>
C###        If the sparsity pattern is row-column the sequence will be:
C###        <UL>
C###        <LI>(number of values) integers indicating the row
C###            numbers for the matrix.
C###        <LI>(number of values) integers indicating the column
C###            numbers for the matrix.
C###        </UL>
C###    <LI>A sequence of double precision numbers indicating the matrix
C###        values. If the matrix is not a sparse matrix the numbers
C###        will be a list of entries stored by columns. If the matrix
C###        has compressed row sparsity there will be (row_offset(row+1)
C###        - row_offset(row) - 1) numbers for each row. If the matrix
C###        has row-column sparsity there will be (number of non-zeros)
C###        numbers. If the matrix has compressed column sparsity there
C###        will be (column_offset(column+1) - column_offset(column)
C###        - 1) numbers for each column.
C###    </UL>
C###    </HTML>
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: MAP_ART_VEIN(nd,nr)
C###  Type: INTEGER
C###  Description:
C###    MAP_ART_VEIN(nd,nr) stores the terminal node (np) associated with each
C###    data point nd in region nr. This is used to map terminal arterials to
C###    terminal venules which have been generated using the bifurcating
C###    distributive algorithm.

C DPN 23/10/97 - add a file name parameter


      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'binf00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'matr00.cmn'
      INCLUDE 'matr00.inc'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISIZE_MFI(3),ISIZE_PHI(2),ISIZE_PHIH(2),
     '  ISIZE_TBH(2),ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  ISR_GQ(NISR_GQM),IUNIT,LD_NP(NDM),MAP_ART_VEIN(0:NDM,NRM),nr,nx
      REAL*8 GK(NZ_GK_M),GKK(NZ_GKK_M),GQ(NZ_GQ_M),
     '  GR(NYROWM),GRR(NOM),MFI(NDM,NTSM,3),
     '  PHI(NY_TRANSFER_M,NTSM),
     '  PHI_H(NY_TRANSFER_M,NTSM),PHI_H_EXACT(NY_TRANSFER_M,NTSM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  T_BH_INV(NY_TRANSFER_M,NY_TRANSFER_M),YP(NYM,NIYM),
     '  YQ(NYQM,NIQM),YQS(NIQSM,NQM),
     '  ZCROSSING(NY_TRANSFER_M,NTSM)
      CHARACTER COMAND*(*),ERROR*(*),FILEFORMAT*(*),FILE_NAME*(*),
     '  INOUTTYPE*(*)

!     Local Variables
      INTEGER IBEG1,IEND1,IBEG2,IEND2,IOTYPEOLD,ISC_DUMMY(1),
     '  ISR_DUMMY(1),FILEID,FILETYPE,matr,
     '  M_DUMMY,N_DUMMY,nc,nj,nomatr,novect,
     '  NUMMATRICES,NUMVECTORS,NYTT,NZ_DUMMY,NZ_MAX,
     '  ROWLIST_DUMMY(1),SPARSITY_TYPE,TAG,vect,VERSION(3)
      CHARACTER FMT*500,
     '  LINE*132,MATRICES(NUMMATRMX)*11,VECTORS(NUMMATRMX)*11
      LOGICAL FOUND

      CALL ENTERS('IOMATR',*9999)

      IOTYPEOLD=IOTYPE
      ROWLIST_DUMMY(1)=0


C LKC 17-JAN-2000 NZ_DUMMY should really be initialised - not
C  even sure if it should be used at all!
      NZ_DUMMY=0

      IF(COMAND(1:4).EQ.'READ') THEN

        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          VERSION(1)=2
          IOTYPE=2
          CALL OPEN_SEQ_FILE(VERSION(1),IUNIT,FILE_NAME,
     '      'matr','OLD',.FALSE.,ERROR,*9999)
C         Read in header information


C LKC 26-JUN-2000 This is invalid FORTRAN - but now works
C
C          READ(IUNIT,'('' Number of matrices: '',I1)') NUMMATRICES
          FMT='('' Number of matrices: '',I1)'
          READ(IUNIT,FMT) NUMMATRICES

          IF(NUMMATRICES.GT.0) THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C            READ(IUNIT,'('' Matrices are:'',7(1X,A11))')
C     '        (MATRICES(matr),matr=1,NUMMATRICES)
            FMT='('' Matrices are:'',7(1X,A11))'
            READ(IUNIT,FMT) (MATRICES(matr),matr=1,NUMMATRICES)
          ENDIF

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C          READ(IUNIT,'(/'' Number of vectors : '',I1)') NUMVECTORS
C          IF(NUMVECTORS.GT.0) THEN
C           READ(IUNIT,'('' Vectors are :'',7(1X,A11))') (VECTORS(vect),
C     '        vect=1,NUMVECTORS)
C          ENDIF
          FMT='(/'' Number of vectors : '',I1)'
          READ(IUNIT,FMT) NUMVECTORS
          IF(NUMVECTORS.GT.0) THEN
            FMT='('' Vectors are :'',7(1X,A11))'
            READ(IUNIT,FMT) (VECTORS(vect),vect=1,NUMVECTORS)
          ENDIF

C LKC 28-NOV-2002 Make sure there is something to read in
          CALL ASSERT(NUMVECTORS+NUMMATRICES.GE.1,
     '      '>> No vectors or matrices to read in',ERROR,*9999)

          IF(NUMMATRICES.GT.0) THEN
C Skip over header AJP 25/9/97
            READ(IUNIT,'(A)') LINE
            READ(IUNIT,'(A)') LINE
          ENDIF
          DO nomatr=1,NUMMATRICES
            READ(IUNIT,'(A)') LINE
            READ(IUNIT,'(A)') LINE
            CALL STRING_TRIM(MATRICES(nomatr),IBEG1,IEND1)
            matr=1
            FOUND=.FALSE.
            DO WHILE(.NOT.FOUND.AND.matr.LE.NUMMATRMX)
              CALL STRING_TRIM(MATRNAME(matr),IBEG2,IEND2)
              IF(MATRICES(nomatr)(IBEG1:IEND1).EQ.
     '          MATRNAME(matr)(IBEG2:IEND2).AND.
     '          MATRTYPE(matr).EQ.1) THEN
                FOUND=.TRUE.
              ELSE
                matr=matr+1
              ENDIF
            ENDDO            
            IF(.NOT.FOUND) THEN
              WRITE(ERROR,'('' >>Matrix '',A,'' is unknown'')')
     '          MATRICES(nomatr)(IBEG1:IEND1)
              GOTO 9999
            ENDIF

            IF(matr.EQ.MATR_GK) THEN
              CALL IO_MATRIX('READ ','DISK',ISC_GK,ISR_GK,IUNIT,
     '          NYT(1,1,nx),NYT(2,1,nx),NYROWM,NZT(1,nx),
     '          NYROWM,NYROWM,NZ_GK_M,ROWLIST_DUMMY,
     '          KTYP24,GK,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GQ) THEN
              CALL IO_MATRIX('READ ','DISK',ISC_GQ,ISR_GQ,IUNIT,
     '          NYT(1,2,nx),NYT(2,2,nx),NYROWM,NZT(2,nx),
     '          NYROWM,NYROWM,NZ_GQ_M,ROWLIST_DUMMY,
     '          KTYP24,GQ,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GD) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_GM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_GKK) THEN
              CALL IO_MATRIX('READ ','DISK',ISC_GKK,ISR_GKK,IUNIT,
     '          NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOM,NZZT(1,nr,nx),
     '          NOM,NOM,NZ_GKK_M,ROWLIST_DUMMY,SPARSEGKK(nx),
     '          GKK,FILEFORMAT,
     '          INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GMM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_YP) THEN
              SPARSITY_TYPE=-1
C 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          M_DUMMY,N_DUMMY,NYM,NZ_DUMMY,NYM,NIYM,NZ_DUMMY,
C     '          ROWLIST_DUMMY,
C     '          SPARSITY_TYPE,YP,FILEFORMAT,INOUTTYPE,.TRUE.,
C     '          ERROR,*9999)
              NZ_MAX=NYM*NIYM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          M_DUMMY,N_DUMMY,NYM,NZ_DUMMY,NYM,NIYM,NZ_MAX,
     '          ROWLIST_DUMMY,
     '          SPARSITY_TYPE,YP,FILEFORMAT,INOUTTYPE,.TRUE.,
     '          ERROR,*9999)

            ELSE IF(matr.EQ.MATR_T_BH) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              TBH_OK=.FALSE.
              SPARSITY_TYPE=-1

C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(1),ISIZE_TBH(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          T_BH,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NY_TRANSFER_M
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(1),ISIZE_TBH(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
C AJP 26-3-97 Checking on the size of the matrix should be done inside
C IO_MATRIX, but there may not be enough information to do that there.
              CALL ASSERT(ISIZE_TBH(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              EVALUATE_TRANSFER=.TRUE.
            ELSE IF(matr.EQ.MATR_T_BH_INV) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(2),ISIZE_TBH(1),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          T_BH_INV,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NY_TRANSFER_M
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_TBH(1),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH_INV,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

C AJP 26-3-97 Checking on the size of the matrix should be done inside
C IO_MATRIX, but there may not be enough information to do that there.
              CALL ASSERT(ISIZE_TBH(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              EVALUATE_INVERSE=.TRUE.

            ELSE IF(matr.EQ.MATR_PHI) THEN

C LKC 16-JUN-2001 USE_TRANSFER does not really need to be set
C
C              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
C     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1


C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_PHI(1),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          PHI,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NTSM

C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_PHI(1),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          PHI,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_PHI(1),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NTSM,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

              CALL ASSERT(ISIZE_PHI(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              NTST=ISIZE_PHI(2)
              EVALUATE_PHI=.TRUE.
            ELSE IF(matr.EQ.MATR_ZCROSSING) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(1),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          ZCROSSING,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

              NZ_MAX=NY_TRANSFER_M*NTSM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          ZCROSSING,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

              CALL ASSERT(ISIZE_TBH(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              NTST=ISIZE_PHI(2)
            ELSE IF(matr.EQ.MATR_PHI_H) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          PHI_H,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NY_TRANSFER_M
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              EVALUATE_INVERSE=.TRUE.
CC AJPe
            ELSE IF(matr.EQ.MATR_YQ) THEN
              SPARSITY_TYPE=-1

C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          M_DUMMY,N_DUMMY,NYQM,NZ_DUMMY,NYQM,NIQM,NZ_DUMMY,
C     '          ROWLIST_DUMMY,
C     '          SPARSITY_TYPE,YQ,FILEFORMAT,INOUTTYPE,.TRUE.,
C     '          ERROR,*9999)
              NZ_MAX=NYQM*NIQM*NAM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          M_DUMMY,N_DUMMY,NYQM,NZ_DUMMY,NYQM,NIQM,NZ_MAX,
     '          ROWLIST_DUMMY,
     '          SPARSITY_TYPE,YQ,FILEFORMAT,INOUTTYPE,.TRUE.,
     '          ERROR,*9999)
            ELSE IF(matr.EQ.MATR_YQS) THEN
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          M_DUMMY,N_DUMMY,NIQSM,NZ_DUMMY,NIQSM,NQM,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,YQS,FILEFORMAT,INOUTTYPE,
C     '          .TRUE.,ERROR,*9999)
              NZ_MAX=NIQSM*NQM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          M_DUMMY,N_DUMMY,NIQSM,NZ_DUMMY,NIQSM,NQM,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,YQS,FILEFORMAT,INOUTTYPE,
     '          .TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_PHI_H_EXACT) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
              NZ_MAX=NY_TRANSFER_M*NY_TRANSFER_M
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H_EXACT,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              EVALUATE_PHI_H_EXACT=.TRUE.
            ELSE IF(matr.EQ.MATR_MFI) THEN
              CALL ASSERT(USE_MAGNETIC.EQ.1,'>> SEt USE_MAGNETIC to 0',
     '          ERROR,*9999)
C LKC 11-JUL-2005 We must always have 3 magnetic components (even if
C dealing with a "2d" problem.
C              ISIZE_MFI(3)=NJT
              ISIZE_MFI(3)=3
              NZ_MAX=NDM*NTSM
              
              DO nj=1,ISIZE_MFI(3)                
                FMT='(/'' Component '',I2)'
                READ(IUNIT,FMT) ISC_DUMMY(1)
                CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '            ISIZE_MFI(1),ISIZE_MFI(2),NDM,
     '            NZ_DUMMY,NDM,NTSM,NZ_MAX,ROWLIST_DUMMY,SPARSITY_TYPE,
     '            MFI(1,1,nj),FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              ENDDO

            ELSE
              ERROR='Unknown matrix type'
              GOTO 9999
            ENDIF

          ENDDO !nc

C         Read in vector values
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
            IF(NUMVECTORS.GT.0) THEN
C            READ(IUNIT,'(A)') LINE
            ENDIF
            DO novect=1,NUMVECTORS
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(/'' Vectors:'')')
              READ(IUNIT,'(A)') LINE
              CALL STRING_TRIM(VECTORS(novect),IBEG1,IEND1)
              vect=1
              FOUND=.FALSE.
              DO WHILE(.NOT.FOUND.AND.vect.LE.NUMMATRMX)
                CALL STRING_TRIM(MATRNAME(vect),IBEG2,IEND2)
                IF(VECTORS(novect)(IBEG1:IEND1).EQ.
     '            MATRNAME(vect)(IBEG2:IEND2).AND.
     '            MATRTYPE(vect).EQ.0) THEN
                  FOUND=.TRUE.
                ELSE
                  vect=vect+1
                ENDIF
              ENDDO

              IF(.NOT.FOUND) THEN
                WRITE(ERROR,'('' >>Vector '',A,'' is unknown'')')
     '            VECTORS(novect)(IBEG1:IEND1)
                GOTO 9999
              ENDIF

              IF(vect.EQ.MATR_GR) THEN
                CALL IO_VECTOR('READ ','DISK',IUNIT,NYT(1,1,nx),
     '            GR,NYROWM,FILEFORMAT,INOUTTYPE,
     '            ERROR,*9999)
              ELSE IF(vect.EQ.MATR_GRR) THEN
                CALL IO_VECTOR('READ ','DISK',IUNIT,NOT(1,1,nr,nx),
     '            GRR,NOM,FILEFORMAT,INOUTTYPE,
     '            ERROR,*9999)
              ELSE IF(vect.EQ.MATR_LD_NP) THEN
C                CALL IO_VECTOR_INT('READ ','DISK',IUNIT,NDT,LD_NP,NDM,
C               &            FILEFORMAT,INOUTTYPE,ERROR,*9999)
 !KSB - this now specific for mapping arteries to veins (15/06/05)
                CALL ASSERT(USE_DATA.EQ.1,'>>Must set USE_DATA to 1!',
     &            ERROR,*9999)
                CALL ASSERT(USE_LUNG.EQ.1,'>>Must set USE_LUNG to 1!',
     &            ERROR,*9999)
                CALL ASSERT(USE_GRID.EQ.1,'>>Must set USE_GRID to 1!',
     &            ERROR,*9999)
                MAP_ART_VEIN(0,nr)=0
                CALL IO_VECTOR_INT('READ ','DISK',IUNIT,
     &            MAP_ART_VEIN(0,nr),MAP_ART_VEIN(1,nr),NDM,FILEFORMAT,
     &            INOUTTYPE,ERROR,*9999)                
              ENDIF
            ENDDO !nc
          ENDIF

C***      Close ASCII input file
          CLOSE(IUNIT)

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

C***      Open binary file

          FILEID=IUNIT !unit number
          FILETYPE=1 !Binary matrix file
          VERSION(1)=1 !File version at 22/10/95
          VERSION(2)=1
          VERSION(3)=0
          CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '      'READ','mat',FILE_NAME,ERROR,*9999)

C***      Read matrices (tags)

          DO TAG=1,NUMBERTAGS

            CALL READ_BIN_TAG_HEADER(FILEID,ERROR,*9999)

            IF(TAGINDEX.EQ.MATR_GK) THEN
              CALL IO_MATRIX('READ ','DISK',ISC_GK,ISR_GK,FILEID,
     '          NYT(1,1,nx),NYT(2,1,nx),NYROWM,NZT(1,nx),
     '          NYROWM,NYROWM,NZ_GK_M,ROWLIST_DUMMY,
     '          KTYP24,GK,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_GQ) THEN
              CALL IO_MATRIX('READ ','DISK',ISC_GQ,ISR_GQ,FILEID,
     '          NYT(1,2,nx),NYT(2,2,nx),NYROWM,NZT(2,nx),
     '          NYROWM,NYROWM,NZ_GQ_M,ROWLIST_DUMMY,
     '          KTYP24,GQ,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_GD) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(TAGINDEX.EQ.MATR_GM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(TAGINDEX.EQ.MATR_GR) THEN
              CALL IO_VECTOR('READ ','DISK',FILEID,NYT(1,1,nx),
     '          GR,NYROWM,FILEFORMAT,INOUTTYPE,
     '          ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_GKK) THEN
              CALL IO_MATRIX('READ ','DISK',ISC_GKK,ISR_GKK,FILEID,
     '          NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOM,NZZT(1,nr,nx),
     '          NOM,NOM,NZ_GKK_M,ROWLIST_DUMMY,SPARSEGKK(nx),
     '          GKK,FILEFORMAT,
     '          INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_GMM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(TAGINDEX.EQ.MATR_GRR) THEN
              CALL IO_VECTOR('READ ','DISK',FILEID,NOT(1,1,nr,nx),
     '          GRR,NOM,FILEFORMAT,INOUTTYPE,
     '          ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_YP) THEN
              SPARSITY_TYPE=-1

C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          M_DUMMY,N_DUMMY,NYM,NZ_DUMMY,NYM,NIYM,NZ_DUMMY,
C     '          ROWLIST_DUMMY,
C     '          SPARSITY_TYPE,YP,FILEFORMAT,INOUTTYPE,.TRUE.,
C     '          ERROR,*9999)
              NZ_MAX=NYM*NIYM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          M_DUMMY,N_DUMMY,NYM,NZ_DUMMY,NYM,NIYM,NZ_MAX,
     '          ROWLIST_DUMMY,
     '          SPARSITY_TYPE,YP,FILEFORMAT,INOUTTYPE,.TRUE.,
     '          ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_T_BH) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(1),ISIZE_TBH(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          T_BH,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NY_TRANSFER_M
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(1),ISIZE_TBH(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
C AJP 26-3-97 Checking on the size of the matrix should be done inside
C IO_MATRIX, but there may not be enough information to do that there.
              CALL ASSERT(ISIZE_TBH(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              EVALUATE_TRANSFER=.TRUE.
            ELSE IF(TAGINDEX.EQ.MATR_T_BH_INV) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(2),ISIZE_TBH(1),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          T_BH_INV,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NY_TRANSFER_M
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_TBH(1),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH_INV,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
C AJP 26-3-97 Checking on the size of the matrix should be done inside
C IO_MATRIX, but there may not be enough information to do that there.
              CALL ASSERT(ISIZE_TBH(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              EVALUATE_INVERSE=.TRUE.
CC AJPs 191297
            ELSE IF(TAGINDEX.EQ.MATR_PHI) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_PHI(1),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          PHI,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NTSM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_PHI(1),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              NTST=ISIZE_PHI(2)
              EVALUATE_PHI=.TRUE.
            ELSE IF(TAGINDEX.EQ.MATR_ZCROSSING) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(1),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          ZCROSSING,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NTSM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          ZCROSSING,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(1).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              NTST=ISIZE_PHI(2)
            ELSE IF(TAGINDEX.EQ.MATR_PHI_H) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          PHI_H,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              NZ_MAX=NY_TRANSFER_M*NTSM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              EVALUATE_INVERSE=.TRUE.
CC AJPe
            ELSE IF(TAGINDEX.EQ.MATR_YQ) THEN
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          M_DUMMY,N_DUMMY,NYQM,NZ_DUMMY,NYQM,NIQM,NZ_DUMMY,
C     '          ROWLIST_DUMMY,
C     '          SPARSITY_TYPE,YQ,FILEFORMAT,INOUTTYPE,.TRUE.,
C     '          ERROR,*9999)
              NZ_MAX=NYQM*NIQM*NAM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          M_DUMMY,N_DUMMY,NYQM,NZ_DUMMY,NYQM,NIQM,NZ_MAX,
     '          ROWLIST_DUMMY,
     '          SPARSITY_TYPE,YQ,FILEFORMAT,INOUTTYPE,.TRUE.,
     '          ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_YQS) THEN
              SPARSITY_TYPE=-1
C LKC 23-JAN-1999 Incorrect passing of NZ_DUMMY
C              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          M_DUMMY,N_DUMMY,NIQSM,NZ_DUMMY,NIQSM,NQM,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,YQS,FILEFORMAT,INOUTTYPE,
C     '          .TRUE.,ERROR,*9999)
              NZ_MAX=NIQSM*NQM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          M_DUMMY,N_DUMMY,NIQSM,NZ_DUMMY,NIQSM,NQM,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,YQS,FILEFORMAT,INOUTTYPE,
     '          .TRUE.,ERROR,*9999)
            ELSE IF(TAGINDEX.EQ.MATR_PHI_H_EXACT) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
              NZ_MAX=NY_TRANSFER_M*NTSM
              CALL IO_MATRIX('READ ','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_MAX,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H_EXACT,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              CALL ASSERT(ISIZE_TBH(2).LE.NY_TRANSFER_M,
     '          '>>NY_TRANSFER_M too small to hold matrix',ERROR,*9999)
              CALL ASSERT(ISIZE_PHI(2).LE.NTSM,
     '          '>>NTSM too small to hold matrix',ERROR,*9999)
              EVALUATE_PHI_H_EXACT=.TRUE.
            ELSE
              WRITE(OP_STRING,'('' >>Warning: Unknown matrix code '','
     '          //'I2,'' ('',A,'')'')') TAGINDEX,
     '          TAGHEADER(1:TAGHEADERBYTES)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL BINSKIPFILE(IUNIT,NUMTAGBYTES,ERROR,*9999)
            ENDIF

            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' Matrix: code '',I2,'' ('',A,'') '
     '          //'read'')') TAGINDEX,TAGHEADER(1:TAGHEADERBYTES)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

          ENDDO !tag

C***      Close binary file
          CALL BINCLOSEFILE(IUNIT,ERROR,*9999)

        ENDIF

      ELSE IF(COMAND(1:5).EQ.'WRITE') THEN

        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
          VERSION(1)=2
          IOTYPE=3
          CALL OPEN_SEQ_FILE(VERSION(1),IUNIT,FILE_NAME,
     '      'matr','NEW',.FALSE.,ERROR,*9999)
C         Write out header information
          WRITE(IUNIT,'('' Number of matrices: '',I1)') MATRLIST(0)
          IF(MATRLIST(0).GT.0) THEN
CCs AJP 191297
C            WRITE(IUNIT,'('' Matrices are:'',9(1X,A8))')
C     '        (MATRNAME(MATRLIST(matr)),matr=1,MATRLIST(0))
            WRITE(IUNIT,'('' Matrices are:'',7(1X,A11))')
     '        (MATRNAME(MATRLIST(matr)),matr=1,MATRLIST(0))
CC AJP
          ENDIF
          WRITE(IUNIT,'(/'' Number of vectors : '',I1)') VECTLIST(0)
          IF(VECTLIST(0).GT.0) THEN
            WRITE(IUNIT,'('' Vectors are :'',7(1X,A11))')
     '        (MATRNAME(VECTLIST(vect)),vect=1,VECTLIST(0))
          ENDIF
          IF(MATRLIST(0).GT.0) THEN
            WRITE(IUNIT,'(/'' Matrices:'')')
          ENDIF
          DO nomatr=1,MATRLIST(0)
            matr=MATRLIST(nomatr)
            WRITE(IUNIT,'(/'' Matrix '',I2,'':'')') nomatr
            IF(matr.EQ.MATR_GK) THEN
              CALL IO_MATRIX('WRITE','DISK',ISC_GK,ISR_GK,IUNIT,
     '          NYT(1,1,nx),NYT(2,1,nx),NYROWM,NZT(1,nx),
     '          NYROWM,NYROWM,NZ_GK_M,ROWLIST_DUMMY,
     '          KTYP24,GK,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GQ) THEN
              CALL IO_MATRIX('WRITE','DISK',ISC_GQ,ISR_GQ,IUNIT,
     '          NYT(1,2,nx),NYT(2,2,nx),NYROWM,NZT(2,nx),
     '          NYROWM,NYROWM,NZ_GQ_M,ROWLIST_DUMMY,
     '          KTYP24,GQ,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GD) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_GM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_GKK) THEN
              CALL IO_MATRIX('WRITE','DISK',ISC_GKK,ISR_GKK,IUNIT,
     '          NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOM,NZZT(1,nr,nx),
     '          NOM,NOM,NZ_GKK_M,ROWLIST_DUMMY,SPARSEGKK(nx),
     '          GKK,FILEFORMAT,
     '          INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GMM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_YP) THEN
              NYTT=0
              DO nc=1,NCT(nr,nx)
                NYTT=NYTT+NYT(1,nc,nx)
              ENDDO !nc
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          NYTT,NIYM,NYM,NZ_DUMMY,NYM,NIYM,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,YP,
     '          FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_T_BH) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(1),ISIZE_TBH(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_T_BH_INV) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
     '          ERROR,*9999)
              CALL ASSERT(ICALC_TRANSFER.EQ.1,
     '          '>>Transfer matrix has not been inverted explicitly',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_TBH(1),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH_INV,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_PHI) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_PHI,'>>Evaluate PHI first',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_PHI(1),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_ZCROSSING) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          ZCROSSING,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_PHI_H) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
C LKC 18-JUL-2002 This is not strictly required
C              CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
C     '          ERROR,*9999)
              SPARSITY_TYPE=-1

C LKC 19-JUl-2002 the rows and cols are wrong .....
C              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
C     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
C     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
C     '          ROWLIST_DUMMY,SPARSITY_TYPE,
C     '          PHI_H,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_PHIH(1),ISIZE_PHIH(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

            ELSE IF(matr.EQ.MATR_YQ) THEN
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          NYQM,NIQM,NYQM,NZ_DUMMY,NYQM,NIQM,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,YQ,
     '          FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

            ELSE IF(matr.EQ.MATR_YQS) THEN
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          NIQSM,NIQSM,NQM,NZ_DUMMY,NIQSM,NQM,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,YQS,FILEFORMAT,INOUTTYPE,
     '          .TRUE.,ERROR,*9999)

            ELSE IF(matr.EQ.MATR_PHI_H_EXACT) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_PHI_H_EXACT,
     '          '>>Evaluate PHI_H_EXACT first',ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H_EXACT,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

            ELSE IF(matr.EQ.MATR_MFI) THEN
              SPARSITY_TYPE=-1

C LKC 11-JUL-2005 We must always have 3 magnetic components (even if
C dealing with a "2d" problem.
C              DO nj=1,NJT
              DO nj=1,3
                WRITE(IUNIT,'(/'' Component '',I2)') nj
                CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '            ISIZE_MFI(1),ISIZE_MFI(2),NDM,
     '            NZ_DUMMY,NDM,NTSM,NZ_DUMMY,
     '            ROWLIST_DUMMY,SPARSITY_TYPE,
     '            MFI(1,1,nj),FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
              ENDDO

            ELSE
              ERROR='Unknown matrix type'
              GOTO 9999
            ENDIF
          ENDDO !nomatr

C         Write out vector values
          IF((INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT').AND.
     '      (VECTLIST(0).GT.0)) THEN
            WRITE(IUNIT,'(/'' Vectors:'')')
            DO novect=1,VECTLIST(0)
              vect=VECTLIST(novect)
              WRITE(IUNIT,'(/'' Vector '',I2,'':'')') novect
              IF(vect.EQ.MATR_GR) THEN
                CALL IO_VECTOR('WRITE','DISK',IUNIT,NYT(1,1,nx),
     '            GR,NYROWM,FILEFORMAT,INOUTTYPE,
     '            ERROR,*9999)
              ELSE IF(vect.EQ.MATR_GRR) THEN
                CALL IO_VECTOR('WRITE','DISK',IUNIT,NOT(1,1,nr,nx),
     '            GRR,NOM,FILEFORMAT,INOUTTYPE,
     '            ERROR,*9999)
              ELSE IF(vect.EQ.MATR_LD_NP) THEN
                CALL IO_VECTOR_INT('WRITE ','DISK',IUNIT,NDT,LD_NP,NDM,
     &            FILEFORMAT,INOUTTYPE,ERROR,*9999)    
              ENDIF
            ENDDO !novect
          ENDIF

C***      Close ASCII output file
          CLOSE(IUNIT)

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

C*** Open binary file

          FILEID=IUNIT !unit number
          FILETYPE=1 !Binary matrix file
          VERSION(1)=1 !Version 22/10/95
          VERSION(2)=1
          VERSION(3)=0
          NUMBERTAGS=MATRLIST(0)+VECTLIST(0)
          CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,
     '      'WRITE','mat',FILE_NAME,ERROR,*9999)

C*** Write tags (matrices)

          DO nomatr=1,MATRLIST(0)
            matr=MATRLIST(nomatr)

            TAGINDEX=matr
            TAGHEADER=MATRNAME(matr)
            NUMBERSUBTAGS=0

            IF(matr.EQ.MATR_GK) THEN
              NUMTAGBYTES=7*INTSIZE
              IF(KTYP24.GT.0) THEN
                IF(INOUTTYPE(1:6).NE.'ARRAYS') THEN
                  IF(KTYP24.EQ.1) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(NYT(1,1,nx)+1+NZT(1,nx))*
     '                INTSIZE
                  ELSE IF(KTYP24.EQ.2) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(1+2*NZT(1,nx))*INTSIZE
                  ELSE IF(KTYP24.EQ.3) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(NYT(2,1,nx)+1+NZT(1,nx))*
     '                INTSIZE
                  ENDIF
                ELSE
                  ERROR='>>Cannot just write the values for a sparse '
     '              //'matrix'
                  GOTO 9999
                ENDIF
              ENDIF
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+NZT(1,nx)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              CALL IO_MATRIX('WRITE','DISK',ISC_GK,ISR_GK,FILEID,
     '          NYT(1,1,nx),NYT(2,1,nx),NYROWM,NZT(1,nx),
     '          NYROWM,NYROWM,NZ_GK_M,ROWLIST_DUMMY,
     '          KTYP24,GK,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GQ) THEN
              NUMTAGBYTES=7*INTSIZE
              IF(KTYP24.GT.0.) THEN
                IF(INOUTTYPE(1:6).NE.'ARRAYS') THEN
                  IF(KTYP24.EQ.1) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(NYT(1,2,nx)+1+NZT(2,nx))*
     '                INTSIZE
                  ELSE IF(KTYP24.EQ.2) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(1+2*NZT(2,nx))*INTSIZE
                  ELSE IF(KTYP24.EQ.1) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(NYT(2,2,nx)+1+NZT(2,nx))*
     '                INTSIZE
                  ENDIF
                ELSE
                  ERROR='>>Cannot just write the values for a sparse '
     '              //'matrix'
                  GOTO 9999
                ENDIF
              ENDIF
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+NZT(2,nx)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              CALL IO_MATRIX('WRITE','DISK',ISC_GQ,ISR_GQ,FILEID,
     '          NYT(1,2,nx),NYT(2,2,nx),NYROWM,NZT(2,nx),
     '          NYROWM,NYROWM,NZ_GQ_M,ROWLIST_DUMMY,
     '          KTYP24,GQ,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GD) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_GM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_GKK) THEN
              NUMTAGBYTES=7*INTSIZE
              IF(SPARSEGKK(nx).GT.0) THEN
                IF(INOUTTYPE(1:6).NE.'ARRAYS') THEN
                  IF(SPARSEGKK(nx).EQ.1) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(NOT(1,1,nr,nx)+1+
     '                NZZT(1,nr,nx))*INTSIZE
                  ELSE IF(SPARSEGKK(nx).EQ.2.OR.SPARSEGKK(nx).EQ.4
     '              .OR.SPARSEGKK(nx).EQ.5) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(1+2*NZZT(1,nr,nx))*INTSIZE
                  ELSE IF(SPARSEGKK(nx).EQ.3) THEN
                    NUMTAGBYTES=NUMTAGBYTES+(NOT(2,1,nr,nx)+1+
     '                NZZT(1,nr,nx))*INTSIZE
                  ENDIF
                ELSE
                  ERROR='>>Cannot just write the values for a sparse '
     '              //'matrix'
                  GOTO 9999
                ENDIF
              ENDIF
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+NZZT(1,nr,nx)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              CALL IO_MATRIX('WRITE','DISK',ISC_GKK,ISR_GKK,FILEID,
     '          NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOM,NZZT(1,nr,nx),
     '          NOM,NOM,NZ_GKK_M,ROWLIST_DUMMY,SPARSEGKK(nx),
     '          GKK,FILEFORMAT,
     '          INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_GMM) THEN
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(matr.EQ.MATR_YP) THEN
              NUMTAGBYTES=7*INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NYTT=0
                DO nc=1,NCT(nr,nx)
                  NYTT=NYTT+NYT(1,nc,nx)
                ENDDO !nc
                NUMTAGBYTES=NUMTAGBYTES+(NYTT*NIYM)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,FILEID,
     '          NYTT,NIYM,NYM,NZ_DUMMY,NYM,NIYM,NZ_DUMMY,
     '          ROWLIST_DUMMY,
     '          SPARSITY_TYPE,YP,FILEFORMAT,INOUTTYPE,.TRUE.,
     '          ERROR,*9999)
            ELSE IF(matr.EQ.MATR_T_BH) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '          ERROR,*9999)
              NUMTAGBYTES=7*INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(ISIZE_TBH(1)*ISIZE_TBH(2))*
     '            DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(1),ISIZE_TBH(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_T_BH_INV) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
     '          ERROR,*9999)
              CALL ASSERT(ICALC_TRANSFER.EQ.1,
     '          '>>Transfer matrix has not been inverted explicitly',
     '          ERROR,*9999)
              NUMTAGBYTES=7*INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(ISIZE_TBH(1)*ISIZE_TBH(2))*
     '            DPSIZE
              ENDIF
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_TBH(1),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          T_BH_INV,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
CC AJPs 191297
            ELSE IF(matr.EQ.MATR_PHI) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_PHI,'>>Evaluate PHI first',
     '          ERROR,*9999)
              NUMTAGBYTES=INTSIZE+2*INTSIZE+INTSIZE+INTSIZE+INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(ISIZE_PHI(1)*ISIZE_PHI(2))*
     '            DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_PHI(1),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_ZCROSSING) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              NUMTAGBYTES=INTSIZE+2*INTSIZE+INTSIZE+INTSIZE+INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(ISIZE_TBH(1)*ISIZE_PHI(2))*
     '            DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          ZCROSSING,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_PHI_H) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
     '          ERROR,*9999)
              NUMTAGBYTES=INTSIZE+2*INTSIZE+INTSIZE+INTSIZE+INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(ISIZE_TBH(2)*ISIZE_PHI(2))*
     '            DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_YQ) THEN
              NUMTAGBYTES=7*INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(NYQM*NIQM)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,FILEID,
     '          NYQM,NIQM,NYQM,NZ_DUMMY,NYQM,NIQM,NZ_DUMMY,
     '          ROWLIST_DUMMY, SPARSITY_TYPE,YQ,FILEFORMAT,
     '          INOUTTYPE,.TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_YQS) THEN
              NUMTAGBYTES=7*INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(NIQSM*NQM)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,FILEID,
     '          NIQSM,NQM,NIQSM,NZ_DUMMY,NIQSM,NQM,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,YQS,FILEFORMAT,INOUTTYPE,
     '          .TRUE.,ERROR,*9999)
            ELSE IF(matr.EQ.MATR_PHI_H_EXACT) THEN
              CALL ASSERT(USE_TRANSFER.EQ.1,'>>USE_TRANSFER is 0',
     '          ERROR,*9999)
              CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME is 0',
     '          ERROR,*9999)
              CALL ASSERT(EVALUATE_PHI_H_EXACT,
     '          '>>Evaluate PHI_H_EXACT first',ERROR,*9999)
              NUMTAGBYTES=INTSIZE+2*INTSIZE+INTSIZE+INTSIZE+INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+(ISIZE_TBH(2)*ISIZE_PHI(2))*
     '            DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              SPARSITY_TYPE=-1
              CALL IO_MATRIX('WRITE','DISK',ISC_DUMMY,ISR_DUMMY,IUNIT,
     '          ISIZE_TBH(2),ISIZE_PHI(2),NY_TRANSFER_M,
     '          NZ_DUMMY,NY_TRANSFER_M,NY_TRANSFER_M,NZ_DUMMY,
     '          ROWLIST_DUMMY,SPARSITY_TYPE,
     '          PHI_H_EXACT,FILEFORMAT,INOUTTYPE,.TRUE.,ERROR,*9999)

            ELSE IF(matr.EQ.MATR_MFI) THEN
              ERROR='Not implemented for binary yet'
              GOTO 9999
            ELSE
              WRITE(OP_STRING,'('' >>Warning: Unknown matrix code '','
     '          //'I2,'' ('',A,'')'')') TAGINDEX,
     '          TAGHEADER(1:TAGHEADERBYTES)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL BINSKIPFILE(IUNIT,NUMTAGBYTES,ERROR,*9999)
            ENDIF

            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' Matrix: code '',I2,'' ('',A,'') '
     '          //'read'')') TAGINDEX,TAGHEADER(1:TAGHEADERBYTES)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

          ENDDO !nomatr (matr)

C*** Write tags (vectors)

          DO novect=1,VECTLIST(0)
            vect=VECTLIST(novect)

            TAGINDEX=vect
            TAGHEADER=MATRNAME(vect)
            NUMBERSUBTAGS=0

            IF(vect.EQ.MATR_GR) THEN
              NUMTAGBYTES=5*INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+NYT(1,1,nx)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              CALL IO_VECTOR('WRITE','DISK',FILEID,NYT(1,1,nx),
     '          GR,NYROWM,FILEFORMAT,INOUTTYPE,
     '          ERROR,*9999)
            ELSE IF(vect.EQ.MATR_GRR) THEN
              NUMTAGBYTES=5*INTSIZE
              IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '          INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
                NUMTAGBYTES=NUMTAGBYTES+NOT(1,1,nr,nx)*DPSIZE
              ENDIF
              CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
              CALL IO_VECTOR('WRITE','DISK',FILEID,NOT(1,1,nr,nx),
     '          GRR,NOM,FILEFORMAT,INOUTTYPE,
     '          ERROR,*9999)
            ELSE
              WRITE(OP_STRING,'('' >>Warning: Unknown matrix code '','
     '          //'I2,'' ('',A,'')'')') TAGINDEX,
     '          TAGHEADER(1:TAGHEADERBYTES)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL BINSKIPFILE(IUNIT,NUMTAGBYTES,ERROR,*9999)
            ENDIF

            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' Vector: code '',I2,'' ('',A,'') '
     '          //'read'')') TAGINDEX,CHARDATA(1:TAGHEADERBYTES)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

          ENDDO !novect (vect)

C***      Close binary file

          CALL BINCLOSEFILE(IUNIT,ERROR,*9999)

        ENDIF

      ELSE
        ERROR=' Command error: COMAND='//COMAND
        GOTO 9999
      ENDIF

      IOTYPE=IOTYPEOLD

      CALL EXITS('IOMATR')
      RETURN
 9999 CALL ERRORS('IOMATR',ERROR)
      CALL EXITS('IOMATR')
      RETURN 1
      END



