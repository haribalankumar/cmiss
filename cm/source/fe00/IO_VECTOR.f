      SUBROUTINE IO_VECTOR(COMAND,DEVICE,IUNIT,VECTOR_M,VECTOR,VMAX,
     '  FILEFORMAT,INOUTTYPE,ERROR,*)

C#### Subroutine: IO_VECTOR
C###  Description:
C###    IO_VECTOR reads and writes a single vector.

C     VECTOR_M = Number of rows in the vector

      IMPLICIT NONE
      INCLUDE 'binf00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER IUNIT,VECTOR_M,VMAX
      REAL*8 VECTOR(*)
      CHARACTER COMAND*(*),DEVICE*(*),ERROR*(*),
     '  FILEFORMAT*(*),INOUTTYPE*(*)
!     Local Variables
      INTEGER CERROR(50),ERR,NUMINDICES(1),
     '  NUMVALUES,NUMVALS(1),row,SPARSITY(1),TIMECODE(1)
      REAL*8 TIME(1)
      LOGICAL ISBINFILEOPEN
      CHARACTER FMT*500

      CALL ENTERS('IO_VECTOR',*9999)

      NUMVALUES=0 !initialising

C***  Reading

      IF(COMAND.EQ.'READ') THEN

C       Must read from disk
        CALL ASSERT(DEVICE(1:4).EQ.'DISK','>>Invalid device for read',
     '    ERROR,*9999)

        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN

C         Read in vector values
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C            READ(IUNIT,'('' Number of rows='',I5)') VECTOR_M
            FMT='('' Number of rows='',I5)'
            READ(IUNIT,FMT) VECTOR_M
            CALL ASSERT(VECTOR_M.LE.VMAX,
     '        '>>Array too small to hold data',ERROR,*9999)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C            READ(IUNIT,'(5(1X,D15.8),:/(5(1X,D15.8)))') (VECTOR(row),
C     '        row=1,VECTOR_M)
            FMT='(5(1X,D15.8),:/(5(1X,D15.8)))'
            READ(IUNIT,FMT) (VECTOR(row),row=1,VECTOR_M)
          ENDIF !array/all_inout

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

C         Read number of indicies
          CALL BINREADFILE(IUNIT,INTTYPE,1,NUMINDICES,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          CALL ASSERT(NUMINDICES(1).EQ.1,
     '      '>>Number of indices <> 1, not implemented',ERROR,*9999)
C         Read the index values
          CALL BINREADFILE(IUNIT,INTTYPE,NUMINDICES(1),INTDATA,
     '      REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          VECTOR_M=INTDATA(1)
C         Read the size (same as M but consistent with a matrix)
          CALL BINREADFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          CALL ASSERT(VECTOR_M.EQ.INTDATA(1),'>>Invalid vector size',
     '      ERROR,*9999)
C         Read the number of values stored
          NUMVALS(1)=NUMVALUES
          CALL BINREADFILE(IUNIT,INTTYPE,1,NUMVALS,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NUMVALUES=NUMVALS(1)
          CALL ASSERT(NUMVALUES.LE.VMAX,
     '      'ERROR: Array too small to hold data',ERROR,*9999)
C         Read the timecode
          CALL BINREADFILE(IUNIT,INTTYPE,1,TIMECODE,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          IF(TIMECODE(1).EQ.1) THEN
C           Read the time
            CALL BINREADFILE(IUNIT,DPTYPE,1,INTDATA,REAL4DATA,
     '        TIME,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          ENDIF
C         Read the sparsity type
          CALL BINREADFILE(IUNIT,INTTYPE,1,SPARSITY,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          CALL ASSERT(SPARSITY(1).EQ.0,'>>Invalid sparisty format',
     '      ERROR,*9999)
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
C           Read the values
            CALL BINREADFILE(IUNIT,DPTYPE,NUMVALUES,INTDATA,
     '        REAL4DATA,VECTOR,CHARDATA,LOGDATA,SINTDATA,
     '        ERROR,*9999)
          ELSE
            IF(NUMVALUES.GT.0) THEN
C             Skip the values
              CALL BINSKIPFILE(IUNIT,NUMVALUES*DPSIZE,ERROR,*9999)
            ENDIF
          ENDIF

        ENDIF

C***  Writing

      ELSE IF(COMAND.EQ.'WRITE') THEN

        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN

          IF(DEVICE(1:4).EQ.'DISK') THEN

C           Write out vector values
            IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '        INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
              WRITE(IUNIT,'('' Number of rows='',I5)') VECTOR_M
              WRITE(IUNIT,'(5(1X,D15.8),:/(5(1X,D15.8)))') (VECTOR(row),
     '          row=1,VECTOR_M)
            ENDIF               !array/all_inout

          ELSE IF(DEVICE(1:4).EQ.'TERM') THEN

            WRITE(OP_STRING,'('' Number of rows='',I5)') VECTOR_M
            CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(1X,5D12.4,:/(1X,5D12.4))') (VECTOR(row),
     '        row=1,VECTOR_M)
            CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

          CALL ASSERT(DEVICE(1:4).EQ.'DISK',
     '      '>>Invalid device for binary write',ERROR,*9999)

C         Write the number of indicies
          INTDATA(1)=1
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the index values
          INTDATA(1)=VECTOR_M
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the size (same as M but consistent with a matrix)
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the number of values stored
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
            NUMVALUES=VECTOR_M
          ELSE
            NUMVALUES=0
          ENDIF
          NUMVALS(1)=NUMVALUES
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,NUMVALS,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NUMVALUES=NUMVALS(1)
C         Write the timecode (not implemented at the moment)
          INTDATA(1)=0
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the time (for future development)
C         Write the sparsity format
          INTDATA(1)=0
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the values stored
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
            CALL BINWRITEFILE(IUNIT,DPTYPE,VECTOR_M,INTDATA,REAL4DATA,
     '        VECTOR,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          ENDIF

        ENDIF
      ENDIF !read/write

      CALL EXITS('IO_VECTOR')
      RETURN
 9999 CALL ERRORS('IO_VECTOR',ERROR)
      CALL EXITS('IO_VECTOR')
      IF(ISBINFILEOPEN(IUNIT)) CALL BINARYCLOSEFILE(IUNIT,ERR,CERROR)
      RETURN 1
      END


