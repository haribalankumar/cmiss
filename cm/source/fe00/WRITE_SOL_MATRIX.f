      SUBROUTINE WRITE_SOL_MATRIX(ISC_GKK,ISR_GKK,nr,nx,GKK,GRR,ERROR,*)

C#### Subroutine: WRITE_SOL_MATRIX
C###  Description:
C###    Writes out a matrix file for the solution matrix and vector.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'binf00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'matr00.inc'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER ISC_GKK(NISC_GKKM),ISR_GKK(NISR_GKKM),nr,nx
      REAL*8 GKK(NZ_GKK_M),GRR(NOM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER FILEID,FILETYPE,IUNIT,OLDIOTYPE,
     '  ROWLIST_DUMMY(1),VERSION(3)

      CALL ENTERS('WRITE_SOL_MATRIX',*9999)

      ROWLIST_DUMMY(1)=0
      IF(KTYP4.GT.0) THEN !Ascii file
        OLDIOTYPE=IOTYPE
        IOTYPE=3
        IUNIT=20
        VERSION(1)=2
        CALL OPEN_SEQ_FILE(VERSION(1),IUNIT,FILE00,'matr',
     '    'NEW',.FALSE.,ERROR,*9999)
        IF(KTYP4.EQ.3) THEN !GKK and GRR
          WRITE(IUNIT,'('' Number of matrices: 1'')')
          WRITE(IUNIT,'('' Matrices are: GKK'')')
          WRITE(IUNIT,'(/'' Number of vectors : 1'')')
          WRITE(IUNIT,'('' Vectors are : GRR'')')
        ELSE IF(KTYP4.EQ.2) THEN !Just GRR
          WRITE(IUNIT,'('' Number of matrices: 0'')')
          WRITE(IUNIT,'(/'' Number of vectors : 1'')')
          WRITE(IUNIT,'('' Vectors are : GRR'')')
        ELSE IF(KTYP4.EQ.1) THEN!Just GKK
          WRITE(IUNIT,'('' Number of matrices: 1'')')
          WRITE(IUNIT,'('' Matrices are: GKK'')')
          WRITE(IUNIT,'(/'' Number of vectors : 0'')')
        ELSE
          ERROR='>>Invalid KTYP4'
          GOTO 9999
        ENDIF
        IF(KTYP4.EQ.1.OR.KTYP4.EQ.3) THEN
          WRITE(IUNIT,'(/'' Matrices:'')')
          WRITE(IUNIT,'(/'' Matrix  1:'')')

          CALL IO_MATRIX('WRITE','DISK',ISC_GKK,ISR_GKK,IUNIT,
     '      NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOM,NZZT(1,nr,nx),
     '      NOM,NOM,NZ_GKK_M,ROWLIST_DUMMY,SPARSEGKK(nx),
     '      GKK,'ASCII ','ALL_INOUT',
     '      .TRUE.,ERROR,*9999)

        ENDIF
        IF(KTYP4.EQ.2.OR.KTYP4.EQ.3) THEN
          WRITE(IUNIT,'(/'' Vectors:'')')
          WRITE(IUNIT,'(/'' Vector  1:'')')

          CALL IO_VECTOR('WRITE','DISK',IUNIT,NOT(1,1,nr,nx),
     '      GRR,NOM,'ASCII ','ALL_INOUT',
     '      ERROR,*9999)

        ENDIF

        CLOSE(IUNIT)
        IOTYPE=OLDIOTYPE

      ELSE IF(KTYP4.LT.0) THEN !Binary file

C*** Open binary file

        FILEID=20 !unit number
        FILETYPE=1 !Binary matrix file
        VERSION(1)=1 !Version 22/5/95
        VERSION(2)=1
        VERSION(3)=0
        IF(KTYP4.EQ.-3) THEN
          NUMBERTAGS=2
        ELSE
          NUMBERTAGS=1
        ENDIF
        CALL OPEN_BIN_FILE(FILEID,FILETYPE,NUMBERTAGS,VERSION,'WRITE',
     '    'mat',FILE00,ERROR,*9999)

        IF(KTYP4.EQ.-1.OR.KTYP4.EQ.-3) THEN
C***      Write out GKK
          TAGINDEX=MATR_GKK
          TAGHEADER='GKK'
          NUMBERSUBTAGS=0
          NUMTAGBYTES=7*INTSIZE
          IF(SPARSEGKK(nx).EQ.1) THEN
            NUMTAGBYTES=NUMTAGBYTES+(NOT(1,1,nr,nx)+1+NZZT(1,nr,nx))*
     '        INTSIZE
          ELSE IF(SPARSEGKK(nx).EQ.2.OR.SPARSEGKK(nx).EQ.4
     '        .OR.SPARSEGKK(nx).EQ.5) THEN
            NUMTAGBYTES=NUMTAGBYTES+(1+2*NZZT(1,nr,nx))*INTSIZE
          ELSE IF(SPARSEGKK(nx).EQ.3) THEN
            NUMTAGBYTES=NUMTAGBYTES+(NOT(2,1,nr,nx)+1+NZZT(1,nr,nx))*
     '        INTSIZE
          ENDIF
          NUMTAGBYTES=NUMTAGBYTES+NZZT(1,nr,nx)*DPSIZE
          CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
          CALL IO_MATRIX('WRITE','DISK',ISC_GKK,ISR_GKK,FILEID,
     '      NOT(1,1,nr,nx),NOT(2,1,nr,nx),NOM,NZZT(1,nr,nx),
     '      NOM,NOM,NZ_GKK_M,ROWLIST_DUMMY,SPARSEGKK(nx),
     '      GKK,'BINARY','ALL_INOUT',
     '      .TRUE.,ERROR,*9999)
        ENDIF
        IF(KTYP4.EQ.-2.OR.KTYP4.EQ.-3) THEN
C***      Write out GRR
          TAGINDEX=MATR_GRR
          TAGHEADER='GRR'
          NUMBERSUBTAGS=0
          NUMTAGBYTES=6*INTSIZE+NOT(1,1,nr,nx)*DPSIZE
          CALL WRITE_BIN_TAG_HEADER(FILEID,ERROR,*9999)
          CALL IO_VECTOR('WRITE','DISK',FILEID,NOT(1,1,nr,nx),
     '      GRR,NOM,'BINARY','ALL_INOUT',
     '      ERROR,*9999)
        ENDIF

C*** Close binary file

        CALL BINCLOSEFILE(FILEID,ERROR,*9999)
C
      ENDIF

      CALL EXITS('WRITE_SOL_MATRIX')
      RETURN
 9999 CALL ERRORS('WRITE_SOL_MATRIX',ERROR)
      CALL EXITS('WRITE_SOL_MATRIX')
      RETURN 1
      END



