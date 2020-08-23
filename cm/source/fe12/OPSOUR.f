      SUBROUTINE OPSOUR(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,NDIPOLES,
     '  nr,DIPOLE_CEN,DIPOLE_DIR,ERROR,*)

C#### Subroutine: OPSOUR
C###  Description:
C###    OPSOUR outputs source parameters for region nr.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),NDIPOLES(NRM),nr
      REAL*8 DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nj,nt

      CALL ENTERS('OPSOUR',*9999)
      WRITE(OP_STRING,
     '  '(/'' The number of dipole sources in region '',I1,'
     '  //''' is '',I2)') nr,NDIPOLES(nr)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO n=1,NDIPOLES(nr)
        WRITE(OP_STRING,'(/'' Dipole '',I2,'': '')') n
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(DIPOLE_CEN_NTIME(n,nr).GT.0) THEN
          WRITE(OP_STRING,'(''   The dipole center is moving'')')
        ELSE
          WRITE(OP_STRING,'(''   The dipole center is static'')')
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(DIPOLE_DIR_NTIME(n,nr).GT.0) THEN
          WRITE(OP_STRING,'(''   The dipole vector is moving'')')
        ELSE
          WRITE(OP_STRING,'(''   The dipole vector is static'')')
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(DIPOLE_CEN_NTIME(n,nr).GT.0) THEN
          WRITE(OP_STRING,'(''   There are '',I3,'' time points for '
     '      //'the dipole centre'')') DIPOLE_CEN_NTIME(n,nr)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''   Initial dipole centre : '',3D12.4)')
     '      (DIPOLE_CEN(nj,0,n,nr),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nt=1,DIPOLE_CEN_NTIME(n,nr)
            WRITE(OP_STRING,'(''   Time point '',I3,'', Time='',D12.4,'
     '        //''' s, Dipole centre : '',3D12.4)') nt,
     '        DIPOLE_CEN(4,nt,n,nr),(DIPOLE_CEN(nj,nt,n,nr),nj=1,NJT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nt
        ELSE
          WRITE(OP_STRING,'(''   Dipole centre: '',3D12.4)')
     '      (DIPOLE_CEN(nj,0,n,nr),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(DIPOLE_DIR_NTIME(n,nr).GT.0) THEN
          WRITE(OP_STRING,'(''   There are '',I3,'' time points for '
     '      //'the dipole vector'')') DIPOLE_DIR_NTIME(n,nr)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''   Initial dipole vector : '',3D12.4)')
     '      (DIPOLE_DIR(nj,0,n,nr),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nt=1,DIPOLE_DIR_NTIME(n,nr)
            WRITE(OP_STRING,'(''   Time point '',I3,'', Time='',D12.4,'
     '        //''' s, Dipole vector : '',3D12.4)') nt,
     '        DIPOLE_DIR(4,nt,n,nr),(DIPOLE_DIR(nj,nt,n,nr),nj=1,NJT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nt
        ELSE
          WRITE(OP_STRING,'(''   Dipole vector: '',3D12.4)')
     '      (DIPOLE_DIR(nj,0,n,nr),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !n

      CALL EXITS('OPSOUR')

      RETURN
 9999 CALL ERRORS('OPSOUR',ERROR)
      CALL EXITS('OPSOUR')
      RETURN 1
      END


