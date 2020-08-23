      SUBROUTINE SGDIPO(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,ISEG,ISDIPO,
     '  ISDIPA,iw,n,nr,SCALE,CSEG,DIPOLE_CEN,DIPOLE_DIR,PATH,
     '  VECT,ERROR,*)

C#### Subroutine: SGDIPO
C###  Description:
C###    SGDIPO creates new dipole segment ISDIPO(iw,n,nr) and dipole
C###    path segment ISDIPA(iw,n,nr) on window iw for dipole n in
C###    region nr.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),ISEG(*),ISDIPO,ISDIPA,iw,n,nr
      REAL*8 DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),SCALE
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL PATH,VECT
!     Local Variables
      INTEGER INDEXDIPOLE,INDEXDIPOLEPATH,INDEX_POLYLINE,
     '  INDEX_OLD,nj,nt
      REAL*8 Z(3,0:NDIPTIMM),CENTER(3),DIRECTION(3)
      CHARACTER CFROMI*2,CLABEL*15

      CALL ENTERS('SGDIPO',*9999)

      IF(VECT) THEN

        INDEXDIPOLE=INDEX_POLYLINE(0,'SOLID','WIDTH2','BLACK')

        CLABEL='DIPOLE'//CFROMI(nr,'(I2)')//CFROMI(n,'(I2)')
        CALL OPEN_SEGMENT(ISDIPO,ISEG,iw,CLABEL,INDEXDIPOLE,INDEX_OLD,
     '    n,1,CSEG,ERROR,*9999)

        DO nj=1,NJT
          CENTER(nj)=DIPOLE_CEN(nj,0,n,nr)
          DIRECTION(nj)=DIPOLE_DIR(nj,0,n,nr)
        ENDDO !nj

        CALL VECTOR(INDEXDIPOLE,iw,CENTER,DIRECTION,SCALE,ERROR,*9999)

        CALL CLOSE_SEGMENT(ISDIPO,iw,ERROR,*9999)

      ELSE IF(PATH) THEN

        INDEXDIPOLEPATH=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')

        CLABEL='DIPOLEPATH'//CFROMI(nr,'(I2)')//CFROMI(n,'(I2)')
        CALL OPEN_SEGMENT(ISDIPA,ISEG,iw,CLABEL,INDEXDIPOLEPATH,
     '    INDEX_OLD,n,1,CSEG,ERROR,*9999)

        IF(DIPOLE_CEN_NTIME(n,nr).GT.0) THEN
          DO nt=0,DIPOLE_CEN_NTIME(n,nr)
            DO nj=1,NJT
              Z(nj,nt)=DIPOLE_CEN(nj,nt,n,nr)
            ENDDO !nj
          ENDDO !nt
          CALL POLYLINE(INDEXDIPOLEPATH,iw,DIPOLE_CEN_NTIME(n,nr)+1,Z,
     '      ERROR,*9999)
        ENDIF
        IF(DIPOLE_DIR_NTIME(n,nr).GT.0) THEN
          IF(DIPOLE_CEN_NTIME(n,nr).GT.0) THEN
            DO nt=0,DIPOLE_DIR_NTIME(n,nr)
              DO nj=1,NJT
                Z(nj,nt)=DIPOLE_CEN(nj,nt,n,nr)+DIPOLE_DIR(nj,nt,n,nr)*
     '            SCALE
              ENDDO !nj
            ENDDO !nt
          ELSE
            DO nt=0,DIPOLE_DIR_NTIME(n,nr)
              DO nj=1,NJT
                Z(nj,nt)=DIPOLE_CEN(nj,0,n,nr)+DIPOLE_DIR(nj,nt,n,nr)*
     '            SCALE
              ENDDO !nj
            ENDDO !nt
          ENDIF
          CALL POLYLINE(INDEXDIPOLEPATH,iw,DIPOLE_DIR_NTIME(n,nr)+1,Z,
     '      ERROR,*9999)
        ENDIF

        CALL CLOSE_SEGMENT(ISDIPA,iw,ERROR,*9999)

      ENDIF

      CALL EXITS('SGDIPO')
      RETURN
 9999 CALL ERRORS('SGDIPO',ERROR)
      CALL EXITS('SGDIPO')
      RETURN 1
      END


