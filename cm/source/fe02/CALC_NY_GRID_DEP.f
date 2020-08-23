      SUBROUTINE CALC_NY_GRID_DEP(NHQ,NQNY,nr,nx,NYNQ,NYQNR,
     '  ERROR,*)

C#### Subroutine: CALC_NY_GRID_DEP
C###  Description:
C###    CALC_NY_GRID_DEP calculates the mapping array NYNQ and
C###    their inverse NQNY for dependent variables.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NHQ(NRM),NQNY(2,NYQM,0:NRCM),nr,nx,
     '  NYNQ(NHM,NQM,0:NRCM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,nc,nh,nq,nrc,nrr,ny,ny_start(0:2)
      CHARACTER CHAR*1
!     External function
      INTEGER IDIGITS

      CALL ENTERS('CALC_NY_GRID_DEP',*9999)

      nc=1
C NOTE: no nc indicies for NYNQ and NQNY nc assumed =1 for all cases
C where ityp4=1 NPS 1/11/96

      CALL ASSERT(NRCM.GE.2,'>>Increase NRCM to be >= 2',
     '  ERROR,*9999)

      DO nrc=0,2
        ny_start(nrc)=0
      ENDDO

C finding starting ny for current region
      DO nrc=0,2
        DO nrr=1,nr-1
          DO nq=NQR(1,nrr),NQR(2,nrr)
            DO nh=1,NHQ(nrr)
              ny_start(nrc)=ny_start(nrc)+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C initalising NYNQ and NQNY for region above and including current

      DO nrc=0,2
        ny=ny_start(nrc)
        DO nrr=nr,NRT
          DO nq=NQR(1,nrr),NQR(2,nrr)
            DO nh=1,NHQ(nrr)
              ny=ny+1
              IF (ny.LE.NYQM) THEN
                DO i=1,2
                  NQNY(i,ny,nrc)=0
                ENDDO
              ENDIF
              NYNQ(nh,nq,nrc)=0
            ENDDO !nh
          ENDDO !nq
          NYQNR(0,nrc,nc,nrr)=0
        ENDDO !nrr
      ENDDO !nrc

C Set up mapping arrays for current region

      DO nrc=0,2
        ny=ny_start(nrc)
        DO nq=NQR(1,nr),NQR(2,nr)
          DO nh=1,NHQ(nr)
            ny=ny+1
            NYQNR(0,nrc,nc,nr)=NYQNR(0,nrc,nc,nr)+1
            IF(nrc.NE.0) THEN
              IF(ny.GT.NYQT(nrc,nc,nx)) NYQT(nrc,nc,nx)=ny
            ENDIF
            IF(NYQNR(0,nrc,nc,nr).LE.NYQM)  THEN
              NYQNR(NYQNR(0,nrc,nc,nr),nrc,nc,nr)=ny
              NQNY(1,ny,nrc)=nh
              NQNY(2,ny,nrc)=nq
            ENDIF
            NYNQ(nh,nq,nrc)=ny
          ENDDO
        ENDDO
      ENDDO
      DO nrc=0,2
        IF(NYQNR(0,nrc,nc,nr).GT.NYQM) THEN
          WRITE(CHAR,'(I1)') IDIGITS(NYQNR(0,nrc,nc,nr))
          WRITE(ERROR,'(''>>Increase NYQM to '',I'//CHAR//')')
     '      NYQNR(0,nrc,nc,nr)
          GO TO 9999 
        ENDIF
C        CALL ASSERT(NYQNR(0,nrc,nc,nr).LE.NYQM,'>>Increase NYQM',
C     '    ERROR,*9999)
      ENDDO !nrc

      CALL EXITS('CALC_NY_GRID_DEP')
      RETURN
 9999 CALL ERRORS('CALC_NY_GRID_DEP',ERROR)
      CALL EXITS('CALC_NY_GRID_DEP')
      RETURN 1
      END


