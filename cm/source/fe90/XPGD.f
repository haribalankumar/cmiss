      SUBROUTINE XPGD(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,NBH,NDIPOLES,
     '  NENP,NKH,NP_INTERFACE,np,nr,NW,nx,NYNP,CE,DIPOLE_CEN,DIPOLE_DIR,
     '  GD,TIME,XG,XP,XPFP,ERROR,*)

C#### Subroutine: XPGD
C###  Description:
C###    XPGD calculates the vector GD which contains domain integrals
C###    or source terms.

C**** Currently (9-9-94) only implemented for dipole sources.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'anal00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),NBH(NHM,NCM,NEM),NDIPOLES(NRM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NKH(NHM,NPM,NCM),NP_INTERFACE(0:NPM,0:3),np,nr,NW(NEM,3),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),GD(NZ_GD_M),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),XPFP(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IMAT,n_dipole,nd,NDTOT,nj,nk,ny1
      REAL*8 CENTRE(3),DIRECTION(3),DSDX(3,3),DXDS(3,3),HYPGREEN,
     '  P(3),SIGMA,SUM,SUMR,TIME,XDR(3),XNO(3,4)
      CHARACTER FORMAT*500

      CALL ENTERS('XPGD',*9999)


      ny1=NYNP(1,1,NH_LOC(1,nx),np,1,1,nr)
      IF(DABS(CE(1,1)).LE.RDELTA) THEN
        IMAT=0
        SIGMA=1.0d0 !using standard laplace
      ELSE
        IMAT=1
        SIGMA=CE(1,1)
      ENDIF
      IF(KTYP92.LT.3) THEN !Normal BIE used

        DO n_dipole=1,NDIPOLES(nr)
          CALL GETDIPOLE(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,n_dipole,
     '      nr,CENTRE,DIPOLE_CEN,DIPOLE_DIR,DIRECTION,TIME,
     '      ERROR,*9999)

          SUM=0.0d0
          SUMR=0.0d0
          DO nj=1,NJT
            XDR(nj)=CENTRE(nj)-XPFP(nj) !r is from source to field point
            SUM=SUM+XDR(nj)*DIRECTION(nj) !r.p
            SUMR=SUMR+XDR(nj)*XDR(nj) !r*r
          ENDDO
          SUMR=DSQRT(SUMR)

C cpb 19/1/95 Adding 2D dipole
          IF(NJT.EQ.2) THEN
C           -r.p/(2*pi*sigma*r)*(1/r)
C           [for dot product r needs to be of unit length]
            SUM=-SUM/(2.0d0*PI*SIGMA*SUMR*SUMR)
          ELSE IF(NJT.EQ.3) THEN
C           r.p/(4*pi*sigma*r^2)*(1/r)
C           [for dot product r needs to be of unit length]
            SUM=SUM/(4.0d0*PI*SIGMA*SUMR*SUMR*SUMR)
          ENDIF

C cpb 19/1/96 Changing from -SUM to +SUM as GD is subtracted from GRR
          GD(ny1)=GD(ny1)+SUM
        ENDDO !n_dipole
      ENDIF !normal BEM

      IF(HYP) THEN
        IF(NKH(NH_LOC(1,nx),np,1).GT.1) THEN !More than one derivative at np
          CALL XPXNO(NBH,NENP,NP_INTERFACE,np,nr,NW,nx,DSDX,
     '      DXDS,XG,XNO,XP,ERROR,*9999)
          NDTOT=MAX(NKH(NH_LOC(1,nx),np,1)-1-KTYP93(1,nr),1)
          DO nd=1,NDTOT
            nk=nd+1
            ny1=NYNP(nk,1,NH_LOC(1,nx),np,1,1,nr)
            DO n_dipole=1,NDIPOLES(nr)
              CALL GETDIPOLE(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,n_dipole,
     '          nr,CENTRE,DIPOLE_CEN,DIPOLE_DIR,DIRECTION,TIME,
     '          ERROR,*9999)
              SUMR=0.0d0
              DO nj=1,NJT
                XDR(nj)=CENTRE(nj)-XPFP(nj)
                SUMR=SUMR+XDR(nj)**2 !r*r
                P(nj)=DIRECTION(nj)
              ENDDO !nj
              SUMR=DSQRT(SUMR) !r
              IF(NJT.EQ.2) THEN
                GD(ny1)=GD(ny1)+HYPGREEN(IGREN(nr),IMAT,CE(1,1),SUMR,
     '            1.0d0,P,XNO(1,nd),XDR)
              ELSE IF(NJT.EQ.3) THEN
                GD(ny1)=GD(ny1)-HYPGREEN(IGREN(nr),IMAT,CE(1,1),SUMR,
     '            1.0d0,P,XNO(1,nd),XDR)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XPGD_1)
        WRITE(OP_STRING,'('' GD vector'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Node '',I5)')np
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        FORMAT='('' GD('',I5,'')='',D12.4)'
        WRITE(OP_STRING,FORMAT)ny1,GD(ny1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(XPGD_1)
      ENDIF

      CALL EXITS('XPGD')
      RETURN
 9999 CALL ERRORS('XPGD',ERROR)
      CALL EXITS('XPGD')
      RETURN 1
      END


