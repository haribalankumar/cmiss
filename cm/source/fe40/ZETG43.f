      SUBROUTINE ZETG43(NBH,ng,AXIAL_STRAIN,CG,DX2,EAU,EAUV,EIV,
     '  PG,XE,ZE,ERROR,*)

C#### Subroutine: ZETG43
C###  Description:
C###    ZETG43 is for beam elements only. Tension is in kN.
C###    Note: geometric basis and axial displacement basis is assumed
C###    linear.  Transverse displacement basis must be cubic Hermite.

C**** CG(1) = Young's modulus (GPa)
C**** CG(2) = Poisson's ratio
C**** CG(3) = Depth of beam (m)
C**** CG(4) = Width of beam (m)
C**** CG(5) = Wall thickness of beam (m) (if tubular x-section only)
C**** EA    is Young's modulus * Xsection area in kN
C**** EI    is Young's modulus * 2nd moment of area in kN.m^2
C**** DUDX  is dU/dX, DVDX is dV/dx, D2VDX2 is d^2/dV^2
C**** EAU   is EA*(dU/dX+.5*(dV/dX)^2) (axial force in kN)
C**** EAUV  is EAU*dV/dX
C**** EIV   is EI*d^2V/dX^2 (bending moment in kNm)
C**** Note: Put DVDX=0 for linear beam theory

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),ng
      REAL*8 AXIAL_STRAIN,CG(NMM),DX2,EAU,EAUV,EIV,
     '  PG(NSM,NUM,NGM,NBM),XE(NSM,NJM),ZE(NSM,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nj,ns,nx
      REAL*8 DVDX,D2VDX2,EA,EI,AREA,DEPTH,SX,SZ,THICK,WIDTH
C     REAL*8 POIS,YMOD
      CALL ENTERS('ZETG43',*9999)

      nx=1 !temporary
C      YMOD =CG(1) !is Young's modulus (GPa)
C      POIS =CG(2) !is Poisson's ratio
      DEPTH=CG(3) !is depth of beam
      WIDTH=CG(4) !is width of beam

      IF(KTYP45.EQ.1) THEN      !rectangular solid
        AREA=WIDTH*DEPTH
C        RIZZ=WIDTH*DEPTH**3/12.d0 !is 2nd mom of area about z-axis
C        RIYY=DEPTH*WIDTH**3/12.d0 !is 2nd mom of area about y-axis
      ELSE IF(KTYP45.EQ.2) THEN !rectangular tube
        THICK=CG(5)              !is wall thickness
        AREA=WIDTH*DEPTH-(WIDTH-2.0d0*THICK)*(DEPTH-2.0d0*THICK)
C        RIZZ=WIDTH*DEPTH**3/12.d0 !is 2nd mom of area about z-axis
C     '    -(WIDTH-2.0d0*THICK)*(DEPTH-2.0d0*THICK)**3/12.d0
C        RIYY=DEPTH*WIDTH**3/12.d0 !is 2nd mom of area about y-axis
C     '    -(DEPTH-2.0d0*THICK)*(WIDTH-2.0d0*THICK)**3/12.d0
      ELSE IF(KTYP45.EQ.3) THEN !ellipsoidal solid
        AREA=PI*WIDTH*DEPTH/4.d0
C        RIZZ=PI*WIDTH*DEPTH**3/64.d0 !is 2nd mom about z-axis
C        RIYY=PI*DEPTH*WIDTH**3/64.d0 !is 2nd mom about y-axis
      ELSE IF(KTYP45.EQ.4) THEN !ellipsoidal tube
        THICK=CG(5)                 !is wall thickness
        AREA=PI*(WIDTH*DEPTH-(WIDTH-2.0d0*THICK)*
     '    (DEPTH-2.0d0*THICK))/4.d0
C        RIZZ=PI*WIDTH*DEPTH**3/64.d0 !is 2nd mom about z-axis
C     '    -PI*(WIDTH-2.0d0*THICK)*(DEPTH-2.0d0*THICK)**3/64.d0
C        RIYY=PI*DEPTH*WIDTH**3/64.d0 !is 2nd mom about y-axis
C     '    -PI*(DEPTH-2.0d0*THICK)*(WIDTH-2.0d0*THICK)**3/64.d0
      ENDIF

C      SMOD=YMOD/(2.d0*(1.d0+POIS)) !is shear modulus (GPa)
C GBS 19/09/95 - not used
C      EAIZZ=YMOD*AREA*RIZZ       !is bending rigidity (z-axis)
C      GAKZZ=SMOD*AREA*RIZZ       !is bending rigidity (z-axis)
C      GAKYY=SMOD*AREA*RIYY       !is bending rigidity (y-axis)

      EA=CG(1)*AREA*1.d2 !(kN)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(''  EA='',E11.3,''(kN)'')') EA
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C     nb=NBH(Nh_LOC(1,nx))
C     DUDX=0.
C     DO ns=1,NST(nb)
C       DUDX=DUDX+PG(ns,2,ng,nb)*DX1*(ZE(ns,1)-XE(ns,1))
C     ENDDO
C     IF(DOP) THEN
C       WRITE(OP_STRING,'('' DUDX='',E11.3)') DUDX
C       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C     ENDIF

C     nb=NBH(NH_LOC(2,nx))
C     DVDX=0.
C     DO ns=1,NST(nb)
C       DVDX=DVDX+PG(ns,2,ng,nb)*DX1*(ZE(ns,2)-XE(ns,2))
C     ENDDO
C     IF(DOP) THEN
C       WRITE(OP_STRING,'('' DVDX='',E11.3)') DVDX
C       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C     ENDIF

C KAT 2001-11-23 EAU is used here before set!
C     This must be wrong so let the user know if they are assuming that
C     this works.
      CALL ASSERT(.FALSE.,'Implementation not complete',ERROR,*9999)

C MLB 16/4/97 Put dvdx to zero as in standard beam problem because
C EAUV has not been defined.
      DVDX=0.0d0
      EAUV=EAU*DVDX
C

C     EAU=EA*(DUDX+0.5d0*DVDX**2)
C     EAUV=EAU*DVDX

      !add displacement in m onto initial coords
      SX=0.0d0
      SZ=0.0d0
      DO nj=1,NJT
        SX=SX+ (XE(2,nj)-XE(1,nj))**2
        SZ=SZ+((XE(2,nj)+ZE(3,nj)*1.D-2)-(XE(1,nj)+ZE(1,nj)*1.D-2))**2
      ENDDO
      SX=DSQRT(SX)
      SZ=DSQRT(SZ)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) ' SX=',SX,' SZ=',SZ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      AXIAL_STRAIN=SZ/SX-1.d0
      EAU=EA*AXIAL_STRAIN !axial force

      IF(NVE(3).GT.1) THEN !transverse deflections defined
        EI=CG(1)*CG(3)*1.D-2 !(kN.m^2)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(''  EI='',E11.3,''(kN.m^2)'')') EI
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        nb=NBH(NH_LOC(2,nx))
        D2VDX2=0.0d0
        DO ns=1,NST(nb)
          D2VDX2=D2VDX2+PG(ns,3,ng,nb)*DX2*(ZE(ns,2)-XE(ns,2))
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' D2VDX2='',E11.3)') D2VDX2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        EIV=EI*D2VDX2 !bending moment

      ELSE
        EIV=0.d0
      ENDIF

      CALL EXITS('ZETG43')
      RETURN
 9999 CALL ERRORS('ZETG43',ERROR)
      CALL EXITS('ZETG43')
      RETURN 1
      END


