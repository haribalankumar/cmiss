      SUBROUTINE CGC11(nr,nw,nx,CC,CG,DENS,ERROR,*)

C#### Subroutine: CGC11
C###  Description:
C###    <HTML>
C###    CGC11 evaluates coefficents for plane stress and plane strain.
C###    CC(i,j) is the 3*3 matrix of elastic material coefficients.
C###    either:
C###    <PRE>
C###       1 isotropic,
C###       2 transversely isotropic wrt Xi(1),
C###       3 transversely isotropic wrt Xi(2),
C###    or:
C###       4 orthotropic.
C###    Note: 2d_stress_vector = CC * 2d_strain_vector.
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nr,nw,nx
      REAL*8 CC(3,*),CG(NMM),DENS
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 POIS,POIS1,POIS2,SMOD1,YMOD,YMOD1,YMOD2,YY

      CALL ENTERS('CGC11',*9999)

      IF(nw.EQ.11.OR.nw.EQ.12) THEN
        IF(IMT(nw).EQ.1) THEN !isotropic
          YMOD=CG(1) !is Young's modulus (kPa)
          POIS=CG(2) !is Poisson's ratio
C cpb 15/4/96 For modal analysis if density is in Mg/m^3 and Youngs
C modulus is in GPa then the frequency is in kHz
          DENS=CG(ILT(nw,nr,nx))*1.0d-3 !is density in Mg/m^3
          IF(nw.EQ.11) THEN
            YY=YMOD/(1.0d0-POIS*POIS)
            CC(1,2)=YY*POIS
            CC(2,1)=CC(1,2)
            CC(3,3)=YY*(1.0d0-POIS)/2.0d0
          ELSE IF(nw.EQ.12) THEN
            YY=YMOD*(1.0d0-POIS)/(1.0d0+POIS)/(1.0d0-2.0d0*POIS)
            CC(1,2)=YY*POIS/(1.d0-POIS)
            CC(2,1)=CC(1,2)
            CC(3,3)=YY*(1.0d0-2.0d0*POIS)/(2.0d0*(1.0d0-POIS))
          ENDIF
          CC(1,1)=YY
          CC(1,3)=0.0d0
          CC(2,2)=YY
          CC(2,3)=0.0d0
          CC(3,1)=0.0d0
          CC(3,2)=0.0d0
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' YMOD='',E11.3,'' POIS='',E11.3,'' YY='','//'E11.3)')
     '        YMOD,POIS,YY
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

        ELSE IF(IMT(nw).EQ.2) THEN !transv isotropy (fibre in 1 dir.n)
          YMOD1=CG(1) !is Young's modulus E1 (GPa)
          YMOD2=CG(2) !is Young's modulus E2 (GPa)
          SMOD1=CG(3) !is shear modulus G1   (GPa)
          POIS1=CG(4) !is Poisson's ratio nu1
          POIS2=CG(5) !is Poisson's ratio nu2
C cpb 15/4/96 For modal analysis if density is in Mg/m^3 and Youngs
C modulus is in GPa then the frequency is in kHz
          DENS=CG(ILT(nw,nr,nx))*1.0d-3 !is density in Mg/m^3
          IF(nw.EQ.11) THEN
            YY=1.0d0/(1.0d0-POIS1*POIS2)
            CC(1,1)=YMOD1*YY
            CC(1,2)=YMOD1*POIS2*YY
            CC(1,3)=0.0d0
c           CC(2,1)=YMOD2*POIS1*YY
            CC(2,1)=CC(1,2)
            CC(2,2)=YMOD2*YY
            CC(2,3)=0.0d0
            CC(3,1)=0.0d0
            CC(3,2)=0.0d0
            CC(3,3)=SMOD1
          ELSE IF(nw.EQ.12) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ENDIF

        ELSE IF(IMT(nw).EQ.3) THEN !transv. isotropic (fibre in 2 dir)
          YMOD1=CG(1) !is Young's modulus E1 (GPa)
          YMOD2=CG(2) !is Young's modulus E2 (GPa)
          SMOD1=CG(3) !is shear modulus G1   (GPa)
          POIS1=CG(4) !is Poisson's ratio nu1
          POIS2=CG(5) !is Poisson's ratio nu2
C cpb 15/4/96 For modal analysis if density is in Mg/m^3 and Youngs
C modulus is in GPa then the frequency is in kHz
          DENS=CG(ILT(nw,nr,nx))*1.0d-3 !is density in Mg/m^3
          IF(nw.EQ.11) THEN
            YY=1.0d0/(1.0d0-POIS1*POIS2)
            CC(1,1)=YMOD1*YY
            CC(1,2)=YMOD1*POIS2*YY
            CC(1,3)=0.0d0
            CC(2,1)=YMOD2*POIS1*YY
            CC(2,2)=YMOD2*YY
            CC(2,3)=0.0d0
            CC(3,1)=0.0d0
            CC(3,2)=0.0d0
            CC(3,3)=SMOD1
          ELSE IF(nw.EQ.12) THEN
            ERROR='>>Not implemented'
            GOTO 9999
          ENDIF

        ELSE IF(IMT(nw).EQ.4) THEN !orthotropic (fibre in 1 dir.n)
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF
      ELSE
        ERROR='>>Invalid NW number'
        GOTO 9999
      ENDIF

      CALL EXITS('CGC11')
      RETURN
 9999 CALL ERRORS('CGC11',ERROR)
      CALL EXITS('CGC11')
      RETURN 1
      END


