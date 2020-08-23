      SUBROUTINE IPVALU(VENTILATION,ERROR,*)

C#### Subroutine: IPVALU
C###  Description:
C###    IPVALU defines parameter values in an input file. Currently set 
C###    up for pulmonary problems only. Specific parameters depend on the 
C###    problem type. This is for non-spatially varying parameters with
C###    variable names that are contained in common block files. This 
C###    enables parameter values to be set in an *.ipvalu file rather
C###    than being hard-coded or entered as material parameters.
C###    Uses ITYP nr=1 and nx=1 -- needs to be extended for multiple
C###    regions and classes.
C###  Created by AJS Nov 2010

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'


!     Parameter List
!       INTEGER 
!       REAL*8 
      LOGICAL VENTILATION
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES,nr,nx
      LOGICAL FILEIP

      CALL ENTERS('IPVALU',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      nr=1
      nx=1
      RMAX=1.7D+37
      RMIN=-RMAX


      IF(VENTILATION)THEN !Ventilation model (currently no ITYP defined for this)
C Gas properties (default = air)
        FORMAT='($,'' Gas viscosity (Pa.s)                           '//
     &    '   [  1.8e-5]: '',D11.4)'
        RDEFLT(1)=1.8d-5 !Pa.s
        IF(IOTYPE.EQ.3) RDATA(1)=GAS_VISCOSITY
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GAS_VISCOSITY=RDATA(1)
        FORMAT='($,'' Gas density (g/mm3)                            '//
     &    '   [1.146e-6]: '',D11.4)'
        RDEFLT(1)=1.146d-6
        IF(IOTYPE.EQ.3) RDATA(1)=GAS_DENSITY
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GAS_DENSITY=RDATA(1)

C Strain energy density function W=c/2*exp(a*J1^2 + b*J2) [Tawhai, 2009, J Appl Physiol v107]
        FORMAT='($,'' Strain energy density function coefficient a   '//
     &    '   [   0.433]: '',D11.4)'
        RDEFLT(1)=0.433d0
        IF(IOTYPE.EQ.3) RDATA(1)=SEDF_COEFFS(1)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SEDF_COEFFS(1)=RDATA(1)
        FORMAT='($,'' Strain energy density function coefficient b   '//
     &    '   [  -0.611]: '',D11.4)'
        RDEFLT(1)=-0.611d0
        IF(IOTYPE.EQ.3) RDATA(1)=SEDF_COEFFS(2)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SEDF_COEFFS(2)=RDATA(1)
        FORMAT='($,'' Strain energy density function coefficient c '//
     &    '(Pa) [    2500]: '',D11.4)'
        RDEFLT(1)=2500.d0
        IF(IOTYPE.EQ.3) RDATA(1)=SEDF_COEFFS(3)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SEDF_COEFFS(3)=RDATA(1)

      ELSEIF((ITYP3(nr,nx).EQ.4).AND.
     &  (COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6)) THEN !Blood flow problem
C Parameters for sheet flow capillary model         
C...  DEFINING THE PARAMETERS THAT INFLUENCE SHEET HEIGHT
        FORMAT='($,'' Unstrained sheet height (m)                   '//
     &    '        [ 3.5d-6]: '',D11.4)'
        RDEFLT(1)=3.5d-6 !m
        IF(IOTYPE.EQ.3) RDATA(1)=H0
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) H0=RDATA(1)
        ! Friction factor
        FORMAT='($,'' Friction factor (K)                           '//
     &    '        [ 12.0d0]: '',D11.4)'
        RDEFLT(1)=12.0d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=K_cap
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) K_cap=RDATA(1)
        FORMAT='($,'' Friction factor (F)                            '//
     &    '       [ 1.8d0]: '',D11.4)'
        RDEFLT(1)=1.8d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=F_cap
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) F_cap=RDATA(1)
        !Zone 2 constants
        FORMAT='($,'' Zone 2 constants (F)                           '//
     &    '       [ 0.104d0]: '',D11.4)'
        RDEFLT(1)=0.104d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=F_sheet
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) F_sheet=RDATA(1)
        FORMAT='($,'' Zone 2 constants (sigma, Pa)                   '//
     &    '       [ 4.45d0*98.06d0]: '',D11.4)'
        RDEFLT(1)=4.45d0*98.06d0 !Pa
        IF(IOTYPE.EQ.3) RDATA(1)=sigma_cap
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) sigma_cap=RDATA(1)
        !Apparent viscosity
        FORMAT='($,'' Apparent viscosity (Pa.s)                   '//
     &    '          [ 1.92d-3]: '',D11.4)'
        RDEFLT(1)=1.92d-3 !Pa.s or cP!? 
        IF(IOTYPE.EQ.3) RDATA(1)=mu_c
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) mu_c=RDATA(1)
        !Recruitment parameters
        FORMAT='($,'' Recruitment parameters (F)                   '//
     &    '         [ 0.6463d0]: '',D11.4)'
        RDEFLT(1)=0.6463d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=F_rec
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) F_rec=RDATA(1)
        FORMAT='($,'' Recruitment parameters (sigma)                 '//
     &    '       [ 2230.d0]: '',D11.4)'
        RDEFLT(1)=2230.d0 !Pa
        IF(IOTYPE.EQ.3) RDATA(1)=sigma_rec
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) sigma_rec=RDATA(1)
        ! DEFINING SHEET AREA AND PATHLENGTH
        FORMAT='($,'' Path length from arteriole to venule (m)       '//
     &    '       [ 1188.d-6]: '',D11.4)'
        RDEFLT(1)=1188.d-6 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=L_c
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) L_c=RDATA(1)
        !Upper and lower bounds for sheet height
        FORMAT='($,'' Lower pressure bound for sheet (cmH2O)       '//
     &    '         [ 0.d0]: '',D11.4)'
        RDEFLT(1)=0.d0 !cmH2O
        IF(IOTYPE.EQ.3) RDATA(1)=Plb_c
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Plb_c=RDATA(1)
        Plb_c=Plb_c*98.07d0 !Converting cmH2O -> Pa
        FORMAT='($,'' Upper pressure bound for sheet (cmH2O)       '//
     &    '         [ 32.d0]: '',D11.4)'
        RDEFLT(1)=32.d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=Pub_c
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Pub_c=RDATA(1)
        Pub_c=Pub_c*98.07d0 !Converting cmH2O -> Pa
        FORMAT='($,'' Upper pressure bound for intra-acinar vessels'//
     &    ' (cmH2O) [ 32.d0]: '',D11.4)'
        RDEFLT(1)=32.d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=Pub_a_v
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Pub_a_v=RDATA(1)
        Pub_a_v=Pub_a_v*98.06d0 !Converting cmH2O->Pa
        !Total capillary surface area in the lung
        FORMAT='($,'' Total capillary surface area (m^2)        '//
     &    '            [ 6.8d1]: '',D11.4)'
       ! - total_cap_area=0.5d0*126.d0/31800.d0 !Area from Gher 1978,approx TLC
        RDEFLT(1)=6.3d1 !m^2
        IF(IOTYPE.EQ.3) RDATA(1)=total_cap_area
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) total_cap_area=RDATA(1)
        IF(NTERMINAL.GT.0)THEN
           total_cap_area=total_cap_area/NTERMINAL
        ELSE       
          write(*,*) 'WARNING: NTERMINAL is zero - assuming 32000 acini'
          write(*,*) 'remember to do fem evaluate order!'
          total_cap_area=total_cap_area/32000d0
        ENDIF
        !Minimum arteriole/venule lengths and radii
        FORMAT='($,'' Terminal arteriole length (m)              '//
     &    '           [ 0.130d-3]: '',D11.4)'
        RDEFLT(1)=0.130d-3 !m
        IF(IOTYPE.EQ.3) RDATA(1)=L_art_terminal
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) L_art_terminal=RDATA(1)
           FORMAT='($,'' Terminal venule length (m)               '//
     &    '             [ 0.130d-3]: '',D11.4)'
        RDEFLT(1)=0.130d-3 !m
        IF(IOTYPE.EQ.3) RDATA(1)=L_vein_terminal
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) L_vein_terminal=RDATA(1)
          FORMAT='($,'' Terminal arteriole radius (m)            '//
     &    '             [ 0.01d-3]: '',D11.4)'
        RDEFLT(1)=0.01d-3 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=R_art_terminal
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) R_art_terminal=RDATA(1)
         FORMAT='($,'' Terminal venule radius (m)                '//
     &    '            [ 0.009d-3]: '',D11.4)'
        RDEFLT(1)=0.009d-3 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=R_vein_terminal
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) R_vein_terminal=RDATA(1)
        FORMAT='($,'' Initial acinar model resistance                '//
     &    '       [ 5000.d0]: '',D11.4)'
        RDEFLT(1)=5000.d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=INITIAL_LPM
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) INITIAL_LPM=RDATA(1)
        FORMAT='($,'' Upper pressure bound for extra-acinar vessels'//
     &    ' (cmH2O) [ 32.d0]: '',D11.4)'
        RDEFLT(1)=32.d0 !no units
        IF(IOTYPE.EQ.3) RDATA(1)=Ptm_max
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Ptm_max=RDATA(1)
        Ptm_max=Ptm_max*98.06d0 !Converting cmH2O->Pa
        FORMAT='($,'' Number of symmetric generations              '//
     &    '         [     9]: '',I3)'
        IDEFLT(1)=9 !no units
        IF(IOTYPE.EQ.3) IDATA(1)=num_symm_gen
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) num_symm_gen=IDATA(1)

C Pulmonary advection-diffusion-exchange problems
      ELSEIF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.11.AND.
     &  ITYP3(nr,nx).EQ.1.AND.ITYP7(nr,nx).GT.1)THEN
        FORMAT='($,'' Mixed venous PO2 (mmHg)                       '//
     &    '   [      40]: '',D11.4)'
        RDEFLT(1)=40.0d0 !mmHg
        IF(IOTYPE.EQ.3) RDATA(1)=MIXED_VENOUS_PO2
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MIXED_VENOUS_PO2=RDATA(1)

        FORMAT='($,'' Air-blood barrier harmonic mean thickness (um)'//
     &    '   [    1.11]: '',D11.4)'
        RDEFLT(1)=1.11d0 !um
        IF(IOTYPE.EQ.3) RDATA(1)=BARRIER_THICKNESS
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) BARRIER_THICKNESS=RDATA(1)
        BARRIER_THICKNESS=BARRIER_THICKNESS*1.0d-6 ! convert um to m

        FORMAT='($,'' Air-blood surface area per acinus (cm2)       '//
     &    '   [      30]: '',D11.4)'
        RDEFLT(1)=30.0d0 !m2
        IF(IOTYPE.EQ.3) RDATA(1)=ACINUS_SURFACE_AREA
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ACINUS_SURFACE_AREA=RDATA(1)
        ACINUS_SURFACE_AREA=ACINUS_SURFACE_AREA*1.0d-4 ! convert cm2 to m2

        FORMAT='($,'' Atmospheric pressure (mmHg)                   '//
     &    '   [     760]: '',D11.4)'
        RDEFLT(1)=760.0d0 !mmHg
        IF(IOTYPE.EQ.3) RDATA(1)=P_ATM
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) P_ATM=RDATA(1)

        FORMAT='($,'' Cardiac output (mm3/s)                        '//
     &    '   [ 0.833E6]: '',D11.4)'
        RDEFLT(1)=0.833d6 !mmHg
        IF(IOTYPE.EQ.3) RDATA(1)=CARDIAC_OUTPUT
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) CARDIAC_OUTPUT=RDATA(1)

      ELSE
C Currently only set up for pulmonary problems. Change this assert 
C statement when other problem types are added. 
        CALL ASSERT((ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.11),
     '    ' >> Define values currently only set up for pulmonary '//
     '    'problems',ERROR,*9999)

      ENDIF !PROBLEM TYPE

      CALL EXITS('IPVALU')
      RETURN
 9999 CALL ERRORS('IPVALU',ERROR)
      CALL EXITS('IPVALU')
      RETURN 1
      END



