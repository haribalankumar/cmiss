C#### Module: CMISS_ARCHIVE4
C###  Description:
C###    Contains archived code from modules FE30 -> FEUSER

CFE30 Subroutine CALFUN     calls modal analysis routine
C     Subroutine CALC_WALL_SS
C     Subroutine FFRONT     activ. pattern with variable time step
C     Subroutine FRONT      activ. pattern with fixed t step-FHN
C     Subroutine REASSTB

CFE40 Subroutine ELAS3D    calc element matrix for 3D elasticity
C     Subroutine IPINI4     input i.c.s and boundary conditions
C     Subroutine OPINI4     output initial and boundary data

CFE50 Function AGP         integrand of Laplace eqtn for perfusion
C     Subroutine IPPRES    input pressure boundary conditions
C     Subroutine PFRF50    calc pressure loading residual terms
C###  Routine: UPVOROFROMDEFORMED update voro mesh after def of fe mesh

CFE70
C###  Routine: IOVSAE   writes output file for Vsaero input
C###  Routine: IPVSA    basic input for VSAERO
C###  Routine: IPVSB    patch geometry input for VSAERO
C###  Routine: IPVSC    wake input for VSAERO
C###  Routine: IPVSD    surface streamline input for VSAERO
C###  Routine: IPVSE    boundary layer input for VSAERO
C###  Routine: IPVSF    off-body velocity scan input for VSAERO
C###  Routine: IPVSG    off-body streamline input for VSAERO
C###  Routine: OPVSA    basic output for VSAERO
C###  Routine: OPVSB    patch geometry output for VSAERO
C###  Routine: OPVSC    wake output for VSAERO
C###  Routine: OPVSD    surface streamline output for VSAERO
C###  Routine: OPVSE    boundary layer output for VSAERO
C###  Routine: OPVSF    off-body velocity scan output for VSAERO
C###  Routine: OPVSG    off-body streamline output for VSAERO

CFE90 Subroutine HERM       generates extra nodes etc for hermite bem.
C     Subroutine UPBE       updates coords. if SOLVE4 used with optimis.
C###  Routine: GRADPHI_I   calcs grad of potential at internal points.
C###  Routine: GRADPHI_N   calcs gradient of the potential at a node.

CFEUSER Subroutine USER_7  helical fibre eqtn: fibre length prescribed
C       Subroutine LSFUN2  Calculates function and first deriv matrix
C                          for Nag routine E04GEF
 
Module FE30
=========== 

C 1/11/00 KAT removed from BLK30

C**** CDMATE(nm,ityp5,ityp2,ityp3) contains typical material values.

      CHARACTER*50
     '  CSTATIC_01_1(20),
     '  CSTATIC_02_1(20),
     '  CSTATIC_03_1(20),CSTATIC_03_2(20),CSTATIC_03_3(20),
     '  CSTATIC_04_1(20),CSTATIC_04_2(20),
     '  CSTATIC_05_1(20),CSTATIC_05_2(20),
     '  CSTATIC_06_1(20),
     '  CSTATIC_07_1(20),
     '  CSTATIC_08_1(20),CSTATIC_08_2(20),
     '  CSTATIC_09_1(20),CSTATIC_09_2(20),
     '  CSTATIC_10_1(20),
     '  CSTATIC_11_1(20),
     '  CSTATIC_12_1(20)

      EQUIVALENCE
     '  (CDMATE(1,1, 1,1),CSTATIC_01_1(1)),
     '  (CDMATE(1,1, 2,1),CSTATIC_02_1(1)),
     '  (CDMATE(1,1, 3,1),CSTATIC_03_1(1)),
     '  (CDMATE(1,1, 3,2),CSTATIC_03_2(1)),
     '  (CDMATE(1,1, 3,3),CSTATIC_03_3(1)),
     '  (CDMATE(1,1, 4,1),CSTATIC_04_1(1)),
     '  (CDMATE(1,1, 4,2),CSTATIC_04_2(1)),
     '  (CDMATE(1,1, 5,1),CSTATIC_05_1(1)),
     '  (CDMATE(1,1, 5,2),CSTATIC_05_2(1)),
     '  (CDMATE(1,1, 6,1),CSTATIC_06_1(1)),
     '  (CDMATE(1,1, 7,1),CSTATIC_07_1(1)),
     '  (CDMATE(1,1, 8,1),CSTATIC_08_1(1)),
     '  (CDMATE(1,1, 8,2),CSTATIC_08_2(1)),
     '  (CDMATE(1,1, 9,1),CSTATIC_09_1(1)),
     '  (CDMATE(1,1, 9,2),CSTATIC_09_2(1)),
     '  (CDMATE(1,1,10,1),CSTATIC_10_1(1)),
     '  (CDMATE(1,1,11,1),CSTATIC_11_1(1)),
     '  (CDMATE(1,1,12,1),CSTATIC_12_1(1))

      DATA CSTATIC_08_1/'(water=1.138e-4,air=0.145e-4)',
     '                  '(water=1000,air=1.226)       ',18*' '/

      CHARACTER*50
     '  CDYNAM_01_1(20),
     '  CDYNAM_02_1(20),
     '  CDYNAM_03_1(20),
     '  CDYNAM_04_1(20),
     '  CDYNAM_05_1(20),CDYNAM_05_2(20),CDYNAM_05_3(20),
     '  CDYNAM_06_1(20),
     '  CDYNAM_07_1(20),
     '  CDYNAM_08_1(20),
     '  CDYNAM_09_1(20),CDYNAM_09_2(20),CDYNAM_09_3(20),CDYNAM_09_4(20),
     '  CDYNAM_09_5(20),CDYNAM_09_6(20),CDYNAM_09_7(20),CDYNAM_09_8(20),
     '  CDYNAM_10_1(20),CDYNAM_10_2(20),
     '  CDYNAM_11_1(20),
     '  CDYNAM_12_1(20)

      EQUIVALENCE
     '  (CDMATE(1,2, 1,1),CDYNAM_01_1(1)),
     '  (CDMATE(1,2, 2,1),CDYNAM_02_1(1)),
     '  (CDMATE(1,2, 3,1),CDYNAM_03_1(1)),
     '  (CDMATE(1,2, 4,1),CDYNAM_04_1(1)),
     '  (CDMATE(1,2, 5,1),CDYNAM_05_1(1)),
     '  (CDMATE(1,2, 5,2),CDYNAM_05_2(1)),
     '  (CDMATE(1,2, 5,3),CDYNAM_05_3(1)),
     '  (CDMATE(1,2, 6,1),CDYNAM_06_1(1)),
     '  (CDMATE(1,2, 7,1),CDYNAM_07_1(1)),
     '  (CDMATE(1,2, 8,1),CDYNAM_08_1(1)),
     '  (CDMATE(1,2, 9,1),CDYNAM_09_1(1)),
     '  (CDMATE(1,2, 9,2),CDYNAM_09_2(1)),
     '  (CDMATE(1,2, 9,3),CDYNAM_09_3(1)),
     '  (CDMATE(1,2, 9,4),CDYNAM_09_4(1)),
     '  (CDMATE(1,2, 9,5),CDYNAM_09_5(1)),
     '  (CDMATE(1,2, 9,6),CDYNAM_09_6(1)),
     '  (CDMATE(1,2, 9,7),CDYNAM_09_7(1)),
     '  (CDMATE(1,2, 9,8),CDYNAM_09_8(1)),
     '  (CDMATE(1,2,10,1),CDYNAM_10_1(1)),
     '  (CDMATE(1,2,10,2),CDYNAM_10_2(1)),
     '  (CDMATE(1,2,11,1),CDYNAM_11_1(1)),
     '  (CDMATE(1,2,12,1),CDYNAM_12_1(1))

      CHARACTER*50
     '  CMODAL_01_1(20),
     '  CMODAL_02_1(20),
     '  CMODAL_03_1(20),
     '  CMODAL_04_1(20),CMODAL_04_2(20),
     '  CMODAL_05_1(20),
     '  CMODAL_06_1(20),
     '  CMODAL_07_1(20),
     '  CMODAL_08_1(20),
     '  CMODAL_09_1(20),
     '  CMODAL_10_1(20),
     '  CMODAL_11_1(20),
     '  CMODAL_12_1(20)

      EQUIVALENCE
     '  (CDMATE(1,3, 1,1),CMODAL_01_1(1)),
     '  (CDMATE(1,3, 2,1),CMODAL_02_1(1)),
     '  (CDMATE(1,3, 3,1),CMODAL_03_1(1)),
     '  (CDMATE(1,3, 4,1),CMODAL_04_1(1)),
     '  (CDMATE(1,3, 4,2),CMODAL_04_2(1)),
     '  (CDMATE(1,3, 5,1),CMODAL_05_1(1)),
     '  (CDMATE(1,3, 6,1),CMODAL_06_1(1)),
     '  (CDMATE(1,3, 7,1),CMODAL_07_1(1)),
     '  (CDMATE(1,3, 8,1),CMODAL_08_1(1)),
     '  (CDMATE(1,3, 9,1),CMODAL_09_1(1)),
     '  (CDMATE(1,3,10,1),CMODAL_10_1(1)),
     '  (CDMATE(1,3,11,1),CMODAL_11_1(1)),
     '  (CDMATE(1,3,12,1),CMODAL_12_1(1))

      CHARACTER*50
     '  CQUASIST_01_1(20),
     '  CQUASIST_02_1(20),
     '  CQUASIST_03_1(20),
     '  CQUASIST_04_1(20),
     '  CQUASIST_05_1(20),
     '  CQUASIST_06_1(20),
     '  CQUASIST_07_1(20),
     '  CQUASIST_08_1(20),
     '  CQUASIST_09_1(20),
     '  CQUASIST_10_1(20),
     '  CQUASIST_11_1(20),
     '  CQUASIST_12_1(20)

      EQUIVALENCE
     '  (CDMATE(1,4, 1,1),CQUASIST_01_1(1)),
     '  (CDMATE(1,4, 2,1),CQUASIST_02_1(1)),
     '  (CDMATE(1,4, 3,1),CQUASIST_03_1(1)),
     '  (CDMATE(1,4, 4,1),CQUASIST_04_1(1)),
     '  (CDMATE(1,4, 5,1),CQUASIST_05_1(1)),
     '  (CDMATE(1,4, 6,1),CQUASIST_06_1(1)),
     '  (CDMATE(1,4, 7,1),CQUASIST_07_1(1)),
     '  (CDMATE(1,4, 8,1),CQUASIST_08_1(1)),
     '  (CDMATE(1,4, 9,1),CQUASIST_09_1(1)),
     '  (CDMATE(1,4,10,1),CQUASIST_10_1(1)),
     '  (CDMATE(1,4,11,1),CQUASIST_11_1(1)),
     '  (CDMATE(1,4,12,1),CQUASIST_12_1(1))

      CHARACTER*50
     '  CFRONT_01_1(20),
     '  CFRONT_02_1(20),
     '  CFRONT_03_1(20),
     '  CFRONT_04_1(20),
     '  CFRONT_05_1(20),
     '  CFRONT_06_1(20),
     '  CFRONT_07_1(20),
     '  CFRONT_08_1(20),
     '  CFRONT_09_1(20),
     '  CFRONT_10_1(20),
     '  CFRONT_11_1(20),
     '  CFRONT_12_1(20)

      EQUIVALENCE
     '  (CDMATE(1,5, 1,1),CFRONT_01_1(1)),
     '  (CDMATE(1,5, 2,1),CFRONT_02_1(1)),
     '  (CDMATE(1,5, 3,1),CFRONT_03_1(1)),
     '  (CDMATE(1,5, 4,1),CFRONT_04_1(1)),
     '  (CDMATE(1,5, 5,1),CFRONT_05_1(1)),
     '  (CDMATE(1,5, 6,1),CFRONT_06_1(1)),
     '  (CDMATE(1,5, 7,1),CFRONT_07_1(1)),
     '  (CDMATE(1,5, 8,1),CFRONT_08_1(1)),
     '  (CDMATE(1,5, 9,1),CFRONT_09_1(1)),
     '  (CDMATE(1,5,10,1),CFRONT_10_1(1)),
     '  (CDMATE(1,5,11,1),CFRONT_11_1(1)),
     '  (CDMATE(1,5,12,1),CFRONT_12_1(1))

      CHARACTER*50
     '  CBUCKLE_01_1(20),
     '  CBUCKLE_02_1(20),
     '  CBUCKLE_03_1(20),
     '  CBUCKLE_04_1(20),
     '  CBUCKLE_05_1(20),
     '  CBUCKLE_06_1(20),
     '  CBUCKLE_07_1(20),
     '  CBUCKLE_08_1(20),
     '  CBUCKLE_09_1(20),
     '  CBUCKLE_10_1(20),
     '  CBUCKLE_11_1(20),
     '  CBUCKLE_12_1(20)

      EQUIVALENCE
     '  (CDMATE(1,6, 1,1),CBUCKLE_01_1(1)),
     '  (CDMATE(1,6, 2,1),CBUCKLE_02_1(1)),
     '  (CDMATE(1,6, 3,1),CBUCKLE_03_1(1)),
     '  (CDMATE(1,6, 4,1),CBUCKLE_04_1(1)),
     '  (CDMATE(1,6, 5,1),CBUCKLE_05_1(1)),
     '  (CDMATE(1,6, 6,1),CBUCKLE_06_1(1)),
     '  (CDMATE(1,6, 7,1),CBUCKLE_07_1(1)),
     '  (CDMATE(1,6, 8,1),CBUCKLE_08_1(1)),
     '  (CDMATE(1,6, 9,1),CBUCKLE_09_1(1)),
     '  (CDMATE(1,6,10,1),CBUCKLE_10_1(1)),
     '  (CDMATE(1,6,11,1),CBUCKLE_11_1(1)),
     '  (CDMATE(1,6,12,1),CBUCKLE_12_1(1))


C MHT 03-05-01 Archiving routine.
      SUBROUTINE CALC_WALL_SS(alpha,beta,Clumen,dH,Kd,L,rm,T_ss,Tlumen,
     '  Twall,YG,TUBE,HUMIDIFIER,ERROR,*)

C#### Subroutine: CALC_WALL_SS
C###  Description:
C###    CALC_WALLTEMP2 finds the temperature at the Mucus-Air-Interface
C###    (MAI) in an airway, given the temperature and water vapour
C###    concentration at the airway centre, and the wall temperature.
C***  Created by Merryn Howatson Tawhai, August 2000

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mesh00.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

!    Parameter List
      REAL*8 alpha,beta,Clumen,dH,Kd,L,rm,T_ss(13),Tlumen,Twall,
     '  YG(NIYGM)
      LOGICAL TUBE, HUMIDIFIER
      CHARACTER ERROR*(*)
!    Local Variables
      INTEGER nnit,NITER
      REAL*8 Cmai,depth,devap,dFdTmai,dheat,F,heat_term,Kl,Kpvc,Kt,Kw,
     '  Kw0,Lpvc,R,rho_w,temp_vol,Tmai,Tmai_old,Vevap,vol_limit,Vsource,
     '  ZEROTOL
      LOGICAL FOUND

      CALL ENTERS('CALC_WALL_SS',*9999)

C  Set parameters
      NITER=1000
      Tmai=YG(6)
      ZEROTOL=ZERO_TOL*1.d6
      Kw0=0.55d-6 !thermal cond. through ASL (kJ/mm/s/K)
      Kl=0.268d-7 !thermal conductivity in lumen (air/water vapour)
      depth=YG(4)     
      R=rm-depth
      Kt=0.55d-6 !thermal conductivity through tissue
      rho_w=0.11d-2 !density of water (ASL)
      IF(TUBE)THEN        
        Kpvc=0.92d-7 !thermal conductivity through pvc (from F&P)
        Lpvc=2.d0 !thickness of ET tube
        IF(HUMIDIFIER) Kt=1.d0
      ELSE
        Kpvc=1.d0
        Lpvc=0.d0
      ENDIF
      vol_limit=PI*(rm**2.d0-R**2.d0)*L
      temp_vol=vol_limit
      FOUND=.FALSE.
      nnit=0
      DO WHILE(.NOT.FOUND)
        nnit=nnit+1
        Cmai=RH_MAI*0.4035584d0*DEXP(-4.97d3/Tmai)
        Vevap=Kd*alpha*(Cmai-Clumen)*2.d0*PI*L*DT/rho_w
        devap=rm-(rm**2.d0-(vol_limit-Vevap)/(PI*L))**0.5d0
        IF(TUBE)THEN
          Vsource=0.d0
          YG(12)=0.d0
        ELSE
          IF(devap.LT.0.01d0)THEN
            Vsource=REPL_ASL*0.002d0/3.d0*PI*rm*L*DT
            devap=rm-(rm**2.d0-(vol_limit-Vevap+Vsource)/(PI*L))**0.5d0
            IF(devap.GT.0.01d0)THEN
              Vsource=(rm**2.d0-R**2.d0)*PI*L+Vevap-vol_limit
            ENDIF
          ELSE IF(devap.GT.0.01d0)THEN
            Vsource=-REPL_ASL*0.002d0/3.d0*PI*rm*L*DT
            devap=rm-(rm**2.d0-(vol_limit-Vevap+Vsource)/(PI*L))**0.5d0
            IF(devap.LT.0.01d0)THEN
              Vsource=(rm**2.d0-R**2.d0)*PI*L+Vevap-vol_limit
            ENDIF
          ELSE
            Vsource=0.d0
          ENDIF
        ENDIF !TUBE
        depth=rm-(rm**2.d0-temp_vol/(PI*L))**0.5d0
        IF(depth.LE.ZERO_TOL*1.d4)THEN
          depth=0.d0
          Kw=1.d0
          heat_term=0.d0
          dheat=0.d0
        ELSE
          Kw=Kw0
          heat_term=dH*Kd*(Clumen-Cmai)*alpha/R
          dheat=dH*Kd*alpha*Cmai*4.97d3/Tmai**2.d0/R
        ENDIF
C        IF(EXPN.AND.heat_term.LE.0.d0) heat_term=0.d0
        Tmai_old=Tmai
        F=Tmai*(Kl*beta*(Kpvc*Kt*depth+Kw*Kpvc*dely+Kw*Kt*Lpvc)
     '    -Kw*Kpvc*Kt*R)/(R*(Kpvc*Kt*depth+Kw*Kpvc*dely+Kw*Kt*Lpvc))
     '    +Twall*Kw*Kt*Kpvc/(Kpvc*Kt*depth+Kw*Kpvc*dely+Kw*Kt*Lpvc)
     '    -Tlumen*Kl*beta/R+heat_term
        dFdTmai=(Kl*beta*(Kpvc*Kt*depth+Kw*Kpvc*dely+Kw*Kt*Lpvc)
     '    -Kw*Kpvc*Kt*R)/(R*(Kpvc*Kt*depth+Kw*Kpvc*dely+Kw*Kt*Lpvc))
     '    -dheat
        Tmai=Tmai-F/dFdTmai
        IF(DABS(Tmai-Tmai_old).LE.ZEROTOL) FOUND=.TRUE.
        IF(nnit.GT.NITER)THEN
          FOUND=.TRUE.
          Tmai=YG(6)
          depth=YG(3)
        ENDIF
      ENDDO
      T_ss(1)=Tmai
      IF(depth.GT.0.d0)THEN
        T_ss(5)=(Tmai*Kw*Kpvc*dely+Tmai*Kw*Kt*Lpvc+Twall*Kpvc*Kt*depth)
     '    /(Kpvc*Kt*depth+Kw*Kpvc*dely+Kw*Kt*Lpvc)
        IF(TUBE)THEN
          T_ss(9)=(Twall*Kt*Lpvc+T_ss(5)*Kpvc*dely)/(Kpvc*dely+Kt*Lpvc)
        ELSE
          T_ss(9)=T_ss(5)
        ENDIF
      ELSE
        IF(TUBE)THEN
          T_ss(5)=(Twall*Kt*Lpvc+Tmai*Kpvc*dely)/(Kpvc*dely+Kt*Lpvc)
          T_ss(9)=T_ss(5)
        ELSE
          T_ss(5)=Tmai
          T_ss(9)=Tmai
        ENDIF
      ENDIF
      T_ss(13)=Twall
      T_ss(2)=(T_ss(1)*0.75d0+T_ss(5)*0.25d0)
      T_ss(3)=(T_ss(1)*0.50d0+T_ss(5)*0.50d0)
      T_ss(4)=(T_ss(1)*0.25d0+T_ss(5)*0.75d0)
      T_ss(6)=(T_ss(5)*0.75d0+T_ss(9)*0.25d0)
      T_ss(7)=(T_ss(5)*0.50d0+T_ss(9)*0.50d0)
      T_ss(8)=(T_ss(5)*0.25d0+T_ss(9)*0.75d0)
      T_ss(10)=(T_ss(9)*0.75d0+Twall*0.25d0)
      T_ss(11)=(T_ss(9)*0.50d0+Twall*0.50d0)
      T_ss(12)=(T_ss(9)*0.25d0+Twall*0.75d0)
      
      CALL EXITS('CALC_WALL_SS')
      RETURN
 9999 CALL ERRORS('CALC_WALL_SS',ERROR)
      CALL EXITS('CALC_WALL_SS')
      RETURN 1
      END



C 25/2/97 LC removed section from : news AJP 25/1/96

C#### Subroutine: IPBAS3
C###  Description:
C###    IPBAS3 inputs basis functions for dependent variables for 
C###    FE30 problems.

C news AJP 25/1/96
c      IF(ITYP4(nr,nx).EQ.4) THEN                      !Collocation
cC GBS 5-DEC-1994   Is this really necessary?
cC        CALL ASSERT(NQT.GT.0,'Collocation pts not defined',ERROR,*9999)
c        DO noelem=1,NEELEM(0,nr)
c          ne=NEELEM(noelem,nr)
c          DO nc=1,2
c            NBH(NH_LOC(1,nx),nc,ne)=NBJ(1,ne) !to allow listing of i.c.s
c          ENDDO
c        ENDDO
c
c      ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.8) THEN !Threshold 
c        DO noelem=1,NEELEM(0,nr)                            !model
c          ne=NEELEM(noelem,nr)
c          DO nc=1,2
c            NBH(NH_LOC(1,nx),nc,ne)=NBJ(1,ne) !to allow listing of i.c.s
c          ENDDO
c        ENDDO
c
c      ELSE IF(ITYP2(nr,nx).EQ.9) THEN                 !Cardiac activ.n
cC GBS 5-DEC-1994   Is this really necessary?
cC        CALL ASSERT(NQT.GT.0,'Global Pts not defined',ERROR,*9999)
c        DO noelem=1,NEELEM(0,nr)
c          ne=NEELEM(noelem,nr)
c          DO nc=1,2
c            NBH(NH_LOC(1,nx),nc,ne)=NBJ(1,ne) !to allow listing of i.c.s
c          ENDDO
c        ENDDO
c
c      ELSE                                            !All other models
C newe ajp 25/1/96

C 25/2/97 LC removed section from : 
C     cpb 11/9/95 Adding zero cross derivative basis functions.

C#### Subroutine: IPBAS3
C###  Description:
C###    IPBAS3 inputs basis functions for dependent variables for 
C###    FE30 problems.

c cpb 11/9/95 Adding zero cross derivative basis functions.
C        IF(BDRY_ELEMENT) THEN
C          HERMITE_3D=.FALSE.
C          noelem=1
C          DO WHILE(.NOT.HERMITE_3D.AND.noelem.LE.NEELEM(0,nr))
C            ne=NEELEM(noelem,nr)
C            nb=NBASEF(NBH(NH_LOC(1,nx),1,ne),1)
C            IF(NJE(ne).EQ.3) THEN
C              IF((IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2).OR.
C     '          (IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4).OR.
C     '          ((IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6).AND.
C     '          (IBT(2,1,nb).EQ.4.AND.IBT(1,2,nb).EQ.2)).OR.
C     '          ((IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6).AND.
C     '          (IBT(2,2,nb).EQ.4.AND.IBT(1,1,nb).EQ.2))) THEN
CC               Hermite in each direction or hermite simplex or 
CC               Hermite sector elements
C                HERMITE_3D=.TRUE.
C              ENDIF
C            ENDIF
C            noelem=noelem+1
C          ENDDO !noelem
C          IF(HERMITE_3D) THEN
C            DO nc=1,2
C              KTYP93(nc,nr)=0
C            ENDDO !nc
C            FORMAT='(/$,'' Do you want to set cross derivatives of'//
C     '        ' the dependent variables to zero [Y]? '',A)'
C            IF(IOTYPE.EQ.3) ADATA(1)='Y'
C            CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '        1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(ADATA(1).EQ.'Y') KTYP93(1,nr)=1
C            IF(KTYP93(1,nr).EQ.1) THEN !Check normal deriv interpolation
C              FOUND=.FALSE.
C              noelem=1
C              DO WHILE(.NOT.FOUND.AND.noelem.LE.NEELEM(0,nr))
C                ne=NEELEM(noelem,nr)
C                nb=NBASEF(NBH(NH_LOC(1,nx),2,ne),1)
C                IF((IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2).OR.
C     '            (IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4).OR.
C     '            ((IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6).AND.
C     '            (IBT(2,1,nb).EQ.4.AND.IBT(1,2,nb).EQ.2)).OR.
C     '            ((IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6).AND.
C     '            (IBT(2,2,nb).EQ.4.AND.IBT(1,1,nb).EQ.2))) THEN
CC                 Hermite in each direction or hermite simplex or 
CC                 Hermite sector elements
C                  FOUND=.TRUE.
C                  KTYP93(2,nr)=1
C                ENDIF
C                noelem=noelem+1
C              ENDDO
C            ENDIF
C          ELSE
C            KTYP93(1,nr)=0
C            KTYP93(2,nr)=0
C          ENDIF
C        ENDIF
C      ENDIF !problem type


      SUBROUTINE CALFUN(IBT,IDO,INP,LEG,LGE,ME,NBH,NBJ,NHE,
     '  NHP,NJE,NKE,NKH,NPB,NPE,NPNODE,NQE,NW,NYCZ,NYRZ,NZD,
     '  A,B,CE,CG,CP,ED,EM,ER,ES,GM,GS,
     '  GSCOPY,HK,PG,RE,RG,SE,VE,WG,
     '  XA,XE,XF,XG,XP,YP,ZA,ZE,ZF,ZG,ZP,
     '  EIGPRM,PRGS,PRGM,EIGVM,EIGV,EIGOLD,
     '  FIX,LNY,M,N,F,X,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),
     '  LEG(*),LGE(*),M,ME(*),N,NBH(NHM,NCM,*),NBJ(NJM,*),NHE(*),
     '  NHP(*),NJE(*),NKE(NKM,NNM,NBM,*),NKH(NHM,NCM,*),NPB(0:NPM,0:*),
     '  NPE(NNM,NBM,*),NPNODE(0:NPM,0:*),
     '  NQE(NSM,NBM,*),NW(NEM,*),NYCZ(*),NYRZ(*),NZD(*)
      REAL*8 A(NSM,*),B(NSM,*),CE(NMM,*),CG(NMM,*),CP(NMM,*),
     '  ED(NVM,*),EIGOLD(*),EIGPRM(NTM,*),EIGV(*),EIGVM(NYM,*),
     '  EM(NVM,*),ER(*),ES(NVM,*),F(*),GM(*),GS(*),GSCOPY(*),HK(*),
     '  PG(NSM,NUM,NGM,*),PRGS(NTM,*),PRGM(NTM,*),RE(NSM,*),RG(*),
     '  SE(NSM,NBM,*),VE(NSM,NKM,*),WG(NGM,*),X(*),XA(NAM,NJM,*),
     '  XE(NSM,*),XF(NSM,*),XG(NJM,*),XP(NKM,NJM,*),YP(NYM,*),
     '  ZA(NAM,NHM,NCM,*),ZE(NSM,*),ZF(NSM,*),ZG(NHM,*),
     '  ZP(NKM,NHM,NCM,*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*),LNY(*)
!     Local Variables
      INTEGER ND,NONODE,NP,NQ,NQR,NQTR,NR
      REAL*8 RL

      CALL ENTERS('CALFUN',*9999)
      NQR=0
      DO NQ=1,NQT
        IF(.NOT.FIX(NQ,4)) THEN
          NQR=NQR+1
          IF(ILP(1,1).EQ.2) CE(1,NQ)=X(NQR)
          IF(ILP(1,1).EQ.3) CP(1,NQ)=X(NQR)
        ENDIF
      ENDDO                  
      NQTR=NQR
      IF(.NOT.FIX(NQT+1,4)) THEN
        NQTR=NQTR+1
        RL=X(NQTR)
        DO NONODE=1,NPNODE(0,1)
          NP=NPNODE(NONODE,1)
          XP(1,1,NP)=RL*DBLE(NP-1)/DBLE(NPT(1)-1)
        ENDDO
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' Current X-vector:''/,(10E11.4))') 
     '      (X(NQ),NQ=1,NQTR)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Current nodal coords:''/,(10E11.4))')
     '      (XP(1,1,NP),NP=1,NPT(1))    
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      NR=1 !should be generalized
      CALL MODAL(IBT,IDO,INP,LEG,LGE,ME,NBH,NBJ,
     '  NHE,NHP,NJE,NKE,NKH,NPB,NPE,NQE,NR,NW,NYCZ,NYRZ,NZD,
     '  A,B,CE,CG,CP,ED,EM,ER,ES,GM,GS,GSCOPY,HK,
     '  PG,RE,RG,SE,VE,WG,XA,XE,XF,XG,XP,YP,ZA,ZE,ZF,ZG,ZP,
     '  EIGPRM,PRGS,PRGM,EIGVM,EIGV,EIGOLD,FIX,LNY)

      DO ND=1,NDT
        F(ND)=YP(ND,4)-EIGV(ND)
      ENDDO
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' Current eigenvalues:''/,(10E11.4))')
     '      (EIGV(ND),ND=1,NDT)   
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Current F-vector:''/,(10E11.4))')
     '      (F(ND),ND=1,NDT)   
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('CALFUN')
      RETURN
 9999 CALL ERRORS('CALFUN',ERROR)
      CALL EXITS('CALFUN')
      RETURN 1
      END


      SUBROUTINE FFRONT(IBT,IME,IMG,INE,ING,INP,ITHRES,LE,LG,MSA,
     '  NBJ,NCNE,NCNP,NEELEM,NEHIS,NEIGM,NEMAX,NGHIS,NGMAX,
     '  NHIS,NJE,NKE,NLET,NLGT,
     '  NMAX,NPE,NPF,NPNODE,NQGE,NRE,NSE1,NSE2,NSE3,NUAT,NUATT,NW,NXI,
     '  CE,CG,CP,PG,SE,SSTA,STA,T,THIS,THRES,TA,VE,XA,XE,XG,XIG,
     '  XP,ZE,ZG,ERROR,*)

C**** Calculates activation pattern by threshold modelling.
C**** NE,NG are element and gauss pt number of current primary (active) point.
C**** ME,MG  "   "  "    "    "   "    "    "  current secondary point.
C**** NGP3(i,j,k),i,j,k=1,3 are the Gauss point numbers ng for 3*3*3
C**** NGP5(i,j,k),i,j,k=1,5 are the Gauss point numbers ng for 5*5*5
C**** MNE(i,j,k),i,j,k=-1,1 are the element numbers surrounding ne.
C**** ITHRES(1,ng,ne) is 0 if Gauss point ng of element ne is not active
C****   "                1  "    "    "    "  "    "     "  " active
C****   "                2  "    "    "   is surrounded by active points
C**** ITHRES(2,ng,ne) is 0  "    "    "   ng is ordinary myocardium
C****   "                1  "    "    "    "  " Purkinje tissue
C**** THRES(1,ng,ne) is time Gauss point became active
C**** NUAT is the number of Gauss points in the active set (ITHRES=1)
C**** NLET is the number of elements with active points in them
C**** LE is the list of elements with active points in them
C**** LG is the list of active Gauss points within an active element
C**** NGHIS (& NEHIS, THIS) are the Gauss points to be activated from history
C****	(ie from data known about the times of activation of epicardial
C****   points)
C**** TA(NG,NE,nnei) is the time that point (ng,ne) would activate
C****   neighbour nnei
C**** STA(NG,NE) is the smallest TA for point (ng,ne)
C**** SSTA(nsa) are the NSA smallest STA's in this iteration
C**** KTYP31 = 1 for forward tracking of model
C****   "    = 2 for backward   "     "    "
C**** The current active Gauss point (and it's element) are designated "primary"
C**** Points or  elements adjacent to primary point or  element are
C****   "secondary" or "neighbours".
C**** MSA is the number of  points allowed to be activated this iteration
C**** NSA is no of points to activate found so far this iteration (NSA.LE.MSA)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:suben00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER NEMAX,NGMAX,NEIGM
      INTEGER IBT(2,NIM,*),IME(2000),IMG(2000),INE(2000),ING(2000),
     '  INP(NNM,NIM,*),ITHRES(3,NGM,*),LE(*),LG(NEMAX,*),MSA,
     '  NBJ(NJM,*),NEELEM(0:NEM,0:*),
     '  NEHIS(*),NGHIS(*),NHIS,NJE(*),
     '  NKE(NKM,NNM,NBM,*),NLET,NLGT(*),NMAX,NPE(NNM,NBM,*),NPF(12,*),
     '  NPNODE(0:NPM,0:*),
     '  NQGE(NGM,*),NRE(*),NSE1,NSE2,NSE3,NUAT,NUATT,NW(NEM,*),
     '  NXI(-NIM:NIM,0:*)
      REAL*8 CE(NMM,*),CG(NMM,*),CP(NMM,*),PG(NSM,NUM,NGM,*),
     '  SE(NSM,NBM,*),SSTA(*),STA(NGMAX,*),
     '  T,TA(NGMAX,NEMAX,NEIGM),THIS(*),THRES(3,NGM,*),VE(NSM,NKM,*),
     '  XA(NAM,NJM,*),XE(NSM,*),XG(NJM,*),XIG(NIM,NGM,*),
     '  XP(NKM,NJM,*),ZE(NSM,*),ZG(NHM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NSAMX
      PARAMETER (NSAMX=50)
      INTEGER I,I1,I2,I3,II1,ISUM,J,J1,J2,J3,JJ1,K,K1,K2,K3,KK1,
     '  ME,MEA(NSAMX),MG,MGA(NSAMX),MNE(-1:1,-1:1,-1:1),
     '  NEA(NSAMX),NGA(NSAMX),NGP3(3,3,3),NGP5(5,5,5),NA,NB,NE,NG,NH,
     '  NITB,NJ1,NK1,NLE1,NLG,NNEIG,NREAC,NSA
      REAL*8 CG1,CG2,CG3,COSALF,COSETA,DUMSQ,DX,DX1,DX2,DX3,DXI1,DXI2,
     '  DXI3,
     '  DXIX(3,3),ETA,GL(3,3),GU(3,3),R,RG,SINALF,SINETA,SQ,STIME,TIME
      LOGICAL EFLAG,EQUALS,AFLAG

      DATA NGP3/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     '  19,20,21,22,23,24,25,26,27/
      DATA NGP5/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     '  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     '  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,
     '  55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,
     '  73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,
     '  91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,
     '  107,108,109,110,111,112,113,114,115,116,117,118,119,120,
     '  121,122,123,124,125/

      CALL ENTERS('FFRONT',*9999)
      CALL NENXI(IBT,INP,NBJ,NEELEM,NPE,NXI,ERROR,*9999)  
                                               !Returns list NXI of elements
      NSA=0                                    !  surrounding NE.
      DO NA=1,NSAMX
        SSTA(NA)=9900.0d0
      ENDDO
C *** For each element in list of active elements,calculate next activated
      DO 500 NLE1=1,NLET
	NE=LE(NLE1)
        CALL XPXE(NBJ(1,NE),NCNE(1,1,ne,NRE(ne)),NCNP,
     '    NJE(NE),NKE(1,1,1,NE),NPE(1,1,NE),NPF(1,1),
     '    NQGE(1,NE),NRE(ne),SE(1,1,NE),XA,XE,XP,ERROR,*9999)
        CALL CPCG(1,NBJ(1,NE),NPE(1,1,NE),CE(1,NE),CG,CP,PG,ERROR,*9999)
        NB=NBJ(1,NE)
        NITB=NIT(NB)
C ***   Element numbers of current primary element and it's secondary elements.
        MNE(-1,-1,-1)=NXI(-1,NXI(-2,NXI(-3,NE)))
        MNE( 0,-1,-1)=NXI( 0,NXI(-2,NXI(-3,NE)))
        MNE( 1,-1,-1)=NXI( 1,NXI(-2,NXI(-3,NE)))
        MNE(-1, 0,-1)=NXI(-1,NXI( 0,NXI(-3,NE)))
        MNE( 0, 0,-1)=NXI( 0,NXI( 0,NXI(-3,NE)))
        MNE( 1, 0,-1)=NXI( 1,NXI( 0,NXI(-3,NE)))
        MNE(-1, 1,-1)=NXI(-1,NXI( 2,NXI(-3,NE)))
        MNE( 0, 1,-1)=NXI( 0,NXI( 2,NXI(-3,NE)))
        MNE( 1, 1,-1)=NXI( 1,NXI( 2,NXI(-3,NE)))
        MNE(-1,-1, 0)=NXI(-1,NXI(-2,NXI( 0,NE)))
        MNE( 0,-1, 0)=NXI( 0,NXI(-2,NXI( 0,NE)))
        MNE( 1,-1, 0)=NXI( 1,NXI(-2,NXI( 0,NE)))
        MNE(-1, 0, 0)=NXI(-1,NXI( 0,NXI( 0,NE)))
        MNE( 0, 0, 0)=NXI( 0,NXI( 0,NXI( 0,NE)))
        MNE( 1, 0, 0)=NXI( 1,NXI( 0,NXI( 0,NE)))
        MNE(-1, 1, 0)=NXI(-1,NXI( 2,NXI( 0,NE)))
        MNE( 0, 1, 0)=NXI( 0,NXI( 2,NXI( 0,NE)))
        MNE( 1, 1, 0)=NXI( 1,NXI( 2,NXI( 0,NE)))
        MNE(-1,-1, 1)=NXI(-1,NXI(-2,NXI( 3,NE)))
        MNE( 0,-1, 1)=NXI( 0,NXI(-2,NXI( 3,NE)))
        MNE( 1,-1, 1)=NXI( 1,NXI(-2,NXI( 3,NE)))
        MNE(-1, 0, 1)=NXI(-1,NXI( 0,NXI( 3,NE)))
        MNE( 0, 0, 1)=NXI( 0,NXI( 0,NXI( 3,NE)))
        MNE( 1, 0, 1)=NXI( 1,NXI( 0,NXI( 3,NE)))
        MNE(-1, 1, 1)=NXI(-1,NXI( 2,NXI( 3,NE)))
        MNE( 0, 1, 1)=NXI( 0,NXI( 2,NXI( 3,NE)))
        MNE( 1, 1, 1)=NXI( 1,NXI( 2,NXI( 3,NE)))
        IF(DOP) THEN
          WRITE(OP_STRING,'(//,'' *ELEMENT* '',I3,'' MNE:'',27I4)')
     '      NE,(((MNE(I,J,K),I=-1,1),J=-1,1),K=-1,1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C ***   For each active gauss point in active element
        DO 400 NLG=1,NLGT(NLE1)
          NG=LG(NLE1,NLG)
C ***     Calculate I, J, K corresponding to ng
	  IF(NMAX.EQ.3)THEN   !3*3*3 Gauss point
	    IF(NG.LT.10)THEN
              K=1
	    ELSE IF(NG.LT.19)THEN
              K=2   
            ELSE
              K=3
	    ENDIF
	    K1=(K-1)*9
	    IF(NG.LT.K1+4)THEN
              J=1
	    ELSE IF(NG.LT.K1+7)THEN
              J=2
            ELSE
              J=3
            ENDIF
	    I=MOD(NG,3)
	    IF(I.EQ.0)I=3
	  ELSE                  !5*5*5 Gauss point
	    IF(NG.LT.26)THEN
	      K=1
	    ELSE IF(NG.LT.51)THEN
              K=2
	    ELSE IF(NG.LT.76)THEN
              K=3
	    ELSE IF(NG.LT.101)THEN
              K=4
	    ELSE
              K=5
	    ENDIF
	    K1=(K-1)*25
	    IF(NG.LT.K1+6)THEN
              J=1
	    ELSE IF(NG.LT.K1+11)THEN
              J=2
	    ELSE IF(NG.LT.K1+16)THEN
              J=3
	    ELSE IF(NG.LT.K1+21)THEN
              J=4
	    ELSE
              J=5
	    ENDIF
	    I=MOD(NG,5)
	    IF(I.EQ.0)I=5
	  ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '(/,'' *Active point: '',''ne='',I3,'' ng='',
     '        I4,3x,'' I='',I2,'' J='',I2,'' K='',I2)') NE,NG,I,J,K
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL XEXG(NBJ(1,NE),NG,NJE(NE),PG,VE(1,1,NE),XE,XG,ERROR,
     '      *9999)
          CALL XGMG(0,NIT(1),1,NJE(NE),DXIX,GL,GU,RG,XG,ERROR,*9999)
          TIME=DABS(TINCR)
	  IF(JTYP9.EQ.1)THEN
            ETA=XG(NJ_LOC(NJL_FIBR,1),1)
	  ELSE
	    ETA=(ATAN2(XG(NJ_LOC(NJL_FIBR,2),1),XG(NJ_LOC(NJL_FIBR,1),
     '        1)))/2.0d0
	  ENDIF
          COSETA=DCOS(ETA)
          SINETA=DSIN(ETA)
          CG1=CG(1,NG)*TIME  !Dist. wavefront moves in dir-1 in this timestep
          CG2=CG(2,NG)*TIME       
C         CG3=CG(3,NG)*TIME
C ***     For endocardial Gauss points, speedup velocity as if Purkinje fibres
          IF(NE.GT.12.AND.NG.LT.NSUBEN)THEN
	    CG1=CG1*VSUBEN
C	    WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' NSUBEN='',I4)')NE,NG,CG1,NSUBEN          
	  ENDIF
          IF(IRVSUBEN.EQ.1)THEN
            IF(NE.EQ.13.AND.NG.GE.NGTOP_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN        
            ELSE IF(NE.EQ.16.AND.NG.GE.NGTOP_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN
            ELSE IF(NE.EQ.17.AND.NG.GE.NGTOP_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN
            ELSE IF(NE.EQ.20.AND.NG.GE.NGTOP_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN
            ELSE IF(NE.EQ.1.AND.NG.LE.NGBOT_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN
            ELSE IF(NE.EQ.4.AND.NG.LE.NGBOT_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN
            ELSE IF(NE.EQ.5.AND.NG.LE.NGBOT_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN
            ELSE IF(NE.EQ.8.AND.NG.LE.NGBOT_LAYER)THEN
              CG1=CG1*VSUBEN
C	      WRITE(59,'('' NE='',I3,'' NG='',I4,'' CG1='',E10.4,
C    '		'' IRVSUBEN='',I4)')NE,NG,CG1,IRVSUBEN
            ENDIF
          ENDIF
          IF(DOP) THEN 
             WRITE(OP_STRING,'('' *TINCR='',E10.3,
     '      '' ETA='',E10.3,/)') TIME,ETA
             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
	  ISUM=0
C ***     Have primary gauss point and local properties.
C ***     For each neighbour of primary gauss point:
	  IF(NITB.LE.2)THEN
	    NK1=0
	  ELSE
	    NK1=1
	  ENDIF                                               
          DO 300 K1=-NK1,NK1
	    KK1=K+K1
            K2=MOD(KK1,NMAX)
	    IF(K2.LE.0)K2=K2+NMAX
	    IF(KK1.LT.1)THEN
	      K3=-1
	    ELSE IF(KK1.GT.NMAX) THEN
	      K3=1
	    ELSE
	      K3=0
	    ENDIF
	    NJ1=NSE1
	    IF(K1.NE.0.AND.NSE3.EQ.1)NJ1=1
          DO 300 J1=-NJ1,NJ1
	    JJ1=J+J1
            J2=MOD(JJ1,NMAX)
	    IF(J2.LE.0) J2 = J2+NMAX
	    IF(JJ1.LT.1)THEN
	      J3=-1
	    ELSE IF(JJ1.GT.NMAX)THEN
	      J3=1
	    ELSE
	      J3=0
	    ENDIF
          DO 300 I1=-NJ1,NJ1
	    II1=I+I1
            I2=MOD(II1,NMAX)
	    IF(I2.LE.0)I2 = I2+NMAX
	    IF(II1.LT.1)THEN
	      I3=-1
	    ELSE IF(II1.GT.NMAX)THEN
	      I3=1
	    ELSE
	      I3=0
	    ENDIF
	    IF(NMAX.EQ.3)THEN
              MG=NGP3(I2,J2,K2)
	    ELSE
              MG=NGP5(I2,J2,K2)
	    ENDIF
            ME=MNE(I3,J3,K3)
C ***       Have primary and secondary point.
C           IF(DOP) THEN
C               WRITE(OP_STRING,
C		  '('' I1='',I2,'' J1='',I2,'' K1='',I2,
C     '       '' I2='',I2,'' J2='',I2,'' K2='',I2,
C     '       '' me='',I3,'' mg='',I3)') I1,J1,K1,I2,J2,K2,ME,MG
C               CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C           ENDIF
            IF(ME.GT.0) THEN
C ***         If gauss point mg is inactive, see whether it should
C ***         be activated in this iteration
              IF(ITHRES(1,MG,ME).EQ.0.
     '          AND.(.NOT.(MG.EQ.NG.AND.ME.EQ.NE))) THEN
		IF(NITB.EQ.3)THEN
		  NNEIG=(I1+NSE1+1)+(J1+NSE1)*(2*NSE1+1)+(K1+1)*
     '		    (2*NSE1+1)*(2*NSE1+1)
	  	ELSE
		  NNEIG=(I1+NSE1+1)+(J1+NSE1)*(2*NSE1+1)
		ENDIF
C ***           If its activation time from this primary point not yet calculated
	        IF(TA(NG,NE,NNEIG).GT.9000.0d0)THEN    ! calculate it
	  	  IF(I3.EQ.0) THEN
                    DXI1=(XIG(1,MG,NB)-XIG(1,NG,NB))
	  	  ELSE
	  	    IF(XIG(1,MG,NB).GT.XIG(1,NG,NB)) THEN
	  	      DXI1=-1.0d0+XIG(1,MG,NB)-XIG(1,NG,NB)
	     	    ELSE
	  	      DXI1=1.0d0-XIG(1,NG,NB)+XIG(1,MG,NB)
	  	    ENDIF
	  	  ENDIF
                  DX1=DXI1*DSQRT(GL(1,1))
 	  	  IF(J3.EQ.0) THEN
                    DXI2=(XIG(2,MG,NB)-XIG(2,NG,NB))
	  	  ELSE
	  	    IF(XIG(2,MG,NB).GT.XIG(2,NG,NB)) THEN
	  	      DXI2=-1.0d0+XIG(2,MG,NB)-XIG(2,NG,NB)
	  	    ELSE
	  	      DXI2=1.0d0-XIG(2,NG,NB)+XIG(2,MG,NB)
	  	    ENDIF
	  	  ENDIF
                  DX2=DXI2*DSQRT(GL(2,2))
                  IF(NITB.EQ.3) THEN
	  	    IF(K3.EQ.0) THEN
                      DXI3=(XIG(3,MG,NB)-XIG(3,NG,NB))
	  	    ELSE
	  	      IF(XIG(3,MG,NB).GT.XIG(3,NG,NB)) THEN
	  	        DXI3=-1.0d0+XIG(3,MG,NB)-XIG(3,NG,NB)
	  	      ELSE
	  	        DXI3=1.0d0-XIG(3,NG,NB)+XIG(3,MG,NB)
	  	      ENDIF
	  	    ENDIF
                    DX3=DXI3*DSQRT(GL(3,3))
                  ENDIF
c                 IF(DOP) 
C                    WRITE(OP_STRING,'('' DXI1='',E10.3,'' DXI2='',
c     '              E10.3,'' DXI3='',E10.3,'' DX1='',E10.3,'' DX2='',
c     '              E10.3,'' DX3='',E10.3)') DXI1,DXI2,DXI3,DX1,DX2,DX3
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                 ENDIF
                  DX=DSQRT(DX1*DX1+DX2*DX2+DX3*DX3)
                  IF(DABS(DX).GT.1.d-6) THEN
                    COSALF=(DX1*COSETA+DX2*SINETA)/DX
                    DUMSQ=COSALF*COSALF
                    IF(DUMSQ.GT.0.999999d0) THEN
                      SINALF=0.0d0
                    ELSE
                      SINALF=DSQRT(1.0d0-DUMSQ)
                    ENDIF
                    SQ=DSQRT((COSALF/CG1)**2+(SINALF/CG2)**2)
                    IF(SQ.GT.1.d-6) THEN
                      R=1.0d0/SQ         
                    ELSE
                      WRITE(OP_STRING,*) ' SQ is zero'
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                      WRITE(IO4,*) ' SQ is zero'
                    ENDIF
	  	    TA(NG,NE,NNEIG)=THRES(1,NG,NE)+DABS(TINCR)*DX/R
	  	    IF(DOP)THEN
                      WRITE(OP_STRING,
     '                   '('' TA-New: NNEIG='',I4,'' NG='',I3,
     '                   '' NE='',I3,'' R='',E10.4,'' DX='',E10.4,
     '                   '' THRES(1,NG,NE)='',E10.4,'' TA='',E10.4)')
     '                  NNEIG,NG,NE,R,DX,THRES(1,NG,NE),TA(NG,NE,NNEIG)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF                       
                  ELSE
                    WRITE(OP_STRING,
     '                '('' DX is zero. NE='',I4,'' NG='',I4)')
     '                NE,NG          
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
	  	  IF(DOP)THEN
                    WRITE(OP_STRING,
     '                 '('' TA-Old: NNEIG='',I4,'' NG='',I3,
     '                 '' NE='',I3,27X,'' THRES(1,NG,NE)='',
     '                 E10.4,'' TA='',E10.4)')
     '              NNEIG,NG,NE,THRES(1,NG,NE),TA(NG,NE,NNEIG)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999) 
                  ENDIF
	        ENDIF
C ***           If this is earliest TA to date: (STA is remembered from previous iterations)
	  	IF(TA(NG,NE,NNEIG).LE.STA(MG,ME))THEN
		  IF(KTYP31.EQ.1)THEN        !Foward model
		    IF(TA(NG,NE,NNEIG).GT.1.0d-4)THEN
                      IF(DOP)THEN
                        IF(STA(MG,ME).LT.9900.0d0) THEN !(Initialized at 10000.0)
                          IF(TA(NG,NE,NNEIG).LT.STA(MG,ME))THEN
                            WRITE(OP_STRING,'(''    New     STA   '',
     '                        ''  -  update    STA  '','' MG='',
     '                        I4,'' ME='',I4,
     '                        '' Previous STA ='',E10.4,'' TA='',
     '                        E10.4)')               
     '                        MG,ME,STA(MG,ME),TA(NG,NE,NNEIG)
                             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ELSE
                            WRITE(OP_STRING,'(''    Old     STA   '',
     '                        ''  -  unchanged STA  '','' MG='',I4,
     '                        '' ME='',I4,
     '                        '' Previous STA ='',E10.4,'' TA='',
     '                        E10.4)')
     '                        MG,ME,STA(MG,ME),TA(NG,NE,NNEIG)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ELSE
                          WRITE(OP_STRING,'(''    First   STA   '',
     '                      ''  -  update    STA  '','' MG='',I4,
     '                      '' ME='',I4,
     '                      '' Previous STA ='',E10.4,'' TA='',E10.4)')
     '                      MG,ME,STA(MG,ME),TA(NG,NE,NNEIG)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF
	  	      STA(MG,ME)=TA(NG,NE,NNEIG)
		    ENDIF
		  ELSE                       !Backward
                    IF(DOP)THEN
                      IF(STA(MG,ME).LT.9900.0d0) THEN !(Initialized at 10000.0)
                        IF(TA(NG,NE,NNEIG).LT.STA(MG,ME))THEN
                          WRITE(OP_STRING,'(''    New     STA   '',
     '                      ''  -  update    STA  '','' MG='',I4,
     '                      '' ME='',I4,
     '                      '' Previous STA ='',E10.4,'' TA='',E10.4)')
     '                      MG,ME,STA(MG,ME),TA(NG,NE,NNEIG)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ELSE
                          WRITE(OP_STRING,'(''    Old     STA   '',
     '                      ''  -  unchanged STA  '','' MG='',I4,
     '                      '' ME='',I4,
     '                      '' Previous STA ='',E10.4,'' TA='',E10.4)')
     '                      MG,ME,STA(MG,ME),TA(NG,NE,NNEIG)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ELSE
                        WRITE(OP_STRING,'(''    First   STA   '',
     '                    ''  -  update    STA  '','' MG='',I4,
     '                    '' ME='',I4,
     '                    '' Previous STA ='',E10.4,'' TA='',E10.4)')
     '                    MG,ME,STA(MG,ME),TA(NG,NE,NNEIG)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999) 
                      ENDIF
                    ENDIF                      
	  	    STA(MG,ME)=TA(NG,NE,NNEIG)
		  ENDIF
                  IF(NSA.LT.MSA)THEN
		    STIME=SSTA(NSA+1)
                  ELSE
		    STIME=SSTA(MSA)
                  ENDIF
 	          IF(STA(MG,ME).LE.STIME)THEN
		    AFLAG=.FALSE.
C ***               Don't activate it if it is to be activated from history
C ***               i.e. on epicardium with known activation times
		    DO NH=1,NHIS
		      IF(MG.EQ.NGHIS(NH).AND.ME.EQ.NEHIS(NH))THEN
     			AFLAG=.TRUE.
		        IF(DOP) THEN
                           WRITE(OP_STRING,
     '                       '(''   Dont activate now From history'')')
                           CALL WRITES(IODI,OP_STRING,ERROR,*9999) 
                        ENDIF  
	   	      ENDIF
                    ENDDO
C ***               Check to see if already activated this iteration
		    DO NA=1,NSA
		      IF(MG.EQ.NGA(NA).AND.ME.EQ.NEA(NA))THEN
		        IF(STA(MG,ME).LT.SSTA(NA))THEN
                          IF(DOP) THEN
                            WRITE(OP_STRING,'(''    Earlier STA   '',
     '                        ''  -  update    SSTA '','' MG='',I4,
     '                        '' ME='',I4,
     '                        '' Previous SSTA='',E10.4,
     '                        '' TA='',E10.4)')
     '                        MG,ME,SSTA(NA),TA(NG,NE,NNEIG)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                          SSTA(NA)=STA(MG,ME)
                          MEA(NA)=NE
                          MGA(NA)=NG
                          AFLAG=.TRUE.
                          CALL SORTAC2(NSA,SSTA,NGA,NEA,MGA,MEA)
                        ELSE
                          AFLAG=.TRUE. !For case of identical TA from another point
	                  IF(DOP) THEN
                            WRITE(OP_STRING,
     '                      '(''    Dont activate now. '',
     '                      '' Identical TA from elsewhere '')')
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
		        ENDIF
		      ENDIF
                    ENDDO
		    IF(DABS(STA(MG,ME)).LT.0.0001d0)THEN
		      AFLAG=.TRUE.
		      IF(DOP) THEN                
                        WRITE(OP_STRING,
     '                    '(''    Dont activate now. '',
     '                   '' STA equals 0. MG='',I5,'' ME='',I3)')MG,ME
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
		    ENDIF
		    IF(.NOT.AFLAG)THEN   !If no reason not to add this point to SSTA list
		      IF(NSA.EQ.MSA)THEN
		        SSTA(MSA)=STA(MG,ME)
		        NGA(MSA)=MG
		        NEA(MSA)=ME
		        MGA(MSA)=NG
		        MEA(MSA)=NE
		        CALL SORTAC1(MSA,SSTA,NGA,NEA,MGA,MEA,ERROR,
     '                    *9999)
		      ELSE
		        NSA=NSA+1
		        SSTA(NSA)=STA(MG,ME)
		        NGA(NSA)=MG
		        NEA(NSA)=ME
		        MGA(NSA)=NG
		        MEA(NSA)=NE
		        CALL SORTAC1(NSA,SSTA,NGA,NEA,MGA,MEA,ERROR,
     '                    *9999)
	 	      ENDIF
		      IF(DOP) THEN
                         WRITE(OP_STRING,
     '                   '(''    *Gauss Point added to SSTA*       '',
     '                   '' MG='',I4,'' ME='',I4,'' STA='',E10.4,
     '                   '' STIME='',
     '                   E10.4,''    NSA='',I3)')MG,ME,STA(MG,ME),
     '                   STIME,NSA
                         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF                                   
		    ENDIF
                  ELSE
	            IF(DOP) THEN 
                      WRITE(OP_STRING,'(''    Dont activate. '',
     '                '' Longer than max. SSTA this itn. & NSA=MSA '')')
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
		  ENDIF
                ELSE
		  IF(DOP) THEN
                    WRITE(OP_STRING,'(''    Dont activate now. '',
     '                '' There is an earlier TA from another point.'')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)       
                  ENDIF
                ENDIF
              ELSE
		IF(ITHRES(1,MG,ME).GE.1.
     '            AND.(.NOT.(MG.EQ.NG.AND.ME.EQ.NE)))THEN
		  ISUM=ISUM+1
		ENDIF
	      ENDIF               
	    ELSE
	      IF(ME.EQ.0)THEN
	        ISUM=ISUM+1
	      ENDIF
            ENDIF
 300      CONTINUE
C ***     If all neighbours are active, delete it from active list
	  IF(ISUM.EQ.NSE2) THEN
	    ITHRES(1,NG,NE)=2
	    EQUALS=.FALSE.
	    DO NA=1,NUAT-1
	      IF(INE(NA).EQ.NE.AND.ING(NA).EQ.NG) EQUALS=.TRUE.
	      IF(EQUALS)THEN
		INE(NA)=INE(NA+1)
		ING(NA)=ING(NA+1)
		IME(NA)=IME(NA+1)
		IMG(NA)=IMG(NA+1)
	      ENDIF
            ENDDO
	    NUAT=NUAT-1
 	    IF(NUAT.EQ.1)THEN
     	      WRITE(OP_STRING,'('' No more active points'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C	      NUAT=0
	      NSA=0
	    ENDIF
          ENDIF
 400    CONTINUE
 500  CONTINUE
C *** Search complete for this iteration,
C *** Sort Gauss points which want to activate this time round
C
C *** History Points  ? to start activation
      IF(NHIS.GT.0)THEN
        NREAC=0
        IF(NSA.EQ.0)THEN
          IF(KTYP31.EQ.1)THEN
	    IF(DABS(THIS(1)).LE.T)THEN
	      NSA=1
	      SSTA(NSA)=THIS(1)
	    ENDIF
          ELSE
	    IF(DABS(THIS(1)).GE.T)THEN
	      NSA=1
	      SSTA(NSA)=THIS(1)
	    ENDIF
          ENDIF
        ENDIF
C ***   Test history points to see if any are to be activated
        DO NH=1,NHIS   !NHIS points are in ascending order
          IF(KTYP31.EQ.1)THEN
            IF(THIS(NH).LE.SSTA(NSA)) THEN
	      NREAC=NREAC+1
	      NG=NGHIS(NH)
	      NE=NEHIS(NH)
C     	      WRITE(59,'('' activated from history : NG='',I3,'' NE='',I3,
C     '	        '' time='',e10.4)')NG,NE,THIS(NH)
 	      IF(DOP) THEN
                WRITE(OP_STRING,
     '          '('' activated from history : NG='',I3,'' NE='',I3,'
     '	        //''' time='',e10.4)') NG,NE,THIS(NH)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
	      SSTA(NSA)=THIS(NH)
	      NGA(NSA)=NG
	      NEA(NSA)=NE
	      MGA(NSA)=0
	      MEA(NSA)=0
	      CALL SORTAC1(NSA,SSTA,NGA,NEA,MGA,MEA,ERROR,*9999)
	    ELSE
	      GOTO 540
 	    ENDIF
          ELSE
	    IF(NSA.EQ.0)THEN
	      STIME=SSTA(1)
	    ELSE
	      STIME=SSTA(NSA)
	    ENDIF
            IF(DABS(THIS(NH)).GE.DABS(STIME)) THEN
	      NREAC=NREAC+1
	      NG=NGHIS(NH)
	      NE=NEHIS(NH)
C	      WRITE(59,549)NG,NE,THIS(NH)
 	      IF(DOP) THEN
                WRITE(OP_STRING,
     '         '('' activated from history : NG='',I3,'' NE='',I3,'
     '	        //''' time='',e10.4)') NG,NE,THIS(NH)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              SSTA(NSA)=THIS(NH)
              NGA(NSA)=NG
              NEA(NSA)=NE
	      MGA(NSA)=0
	      MEA(NSA)=0
	      CALL SORTAC1(NSA,SSTA,NGA,NEA,MGA,MEA,ERROR,*9999)
	    ELSE
	      GOTO 540
 	    ENDIF
          ENDIF
        ENDDO
 540    CONTINUE
        NHIS=NHIS-NREAC
        DO NH=1,NHIS
          NEHIS(NH)=NEHIS(NH+NREAC)
          NGHIS(NH)=NGHIS(NH+NREAC)
          THIS(NH)=THIS(NH+NREAC)
        ENDDO
        IF(DOP) THEN 
          WRITE(OP_STRING,'('' NHIS='',I4,'' NREAC='',I4)')
     '      NHIS,NREAC
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF
C *** Activate those on list to be activated
      WRITE(OP_STRING,'('' Points for activation this iteration'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C      WRITE(IO4,'('' Points for activation this iteration'')')
      DO NA=1,NSA
C	IF(DOP)THEN
	  WRITE(OP_STRING,'('' NG='',I3,'' NE='',I3,'' SSTA='',
     '           E10.4,'' ACTIVATED BY'',
     '           '' MG='',I3,'' ME='',I3)')
     '        NGA(NA),NEA(NA),SSTA(NA),MGA(NA),MEA(NA)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C	  WRITE(59,535)NGA(NA),NEA(NA),SSTA(NA),MGA(NA),MEA(NA)
C	ENDIF
        ITHRES(1,NGA(NA),NEA(NA))=1
        THRES(1,NGA(NA),NEA(NA))=SSTA(NA)
	ING(NUAT+NA)=NGA(NA)
	INE(NUAT+NA)=NEA(NA)
	IMG(NUAT+NA)=MGA(NA)
	IME(NUAT+NA)=MEA(NA)
      ENDDO
      NUAT=NUAT+NSA
      NUATT=NUATT+NSA
      IF(NSA.GT.0)THEN
        DT=DABS(SSTA(NSA))-T
      ELSE
	DT=TINCR
      ENDIF

C *** Create new list of active elements and gauss points
      NLET=0
      DO NA=1,NUAT
	EFLAG=.FALSE.
	DO NLE1=1,NLET
	  IF(INE(NA).EQ.LE(NLE1))THEN
	    EFLAG=.TRUE.
	    LG(NLE1,NLGT(NLE1)+1)=ING(NA)
	    NLGT(NLE1)=NLGT(NLE1)+1
	  ENDIF
        ENDDO
	IF(.NOT.EFLAG)THEN
	  NLET=NLET+1
	  LE(NLET)=INE(NA)
	  NLGT(NLET)=1
	  LG(NLET,1)=ING(NA)
	ENDIF
      ENDDO

      IF(DOP)THEN
        WRITE(OP_STRING,'(/,'' List of current active points'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
	DO NA=1,NUAT
	  NG=ING(NA)
	  NE=INE(NA)
	  MG=IMG(NA)
	  ME=IME(NA)
	  WRITE(OP_STRING,'('' NA='',I3,'' NG='',I3,'' NE='',I3,
     '       '' TIME='',E10.4,'' MG='',I3,'' ME='',I3)')NA,NG,NE,
     '       THRES(1,NG,NE),MG,ME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(/)')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
	DO NLE1=1,NLET
	  WRITE(OP_STRING,'('' LIST='',I3,'' ELEM='',I3)')
     '      NLE1,LE(NLE1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
	  DO NLG=1,NLGT(NLE1)
	    WRITE(OP_STRING,
     '        '(10X,'' LISTG='',I4,'' NG='',I4,'' NLE1='',
     '        I5)') NLG,LG(NLE1,NLG),NLE1
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      WRITE(OP_STRING,
     '  '(/,'' NUAT='',I8,'' NUATT='',I8,'' MSA='',I8)') NUAT,
     '  NUATT,MSA
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C     WRITE(59,'('' NUAT='',I8,'' NUATT='',I8,'' MSA='',I8)') NUAT,
C    '     NUATT,MSA

      CALL EXITS('FFRONT')
      RETURN
 9999 CALL ERRORS('FFRONT',ERROR)
      CALL EXITS('FFRONT')
      RETURN 1
      END


      SUBROUTINE FRONT(NBJ,NGAP,NKE,NQGE,NWQ,NXQ,CQ,DNUDXQ,DXDXIQ,
     '  GUQ,PROPQ,T,XA,XQ,ZA,ERROR,*)

C**** OLD ****
C****
C**** NO LONGER USED  ----  GBS  11 March 1994
C****
C**** ITYP3(nr)=3 constant wavespeed FitzHugh-Nagumo model
C**** Calculates activation pattern by using FitzHugh-Nagumo equations
C**** with constant diffusion term.
C**** The equation is modelled explicitly, with explicit finite
C**** differences for diffusion.
C**** ZA(1,1,1,nq) is membrane potential at global point nq
C**** ZA(1,2,1,nq) is recovery variable  at global point nq
C**** ZA(2,1,1,nq) is activation   time of first  pulse at nq
C**** ZA(2,2,1,nq) is deactivation time "    "      "    "  "
C**** ZA(2,3,1,nq) is u                 "    "      "    "  "
C**** ZA(2,4,1,nq) is activation   time of second pulse at nq
C****      ...
C**** ZA(3,1,1,nq) is time since global point became active
C**** ZA(3,3,1,nq) used for temporary storage of u at nq
C**** NWQ(1,nq) is 0 if nq is internal
C****    "          mq1 if nq on external bdy (for calc no-flux b.c.)
C**** NWQ(2,nq) is 0 if nq is internal
C****    "          mq2 if nq on external bdy
C**** NWQ(4,nq) is 0 if global point nq is not active
C****    "            1  "    "     "    "  " active
C****    "            2  "    "     "    " surrounded by active points
C**** KTYP32 = 1 for standard FHN equation (C2*V)
C****    "   = 2 for Rogers/McCulloch form of FHN (C2*U*V)
C**** UMQ contains the value of U at the nodes of the local quadratic
C**** element about nq.

C**** 7-SEP-1992
C**** Currently we are forming a local quadratic element about each
C**** internal g.p. nq, and using finite differences to calculate
C**** D*del2u.  This is currently approximated as D*u,ij*g(ij).
C**** Boundary g.p. require no-flux b.c's to be applied after computing
C**** the values of all internal points.  

C**** 14-SEP-1992
C**** Changing now to full anisotropic diffusion, using values in XQ
C**** and DXQ for solving equation.

C**** 15-SEP-1992
C**** Dij and Dij,k are now calculated in UPMATE and stored in PROPQ.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
!     Parameter List
      INTEGER NBJ(NJM,*),NGAP(NIM,*),NKE(NKM,NNM,NBM,*),NQGE(NGM,*),
     '  NWQ(6,0:*),NXQ(-NIM:NIM,0:*)
      REAL*8 CQ(NMM,*),DNUDXQ(3,3,*),DXDXIQ(3,3,*),GUQ(3,3,*),
     '	PROPQ(3,3,4,2,*),T,XA(NAM,NJM,*),XQ(NJM,*),ZA(NAM,NHM,NCM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IB,IJ,IK,J,JB,K,MQ,NA,NB,NE,NELEM,NI,NI1,NI2,
     '  NII,NIJ,NIK,NITB,NJ,NQ,NR
      REAL*8 A,B,C1,C2,D,DEL2U,DET,DI,DIFFU,DU,DUDX(3),D2UDX2(3,3),DV,
     '  F1,F2,GL(3,3),GU(3,3),UMQ(-1:1,-1:1,-1:1),SUM,SUM1,SUM2,TOTAL,
     '  U,U1,U2,UP,V,VP
      LOGICAL INTERNAL,SURROUND

C**** FitzHugh-Nagumo model
C**** U and V at current position
C**** DIFFU is diffusion term (D*del^2 u)
C**** F1 gives du/dt, F2 gives dv/dt

      F1(U,V,A,C1,C2,DIFFU)=DIFFU+C1*U*(U-A)*(1.0-U)
C *** Uses first term if ktyp32=1 (standard FHN)     
C ***  or second term if ktyp32=2 (Rogers/McCulloch FHN)
     '   -(2-KTYP32)*C2*V
     '   -(KTYP32-1)*C2*U*V
      F2(U,V,B,D)=B*(U-D*V)

      CALL ENTERS('FRONT',*9999)

      NR=1 !Need to update this maybe
      NB=NBJ(1,1)
      NITB=NIT(NB)

C *** For each g.p. nq, put the value of u to be used for this time
C *** step into ZA(3,3,1,nq) - temporary storage
      DO NQ=1,NQT   
        IF(ZA(2,3,1,NQ).GT.0.0.AND.
     '    T.GE.ZA(2,1,1,NQ).AND.T.LT.ZA(2,2,1,NQ)) THEN
          ZA(3,3,1,NQ)=ZA(2,3,1,NQ)  !Inside time for pulse 1 => use u1
        ELSE IF(ZA(2,6,1,NQ).GT.0.0.AND.
     '    T.GE.ZA(2,4,1,NQ).AND.T.LT.ZA(2,5,1,NQ)) THEN
          ZA(3,3,1,NQ)=ZA(2,6,1,NQ)  !Inside time for pulse 2 => use u2
        ELSE
          ZA(3,3,1,NQ)=ZA(1,1,1,NQ)  !u(nq)(t) = u(nq)(t-1) - previous u
        ENDIF
      ENDDO

C *** Loop now by grid point & solve for all internal g.p.
      DO NQ=1,NQT
      	if(dop) then
      	  WRITE(OP_STRING,'('' nq: '',I5)') nq
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      	endif

C *** If we have an internal g.p. continue
        IF(NWQ(1,NQ).EQ.0) THEN    !Not on boundary
      	  IK=MAX(0,NITB-2)         !zero for 1,2-D, one for 3-D
      	  IJ=MIN(NITB-1,1)         !zero for 1-D, one for 2,3-D
C *** Formulate a local quadratic element about nq defined by UMQ,
C *** storing the value of u at each node point mq.
      	  DO NIK=-IK,IK
      	    DO NIJ=-IJ,IJ
      	      DO NII=-1,1
C *** MQ is the neighbouring g.p.
      		MQ=NXQ(NII,NXQ(NIJ*2,NXQ(NIK*3,NQ)))
C *** Value of u at MQ
      		UMQ(NII,NIJ,NIK)=ZA(3,3,1,MQ)
      	      ENDDO
      	    ENDDO
      	  ENDDO

C *** Compute u,k by first order finite differences about nq.
          DUDX(1)=UMQ(1,0,0)-UMQ(-1,0,0)
          DUDX(2)=UMQ(0,1,0)-UMQ(0,-1,0)
          DUDX(3)=UMQ(0,0,1)-UMQ(0,0,-1)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' du/dx = '',3F10.5)')
     '        (dudx(ni),ni=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C *** Compute u,ij  = d2u/dxi(i)dxi(j) by taking second order finite
C *** differences about nq.
      	  DO IB=1,NITB
      	    DO JB=1,NITB
      	      IF(IB.EQ.JB) THEN
                IF(IB.EQ.1) THEN         ! u,11
      	      	  D2UDX2(IB,JB)=(UMQ(1,0,0)-2.0*UMQ(0,0,0)+
     '              UMQ(-1,0,0))*4.0
      	        ELSE IF(IB.EQ.2) THEN    ! u,22
      	       	  D2UDX2(IB,JB)=(UMQ(0,1,0)-2.0*UMQ(0,0,0)+
     '              UMQ(0,-1,0))*4.0
      	        ELSE IF(IB.EQ.3) THEN    ! u,33
      	    	  D2UDX2(IB,JB)=(UMQ(0,0,1)-2.0*UMQ(0,0,0)+
     '	     	    UMQ(0,0,-1))*4.0
      	        ENDIF
       	      ELSE
      	        IF(IB+JB.EQ.3) THEN       ! u,12 and u,21
      	 	  D2UDX2(IB,JB)=(UMQ(1,1,0)-UMQ(1,-1,0)-
     '	  	    UMQ(-1,1,0)+UMQ(-1,-1,0))
      	        ELSE IF(IB+JB.EQ.4) THEN  ! u,13 and u,31
      	    	  D2UDX2(IB,JB)=(UMQ(1,0,1)-UMQ(1,0,-1)-
     '	     	    UMQ(-1,0,1)+UMQ(-1,0,-1))
      	        ELSE IF(IB+JB.EQ.5) THEN  ! u,23 and u,32
      	       	  D2UDX2(IB,JB)=(UMQ(0,1,1)-UMQ(0,-1,1)-
     '	       	    UMQ(0,1,-1)+UMQ(0,-1,-1))
      	        ENDIF
      	      ENDIF
      	    ENDDO
      	  ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'('' d2u/dx2'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nj=1,nitb
      	      WRITE(OP_STRING,'(3F10.5)') (d2udx2(ni,nj),ni=1,nitb)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
      
C *** Compute diffusion term = (Dki,j*u,k+Dki*u,jk)*gij
          DIFFU=0.0
          DO I=1,NITB
            DO J=1,NITB
              SUM1=0.0
              SUM2=0.0
              DO K=1,NITB
                SUM1=SUM1+PROPQ(K,I,J+1,1,NQ)*DUDX(K)
                SUM2=SUM2+PROPQ(K,I,1,1,NQ)*D2UDX2(J,K)
              ENDDO
              DIFFU=DIFFU+(SUM1+SUM2)*GUQ(I,J,NQ)
            ENDDO
      	  ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Diffusion = '',F10.5)') diffu
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C *** FitzHugh-Nagumo model
      	  A =CQ(3,NQ)    !threshold potential
      	  B =CQ(4,NQ)    !recovery rate const
      	  C1=CQ(1,NQ)    !excitation rate const
      	  C2=CQ(2,NQ)    !excitation decay const
      	  D =CQ(5,NQ)    !recovery decay const
      	  U =ZA(3,3,1,NQ)  !is previous time step excitation var 
      	  V =ZA(1,2,1,NQ)  !is previous time step recovery variable
      	  IF(DOP) THEN
      	    WRITE(OP_STRING,
     '        '('' FHN constants (a,b,c1,c2,d) are: '','
     '	      //'5F12.6)') a,b,c1,c2,d
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      	  ENDIF

C *** Using Euler Predictor-Corrector method to solve for new u and v
      	  DU=F1(U,V,A,C1,C2,DIFFU)  !is initial slope for u
      	  DV=F2(U,V,B,D)            !is initial slope for v
      	  UP=U+DT*DU       !is Euler prediction for u
      	  VP=V+DT*DV       !is Euler prediction for v
      	  ZA(1,1,1,NQ)=U+DT*0.5*(DU+F1(UP,VP,A,C1,C2,DIFFU)) !u(nq,t)
      	  ZA(1,2,1,NQ)=V+DT*0.5*(DV+F2(UP,VP,B,D))           !v(nq,t)
      	  IF(DOP) THEN
      	    WRITE(OP_STRING,'('' du,dv,up,vp: '',4F12.6)')
     '        du,dv,up,vp
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      	    WRITE(OP_STRING,
     '        '('' Old u,v: '',2E14.5,'' New u,v: '',2E14.5)')
     '	      u,v,za(1,1,1,nq),za(1,2,1,nq)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      	    WRITE(OP_STRING,'('' ============='')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      	  endif
      	ENDIF   !Internal point
      ENDDO     !Global grid point

C *** Update bdy grid points
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Updating bdry/initial conditions'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO NQ=1,NQT
        IF(ZA(2,3,1,NQ).GT.0.0.AND.
     '    T.GE.ZA(2,1,1,NQ).AND.T.LT.ZA(2,2,1,NQ)) THEN
          ZA(1,1,1,NQ)=ZA(2,3,1,NQ)  !Inside time for pulse 1 => use u1
      	  IF(DOP) THEN
      	    WRITE(OP_STRING,'('' nq: '', I5,'' Defined u: '',F12.5)')
     '	      nq,za(1,1,1,nq)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      	  ENDIF
        ELSE IF(ZA(2,6,1,NQ).GT.0.0.AND.
     '    T.GE.ZA(2,4,1,NQ).AND.T.LT.ZA(2,5,1,NQ)) THEN
          ZA(1,1,1,NQ)=ZA(2,6,1,NQ)  !Inside time for pulse 2 => use u2
      	  IF(DOP) THEN
      	    WRITE(OP_STRING,'('' nq: '', I5,'' Defined u: '',F12.5)')
     '	      nq,za(1,1,1,nq)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      	  ENDIF
        ELSE
          IF(NWQ(1,NQ).NE.0) THEN    !On external bdy
C *** Generic no-flux bdy cond.n applied at external bdy point
      	    U1=ZA(1,1,1,NWQ(1,NQ))   !Pick up values for g.p. as
      	    U2=ZA(1,1,1,NWQ(2,NQ))   !defined in DEGRID
      	    ZA(1,1,1,NQ)=(4.0*U1-U2)/3.0
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' nq: '', I5,'' u1,u2,u: '',3F12.5)')
     '          nq,u1,u2,za(1,1,1,nq)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
      ENDDO

C *** Update flags on each grid point.
      DO NQ=1,NQT
        IF(NWQ(4,NQ).GT.0) THEN            !nq active already
          ZA(3,1,1,NQ)=ZA(3,1,1,NQ)+DABS(TINCR)  !Time since activation
          IF(NWQ(4,NQ).EQ.1) THEN  !Check if surrounded by active g.p.
            SURROUND=.TRUE.
            DO NIK=-1,1
              DO NIJ=-1,1
                DO NII=-1,1
                  MQ=NXQ(NII,NXQ(NIJ*2,NXQ(NIK*3,NQ)))
                  IF(MQ.NE.0.AND.NWQ(4,MQ).EQ.0) SURROUND=.FALSE.
                ENDDO
              ENDDO
            ENDDO
            IF(SURROUND) NWQ(4,NQ)=2
          ENDIF
        ELSE IF(ZA(1,1,1,NQ).GT.1.E-5) THEN !If U>0 then active
          NWQ(4,NQ)=1                     !Set state active
        ENDIF
        IF(ZA(1,1,1,NQ).LT.-1.E-5) THEN     !If U<0 then inactive
          NWQ(4,NQ)=0                     !Set state inactive
        ENDIF
      ENDDO

      CALL EXITS('FRONT')
      RETURN
 9999 CALL ERRORS('FRONT',ERROR)
      CALL EXITS('FRONT')
      RETURN 1
      END
      

Module FE40
=========== 


      SUBROUTINE ELAS3D(NBH,NBJ,nw,CG,
     '  ED,EM,ER,ES,PG,RG,SE,VE,WG,XE,XG,ERROR,*)

C#### Subroutine: ELAS3D
C###  Description:
C###    ELAS3D calculates the element load vector ER and the stiffness 
C###    matrix ES for 3D elasticity stress (nw=9).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:b15.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp40.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),nw
      REAL*8 CG(NMM,NGM),ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),
     '  ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM),VE(NSM,NKM),WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,m,ms,ms1,ms2,nb,ng,NGTB,n,nr,ns,
     '  ns1,ns2,ns_i,ns_j,ns_m,ns_n,NSTB
      REAL*8 CC(6,6),C1,C2,C3,DENS,DXIX(3,3),DXIXN(3,3),
     '  GL(3,3),GLN(3,3),GU(3,3),GUN(3,3),
     '  PM,PMX,PMY,PMZ,PN,PNX,PNY,PNZ,RT(3,3),RWG,S_thermal(3)
      LOGICAL ROTATE_FIBRE

      CALL ENTERS('ELAS3D',*9999)
      nr=1   !Temporary
      IF(DOP) THEN
        WRITE(OP_STRING,'('' 3D elasticity'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      nb=NBH(1) !is dependent variable basis
      NGTB=NGT(nb)
      NSTB=NST(nb)

      ROTATE_FIBRE=.FALSE.
      IF(IMT(nw).EQ.5) THEN !Special case anisotropic
        IF(NJ_LOC(NJL_FIBR,0).GT.0) THEN
          ROTATE_FIBRE=.TRUE.
        ENDIF
      ENDIF

      DO ng=1,NGTB
        CALL XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*9999)
        !get material coord info in DXIXN etc
        !if not fibre, then DXIXN==DXIX
        IF(ROTATE_FIBRE) THEN !Calculate wrt nu coords
          CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXIXN,GLN,GUN,
     '      RG(ng),XG,ERROR,*9999)
        ELSE
          CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIXN,GLN,GUN,
     '      RG(ng),XG,ERROR,*9999)
        ENDIF
        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
        RWG=RG(ng)*WG(ng,nb)
        CALL CGC9(nw,CC,C1,C2,C3,CG(1,ng),DENS,S_thermal,ERROR,*9999)
! Isotropic or Anisotropic (special case)
        IF(IMT(nw).EQ.1.OR.IMT(nw).EQ.5) THEN 
          IF(ROTATE_FIBRE) THEN !Calculate rotation matrix at gauss pt
            CALL MAT_VEC_NG(3,RT(1,1),RT(1,2),RT(1,3),XG,
     '        ERROR,*9999)
          ENDIF
          C3=0.25d0*C3 !since c3 always appears multiplied by 1/4
          DO ms=1,NSTB
            ms1=ms+NSTB
            ms2=ms+2*NSTB
            PM =PG(ms,1,ng,nb)
            PMX=PG(ms,2,ng,nb)*DXIXN(1,1)
     '         +PG(ms,4,ng,nb)*DXIXN(2,1)
     '         +PG(ms,7,ng,nb)*DXIXN(3,1)
            PMY=PG(ms,2,ng,nb)*DXIXN(1,2)
     '         +PG(ms,4,ng,nb)*DXIXN(2,2)
     '         +PG(ms,7,ng,nb)*DXIXN(3,2)
            PMZ=PG(ms,2,ng,nb)*DXIXN(1,3)
     '         +PG(ms,4,ng,nb)*DXIXN(2,3)
     '         +PG(ms,7,ng,nb)*DXIXN(3,3)
            !Thermal is only implemented for isotropic, so it should
            !be wrt element coords
            IF(KTYP43.GT.0) THEN !include thermal stress effects
C              ER(ms ) = ER(ms ) - S_thermal(1)*PMX*RWG 
C              ER(ms1) = ER(ms1) - S_thermal(2)*PMY*RWG 
C              ER(ms2) = ER(ms2) - S_thermal(3)*PMZ*RWG
C GMH 26/12/95 Use PMX etc from DXIX
              ER(ms ) = ER(ms ) - S_thermal(1)*
     '          (PG(ms,2,ng,nb)*DXIX(1,1)
     '          +PG(ms,4,ng,nb)*DXIX(2,1)
     '          +PG(ms,7,ng,nb)*DXIX(3,1))*RWG 
              ER(ms1) = ER(ms1) - S_thermal(2)*
     '          (PG(ms,2,ng,nb)*DXIX(1,2)
     '          +PG(ms,4,ng,nb)*DXIX(2,2)
     '          +PG(ms,7,ng,nb)*DXIX(3,2))*RWG 
              ER(ms2) = ER(ms2) - S_thermal(3)*
     '          (PG(ms,2,ng,nb)*DXIX(1,3)
     '          +PG(ms,4,ng,nb)*DXIX(2,3)
     '          +PG(ms,7,ng,nb)*DXIX(3,3))*RWG
            ENDIF !ktyp43
         
            DO ns=1,NSTB
              ns1=ns+NSTB
              ns2=ns+2*NSTB
              PN =PG(ns,1,ng,nb)
              PNX=PG(ns,2,ng,nb)*DXIXN(1,1)
     '           +PG(ns,4,ng,nb)*DXIXN(2,1)
     '           +PG(ns,7,ng,nb)*DXIXN(3,1)
              PNY=PG(ns,2,ng,nb)*DXIXN(1,2)
     '           +PG(ns,4,ng,nb)*DXIXN(2,2)
     '           +PG(ns,7,ng,nb)*DXIXN(3,2)
              PNZ=PG(ns,2,ng,nb)*DXIXN(1,3)
     '           +PG(ns,4,ng,nb)*DXIXN(2,3)
     '           +PG(ns,7,ng,nb)*DXIXN(3,3)
         
C GMH 26/12/95 Should this be wrt element or material coords???
              EM(ms ,ns) =EM(ms ,ns) +DENS*PM*PN*RWG
              EM(ms1,ns1)=EM(ms1,ns1)+DENS*PM*PN*RWG
              EM(ms2,ns2)=EM(ms2,ns2)+DENS*PM*PN*RWG
C GMH 22/12/95 We must rotate this component of the ES matrix by
C              the rotation matrix at this point.  We can store
C              the component in ED as it gets overwritten at the
C              end of this routine.
C              ES( ms,ns ) =ES(ms,ns)
C     '          +(C1*PMX*PNX+C3*(PMY*PNY+PMZ*PNZ))*RWG    !(u,u) term
C              ES( ms,ns1) =ES(ms,ns1)
C     '          +(C2*PMX*PNY+C3*PMY*PNX)*RWG              !(u,v) term
C              ES( ms,ns2) =ES(ms,ns2)
C     '          +(C2*PMX*PNZ+C3*PMZ*PNX)*RWG              !(u,w) term
C         
C              ES(ms1,ns ) =ES(ms1,ns)
C     '          +(C2*PMY*PNX+C3*PMX*PNY)*RWG              !(v,u) term
C              ES(ms1,ns1) =ES(ms1,ns1)
C     '          +(C1*PMY*PNY+C3*(PMZ*PNZ+PMX*PNX))*RWG    !(v,v) term
C              ES(ms1,ns2) =ES(ms1,ns2)
C     '          +(C2*PMY*PNZ+C3*PMZ*PNY)*RWG              !(v,w) term
C         
C              ES(ms2,ns ) =ES(ms2,ns)
C     '          +(C2*PMZ*PNX+C3*PMX*PNZ)*RWG              !(w,u) term
C              ES(ms2,ns1) =ES(ms2,ns1)
C     '          +(C2*PMZ*PNY+C3*PMY*PNZ)*RWG              !(w,v) term
C              ES(ms2,ns2) =ES(ms2,ns2)
C     '          +(C1*PMZ*PNZ+C3*(PMX*PNX+PMY*PNY))*RWG    !(w,w) term
              ED( ms,ns ) =
     '          (C1*PMX*PNX+C3*(PMY*PNY+PMZ*PNZ))*RWG    !(u,u) term
              ED( ms,ns1) =
     '          (C2*PMX*PNY+C3*PMY*PNX)*RWG              !(u,v) term
              ED( ms,ns2) =
     '          (C2*PMX*PNZ+C3*PMZ*PNX)*RWG              !(u,w) term
              ED(ms1,ns ) =
     '          (C2*PMY*PNX+C3*PMX*PNY)*RWG              !(v,u) term
              ED(ms1,ns1) =
     '          (C1*PMY*PNY+C3*(PMZ*PNZ+PMX*PNX))*RWG    !(v,v) term
              ED(ms1,ns2) =
     '          (C2*PMY*PNZ+C3*PMZ*PNY)*RWG              !(v,w) term 
              ED(ms2,ns ) =
     '          (C2*PMZ*PNX+C3*PMX*PNZ)*RWG              !(w,u) term
              ED(ms2,ns1) =
     '          (C2*PMZ*PNY+C3*PMY*PNZ)*RWG              !(w,v) term
              ED(ms2,ns2) =
     '          (C1*PMZ*PNZ+C3*(PMX*PNX+PMY*PNY))*RWG    !(w,w) term

            ENDDO !ns
          ENDDO !ms
          
C         Now rotate ED, and add to ES
          DO i=1,3
            ns_i=(i-1)*NSTB
            DO j=1,3
              ns_j=(j-1)*NSTB
              DO ms=1,NSTB
                DO ns=1,NSTB
                  IF(ROTATE_FIBRE) THEN
                    !ES=RT(transposed)*ES*RT
                    !This is the rotation of one block of ED
                    !The i,j th element of ES gets the m,n th 
                    !element of ED added to it.
                    DO m=1,3
                      ns_m=(m-1)*NSTB
                      DO n=1,3
                        ns_n=(n-1)*NSTB
                        ES(ns_i+ms,ns_j+ns)=ES(ns_i+ms,ns_j+ns)+
     '                    RT(i,m)*ED(ns_m+ms,ns_n+ns)*RT(j,n)
                      ENDDO !n
                    ENDDO !m
                  ELSE
                    ES(ns_i+ms,ns_j+ns)=ES(ns_i+ms,ns_j+ns)+
     '                ED(ns_i+ms,ns_j+ns)
                  ENDIF
                ENDDO !ns
              ENDDO !ms
            ENDDO !j
          ENDDO !i
! Transversly isotropic (fibre in 1 dir.n)
        ELSE IF(IMT(nw).EQ.2) THEN 

! Transversly isotropic (fibre in 2 dir.n)
        ELSE IF(IMT(nw).EQ.3) THEN 

! Orthotropic (fibre in 1 dir.n) 
        ELSE IF(IMT(nw).EQ.4) THEN

        ENDIF !anisotropy
      ENDDO !ng

      DO ms=1,NSTB
        ms1=ms+NSTB
        ms2=ms+2*NSTB
        ER(ms ) = ER(ms ) *SE(ms,nb)
        ER(ms1) = ER(ms1) *SE(ms,nb)
        ER(ms2) = ER(ms2) *SE(ms,nb)
        DO ns=1,NSTB
          ns1=ns+NSTB
          ns2=ns+2*NSTB
          ES(ms ,ns ) = ES(ms ,ns ) *SE(ms,nb)*SE(ns,nb)
          ES(ms ,ns1) = ES(ms ,ns1) *SE(ms,nb)*SE(ns,nb)
          ES(ms ,ns2) = ES(ms ,ns2) *SE(ms,nb)*SE(ns,nb)
          ES(ms1,ns ) = ES(ms1,ns ) *SE(ms,nb)*SE(ns,nb)
          ES(ms1,ns1) = ES(ms1,ns1) *SE(ms,nb)*SE(ns,nb)
          ES(ms1,ns2) = ES(ms1,ns2) *SE(ms,nb)*SE(ns,nb)
          ES(ms2,ns ) = ES(ms2,ns ) *SE(ms,nb)*SE(ns,nb)
          ES(ms2,ns1) = ES(ms2,ns1) *SE(ms,nb)*SE(ns,nb)
          ES(ms2,ns2) = ES(ms2,ns2) *SE(ms,nb)*SE(ns,nb)
        ENDDO !ns
      ENDDO !ms

C GMH 22/12/95 Code above assumes that ED is being overwritten -
C              do not rely on ED being zero.
      IF(DAMPING_FACTOR1.GT.0.d0.OR.DAMPING_FACTOR2.GT.0.d0) THEN
        DO ms=1,2*NSTB
          DO ns=1,2*NSTB
            ED(ms,ns)=DAMPING_FACTOR1*DAMPING_FREQUENCY*EM(ms,ns)
     '               +DAMPING_FACTOR2*DAMPING_FREQUENCY*ES(ms,ns)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('ELAS3D')
      RETURN
 9999 CALL ERRORS('ELAS3D',ERROR)
      CALL EXITS('ELAS3D')
      RETURN 1
      END



      SUBROUTINE IPINI4(IBT,IDO,INP,NBH,NBJ,NEELEM,
     '  NGAP,NHE,NHP,NJE,NJP,NKH,NKJ,NP_INTERFACE,
     '  NPNE,NPNODE,NPNY,NQE,nr,NVHP,NVJE,NW,nx,NXI,NYNE,NYNP,NYNR,
     '  XP,YP,ZA,ZP,FIX,ERROR,*)

C#### Subroutine: IPINI4
C###  Description:
C###    Inputs initial conditions and boundary conditions
C**** YP(ny,1) contains essential b.c.s,defined by FIX(ny,1).
C**** YP(ny,2)    "     natural           "        FIX(ny,2).
C**** YP(ny,4)    "     initial solution.
C**** Note: Load params defined by elements or nodes are assumed to be
C****       aligned with the element or global coord system,respect.ly,
C****       (in the latter case a rotation is made for beams & plates).
C**** Note: Direction cosines of normals to principal axes of beam
C****       elements in 3D space are carried in CE(ILT(2)+1...,ne)
C****       but  are not interpolated at the element Gauss points.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:dx00.cmn'
      INCLUDE 'cmiss$reference:fluid00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:titl40.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),
     '  NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NEELEM(0:NEM,0:*),NGAP(NIM,*),NHE(*),
     '  NHP(*),NJE(*),NJP(*),NKH(NHM,NPM,*),
     '  NKJ(NJM,*),
     '  NP_INTERFACE(0:NPM,0:*),NPNE(NNM,NBFM,*),
     '  NPNODE(0:NPM,0:*),NPNY(0:6,NYM,0:*),NQE(NSM,NBFM,*),nr,
     '  NVHP(NHM,NPM,*),NVJE(NNM,NBFM,NJM,*),
     '  NW(NEM,*),nx,NXI(-NIM:NIM,0:*),
     '  NYNE(NAM,NHM,0:NRCM,NCM,*),NYNP(NKM,NVM,NHM,NPM,0:NRCM,*),
     '  NYNR(0:NYM,0:NRCM,*)
      REAL*8 XP(NKM,NVM,NJM,*),YP(NYM,*),
     '  ZA(NAM,NHM,NCM,*),ZP(NKM,NVM,NHM,NPM,*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*)
!     Local Variables
      INTEGER GETNYP,ICHAR,INFO,loop,LOOPT,
     '  mh,N1CHAR,N2CHAR,nb,ne,nh,nk,noelem,
     '  nonode,NOQUES,np,nv,ny,NYTOT
      CHARACTER CFROMI*2,CHAR1*2
      LOGICAL FILEIP

      CALL ENTERS('IPINI4',*9999)
      nv=1 ! temporary cpb 22/11/94
      FILEIP=.FALSE.
      NOQUES=0

C *** CE(14,ne) used for total pressure load on membranes
C     DO noelem=1,NEELEM(0,nr)
C       ne=NEELEM(noelem,nr)
C       CE(14,ne)=0.d0
C     ENDDO

      IF(ETYP(8).OR.ETYP(10)) THEN !fluid-coupled
        FORMAT='($,'' Specify the fluid mass density'//
     '    ' [1000.0]: '',E12.5)'
        RDEFLT(1)=1000.0d0
        IF(IOTYPE.EQ.3) RDATA(1)=RHOL
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) RHOL=RDATA(1)

        FORMAT='($,'' Specify the g-acceleration [9.81]: '',E12.5)'
        RDEFLT(1)=9.81d0
        IF(IOTYPE.EQ.3) RDATA(1)=GACN
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GACN=RDATA(1)
      ENDIF

      IF(ETYP(7).OR.ETYP(8)) THEN !shell
        !To determine the 3rd basis functions (GEOM) derivs for shells.
        !Only the hermite basis is considered, since hermite is used for
        !shell geometry.
        nb=NBJ(1,NEELEM(1,nr))
        CALL ASSERT(NKT(0,nb).EQ.4,'>>Geometry must be cubic Hermite',
     '    ERROR,*9999)
        CALL GAUS20(IDO(1,0,nb),INP(1,1,nb),
     '    nb,NGAP(1,nb),D3PG(1,1,1,nb),ERROR,*9999)
      ENDIF

C     Check basis of beam element variables
      IF(ETYP(3)) THEN !some beam elements
        !Note: geom basis is assumed linear & axial displ. is cubic 
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).EQ.3) THEN !beam element
            CALL ASSERT(NKT(0,NBJ(1,ne)).EQ.1,
     '        ' >>Geometry basis sb linear',ERROR,*9999)
            nb=NBH(1,1,ne)
            IF(nb.GT.0) CALL ASSERT(NKT(0,nb).EQ.2,
     '        ' >>Axial displ. basis sb cubic',ERROR,*9999)
            nb=NBH(2,1,ne)
            IF(nb.GT.0) CALL ASSERT(NKT(0,nb).EQ.2,
     '        ' >>Transverse displ. basis sb cubic',ERROR,*9999)
            nb=NBH(3,1,ne)
            IF(nb.GT.0) CALL ASSERT(NKT(0,nb).EQ.2,
     '        ' >>Transverse displ. basis sb cubic',ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO nh=1,NHE(ne)
          IF(NBH(nh,1,ne).EQ.0) THEN
            DO mh=1,NHE(ne)
              IF(NBH(mh,1,ne).NE.0) THEN
                NBH(nh,1,ne)=NBH(mh,1,ne)
                GO TO 39
              ENDIF
            ENDDO
 39         CONTINUE
          ENDIF
        ENDDO
      ENDDO

c     DO nonode=1,NPNODE(0,nr)
c       np=NPNODE(nonode,nr)
c       NHP(np)=0
c       DO nh=1,NHM
c         NKH(nh,np,1)=0
c       ENDDO
c       DO noelem=1,NEELEM(0,nr)
c         ne=NEELEM(noelem,nr)
c         nb=NBH(1,1,ne)
c         ie=NW(ne,1)
c         DO nn=1,NNT(nb)
c           IF(NPNE(nn,nb,ne).EQ.NP) THEN
c             IF(NHE(ne).GT.NHP(np)) NHP(np)=NHE(ne)
c             DO nvar=1,NVE(ie)
c               nh=NHV(nvar,ie)
c               nb=NBH(nh,1,ne)
c               IF(NKT(0,nb).GT.NKH(nh,np,1)) NKH(nh,np,1)=NKT(0,nb)
c             ENDDO
c             GO TO 44
c           ENDIF
c         ENDDO
c44       CONTINUE
c       ENDDO
c       DO nh=1,NHP(np)
c         IF(NKH(nh,np,1).EQ.0) THEN
c           DO mh=1,NHP(np)
c             IF(NKH(mh,np,1).NE.0) THEN
c               NKH(nh,np,1)=NKH(mh,np,1)
c               GO TO 452
c             ENDIF
c           ENDDO
c452        CONTINUE
c         ENDIF
c       ENDDO
c     ENDDO

      IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.4.OR.
     '  (ITYP5(nr,nx).EQ.1.AND.ITYP6(nr,nx).EQ.2)) THEN
 590    FORMAT='('' Specify whether initial solution is [1]:'''//
     '    '/''   (1) Zero or undeformed (if nonlinear)'''//
     '    '/''   (2) Read in'''//
     '    '/''   (3) Restart from previous solution'''//
     '    '/,$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP5
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP5=IDATA(1)
        IF(KTYP5.EQ.3) THEN
          RESTAR=.TRUE.
        ELSE
          RESTAR=.FALSE.
        ENDIF
      ELSE
        KTYP5=1
      ENDIF

      IF(KTYP5.EQ.1) THEN
        DO ny=1,NYT(2,1,nx) !should use NYNR here
          YP(ny,4)=0.d0
          YP(ny,6)=0.d0
          YP(ny,7)=0.d0
        ENDDO
        IF(ITYP6(nr,nx).EQ.2) THEN
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nh=1,NJP(np)
              DO nv=1,NVHP(nh,np,1)
                DO nk=1,NKH(nh,np,1)
                  ZP(nk,nv,nh,np,1)=XP(nk,nv,nh,np)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL ZPYP(4,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '      nr,NVHP,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' YP(ny,4): '',10E11.3)')
     '        (YP(ny,4),NY=1,NYT(2,1,nx)) !should use NYNR here
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ELSE IF(KTYP5.EQ.3) THEN
        !Restart from previous solution
      ENDIF

      IF(KTYP5.EQ.2) THEN
        LOOPT=4
      ELSE IF(ITYP6(nr,nx).EQ.1) THEN !linear
        LOOPT=2
      ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
        LOOPT=3
      ENDIF
      DO loop=1,LOOPT
        nonode=0
        IF(loop.EQ.1) FORMAT='('' Essential boundary conditions'//
     '    ' defined at nodes:'')'
        IF(loop.EQ.2) FORMAT='(/'' Force boundary conditions defined'//
     '    ' at nodes:'')'
        IF(loop.EQ.3) FORMAT='(/'' Incremented boundary conditions:'')'
        IF(loop.EQ.4) FORMAT='(/'' Initial conditions defined'//
     '    ' at nodes:'')'
        CALL GINOUT(IOTYPE,0,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
 6100   FORMAT='(/$,'' Enter node number [EXIT]: '',I3)'
        IF(IOTYPE.EQ.3)THEN
          nonode=nonode+1
          IF(nonode.LE.NPNODE(0,nr)) THEN
            np=NPNODE(nonode,nr)
            IDATA(1)=NP
          ELSE
            IDATA(1)=0
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NPM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

        IF(IDATA(1).NE.0) THEN
          IF(IOTYPE.NE.3)THEN
            np=IDATA(1)
          ENDIF
          DO nh=1,NHP(np)
            CHAR1=CFROMI(nh,'(I1)')
            FORMAT='('' Dependent variable number '//CHAR1(1:1)//' :'')'
            CALL GINOUT(IOTYPE,0,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '        ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            DO nv=1,NVHP(nh,np,1)
              IF(NVHP(nh,np,1).GT.1) THEN
                CHAR1=CFROMI(nv,'(I2)')
                FORMAT='('' For version number '//CHAR1(1:2)//':'')'
                CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
              ENDIF
              DO nk=1,NKH(nh,np,1)
                RDEFLT(1)=1.1d6
                NYTOT=GETNYP(nk,nv,nh,np,0,1,nr,NP_INTERFACE,NPNY,
     '            NYNP,NYNR)
                IF(nk.EQ.1) THEN
                  FORMAT='($,'' The value of the dependent variable is '
     '              //'[no b.c.]: '',G12.5)'
                ELSE IF(nk.GT.1) THEN
                  CHAR1=CFROMI(nk,'(I1)')
                  FORMAT='($,'' The value of derivative number '//
     '              CHAR1(1:1)//' is [no b.c.]: '',G12.5)'
                ENDIF
                IF(IOTYPE.EQ.3) THEN
                  IF(loop.LE.3) THEN
                    IF(FIX(NYTOT,loop))THEN
                      RDATA(1)=YP(NYTOT,loop)
                    ELSE
                      RDATA(1)=RDEFLT(1)
                    ENDIF
                  ELSE IF(loop.EQ.4) THEN
                    IF(FIX(NYTOT,loop))THEN
                      RDATA(1)=YP(NYTOT,4)
                    ELSE
                      RDATA(1)=RDEFLT(1)
                    ENDIF
                  ENDIF
                ENDIF
                CALL GINOUT(IOTYPE,5,IVDU,IFILE,N1CHAR,N2CHAR,NOQUES,
     '            FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '            ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '            LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                IF(RDATA(1).LE.1.E6) THEN
                  IF(IOTYPE.NE.3) THEN
                    IF(loop.LE.3) THEN
                      YP(NYTOT,loop)=RDATA(1)
                      FIX(NYTOT,loop)=.TRUE.
                    ELSE IF(loop.EQ.4) THEN
                      YP(NYTOT,4)=RDATA(1)
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO  ! nk
            ENDDO  ! nv
          ENDDO  ! nh
          GO TO 6100
        ENDIF
      ENDDO

      CALL EXITS('IPINI4')
      RETURN
 9999 CALL ERRORS('IPINI4',ERROR)
      CALL EXITS('IPINI4')
      RETURN 1
      END

      SUBROUTINE OPINI4(NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NW,nx,
     '  YP,FIX,ERROR,*)

C#### Subroutine: OPINI4
C###  Description:
C###    Output of initial and boundary data.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:titl40.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,*),NEELEM(0:*),NHE(*),NHP(*),NKH(NHM,NPM,*),
     '  NPNODE(0:*),nr,NW(NEM,*),nx
      REAL*8 YP(NYM,*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*)
!     Local Variables
      INTEGER IBEG,IEND,nc,nh,nhx,NHK,NHKT,nk,nonode,np,NY
      CHARACTER CFROMI*2,CHAR*2,FORMAT*200

      CALL ENTERS('OPINI4',*9999)
      NY=0
      nc=1 !Temporary AJP 18-12-91

      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        NHKT=0
        DO nhx=1,NHP(np)
          nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
          DO nk=1,NKH(nh,np,nc)
            NHKT=NHKT+1
          ENDDO
        ENDDO
        CHAR=CFROMI(NHKT,'(I2)')
        CALL TRIM(CHAR,IBEG,IEND)
        FORMAT='('' Node'',I4,'' Essential   b.c.s: '','
     '    //CHAR(IBEG:IEND)//'L1,9E12.4,/(33X,9E12.4))'
        WRITE(OP_STRING,FORMAT) np,(FIX(NHK+ny,1),NHK=1,NHKT),
     '    (YP(NHK+ny,1),NHK=1,NHKT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(10X,''Natural     b.c.s: '','
     '    //CHAR(IBEG:IEND)//'L1,9E12.4,/(33X,9E12.4))'
        WRITE(OP_STRING,FORMAT) (FIX(NHK+ny,2),NHK=1,NHKT),
     '    (YP(NHK+ny,2),NHK=1,NHKT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ITYP6(nr,nx).EQ.2) THEN
          FORMAT='(10X,''Incremental b.c.s: '','
     '    //CHAR(IBEG:IEND)//'L1,9E12.4,/(33X,9E12.4))'
          WRITE(OP_STRING,FORMAT) (FIX(NHK+ny,3),NHK=1,NHKT),
     '      (YP(NHK+ny,3),NHK=1,NHKT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        FORMAT='(10X,''Initial conditions:'','
     '    //CHAR(IBEG:IEND)//'X,9E12.4,/(33X,9E12.4))'
        WRITE(OP_STRING,FORMAT) (YP(NHK+ny,4),NHK=1,NHKT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ny=ny+NHKT
      ENDDO

      WRITE(OP_STRING,'(/12X,''   No. of       No. of  '','
     '  //'/12X,''  variables   derivs/var'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        WRITE(OP_STRING,'(3X,''Node'',I4,8X,I1,7X,6I3)')
     '    np,NHP(np),(NKH(nh,np,nc),nh=1,NHP(np))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      CALL EXITS('OPINI4')
      RETURN           
 9999 CALL ERRORS('OPINI4',ERROR)
      CALL EXITS('OPINI4')
      RETURN 1
      END


Module FE50
=========== 

C 25/2/97 LC removed section from :
C
C#### Subroutine: D_PFRF
C###  Description:
C###    D_PFRF evaluates derivatives wrt geometric variables of 
C###    contributions to element residuals RE from the pressure PF(i) 
C###    acting on the Xi(3)=IXF face (IXF=iface-1).
C
C old MPN 28-Jun-1995: integrals now wrt deformed coords
C        IF(ITYP10(nr).EQ.1) THEN
C          DO nhx=1,NJ_LOC(NJL_GEOM,0)
C            nh=NH_LOC(nhx,nx)
C            DO ni=1,NITB
C              DZI(nhx,ni)=ZG(nhx,NU1(ni))
C              D_DZI(nhx,ni)=D_ZG(nhx,NU1(ni))
C            ENDDO
C          ENDDO
C        ELSE IF(ITYP10(nr).EQ.2) THEN
C          RX=XG(1,1)
C          RZ=ZG(1,1)
C          D_RZ=D_ZG(1,1)
C          DT=ZG(2,1)-XG(2,1)
C          D_DT=D_ZG(2,1)
C          CT=DCOS(DT)
C          D_CT=-DSIN(DT)*D_DT
C          ST=DSIN(DT)
C          D_ST= DCOS(DT)*D_DT
C          DO ni=1,NITB
C            DZI(1,ni)= ZG(1,NU1(ni))*CT-RZ*ST*ZG(2,NU1(ni))
C            D_DZI(1,ni)= D_ZG(1,NU1(ni))*CT +ZG(1,NU1(ni))*D_CT
C     '        -D_RZ*ST*ZG(2,NU1(ni)) -RZ*D_ST*ZG(2,NU1(ni)) 
C     '                               -RZ*ST*D_ZG(2,NU1(ni))
C            DZI(2,ni)=(ZG(1,NU1(ni))*ST+RZ*CT*ZG(2,NU1(ni)))/RX
C            D_DZI(2,ni)=(D_ZG(1,NU1(ni))*ST +ZG(1,NU1(ni))*D_ST
C     '        +D_RZ*CT*ZG(2,NU1(ni)) +RZ*D_CT*ZG(2,NU1(ni)) 
C     '                               +RZ*CT*D_ZG(2,NU1(ni)))/RX
C            DZI(3,ni)= ZG(3,NU1(ni))
C            D_DZI(3,ni)= D_ZG(3,NU1(ni))
C          ENDDO
C        ELSE IF(ITYP10(nr).EQ.3) THEN
C          RX=XG(1,1)
C          RZ=ZG(1,1)
C          D_RZ=D_ZG(1,1)
C          DT=ZG(2,1)-XG(2,1)
C          D_DT=D_ZG(2,1)
C          CT=DCOS(DT)
C          D_CT=-DSIN(DT)*D_DT
C          ST=DSIN(DT)
C          D_ST= DCOS(DT)*D_DT
C          CX=DCOS(XG(3,1))
C          SX=DSIN(XG(3,1))
C          CZ=DCOS(ZG(3,1))
C          D_CZ=-DSIN(ZG(3,1))*D_ZG(3,1)
C          SZ=DSIN(ZG(3,1))
C          D_SZ= DCOS(ZG(3,1))*D_ZG(3,1)
C          CC=CX*CZ
C          D_CC=CX*D_CZ
C          SS=SX*SZ
C          D_SS=SX*D_SZ
C          CS=CX*SZ
C          D_CS=CX*D_SZ
C          SC=SX*CZ
C          D_SC=SX*D_CZ
C          DO ni=1,NITB
C            RB=ZG(1,NU1(ni))
C            D_RB=D_ZG(1,NU1(ni))
C            TB=ZG(2,NU1(ni))
C            D_TB=D_ZG(2,NU1(ni))
C            PB=ZG(3,NU1(ni))
C            D_PB=D_ZG(3,NU1(ni))
C            DZI(1,ni)=CC*(RB*CT-RZ*ST*TB)-CS*RZ*CT*PB+SC*RZ*PB+SS*RB
C            D_DZI(1,ni)=D_CC*RB*CT +CC*D_RB*CT +CC*RB*D_CT
C     '        -D_CC*RZ*ST*TB-CC*D_RZ*ST*TB-CC*RZ*D_ST*TB-CC*RZ*ST*D_TB 
C     '        -D_CS*RZ*CT*PB-CS*D_RZ*CT*PB-CS*RZ*D_CT*PB-CS*RZ*CT*D_PB
C     '        +D_SC*RZ*PB+SC*D_RZ*PB+SC*RZ*D_PB  +  D_SS*RB+SS*D_RB 
C            DZI(2,ni)=(RZ*CZ*CT*TB+(RB*CZ-RZ*SZ*PB)*ST)/(RX*CX)
C            D_DZI(2,ni)=(D_RZ*CZ*CT*TB +RZ*D_CZ*CT*TB +RZ*CZ*D_CT*TB
C     '                                                +RZ*CZ*CT*D_TB
C     '        +D_RB*CZ*ST +RB*D_CZ*ST +RB*CZ*D_ST 
C     '        -D_RZ*SZ*PB*ST -RZ*D_SZ*PB*ST -RZ*SZ*D_PB*ST
C     '                                      -RZ*SZ*PB*D_ST)/(RX*CX)
C            DZI(3,ni)=(CC*RZ*PB+CS*RB+SC*(RZ*ST*TB-RB*CT)+SS*RZ*CT*PB)
C     '                                                 /RX
C            D_DZI(3,ni)=(D_CC*RZ*PB +CC*D_RZ*PB +CC*RZ*D_PB 
C     '        +D_CS*RB +CS*D_RB 
C     '        +D_SC*RZ*ST*TB+SC*D_RZ*ST*TB+SC*RZ*D_ST*TB+SC*RZ*ST*D_TB 
C     '        -D_SC*RB*CT -SC*D_RB*CT -SC*RB*D_CT 
C     '        +D_SS*RZ*CT*PB+SS*D_RZ*CT*PB+SS*RZ*D_CT*PB+SS*RZ*CT*D_PB)
C     '                                                 /RX
C          ENDDO
C        ELSE IF(ITYP10(nr).EQ.4) THEN
C          SLX=DSINH(XG(1,1))
C          SLZ=DSINH(ZG(1,1))
C          D_SLZ=DCOSH(ZG(1,1))*D_ZG(1,1)
C          SMX=DSIN (XG(2,1))
C          SMZ=DSIN (ZG(2,1))
C          D_SMZ=DCOS (ZG(2,1))*D_ZG(2,1)
C          CLX=DSQRT(1.0d0+SLX*SLX)
C          CLZ=DSQRT(1.0d0+SLZ*SLZ)
C          D_CLZ= SLZ*D_SLZ/DSQRT(1.0d0+SLZ*SLZ)
C          CMX=DSQRT(1.0d0-SMX*SMX)
C          CMZ=DSQRT(1.0d0-SMZ*SMZ)
C          D_CMZ=-SMZ*D_SMZ/DSQRT(1.0d0-SMZ*SMZ)
C          CSLX=CLX/SLX
C          CSMX=CMX/SMX
C          DT=ZG(3,1)-XG(3,1)
C          D_DT=D_ZG(3,1)
C          CT=DCOS(DT)
C          D_CT=-DSIN(DT)*D_DT
C          ST=DSIN(DT)
C          D_ST= DCOS(DT)*D_DT
C          CCL=CLX*CLZ
C          D_CCL=CLX*D_CLZ
C          CSL=CLX*SLZ
C          D_CSL=CLX*D_SLZ
C          SCL=SLX*CLZ
C          D_SCL=SLX*D_CLZ
C          SSL=SLX*SLZ
C          D_SSL=SLX*D_SLZ
C          CC=CMX*CMZ
C          D_CC=CMX*D_CMZ
C          CS=CMX*SMZ
C          D_CS=CMX*D_SMZ
C          SC=SMX*CMZ
C          D_SC=SMX*D_CMZ
C          SS=SMX*SMZ
C          D_SS=SMX*D_SMZ
C          G1=SLX*SLX+SMX*SMX
C          G3=SLX*SLX*SMX*SMX
C          DO ni=1,NITB
C            DLB=ZG(1,NU1(ni))
C            D_DLB=D_ZG(1,NU1(ni))
C            DMB=ZG(2,NU1(ni))
C            D_DMB=D_ZG(2,NU1(ni))
C            DTB=ZG(3,NU1(ni))
C            D_DTB=D_ZG(3,NU1(ni))
C            DZI(1,ni)=(( SSL*CC+CCL*SS*CT)*DLB+(-SCL*CS+CSL*SC*CT)*DMB
C     '                                              -CSL*SS*ST*DTB)/G1
C            D_DZI(1,ni)=(D_SSL*CC*DLB +SSL*D_CC*DLB +SSL*CC*D_DLB
C     '        +D_CCL*SS*CT*DLB +CCL*D_SS*CT*DLB +CCL*SS*D_CT*DLB
C     '                                          +CCL*SS*CT*D_DLB 
C     '        -D_SCL*CS*DMB -SCL*D_CS*DMB -SCL*CS*D_DMB 
C     '        +D_CSL*SC*CT*DMB +CSL*D_SC*CT*DMB +CSL*SC*D_CT*DMB 
C     '                                          +CSL*SC*CT*D_DMB 
C     '        -D_CSL*SS*ST*DTB -CSL*D_SS*ST*DTB -CSL*SS*D_ST*DTB 
C     '                                          -CSL*SS*ST*D_DTB)/G1
C            DZI(2,ni)=((-CSL*SC+SCL*CS*CT)*DLB+( CCL*SS+SSL*CC*CT)*DMB
C     '                                              -SSL*CS*ST*DTB)/G1
C            D_DZI(2,ni)=(-D_CSL*SC*DLB -CSL*D_SC*DLB -CSL*SC*D_DLB
C     '        +D_SCL*CS*CT*DLB +SCL*D_CS*CT*DLB +SCL*CS*D_CT*DLB
C     '                                          +SCL*CS*CT*D_DLB 
C     '        +D_CCL*SS*DMB +CCL*D_SS*DMB +CCL*SS*D_DMB 
C     '        +D_SSL*CC*CT*DMB +SSL*D_CC*CT*DMB +SSL*CC*D_CT*DMB 
C     '                                          +SSL*CC*CT*D_DMB 
C     '        -D_SSL*CS*ST*DTB -SSL*D_CS*ST*DTB -SSL*CS*D_ST*DTB
C     '                                          -SSL*CS*ST*D_DTB)/G1
C            DZI(3,ni)=(SCL*SS*ST*DLB+SSL*SC*ST*DMB+SSL*SS*CT*DTB)/G3
C            D_DZI(3,ni)=(D_SCL*SS*ST*DLB +SCL*D_SS*ST*DLB 
C     '                  +SCL*SS*D_ST*DLB +SCL*SS*ST*D_DLB 
C     '        +D_SSL*SC*ST*DMB +SSL*D_SC*ST*DMB +SSL*SC*D_ST*DMB
C     '                                          +SSL*SC*ST*D_DMB 
C     '        +D_SSL*SS*CT*DTB +SSL*D_SS*CT*DTB +SSL*SS*D_CT*DTB 
C     '                                          +SSL*SS*CT*D_DTB)/G3
C          ENDDO
C        ENDIF

C 25/2/97 LC removed section from :
C
C#### Subroutine: ENERGY
C###  Description:
C###    ENERGY calculates derivatives of strain energy function wrt 
C###    either principal strain invariants (KTYP55(nr)=1), or 
C###    principal extensions (KTYP55(nr)=2), or physical strains
C
C old MPN 13-Apr-96: init extns handled by 'growth' defm tensor
C                    see ZGTG53.
C          L0_fibre=CG(28)               !initial fibre ext ratio
C          L0_sheet=CG(29)               !initial sheet ext ratio
C          L0_sheetnormal=CG(30)         !initial sheetnormal ext ratio
C          E0_fibre=0.5d0*(L0_fibre*L0_fibre-1.d0) !initial fibre strain
C          E0_sheet=0.5d0*(L0_sheet*L0_sheet-1.d0) !initial cross strain
C          E0_sheetnormal=0.5d0*(L0_sheetnormal*L0_sheetnormal-1.d0) !initial sheet-normal strain
C          PP1=P1+E0_fibre
C          PP2=P2+E0_sheet
C          PP3=P3+E0_sheetnormal



C LC 25/2/97 archived section :
C       MPN 15-Sep-95: not used anymore
C       MPN 16-Aug-94: KTYP57(nr)=5 now used for entering cavity pressures 
C
C#### Subroutine: IPINI5
C###  Description:
C###    IPINI5 inputs initial conditions and boundary conditions.
C
C!!!! OLD!
C!!!! MPN 15-Sep-95: not used anymore
C!!!!
C!news MPN 16-Aug-94: KTYP57(nr)=5 now used for entering cavity pressures
C      IF(KTYP57(nr).EQ.4) THEN
C!old      IF(KTYP57(nr).EQ.4.OR.KTYP57(nr).EQ.5) THEN
C        IF(IOTYPE.NE.3) THEN
CC!!! should use NYNR here instead of NYT
C          NYT(2,1,nx)=NYT(2,1,nx)+1  !extra global cavity
C          YP(NYT(2,1,nx),1)=0.0d0     !pressure equation
C          YP(NYT(2,1,nx),2)=0.0d0
C          YP(NYT(2,1,nx),3)=0.0d0
C          FIX(NYT(2,1,nx),1)=.FALSE.
C          FIX(NYT(2,1,nx),2)=.FALSE.
C          FIX(NYT(2,1,nx),3)=.TRUE.
C        ENDIF
C        IF(KTYP5.EQ.1) THEN
C          IF(IOTYPE.NE.3) THEN
C            YP(NYT(2,1,nx),4)=0.0d0
C            PCAVITY=0.0d0
C          ENDIF
C        ELSE IF(KTYP5.EQ.2) THEN
C          FORMAT='($,'' Enter the initial ventricular pressure'//
C     '      ' [0.0]: '',G25.17)'
C          RDEFLT(1)=0.0d0
C          IF(IOTYPE.EQ.3) RDATA(1)=PCAVITY
C          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '      FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            PCAVITY=RDATA(1)
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              DO iface=1,2
C                IF(FIXP(iface,ne)) PF(iface,ne)=PCAVITY
C              ENDDO !iface
C            ENDDO !noelem (ne)
C            YP(NYT(2,1,nx),4)=PCAVITY
C          ENDIF
C        ELSE IF(KTYP5.EQ.3) THEN
C          IF(IOTYPE.NE.3) THEN
C            DO iface=1,2
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                IF(FIXP(iface,ne)) THEN
C                  PCAVITY=PF(iface,ne)
C                  YP(NYT(2,1,nx),4)=PCAVITY
C                  GOTO 9998
C                ENDIF
C              ENDDO !noelem (ne)
C            ENDDO !iface
C          ENDIF
C        ENDIF
C      ENDIF

C 25/2/97 LC removed from : MPN 15-Sep-95: not used anymore
C#### Subroutine: OPINI5
C###  Description:
C###    OPINI5 outputs initial and boundary data.
C
C!!!! MPN 15-Sep-95: not used anymore
C!!!!
C      IF(KTYP57(nr).GT.1) THEN
C        WRITE(OP_STRING,
C     '      '(/'' Element pressures (applied to Xi(3) faces)'','':''/)')
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        DO noelem=1,NEELEM(0)
C          ne=NEELEM(noelem)
C          WRITE(OP_STRING,
C     '      '('' Element'',I3,'' NW='',I1,''  Increments:       '''
C     '      //'2('' PE('',I1,'')='',D11.4))') 
C     '      ne,NW(ne,1),(i,PE(i,ne),i=1,2)
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C          WRITE(OP_STRING,'(16X,''  Current pressures:'''
C     '      //'2('' PF('',I1,'')='',D11.4))') 
C     '      (i,PF(i,ne),i=1,2)
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        ENDDO
C      ENDIF


C LC 25/2/97 removed section : 
C MPN 28-Jun-1995: integrals now wrt deformed coords
C#### Subroutine: PFRF
C###  Description:
C###    PFRF evaluates contribution RF(ns,nj) to element residuals RE 
C###    from the pressure PF(i) acting on the Xi(3)=IXF face 
C###    (IXF=iface-1).
C
C old MPN 28-Jun-1995: integrals now wrt deformed coords
C          IF(ITYP10(nr).EQ.1) THEN
C            DO nhx=1,NJ_LOC(NJL_GEOM,0)
C              nx=NH_LOC(nhx,nx)
C              DO ni=1,NITB
C                DZI(nhx,ni)=ZG(nhx,NU1(ni))
C              ENDDO !ni
C            ENDDO !nhx
C          ELSE IF(ITYP10(nr).EQ.2) THEN
C            RX=XG(1,1)
C            RZ=ZG(1,1)
C            DT=ZG(2,1)-XG(2,1)
C            CT=DCOS(DT)
C            ST=DSIN(DT)
C            DO ni=1,NITB
C              DZI(1,ni)= ZG(1,NU1(ni))*CT-RZ*ST*ZG(2,NU1(ni))
C              DZI(2,ni)=(ZG(1,NU1(ni))*ST+RZ*CT*ZG(2,NU1(ni)))/RX
C              DZI(3,ni)= ZG(3,NU1(ni))
C            ENDDO !ni
C          ELSE IF(ITYP10(nr).EQ.3) THEN
C            RX=XG(1,1)
C            RZ=ZG(1,1)
C            DT=ZG(2,1)-XG(2,1)
C            CT=DCOS(DT)
C            ST=DSIN(DT)
C            CX=DCOS(XG(3,1))
C            SX=DSIN(XG(3,1))
C            CZ=DCOS(ZG(3,1))
C            SZ=DSIN(ZG(3,1))
C            CC=CX*CZ
C            SS=SX*SZ
C            CS=CX*SZ
C            SC=SX*CZ
C            DO ni=1,NITB
C              RB=ZG(1,NU1(ni))
C              TB=ZG(2,NU1(ni))
C              PB=ZG(3,NU1(ni))
C              DZI(1,ni)=CC*(RB*CT-RZ*ST*TB)-CS*RZ*CT*PB+SC*RZ*PB+SS*RB
C              DZI(2,ni)=(RZ*CZ*CT*TB+(RB*CZ-RZ*SZ*PB)*ST)/(RX*CX)
C              DZI(3,ni)=(CC*RZ*PB+CS*RB+SC*(RZ*ST*TB-RB*CT)+SS*RZ*CT*PB)
C     '                                                   /RX
C            ENDDO !ni
C          ELSE IF(ITYP10(nr).EQ.4) THEN
C            SLX=DSINH(XG(1,1))
C            SLZ=DSINH(ZG(1,1))
C            SMX=DSIN (XG(2,1))
C            SMZ=DSIN (ZG(2,1))
C            CLX=DSQRT(1.0d0+SLX*SLX)
C            CLZ=DSQRT(1.0d0+SLZ*SLZ)
C            CMX=DSQRT(1.0d0-SMX*SMX)
C            CMZ=DSQRT(1.0d0-SMZ*SMZ)
C            CSLX=CLX/SLX
C            CSMX=CMX/SMX
C            DT=ZG(3,1)-XG(3,1)
C            CT=DCOS(DT)
C            ST=DSIN(DT)
C            CCL=CLX*CLZ
C            CSL=CLX*SLZ
C            SCL=SLX*CLZ
C            SSL=SLX*SLZ
C            CC=CMX*CMZ
C            CS=CMX*SMZ
C            SC=SMX*CMZ
C            SS=SMX*SMZ
C            G1=SLX*SLX+SMX*SMX
C            G3=SLX*SLX*SMX*SMX
C            DO ni=1,NITB
C              DLB=ZG(1,NU1(ni))
C              DMB=ZG(2,NU1(ni))
C              DTB=ZG(3,NU1(ni))
C              DZI(1,ni)=(( SSL*CC+CCL*SS*CT)*DLB+(-SCL*CS+CSL*SC*CT)*DMB
C     '                                                -CSL*SS*ST*DTB)/G1
C              DZI(2,ni)=((-CSL*SC+SCL*CS*CT)*DLB+( CCL*SS+SSL*CC*CT)*DMB
C     '                                                -SSL*CS*ST*DTB)/G1
C              DZI(3,ni)=(SCL*SS*ST*DLB+SSL*SC*ST*DMB+SSL*SS*CT*DTB)/G3
C            ENDDO !ni
C          ENDIF


C LC 25/2/97 removed section from : Temporary MPN 12-Nov-94
C
C#### Subroutine: OPSAIL
C###  Description:
C###    OPSAIL performs sail parameter output.
C 
C      nv=1 ! Temporary MPN 12-Nov-94
C      DO nj=1,3
C        ZTACK(nj) =XP(1,nv,nj,NPTACK)
C        ZCLEW(nj) =XP(1,nv,nj,NPCLEW)
C        ZHEADL(nj)=XP(1,nv,nj,NPHEADL)
C        ZHEADT(nj)=XP(1,nv,nj,NPHEADT)
C      ENDDO
C      CLUFF=SQRT((ZHEADL(1)-ZTACK(1))**2+(ZHEADL(3)-ZTACK(3))**2)
C      CLECH=SQRT((ZHEADT(1)-ZCLEW(1))**2+(ZHEADT(3)-ZCLEW(3))**2)
C      DO j=1,2
C        XYLUF1(1,j,1)=0.D0
C        XZLUF1(1,j,1)=0.D0
C        XYLUF2(1,j,2)=0.D0
C        XZLUF2(1,j,2)=0.D0
C        XYLEC1(1,j,1)=0.D0
C        XZLEC1(1,j,1)=0.D0
C        XYLEC2(1,j,2)=0.D0
C        XZLEC2(1,j,2)=0.D0
C      ENDDO
C      XYLUF2(1,1,2)=CLUFF
C      XZLUF2(1,1,2)=CLUFF
C      XYLEC2(1,1,2)=CLECH
C      XZLEC2(1,1,2)=CLECH
C      DXLUFF=CLUFF/DBLE(NTSECT-1)
C      DXLECH=CLECH/DBLE(NTSECT-1)
C      np1=NPTACK
C      DO nosect=1,NTSECT
C        np2=np1+1
C        np3=np1+2
C        COS=(ZHEADL(1)-ZTACK(1))/CLUFF
C        SIN=(ZHEADL(3)-ZTACK(3))/CLUFF
C        XLUFF=DXLUFF*DBLE(nosect-1)
C        XLECH=DXLECH*DBLE(nosect-1)
C        COS=(ZHEADL(1)-ZTACK(1))/CLUFF
C        SIN=(ZHEADL(3)-ZTACK(3))/CLUFF
C        DRAFT=(XP(1,nv,1,np2)-XP(1,nv,1,np1))/
C     '    (XP(1,nv,1,np3)-XP(1,nv,1,np1))
C        XYCHO1(1,1,2)= (XP(1,nv,1,np2)-XP(1,nv,1,np1))*COS
C     '                +(XP(1,nv,2,np2)-XP(1,nv,2,np1))*SIN
C        XYCHO1(1,2,2)=-(XP(1,nv,1,np2)-XP(1,nv,1,np1))*SIN
C     '                +(XP(1,nv,2,np2)-XP(1,nv,2,np1))*COS
C        XYCHO1(2,1,1)= XP(1,nv,1,np1)*COS+XP(1,nv,2,np1)*SIN
C        XYCHO1(2,2,1)=-XP(1,nv,1,np1)*SIN+XP(1,nv,2,np1)*COS
C        np1=np1+3
C      ENDDO

C      CALL EXITS('OPSAIL')
C 9998  RETURN


C MPN 17Jul2000 removed section from :
C#### Subroutine: ZGTG5A
C###  Description:
C###    ZGTG5A adds active fibre stress to TG array.

C**** OLD MPN 5Jun2000
C**** FEXT(1,ng,ne) is current muscle fibre extension ratio
C****   "  2    "   "  previous   "     "       "       "
C****   "  3    "   "  muscle fibre ext. ratio at time of activation
C****   "  4    "   "  zero when gauss point is inactive else 1
C****   "  5    "   " previous hereditary integral for 1st time constant
C****   "  6    "   "     "         "         "     "  2nd  "      "
C****   "  7    "   "     "         "         "     "  3rd  "      "
C**** KTYP59(nr) is elastance/Hill-type/fading-memory formulation
C**** DEL_T  is the time step used (taken from load step loop in FE07)
C**** TV_SLO is slope of tension/vel. in lengthening (before yield)
C**** YIELDR is ratio of yield tension to isometric tension
C**** SNLPA  is static nonlinearity parameter "a"
C**** NTACTV is #dynamic terms in the material response function
C**** ACOEFF(nactv), nactv=1,NTACTV are coeffs   for lin dynamic terms
C**** ALFA(nactv),     "       "     "time constants  "     "      "
C****
C**** KTYP59(nr)=3 Fading memory model
C**** ========
C**** PARAMS FOR THE EXTRACELLULAR CA CONC RELATION Cao(t)****
C**** time1  is the time (s)    at junction of 1st & 2nd cubic elements.
C**** Cao1   is the Ca conc. (mM)     "     "   "  "  "    "      "
C**** Caslo1 is the slope of the fn   "     "   "  "  "    "      " 
C**** time2  is the time (s)    at junction of 2nd & 3rd cubic elements.
C**** Cao2   is the Ca conc. (mM)     "     "   "  "  "   "  "    "
C**** Caslo2 is the slope of the fn   "     "   "  "  "   "  "    "
C**** t_end  is the time (s) after onset of contraction, 
C****                        when Cao drops to zero
C****
C**** PARAMETERS FOR THE ISOMETRIC TENSION RELATION, T0(Cao,FEXT) ****
C**** FEXTo  is the exten ratio at which there is zero active tension
C**** FEXTmx is the max permitted exten ratio (limited by expt'l data)
C**** T2_k   is the constant of the calcium dep nodal value @ XI2=1
C**** T2_a   is the pole     "   "     "     "    "     "   @ XI2=1
C**** dT1_k  is the constant of the calcium dep nodal slope @ XI2=0
C**** dT1_a  is the pole     "   "     "     "    "     "   @ XI2=0
C**** dT2_k  is the constant "   "     "     "    "     "   @ XI2=1
C**** dT2_a  is the pole     "   "     "     "    "     "   @ XI2=1
C****
C**** VARIABLES ****
C**** XI1    is the local [0,1] variable varying with time
C**** XI2    is the local [0,1] variable varying with FEXT
C**** Ca(i), i=1,4 are the time dep nodal params for the Cao function
C**** T(i),  i=1,4 are the Ca and length dep nodal params for the T0 fn

      REAL*8 Ca(4),Cai,Cao,Caosqr,DFEXT,HMT_C50,
     '  HMT_n,HMT_pC50,HMT_ZSS,lambda,PSI10,PSI11,PSI20,PSI21,
     '  SUM,T(4),T0,TIMEAA,Q,VEL,XI,XI1,XI2
      REAL*8 Cao1,Cao2,Caslo1,Caslo2,dT1_a,dT1_k,dT2_k,dT2_a,
     '  FEXTo,FEXTmx,T2_a,T2_k,t_end,time1,time2
      DATA Cao1 /1.347d0/, Cao2 /0.6918d0/, Caslo1 /22.48d0/, 
     '  Caslo2 /-1.916d0/,
     '  dT1_a /0.77d0/, dT1_k /449.d0/, dT2_k /31.d0/, dT2_a /0.13d0/,
     '  FEXTo /0.865d0/, FEXTmx /1.189d0/, T2_a /0.13d0/, T2_k /137.d0/,
     '  t_end /1.023d0/, time1 /0.1d0/, time2 /0.4173d0/

      PSI10(XI) = 1.0d0 - 3.0d0*XI*XI + 2.0d0*XI*XI*XI !Cubic H basis fn
      PSI11(XI) = XI*(XI-1.0d0)*(XI-1.0d0)
      PSI20(XI) = XI*XI*(3.0d0-2.0d0*XI)
      PSI21(XI) = XI*XI*(XI-1.0d0)



C news MPN 23May2000: adding Steady State HMT (several changes below)
      IF(KTYP59(nr).EQ.1.OR.      !SS tension-length-Ca relation
     '   KTYP59(nr).EQ.2) THEN    !Steady state HMT

        lambda = FEXT(1)
        Cai = FEXT(4)
        IF(Cai.LT.0.0d0) THEN
          WRITE(OP_STRING,'('' Cannot handle negative Cai='','
     '      //'D12.3,''; continuing using zero active stress'')') Cai
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GOTO 9998
        ENDIF
        IF(KTYP59(nr).EQ.1) THEN         !SS tension-length-Ca relation
          TNA = Ca_max*Cai**Ca_h/(Cai**Ca_h + Ca_c50**Ca_h)
     '          *Tref*(1.d0+T0_beta*(lambda-1.d0))
        ELSE IF(KTYP59(nr).EQ.2) THEN    !Steady State HMT
C         Taken from Hunter, McCulloch + ter Keurs, 
C         Prog Biophys Mol Biol, 69 (1998) pp. 300-301
          HMT_pC50 = HMT_pC50_ref*(1.d0+HMT_pC50_beta*(lambda-1.d0))
          HMT_C50 = 1.d1**(3.d0-HMT_pC50) ! changed for milliMolar !!!
          HMT_n = HMT_n_ref*(1.d0+HMT_n_beta*(lambda-1.d0))
          HMT_ZSS = (Cai**HMT_n)/(Cai**HMT_n + HMT_C50**HMT_n)
          TNA = Tref*(1.d0+T0_beta*(lambda-1.d0))*HMT_ZSS
        ENDIF
        IF(DOP) THEN
C$        call mp_setlock()
          WRITE(OP_STRING,'('' Cai='',D12.3,'' TNA='',D12.3)') Cai,TNA
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$        call mp_unsetlock()
        ENDIF
C newe MPN 23May2000

      ELSE IF(KTYP59(nr).EQ.3) THEN !fading memory model
CC****   Determine the time after activation.
C        TIMEAA=DEL_T*DBLE(NOSTEP)-TDELAY
C        IF(DOP) THEN
CC$        call mp_setlock()
C          WRITE(OP_STRING,'('' TIMEAA='',D13.6)') TIMEAA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C        ENDIF
        
C        IF(TIMEAA.GE.0.0d0)THEN
C          IF(FEXT(4).EQ.0.0d0) THEN
C            FEXT(3)=FEXT(1)
C            FEXT(4)=1.0d0
C          ENDIF
CC****     Determine change in sarcomere length and current estimate of 
CC****     the velocity. Note: positive velocity for lengthening.
C          DFEXT = (FEXT(1)-FEXT(2))
C          VEL = DFEXT/DEL_T
C          IF(DOP) THEN
CC$          call mp_setlock()
C            WRITE(OP_STRING,'('' FEXT(1..7)='',7D10.3/'' DFEXT='','
C     '        //'D11.4,'' DFEXT='',D11.4)') (FEXT(i),I=1,7),DFEXT,VEL
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
C          ENDIF
        
CC****     Calculate Cao as a function of time after activation.
CC****     Uses three cubic hermite elements.
C          XI1 = 0.0d0
C          DO j=1,4
C            Ca(j) = 0.0d0
C          ENDDO
C          IF((TIMEAA.GE.0.0d0).AND.(TIMEAA.LT.time1))THEN
C            XI1 = TIMEAA/time1
C            Ca(2) = Cao1
C            Ca(4) = Caslo1*time1
C          ELSE IF((TIMEAA.GE.time1).AND.(TIMEAA.LT.time2))THEN
C            XI1 = (TIMEAA-time1)/(time2-time1)
C            Ca(1) = Cao1
C            Ca(2) = Cao2
C            Ca(3) = Caslo1*(time2-time1)
C            Ca(4) = Caslo2*(time2-time1)
C          ELSE IF((TIMEAA.GE.time2).AND.(TIMEAA.LT.t_end))THEN
C            XI1 = (TIMEAA-time2)/(t_end-time2)
C            Ca(1) = Cao2
C            Ca(3) = Caslo2*(t_end-time2)
C          ENDIF
C          Cao = PSI10(XI1)*Ca(1) + PSI20(XI1)*Ca(2) + PSI11(XI1)*Ca(3) +
C     '          PSI21(XI1)*Ca(4)
C          IF(DOP) THEN
CC$          call mp_setlock()
C            WRITE(OP_STRING,'('' Cao='',D13.6)') Cao
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
C          ENDIF
        
CC****     T0 calculated as a fn of FEXT(1) and Cao. 
CC****     Uses 1 cubic hermite element.
C          XI2 = 0.0d0
C          DO j=1,4
C            T(j) = 0.0d0
C          ENDDO
C          IF(FEXT(1).GE.FEXTo)THEN
C            CALL ASSERT(FEXT(1).LE.FEXTmx,'>>Extension Ratio greater'
C     '        //' than maximum FEXT ',ERROR,*9999)
C            XI2 = (FEXT(1)-FEXTo)/(FEXTmx-FEXTo)  !Local Variable XI2.
C            Caosqr=Cao*Cao
C            T(1) = 0.0d0                          !Ca dep parameters
C            T(2) = T2_k  * Caosqr/(Caosqr+T2_a )  !(Hill type relations).
C            T(3) = dT1_k * Caosqr/(Caosqr+dT1_a)
C            T(4) = dT2_k * Caosqr/(Caosqr+dT2_a)
C          ENDIF
C          T0 = PSI10(XI2)*T(1) + PSI20(XI2)*T(2) +PSI11(XI2)*T(3) 
C     '       +PSI21(XI2)*T(4)
C          IF(DOP) THEN
CC$          call mp_setlock()
C            WRITE(OP_STRING,'('' T0='',D13.6)') T0
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
C          ENDIF
        
CC****     Determine the magnitude of the active muscle fibre stress.
C          IF(VEL.GE.0.0d0) THEN        !lengthing
C            TNA = TV_SLO*VEL + T0
C            IF(TNA.GT.(YIELDR*T0)) TNA = YIELDR*T0
        
C          ELSE IF(VEL.LT.0.0d0) THEN   !shortening
CC****       Calculate hereditary integral, summing over each rate constant
C            Q=0.0d0
C            DO nactv=1,NTACTV
C              CALL ASSERT(ALFA(nactv)*DEL_T.LT.80.0d0, !to avoid underflow
C     '          '>>Reduce time step using the DEFINE ACTIVE command',
C     '          ERROR,*9999)
C              Q = Q + ACOEFF(nactv)*(DEXP(-ALFA(nactv)*DEL_T)
C     '          *FEXT(nactv+4)+DEXP(-ALFA(nactv)*DEL_T/2.0d0)*DFEXT)
C            ENDDO
C            IF(DOP) THEN
CC$            call mp_setlock()
C              WRITE(OP_STRING,'('' Q='',D13.6)') Q
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
C            ENDIF
C            TNA=T0*(1.0d0+SNLPA*Q)/(1.0d0-Q)
C          ENDIF
C        ENDIF !timeaa

C**** END OLD MPN 5Jun2000


      REAL*8 FUNCTION AGP(nr,ns,PPG,XG,ZG)

C**** Evaluates integrand of the domain integral in the Laplace eqtn
C**** for perfusion.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER nr,ns
      REAL*8 PPG(64,*),XG(NJM,*),ZG(NHM,*)
!     Local Variables
      INTEGER nh,nj,njj,NU1(0:3)
      REAL*8 AA,G1,G3,RAD,SLX,SMX
      DATA NU1/1,2,4,7/

      nh=NJ_LOC(NJL_GEOM,0)+2
      IF(ITYP10(nr).EQ.1) THEN
        AGP=0.0D0
        DO njj=1,NJ_LOC(NJL_GEOM,0)
          nj=NJ_LOC(NJL_GEOM,njj)
          AGP=AGP+ZG(nh,NU1(nj))*PPG(ns,1+nj)
        ENDDO
      ELSE IF(ITYP10(nr).EQ.2) THEN
        AGP=ZG(nh,2)*PPG(ns,2)
        IF(NJ_LOC(NJL_GEOM,0).GE.2) 
     '    AGP=AGP+ZG(nh,4)*PPG(ns,3)*XG(1,1)**2
        IF(NJ_LOC(NJL_GEOM,0).EQ.3) AGP=AGP+ZG(nh,7)*PPG(ns,4)
      ELSE IF(ITYP10(nr).EQ.3) THEN
        RAD=XG(1,1)
        AGP=ZG(nh,2)*PPG(ns,2)
        IF(NJ_LOC(NJL_GEOM,0).GE.2) 
     '    AGP=AGP+ZG(nh,4)*PPG(ns,3)*(RAD*DCOS(XG(3,1)))**2
        IF(NJ_LOC(NJL_GEOM,0).EQ.3) AGP=AGP+ZG(nh,7)*PPG(ns,4)* RAD**2
      ELSE IF(ITYP10(nr).EQ.4) THEN
        AA=FOCUS**2
        SLX=DSINH(XG(1,1))
        SMX=DSIN (XG(2,1))
        G1 =AA*(SLX*SLX+SMX*SMX)
        G3 =AA*(SLX*SLX*SMX*SMX)
        AGP=G1*(ZG(nh,2)*PPG(ns,2)+ZG(nh,4)*PPG(ns,3))
        IF(NJ_LOC(NJL_GEOM,0).EQ.3) AGP=AGP+G3*ZG(nh,7)*PPG(ns,4)
      ENDIF

      RETURN
      END


      SUBROUTINE IPPRES(NBH,NEELEM,NELIST,NKH,NPNODE,nr,NVHP,NW,
     '  PE,PF,XP,ZP,FIXP,ERROR,*)

C#### Subroutine: IPPRES
C###  Description:
C###    Inputs pressure boundary conditions
C#### FIXP(1..2,ne) specifies pressure bcs applied to Xi3=0,1 faces
C#### PE(1..2,ne) is the pressure increment applied to the Xi3=0,1 face
C#### PF(1..2,ne) is the current pressure on the Xi3=0,1 face
C**** If KTYP57(nr)>1 and KTYP51(nr) = 3 or 4 (3D or membrane)
C****   a pressure b.c. is applied to the Xi(3)=0,1 faces of
C****   elements with NW(ne,1)=2,3, respectively, or both if NW(ne,1)=4.
C**** If KTYP57(nr)>1 and KTYP51(nr) = 5 (string)
C****   a pressure b.c. is applied to the string normal.

C!!!! OLD!
C!!!! MPN 15-Sep-95: pressure bc's now element variables, so
C!!!! this routine is no longer used or called from IPINI5.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:acti01.cmn'
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grou00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:lvpr00.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NE_R_M),
     '  NKH(NHM,NPM,NCM),NPNODE(0:NP_R_M,0:NRM),nr,NVHP(NHM,NPM,NCM),
     '  NW(NEM,2)
      REAL*8 PE(2,NEM),PF(2,*),XP(NKM,NVM,NJM,NPM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIXP(2,*)
!     Local Variables
      INTEGER MAX_LIST
      PARAMETER(MAX_LIST=30)
      INTEGER i,I1FACE,I2FACE,IBEG,ICHAR,IEND,iface,INFO,
     '  IOSTAT,IPF07,IXFACE,n,N1char,N2char,
     '  N1ELEM,N1GREL,N2ELEM,NBPRES,ne,NE1,
     '  nh,nhh2,NIFACE,nk,noelem,nogrel,nolist,nonode,NOQUES,
     '  np,NPF07,NPRBC,nv
      REAL*8 CAV_PRESS_CUR,CAV_PRESS_INCR,PRESS(100)
      CHARACTER CHAR*100,CHAR1*100,CUPPER*100,LABEL*100
      LOGICAL DEFINED,FILEIP

      CALL ENTERS('IPPRES',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=0 !temporary (needs deleting later)

      FORMAT='('' Specify option for Xi(3) faces [1]:'''//
     '  '/''   (1) No pressure boundary conditions'''//
     '  '/''   (2) Boundary pressure increments entered'''//
     '  '/''   (3) Boundary pressures read from file'''//
     '  '/''   (4) Bdry pressures computed from prescribed volume'''//
     '  '/''   (5) Bdry pressures computed from cavity pressure(s)'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP57(nr)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,5,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        KTYP57(nr)=IDATA(1)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO iface=1,2
            PF(iface,ne)=0.0D0
          ENDDO
        ENDDO
      ENDIF

!news MPN 16-Aug-94: KTYP57(nr)=5 now used for entering cavity pressures
      IF(KTYP57(nr).EQ.2.OR.KTYP57(nr).EQ.4) THEN !Press increm
!old  IF(KTYP57(nr).EQ.2.OR.KTYP57(nr).EQ.4.OR.KTYP57.EQ.5) THEN 
!old    !Press incr
        IF(IOTYPE.EQ.3) N1ELEM=1

 5820   CDATA(1)='ELEMENTS' !for use with group input
        FORMAT='($,'' Enter element #/name [EXIT]: '',I5)'
        IF(IOTYPE.EQ.3) THEN
          DO noelem=N1ELEM,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            N2ELEM=noelem
            IF(FIXP(1,ne).OR.FIXP(2,ne)) GOTO 5825
          ENDDO
          N2ELEM=0
 5825     IF(N2ELEM.EQ.0) THEN
            IDATA(1)=0
          ELSE
            IDATA(1)=ne
            N1ELEM=N2ELEM+1
          ENDIF
          IDATA(0)=1 !write out one element at a time
        ENDIF !iotype=3
        N1char=3  !start of format string after $
        N2char=32 !end of string (excluding data format) from $
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,N1char,N2char,NOQUES,
     '    FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NEM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IDATA(1).NE.0) THEN !not default
          NELIST(0)=IDATA(0)
          DO n=1,IDATA(0)
            NELIST(n)=IDATA(n)
          ENDDO !n
          IF(KTYP51(nr).EQ.3) THEN !3D
            NIFACE=2
          ELSE IF(KTYP51(nr).GE.4) THEN !membrane,string or shell
            NIFACE=1
          ENDIF
          IF(IOTYPE.EQ.3) I1FACE=1

C         Enter extra information for first element in group
          ne=NELIST(1) !rest of group filled at end
          IF(KTYP51(nr).EQ.3.OR.KTYP51(nr).EQ.4) THEN !3D or membrane
 5830       FORMAT='($,'' Enter Xi(3) face number [EXIT]: '',I3)'
            IF(IOTYPE.EQ.3) THEN
              DO iface=I1FACE,NIFACE
                I2FACE=iface
                IF(FIXP(iface,ne)) GOTO 5835
              ENDDO !iface
              I2FACE=0
 5835         IF(I2FACE.EQ.0) THEN
                IDATA(1)=0
              ELSE
                IDATA(1)=I2FACE
                I1FACE=I2FACE+1
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NIFACE,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IDATA(1).NE.0) THEN
              iface=IDATA(1)
              IXFACE=iface-1
              WRITE(CHAR1,'(I1)') IXFACE
              FORMAT='($,'' Enter the current pressure   applied to'//
     '          ' the Xi(3)='//CHAR1(1:1)//' face [0.0]: '',G25.17)'
              RDEFLT(1)=0.0D0
              IF(IOTYPE.EQ.3) RDATA(1)=PF(iface,ne)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                NW(ne,1)=NW(ne,1)+iface
                PF(iface,ne)=RDATA(1)
                FIXP(iface,ne)=.TRUE.
              ENDIF
              FORMAT='($,'' Enter the pressure increment applied to'//
     '          ' the Xi(3)='//CHAR1(1:1)//' face [0.0]: '',G25.17)'
              RDEFLT(1)=0.0D0
              IF(IOTYPE.EQ.3) RDATA(1)=PE(iface,ne)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                PE(iface,ne)=RDATA(1)
              ENDIF
              GO TO 5830
            ENDIF

          ELSE !string
            FIXP(1,ne)=.TRUE.
            NW(ne,1)=2
            FORMAT='($,'' Enter the pressure increment'
     '        //' [0.0]: '',G25.17)'
            RDEFLT(1)=0.0D0
            IF(IOTYPE.EQ.3) RDATA(1)=PE(1,ne)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              PE(1,ne)=RDATA(1)
            ENDIF
          ENDIF !KTYP51(nr)

!news 	  MPN 21/12/93
      	  NBPRES=NBH(NJ_LOC(NJL_GEOM,0)+1,1,ne) !basis fn of hyd. press.
          IF(NNT(NBPRES).EQ.0) THEN  !element based pressure variables
            IF(NW(ne,1).EQ.1) THEN                       !no press bc's
	      NPRBC=0
      	    ELSE IF(NW(ne,1).EQ.2.OR.NW(ne,1).EQ.3) THEN  !1 press bc
    	      NPRBC=1
            ELSE IF(NW(ne,1).EQ.4) THEN                   !2 press bc's
              NPRBC=2
            ENDIF
            IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.2) THEN  !3D,Incompress
              IF(NPRBC.GT.NAT(NBPRES)) THEN
                WRITE(OP_STRING,'('' >>WARNING: There are more '
     '            //'press bc.s prescribed than there are aux params '
     '            //'for this element'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999) 
              ENDIF
            ELSE IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).EQ.3) THEN 
              !3D,Incomp+fluid
              IF(NPRBC.GT.NAT(NBPRES)-1) THEN
                WRITE(OP_STRING,'('' >>WARNING: There are more '
     '            //'than (# aux params-1) press bc.s prescribed '
     '            //'for this element'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999) 
              ENDIF
            ENDIF
          ENDIF
!newe	

C         Apply same information to rest of element group
          DO n=2,NELIST(0)
            ne=NELIST(n)
            IF(KTYP51(nr).EQ.3) THEN !3D
              NIFACE=2
            ELSE IF(KTYP51(nr).GE.4) THEN !membrane,string or shell
              NIFACE=1
            ENDIF
            NW(ne,1)=NW(NELIST(1),1)
            DO iface=1,NIFACE
              PF(iface,ne)=PF(iface,NELIST(1))
              PE(iface,ne)=PE(iface,NELIST(1))
              FIXP(iface,ne)=FIXP(iface,NELIST(1))
            ENDDO !iface
          ENDDO !n

          GO TO 5820
        ENDIF !input of element#s

!news MPN 16-Aug-94: KTYP57(nr)=5 now used for entering cavity pressures
        IF(KTYP57(nr).EQ.4) THEN
!old        IF(KTYP57(nr).EQ.4.OR.KTYP57(nr).EQ.5) THEN
C**       This should be generalized to account for two ventricles
C**       by looping over the next section and adding an index for
C**       and right...
C         CSIDE(1)='left'
C         CSIDE(2)='right'
C         DO iside=1,2
          FORMAT='($,'' Enter the region number of the'//
     '      ' ventricular cavity [2]: '',I1)'
          IDEFLT(1)=2
          IF(IOTYPE.EQ.3) IDATA(1)=NRCAVITY
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NRCAVITY=IDATA(1)
            DO nonode=1,NPNODE(0,NRCAVITY)
              np=NPNODE(nonode,NRCAVITY)
              DO nhh2=1,NJ_LOC(NJL_GEOM,0)
                nh=NJ_LOC(NJL_GEOM,nhh2)
                DO nv=1,NVHP(nh,np,1)
                  DO nk=1,NKH(nh,np,1)
                    ZP(nk,nv,nh,np,1)=XP(nk,nv,nh,np)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF(KTYP57(nr).EQ.4) THEN
            FORMAT='($,'' Enter the volume increment'//
     '        ' applied to this region [0.0]: '',G25.17)'
            RDEFLT(1)=0.0D0
            IF(IOTYPE.EQ.3) RDATA(1)=VINCR
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) VINCR=RDATA(1)
C**         May want to add this bit later, if reference volume is
C**         different from the current volume.
C           FORMAT='($,'' Specify option for reference volume [1]:'''//
C    '        '/''   (1) Undeformed cavity volume'''//
C    '        '/''   (2) Current volume'''//
C    '        '/''   (3) Prescribed volume'''//
C    '        '/$,''    '',I1)'
C           IF(IOTYPE.EQ.3) IDATA(1)=IVREF
C           CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C    '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
C    '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C           IF(IOTYPE.NE.3) IVREF=IDATA(1)
            FORMAT='($,'' Enter the fluid compliance (penalty factor)'//
     '        ' [1.E-5]: '',G25.17)'
            RDEFLT(1)=1.0D-5
            IF(IOTYPE.EQ.3) RDATA(1)=CARTERY
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              CARTERY=RDATA(1)
              RARTERY=0.0D0 !Not used for the prescribed volume case
            ENDIF
!news MPN 16-Aug-94: KTYP57(nr)=5 now used for entering cavity pressures
!old
!          ELSE IF(KTYP57(nr).EQ.5) THEN
!            FORMAT='($,'' Enter the time of outflow valve opening'//
!     '        ' [0.0]: '',G25.17)'
!            RDEFLT(1)=0.0D0
!            IF(IOTYPE.EQ.3) RDATA(1)=TOPEN
!            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
!     '        FORMAT,1,
!     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
!     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
!            IF(IOTYPE.NE.3) TOPEN=RDATA(1)
!            FORMAT='($,'' Enter the arterial compliance'//
!     '        ' [0.0]: '',G25.17)'
!            RDEFLT(1)=0.0D0      ! **** Needs default value here
!            IF(IOTYPE.EQ.3) RDATA(1)=CARTERY
!            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
!     '        FORMAT,1,
!     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
!     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
!            IF(IOTYPE.NE.3) THEN
!              CARTERY=RDATA(1)
!            ENDIF
!            FORMAT='($,'' Enter the total peripheral resistance'//
!     '        ' [0.0]: '',G25.17)'
!            RDEFLT(1)=0.0D0      ! **** Needs default value here
!            IF(IOTYPE.EQ.3) RDATA(1)=RARTERY
!            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
!     '        FORMAT,1,
!     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
!     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
!            IF(IOTYPE.NE.3) THEN
!              RARTERY=RDATA(1)
!            ENDIF
          ENDIF
        ENDIF

      ELSE IF(KTYP57(nr).EQ.3) THEN !Boundary pressures input from file
        IF(IOTYPE.EQ.3) THEN
          CALL TRIM(FILE07,IBEG,IEND)
          CDATA(1)=FILE07(IBEG:IEND)
        ENDIF
        CALL TRIM(FILE00,IBEG,IEND)
        FORMAT='($,'' Enter file name ['//FILE00(IBEG:IEND)//']'
     '    //' (file extension is .PRESSURE): '',A)'
        CDEFLT(1)='FILE'
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          CALL TRIM(CDATA(1),IBEG,IEND)
          FILE07=CDATA(1)(IBEG:IEND)
          CALL TRIM(FILE07,IBEG,IEND)
          CALL OPENF(7,'DISK',FILE07(IBEG:IEND)//'.pressure','OLD',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          IF(KTYP51(nr).EQ.3) THEN
            NIFACE=2
          ELSE IF(KTYP51(nr).EQ.4) THEN
            NIFACE=1
          ENDIF
          READ(7,*) ((FIXP(i,NEELEM(noelem,nr)),i=1,NIFACE),
     '      noelem=1,NEELEM(0,nr))
          NPF07=0  !# of prescribed pressures in file
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO iface=1,NIFACE
              IF(FIXP(iface,ne)) NPF07=NPF07+1
            ENDDO
          ENDDO
          READ(7,*,IOSTAT=IOSTAT) TIME,(PRESS(IPF07),IPF07=1,NPF07)
          IF(IOSTAT.NE.0) THEN
            IF(IOSTAT.EQ.36) THEN
              ERROR=' End of pressure file'
            ELSE
              ERROR=' Real data read error in pressure file'
            ENDIF
            CALL CLOSEF(7,ERROR,*9999)
            GO TO 9999
          ENDIF
          IPF07=0
          DO noelem=1,NEELEM(0,nr)
            IPF07=IPF07+1
            ne=NEELEM(noelem,nr)
            DO iface=1,NIFACE
              IF(FIXP(iface,ne)) PF(iface,ne)=PRESS(IPF07)
            ENDDO
          ENDDO
          CALL CLOSEF(7,ERROR,*9999)
        ENDIF

!news MPN 16-Aug-94: KTYP57(nr)=5 now used for entering cavity pressures
      ELSE IF(KTYP57(nr).EQ.5) THEN !Bdry press computed from cavity press
C ***   Only used for entering from file or prompt
        IF(IOTYPE.NE.3) THEN
 6000     FORMAT='($,'' Enter cavity number [EXIT]: '',I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

          IF(IDATA(1).NE.0) THEN
 6010       FORMAT='($,'' Enter the element list (sep by ;) '
     '        //'or group name: '',A)'
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C ***       Replace all ;'s with commas in CDATA(1)
            CALL INSERT(CDATA(1),';',',')
            IF(DOP) THEN
              CALL TRIM(CDATA(1),IBEG,IEND)
              WRITE(OP_STRING,'('' Character elem list: '',A)')
     '          CDATA(1)(IBEG:IEND)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999) 
            ENDIF
C ***       Parse list of elements
            CALL PARSIL(CDATA(1),MAX_LIST,NELIST(0),NELIST(1),
     '        ERROR,*6060)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' # elems in list: '',I3)') NELIST(0)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999) 
              WRITE(OP_STRING,'('' Elem list: '',20I4)') 
     '          (NELIST(nolist),nolist=1,NELIST(0))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999) 
            ENDIF
C ***       Check that elements in list are defined in current region
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              DEFINED=.FALSE.
              DO noelem=1,NEELEM(0,nr)
                NE1=NEELEM(noelem,nr)
                IF(ne.EQ.NE1) THEN
                  DEFINED=.TRUE.
                  GO TO 6055
                ENDIF
              ENDDO
 6055         IF(.NOT.DEFINED) THEN
                WRITE(CHAR,'(I5)') ne
                WRITE(OP_STRING,'('' >>Element '//CHAR(1:4)
     '            //' is not defined'')') 
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999) 
                GO TO 6010  !reenter element list
              ENDIF
            ENDDO

            GO TO 6200 !elem list determined
 6060       CONTINUE
C ***       Parse group name
            CALL TRIM(CDATA(1),IBEG,IEND)
            CHAR=CUPPER(CDATA(1)(IBEG:IEND))
            N1GREL=0
            DO nogrel=1,NTGREL
              LABEL=CUPPER(LAGREL(nogrel))
              IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG:IEND)) THEN
                N1GREL=nogrel
                GO TO 6100
              ENDIF
            ENDDO
 6100       IF(N1GREL.GT.0) THEN
              NELIST(0)=LIGREL(0,N1GREL)
              DO nolist=1,NELIST(0)
                NELIST(nolist)=LIGREL(nolist,N1GREL)
              ENDDO
            ELSE
              WRITE(OP_STRING,'('' >>Group name '//CHAR(IBEG:IEND)
     '          //' is not defined'')') 
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999) 
              GO TO 6010 !reenter group name
            ENDIF
 6200       CONTINUE !Group name and element list determined

 6300       FORMAT='($,'' Enter Xi(3) face number [EXIT]: '',I3)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NIFACE,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IDATA(1).NE.0) THEN
              iface=IDATA(1)
              IXFACE=iface-1
              WRITE(CHAR1,'(I1)') IXFACE
              FORMAT='($,'' Enter the current pressure   applied to'//
     '          ' the Xi(3)='//CHAR1(1:1)//' face [0.0]: '',G25.17)'
              RDEFLT(1)=0.0D0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              CAV_PRESS_CUR=RDATA(1)

              FORMAT='($,'' Enter the pressure increment applied to'//
     '          ' the Xi(3)='//CHAR1(1:1)//' face [0.0]: '',G25.17)'
              RDEFLT(1)=0.0D0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              CAV_PRESS_INCR=RDATA(1)

              DO nolist=1,NELIST(0)
                ne=NELIST(nolist)
                PF(iface,ne)=CAV_PRESS_CUR
                PE(iface,ne)=CAV_PRESS_INCR
                FIXP(iface,ne)=.TRUE.
                IF(NW(ne,1).EQ.1) THEN
                  NW(ne,1)=NW(ne,1)+iface
                ELSE IF(NW(ne,1).EQ.2) THEN
                  NW(ne,1)=NW(ne,1)+2*(iface-1)
                ELSE IF(NW(ne,1).EQ.3) THEN
                  NW(ne,1)=NW(ne,1)+2-iface
                ENDIF
              ENDDO
              GO TO 6300 !next Xi3 face
            ENDIF
            GO TO 6000 !next cavity
          ENDIF

C ***     This is only used for entering cavity pressures and then
C ***     setting face pressures for each of the elements in given list.
C ***     For future IPINIT file writes and then reads of the same data 
C ***     use KTYP57(nr)=2 code as cavity element list is not stored.
          KTYP57(nr)=2
        ENDIF
!newe
      ENDIF

      CALL EXITS('IPPRES')
      RETURN
 9999 CALL ERRORS('IPPRES',ERROR)
      CALL EXITS('IPPRES')
      RETURN 1
      END


      SUBROUTINE PFRF50(IBT,IDO,INP,IXF,NAN,NBH,NBJ,NGAP,NJE,NPF,nr,
     '  PF,PG,RF,WG,XE,XG,XG1,ZG,ERROR,*)

C#### Subroutine: PFRF50
C###  Description:
C###    Evaluates contribution RF(ns,nj) to element residuals RE from
C###    the pressure PF(i) acting on the Xi(3)=IXF face (IXF=iface-1).

C**** Note: AZL & AZU are deformed state metric tensors wrt Nu.
C**** 21Sep88: This version of PFRF differs from earlier version only in
C**** that calculations are wrt Nu instead of Xi-coords. This way is
C**** less efficient. Also, this version calls the subroutine DLZJDX

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IXF,
     '  NAN(NIM,NAM,*),NBH(*),NBJ(*),NGAP(NIM,*),NJE,NPF(*),nr
      REAL*8 PF,PG(NSM,NUM,NGM,*),RF(32,*),WG(NGM,*),
     '  XE(NSM,*),XG(NJM,*),XG1(NJM,*),ZG(NHM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,NBFF,ng,ng1,ng2,NGI1,NGI2,ni,NITB,nj,njj,ns,NU1(0:3)
      REAL*8 AZ,AZL(3,3),AZU(3,3),D(5,5),DXIX(3,3),DXIXJ(3,3),DZI(3,3),
     '  GXL(3,3),GXU(3,3),RAZ,RWG,RWX,SUM,XI(3)
      DATA NU1/1,2,4,7/
      DATA D/5*0.0D0,-0.288675134594813D0,0.288675134594813D0,3*0.0D0,
     '       -0.387298334620741D0,0.0D0,0.387298334620741D0,2*0.0D0,
     '       -0.430568155797026D0,    -0.169990521792428D0,
     '        0.169990521792428D0,     0.430568155797026D0,  0.0D0,
     '       -0.453089922969332D0,    -0.269234655052841D0,  0.0D0,
     '        0.269234655052841D0,     0.453089922969332D0/

      CALL ENTERS('PFRF50',*9999)

      NITB=NIT(NBH(1))

      DO njj=1,NJ_LOC(NJL_GEOM,0)
        nj=NJ_LOC(NJL_GEOM,njj)
        NBFF=NPF(9+nj)
        DO ns=1,NST(NBFF)+NAT(NBFF)
          RF(ns,nj)=0.0D0
        ENDDO !ns
      ENDDO !nj

      XI(3)=DBLE(IXF)
      NGI1=NGAP(1,NPF(10))
      NGI2=NGAP(2,NPF(10))
      ng=0
      DO ng2=1,NGI2
        DO ng1=1,NGI1
          ng=ng+1
          XI(1)=0.5d0+D(ng1,NGI1)
          XI(2)=0.5d0+D(ng2,NGI2)
          CALL ZEZI(0,IBT,IDO,INP,NAN,NBJ,NJM,NJE,
     '      DXIX,XE,XG1,XI,ERROR,*9999)
C ***     Calculate DXIX derivatives of Xi wrt NU coords (IP=1)
          CALL XGMG(1,NITB,NBJ(1),NJE,nr,DXIXJ,GXL,GXU,RWX,XG,
     '      ERROR,*9999)
          CALL ZEZI(1,IBT,IDO,INP,NAN,NBH,NHM,NJE,
     '      DXIX,XE,XG1,XI,ERROR,*9999)
          CALL ZGMG(NBH(1),NJE,AZ,AZL,AZU,ZG,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '(/'' >>>PFRF50 diagnostic op at Gauss pt '',I2)') ng
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO njj=1,NJ_LOC(NJL_GEOM,0)
              nj=NJ_LOC(NJL_GEOM,njj)
              WRITE(OP_STRING,'(''  XG('',I1,'',ni): '',4E12.4)') 
     '          nj,(XG(nj,NU1(ni)),ni=0,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''  ZG('',I1,'',ni): '',4E12.4)') 
     '          nj,(ZG(nj,NU1(ni)),ni=0,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nj
            DO mi=1,NITB
              WRITE(OP_STRING,'('' AZU('',I1,'',ni): '',3E12.4)') 
     '          MI,(AZU(mi,ni),ni=1,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !mi
          ENDIF
          RAZ=DSQRT(AZ)
          RWG=RAZ*WG(ng,NPF(10))

C ***     Calculate covariant derivs of deformed theta wrt undef Xj
          CALL DLZJDX(ITYP10(1),NBFF,NJE,DZI,XG,ZG,ERROR,*9999)

          DO njj=1,NJ_LOC(NJL_GEOM,0)
            nj=NJ_LOC(NJL_GEOM,njj)
            NBFF=NPF(9+nj)
            SUM=0.0D0
            DO ni=1,NITB
              SUM=SUM+AZU(ni,3)*DZI(nj,ni)
            ENDDO !ni
            DO ns=1,NST(NBFF)+NAT(NBFF)
              IF(IXF.EQ.1)THEN           !positive face
                RF(ns,nj)=RF(ns,nj)-PF*SUM*PG(ns,1,ng,NBFF)*RWG
              ELSE                       !negative face
                RF(ns,nj)=RF(ns,nj)+PF*SUM*PG(ns,1,ng,NBFF)*RWG
              ENDIF
            ENDDO !ns
          ENDDO !nj
        ENDDO !ng1
      ENDDO !ng2

      CALL EXITS('PFRF50')
      RETURN
 9999 CALL ERRORS('PFRF50',ERROR)
      CALL EXITS('PFRF50')
      RETURN 1
      END


      SUBROUTINE UPVOROFROMDEFORMED(NBH,NHE,NKE,NPF,NPNE,nc,nx,
     '  NVHE,NW,CURVCORRECT,SE,ZA,ZP,IBT,IDO,INP,NPNODE,
     '  XIP,XP,NEP,ERROR,*)


C#### Subroutine: UPVOROFROMDEFORMED
C###  Description:
C###   UPVOROFROMDEFORMED  uses ZPZE and PXI to update the position of
C###   nodes in region 2 from the deformed geometry of region 1.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBH(NHM,NCM,NEM),nc,NEP(NPM), 
     '  NHE(NEM,NXM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NW(NEM,3,NXM),nx,NPNODE(0:NP_R_M,0:NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),XIP(NIM,NPM),
     ' ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ni,NITB,nonode,np,nj,njj
      REAL*8 XI(3),ZE(NSM,NHM)
!     Functions
      REAL*8 PXI

      CALL ENTERS('UPVOROFROMDEFORMED',*9999)

!Temp ne on next line
      nb=NBH(1,1,1)
      NITB=NIT(nb)
      DO nonode=1,NPNODE(0,2)
        np=NPNODE(nonode,2)
        ne=NEP(np)
        nb=NBH(1,1,ne)
        CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),1,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '    CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '    ZE,ZP,ERROR,*9999)
        DO ni=1,NITB
          XI(ni)=XIP(ni,np)
        ENDDO
        DO njj=1,NJ_LOC(NJL_GEOM,0,2)
          nj=NJ_LOC(NJL_GEOM,njj,2)
          XP(1,1,nj,np)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      1,XI,ZE(1,nj))
        ENDDO
      ENDDO

      CALL EXITS('UPVOROFROMDEFORMED')
      RETURN
 9999 CALL ERRORS('UPVOROFROMDEFORMED',ERROR)
      CALL EXITS('UPVOROFROMDEFORMED')
      RETURN 1
      END


Module FE70
=========== 

      SUBROUTINE IOVSAE(IBT,IDO,INP,NBJ,
     '  NKE,NPF,NPNE,NRE,NVJE,
     '  SE,XA,XE,XP,ERROR,*)

C#### Subroutine: IOVSAE
C###  Description:
C###    IOVSAE writes output file for Vsaero input.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:vsa00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vsc00.cmn'
      INCLUDE 'cmiss$reference:vsd00.cmn'
      INCLUDE 'cmiss$reference:vse00.cmn'
      INCLUDE 'cmiss$reference:vsf00.cmn'
      INCLUDE 'cmiss$reference:vsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NKE(NKM,NNM,NBFM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XI(2),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IPRNAB,ITOT1,ITOT2,ix1,ix2,j,nb,nc,ne,nee,nj,
     '  noele,nopat,nosect,noseps,NTELE
      REAL*8 BX,BY,BZ,PXI,XYZ(3)

      CALL ENTERS('IOVSAE',*9999)
      nc=1 ! Temporary MPN 12-Nov-94

C *** Card 1
      WRITE(IOFILE2,*) 'VSaero data from CMISS'

C *** Card 2
      WRITE(IOFILE2,'(6I5)') IPRI,IPRLEV,IPRESS,MSTOP,MSTART,MODIFY
      IF(IPRLEV.EQ.5) THEN
C ***   Card 2a
        WRITE(IOFILE2,'(5I5)') IPRGOM,IPRNAB,IPRWAK,IPRCPV,IPRPPI
      ENDIF

C *** Card 3
      WRITE(IOFILE2,'(8I5)') MODE,NPNMAX,NRBMAX,ITGSMX,IMERGE,NSUB,
     '  NSPMAX,NPCMAX

C *** Card 3a
C     IF(NRBMAX.LT.0) THEN
C       WRITE(IOFILE2,'(16I5)') (NROWB(i),i=1,IABS(NRBMAX))
C     ENDIF

C *** Card 4
      IF(MODE.EQ.1) THEN
        WRITE(IOFILE2,'(3I5)') NWIT,NVPI,IBLTYP
      ELSE IF(MODE.EQ.2) THEN
        WRITE(IOFILE2,'(2I5)') nt,NHC
      ENDIF
C *** Card set 4a
C     IF(NVPI.GT.0.AND.IBLTYP.EQ.0) THEN
C       WRITE(9,'(I5)') NPSETS
C       DO 100 n=1,NPSETS
C         WRITE(9,'(16I5)') NPCHBL,NBCOL,(KOL(i),i=1,NBCOL)
C100    CONTINUE
C     ENDIF

C *** Card 5
      WRITE(IOFILE2,'(7F10.2)') RSYM,RGPR,RNF,RFF,RCORE,SOLRES,TOL

C *** Card 6
      WRITE(IOFILE2,'(5F10.2)') ALDEG,YAWDEG,RMACH,VMOD,COMFAC
C *** Card 6a
      IF(MODE.EQ.2) THEN
        WRITE(IOFILE2,'(5F10.2)') ALBAR,RFREQ,HX,HY,HZ
      ENDIF

C *** Card 7
      WRITE(IOFILE2,'(6F10.2)') CBAR,SREF,SSPAN,RMPX,RMPY,RMPZ

C *** Card 8
      WRITE(IOFILE2,'(5I5)') NORSET,NVORT,NPASUM,JETPAN,NBCHGE
C *** Card set 8a
C     IF(NORSET.GT.0) THEN
C       WRITE(IOFILE2,'(5I5,2F10.2)') (NORPCH(i),NORF(i),NORL(i),
C    '    NOCF(i),NOCL(i),VNORM(i),ADUB(i),i=1,NORSET)
C     ENDIF
C *** Card set 8b
C     IF(NVORT.GT.0) THEN
C       WRITE(IOFILE2,'(F10.2)') VORT
C       WRITE(IOFILE2,'(3F10.2)') (RXV(i),RYV(i),RZV(i),i=1,NVORT+1)
C     ENDIF
C *** Card set 8c
C     IF(NPASUM.GT.0) THEN
C       WRITE(IOFILE2,'(5I5)') (NPSPCH(i),NPSRF(i),NPSRL(i),NPSCF(i),
C    '    NPSCL(i),i=1,NPASUM)
C     ENDIF
C *** Card set 8d
C     IF(JETPAN.GT.0) THEN
C       WRITE(IOFILE2,'(5I5,2F10.2)') (JETPCH(i),JETRF(i),JETRL(i),
C    '    JETCF(i),JETCL(i),VIN(i),VOUT(i),i=1,JETPAN)
C     ENDIF
C *** Card set 8e
C     IF(NBCHGE.GT.0) THEN
C       WRITE(IOFILE2,'(4I5)') (KPAN(i),KSIDE(i),NEWNAB(i),NEWSID(i),
C    '    I=1,NBCHGE)
C     ENDIF

      DO 600 nopat=1,NTPAT
C ***   Card 9 (Component card)
        WRITE(IOFILE2,'(5F10.2)') CTX(nopat),CTY(nopat),CTZ(nopat),
     '    SCAL(nopat),THET(nopat)
C ***   Card IOFILE2a
        IF(SCAL(nopat).LT.0.0) THEN
          WRITE(IOFILE2,'(6F10.2)') CPX(nopat),CPY(nopat),CPZ(nopat),
     '      CHX(nopat),CHY(nopat),CHZ(nopat)
        ENDIF
C ***   Card 10
        WRITE(IOFILE2,'(4I5)') IDENT(nopat),MAKE(nopat),KOMP(nopat),
     '    KLASS(nopat)
C ***   Card sets 11(section card),12(section def.n) & 14(node card)
C ***   Card sets 11(section card),12(section def.n) & 14(node card)
C ***   Define panels by subdividing each element so that every gauss point
C ***   in the element is the centre of a panel             FMM 12/10/88
C ***   Note Y and Z coords reversed so sail is defined as wing in VSAERO
      	NTELE=NTNODE(nopat)/3
        ITOT2=2
        nosect=0
        DO nee=1,(NET(1)-NTELE+1),NTELE
          IF(nee.EQ.(NET(1)-NTELE+1)) ITOT2=3
      	  DO ix2=0,ITOT2
      	    nosect=nosect+1
            WRITE(IOFILE2,'(6F10.2,4I5)') STX(nopat),STY(nopat),
     '        STZ(nopat),SCALE(nopat),ALF(nopat),THETA(nopat),
     '        INMODE(nopat),NODES(nosect,nopat),NPS(nopat)
            XI(2)=DBLE(ix2)*1.0/3.0
            DO noele=1,NTELE
              ne=(nee-1)+noele
              CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              ITOT1=2
              IF(noele.EQ.NTELE)ITOT1=3
      	      DO ix1=0,ITOT1
      	        XI(1)=DBLE(ix1)*1.D0/3.D0
                DO nj=1,3
                  nb=NBJ(nj,ne)
	          XYZ(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XI,XE(1,nj))
        	ENDDO
                BX=XYZ(1)-CTX(nopat)
       	        BY=XYZ(3)-CTY(nopat)
            	BZ=XYZ(2)-CTZ(nopat)
                WRITE(IOFILE2,'(3F10.2)') BX,BY,BZ
              ENDDO
            ENDDO
            WRITE(IOFILE2,'(30X,4I5)') NODEC(nopat),NPC(nopat),
     '        INTC(nopat),MOVE(nopat)
          ENDDO
        ENDDO
 600  CONTINUE

C *** Card set 17 (wake grid plane)
      WRITE(IOFILE2,'(F10.2)') XMIN
      WRITE(IOFILE2,'(F10.2)') XMAX
C *** Card 18 (node card)
      WRITE(IOFILE2,'(30X,4I5)') NODE,NPCW,INTCW,MARK

      DO 700 nopat=1,NTPAT
C ***   Card 19 (wake)
        WRITE(IOFILE2,'(3I5)') IDENTW(nopat),IFLEXW(nopat),IDEFW(nopat)
C ***   Card 20 (separation line)
        WRITE(IOFILE2,'(7I5)') KWPACH(nopat),KWSIDE(nopat),
     '    KWLINE(nopat),KWPAN1(nopat),KWPAN2(nopat),INPUT(nopat),
     '    NODEWS(1,nopat)
C ***   Card 21 (streamwise wake-line geometry)
        DO 650 noseps=1,NTSEPS(nopat)
          WRITE(IOFILE2,'(3F10.2)') SWPX(noseps,nopat),
     '      SWPY(noseps,nopat),DELTAZ(noseps,nopat)
 650    CONTINUE
C ***   Card 22 (wake node card)
        WRITE(IOFILE2,'(30X,3I5)') NODEWC(nopat),NPCP(nopat),
     '    INTCP(nopat)
        WRITE(IOFILE2,'(7I5)') 0,0,0,0,0,0,NODEWS(2,nopat)
 700  CONTINUE

C *** Card set 24 and 25 (surface streamline data)
      DO 800 i = 1,NSTRM
        WRITE(IOFILE2,'(F10.0,2I5)')F(i),KP(i)
800   CONTINUE

      IF(NVPI.GT.0) THEN
C ***   Card 26 (boundary-layer data)
        WRITE(IOFILE2,'(5F10.0)')RND,TRIPUP,TRIPOP,XPRINT,XSKIP
      ENDIF

C *** Off body velocity scan
C *** Card 27 (Scan Box)
      DO 850 i = 1,NBOX
    	  WRITE(IOFILE2,'(6I5)')MOLD(i),MEET(i),NEAR(i),INCPRI(i),
     '    INCPRJ(i),INCPRK(i)
C *** Card 28 (First corner of skewed box)
        WRITE(IOFILE2,'(3F10.0,I5)')X0(i),Y0(i),Z0(i),np(i)
        IF(np(i).GT.1)THEN
C *** Card 29 (Second Corner of skewed box)
          WRITE(IOFILE2,'(3F10.0,I5)')X1(i),Y1(i),Z1(i),NP1(i)
        	IF(NP1(i).LT.0)THEN
C *** Card 29a ( Point location along first edge) 	
            WRITE(IOFILE2,'(8F10.0)')(AL1(I,J),j=1,ABS(NP1(i)))
          ENDIF
        ENDIF
        IF(np(i).GT.2)THEN
C *** Card 30 (Third corner point of skewed box)
          WRITE(IOFILE2,'(3F10.0,I5)')X2(i),Y2(i),Z2(i),NP2(i)
          IF(NP2(i).LT.0)THEN
C *** Card 30a (Location of scan lines)
            WRITE(IOFILE2,'(8F10.0)')(AL2(i,j),j=1,ABS(NP2(i)))
          ENDIF
        ENDIF
        IF(np(i).EQ.4)THEN
C *** Card 31 (Fourth corner point of skewed box)
          WRITE(IOFILE2,'(3F10.0,I5)')X3(i),Y3(i),Z3(i),NP3(i)
          IF(NP3(i).LT.0)THEN
C *** Card 31a (Locations of Scan lines)
            WRITE(IOFILE2,'(8F10.0)')(AL3(i,j),j=1,ABS(NP3(i)))
          ENDIF
        ENDIF
850   CONTINUE

C *** Cards 32-34 Cylindrical Scan volumes

C ***  Off-body streamline data
C *** Card 35 (Location of staring point for streamline calc)
      DO 900 i=1,NPOINT
        WRITE(IOFILE2,'(6F10.0,I5)')RSX(i),RSY(i),RSZ(i),SU(i),SD(i),
     '    DELS(i),NEAR2(i)
900   CONTINUE
       	
      CALL EXITS('IOVSAE')
      RETURN
 9999 CALL ERRORS('IOVSAE',ERROR)
      CALL EXITS('IOVSAE')
      RETURN 1
      END


      SUBROUTINE IPVSA(ERROR,*)

C#### Subroutine: IPVSA
C###  Description:
C###    IPVSA performs basic input for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsa00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPVSA',*9999)
      FILEIP=.FALSE.
      NOQUES=0

C *** Card 2
      FORMAT='('' Enter input print control [0]:'''//
     '  '/''   (0) Print all i/p data except patch geometry'''//
     '  '/''   (1) Print all i/p data except detailed coordinates of'//
     '             ' patch geometry'''//
     '  '/''   (2) Print all i/p data'''//
     '  '/$,''    '',I1)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,0,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,2,
c    '  INFO,ERROR,*9999)
      IPRI=IDATA(1)

      FORMAT='('' Enter output print control [0]:'''//
     '  '/''   (0) Basic print level'''//
     '  '/''   (1) + panel corner points'''//
     '  '/''   (2) + doublet solution'''//
     '  '/''   (3) + p/p(infinity)'''//
     '  '/''   (4) + corner point analysis'''//
     '  '/''   (5) + further options'''//
     '  '/$,''    '',I1)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,0,5,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,5,
c    '  INFO,ERROR,*9999)
      IPRLEV=IDATA(1)

      FORMAT='('' Enter printout frequency in wake iteration or'//
     '          ' time-stepping [0]:'''//
     '  '/''   (0) Prints last step only'''//
     '  '/''   (1) Prints every step'''//
     '  '/''   (n) Prints at every nth step'''//
     '  '/$,''   '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,0,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,99,
c    '  INFO,ERROR,*9999)
      IPRESS=IDATA(1)

      FORMAT='('' Enter extent of calculations [0]:'''//
     '  '/''   (0) Complete run'''//
     '  '/''   (1) Stop after GEOMIN (i/p file only requires basic'//
     '           ' data & patch geometry)'''//
     '  '/''   (2) Stop after SURPAN (i/p file only requires basic'//
     '           ' data & patch geometry)'''//
     '  '/''   (3) Stop after WAKPAN (i/p file includes wake'//
     '           ' geometry)'''//
     '  '/''   (4) Complete run with restart file formed after'//
     '           ' ANALIZ'''//
     '  '/$,''    '',I1)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,0,4,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,4,
c    '  INFO,ERROR,*9999)
      MSTOP=IDATA(1)

      FORMAT='('' Enter restart control [0]:'''//
     '  '/''   (0) Regular first run'''//
     '  '/''   (3) Program restart for further solns & wake shape'//
     '           ' iterations &/or bdry layer calculations'''//
     '  '/''   (4) Program restart for off-body velocities'//
     '           ' & streamlines using earlier solution'''//
     '  '/''  *(5) Program restart for more surface streamlines'//
     '           ' & bdry layer calcs.'''//
     '  '/$,''    '',I1)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,0,5,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,5,
c    '  INFO,ERROR,*9999)
      MSTART=IDATA(1)

      IF(MSTART.GT.0) THEN
        FORMAT='('' Enter details of restart run [0]:'''//
     '    '/''   (0) no change in basic conditions'''//
     '    '/''  (>0) Original basic i/p + restart changes'''//
     '    '/''   (2) Wake input included in changes'''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,2,
c    '    INFO,ERROR,*9999)
        MODIFY=IDATA(1)
      ENDIF

      IF(IPRLEV.EQ.5) THEN
C ***   Card 2a
        FORMAT='('' Enter additional print control [0]:'''//
     '    '/''      (0) Print off'''//
     '    '/''      (1) Panel corner points printed for all panels'''//
     '    '/''  (200+n) Panel corner points printed for panels on'//
     '                ' patch n'''//
     '    '/''      (2) Panel control points & unit normal vectors'//
     '                ' printed for all panels'''//
     '    '/''  (400+n) Panel control points & unit normal vectors'//
     '                ' printed for panels on patch n'''//
     '    '/''     (-1) Prints basic points (input & generated) for'//
     '                ' defined sections on all patches'''//
     '    '/'' -(200+n) Prints basic points (input & generated) for'//
     '                ' patch n'''//
     '    '/$,''    '',I4)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,-IMAX,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,-IMAX,IMAX,
c    '    INFO,ERROR,*9999)
        IPRGOM=IDATA(1)

        FORMAT='('' Enter printout control for panel neighbour info'//
     '            ' [0]:'''//
     '    '/''   (0) Off'''//
     '    '/''   (1) Print info for panels on patch edges, panels at'//
     '             ' wake shedding lines'//
     '    '/''       & panels which have failed to find neighbours'''//
     '    '/''   (2) Print neighbour info for all panels'''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,2,
c    '    INFO,ERROR,*9999)
        IPRNAD=IDATA(1)

        FORMAT='('' Enter printout control for wake data [0]:'''//
     '    '/''   (0) Off'''//
     '    '/''   (1) Print wake-shedding info for each wake column'''//
     '    '/''   (2) + details of wake line geometry'''//
     '    '/''   (3) + wake panel doublet values'''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,3,
c    '    INFO,ERROR,*9999)
        IPRWAK=IDATA(1)

        FORMAT='('' Enter printout control for panel corner point'//
     '            ' analysis [0]:'''//
     '    '/''   (0) Off'''//
     '    '/''   (1) Print x,y,z,vx,vy,vz,V & CP for each panel'//
     '             ' corner point'''//
     '    '/''   (2) + panel corner point doublet & source values'''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,2,
c    '    INFO,ERROR,*9999)
        IPRCPV=IDATA(1)

        FORMAT='('' Enter printout control for p/p(infinity) values'//
     '            ' & veloc & press @ panel centres [0]:'''//
     '    '/''   (0) Off'''//
     '    '/''   (1) On'''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,1,
c    '    INFO,ERROR,*9999)
        IPRPPI=IDATA(1)
      ENDIF

C *** Card 3
      FORMAT='('' Enter mode [1]:'''//
     '  '/''   (1) Steady calculation'''//
     '  '/''   (2) Unsteady calculation'''//
     '  '/$,''    '',I1)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IONE,1,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,2,
c    '  INFO,ERROR,*9999)
      MODE=IDATA(1)

C     IDEFLT(1)=1000
C     FORMAT='($,'' Enter upper limit on number of panels [1000]:'
C    '  //' '',I4)'
C     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,1000,
C    '  INFO,ERROR,*9999)
C     NPNMAX=IDATA(1)
      NPNMAX=9*NET(1)

      IDEFLT(1)=140
      FORMAT='($,'' Enter limit for block size in Gauss-Seidel [140]:'
     '  //' '',I3)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IDEFLT,1,140,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,140,
c    '  INFO,ERROR,*9999)
      NRBMAX=IDATA(1)

      IDEFLT(1)=20
      FORMAT='($,'' Enter limit on no of Gauss-Seidel iterations [20]:'
     '  //' '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IDEFLT,1,20,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,20,
c    '  INFO,ERROR,*9999)
      ITGSMX=IDATA(1)

      IMERGE=0
      IDEFLT(1)=10
      FORMAT='($,'' Enter no of subpanel intervals used on a'
     '  //' near-field wake panel [10]: '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IDEFLT,1,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,100,
c    '  INFO,ERROR,*9999)
      NSUB=IDATA(1)

      IDEFLT(1)=25
      FORMAT='($,'' Enter limit on no of subpanels used on a'
     '  //' near-field wake panel [25]: '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IDEFLT,1,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,100,
c    '  INFO,ERROR,*9999)
      NSPMAX=IDATA(1)

      IF(MODE.EQ.1) THEN
        IDEFLT(1)=2
        FORMAT='($,'' Enter limit on no of pred./corr. cycles'
     '    //' in steady wake relax.n [2]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,1,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,100,
c    '    INFO,ERROR,*9999)
        NPCMAX=IDATA(1)
      ELSE IF(MODE.EQ.2) THEN
        NPCMAX=0
      ENDIF

C *** Card 4
      IF(MODE.EQ.1) THEN
        IDEFLT(1)=0
        FORMAT='($,'' Enter no of wake shape iter.ns/pot.flow iter.n'
     '    //' [0=rigid wake]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,100,
c    '    INFO,ERROR,*9999)
        NWIT=IDATA(1)

        IDEFLT(1)=0
        FORMAT='($,'' Enter no of viscous iter.ns/pot.flow iter.n'
     '    //' [0=pot.flow only]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,100,
c    '    INFO,ERROR,*9999)
        NVPI=IDATA(1)

        IF(NVPI.GT.0) THEN
          IBLTYP=1
        ENDIF
      ELSE IF(MODE.EQ.2) THEN
        IDEFLT(1)=0
        FORMAT='($,'' Enter no of time steps [0]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,99,
c    '    INFO,ERROR,*9999)
        nt=IDATA(1)

        IDEFLT(1)=0
        FORMAT='($,'' Enter no of half cycles [0]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,100,
c    '    INFO,ERROR,*9999)
        NHC=IDATA(1)
      ENDIF

C *** Card 5
      IDEFLT(1)=1
      FORMAT='($,'' Enter (0)symmetry about y=0, (1)asymmetry [1]:'
     '  //' '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IDEFLT,0,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,1,
c    '  INFO,ERROR,*9999)
      RSYM=DBLE(IDATA(1))

      IDEFLT(1)=1
      FORMAT='($,'' Enter (0)free air, (1)ground plane @ z=0 [1]:'
     '  //' '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IDEFLT,0,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,1,
c    '  INFO,ERROR,*9999)
      RGPR=DBLE(IDATA(1))

      RDEFLT(1)=2.5D0
      FORMAT='($,'' Enter radius of near-field factor [2.5]:'
     '  //' '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.d0,RMAX,
c    '  INFO,ERROR,*9999)
      RNF=RDATA(1)

      RDEFLT(1)=5.D0
      FORMAT='($,'' Enter radius of far-field factor [5.0]:'
     '  //' '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      RFF=RDATA(1)

      RDEFLT(1)=0.05D0
      FORMAT='($,'' Enter core radius on vortex filaments [0.05]:'
     '  //' '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      RCORE=RDATA(1)

      RDEFLT(1)=0.2D0
      FORMAT='($,'' Enter Gauss-Seidel residual limit as % of max'
     '  //' doublet size [0.2]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      SOLRES=RDATA(1)

      RDEFLT(1)=0.2D0
      FORMAT='($,'' Enter tolerance limit for test of proximity of a'
     '  //' point to a panel edge [0.2]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      TOL=RDATA(1)

C *** Card 6
      RDEFLT(1)=0.D0
      FORMAT='($,'' Enter incidence of x-axis in degrees [0.]: '','
     '  //'E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      ALDEG=RDATA(1)

      RDEFLT(1)=0.D0
      FORMAT='($,'' Enter yaw of x-axis in degrees [0.]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      YAWDEG=RDATA(1)

      RDEFLT(1)=0.D0
      FORMAT='($,'' Enter mach number [0.]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      RMACH=RDATA(1)

      RDEFLT(1)=1.D0
      FORMAT='($,'' Enter onset flow velocity magnitude [1.]:'
     '  //' '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      VMOD=RDATA(1)

      RDEFLT(1)=0.D0
      FORMAT='($,'' Enter compressibility algorithm factor [0.]:'
     '  //' '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      COMFAC=RDATA(1)

C *** Card 6a
      IF(MODE.EQ.2) THEN
        RDEFLT(1)=0.D0
        FORMAT='($,'' Enter amplitude of motion in degrees [0.]:'
     '    //' '',E11.4)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '    INFO,ERROR,*9999)
        ALBAR=RDATA(1)

        RDEFLT(1)=0.D0
        FORMAT='($,'' Enter reduced frequency [0.]: '',E11.4)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '    INFO,ERROR,*9999)
        RFREQ=RDATA(1)

        RDEFLT(1)=0.D0
        RDEFLT(2)=1.D0
        RDEFLT(3)=0.D0
        FORMAT='($,'' Enter pivot axis unit vector thru ref. moment'
     '    //' point [0.,1.,0.]: '',3E11.4)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.d0,1.d0,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RDEFLT,0.,1.,
c    '    INFO,ERROR,*9999)
        HX=RDATA(1)
        HY=RDATA(2)
        HZ=RDATA(3)
      ENDIF

C *** Card 7
      RDEFLT(1)=1.D0
      FORMAT='($,'' Enter reference chord used for normalizing'
     '  //' pitching moment [1.]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      CBAR=RDATA(1)

      RDEFLT(1)=1.D0
      FORMAT='($,'' Enter reference area [1.]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      SREF=RDATA(1)

      RDEFLT(1)=1.D0
      FORMAT='($,'' Enter semi-span for normalizing rolling & yawing'
     '  //' moments [1.]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      SSPAN=RDATA(1)

      RDEFLT(1)=0.D0
      RDEFLT(2)=0.D0
      RDEFLT(3)=0.D0
      FORMAT='($,'' Enter coords of ref. moment point'
     '  //' (if unsteady must be on pivot axis)[0,0,0]: '',3E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RDEFLT,0.,RMAX,
c    '  INFO,ERROR,*9999)
      RMPX=RDATA(1)
      RMPY=RDATA(2)
      RMPZ=RDATA(3)

C *** Card 8
      NORSET=0
      NVORT =0
      NPASUM=0
      JETPAN=0
      NBCHGE=0

      CALL EXITS('IPVSA')
      RETURN
 9999 CALL ERRORS('IPVSA',ERROR)
      CALL EXITS('IPVSA')
      RETURN 1
      END


      SUBROUTINE IPVSB(ERROR,*)

C#### Subroutine: IPVSB
C###  Description:
C###    IPVSB performs patch geometry input for VSAERO. Note: Currently
C###    dimensioned for 9 patches (sails),10 sections per patch and 20 
C###    nodes per section.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,INFO,nonode,nopat,NOQUES,nosect
      LOGICAL FILEIP
      CHARACTER CHAR*11

      CALL ENTERS('IPVSB',*9999)
      FILEIP=.FALSE.
      NOQUES=0

      IDEFLT(1)=1
      FORMAT='($,'' Enter the number of patches [1]: '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IDEFLT,1,10,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,10,
c    '  INFO,ERROR,*9999)
      NTPAT=IDATA(1)

      DO 600 nopat=1,NTPAT
        WRITE(CHAR,'(I1)') nopat    

C ***   Card 9 (Component card)
        RDEFLT(1)=0.D0
        RDEFLT(2)=0.D0
        RDEFLT(3)=0.D0
        FORMAT='($,'' Enter coordinates of patch '//CHAR(1:1)//
     '    ' origin [0,0,0]: '',3E11.4)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RDEFLT,-RMAX,RMAX,
c    '    INFO,ERROR,*9999)
        CTX(nopat)=RDATA(1)
        CTY(nopat)=RDATA(2)
        CTZ(nopat)=RDATA(3)

        RDEFLT(1)=1.D0
        FORMAT='($,'' Enter scale factor on patch (-ve if rotation'
     '    //' not about default y-axis)[1.]: '',E11.4)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,-RMAX,RMAX,
c    '    INFO,ERROR,*9999)
        SCAL(nopat)=RDATA(1)

        RDEFLT(1)=0.0
        FORMAT='($,'' Enter rotation angle in degrees [0.]: '',E11.4)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,-RMAX,RMAX,
c    '    INFO,ERROR,*9999)
        THET(nopat)=RDATA(1)

C ***   Card 9a
        IF(SCAL(nopat).LT.0.D0) THEN
          CPX(nopat)=0.D0
          CPY(nopat)=0.D0
          CPZ(nopat)=0.D0

          RDEFLT(1)=0.D0
          RDEFLT(2)=0.D0
          RDEFLT(3)=1.D0
          FORMAT='($,'' Enter coordinates of rotation axis'
     '      //' [0.,0.,1.]: '',3E11.4)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
c         CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RDEFLT,-RMAX,
c    '      RMAX,INFO,ERROR,*9999)
          CHX(nopat)=RDATA(1)
          CHY(nopat)=RDATA(2)
          CHZ(nopat)=RDATA(3)

        ELSE IF(SCAL(nopat).GE.0.D0) THEN
          CHX(nopat)=0.D0
          CHY(nopat)=1.D0
          CHZ(nopat)=0.D0
        ENDIF

C ***   Card 10
        IDEFLT(1)=3
        FORMAT='('' Enter patch option [3]:'''//
     '    '/''   (1) Wing-type patch   '''//
     '    '/''   (2) Body-type patch   '''//
     '    '/''   (3) Single sheet patch'''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,3,
c     '   INFO,ERROR,*9999)
        IDENT(nopat)=IDATA(1)
        MAKE(nopat)=0
        KOMP(nopat)=nopat
        KLASS(nopat)=1

        FORMAT='($,'' Enter number of sections in patch [1]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IONE,1,IOIMX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,IOIMX,
c    '    INFO,ERROR,*9999)
        NTSECT(nopat)=IDATA(1)

C ***   Card sets 11(section card),12(section def.n) & 14(node card)
        FORMAT='($,'' Enter number of nodes per section [1]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IONE,1,IOIMX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,IOIMX,
c    '    INFO,ERROR,*9999)
        NTNODE(nopat)=IDATA(1)

        STX(nopat)=0.D0
        STY(nopat)=0.D0
        STZ(nopat)=0.D0
        SCALE(nopat)=1.D0
        ALF(nopat)  =0.D0
        THETA(nopat)=0.D0
        INMODE(nopat)=4

        DO 570 nosect=1,NTSECT(nopat)
          WRITE(CHAR,'(I2)') nosect
          CALL TRIM(CHAR,IBEG,IEND)
          DO 560 nonode=1,NTNODE(nopat)
            IDEFLT(nonode)=0
 560      CONTINUE
          FORMAT='($,'' Enter nodes in section '//CHAR(IBEG:IEND)//
     '      ': '',20I3)'
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,NTNODE(nopat),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
c         CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,NTNODE(nopat),IDATA,
c    '      IDEFLT,1,NPM,INFO,ERROR,*9999)
          DO 565 nonode=1,NTNODE(nopat)
            NODS(nonode,nosect,nopat)=IDATA(nonode)
 565      CONTINUE
          IF(nosect.LT.NTSECT(nopat)) THEN
            NODES(nosect,nopat)=0
          ELSE
            IF(nopat.LT.NTPAT) THEN
              NODES(nosect,nopat)=4
            ELSE IF(nopat.EQ.NTPAT) THEN
              NODES(nosect,nopat)=5
            ENDIF
          ENDIF
          NPS(nopat)=0

          NODEC(nopat)=3
          NPC(nopat)  =0
          INTC(nopat) =3
          MOVE(nopat) =0
 570    CONTINUE
 600  CONTINUE

      CALL EXITS('IPVSB')
      RETURN
 9999 CALL ERRORS('IPVSB',ERROR)
      CALL EXITS('IPVSB')
      RETURN 1
      END


      SUBROUTINE IPVSC(XP,ERROR,*)

C#### Subroutine: IPVSC
C###  Description:
C###    IPVSC performw wake input for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vsc00.cmn'
!     Parameter List
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INDENTW,INFO,nopat,np,NOQUES,noseps,nv
      LOGICAL FILEIP
      CHARACTER CHAR*11

      CALL ENTERS('IPVSC',*9999)
      FILEIP=.FALSE.
      NOQUES=0

      nv=1 ! Temporary MPN 12-Nov-94
C *** Card set 17 (wake grid plane)
      XMIN=XP(1,nv,1,1)
      DO np=2,NPT(1)
        IF(XP(1,nv,1,np).LT.XMIN) XMIN=XP(1,nv,1,np)
      ENDDO
      RDEFLT(1)=XMIN
      WRITE(CHAR,'(E11.4)') XMIN   
      FORMAT='($,'' Enter x-coord of upstream wake grid plane'''
     '  //',/''  (must be upstream of all separation lines)'''
     '  //',/$,''  ['//CHAR(1:11)//'=lowest x-coord in mesh]: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,-RMAX,RMAX,
c    '  INFO,ERROR,*9999)
      XMIN=RDATA(1)

      XMAX=XP(1,nv,1,1)
      DO np=2,NPT(1)
        IF(XP(1,nv,1,np).GT.XMAX) XMAX=XP(1,nv,1,np)
      ENDDO
      RDEFLT(1)=XMAX+4.0*(XMAX-XMIN)
      WRITE(CHAR,'(E11.4)') RDEFLT(1)    
      FORMAT='($,'' Enter x-coord of downstream wake grid plane'''
     '  //',/$,''  ['//CHAR(1:11)//']: '',E11.4)'
      CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
c     CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,-RMAX,RMAX,
c    '  INFO,ERROR,*9999)
      XMAX=RDATA(1)

C *** Card 18 (node card)
      NODE=3
      FORMAT='($,'' Enter no of additional wake grid planes'
     '  //' (max 29)[0]: '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,1,29,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,1,29,
c    '  INFO,ERROR,*9999)
      IF(IDATA(1).GT.0) THEN
        NPCW=IDATA(1)+1
      ELSE
        NPCW=0
      ENDIF

      IF(NPCW.GT.0) THEN
        IDEFLT(1)=3
        FORMAT='('' Enter type of spacing [3]:'''//
     '    '/''   (0) Full cosine                            '''//
     '    '/''   (1) Half-cosine (smaller panels upstream)  '''//
     '    '/''   (2) Half-cosine (smaller panels downstream)'''//
     '    '/''   (3) Equal                                  '''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,3,
c    '    INFO,ERROR,*9999)
        INTCW=IDATA(1)
      ELSE
        INTCW=3
      ENDIF
      MARK=0

      DO 600 nopat=1,NTPAT

C ***   Card 19 (wake)
        WRITE(CHAR,'(I1)') nopat
        IDEFLT(1)=1
        FORMAT='('' Enter type of wake for patch '//CHAR(1:1)//
     '    ' [1]:'''//
     '    '/''   (0) no wake'''//
     '    '/''   (1) Regular wake'''//
     '    '/''   (2) Unsteady wake'''//
     '    '/''   (3) Separated wake'''//
     '    '/''  *(4) Jet model'''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,4,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,4,
c    '    INFO,ERROR,*9999)
        IDENTW(nopat)=IDATA(1)

        FORMAT='($,'' Specify whether wake is (0)flexible or (1)rigid'
     '    //' [0]: '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,1,
c    '    INFO,ERROR,*9999)
        IFLEXW(nopat)=IDATA(1)
        IDEFW(nopat)=0

C ***   Card 20 (separation line)
        KWPACH(nopat)=nopat
        IDEFLT(1)=2
        FORMAT='('' Enter patch side parallel to separation'//
     '    ' line [2]:'''//
     '    '/''   (1) Foot  (inside)       '''//
     '    '/''   (2) Leech (trailing edge)'''//
     '    '/''   (3) Head  (tip)          '''//
     '    '/''   (4) Luff  (leading edge) '''//
     '    '/$,''    '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,1,4,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,4,
c    '    INFO,ERROR,*9999)
        KWSIDE(nopat)=IDATA(1)

        FORMAT='($,'' Enter no of element sides between patch side'//
     '    ' & separ.n line [0]: '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,20,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,20,
c    '    INFO,ERROR,*9999)
        KWLINE(nopat)=IDATA(1)
        KWPAN1(nopat)=0
        KWPAN2(nopat)=0
        INPUT(nopat) =2
        NODEWS(1,nopat)=0

C ***   Card 21 (streamwise wake-line geometry)
        noseps=0
 500    CONTINUE
        RDEFLT(1)=-1.D6
        FORMAT='($,'' Enter x,z coords of a point on wake line'','//
     '    '/$,''  (but not separation point x=0,z=0)[exit]: '',2E11.3)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,2,RDATA,RDEFLT,-RMAX,RMAX,
c    '    INFO,ERROR,*9999)
        IF(RDATA(1).GT.-1.D5) THEN
          noseps=noseps+1
          SWPX(noseps,nopat)=RDATA(1)
          SWPY(noseps,nopat)=RDATA(2)
          DELTAZ(noseps,nopat)=0
          GO TO 500
        ENDIF
        NTSEPS(nopat)=noseps

C ***   Card 22 (wake node card)
        NODEWC(nopat)=3
        NPCP(nopat)   =0
        INTCP(nopat)  =0

        IF(nopat.LT.NTPAT) THEN
          NODEWS(2,nopat)=3
        ELSE IF(nopat.EQ.NTPAT) THEN
          NODEWS(2,nopat)=5
        ENDIF

        IF(INDENTW.EQ.4) THEN
        ENDIF

 600  CONTINUE

      CALL EXITS('IPVSC')
      RETURN
 9999 CALL ERRORS('IPVSC',ERROR)
      CALL EXITS('IPVSC')
      RETURN 1
      END


      SUBROUTINE IPVSD(ERROR,*)

C#### Subroutine: IPVSD
C###  Description:
C###    IPVSD performs surface streamline input for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsd00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,j,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPVSD',*9999)
      FILEIP=.FALSE.
      NOQUES=0

C *** Card set 24 (starting points of each streamline)
      j = 0
 100  CONTINUE
      FORMAT='($,'' Enter panel number for streamline [exit]: '',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,1,99,
c    '  INFO,ERROR,*9999)
      IF(IDATA(1).GT.0) THEN
        j = j+1
        KP(j)=IDATA(1)
        RDEFLT(1)=0.5D0
        FORMAT='($,'' Enter location of streamline in panel'//
     '    ' (.05<.<.95) [0.5]: '',F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.05d0,0.95d0,INFO,ERROR,*9999)
c       CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.05,0.95,
c    '    INFO,ERROR,*9999)
        F(j)=RDATA(1)
        GO TO 100
      ELSE
C ***   Card 25 (last card of streamline input)
        NSTRM = j+1
        F(NSTRM)  = 2.D0
        KP(NSTRM) = 0.D0
      ENDIF

      CALL EXITS('IPVSD')
      RETURN
 9999 CALL ERRORS('IPVSD',ERROR)
      CALL EXITS('IPVSD')
      RETURN 1
      END


      SUBROUTINE IPVSE(ERROR,*)

C#### Subroutine: IPVSE
C###  Description:
C###    IPVSE performs boundary layer input for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsa00.cmn'
      INCLUDE 'cmiss$reference:vse00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPVSE',*9999)
      FILEIP=.FALSE.
      NOQUES=0

      IF(NVPI.GT.0) THEN
C ***   Card 26 (boundary-layer data)
	RDEFLT(1)=1.D0
	FORMAT='($,'' Enter Reynolds number in millions [1]:'',F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.d0,10.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.0,10.0,
c    '    INFO,ERROR,*9999)
        RND=RDATA(1)
	TRIPUP=1.D0
	TRIPOP=0.D0
	FORMAT='($,'' Boundary layer printout, cross flow parameters'//
     '    '/'' [0.0]:      (0.0) Not printed        '''//
     '    '/      ''       (1.0) Printed            '''//
     '    '/$,''             '',F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RZERO,0.0,1.0,
c    '    INFO,ERROR,*9999)
	XPRINT=RDATA(1)
	RDEFLT(1)=1.D0
	FORMAT='($,'' Number of integral intervals to be skipped'//
     '    ' between printout [1.]:'',F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,1.d0,200.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,1.0,200.0,
c    '    INFO,ERROR,*9999)
	XSKIP=RDATA(1)
      ENDIF

      CALL EXITS('IPVSE')
      RETURN
 9999 CALL ERRORS('IPVSE',ERROR)
      CALL EXITS('IPVSE')
      RETURN 1
      END


      SUBROUTINE IPVSF(ERROR,*)

C#### Subroutine: IPVSF
C###  Description:
C###    IPVSF performs off-body velocity scan input for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsf00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ICHAR,INFO,j,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPVSF',*9999)
      FILEIP=.FALSE.
      NOQUES=0

C *** Card Set 27 (Scan Box)
      j = 0
100   CONTINUE
      FORMAT= '($,'' Enter Scan box number [exit]:'',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,1,99,
c    '  INFO,ERROR,*9999)
      IF(IDATA(1).GT.0)THEN
	j=j+1
	MOLD(j)=1
	IDEFLT(1)=1
	FORMAT='($,'' Control of line intersection routine'//
     '    ' (0)Active (1)Off [1]:'',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,0,1,
c    '    INFO,ERROR,*9999)
	MEET(j)=IDATA(1)
	IDEFLT(1)=-1
	FORMAT='($,'' Control of near field routine in velocity'//
     '    ' calculation [-1]:'''//
     '    '/''           (0) Active                         '''//
     '    '/''           (1) Active only for surface panels '''//
     '    '/''          (-1) Inactive                       '''//
     '    '/$,''             '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,-1,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c	CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,-1,1,
c    '    INFO,ERROR,*9999)
 	NEAR(j)=IDATA(1)
	FORMAT='($,'' Print results for every nth plane,n=0 all'//
     '    ' results [0]:'',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,20,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c	CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,20,
c    '    INFO,ERROR,*9999)
	INCPRI(j)=IDATA(1)
	FORMAT='($,'' Print results for every nth line, n=0 all'//
     '    ' results [0]:'',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,20,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c	CALL IINOUT(IOTYPE,IVDU,FORMAT,1,IDATA,IZERO,0,20,
c    '    INFO,ERROR,*9999)
	INCPRJ(j)=IDATA(1)
	FORMAT='($,'' Print results for every nth point, n=0 all'//
     '    ' results [0]:'',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,20,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c	CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,20,
c    '    INFO,ERROR,*9999)
	INCPRK(j)=IDATA(1)
C ***Card 28 (First corner of skewed box)
	FORMAT='($,'' Enter coordinates,x0,y0,z0, of second corner'//
     '     ' of box [0.0,0.0,0.0]:'',3F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RZERO,0.0,1.0,
c    '    INFO,ERROR,*9999)
	X0(j)=RDATA(1)
	Y0(j)=RDATA(2)
	Z0(j)=RDATA(3)
	IDEFLT(1)=1
	FORMAT='($,'' Specify box type [1]:'''//
     '    '/''         (1) Single point                 '''//
     '    '/''         (2) Single line                  '''//
     '    '/''         (3) Lines within a parallelogram '''//
     '    '/''         (4) Lines within a parallelpiped '''//
     '    '$,''            '',I1)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,1,4,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c	CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,4,
c    '    INFO,ERROR,*9999)
	np(j)=IDATA(1)
	IF(np(j).GT.1)THEN
C *** Card 29 (Second corner of skewed box)
	  FORMAT='($,'' Enter coordinates,x1,y1,z1, of second corner'//
     '      ' of box [0.0,0.0,0.0]:'',3F10.0)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	  CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RZERO,0.0,1.0,
c    '      INFO,ERROR,*9999)
	  X1(j)=RDATA(1)
	  Y1(j)=RDATA(2)
	  Z1(j)=RDATA(3)
	  FORMAT='($,'' Enter number of points,n  along the straight'//
     '      ' line (x0,y0,z0) (x1,y1,z1), (0<n<8) [0]:'''//
     '      '/''      (n) Equally spaced    '''//
     '      '/''      (-n) Specify location '''//
     '      '$,''         '',I2)'
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IZERO,0,10,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
c 	  CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,10,
c    '      INFO,ERROR,*9999)
	  NP1(j)=IDATA(1)
	  IF(NP1(i).LT.0)THEN
C ***  Card 29a (Point locations along the first edge of box)
	  FORMAT='($,'' Enter the normalized locations of each point'//
     '      ' along the first edge of the box'//
     '      '  [0.,0.,0.,0.,0.,0.,0.,0.]:'',8F10.0)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,8,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	  CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,8,RDATA,RZERO,0.0,1.0,
c    '      INFO,ERROR,*9999)
	  DO 200 i = 1,ABS(NP1(j))
             AL1(i,j)=RDATA(i)
200       CONTINUE
	  ENDIF
         ENDIF
	IF(np(j).GT.2)THEN
C ***  Card 30 (Third corner point of skewed box)
	 FORMAT='($,''Enter coordinates, x2,y,z2, of third corner of'//
     '     ' box [0.0,0.0,0.0]:'',3F10.0)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,3,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	 CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RZERO,0.0,1.0,
c    '     INFO,ERROR,*9999)
         X2(j)=RDATA(1)
	 Y2(j)=RDATA(2)
	 Z2(j)=RDATA(3)
	 FORMAT='($,'' Enter the number of points, n along the'//
     '     ' straight line (x0,y0,z0) (x2,y2,z2), (0<n<8) [0]:'''//
     '     '/''      (n) Equally spaced    '''//
     '     '/''     (-n) Specify location '''//
     '     '$,''           '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IZERO,0,10,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c	 CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,10,
c    '     INFO,ERROR,*9999)
	 NP2(j)=IDATA(1)
	 IF(NP2(j).LT.0)THEN
C ***  Card 30a (Specify location of scan lines)
	  FORMAT='($,'' Enter the normalized locations of each point'//
     '      ' along the second edge of the box'//
     '      '  [0.,0.,0.,0.,0.,0.,0.,0.]:'',8F10.0)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,8,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	  CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,8,RDATA,RZERO,0.0,1.0,
c    '      INFO,ERROR,*9999)
	  DO 210 i = 1,ABS(NP2(j))
             AL2(i,j)=RDATA(i)
210       CONTINUE
	  ENDIF
         ENDIF
	IF(np(j).EQ.4)THEN
C ***  Card 31 (Fourth corner point of skewed box)
	 FORMAT='($,''Enter coordinates, x3,y3,z3, of fourth corner'//
     '     ' of box [0.0,0.0,0.0]:'',3F10.0)'
         CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,3,
     '     ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '     LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	 CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RZERO,0.0,1.0,
c    '     INFO,ERROR,*9999)
         X3(j)=RDATA(1)
	 Y3(j)=RDATA(2)
	 Z3(j)=RDATA(3)
	 FORMAT='($,'' Enter the number of points, n along the'//
     '     ' straight line (x0,y0,z0) (x3,y3,z3), (0<n<8) [0]:'''//
     '     '/''      (n) Equally spaced    '''//
     '     '/''     (-n) Specify location '''//
     '     '$,''          '',I2)'
         CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '     FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '     IZERO,0,10,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '     INFO,ERROR,*9999)
c	 CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,10,
c    '     INFO,ERROR,*9999)
	 NP3(j)=IDATA(1)
	 IF(NP3(j).LT.0)THEN
C ***  Card 31a (Specify locations of Scan planes)
	  FORMAT='($,'' Enter the normalized locations of each point'//
     '      ' along the third edge of the box '//
     '      ' [0.,0.,0.,0.,0.,0.,0.,0.]:'',8F10.0)'
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,8,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	  CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,8,RDATA,RZERO,0.0,1.0,
c    '      INFO,ERROR,*9999)
	  DO 220 i = 1,ABS(NP3(j))
             AL3(i,j)=RDATA(i)
220       CONTINUE
	  ENDIF
         ENDIF
         GOTO 100
       ENDIF
       NBOX=j+1
       MOLD(NBOX)=0
      CALL EXITS('IPVSF')
      RETURN
 9999 CALL ERRORS('IPVSF',ERROR)
      CALL EXITS('IPVSF')
      RETURN 1
      END


      SUBROUTINE IPVSG(ERROR,*)

C#### Subroutine: IPVSG
C###  Description:
C###    IPVSG performs off-body streamline input for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vsg00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ICHAR,INFO,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPVSG',*9999)
      FILEIP=.FALSE.
      NOQUES=0

C *** Card 35 (Location of starting point for streamline conditions)
      i=0
100   CONTINUE
      FORMAT='($,'' Enter streamline number, return to finish'//
     '  ' [exit]:'',I2)'
      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '  IZERO,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
c     CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,1,99,
c    '  INFO,ERROR,*9999)
      IF(IDATA(1).GT.0)THEN
	i=i+1
	FORMAT='($,'' Enter x,y,z coordinates of starting point'//
     '    ' [0.0,0.0,0.0]:'',3F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,3,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,3,RDATA,RZERO,0.0,1.0,
c    '    INFO,ERROR,*9999)
	RSX(i)=RDATA(1)
	RSY(i)=RDATA(2)
	RSZ(i)=RDATA(3)
	FORMAT='($,'' Enter distance upstream from starting point'//
     '    ' [0.0]:'',F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.0,1.0,
c    '    INFO,ERROR,*9999)
	SU(i)=RDATA(1)
	FORMAT='($,'' Enter distance downstream from starting point'//
     '    ' [0.0]:'',F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.0,1.0,
c    '    INFO,ERROR,*9999)
	SD(i)=RDATA(1)
        FORMAT='($,'' Enter basic length increment along streamline'//
     '    ' [0.0]:'',F10.0)'
        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,0.d0,1.d0,INFO,ERROR,*9999)
c	CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,0.0,1.0,
c    '    INFO,ERROR,*9999)
	DELS(i)=RDATA(1)
        IDEFLT(1)=-1
        FORMAT='($,'' Controls on near-field routine in velocity'//
     '    ' calculation [-1]:'''//
     '    '/''  (0) Active                   '''//
     '    '/''  (1) Active on surface panels '''//
     '    '/'' (-1) Off                      '''//
     '    '/$,''   '',I2)'
        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,-1,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
c       CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,-1,1,
c    '    INFO,ERROR,*9999)
        NEAR2(i)= IDATA(1)
	GOTO 100
      ENDIF
      NPOINT=i+1
      RSX(NPOINT)  =0.D0
      RSY(NPOINT)  =0.D0
      RSZ(NPOINT)  =0.D0
      SU(NPOINT)   =0.D0
      SD(NPOINT)   =0.D0
      DELS(NPOINT) =0.D0
      NEAR2(NPOINT)=0.D0
      CALL EXITS('IPVSG')
      RETURN
 9999 CALL ERRORS('IPVSG',ERROR)
      CALL EXITS('IPVSG')
      RETURN 1
      END


      SUBROUTINE OPVSA(ERROR,*)

C#### Subroutine: OPVSA
C###  Description:
C###    OPVSA performs basic output for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsa00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('OPVSA',*9999)
      WRITE(OP_STRING,'(''      IPRI = '',I1)') IPRI
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IPRLEV = '',I1)') IPRLEV
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IPRESS = '',I2)') IPRESS
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     MSTOP = '',I1)') MSTOP
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    MSTART = '',I1)') MSTART
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    MODIFY = '',I1)') MODIFY
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IPRGOM = '',I4)') IPRGOM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IPRNAD = '',I1)') IPRNAD
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IPRWAK = '',I1)') IPRWAK
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IPRCPV = '',I1)') IPRCPV
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IPRPPI = '',I1)') IPRPPI
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      MODE = '',I1)') MODE
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    NPNMAX = '',I4)') NPNMAX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    NRBMAX = '',I3)') NRBMAX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    ITGSMX = '',I2)') ITGSMX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IMERGE = '',I1)') IMERGE
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      NSUB = '',I2)') NSUB
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    NSPMAX = '',I2)') NSPMAX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    NPCMAX = '',I2)') NPCMAX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      NWIT = '',I2)') NWIT
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      NVPI = '',I2)') NVPI
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    IBLTYP = '',I1)') IBLTYP
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''        nt = '',I2)') nt
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''       NHC = '',I2)') NHC
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      RSYM = '',E11.4)') RSYM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      RGPR = '',E11.4)') RGPR
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''       RNF = '',E11.4)') RNF
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''       RFF = '',E11.4)') RFF
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     RCORE = '',E11.4)') RCORE
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    SOLRES = '',E11.4)') SOLRES
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''       TOL = '',E11.4)') TOL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     ALDEG = '',E11.4)') ALDEG
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    YAWDEG = '',E11.4)') YAWDEG
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     RMACH = '',E11.4)') RMACH
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      VMOD = '',E11.4)') VMOD
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    COMFAC = '',E11.4)') COMFAC
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     ALBAR = '',E11.4)') ALBAR
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     RFREQ = '',E11.4)') RFREQ
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''        HX = '',E11.4)') HX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''        HY = '',E11.4)') HY
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''        HZ = '',E11.4)') HZ
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      CBAR = '',E11.4)') CBAR
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      SREF = '',E11.4)') SREF
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     SSPAN = '',E11.4)') SSPAN
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      RMPX = '',E11.4)') RMPX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      RMPY = '',E11.4)') RMPY
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      RMPZ = '',E11.4)') RMPZ
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    NORSET = '',I1)') NORSET
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     NVORT = '',I1)') NVORT
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    NPASUM = '',I1)') NPASUM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    JETPAN = '',I1)') JETPAN
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    NBCHGE = '',I1)') NBCHGE
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('OPVSA')
      RETURN
 9999 CALL ERRORS('OPVSA',ERROR)
      CALL EXITS('OPVSA')
      RETURN 1
      END


      SUBROUTINE OPVSB(ERROR,*)

C#### Subroutine: OPVSB
C###  Description:
C###    OPVSB performs patch geometry output for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nonode,nopat,nosect

      CALL ENTERS('OPVSB',*9999)
      WRITE(OP_STRING,'('' NTPAT='',I2)')NTPAT
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      DO 300 nopat = 1,NTPAT
        WRITE(OP_STRING,'('' Patch Number:'',I2)') nopat
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CTX = '',E11.4)') CTX(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CTY = '',E11.4)') CTY(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CTZ = '',E11.4)') CTZ(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      SCAL = '',E11.4)') SCAL(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      THET = '',E11.4)') THET(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CPX = '',E11.4)') CPX(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CPY = '',E11.4)') CPY(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CPZ = '',E11.4)') CPZ(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CHX = '',E11.4)') CHX(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CHY = '',E11.4)') CHY(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       CHZ = '',E11.4)') CHZ(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     IDENT = '',I1)') IDENT(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      MAKE = '',I1)') MAKE(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      KOMP = '',I2)') KOMP(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     KLASS = '',I1)') KLASS(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    NTSECT = '',I2)') NTSECT(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    NTNODE = '',I2)') NTNODE(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       STX = '',E11.4)') STX(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       STY = '',E11.4)') STY(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       STZ = '',E11.4)') STZ(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     SCALE = '',E11.4)') SCALE(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       ALF = '',E11.4)') ALF(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     THETA = '',E11.4)') THETA(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    INMODE = '',I1)') INMODE(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO 200 nosect=1,NTSECT(nopat)
          WRITE(OP_STRING,'(''  SECTION NUMBER='',I2)') nosect
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  NODS = '',20I3)')
     '      (NODS(nonode,nosect,nopat),nonode=1,NTNODE(nopat))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NODES = '',I3)') NODES(nosect,nopat)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 200    CONTINUE
        WRITE(OP_STRING,'(''       NPS = '',I1)') NPS(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     NODEC = '',I1)') NODEC(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''       NPC = '',I1)') NPC(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      INTC = '',I1)') INTC(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      MOVE = '',I1)') MOVE(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 300  CONTINUE

      CALL EXITS('OPVSB')
      RETURN
 9999 CALL ERRORS('OPVSB',ERROR)
      CALL EXITS('OPVSB')
      RETURN 1
      END


      SUBROUTINE OPVSC(ERROR,*)

C#### Subroutine: OPVSC
C###  Description:
C###    OPVSC perrforms wake output for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vsc00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nopat,noseps

      CALL ENTERS('OPVSC',*9999)
      WRITE(OP_STRING,'(''      XMIN = '',E11.4)') XMIN
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      XMAX = '',E11.4)') XMAX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      NODE = '',I1)') NODE
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      NPCW = '',I2)') NPCW
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     INTCW = '',I1)') INTCW
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      MARK = '',I1)') MARK
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      DO 300 nopat=1,NTPAT
        WRITE(OP_STRING,'('' Patch Number:'',I2)') nopat
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    IDENTW = '',I1)') IDENTW(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    IFLEXW = '',I1)') IFLEXW(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     IDEFW = '',I1)') IDEFW(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    KWPACH = '',I2)') KWPACH(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    KWSIDE = '',I1)') KWSIDE(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    KWLINE = '',I2)') KWLINE(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    KWPAN1 = '',I1)') KWPAN1(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
       	WRITE(OP_STRING,'(''    KWPAN2 = '',I1)') KWPAN2(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     INPUT = '',I1)') INPUT(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    NODEWS = '',I1)') NODEWS(1,nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'(''    NTSEPS = '',I1)') NTSEPS(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO 200 noseps=1,NTSEPS(NTPAT)
          WRITE(OP_STRING,'(''       SWPX = '',E11.3)')
     '	    SWPX(noseps,NTPAT)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''       SWPY = '',E11.3)')
     '	    SWPY(noseps,NTPAT)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''     DELTAZ = '',I1)')
     '	    DELTAZ(noseps,NTPAT)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 200    CONTINUE
        WRITE(OP_STRING,'(''    NODEWC = '',I1)') NODEWC(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      NPCP = '',I1)') NPCP(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''     INTCP = '',I1)') INTCP(nopat)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 300  CONTINUE

      CALL EXITS('OPVSC')
      RETURN
 9999 CALL ERRORS('OPVSC',ERROR)
      CALL EXITS('OPVSC')
      RETURN 1
      END


      SUBROUTINE OPVSD(ERROR,*)

C#### Subroutine: OPVSD
C###  Description:
C###    OPVSD performs surface streamline output for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vsd00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER j

      CALL ENTERS('OPVSD',*9999)

      DO 100 j=1,NSTRM
        WRITE(OP_STRING,'('' Panel Number:'',I2)') j
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''        KP = '',I2)') KP(j)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''         F = '',I2)') F(j)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
100   CONTINUE
	
      CALL EXITS('OPVSD')
      RETURN
 9999 CALL ERRORS('OPVSD',ERROR)
      CALL EXITS('OPVSD')
      RETURN 1
      END


      SUBROUTINE OPVSE(ERROR,*)

C#### Subroutine: OPVSE
C###  Description:
C###    OPVSE performs boundary layer output for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vse00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('OPVSE',*9999)
      WRITE(OP_STRING,'(''       RND = '',F10.0)') RND
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    TRIPUP = '',F10.0)') TRIPUP
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    TRIPOP = '',F10.0)') TRIPUP
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''    XPRINT = '',F10.0)') XPRINT
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''     XSKIP = '',F10.0)') XSKIP
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
				
      CALL EXITS('OPVSE')
      RETURN
 9999 CALL ERRORS('OPVSE',ERROR)
      CALL EXITS('OPVSE')
      RETURN 1
      END


      SUBROUTINE OPVSF(ERROR,*)

C#### Subroutine: OPVSF
C###  Description:
C###    OPVSF performs off-body velocity scan output for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vsf00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j

      CALL ENTERS('OPVSF',*9999)

      DO i=1,NBOX
      	WRITE(OP_STRING,'('' Scan box number:'',I2)') i
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      MOLD = '',I1)') MOLD(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      MEET = '',I1)') MEET(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''      NEAR = '',I2)') NEAR(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    INCPRI = '',I2)') INCPRI(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    INCPRJ = '',I2)') INCPRJ(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    INCPRK = '',I2)') INCPRK(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Coords of first corner:'',3F10.3)')
     '    X0(i),Y0(i),Z0(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Box type np = '',I2)') np(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        IF(np(i).GT.1) THEN
          WRITE(OP_STRING,'('' Coords of second corner:'',3F10.3)')
     '      X1(i),Y1(i),Z1(i)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''       NP1 = '',I3)') NP1(i)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(NP1(i).LT.0)THEN
            WRITE(OP_STRING,'('' Point locations AL1'',8F10.3)')
     '        (AL1(j,i),j=1,10)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(np(i).GT.2)THEN
          WRITE(OP_STRING,'('' Coords of third corner:'',3F10.3)')
     '      X2(i),Y2(i),Z2(i)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''       NP2 = '',I3)') NP2(i)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(NP2(i).LT.0)THEN
            WRITE(OP_STRING,'('' Point locations AL2'',8F10.3)')
     '        (AL2(j,i),j=1,10)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(np(i).EQ.4)THEN
          WRITE(OP_STRING,'('' Coords of fourth corner:'',3F10.3)')
     '      X3(i),Y3(i),Z3(i)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''       NP3 = '',I3)') NP3(i)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(NP3(i).LT.0)THEN
            WRITE(OP_STRING,'('' Point locations AL3'',8F10.3)')
     '        (AL3(j,i),j=1,10)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDDO !i

      CALL EXITS('OPVSF')
      RETURN
 9999 CALL ERRORS('OPVSF',ERROR)
      CALL EXITS('OPVSF')
      RETURN 1
      END


      SUBROUTINE OPVSG(ERROR,*)

C#### Subroutine: OPVSG
C###  Description:
C###    OPVSG performs off-body streamline output for VSAERO.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:vsb00.cmn'
      INCLUDE 'cmiss$reference:vsg00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i

      CALL ENTERS('OPVSG',*9999)
      WRITE(OP_STRING,
     '   '(''  Starting point for off-body streamlines'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      DO i=1,NPOINT
      	WRITE(OP_STRING,'('' Streamline number:'',I2)') i
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'('' Coord of starting point:'',3F10.0)')
     '    RSX(i),RSY(i),RSZ(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'('' Distance downstream from start SU = '','
     '	  //'F10.0)') SU(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'('' Distance upstream from start SD = '','
     '	  //'F10.0)') SD(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,'('' Length increment DELS = '',F10.0)')
     '    DELS(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	WRITE(OP_STRING,
     '    '('' Control nearfield routine NEAR2 = '',I2)') NEAR2(i)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDDO
	
      CALL EXITS('OPVSG')
      RETURN
 9999 CALL ERRORS('OPVSG',ERROR)
      CALL EXITS('OPVSG')
      RETURN 1
      END


Module FE90
=========== 


      SUBROUTINE HERM(NBH,NCNP,NEELEM,NHP,NJP,NKH,NKJ,NPE,
     '                NPNODE,NPTNEW,NR,NW,
     '                NYNP,SE,XP,ERROR,*)

C**** Generates extra nodes for hermite BEM problems.
C**** March 1992          
C**** If hermite interpolation is used then an equation will be
C**** generated at each equivalent Lagrange cubic node.  The system of
C**** equations can then be solved in least squares sense or reduced
C**** back to a smaller system of equations by using continuity
C**** constraints.  Need to generate equivalent Lagrange cubic nodal
C**** positions.
C**** 15-5-92
C**** Can in fact generate a far smaller set of equations.  In 2d only
C**** need to have a collocation point in the middle of each element
C**** (as long as boundary is simply connected).  In 3d, depending on
C**** the connectivity of the surface it may be possible to get away
C**** with extra collocation points only at the middle of each side
C**** and possibly one in the middle of each element.
C**** When the problem is coupled to a FE region then no extra
C**** collocation points are needed on any elements which are shared
C**** between the FE and BE regions.
C**** 15-5-92 NONE OF THE ABOVE HAS YET BEEN IMPLEMENTED
C**** 17-5-92 The 2d case consistent with the above comment is being
C**** implemented.
C**** 25-6-92 The minimum number of equations does not seem to be
C**** sufficient (they are dependent). This is being investigated 
C**** further.
C**** ********* 
C**** 10-12-92 
C**** The hypersingular integral equation approach is being used to
C**** generate the extra equations.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,*),NCNP(NHM,0:NCM,NPM,*),
     '        NEELEM(0:NEM,0:*),
     '        NHP(*),NJP(*),NKH(NHM,NCM,*),NKJ(NJM,*),NPE(NNM,NBM,*),
     '        NPNODE(0:NPM,0:*),NPTNEW,NR,NW(*),
     '        NYNP(NKM,NHM,0:NCM,NPM,0:*)
      REAL*8 SE(NSM,NBM,*),XP(NKM,NJM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NB,NE,NE2,NH,NJ,NK,NN1,NN2,NOELEM,NOELEM2,
     '  NONODE,NONODE1,NONODE2,NP1,NP2,
     '  NP3,NP4,NPNEW,nx,NYC(4)
      LOGICAL DEBUG,EXISTS,HYP,NEW,SIMPLE,TESTADAP

      COMMON /TEST/HYP,NEW,DEBUG,SIMPLE,TESTADAP

      CALL ENTERS('HERM',*9999)
      nx=1 !Temporary

      IF(.NOT.HYP)THEN !Old collocation loop
        NPTNEW=0
        CALL ASSERT(NBT+2.LE.NBM,'  >>NBM too small',ERROR,*9999)
        DO NOELEM=1,NEELEM(0,NR) 
          NE=NEELEM(NOELEM,NR)
          IF(NW(NE).EQ.5)THEN !1D cubic hermite dependent var
        			    !interpolation
            IF(NEW)THEN !Use new code (i.e. one extra node per element)
              WRITE(OP_STRING,*)'  Using new code - generating one ',
     '          'extra eqtn'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              NNT(NBT+1)=3 !Generate 1 new node per element
              NKT(0,NBT+1)=1
              NB=NBH(1,1,NE)
              NP1=NPE(1,NB,NE)
              NP2=NPE(2,NB,NE)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
              	  NKH(NH,2,NPNEW)=1
                ELSE
                  !If Hermite interpolation isn't used for the norm der
                  !then nothing gets added below nc=2.
              	  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for xi=1/2
                XP(1,NJ,NPNEW)= 0.5D0*XP(1,NJ,NP1)+
     '      		      1.0D0/8.0D0*XP(2,NJ,NP1)*
     '                        SE(2,NB,NE)+0.5D0*XP(1,NJ,NP2)-
     '      	              1.0D0/8.0D0*XP(2,NJ,NP2)*
     '                                 SE(4,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(1,NBT+1,NE)=NPE(1,NB,NE)
              NPE(2,NBT+1,NE)=NPNEW
              NPE(3,NBT+1,NE)=NPE(2,NB,NE)
            ELSE !Use old code which generates 2 new nodes per element
              NNT(NBT+1)=4 
              NKT(0,NBT+1)=1
              NB=NBH(1,1,NE)
              NP1=NPE(1,NB,NE)
              NP2=NPE(2,NB,NE)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
              	  NKH(NH,2,NPNEW)=1
                ELSE
                  !If Hermite interpolation isn't used for the norm der
                  !then nothing gets added below nc=2.
              	  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for xi=1/3
                XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP1)+
     '      		      4.0D0/27.0D0*XP(2,NJ,NP1)*
     '                                 SE(2,NB,NE)+7.0D0/27.0D0*
     '                                 XP(1,NJ,NP2)-2.0D0/27.0D0*
     '                                 XP(2,NJ,NP2)*SE(4,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
              	  NKH(NH,2,NPNEW)=1
                ELSE
                  !If Hermite interpolation isn't used for the normal 
                  !der then nothing gets added below nc=2.
              	  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for xi=2/3
                XP(1,NJ,NPNEW)=  7.0D0/27.0D0*XP(1,NJ,NP1)+
     '  			       2.0D0/27.0D0*XP(2,NJ,NP1)*
     '                                 SE(2,NB,NE)+20.0D0/27.0D0*
     '                                 XP(1,NJ,NP2)-4.0D0/27.0D0*
     '                                 XP(2,NJ,NP2)*SE(4,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(1,NBT+1,NE)=NPE(1,NB,NE)
              NPE(2,NBT+1,NE)=NPNEW-1
              NPE(3,NBT+1,NE)=NPNEW
              NPE(4,NBT+1,NE)=NPE(2,NB,NE)
            ENDIF !End of new or old code choice
          ELSEIF(NW(NE).EQ.9)THEN !2D bicubic hermite
            !Generate up to 12 new nodes per elem
            NNT(NBT+1)=16
            NKT(0,NBT+1)=1
            NB=NBH(1,1,NE)
            NP1=NPE(1,NB,NE)
            NP2=NPE(2,NB,NE)
            NP3=NPE(3,NB,NE)
            NP4=NPE(4,NB,NE)
            !Check whether the extra nodes along the line between NP1
        	  !and NP2 have already been created
            EXISTS=.FALSE.
            DO NOELEM2=1,NOELEM-1    !Check the elements already covered
              NE2=NEELEM(NOELEM2,NR) !to see if both NP1 and NP2 
              DO NN1=1,NNT(NB)       !belong to one of these
                IF(NP1.EQ.NPE(NN1,NB,NE2))THEN !NP1 in element NE2
                  DO NN2=1,NNT(NB)
                    IF(NP2.EQ.NPE(NN2,NB,NE2))THEN !NP2 also in elem NE2
                      EXISTS=.TRUE.
                      GOTO 1000
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
 1000       CONTINUE
            IF (EXISTS)THEN !extra nodes already created in element NE2
              !Need to check if NE2 is a bicubic hermite element or a
              !linear-cubic hermite element (affects nodal numbering)
              NPE(1,NBT+1,NE)=NP1 
              NPE(4,NBT+1,NE)=NP2
              IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(2,NBT+1,NE)=NPE(2,NBT+1,NE2)
                    NPE(3,NBT+1,NE)=NPE(3,NBT+1,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(2,NBT+1,NE)=NPE(5,NBT+1,NE2)
                    NPE(3,NBT+1,NE)=NPE(9,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP2 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+1,NE)=NPE(8,NBT+1,NE2)
                    NPE(3,NBT+1,NE)=NPE(12,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=2, NN2=1 ',
     '                'OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+1,NE)=NPE(14,NBT+1,NE2)
                    NPE(3,NBT+1,NE)=NPE(15,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=3, NN2=1 ',
     '                'OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
               !NE2 is a cubic hermite-linear element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(2,NBT+1,NE)=NPE(2,NBT+2,NE2)
                    NPE(3,NBT+1,NE)=NPE(3,NBT+2,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(2,NBT+1,NE)=NPE(3,NBT+2,NE2)
                    NPE(3,NBT+1,NE)=NPE(5,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP2 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+1,NE)=NPE(4,NBT+2,NE2)
                    NPE(3,NBT+1,NE)=NPE(6,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=2, NN2=1 ',
     '                'OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+1,NE)=NPE(6,NBT+2,NE2)
                    NPE(3,NBT+1,NE)=NPE(7,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=3, NN2=1 ',
     '                'OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE !Create extra nodes at (1/3,0) and (2/3,0)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (1/3,0)
                XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP1)+
     '                           7.0D0/27.0D0*XP(1,NJ,NP2)
     '                          +4.0D0/27.0D0*XP(2,NJ,NP1)*SE(2,NB,NE)
     '                          -2.0D0/27.0D0*XP(2,NJ,NP2)*SE(6,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (2/3,0)
                XP(1,NJ,NPNEW)= 7.0D0/27.0D0*XP(1,NJ,NP1)+
     '                         20.0D0/27.0D0*XP(1,NJ,NP2)
     '                         +2.0D0/27.0D0*XP(2,NJ,NP1)*SE(2,NB,NE)
     '                         -4.0D0/27.0D0*XP(2,NJ,NP2)*SE(6,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(1,NBT+1,NE)=NP1
              NPE(2,NBT+1,NE)=NPNEW-1
              NPE(3,NBT+1,NE)=NPNEW
              NPE(4,NBT+1,NE)=NP2
            ENDIF !End of NP1-NP2 side of element
            !Check whether the extra nodes along the line between NP1 
            !and NP3 have already been created
            EXISTS=.FALSE.
            DO NOELEM2=1,NOELEM-1    !Check the elements already covered
              NE2=NEELEM(NOELEM2,NR) !to see if both NP1 and NP3 
              DO NN1=1,NNT(NB)       !belong to one of these
                IF(NP1.EQ.NPE(NN1,NB,NE2))THEN !NP1 in element NE2
                  DO NN2=1,NNT(NB)
                    IF(NP3.EQ.NPE(NN2,NB,NE2))THEN !NP3 also in elem NE2
                      EXISTS=.TRUE.
                      GOTO 1100
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
 1100       CONTINUE
            IF (EXISTS)THEN  !extra nodes already created in elem NE2
              !Need to check if NE2 is a bicubic hermite element or a
              !linear-cubic hermite element (affects nodal numbering)
              NPE(1,NBT+1,NE)=NP1
              NPE(13,NBT+1,NE)=NP3
              IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(5,NBT+1,NE)=NPE(2,NBT+1,NE2)
                    NPE(9,NBT+1,NE)=NPE(3,NBT+1,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(5,NBT+1,NE)=NPE(5,NBT+1,NE2)
                    NPE(9,NBT+1,NE)=NPE(9,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP3 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(5,NBT+1,NE)=NPE(8,NBT+1,NE2)
                    NPE(9,NBT+1,NE)=NPE(12,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)' NP1-NP3 SIDE: NN1=2, NN2=1 ',
     '                'OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(5,NBT+1,NE)=NPE(14,NBT+1,NE2)
                    NPE(9,NBT+1,NE)=NPE(15,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)' NP1-NP3 SIDE: NN1=3, NN2=1 ',
     '                'OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP1-NP3 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
               !NE2 is a cubic hermite-linear element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(5,NBT+1,NE)=NPE(2,NBT+2,NE2)
                    NPE(9,NBT+1,NE)=NPE(3,NBT+2,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(5,NBT+1,NE)=NPE(3,NBT+2,NE2)
                    NPE(9,NBT+1,NE)=NPE(5,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP3 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(5,NBT+1,NE)=NPE(4,NBT+2,NE2)
                    NPE(9,NBT+1,NE)=NPE(6,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)' NP1-NP3 SIDE: NN1=2, NN2=1 ',
     '                'OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(5,NBT+1,NE)=NPE(6,NBT+2,NE2)
                    NPE(9,NBT+1,NE)=NPE(7,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP3 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP1-NP3 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE !Create extra nodes at (0,1/3) and (0,2/3)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (0,1/3)
                XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP1)+
     '                           7.0D0/27.0D0*XP(1,NJ,NP3)
     '                          +4.0D0/27.0D0*XP(3,NJ,NP1)*SE(3,NB,NE)
     '                          -2.0D0/27.0D0*XP(3,NJ,NP3)*SE(11,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (0,2/3)
                XP(1,NJ,NPNEW)= 7.0D0/27.0D0*XP(1,NJ,NP1)+
     '                         20.0D0/27.0D0*XP(1,NJ,NP3)
     '                         +2.0D0/27.0D0*XP(3,NJ,NP1)*SE(3,NB,NE)
     '                         -4.0D0/27.0D0*XP(3,NJ,NP3)*SE(11,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(1,NBT+1,NE)=NP1
              NPE(5,NBT+1,NE)=NPNEW-1
              NPE(9,NBT+1,NE)=NPNEW
              NPE(13,NBT+1,NE)=NP3
            ENDIF !End of NP1-NP3 side of element
            !Check whether the extra nodes along the line between NP2 
            !and NP4 have already been created
            EXISTS=.FALSE.
            DO NOELEM2=1,NOELEM-1    !Check the elems already covered 
              NE2=NEELEM(NOELEM2,NR) !to see if both NP2 and NP4 
              DO NN1=1,NNT(NB)       !belong to one of these
                IF(NP2.EQ.NPE(NN1,NB,NE2))THEN !NP1 in element NE2
                  DO NN2=1,NNT(NB)
                    IF(NP4.EQ.NPE(NN2,NB,NE2))THEN !NP4 also in elem NE2
                      EXISTS=.TRUE.
                      GOTO 1200
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
 1200       CONTINUE
            IF (EXISTS)THEN !extra nodes already created in element NE2
              !Need to check if NE2 is a bicubic hermite element or a
              !linear-cubic hermite element (affects nodal numbering)
              NPE(4,NBT+1,NE)=NP2
              NPE(16,NBT+1,NE)=NP4
              IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(8,NBT+1,NE)=NPE(2,NBT+1,NE2)
                    NPE(12,NBT+1,NE)=NPE(3,NBT+1,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(8,NBT+1,NE)=NPE(5,NBT+1,NE2)
                    NPE(12,NBT+1,NE)=NPE(9,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP2-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(8,NBT+1,NE)=NPE(8,NBT+1,NE2)
                    NPE(12,NBT+1,NE)=NPE(12,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP2-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(8,NBT+1,NE)=NPE(14,NBT+1,NE2)
                    NPE(12,NBT+1,NE)=NPE(15,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP2-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP2-NP4 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
               !NE2 is a cubic hermite-linear element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(8,NBT+1,NE)=NPE(2,NBT+2,NE2)
                    NPE(12,NBT+1,NE)=NPE(3,NBT+2,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(8,NBT+1,NE)=NPE(3,NBT+2,NE2)
                    NPE(12,NBT+1,NE)=NPE(5,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP2-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(8,NBT+1,NE)=NPE(4,NBT+2,NE2)
                    NPE(12,NBT+1,NE)=NPE(6,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP2-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(8,NBT+1,NE)=NPE(6,NBT+2,NE2)
                    NPE(12,NBT+1,NE)=NPE(7,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP2-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP2-NP4 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE !Create extra nodes at (1,1/3) and (1,2/3)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (1,1/3)
                XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP2)+
     '                           7.0D0/27.0D0*XP(1,NJ,NP4)
     '                          +4.0D0/27.0D0*XP(3,NJ,NP2)*SE(7,NB,NE)
     '                          -2.0D0/27.0D0*XP(3,NJ,NP4)*SE(15,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (1,2/3)
                XP(1,NJ,NPNEW)= 7.0D0/27.0D0*XP(1,NJ,NP2)+
     '                         20.0D0/27.0D0*XP(1,NJ,NP4)
     '                         +2.0D0/27.0D0*XP(3,NJ,NP2)*SE(7,NB,NE)
     '                         -4.0D0/27.0D0*XP(3,NJ,NP4)*SE(15,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(4,NBT+1,NE)=NP2
              NPE(8,NBT+1,NE)=NPNEW-1
              NPE(12,NBT+1,NE)=NPNEW
              NPE(16,NBT+1,NE)=NP4
            ENDIF !End of NP2-NP4 side of element
            !Check whether the extra nodes along the line between NP3 
            !and NP4 have already been created
            EXISTS=.FALSE.
            DO NOELEM2=1,NOELEM-1    !Check the elems already covered 
              NE2=NEELEM(NOELEM2,NR) !to see if both NP2 and NP4 belong 
              DO NN1=1,NNT(NB)       !to one of these
                IF(NP3.EQ.NPE(NN1,NB,NE2))THEN !NP3 in element NE2
                  DO NN2=1,NNT(NB)
                    IF(NP4.EQ.NPE(NN2,NB,NE2))THEN !NP4 also in elem NE2
                      EXISTS=.TRUE.
                      GOTO 1300
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
 1300       CONTINUE
            IF (EXISTS)THEN !extra nodes already created in element NE2
              !Need to check if NE2 is a bicubic hermite element or a
              !linear-cubic hermite element (affects nodal numbering)
              NPE(13,NBT+1,NE)=NP3
              NPE(16,NBT+1,NE)=NP4
              IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(14,NBT+1,NE)=NPE(2,NBT+1,NE2)
                    NPE(15,NBT+1,NE)=NPE(3,NBT+1,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(14,NBT+1,NE)=NPE(5,NBT+1,NE2)
                    NPE(15,NBT+1,NE)=NPE(9,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(14,NBT+1,NE)=NPE(8,NBT+1,NE2)
                    NPE(15,NBT+1,NE)=NPE(12,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(14,NBT+1,NE)=NPE(14,NBT+1,NE2)
                    NPE(15,NBT+1,NE)=NPE(15,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP3-NP4 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
               !NE2 is a cubic hermite-linear element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(14,NBT+1,NE)=NPE(2,NBT+2,NE2)
                    NPE(15,NBT+1,NE)=NPE(3,NBT+2,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(14,NBT+1,NE)=NPE(3,NBT+2,NE2)
                    NPE(15,NBT+1,NE)=NPE(5,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(14,NBT+1,NE)=NPE(4,NBT+2,NE2)
                    NPE(15,NBT+1,NE)=NPE(6,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(14,NBT+1,NE)=NPE(6,NBT+2,NE2)
                    NPE(15,NBT+1,NE)=NPE(7,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP3-NP4 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE !Create extra nodes at (1/3,1) and (2/3,1)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
        
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite interp isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (1/3,1)
                XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP3)+
     '                           7.0D0/27.0D0*XP(1,NJ,NP4)
     '                          +4.0D0/27.0D0*XP(2,NJ,NP2)*SE(10,NB,NE)
     '                          -2.0D0/27.0D0*XP(2,NJ,NP4)*SE(14,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite interp isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (2/3,1)
                XP(1,NJ,NPNEW)= 7.0D0/27.0D0*XP(1,NJ,NP3)+
     '                         20.0D0/27.0D0*XP(1,NJ,NP4)
     '                         +2.0D0/27.0D0*XP(2,NJ,NP2)*SE(10,NB,NE)
     '                         -4.0D0/27.0D0*XP(2,NJ,NP4)*SE(14,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(13,NBT+1,NE)=NP3
              NPE(14,NBT+1,NE)=NPNEW-1
              NPE(15,NBT+1,NE)=NPNEW
              NPE(16,NBT+1,NE)=NP4
            ENDIF !End of NP3-NP4 side of element
        
            !Create the 4 extra middle nodes
            NPTNEW=NPTNEW+1
            NPNEW=NPNODE(0,NR)+NPTNEW
            NHP(NPNEW)=NHP(NP1)
            NJP(NPNEW)=NJP(NP1)
            DO NJ=1,NJT
              NKJ(NJ,NPNEW)=1
            ENDDO
            DO NH=1,NHP(NPNEW)
              NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
              NKH(NH,1,NPNEW)=1
              IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                NKH(NH,2,NPNEW)=1
              ELSE
               !If Hermite interp isn't used for the normal derivative
               !then nothing gets added below nc=2.
                NKH(NH,2,NPNEW)=0
              ENDIF
            ENDDO
            DO NJ=1,NJT !Node for (1/3,1/3)
              XP(1,NJ,NPNEW)= 20.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP1)
     '                       + 7.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP2)
     '                       +20.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP3)
     '                       + 7.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP4)
     '                       + 4.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP1)*
     '                         SE(2,NB,NE)
     '                       - 2.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP2)*
     '                         SE(6,NB,NE)
     '                       + 4.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP3)*
     '                         SE(10,NB,NE)
     '                       - 2.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP4)*
     '                         SE(14,NB,NE)
     '                       +20.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP1)*
     '                         SE(3,NB,NE)
     '                       + 7.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP2)*
     '                         SE(7,NB,NE)
     '                       -20.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP3)*
     '                         SE(11,NB,NE)
     '                       - 7.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP4)*
     '                         SE(15,NB,NE)
     '                       + 4.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP1)*
     '                         SE(4,NB,NE)
     '                       - 2.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP2)*
     '                         SE(8,NB,NE)
     '                       - 4.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP3)*
     '                         SE(12,NB,NE)
     '                       + 2.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP4)*
     '                         SE(16,NB,NE)
            ENDDO
            NPNODE(NPNEW,NR)=NPNEW
            NPTNEW=NPTNEW+1
            NPNEW=NPNODE(0,NR)+NPTNEW
            NHP(NPNEW)=NHP(NP1)
            NJP(NPNEW)=NJP(NP1)
            DO NJ=1,NJT
              NKJ(NJ,NPNEW)=1
            ENDDO
            DO NH=1,NHP(NPNEW)
              NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
              NKH(NH,1,NPNEW)=1
              IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                NKH(NH,2,NPNEW)=1
              ELSE
               !If Hermite inter isn't used for the normal derivative
               !then nothing gets added below nc=2.
                NKH(NH,2,NPNEW)=0
              ENDIF
            ENDDO
            DO NJ=1,NJT !Node for (1/3,2/3)
              XP(1,NJ,NPNEW)= 20.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP1)
     '                       + 7.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP2)
     '                       +20.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP3)
     '                       + 7.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP4)
     '                       + 4.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP1)*
     '                         SE(2,NB,NE)
     '                       - 2.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP2)*
     '                         SE(6,NB,NE)
     '                       + 4.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP3)*
     '                         SE(10,NB,NE)
     '                       -2.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP4)*
     '                         SE(14,NB,NE)
     '                       +20.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP1)*
     '                         SE(3,NB,NE)
     '                       + 7.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP2)*
     '                         SE(7,NB,NE)
     '                       -20.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP3)*
     '                         SE(11,NB,NE)
     '                       - 7.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP4)*
     '                         SE(15,NB,NE)
     '                       + 4.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP1)*
     '                         SE(4,NB,NE)
     '                       - 2.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP2)*
     '                         SE(8,NB,NE)
     '                       - 4.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP3)*
     '                         SE(12,NB,NE)
     '                       + 2.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP4)*
     '                         SE(16,NB,NE)
            ENDDO
            NPNODE(NPNEW,NR)=NPNEW
            NPTNEW=NPTNEW+1
            NPNEW=NPNODE(0,NR)+NPTNEW
            NHP(NPNEW)=NHP(NP1)
            NJP(NPNEW)=NJP(NP1)
            DO NJ=1,NJT
              NKJ(NJ,NPNEW)=1
            ENDDO
            DO NH=1,NHP(NPNEW)
              NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
              NKH(NH,1,NPNEW)=1
              IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                NKH(NH,2,NPNEW)=1
              ELSE
               !If Hermite inter isn't used for the normal derivative
               !then nothing gets added below nc=2.
                NKH(NH,2,NPNEW)=0
              ENDIF
            ENDDO
            DO NJ=1,NJT !Node for (2/3,1/3)
              XP(1,NJ,NPNEW)=  7.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP1)
     '                       +20.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP2)
     '                       + 7.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP3)
     '                       +20.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP4)
     '                       + 2.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP1)*
     '                         SE(2,NB,NE)
     '                       - 4.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP2)*
     '                         SE(6,NB,NE)
     '                       + 2.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP3)*
     '                         SE(10,NB,NE)
     '                       - 4.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP4)*
     '                         SE(14,NB,NE)
     '                       + 7.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP1)*
     '                         SE(3,NB,NE)
     '                       +20.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP2)*
     '                         SE(7,NB,NE)
     '                       - 7.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP3)*
     '                         SE(11,NB,NE)
     '                       -20.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP4)*
     '                         SE(15,NB,NE)
     '                       + 2.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP1)*
     '                         SE(4,NB,NE)
     '                       - 4.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP2)*
     '                         SE(8,NB,NE)
     '                       - 2.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP3)*
     '                         SE(12,NB,NE)
     '                       + 4.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP4)*
     '                         SE(16,NB,NE)
            ENDDO
            NPNODE(NPNEW,NR)=NPNEW
            NPTNEW=NPTNEW+1
            NPNEW=NPNODE(0,NR)+NPTNEW
            NHP(NPNEW)=NHP(NP1)
            NJP(NPNEW)=NJP(NP1)
            DO NJ=1,NJT
              NKJ(NJ,NPNEW)=1
            ENDDO
            DO NH=1,NHP(NPNEW)
              NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
              NKH(NH,1,NPNEW)=1
              IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                NKH(NH,2,NPNEW)=1
              ELSE
               !If Hermite inter isn't used for the normal derivative
               !then nothing gets added below nc=2.
                NKH(NH,2,NPNEW)=0
              ENDIF
            ENDDO
            DO NJ=1,NJT !Node for (2/3,2/3)
              XP(1,NJ,NPNEW)=  7.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP1)
     '                       +20.0D0/27.0D0*7.0D0/27.0D0*XP(1,NJ,NP2)
     '                       + 7.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP3)
     '                       +20.0D0/27.0D0*20.0D0/27.0D0*XP(1,NJ,NP4)
     '                       + 2.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP1)*
     '                         SE(2,NB,NE)
     '                       - 4.0D0/27.0D0*7.0D0/27.0D0*XP(2,NJ,NP2)*
     '                         SE(6,NB,NE)
     '                       + 2.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP3)*
     '                         SE(10,NB,NE)
     '                       - 4.0D0/27.0D0*20.0D0/27.0D0*XP(2,NJ,NP4)*
     '                         SE(14,NB,NE)
     '                       + 7.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP1)*
     '                         SE(3,NB,NE)
     '                       +20.0D0/27.0D0*2.0D0/27.0D0*XP(3,NJ,NP2)*
     '                         SE(7,NB,NE)
     '                       - 7.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP3)*
     '                         SE(11,NB,NE)
     '                       -20.0D0/27.0D0*4.0D0/27.0D0*XP(3,NJ,NP4)*
     '                         SE(15,NB,NE)
     '                       + 2.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP1)*
     '                         SE(4,NB,NE)
     '                       - 4.0D0/27.0D0*2.0D0/27.0D0*XP(4,NJ,NP2)*
     '                         SE(8,NB,NE)
     '                       - 2.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP3)*
     '                         SE(12,NB,NE)
     '                       + 4.0D0/27.0D0*4.0D0/27.0D0*XP(4,NJ,NP4)*
     '                         SE(16,NB,NE)
            ENDDO
            NPNODE(NPNEW,NR)=NPNEW
        
            NPE(6,NBT+1,NE)=NPNEW-3
            NPE(7,NBT+1,NE)=NPNEW-1
            NPE(10,NBT+1,NE)=NPNEW-2
            NPE(11,NBT+1,NE)=NPNEW
        
            IF(DOP)THEN
              WRITE(OP_STRING,'('' npe(nn,'',I2,'','',I2,'')='')')
     '          NBT+1,NE    
        	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(16(I3,2X))')
     '          (NPE(NN1,NBT+1,NE),NN1=1,NNT(NBT+1))
        	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO NN1=1,NNT(NBT+1)
                WRITE(OP_STRING,'('' xp(1,nj,'',I3,'')='')')
     '            NPE(NN1,NBT+1,NE)
        	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' '',3(E10.4,3X))')
     '            (XP(1,NJ,NPE(NN1,NBT+1,NE)),NJ=1,NJT)
        	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,*)
        	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
        
          ELSEIF(NW(NE).EQ.11)THEN !2d linear-cubic hermite element
        			      !Only need to generate extra nodes
        			      !in the cubic direction
            NNT(NBT+2)=8 !Need to use NBT+2 since NBT+1 is allocated to
            NKT(0,NBT+2)=1 !2d bicubic lagrange basis function 
            !Generate up to 4 new nodes per element
            NB=NBH(1,1,NE)
            NP1=NPE(1,NB,NE)
            NP2=NPE(2,NB,NE)
            NP3=NPE(3,NB,NE)
            NP4=NPE(4,NB,NE)
            !If the element is not a collapsed node element then only
            !need to create extra nodes along the NP1-NP3 and
            !NP2-NP4 sides.
            !Check whether the extra nodes along the line between NP1 
            !and NP3 have already been created
            EXISTS=.FALSE.
            IF(NP1.NE.NP3)THEN !Nodes are not equal
        	    DO NOELEM2=1,NOELEM-1    !Check the elements already covered
        	      NE2=NEELEM(NOELEM2,NR) !to see if both NP1 and NP3 
        	      DO NN1=1,NNT(NB)       !belong to one of these
        		IF(NP1.EQ.NPE(NN1,NB,NE2))THEN !NP1 in element NE2
        		  DO NN2=1,NNT(NB)
        		    IF(NP3.EQ.NPE(NN2,NB,NE2))THEN !NP3 also in elem NE2
        		      EXISTS=.TRUE.
        		      GOTO 2100
        		    ENDIF
        		  ENDDO
        		ENDIF
        	      ENDDO
        	    ENDDO
 2100   	    CONTINUE
        	    IF (EXISTS)THEN  !extra nodes already created in elem NE2
        	      !Need to check if NE2 is a bicubic hermite element or a
        	      !linear-cubic hermite element (affects nodal numbering)
        	      NPE(1,NBT+2,NE)=NP1
        	      NPE(7,NBT+2,NE)=NP3
        	      IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
        		IF(NN1.EQ.1)THEN
        		  IF(NN2.EQ.2)THEN
        		    NPE(3,NBT+2,NE)=NPE(2,NBT+1,NE2)
        		    NPE(5,NBT+2,NE)=NPE(3,NBT+1,NE2)
        		  ELSEIF(NN2.EQ.3)THEN
        		    NPE(3,NBT+2,NE)=NPE(5,NBT+1,NE2)
        		    NPE(5,NBT+2,NE)=NPE(9,NBT+1,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP1-NP3 SIDE: NN1=1 AND NN2=4 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.2)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(3,NBT+2,NE)=NPE(8,NBT+1,NE2)
        		    NPE(5,NBT+2,NE)=NPE(12,NBT+1,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP1-NP3 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.3)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(3,NBT+2,NE)=NPE(14,NBT+1,NE2)
        		    NPE(5,NBT+2,NE)=NPE(15,NBT+1,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP1-NP3 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSE
        		  WRITE(OP_STRING,*)
     '                      ' NP1-NP3 SIDE: NN1=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		  GOTO 9999
        		ENDIF
        	      ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
        	       !NE2 is a cubic hermite-linear element
        		IF(NN1.EQ.1)THEN
        		  IF(NN2.EQ.2)THEN
        		    NPE(3,NBT+2,NE)=NPE(2,NBT+2,NE2)
        		    NPE(5,NBT+2,NE)=NPE(3,NBT+2,NE2)
        		  ELSEIF(NN2.EQ.3)THEN
        		    NPE(3,NBT+2,NE)=NPE(3,NBT+2,NE2)
        		    NPE(5,NBT+2,NE)=NPE(5,NBT+2,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP1-NP3 SIDE: NN1=1 AND NN2=4 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.2)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(3,NBT+2,NE)=NPE(4,NBT+2,NE2)
        		    NPE(5,NBT+2,NE)=NPE(6,NBT+2,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP1-NP3 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.3)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(3,NBT+2,NE)=NPE(6,NBT+2,NE2)
        		    NPE(5,NBT+2,NE)=NPE(7,NBT+2,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP1-NP3 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSE
        		  WRITE(OP_STRING,*)
     '                      ' NP1-NP3 SIDE: NN1=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		  GOTO 9999
        		ENDIF
        	      ENDIF
        	    ELSE !Create extra nodes at (0,1/3) and (0,2/3)
        	      NPTNEW=NPTNEW+1
        	      NPNEW=NPNODE(0,NR)+NPTNEW
        	      NHP(NPNEW)=NHP(NP1)
        	      NJP(NPNEW)=NJP(NP1)
        	      DO NJ=1,NJT
        		NKJ(NJ,NPNEW)=1
        	      ENDDO
        	      DO NH=1,NHP(NPNEW)
        		NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
        		NKH(NH,1,NPNEW)=1
        		IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
        		  NKH(NH,2,NPNEW)=1
        		ELSE
        		 !If Hermite inter isn't used for the normal derivative
        		 !then nothing gets added below nc=2.
        		  NKH(NH,2,NPNEW)=0
        		ENDIF
        	      ENDDO
        	      DO NJ=1,NJT !Node for (0,1/3)
        		XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP1)
     '                                  +7.0D0/27.0D0*XP(1,NJ,NP3)
     '  			        +4.0D0/27.0D0*XP(2,NJ,NP1)
     '                                   *SE(2,NB,NE)
     '  				-2.0D0/27.0D0*XP(2,NJ,NP3)
     '                                   *SE(6,NB,NE)
        	      ENDDO
        	      NPNODE(NPNEW,NR)=NPNEW
        	      NPTNEW=NPTNEW+1
        	      NPNEW=NPNODE(0,NR)+NPTNEW
        	      NHP(NPNEW)=NHP(NP1)
        	      NJP(NPNEW)=NJP(NP1)
        	      DO NJ=1,NJT
        		NKJ(NJ,NPNEW)=1
        	      ENDDO
        	      DO NH=1,NHP(NPNEW)
        		NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
        		NKH(NH,1,NPNEW)=1
        		IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
        		  NKH(NH,2,NPNEW)=1
        		ELSE
        		 !If Hermite inter isn't used for the normal derivative
        		 !then nothing gets added below nc=2.
        		  NKH(NH,2,NPNEW)=0
        		ENDIF
        	      ENDDO
        	      DO NJ=1,NJT !Node for (0,2/3)
        		XP(1,NJ,NPNEW)= 7.0D0/27.0D0*XP(1,NJ,NP1)
     '                                 +20.0D0/27.0D0*XP(1,NJ,NP3)
     '  			        +2.0D0/27.0D0*XP(2,NJ,NP1)
     '                                   *SE(2,NB,NE)
     '  				-4.0D0/27.0D0*XP(2,NJ,NP3)
     '                                   *SE(6,NB,NE)
        	      ENDDO
        	      NPNODE(NPNEW,NR)=NPNEW
        	      NPE(1,NBT+2,NE)=NP1
        	      NPE(3,NBT+2,NE)=NPNEW-1
        	      NPE(5,NBT+2,NE)=NPNEW
        	      NPE(7,NBT+2,NE)=NP3
        	    ENDIF !End of NP1-NP3 side of element
            ELSE !Element has some collapsed nodes
        	    NPE(1,NBT+2,NE)=NP1
        	    NPE(3,NBT+2,NE)=NP1
        	    NPE(5,NBT+2,NE)=NP1
        	    NPE(7,NBT+2,NE)=NP1
            ENDIF
            !Check whether the extra nodes along the line between NP2 and
        	  !NP4 have already been created
            EXISTS=.FALSE.
            IF(NP2.NE.NP4)THEN
        	    DO NOELEM2=1,NOELEM-1    !Check the elements already covered 
        	      NE2=NEELEM(NOELEM2,NR) !to see if both NP2 and NP4 
        	      DO NN1=1,NNT(NB)       !belong to one of these
        		IF(NP2.EQ.NPE(NN1,NB,NE2))THEN !NP1 in element NE2
        		  DO NN2=1,NNT(NB)
        		    IF(NP4.EQ.NPE(NN2,NB,NE2))THEN !NP4 also in elem NE2
        		      EXISTS=.TRUE.
        		      GOTO 2200
        		    ENDIF
        		  ENDDO
        		ENDIF
        	      ENDDO
        	    ENDDO
 2200   	    CONTINUE
        	    IF (EXISTS)THEN !extra nodes already created in element NE2
        	      !Need to check if NE2 is a bicubic hermite element or a
        	      !linear-cubic hermite element (affects nodal numbering)
        	      NPE(2,NBT+2,NE)=NP2
        	      NPE(8,NBT+2,NE)=NP4
        	      IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
        		IF(NN1.EQ.1)THEN
        		  IF(NN2.EQ.2)THEN
        		    NPE(4,NBT+2,NE)=NPE(2,NBT+1,NE2)
        		    NPE(6,NBT+2,NE)=NPE(3,NBT+1,NE2)
        		  ELSEIF(NN2.EQ.3)THEN
        		    NPE(4,NBT+2,NE)=NPE(5,NBT+1,NE2)
        		    NPE(6,NBT+2,NE)=NPE(9,NBT+1,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP2-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.2)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(4,NBT+2,NE)=NPE(8,NBT+1,NE2)
        		    NPE(6,NBT+2,NE)=NPE(12,NBT+1,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP2-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.3)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(4,NBT+2,NE)=NPE(14,NBT+1,NE2)
        		    NPE(6,NBT+2,NE)=NPE(15,NBT+1,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP2-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSE
        		  WRITE(OP_STRING,*)
     '                      ' NP2-NP4 SIDE: NN1=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		  GOTO 9999
        		ENDIF
        	      ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
        	       !NE2 is a cubic hermite-linear element
        		IF(NN1.EQ.1)THEN
        		  IF(NN2.EQ.2)THEN
        		    NPE(4,NBT+2,NE)=NPE(2,NBT+2,NE2)
        		    NPE(6,NBT+2,NE)=NPE(3,NBT+2,NE2)
        		  ELSEIF(NN2.EQ.3)THEN
        		    NPE(4,NBT+2,NE)=NPE(3,NBT+2,NE2)
        		    NPE(6,NBT+2,NE)=NPE(5,NBT+2,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP2-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.2)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(4,NBT+2,NE)=NPE(4,NBT+2,NE2)
        		    NPE(6,NBT+2,NE)=NPE(6,NBT+2,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP2-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSEIF(NN1.EQ.3)THEN
        		  IF(NN2.EQ.4)THEN
        		    NPE(4,NBT+2,NE)=NPE(6,NBT+2,NE2)
        		    NPE(6,NBT+2,NE)=NPE(7,NBT+2,NE2)
        		  ELSE
        		    WRITE(OP_STRING,*)
     '                        ' NP2-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
        		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		    GOTO 9999
        		  ENDIF
        		ELSE
        		  WRITE(OP_STRING,*)
     '                      ' NP2-NP4 SIDE: NN1=4 - ERROR'
        		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        		  GOTO 9999
        		ENDIF
        	      ENDIF
        	    ELSE !Create extra nodes at (1,1/3) and (1,2/3)
        	      NPTNEW=NPTNEW+1
        	      NPNEW=NPNODE(0,NR)+NPTNEW
        	      NHP(NPNEW)=NHP(NP1)
        	      NJP(NPNEW)=NJP(NP1)
        	      DO NJ=1,NJT
        		NKJ(NJ,NPNEW)=1
        	      ENDDO
        	      DO NH=1,NHP(NPNEW)
        		NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
        		NKH(NH,1,NPNEW)=1
        		IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
        		  NKH(NH,2,NPNEW)=1
        		ELSE
        		 !If Hermite inter isn't used for the normal derivative
        		 !then nothing gets added below nc=2.
        		  NKH(NH,2,NPNEW)=0
        		ENDIF
        	      ENDDO
        	      DO NJ=1,NJT !Node for (1,1/3)
        		XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP2)
     '                                  +7.0D0/27.0D0*XP(1,NJ,NP4)
     '  				+4.0D0/27.0D0*XP(2,NJ,NP2)
     '                                   *SE(4,NB,NE)
     '  				-2.0D0/27.0D0*XP(2,NJ,NP4)
     '                                   *SE(8,NB,NE)
        	      ENDDO
        	      NPNODE(NPNEW,NR)=NPNEW
        	      NPTNEW=NPTNEW+1
        	      NPNEW=NPNODE(0,NR)+NPTNEW
        	      NHP(NPNEW)=NHP(NP1)
        	      NJP(NPNEW)=NJP(NP1)
        	      DO NJ=1,NJT
        		NKJ(NJ,NPNEW)=1
        	      ENDDO
        	      DO NH=1,NHP(NPNEW)
        		NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
        		NKH(NH,1,NPNEW)=1
        		IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
        		  NKH(NH,2,NPNEW)=1
        		ELSE
        		 !If Hermite inter isn't used for the normal derivative
        		 !then nothing gets added below nc=2.
        		  NKH(NH,2,NPNEW)=0
        		ENDIF
        	      ENDDO
        	      DO NJ=1,NJT !Node for (1,2/3)
        		XP(1,NJ,NPNEW)=  7.0D0/27.0D0*XP(1,NJ,NP2)
     '                                 +20.0D0/27.0D0*XP(1,NJ,NP4)
     '  				+2.0D0/27.0D0*XP(2,NJ,NP2)
     '                                   *SE(4,NB,NE)
     '  				-4.0D0/27.0D0*XP(2,NJ,NP4)
     '                                   *SE(8,NB,NE)
        	      ENDDO
        	      NPNODE(NPNEW,NR)=NPNEW
        	      NPE(2,NBT+2,NE)=NP2
        	      NPE(4,NBT+2,NE)=NPNEW-1
        	      NPE(6,NBT+2,NE)=NPNEW
        	      NPE(8,NBT+2,NE)=NP4
        	    ENDIF !End of NP2-NP4 side of element
            ELSE !Element has some collapsed nodes
        	    NPE(2,NBT+2,NE)=NP2
        	    NPE(4,NBT+2,NE)=NP2
        	    NPE(6,NBT+2,NE)=NP2
        	    NPE(8,NBT+2,NE)=NP2
            ENDIF
              
            IF(DOP)THEN
              WRITE(OP_STRING,'('' npe(nn,'',I2,'','',I2,'')='')')
     '          NBT+2,NE    
        	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(16(I3,2X))')
     '          (NPE(NN1,NBT+2,NE),NN1=1,NNT(NBT+2))
        	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO NN1=1,NNT(NBT+2)
                WRITE(OP_STRING,'('' xp(1,nj,'',I3,'')='')')
     '            NPE(NN1,NBT+2,NE)   
        	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' '',3(E10.4,3X))')
     '            (XP(1,NJ,NPE(NN1,NBT+2,NE)),NJ=1,NJT)
        	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        	      WRITE(OP_STRING,*)
        	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
        
          ELSE IF(NW(NE).EQ.13)THEN !2d cubic hermite-linear element
        			      !Only need to generate extra nodes
        			      !in the cubic direction
            NNT(NBT+2)=8 !Need to use NBT+2 since NBT+1 is allocated to
            NKT(0,NBT+2)=1 !2d bicubic lagrange basis function 
            !Generate up to 4 new nodes per element
            NB=NBH(1,1,NE)
            NP1=NPE(1,NB,NE)
            NP2=NPE(2,NB,NE)
            NP3=NPE(3,NB,NE)
            NP4=NPE(4,NB,NE)
            !Only need to create extra nodes along the NP1-NP2 and
        	  !NP3-NP4 sides
            !Check whether the extra nodes along the line between NP1 and
        	  !NP2 have already been created
            EXISTS=.FALSE.
            DO NOELEM2=1,NOELEM-1    !Check the elements already covered
              NE2=NEELEM(NOELEM2,NR) !to see if both NP1 and NP3 
              DO NN1=1,NNT(NB)       !belong to one of these
                IF(NP1.EQ.NPE(NN1,NB,NE2))THEN !NP1 in element NE2
                  DO NN2=1,NNT(NB)
                    IF(NP2.EQ.NPE(NN2,NB,NE2))THEN !NP2 also in elem NE2
                      EXISTS=.TRUE.
                      GOTO 3100
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
 3100       CONTINUE
            IF (EXISTS)THEN  !extra nodes already created in elem NE2
              !Need to check if NE2 is a bicubic hermite element or a
        	    !linear-cubic hermite element (affects nodal numbering)
              NPE(1,NBT+2,NE)=NP1
              NPE(4,NBT+2,NE)=NP2
              IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(2,NBT+2,NE)=NPE(2,NBT+1,NE2)
                    NPE(3,NBT+2,NE)=NPE(3,NBT+1,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(2,NBT+2,NE)=NPE(5,NBT+1,NE2)
                    NPE(3,NBT+2,NE)=NPE(9,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP2 SIDE: NN1=1 AND NN2=4 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+2,NE)=NPE(8,NBT+1,NE2)
                    NPE(3,NBT+2,NE)=NPE(12,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP2 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+2,NE)=NPE(14,NBT+1,NE2)
                    NPE(3,NBT+2,NE)=NPE(15,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP3 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
               !NE2 is a cubic hermite-linear element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(2,NBT+2,NE)=NPE(2,NBT+2,NE2)
                    NPE(3,NBT+2,NE)=NPE(3,NBT+2,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(2,NBT+2,NE)=NPE(3,NBT+2,NE2)
                    NPE(3,NBT+2,NE)=NPE(5,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP2 SIDE: NN1=1 AND NN2=4 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+2,NE)=NPE(4,NBT+2,NE2)
                    NPE(3,NBT+2,NE)=NPE(6,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP2 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(2,NBT+2,NE)=NPE(6,NBT+2,NE2)
                    NPE(3,NBT+2,NE)=NPE(7,NBT+2,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP1-NP2 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP1-NP2 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE !Create extra nodes at (1/3,0) and (2/3,0)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (1/3,0)
                XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP1)+
     '                           7.0D0/27.0D0*XP(1,NJ,NP2)
     '                          +4.0D0/27.0D0*XP(2,NJ,NP1)*SE(2,NB,NE)
     '                          -2.0D0/27.0D0*XP(2,NJ,NP2)*SE(6,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite inter isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (2/3,0)
                XP(1,NJ,NPNEW)= 7.0D0/27.0D0*XP(1,NJ,NP1)+
     '                         20.0D0/27.0D0*XP(1,NJ,NP2)
     '                         +2.0D0/27.0D0*XP(2,NJ,NP1)*SE(2,NB,NE)
     '                         -4.0D0/27.0D0*XP(2,NJ,NP2)*SE(6,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(1,NBT+2,NE)=NP1
              NPE(2,NBT+2,NE)=NPNEW-1
              NPE(3,NBT+2,NE)=NPNEW
              NPE(4,NBT+2,NE)=NP2
            ENDIF !End of NP1-NP2 side of element
        
            !Check whether the extra nodes along the line between NP3 and
        	  !NP4 have already been created
            EXISTS=.FALSE.
            DO NOELEM2=1,NOELEM-1    !Check the elements already covered 
              NE2=NEELEM(NOELEM2,NR) !to see if both NP2 and NP4 
              DO NN1=1,NNT(NB)       !belong to one of these
                IF(NP3.EQ.NPE(NN1,NB,NE2))THEN !NP3 in element NE2
                  DO NN2=1,NNT(NB)
                    IF(NP4.EQ.NPE(NN2,NB,NE2))THEN !NP4 also in elem NE2
                      EXISTS=.TRUE.
                      GOTO 3200
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
 3200       CONTINUE
            IF (EXISTS)THEN !extra nodes already created in element NE2
              !Need to check if NE2 is a bicubic hermite element or a
        	    !linear-cubic hermite element (affects nodal numbering)
              NPE(5,NBT+2,NE)=NP3
              NPE(8,NBT+2,NE)=NP4
              IF(NW(NE2).EQ.9)THEN !NE2 is a bicubic hermite element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(5,NBT+2,NE)=NPE(2,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(3,NBT+1,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(5,NBT+2,NE)=NPE(5,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(9,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(5,NBT+2,NE)=NPE(8,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(12,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999                  
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(5,NBT+2,NE)=NPE(14,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(15,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP3-NP4 SIDE: NN1=4 - ERROR'
        		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ELSEIF(NW(NE2).EQ.11.OR.NW(NE2).EQ.13)THEN
               !NE2 is a cubic hermite-linear element
                IF(NN1.EQ.1)THEN
                  IF(NN2.EQ.2)THEN
                    NPE(4,NBT+2,NE)=NPE(2,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(3,NBT+1,NE2)
                  ELSEIF(NN2.EQ.3)THEN
                    NPE(4,NBT+2,NE)=NPE(3,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(5,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=1 AND NN2=4 - ERROR'
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.2)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(4,NBT+2,NE)=NPE(4,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(6,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=2, NN2=1 OR 3 - ERROR'
                    GOTO 9999
                  ENDIF
                ELSEIF(NN1.EQ.3)THEN
                  IF(NN2.EQ.4)THEN
                    NPE(4,NBT+2,NE)=NPE(6,NBT+1,NE2)
                    NPE(6,NBT+2,NE)=NPE(7,NBT+1,NE2)
                  ELSE
                    WRITE(OP_STRING,*)
     '                ' NP3-NP4 SIDE: NN1=3, NN2=1 OR 2 - ERROR'
                    GOTO 9999
                  ENDIF
                ELSE
                  WRITE(OP_STRING,*)' NP3-NP4 SIDE: NN1=4 - ERROR'
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE !Create extra nodes at (1/3,1) and (2/3,1)
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
        
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite interp isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (1/3,1)
                XP(1,NJ,NPNEW)= 20.0D0/27.0D0*XP(1,NJ,NP3)+
     '                           7.0D0/27.0D0*XP(1,NJ,NP4)
     '                          +4.0D0/27.0D0*XP(2,NJ,NP2)*SE(10,NB,NE)
     '                          -2.0D0/27.0D0*XP(2,NJ,NP4)*SE(14,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPTNEW=NPTNEW+1
              NPNEW=NPNODE(0,NR)+NPTNEW
              NHP(NPNEW)=NHP(NP1)
              NJP(NPNEW)=NJP(NP1)
              DO NJ=1,NJT
                NKJ(NJ,NPNEW)=1
              ENDDO
              DO NH=1,NHP(NPNEW)
                NCNP(NH,0,NPNEW,NR)=NCNP(NH,0,NP1,NR)
                NKH(NH,1,NPNEW)=1
                IF(NBH(NH,1,NE).EQ.NBH(NH,2,NE))THEN
                  NKH(NH,2,NPNEW)=1
                ELSE
                 !If Hermite interp isn't used for the normal derivative
                 !then nothing gets added below nc=2.
                  NKH(NH,2,NPNEW)=0
                ENDIF
              ENDDO
              DO NJ=1,NJT !Node for (2/3,1)
                XP(1,NJ,NPNEW)= 7.0D0/27.0D0*XP(1,NJ,NP3)+
     '                         20.0D0/27.0D0*XP(1,NJ,NP4)
     '                         +2.0D0/27.0D0*XP(2,NJ,NP2)*SE(10,NB,NE)
     '                         -4.0D0/27.0D0*XP(2,NJ,NP4)*SE(14,NB,NE)
              ENDDO
              NPNODE(NPNEW,NR)=NPNEW
              NPE(5,NBT+2,NE)=NP3
              NPE(6,NBT+2,NE)=NPNEW-1
              NPE(7,NBT+2,NE)=NPNEW
              NPE(8,NBT+2,NE)=NP4
            ENDIF !End of NP3-NP4 side of element
          ENDIF
        ENDDO !End of loop over elements
        
C ***   Need to adjust NYNP for the extra nodes NP
C ***   Note that there will be NPNODE(0,NR)+NPTNEW equations generated.
C ***   Need to ensure that the extra node equations contribute to the
C ***   global matrices in the same way as an equation for a derivative
C ***   would have.
        NONODE1=NPNODE(0,NR)+1 !First new node number
        NP1=NPNODE(NONODE1,NR) !New node generated above
        DO NONODE2=1,NPNODE(0,NR) !Loop over derivs for existing nodes
          NP2=NPNODE(NONODE2,NR)
          DO NH=1,NHP(NP2)
            DO NC=1,NCNP(NH,0,NP2,NR)
              DO NK=2,NKH(NH,NC,NP2) !Fill in places for deriv equations
                NYNP(1,NH,NC,NP1,NR)=NYNP(NK,NH,NC,NP2,NR)
                NONODE1=NONODE1+1
                NP1=NPNODE(NONODE1,NR)
              ENDDO !End of NK loop
            ENDDO !End of NC loop
          ENDDO !End of NH loop
        ENDDO !End of loop over existing nodes
        DO NC=1,NCT(nr)
          NYC(NC)=NYT(NC,NR,nx)
        ENDDO
        DO NONODE=NONODE1,NPNODE(0,NR)+NPTNEW !Fill in extra places
          !If Hermite interpolation is not used for the normal derivative
          !then nothing gets added for nc=2.
          NP1=NPNODE(NONODE,NR)
          DO NH=1,NHP(NP1)
            DO NC=1,NCNP(NH,0,NP1,NR)
              DO NK=1,NKH(NH,NC,NP1)
                IF(NC.LE.2)THEN
                  NYC(NC)=NYC(NC)+1
                  NYNP(NK,NH,NC,NP1,NR)=NYC(NC)
                ELSEIF(NC.GT.2)THEN
                  NYC(2)=NYC(2)+1
                  NYNP(NK,NH,NC,NP1,NR)=NYC(2)
                ENDIF
              ENDDO !End NK loop
            ENDDO !End NC loop
          ENDDO !End NH loop
        ENDDO !End of nonode loop (i.e. NP loop)
      ELSE !End of old collocation loop
        NPTNEW=0
      ENDIF 

 9998 CALL EXITS('HERM')
      RETURN
 9999 CALL ERRORS('HERM',ERROR)
      CALL EXITS('HERM')
      RETURN 1
      END



      SUBROUTINE UPBE(NHP,NKH,NYNP,CE,XP,YP,ERROR,*)

C**** Update appropriate coordinates and boundary conditions if used in
C**** conjuction with optimisation routine

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp90.cmn'
      INCLUDE 'cmiss$reference:ktyp100.cmn'
      INCLUDE 'cmiss$reference:opti00.cmn'
!     Parameter List
      INTEGER NHP(*),NKH(NHM,NCM,*),NYNP(NKM,NHM,0:NCM,*)
      REAL*8 CE(*),XP(NKM,NJM,*),YP(NYM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NC,NH,NK,NONODE,NP,NY
      REAL*8 S

      CALL ENTERS('UPBE',*9999)
      NC=1 !Temporary AJP 19-12-91
      IF(KTYP27.EQ.3) THEN !Coupled saturated-unsaturated problem
        IF(KTYP100.EQ.1)THEN !Flow from cavities
C ***     Nodal coordinates have been updated in SOLVE1.
C ***     We need to update the bcs used in the BE region.
          S=CE(1) !Value of s for unsaturated flow
          DO NONODE=1,NPJOIN(0,1)
            NP=NPJOIN(NONODE,1)
            DO NH=1,NHP(NP)
              DO NK=1,MAX(NKH(NH,NC,NP)-KTYP93(NC),1)
                NY=NYNP(NK,NH,NC,NP)
                IF(ITYP10(1).EQ.1)THEN !Cartesians
                  YP(NY,1)=DEXP(-S*XP(1,NJT,NP))
                  !Update essential bcs (no corners)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSEIF(KTYP100.EQ.2)THEN !Flow around cavities
C ***     Nodal coordinates have been updated in SOLVE1.
C ***     We need to update the bcs used in the BE region.
          S=CE(1) !Value of s for unsaturated flow
          DO NONODE=1,NPJOIN(0,1)
            NP=NPJOIN(NONODE,1)
            DO NH=1,NHP(NP)
              DO NK=1,MAX(NKH(NH,NC,NP)-KTYP93(NC),1)
                NY=NYNP(NK,NH,NC,NP)
                IF(ITYP10(1).EQ.1)THEN !Cartesians
                  YP(NY,1)=(THETA_SAT-1.0D0)*DEXP(-S*XP(1,NJT,NP))
                  !Update essential bcs (no corners)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
C ***     Since Theta_sat has been changed need to update the bc at
C ***     any nodes that are on the boundary of both regions
          DO NONODE=1,NPSHARE(0)
            NP=NPSHARE(NONODE)
            DO NH=1,NHP(NP)
              DO NK=1,MAX(NKH(NH,NC,NP)-KTYP93(NC),1)
                NY=NYNP(NK,NH,NC,NP)
                IF(ITYP10(1).EQ.1)THEN !Cartesians
                  YP(NY,1)=(THETA_SAT-1.0D0)*DEXP(-S*XP(1,NJT,NP))
                  !Update essential bcs (no corners)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

CC ***     NNODES = TOTAL NUMBER OF NODES IN THE COUPLED PROBLEM.
C        NNODES=(2*NRAELE+1)*2*NTHELE+2*NTHELE
C        S=CE(1,1) !Change this to allow for the 2 regions.
C        DO NOOPTI=1,NTOPTI
C          NP1=NP1OPT(NOOPTI)
C          NP2=NP2OPT(NOOPTI)
C          XP(1,1,NP1)=PAOPTI(NOOPTI)/PBOPTI(NOOPTI)*XP(1,1,NP1)
C          XP(1,2,NP1)=PAOPTI(NOOPTI)/PBOPTI(NOOPTI)*XP(1,2,NP1)
C          XP(1,1,NP2)=XP(1,1,NP1)
C          XP(1,2,NP2)=XP(1,2,NP1)
C          PBOPTI(NOOPTI)=PAOPTI(NOOPTI)
CC Updating the coord. of the F.E. nodes and the outer circle of 
CC B.E nodes
C          IF(NOOPTI.NE.1.AND.NOOPTI.NE.NTOPTI)THEN
C            XP(1,1,NNODES-(2*NTHELE)-(NOOPTI-2))=-XP(1,1,NP1)
C            XP(1,2,NNODES-(2*NTHELE)-(NOOPTI-2))=XP(1,2,NP1)
C            XP(1,1,NNODES-(NOOPTI-2))=-XP(1,1,NP2)
C            XP(1,2,NNODES-(NOOPTI-2))=XP(1,2,NP2)
C          ENDIF
CC ***     Updating the x coord. of the F.E. nodes
C          DIF=(XP(1,1,NP1)-XP(1,1,NOOPTI))/(2*NRAELE)
C          DO NNRAELE=1,2*NRAELE-1
C            XP(1,1,2*NTHELE*NNRAELE+NOOPTI)=NNRAELE*DIF+XP(1,1,NOOPTI)
C            IF(NOOPTI.NE.1.AND.NOOPTI.NE.NTOPTI)THEN
C              XP(1,1,2*NTHELE*(NNRAELE+1)-(NOOPTI-2))=
C     '          -XP(1,1,2*NTHELE*NNRAELE+NOOPTI)
C            ENDIF
C          ENDDO
CC ***     Updating the y coord. of the F.E. nodes
C          DIF=(XP(1,2,NP1)-XP(1,2,NOOPTI))/(2*NRAELE)
C          DO NNRAELE=1,2*NRAELE-1
C            XP(1,2,2*NTHELE*NNRAELE+NOOPTI)=NNRAELE*DIF+XP(1,2,NOOPTI)
C            IF(NOOPTI.NE.1.AND.NOOPTI.NE.NTOPTI)THEN
C              XP(1,2,2*NTHELE*(NNRAELE+1)-(NOOPTI-2))=
C     '          XP(1,2,2*NTHELE*NNRAELE+NOOPTI)
C            ENDIF
C          ENDDO
C        ENDDO
CC ***     Updating the boundary conditions for the B.E. nodes.
CC ***     NPI=LAST FINITE ELEMENT NODE.
C        NPI=(2*NRAELE+1)*2*NTHELE
C        NY=0
C        DO NP=1,NPI           !  Looping over nodes
C          DO NH=1,NHP(NP)     !  Looping over dependent variables
C            DO NK=1,NKH(NH,NP)!  Looping over the derivatives
C              NY=NY+1
C            ENDDO
C          ENDDO
C        ENDDO
C        NY=NY+1
C        DO NP=NPI+1,NNODES
C          YP(NY,1)=DEXP(-S*XP(1,2,NP))
CC ***    This step is only valid for 2-D
C          NY=NY+2
C        ENDDO
C      ENDIF
C
C
 9998 CALL EXITS('UPBE')
      RETURN
 9999 CALL ERRORS('UPBE',ERROR)
      CALL EXITS('UPBE')
      RETURN 1
      END

      SUBROUTINE GRADPHI_I(IBT,IDO,INP,NBH,NBJ,
     '  NEELEM,NGAP,NHE,NHP,NJE,NKE,NKH,NLL,NP_INTERFACE,NPF,NPNE,
     '  NPNODE,NQE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,
     '  CE,DET,DL,DRDN,GRADPHI,PG,RAD,RG,SE,STEPSIZE,WG,
     '  XA,XE,XG1,XIG,XN,XP,XPFP,XR,YP,ZA,ZE,ZF,ZP,STOP,ERROR,*)

C#### Subroutine: GRADPHI_I
C###  Description:
C###    GRADPHI_I calculates the gradient of the potential at the 
C###    internal point XPFP(nj).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:bem000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:ktyp90.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM),NJE(NEM),NKE(NKM,NNM,NBFM,NEFM),
     '  NKH(NHM,NPM,NCM),NLL(12,NEM),NP_INTERFACE(0:NPM,0:3),
     '  NPF(15,NFM),NPNE(NNM,NBFM,NEFM),NPNODE(0:NP_R_M,0:NRM),
     '  NQE(NSM,NBFM,NEFM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEFM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEFM),NW(NEM,2),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),DET(NBFM,0:NNM,NGM,6),DL(3,NLM),DRDN(NGM),
     '  GRADPHI(*),PG(NSM,NUM,NGM,NBM),RAD(NGM),RG(NGM),
     '  SE(NSM,NBFM,NEFM),STEPSIZE,WG(NGM,NBM),XA(NAM,NJM,NQM),
     '  XE(NSM,NJM),XG1(NJM,NUM,NGM),XIG(NIM,NGM,NBM),XN(NJM,NGM),
     '  XP(NKM,NVM,NJM,NPM),XPFP(*),XR(NJM,NGM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZF(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL STOP
!     Local Variables
      INTEGER IADAPTIVE,intscheme,MK,nbuh,nbuhp,nb1j,nb1jp,nbbem,
     ' nbqh,nbqhp,ne,ng,nh,ni,NITB,nj,nj_grad,nk,nn,NNMIN,
     ' noelem,np,ns,nu,NU1(0:3)
      REAL*8 D(2),DGREEN,DXXI(3,3),HYPGREEN,MINDIST,
     ' RWG,SUM,SUMQ,SUMR,SUMU,SUMXG,XGC(3),XIMIN(2),
     ' XNO(3),XPC(3)
      LOGICAL ADAPTIVE,INTERFACE

      DATA NU1/1,2,4,7/

      CALL ENTERS('GRADPHI_I',*9999)

c cpb 21/6/95 This needs to be fixed for the new bem integration scheme.

      nh=NH_LOC(1,nx) !Temporary GMH 1-May-94
      DO nj=1,NJT
        GRADPHI(nj)=0.0d0
      ENDDO
      CALL EQTYPE(IBT,NBH,NEELEM,NJE,nr,NW,nx,ERROR,*9999)
      !Transfer global solution vector to element solution vector
      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,NYNP,
     '  YP,ZA,ZP,ERROR,*9999) !Transfer dependent variable
      CALL YPZP(2,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,NYNP,
     '    YP,ZA,ZP,ERROR,*9999)  !Transfer normal derivative
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        INTERFACE=.FALSE. !INTERFACE is true if element ne is on the
                          !interface between regions in a coupled 
                          !problem and the
                          !region to which ne belongs is NOT the 
                          !region with the smallest region number.
        np=NPNE(1,NBJ(1,ne),ne)
        IF((NP_INTERFACE(np,0).GT.1).AND.(NP_INTERFACE(np,1).NE.nr))
     '         INTERFACE=.TRUE.
        !Transfer global parameters to local parameters
        CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),
     '    NPF(1,1),NPNE(1,1,ne),NQE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA,XE,XP,ERROR,*9999)
        !Calculate the minimum distance from the node to the nodes 
        !in element ne
        CALL DIST(intscheme,NBJ(1,ne),NJE(ne),NLL(1,ne),nnmin,
     '      NPNE(1,1,ne),nr,NW(ne,2),DL,MINDIST,XP,XPFP,ERROR,*9999)
        IF(MINDIST.GT.STEPSIZE) THEN
          NNMIN=0  !Singularity is outside element.
        ELSE
          STOP=.TRUE.
          NNMIN=0 !Should never actually be in the element.
          !We are close to an element so stop tracking current line.
        ENDIF
        CALL QUADBE(intscheme,NBJ(1,ne),nbbem,ERROR,*9999)
        nb1jp=NFBASE(1,NBASEF(NBJ(1,ne),nbbem)) 
        nbuhp=NFBASE(1,NBASEF(NBH(NH_LOC(1,nx),1,ne),nbbem)) 
        nbqhp=NFBASE(1,NBASEF(NBH(NH_LOC(1,nx),2,ne),nbbem)) 
        ADAPTIVE=.FALSE. !Temporary
        IF(ADAPTIVE.AND.(NW(ne,2).NE.11).AND.(NW(ne,2).NE.13)) THEN
          !Don't use adaptive integration on distorted cubic linear
          !elements
          !Set up adaptive basis function based on Telles's rule.
          CALL GAUSS11(IBT,IDO,INP,
     '                 nb1jp,NGAP,D,DET,PG,XIG,XIMIN,ERROR,*9999)
          IADAPTIVE=1
        ELSE
          IADAPTIVE=0
        ENDIF

        nbuh=NBASEF(NBH(NH_LOC(1,nx),1,ne),nbbem) 
        nbqh=NBASEF(NBH(NH_LOC(1,nx),2,ne),nbbem) 
        nb1j=NBASEF(NBJ(1,ne),nbbem) 
        !Transfer global solution vector to element solution vector
        CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,2),nx,SE(1,1,ne),
     '    ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
        CALL ZPZE(NBH(1,2,ne),2,NHE(ne),NKE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,2),nx,SE(1,1,ne),
     '    ZA(1,1,1,ne),ZF,ZP,ERROR,*9999)
          !ZE(ns,nhx) contains the dependent variable values
          !ZF(ns,nh) contains the normal derivative values

        DO ng=1,NGT(nb1j)
          DO nj=1,NJE(ne)
            DO nu=1,NUT(nb1jp)
              SUMXG=0.0d0
              DO ns=1,NST(nb1jp)
                SUMXG=SUMXG+PG(ns,nu,ng,nb1j)*XE(ns,nj)
              ENDDO
              XG1(nj,nu,ng)=SUMXG
            ENDDO
          ENDDO
          !DXXI(nj,ni) contains the values of dXj/dXni at the
          !Gauss point.
          NITB=NIT(nb1jp)
          DO ni=1,NITB
            DO nj=1,NJE(ne)
              DXXI(nj,ni)=XG1(nj,NU1(ni),ng)
            ENDDO
          ENDDO
          IF(NITB.EQ.1) THEN
            RG(ng)=SQRT(DXXI(1,1)*DXXI(1,1)+DXXI(2,1)*DXXI(2,1))
            !Jacobian for 1D integral in 2D space
          ELSE !Calculate cross product of 2 tangent vectors
            RG(ng)=SQRT((DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))**2+
     '                  (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))**2+
     '                  (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))**2)
            !Jacobian for 2D integral in 3D space
          ENDIF
          IF(ITYP10(nr).GT.1) THEN
            !Trans into cartesian to find distance (RAD).
            DO nj=1,NJE(ne)
              XR(nj,ng)=XPFP(nj)
            ENDDO
            CALL COORD(ITYP10(nr),1,XR(1,ng),XPC,ERROR,*9999)
            DO nj=1,NJE(ne)
              XR(nj,ng)=XG1(nj,1,ng)
            ENDDO
            CALL COORD(ITYP10(nr),1,XR(1,ng),XGC,ERROR,*9999)
          ENDIF
          CALL NORMAL(NJE(ne),nr,XG1(1,1,ng),XN(1,ng),
     '      INTERFACE,ERROR,*9999)
          IF(ITYP10(nr).EQ.1) THEN
            DO nj=1,NJE(ne)
              XR(nj,ng)=XG1(nj,1,ng)-XPFP(nj)
            ENDDO
          ELSE
            DO nj=1,NJE(ne)
              XR(nj,ng)=XGC(nj)-XPC(nj)
            ENDDO
          ENDIF
        ENDDO !End of ng loop for variables not depending on which
              !direction the equation was differentiated.
        DO nj_grad=1,NJT !Loop over each derivative direction.
          DO nj=1,NJT
            XNO(nj)=0.0d0
          ENDDO
          XNO(nj_grad)=1.0d0
          DO ng=1,NGT(nb1j)
            SUM=0.0d0
            DO nj=1,NJE(ne)
              SUM=SUM+XR(nj,ng)**2
            ENDDO
            SUMR=XR(nj_grad,ng) !the njth component of XR
            RAD(ng)=DSQRT(SUM)
            DRDN(ng)=-SUMR/RAD(ng) !dR/dn0 where n0 is a unit vector 
                                   !in the nj_grad direction.
          ENDDO
          ns=0
          DO nn=1,NNT(nbuhp)     !Dependent variable loop
            DO nk=1,NKT(nn,nbuhp)
              MK=NKE(nk,nn,nbuhp,ne)
              ns=ns+1
              IF(mk.gt.0.and.
     '          nk.LE.MAX(NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr),1)) THEN
                SUMU=0.0d0
                DO ng=1,NGT(nbuh)
                  RWG=RG(ng)*WG(ng,nbuh)
                  SUMU=SUMU-RWG*PG(ns,1,ng,nbuh)*
     '             HYPGREEN(IGREN(nr),0,CE,RAD(ng),
     '             DET(nbuhp,NNMIN,ng,1+IADAPTIVE),
     '             XN(1,ng),XNO,XR(1,ng))
                ENDDO !End of ng loop
                GRADPHI(nj_grad)=SUMU*ZE(ns,1)+GRADPHI(nj_grad)
              ENDIF !End of MK > 0 loop
            ENDDO !End of nk loop
          ENDDO !End of nn loop
          ns=0
          DO nn=1,NNT(nbqhp) !normal derivative loop
            DO nk=1,NKT(nn,nbqhp)
              MK=NKE(nk,nn,nbqhp,ne)
              ns=ns+1
              IF(mk.gt.0.and.
     '          nk.LE.MAX(NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr),1)) THEN
                SUMQ=0.0d0
                DO ng=1,NGT(nbqh)
                  SUMQ=SUMQ+PG(ns,1,ng,nbqh)*
     '              DGREEN(IGREN(nr),1,nh,nh,CE,DRDN(ng),RAD(ng),
     '              DET(nbqhp,NNMIN,ng,1+IADAPTIVE),XN(1,ng),XR(1,NG))
     '              *RG(ng)*WG(ng,nbqh)
                ENDDO !End of ng loop
                GRADPHI(nj_grad)=SUMQ*ZF(ns,1)+GRADPHI(nj_grad)
              ENDIF !End of MK > 0 loop.
            ENDDO !End of nk loop
          ENDDO !End of nn loop
        ENDDO !End of nj_grad loop (i.e. derivative direction loop).
      ENDDO !End of element loop

      CALL EXITS('GRADPHI_I')
      RETURN
9999  CALL ERRORS('GRADPHI_I',ERROR)
      CALL EXITS('GRADPHI_I')
      RETURN 1
      END


      SUBROUTINE GRADPHI_N(NJP,NKH,NP_INTERFACE,np,nr,NYNP,
     '  GRADPHI,XG,XN,XP,YP,ERROR,*)

C#### Subroutine: GRADPHI_N
C###  Description:
C###    GRADPHI_N calculates the gradient of the potential at node np.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
!     Parameter List
      INTEGER NJP(NPM),NKH(NHM,NPM,NCM),NP_INTERFACE(0:NPM,0:3),np,
     '  nr,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 GRADPHI(*),XG(NJM,NUM),XN(NJM,NGM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,nj,nj_grad,NJE2,nk,M(5),nv,nx,ny1,ny2
      REAL*8 D,DSDX(3,3),DXDS(3,3),SUM
      LOGICAL INTERFACE2
      DATA M/1,2,3,1,2/

      CALL ENTERS('GRADPHI_N',*9999)

      nv=1                      !Temporary MPN 12-Nov-94
      nx=1                      !Temporary 
      DO nj=1,NJT
        GRADPHI(nj)=0.0D0
      ENDDO
      IF(NKH(NH_LOC(1,nx),np,1).GT.1) THEN  !need derivatives defined
                                !Calculate dSi/dXj at the point np.
                                !Need to firstly find dXj/dSi and then invert.
        NJE2=NJP(np)
        DO nj=1,NJE2
          DXDS(nj,1)=XP(2,nv,nj,np)
          DXDS(nj,2)=XP(3,nv,nj,np)
        ENDDO
        DO nj=1,NJE2
          XG(nj,2)=DXDS(nj,1)
          XG(nj,4)=DXDS(nj,2)
        ENDDO
        INTERFACE2=.FALSE.      !INTERFACE2 is true if node np is on the
                                !interface between regions in a coupled 
                                !problem and the region to which
                                !NE2 belongs is NOT the region with
                                !the smallest region number.
        IF((NP_INTERFACE(np,0).GT.1).AND.(NP_INTERFACE(np,1).NE.nr))
     '    INTERFACE2=.TRUE.
        CALL NORMAL(NJE2,nr,XG,XN(1,1),INTERFACE2,ERROR,*9999) 
                                !Find normal vector at np
        DO nj=1,NJE2
          DXDS(nj,NJE2)=XN(nj,1) !dXj/dn
        ENDDO
                                !Calculate DSDX(ni,nj)=dSi/dXj
        IF(NJE2.EQ.2) THEN       !Inverse of 2x2 matrix
          D=DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1)
          CALL ASSERT(DABS(D).GE.RDELTA,'  >>Zero det ?',
     '      ERROR,*9999)
          DSDX(1,1)=DXDS(2,2)/D
          DSDX(1,2)=-DXDS(1,2)/D
          DSDX(2,1)=-DXDS(2,1)/D
          DSDX(2,2)=DXDS(1,1)/D
        ELSE IF(NJE2.EQ.3) THEN   !Inverse of 3x3 matrix
          D=DXDS(1,1)*
     '      (DXDS(2,2)*DXDS(3,3)-DXDS(3,2)*DXDS(2,3))
     '      +DXDS(1,2)*
     '      (DXDS(2,3)*DXDS(3,1)-DXDS(3,3)*DXDS(2,1))
     '      +DXDS(1,3)*
     '      (DXDS(2,1)*DXDS(3,2)-DXDS(3,1)*DXDS(2,2))
          DO ni=1,3
            DO nj=1,3
              DSDX(ni,nj)=(DXDS(M(nj+1),M(ni+1))*DXDS(M(nj+2),
     '          M(ni+2))-DXDS(M(nj+2),M(ni+1))*DXDS(M(nj+1),
     '          M(ni+2)))/D
            ENDDO
          ENDDO
        ENDIF      
        
        DO nj_grad=1,NJE2
          SUM=0.0d0
          DO ni=1,NJE2-1
            nk=ni+1             !Not general
            ny1=NYNP(nk,1,NH_LOC(1,nx),np,0,1,nr) ! global variable number
            SUM=SUM+DSDX(ni,nj_grad)*YP(ny1,1)
          ENDDO
          ny2=NYNP(1,1,NH_LOC(1,nx),np,0,2,nr) !normal derivative contribution
          SUM=SUM+DSDX(NJE2,nj_grad)*YP(ny2,2)
          GRADPHI(nj_grad)=SUM
        ENDDO
      ENDIF                     !nkh loop

      CALL EXITS('GRADPHI_N')
      RETURN
 9999 CALL ERRORS('GRADPHI_N',ERROR)
      CALL EXITS('GRADPHI_N')
      RETURN 1
      END

Module feuser.f
===============

      SUBROUTINE USER_7(ERROR,*)

C#### Subroutine: USER_7
C###  Description:
C###    Evaluates helically wound cylinder formula: 
C###    fibre length prescribed

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'  
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IFAIL,N,NO_DAT
      REAL*8 AA,ANGLE,FIBRE,PITCH,PP,RADIUS,RLENGTH,TOL,TURNS,VOLUME,
     '  X,X02AJF,X_LARGEST,X_SMALLEST,Y
      REAL*8 A(4),X_IMAG(4),X_REAL(4)

      CALL ENTERS('USER_7',*9999)
      PARAM_NAME(1)='Volume'
      PARAM_NAME(2)='Fibre turns'
      PARAM_NAME(3)='Init.radius'
      NAME_COL(1)='FIBRE LENGTH'
      NAME_COL(2)='FIBRE PITCH'
      NAME_COL(3)='FIBRE ANGLE'
      NAME_COL(4)='CELL RADIUS'
      NAME_COL(5)='CELL LENGTH'
      NT_COL=5

      VOLUME=PARAM(1) !Cell volume V
      TURNS =PARAM(2) !No of fibre turns n
      RADIUS=PARAM(3) !Initial cell radius r

      RLENGTH=VOLUME/(PI*RADIUS**2)                  !Init cell length l
      PITCH=RLENGTH/(2.0d0*PI*RADIUS*TURNS)          !Initial fib pitch
      ANGLE=DATAN(PITCH)*180.0d0/PI                  !Initial fib angle
      FIBRE=2.0d0*PI*RADIUS*DSQRT(1.0d0+PITCH*PITCH) !Fib len per turn s
      DATA(1,1)=FIBRE
      DATA(1,2)=PITCH                             
      DATA(1,3)=ANGLE                              
      DATA(1,4)=RADIUS                     
      DATA(1,5)=RLENGTH          

      TOL=X02AJF()
      DO NO_DAT=2,NT_DAT   

C **    Fibre length s=FIBRE comes from 1st column
        FIBRE=DATA(NO_DAT,1)
        WRITE(OP_STRING,'('' fibre  ='',e12.4)') FIBRE
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C **    Solve cubic to find pitch p from fibre length s
        A(1)= 1.0d0
        A(2)= 0.0d0
        A(3)=-1.0d0
        IF(TURNS.GT.0.AND.FIBRE.GT.0) THEN
          AA=4.0d0*PI*VOLUME/(TURNS*FIBRE**3)
        ELSE
          AA=0.0d0
        ENDIF
        A(4)=AA
        N=4
        IFAIL=1
        CALL C02AEF(A,N,X_REAL,X_IMAG,TOL,IFAIL)
        IF(IFAIL.EQ.0) THEN
          WRITE(OP_STRING,'('' Roots:'',3(5X,D12.4,'','',D12.4),'
     '      //''' at a='',E12.4)') (X_REAL(I),X_IMAG(I),I=1,3),AA
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          X_SMALLEST=1.0d0
          X_LARGEST =0.0d0
          DO I=1,3
            X=DBLE(X_REAL(I))
            Y=DBLE(X_IMAG(I))
            IF(X.GT.0.0d0.AND.X.LT.1.0d0.AND.Y.EQ.0.0d0) THEN
              WRITE(OP_STRING,'('' Real root found in 0<x<1 at x='','
     '          //'F6.4)') X
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(X.LT.X_SMALLEST) X_SMALLEST=X
              IF(X.GT.X_LARGEST ) X_LARGEST =X
            ENDIF
          ENDDO
          WRITE(OP_STRING,'('' Smallest & largest real roots are '','
     '      //'2F10.4)') X_SMALLEST,X_LARGEST
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(AA.LE.0.384900d0.AND.X_SMALLEST.GT.0.0d0.AND.
     '      X_SMALLEST.LT.1.0d0) THEN
            X=X_SMALLEST
          ELSE IF(AA.GT.0.384900d0.AND.X_LARGEST.GT.0.0d0.AND.
     '      X_LARGEST.LT.1.0d0) THEN
            X=X_LARGEST
          ELSE
            WRITE(OP_STRING,
     '        '('' Set x=0 since no real x in 0<x<1'')')
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            X=0.0d0
          ENDIF
          WRITE(OP_STRING,'('' x      ='',e12.4)') X
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          PITCH=X/DSQRT(1.0d0-X*X)
          WRITE(OP_STRING,'('' pitch  ='',e12.4)') PITCH
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DATA(NO_DAT,2)=PITCH
          ANGLE=DATAN(PITCH)*180.0d0/PI
          DATA(NO_DAT,3)=ANGLE

C **      Calculate cell radius r from s and p
          PP=DSQRT(1.0d0+PITCH*PITCH)
          RADIUS=FIBRE/(2.0d0*PI*PP)
          WRITE(OP_STRING,'('' radius ='',e12.4)') RADIUS
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DATA(NO_DAT,4)=RADIUS

C **      Calculate cell length l from V and r
          RLENGTH=VOLUME/(PI*RADIUS**2)
          WRITE(OP_STRING,'('' rlength='',e12.4)') RLENGTH
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DATA(NO_DAT,5)=RLENGTH
        ELSE 
          WRITE(OP_STRING,
     '      '('' C02AEF failed with IFAIL='',I3)') IFAIL
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DATA(NO_DAT,2)=0.0d0
          DATA(NO_DAT,3)=0.0d0
          DATA(NO_DAT,4)=0.0d0
          DATA(NO_DAT,5)=0.0d0
        ENDIF
      ENDDO

      CALL EXITS('USER_7')
      RETURN
 9999 CALL ERRORS('USER_7',ERROR)
      CALL EXITS('USER_7')
      RETURN 1
      END


      SUBROUTINE LSFUN2(M,N,X_FIT,FVECC,FJACC,LJC)

C#### Subroutine: LSFUN2
C###  Description:
C###    Calculates function and first deriv matrix for Nag routine E04GEF

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
!     Parameter List
      INTEGER LJC,M,N
      REAL*8 FJACC(LJC,*),FVECC(*),X_FIT(*)
!     Local Variables
      CHARACTER ERROR*10
      INTEGER NO_DAT,NO_FIT,NO_PAR
      REAL*8 C1,C2,E11,RL1,S11

      DO NO_FIT=1,NT_FIT
        NO_PAR=NP_FIT(NO_FIT)
        PARAM(NO_PAR)=X_FIT(NO_FIT)
      ENDDO
      C1=PARAM(1)
      C2=PARAM(2)
C      C3=PARAM(3)
C      C4=PARAM(4)
C      C5=PARAM(5)
C      C6=PARAM(6)
      DO NO_DAT=1,NT_DAT
        RL1=DATA(NO_DAT,1)
        E11=0.5D0*(RL1-1.D0)
        S11=DATA(NO_DAT,3)
        FVECC(NO_DAT)  =C1*DEXP(C2*E11*E11)-S11
        FJACC(NO_DAT,1)=DEXP(C2*E11*E11)
        FJACC(NO_DAT,2)=C1*E11*E11*DEXP(C2*E11*E11)
        WRITE(OP_STRING,'(6E12.3)')
     '    RL1,E11,S11,FVECC(NO_DAT),FJACC(NO_DAT,1),FJACC(NO_DAT,2)
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDDO

 9999 RETURN
      END

C MHT 03-05-01 Archiving routine.
C      SUBROUTINE REASSTB(LD,LDR,nb,NE_OLD,NE_TERM,NPC_MIN,NPNE,N_ELM,
C     '  N_TRM,G_LL,MIN_DIST,XP,ERROR,*)
C
CC#### Subroutine: REASSTB
CC###  Description:
CC###    REASSTB reassigns random points that belong to a set for a
CC###    terminal branch.  Used for generating a Monte-Carlo lung tree.
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:geom00.cmn'
C
C      !Parameter list
C      INTEGER LD(NDM),LDR(0:NDM),nb,NE_OLD(2*NORM),NE_TERM(NE_R_M),
C     '  NPC_MIN(NP_R_M),NPNE(NNM,NBFM,NEM),N_ELM,N_TRM
C      REAL*8 G_LL,MIN_DIST(NE_R_M),XP(NKM,NVM,NJM,NPM)
C      CHARACTER ERROR*(*)
C      !Local variables
C      INTEGER nd,nj,nob,nonob,nontb,nop,nsp,ntb,ntp
C      REAL*8 DIST,PNT1(3),PNT2(3)
C
C      CALL ENTERS('REASSTB',*9999)
C
C      DO nontb=1,N_TRM
C        ntb=NE_TERM(nontb) !global element # of terminal branch
C        ntp=NPNE(2,nb,ntb) !end node # of terminal branch
C        DO nj=1,NJT
C          PNT1(nj)=XP(1,1,nj,ntp) !coordinates at end of TB
C        ENDDO !nj
C        MIN_DIST(ntb)=1000.0d0
C        DO nonob=1,N_ELM !find the closest neighbouring branch
C          DIST=0.0d0
C          nob=NE_OLD(nonob)
C          nop=NPNE(2,nb,nob)
C          DO nj=1,NJT
C            PNT2(nj)=XP(1,1,nj,nop)
C            DIST=DIST+(PNT1(nj)-PNT2(nj))**2.0d0
C          ENDDO !nj
C          DIST=DSQRT(DIST)
C          IF(DIST.LT.MIN_DIST(ntb))THEN
C            MIN_DIST(ntb)=DIST
C            NPC_MIN(ntb)=nob
C          ENDIF !DIST
C        ENDDO !nobranch
C      ENDDO !ntb
C      DO nd=1,NDT !reassign random points to closest branch
C        nsp=LDR(nd)
C        DO nontb=1,N_TRM
C          ntb=NE_TERM(nontb)
C          IF(nsp.EQ.ntb)THEN
C            IF(MIN_DIST(ntb).LT.G_LL)THEN
C              LD(nd)=NPC_MIN(ntb)
C            ENDIF !MIN_DIST
C          ENDIF !nsp
C        ENDDO !nontb (ntb)
C      ENDDO !nd
C
C      CALL EXITS('REASSTB')
C      RETURN
C 9999 CALL ERRORS('REASSTB',ERROR)
C      CALL EXITS('REASSTB')
C      RETURN 1
C      END



