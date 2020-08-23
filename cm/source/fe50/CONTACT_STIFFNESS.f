      SUBROUTINE CONTACT_STIFFNESS(IBT,IDO,INP,ISC_GK,ISR_GK,NBH,
     '  NBHF,NFF,NHE,NKEF,NKHE,NNF,NPF,NPNE,NRE,NVHE,nx,NYNP,
     '  Z_CONT_LIST,GK,SE,Z_CONT,FIX,ERROR,*)

C#### Subroutine: CONTACT_STIFFNESS
C###  Description:
C###    Calculates contact stiffness component and modifes global
C###    stiffness matrix GK.
C### 25/02/08 JHC For frictionless contact, K, the contact stiffness component is analytic linearisation
C###    of the contact residual, which results in
C###    K=increment_of_normal_contact_force + increment_of_variational_gap or
C###    K=inc_t_n + inc_var_gap
C### 13/03/08 JHC For frictional contact, K, the contact stiffness component is analytic linearisation
C###    of the frictioanl contact residual, which results in
C###    K=increment_of_tangential_contact_force * variational_XI + tangential_force * increment_of_variational_XI or
C###    K=inc_t_trial_1 * D(1) + inc_t_trial_2 * D(2) + t_1 * K_ct(1) + t_2 * K_ct(2)
C###    where inc_t_trial * D(1) + inc_t_trial_2 * D(2) is denoted as K_direct
C###    For more details, refer to the Computational Contact and Impact Mechanics by Tod Laursen

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'gen000.cmn'      
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'      
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nonl00.cmn'      
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     &  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),
     &  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NFF(6,NEM),
     &  NHE(NEM),NKEF(0:4,16,6,NBFM),NKHE(NKM,NNM,NHM,NEM),
     &  NNF(0:17,6,NBFM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     &  NVHE(NNM,NBFM,NHM,NEM),nx,
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),Z_CONT_LIST(NDM,2,7)
      REAL*8 GK(NZ_GK_M),SE(NSM,NBFM,NEM),
     &  Z_CONT(NDM,2,67)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)

!     Local Variables
      REAL*8 basem,basem_deri(2),basem_double_deri(2,2),
     &  basen,basen_deri(2),basen_double_deri(2,2),cov_m(2,2),
     &  D_m(2),D_n(2),DOT,EXI_contm(3),EXI_contn(3),
     &  G_n(2),inc_t_n,inc_t_slip(2),inc_t_trial(2),inc_var_gap,
     &  inv_A(2,2),inv_m(2,2),jacobian,
     &  K,K_ct(2),K_direct,kappa(2,2),lamda_n,
     &  N_m(2),N_n(2),Ntang_m(2,2),Ntang_n(2,2),
     &  norm_gap,normal_comp(3,2),
     &  omega(2,2,2),P_n(2),Pt(3),Pt_contra(2),Pt_covar(2),
     &  scalem,scalen,slip(2),T_m(2),T_n(2),
     &  Ttang_m(2,2),Ttang_n(2,2),
     &  tang_comp(3,2,2),tang_deri_comp(3,2,2),
     &  tmp_a,tmp_b,tmp_c,tmp_d,
     &  weight,XI(3,2),XI_contm(2),XI_contn(2)
      INTEGER i,j,m,n,nvr,nvc,nb_cont(2),ne_cont(2),nh,
     &  nr_cont(2),nkr,nkc,nnr,nnc,npr,npc,nyr,nyc,nf_cont(2),
     &  nsfr,nsfc,nsr,nsc,nef_cont(2),nnkr,nnkc,nh_y,nh_z,
     &  nz,STATUS,v,x,y,z
     
!     Storage of face information for master/slave setup in CALC_FACE_INFORMATION_DEP
      REAL*8 SF_cont(NSM,NBFM,2)
      INTEGER NKHF_cont(NKM,NNM,NHM,2),NPNF_cont(NNM,NBFM,2),
     &  NVHF_cont(NNM,NBFM,NHM,2)     

!     Functions
      REAL*8 PSI1

      CALL ENTERS('CONTACT_STIFFNESS',*9999)

C*** 22/02/08 JHC CONT_WT removed from common block nonl00.cmn
C      weight=1.0d0/DBLE(CONT_WT**2) ! integration point weight

      DO j=1,NDT   ! Loop over contact points
C*** 22/02/08 JHC allow different number of contact points at different face
        weight=1.0d0/DBLE(Z_CONT_LIST(j,1,6)) ! integration point weight
      
        IF ((Z_CONT_LIST(j,1,3).EQ.1).OR.
     '    (Z_CONT_LIST(j,1,4).EQ.2)) THEN
     
          DO i=1,2  ! Loop over Slave/Master                                                            
            ne_cont(i)=Z_CONT_LIST(j,i,1) ! global elem # for contact pt j
            nr_cont(i)=NRE(ne_cont(i)) ! region # for elem ne
            IF (KTYP5H(nr_cont(i)).EQ.0) THEN
              nef_cont(i)=Z_CONT_LIST(j,i,2) ! local face # for contact pt j
              nf_cont(i)=NFF(nef_cont(i),ne_cont(i)) ! global face #
            ENDIF
            DO nh=1,3   ! Loop over Z_CONT fields
              XI(nh,i)=Z_CONT(j,i,nh)
              normal_comp(nh,i)=Z_CONT(j,i,nh+3)
              tang_comp(nh,i,1)=Z_CONT(j,i,nh+6)
              tang_comp(nh,i,2)=Z_CONT(j,i,nh+9)
            ENDDO !nh
            norm_gap=Z_CONT(j,i,13)
C*** 09/10/08 JHC No need to re-initialise lamda_n here. Already done in Contact_Residual.f
C            IF (CONT_IT.EQ.1) THEN
C              Z_CONT(j,i,17)=0.0d0 ! Initialise normal contact force.
C            ENDIF
            lamda_n=Z_CONT(j,i,17) ! Current estimate of normal contact force.                          
            STATUS=Z_CONT_LIST(j,i,4) ! frictionless or tied contact
                                                       
          ENDDO !i

C*** 25/02/08 JHC tang_deri_comp(nh,1,2) is the derivative of tangent_1 vector
C                 wrt to the xi_2 direction                     
          DO nh=1,3  
            tang_deri_comp(nh,1,1)=Z_CONT(j,1,nh+30)       
            tang_deri_comp(nh,1,2)=Z_CONT(j,1,nh+33)       
            tang_deri_comp(nh,2,1)=Z_CONT(j,1,nh+36)       
            tang_deri_comp(nh,2,2)=Z_CONT(j,1,nh+39)
          ENDDO

C*** 25/02/08 JHC cov_m is the covariant metric tensor m calculated in UPDATA.f      
          cov_m(1,1)=Z_CONT(j,1,47)
          cov_m(1,2)=Z_CONT(j,1,48)
          cov_m(2,1)=Z_CONT(j,1,49)
          cov_m(2,2)=Z_CONT(j,1,50)

          inv_m(1,1)=Z_CONT(j,1,51)
          inv_m(1,2)=Z_CONT(j,1,52)
          inv_m(2,1)=Z_CONT(j,1,53)
          inv_m(2,2)=Z_CONT(j,1,54)

          inv_A(1,1)=Z_CONT(j,1,55)
          inv_A(1,2)=Z_CONT(j,1,56)
          inv_A(2,1)=Z_CONT(j,1,57)
          inv_A(2,2)=Z_CONT(j,1,58) 

C*** 18/03/08 JHC Additional variables for frictional contact
C*** 29/09/08 JHC These variables now needed for tied contact
          IF (STATUS.EQ.3.OR.STATUS.EQ.2) THEN ! frictional or tied contact
            DO nh=1,2
              Pt_contra(nh)=Z_CONT(j,1,nh+44) 
              Pt_covar(nh)=Z_CONT(j,1,nh+63) 
            ENDDO

            DO nh=1,3 
              Pt(nh)=Z_CONT(j,1,nh+60)
            ENDDO 

C***        slip is the relative movement of xi points between each Newton step
C***        computed in UPDATA
            slip(1)=Z_CONT(j,1,43)
            slip(2)=Z_CONT(j,1,44)
          ENDIF

C  ***    Calculate face info for slave and master surfaces
          DO i=1,2

            IF (KTYP5H(nr_cont(i)).EQ.0) THEN
                                          
              CALL CALC_FACE_INFORMATION_DEP(NBH(1,1,ne_cont(i)),
     '          NBHF(1,1,nf_cont(i)),nef_cont(i),NHE(ne_cont(i)),
     '          NKHE(1,1,1,ne_cont(i)),NKEF,NKHF_cont(1,1,1,i),NNF,
     '          NPNE(1,1,ne_cont(i)),NPNF_cont(1,1,i),
     '          NVHE(1,1,1,ne_cont(i)),NVHF_cont(1,1,1,i),nx,
     '          SE(1,1,ne_cont(i)),SF_cont(1,1,i),ERROR,*9999) 
            ENDIF !KTYP5H

          ENDDO !i    

C ***     Modify GK with contact contribution
          DO m=1,2
            DO n=1,2
              DO nh_y=1,3
                DO nh_z=1,3
                                               
C*** 25/02/08 JHC Calculation of tangent stiffness components moved inside the nodal loop
C                  IF (STATUS.EQ.1) THEN ! frictionless contact
C                    IF (norm_gap.LE.0.0d0) THEN !outside volume
C                      K=(CONT_STIFF*exp((CONT_STIFF/F_hat)*norm_gap))
C     '                  *normal_comp(nh_y,m)*normal_comp(nh_z,n)
C                    ELSE !inside volume
C                      K=CONT_STIFF*normal_comp(nh_y,m)
C     '                  *normal_comp(nh_z,n)
C                    ENDIF 
C                  ELSEIF (STATUS.EQ.2) THEN ! tied contact
C                    K=TIED_STIFF*normal_comp(nh_y,m)*
C     '                normal_comp(nh_z,n)+TIED_STIFF*
C     '                tang_comp(nh_y,m,1)*tang_comp(nh_z,n,1)+
C     '                TIED_STIFF*tang_comp(nh_y,m,2)*tang_comp(nh_z,n,2)
C                  ELSEIF (STATUS.EQ.3) THEN ! contact with friction   
C                    IF (norm_gap.LE.0.0d0) THEN !outside volume
C                      K=(CONT_STIFF*exp((CONT_STIFF/F_hat)*norm_gap))*
C     '                  normal_comp(nh_y,m)*normal_comp(nh_z,n)+
C     '                  (CONT_STIFF*exp((CONT_STIFF/F_hat)*norm_gap))*
C     '                  FRIC_COEFF*tang_comp(nh_y,m,1)*
C     '                  tang_comp(nh_z,n,1)+
C     '                  (CONT_STIFF*exp((CONT_STIFF/F_hat)*norm_gap))*
C     '                  FRIC_COEFF*tang_comp(nh_y,m,2)*
C     '                  tang_comp(nh_z,n,2)     
C                    ELSE !inside volume
C                      K=CONT_STIFF*
C     '                  normal_comp(nh_y,m)*normal_comp(nh_z,n)+
C     '                  CONT_STIFF*FRIC_COEFF*
C     '                  tang_comp(nh_y,m,1)*tang_comp(nh_z,n,1)+
C     '                  CONT_STIFF*FRIC_COEFF*
C     '                  tang_comp(nh_y,m,2)*tang_comp(nh_z,n,2) 
C                    ENDIF           
C                  ENDIF !STATUS

                  IF (KTYP5H(nr_cont(m)).EQ.1) THEN
                    nb_cont(m)=NBH(nh_y,1,ne_cont(m)) ! element bases function #
                  ELSE
                    nb_cont(m)=NBHF(nh_y,1,nf_cont(m)) ! face bases function #
                  ENDIF

                  IF (KTYP5H(nr_cont(n)).EQ.1) THEN
                    nb_cont(n)=NBH(nh_z,1,ne_cont(n)) ! element bases function #
                  ELSE
                    nb_cont(n)=NBHF(nh_z,1,nf_cont(n)) ! face bases function #
                  ENDIF

                  DO nnr=1,NNT(nb_cont(m))   ! weight contact pt out to global nodes
                    DO nnkr=1,NKT(nnr,nb_cont(m))   ! loop over nk

                      IF (KTYP5H(nr_cont(m)).EQ.1) THEN
                        nkr=NKHE(nnkr,nnr,nh_y,ne_cont(m))
                        nvr=NVHE(nnr,nb_cont(m),nh_y,ne_cont(m)) 
                      ELSE
                        nkr=NKHF_cont(nnkr,nnr,nh_y,m)
                        nvr=NVHF_cont(nnr,nb_cont(m),nh_y,m)  ! version #                                                                                                       
                      ENDIF

                      DO nnc=1,NNT(nb_cont(n))
                        DO nnkc=1,NKT(nnc,nb_cont(n))

                          IF (KTYP5H(nr_cont(n)).EQ.1) THEN
                            nkc=NKHE(nnkc,nnc,nh_z,ne_cont(n))
                            nvc=NVHE(nnc,nb_cont(n),nh_z,ne_cont(n))  ! version #          
                          ELSE
                            nkc=NKHF_cont(nnkc,nnc,nh_z,n)
                            nvc=NVHF_cont(nnc,nb_cont(n),nh_z,n)  ! version #                                                                                                        
                          ENDIF

                          jacobian=Z_CONT(j,2,16)

                          IF (KTYP5H(nr_cont(m)).EQ.1) THEN
                            EXI_contm(1)=XI(1,m)
                            EXI_contm(2)=XI(2,m)
                            EXI_contm(3)=XI(3,m)
C*** 01/09/2008 JHC should be using correct basis function number
C                            basem=PSI1(IBT,IDO,INP,nb_cont(m),
C     '                        1,nnkr,nnr,EXI_contm)
                            basem=PSI1(IBT(1,1,nb_cont(m)),
     &                        IDO(1,1,0,nb_cont(m)),INP(1,1,nb_cont(m)),
     &                        nb_cont(m),1,nnkr,nnr,EXI_contm)
                          ELSE
                            XI_contm(1)=XI(NPF(1,nf_cont(m)),m)
                            XI_contm(2)=XI(NPF(3,nf_cont(m)),m)
C*** 01/09/2008 JHC should be using correct basis function number
C                            basem=PSI1(IBT,IDO,INP,nb_cont(m),
C     '                        1,nnkr,nnr,XI_contm)
                            basem=PSI1(IBT(1,1,nb_cont(m)),
     &                        IDO(1,1,0,nb_cont(m)),
     &                        INP(1,1,nb_cont(m)),nb_cont(m),
     &                        1,nnkr,nnr,XI_contm)

C*** 25/02/08 JHC derivatives of basis functions wrt xi directions in the face and the double derivatives
C*** 01/09/2008 JHC should be using correct basis function number
C                            basem_deri(1)=PSI1(IBT,IDO,INP,nb_cont(m),
C     &                        2,nnkr,nnr,XI_contm)
C                            basem_deri(2)=PSI1(IBT,IDO,INP,nb_cont(m),
C     &                        4,nnkr,nnr,XI_contm)
                            basem_deri(1)=PSI1(IBT(1,1,nb_cont(m)),
     &                        IDO(1,1,0,nb_cont(m)),INP(1,1,nb_cont(m)),
     &                        nb_cont(m),2,nnkr,nnr,XI_contm)
                            basem_deri(2)=PSI1(IBT(1,1,nb_cont(m)),
     &                        IDO(1,1,0,nb_cont(m)),INP(1,1,nb_cont(m)),
     &                        nb_cont(m),4,nnkr,nnr,XI_contm)

C*** 01/09/2008 JHC should be using correct basis function number
C                            basem_double_deri(1,1)=PSI1(IBT,IDO,INP, ! d^2base/dxi1^2
C     &                        nb_cont(m),3,nnkr,nnr,XI_contm)
C                            basem_double_deri(1,2)=PSI1(IBT,IDO,INP, ! d^2base/dxi1dxi2
C     &                        nb_cont(m),6,nnkr,nnr,XI_contm)
C                            basem_double_deri(2,1)=
C     &                        basem_double_deri(1,2)
C                            basem_double_deri(2,2)=PSI1(IBT,IDO,INP,! d^2base/dxi2^2
C     &                        nb_cont(m),5,nnkr,nnr,XI_contm)
                            basem_double_deri(1,1)=PSI1(
     &                        IBT(1,1,nb_cont(m)),IDO(1,1,0,nb_cont(m)),
     &                        INP(1,1,nb_cont(m)), 
     &                        nb_cont(m),3,nnkr,nnr,XI_contm) ! d^2base/dxi1^2
                            basem_double_deri(1,2)=PSI1(
     &                        IBT(1,1,nb_cont(m)),IDO(1,1,0,nb_cont(m)),
     &                        INP(1,1,nb_cont(m)), 
     &                        nb_cont(m),6,nnkr,nnr,XI_contm) ! d^2base/dxi1dxi2
                            basem_double_deri(2,1)=
     &                        basem_double_deri(1,2)
                            basem_double_deri(2,2)=PSI1(
     &                        IBT(1,1,nb_cont(m)),IDO(1,1,0,nb_cont(m)),
     &                        INP(1,1,nb_cont(m)),
     &                        nb_cont(m),5,nnkr,nnr,XI_contm) ! d^2base/dxi2^2
                          ENDIF

                          IF (KTYP5H(nr_cont(n)).EQ.1) THEN
                            EXI_contn(1)=XI(1,n)
                            EXI_contn(2)=XI(2,n)
                            EXI_contn(3)=XI(3,n)
C*** 01/09/2008 JHC should be using correct basis function number
C                            basen=PSI1(IBT,IDO,INP,nb_cont(n),
C     '                        1,nnkc,nnc,EXI_contn)
                            basen=PSI1(IBT(1,1,nb_cont(n)),
     &                        IDO(1,1,0,nb_cont(n)),INP(1,1,nb_cont(n)),
     &                        nb_cont(n),1,nnkc,nnc,EXI_contn)
                          ELSE
                            XI_contn(1)=XI(NPF(1,nf_cont(n)),n)
                            XI_contn(2)=XI(NPF(3,nf_cont(n)),n)
C*** 01/09/2008 JHC should be using correct basis function number
C                            basen=PSI1(IBT,IDO,INP,nb_cont(n),1,
C     '                        nnkc,nnc,XI_contn)
                            basen=PSI1(IBT(1,1,nb_cont(n)),
     &                        IDO(1,1,0,nb_cont(n)),INP(1,1,nb_cont(n)),
     &                        nb_cont(n),1,nnkc,nnc,XI_contn)

C*** 25/02/08 JHC derivatives of basis functions wrt xi directions in the face
C*** 01/09/2008 JHC should be using correct basis function number
C                            basen_deri(1)=PSI1(IBT,IDO,INP,nb_cont(n),
C     &                        2,nnkc,nnc,XI_contn)
C                            basen_deri(2)=PSI1(IBT,IDO,INP,nb_cont(n),
C     &                        4,nnkc,nnc,XI_contn)
                            basen_deri(1)=PSI1(IBT(1,1,nb_cont(n)),
     &                        IDO(1,1,0,nb_cont(n)),INP(1,1,nb_cont(n)),
     &                        nb_cont(n),2,nnkc,nnc,XI_contn)
                            basen_deri(2)=PSI1(IBT(1,1,nb_cont(n)),
     &                        IDO(1,1,0,nb_cont(n)),INP(1,1,nb_cont(n)),
     &                        nb_cont(n),4,nnkc,nnc,XI_contn)

C*** 01/09/2008 JHC should be using correct basis function number
C                            basen_double_deri(1,1)=PSI1(IBT,IDO,INP, ! d^2base/dxi1^2
C     &                        nb_cont(n),3,nnkc,nnc,XI_contn)
C                            basen_double_deri(1,2)=PSI1(IBT,IDO,INP, ! d^2base/dxi1dxi2
C     &                        nb_cont(n),6,nnkc,nnc,XI_contn)
C                            basen_double_deri(2,1)=
C     &                        basen_double_deri(1,2)
C                            basen_double_deri(2,2)=PSI1(IBT,IDO,INP,! d^2base/dxi2^2
C     &                        nb_cont(n),5,nnkc,nnc,XI_contn)
                            basen_double_deri(1,1)=PSI1(
     &                        IBT(1,1,nb_cont(n)),IDO(1,1,0,nb_cont(n)),
     &                        INP(1,1,nb_cont(n)),
     &                        nb_cont(n),3,nnkc,nnc,XI_contn) ! d^2base/dxi1^2
                            basen_double_deri(1,2)=PSI1(
     &                        IBT(1,1,nb_cont(n)),IDO(1,1,0,nb_cont(n)),
     &                        INP(1,1,nb_cont(n)),
     &                        nb_cont(n),6,nnkc,nnc,XI_contn) ! d^2base/dxi1dxi2
                            basen_double_deri(2,1)=
     &                        basen_double_deri(1,2)
                            basen_double_deri(2,2)=PSI1(
     &                        IBT(1,1,nb_cont(n)),IDO(1,1,0,nb_cont(n)),
     &                        INP(1,1,nb_cont(n)),
     &                        nb_cont(n),5,nnkc,nnc,XI_contn) ! d^2base/dxi2^2
                          ENDIF

                          IF (KTYP5H(nr_cont(m)).EQ.1) THEN
                            npr=NPNE(nnr,nb_cont(m),ne_cont(m))   ! global node #                                                                                                 
                          ELSE
                            npr=NPNF_cont(nnr,nb_cont(m),m)   ! global node #                                                                                                        
                          ENDIF

                          IF (KTYP5H(nr_cont(n)).EQ.1) THEN
                            npc=NPNE(nnc,nb_cont(n),ne_cont(n))   ! global node #    
                          ELSE
                            npc=NPNF_cont(nnc,nb_cont(n),n)   ! global node #    
                          ENDIF 
                          nyr=NYNP(nkr,nvr,nh_y,npr,1,1,nr_cont(m))
                          nyc=NYNP(nkc,nvc,nh_z,npc,1,1,nr_cont(n))

                          IF (KTYP5H(nr_cont(m)).EQ.1) THEN
                            nsr=nnkr+(nnr-1)*NKT(nnr,nb_cont(m))
                            scalem=SE(nsr,nb_cont(m),ne_cont(m))    
                          ELSE
                            nsfr=nnkr+(nnr-1)*NKT(nnr,nb_cont(m))
                            scalem=SF_cont(nsfr,nb_cont(m),m)    
                          ENDIF

                          IF (KTYP5H(nr_cont(n)).EQ.1) THEN
                            nsc=nnkc+(nnc-1)*NKT(nnc,nb_cont(n))
                            scalen=SE(nsc,nb_cont(n),ne_cont(n))    
                          ELSE
                            nsfc=nnkc+(nnc-1)*NKT(nnc,nb_cont(n))
                            scalen=SF_cont(nsfc,nb_cont(n),n) 
                          ENDIF

C*** 25/02/08 JHC Refer to Computational Contact and Impact Mechanics by Tod Laursen (chapter 5) for implementations
C                 Variable names are made to be consistent with variables used in the book             
                          IF((STATUS.EQ.1).OR.(STATUS.EQ.3)) THEN ! for frictionless or frictional

C***                        **** FRICTIONLESS PART ****      
                            DO x=1,2
                              DO y=1,2
                                DOT=0.0d0
                                DO nh=1,3
                                  DOT=DOT+
     &                             tang_deri_comp(nh,x,y)*
     &                             normal_comp(nh,1)
                                ENDDO
                                kappa(x,y)=DOT
                              ENDDO
                            ENDDO

C*** 28/02/08 JHC Penalty method for frictionless contact mechanics
                            IF(PENALTY.EQ.1) THEN ! penalty method
                              IF(AUGMENT.EQ.0) THEN ! AUGMENT must be used instead of AUG_IT
                                IF(norm_gap.LT.0.0d0) THEN
                                  inc_t_n=0.0d0
                                ELSE
                                  inc_t_n=CONT_STIFF*
     &                              normal_comp(nh_y,m)*scalem*basem*
     &                              normal_comp(nh_z,n)*scalen*basen
                                ENDIF
                              ELSE ! Augmented Lagrangian Method
                                IF(Z_CONT(j,1,46).LT.0.0d0) THEN 
                                  inc_t_n=0.0d0
                                ELSE
                                  inc_t_n=CONT_STIFF*
     &                              normal_comp(nh_y,m)*scalem*basem*
     &                              normal_comp(nh_z,n)*scalen*basen
                                ENDIF
                              ENDIF

                            ELSE ! X-constraint method 

                              IF(norm_gap.LT.0.0d0) THEN !outside volume
                                IF(DABS(lamda_n).LE.ZERO_TOL) THEN ! division by zero
                                  inc_t_n=0.0d0
                                ELSE
                                  inc_t_n=(CONT_STIFF
     '                              *exp((CONT_STIFF/lamda_n)*norm_gap))
     '                              *normal_comp(nh_y,m)*scalem*basem
     '                              *normal_comp(nh_z,n)*scalen*basen
                                ENDIF
                              ELSE !inside volume
                                inc_t_n=CONT_STIFF
     '                            *normal_comp(nh_y,m)*scalem*basem
     '                            *normal_comp(nh_z,n)*scalen*basen
                              ENDIF 
                            ENDIF

                            DO x=1,2
                              IF(m.EQ.1) THEN
                                N_m(x)=0.0d0
                              ELSE
                                N_m(x)=scalem*basem_deri(x)*
     &                            normal_comp(nh_y,m)
                              ENDIF
                              IF(n.EQ.1) THEN
                                N_n(x)=0.0d0
                              ELSE
                                N_n(x)=scalen*basen_deri(x)*
     &                            normal_comp(nh_z,n)
                              ENDIF
                              T_m(x)=scalem*basem*tang_comp(nh_y,m,x)
                              T_n(x)=scalen*basen*tang_comp(nh_z,n,x)
                            ENDDO
                            DO x=1,2 
                              D_m(x)=0.0d0
                              D_n(x)=0.0d0
                              DO y=1,2  
                                D_m(x)=D_m(x)+inv_A(x,y)
     &                            *(T_m(y)+norm_gap*N_m(y))
                                D_n(x)=D_n(x)+inv_A(x,y)
     &                            *(T_n(y)+norm_gap*N_n(y))
                              ENDDO
                            ENDDO

                            inc_var_gap=0.0d0

C*** 03/03/08 JHC Only add this term for non-host mesh problems
C*** 25/09/08 JHC Add the geometric term after user-specified Newton step
                            IF((KTYP5H(nr_cont(m)).EQ.0).AND.
     &                        (KTYP5H(nr_cont(n)).EQ.0).AND.
     &                        (CONV_IT.GE.ADD_GEOM)) THEN
                              DO x=1,2
                                DO y=1,2
                                  tmp_a=0.0d0
                                  tmp_b=0.0d0
                                  DO z=1,2
                                    tmp_a=tmp_a+kappa(x,z)*D_m(z)
                                    tmp_b=tmp_b+kappa(y,z)*D_n(z)
                                  ENDDO ! z

                                  inc_var_gap=inc_var_gap+
     &                              kappa(y,x)*D_m(y)*D_n(x)+
     &                              norm_gap*inv_m(x,y)*
     &                              (N_m(x)-tmp_a)*
     &                              (N_n(y)-tmp_b)  
                                ENDDO ! y

                                inc_var_gap=inc_var_gap-
     &                            D_m(x)*N_n(x)-N_m(x)*D_n(x)
                              ENDDO ! x
                            ENDIF ! non-host mesh

                            K=inc_t_n+inc_var_gap*Z_CONT(j,1,18)

C*** 18/03/08 JHC Added Coulomb's frictional contact
                            IF(STATUS.EQ.3) THEN
C*** Not implemented (at least not tested) with host-mesh problems
                              IF((KTYP5H(nr_cont(m)).EQ.1).OR.
     &                          (KTYP5H(nr_cont(n)).EQ.1)) THEN
                                CALL ASSERT(.FALSE.,
     &                            '>>Frictional contact not tested for '
     &                            //'host meshes',ERROR,*9999)
                              ENDIF

                              DO x=1,2
                                DO y=1,2
                                  DO z=1,2
                                    DOT=0.0d0
                                    DO nh=1,3
                                      DOT=DOT+
     &                                 tang_deri_comp(nh,x,y)*
     &                                 tang_comp(nh,1,z)
                                    ENDDO
                                    omega(x,y,z)=DOT
                                  ENDDO
                                ENDDO
                              ENDDO
                              
                              DO x=1,2
                                DO y=1,2
                                  IF (m.EQ.1) THEN
                                    Ttang_m(x,y)=0.0d0
                                    Ntang_m(x,y)=0.0d0
                                  ELSE
                                    Ttang_m(x,y)=basem_deri(y)*
     &                                scalem*tang_comp(nh_y,m,x)
                                    Ntang_m(x,y)=basem_double_deri(x,y)*
     &                                scalem*normal_comp(nh_y,m)
                                  ENDIF

                                  IF (n.EQ.1) THEN
                                    Ttang_n(x,y)=0.0d0
                                    Ntang_n(x,y)=0.0d0
                                  ELSE
                                    Ttang_n(x,y)=basen_deri(y)*
     &                                scalen*tang_comp(nh_z,n,x)
                                    Ntang_n(x,y)=basen_double_deri(x,y)*
     &                                scalen*normal_comp(nh_z,n)
                                  ENDIF
                                ENDDO   

                                IF (n.EQ.1) THEN
                                  P_n(x)=0.0d0
                                ELSE
                                  P_n(x)=basen_deri(x)*scalen*Pt(nh_z)
                                ENDIF
                              ENDDO

                              DO x=1,2
                                K_ct(x)=0.0d0
                                tmp_a=0.0d0
                                tmp_b=0.0d0
                                DO y=1,2
                                  K_ct(x)=K_ct(x)+
     &                              Ttang_m(y,x)*D_n(y)+
     &                              D_m(y)*Ttang_n(y,x)+
     &                              Ttang_m(x,y)*D_n(y)+
     &                              D_m(y)*Ttang_n(x,y)+
     &                              norm_gap*
     &                              (Ntang_m(x,y)*D_n(y)+
     &                              D_m(y)*Ntang_n(x,y))
                                  tmp_a=tmp_a+kappa(x,y)*D_m(y)
                                  tmp_b=tmp_b+kappa(x,y)*D_n(y)
                                  DO z=1,2
                                    K_ct(x)=K_ct(x)-
     &                                omega(y,z,x)*D_m(y)*D_n(z)-
     &                                omega(x,z,y)*
     &                                (D_m(y)*D_n(z)+D_m(z)*D_n(y))
                                    tmp_c=0.0d0
                                    tmp_d=0.0d0
                                    DO v=1,2
                                      tmp_c=tmp_c+omega(v,z,z)*D_m(v)
                                      tmp_d=tmp_d+omega(v,z,z)*D_n(v)
                                    ENDDO
                                    K_ct(x)=K_ct(x)-inv_m(y,z)*T_m(y)*
     &                                (Ttang_n(z,x)-tmp_d)-inv_m(y,z)*
     &                                (Ttang_m(z,x)-tmp_c)*T_n(y) 
                                  ENDDO
                                ENDDO
                                K_ct(x)=K_ct(x)-normal_comp(nh_y,m)*
     &                            scalem*basem*(N_n(x)-tmp_b)-
     &                            (N_m(x)-tmp_a)*normal_comp(nh_z,n)*
     &                            scalen*basen
                              ENDDO

                              K_direct=0.0d0

                              DO x=1,2
                                G_n(x)=0.0d0
                                DO y=1,2
                                  G_n(x)=G_n(x)-slip(y)*(
     &                              Ttang_n(x,y)+Ttang_n(y,x))
                                  DO z=1,2
                                    G_n(x)=G_n(x)+slip(y)*
     &                                (omega(x,z,y)+omega(y,z,x))*
     &                                D_n(z)
                                  ENDDO ! z
                                ENDDO ! y
                              ENDDO ! x

                              DO x=1,2
                                inc_t_trial(x)=G_n(x)
                                DO y=1,2
                                  inc_t_trial(x)=inc_t_trial(x)+
     &                              (cov_m(x,y)*D_n(y))
                                ENDDO
                                inc_t_trial(x)=FRIC_STIFF*
     &                            inc_t_trial(x)
                              ENDDO

                              IF (Z_CONT(j,1,60).GT.0.0d0) THEN ! stick
                                DO x=1,2
                                  K_direct=K_direct+
     &                              D_m(x)*inc_t_trial(x)
                                ENDDO
                              ELSE                              ! slip
                                DO x=1,2
                                  IF (AUGMENT.EQ.0) THEN ! penalty method
                                    IF (norm_gap.LT.0.0d0) THEN
                                      inc_t_slip(x)=0.0d0
                                    ELSE
                                      inc_t_slip(x)=-CONT_STIFF*
     &                                  FRIC_COEFF*Pt_covar(x)*
     &                                  normal_comp(nh_z,n)*basen*
     &                                  scalen
                                    ENDIF
                                  ELSE ! Augmented Lagrangian Method    
                                    IF (Z_CONT(j,1,46).LT.0.0d0) THEN
                                      inc_t_slip(x)=0.0d0
                                    ELSE
                                      inc_t_slip(x)=-CONT_STIFF*
     &                                  FRIC_COEFF*Pt_covar(x)*
     &                                  normal_comp(nh_z,n)*basen*
     &                                  scalen
                                    ENDIF
                                  ENDIF

                                  DO y=1,2
                                    IF (x.EQ.y) THEN   
                                      inc_t_slip(x)=inc_t_slip(x)+
     &                                  FRIC_COEFF*(Z_CONT(j,1,18)/ 
     &                                  Z_CONT(j,1,59))*             ! Z_CONT(j,1,59) is the magnitude of trial tang. forces
     &                                  inc_t_trial(y)*(1.0d0-
     &                                  Pt_contra(y)*Pt_covar(x))
                                    ELSE
                                      inc_t_slip(x)=inc_t_slip(x)-
     &                                  FRIC_COEFF*(Z_CONT(j,1,18)/ 
     &                                  Z_CONT(j,1,59))*             ! Z_CONT(j,1,59) is the magnitude of trial tang. forces
     &                                  inc_t_trial(y)*
     &                                  Pt_contra(y)*Pt_covar(x)
                                    ENDIF

                                    tmp_a=0.0d0
                                    DO z=1,2
                                      DOT=0.0d0
                                      DO nh=1,3
                                        DOT=DOT+Pt(nh)*
     &                                    tang_deri_comp(nh,y,z) 
                                      ENDDO
                                      tmp_a=tmp_a+DOT*D_n(z)
                                    ENDDO
                                    inc_t_slip(x)=inc_t_slip(x)+
     &                                FRIC_COEFF*Z_CONT(j,1,18)*
     &                                Pt_contra(y)*Pt_covar(x)*
     &                                (P_n(y)+tmp_a)
                                  ENDDO ! y
                                ENDDO ! x

                                DO x=1,2
                                  K_direct=K_direct+
     &                              D_m(x)*inc_t_slip(x)
                                ENDDO   
                              ENDIF ! Stick or Slip

                              K=K+K_direct+( 
     &                          (Z_CONT(j,1,20)*inv_A(1,1)+              ! (Z_CONT(j,1,20) is t_1
     &                          Z_CONT(j,1,22)*inv_A(2,1))*K_ct(1)+      ! (Z_CONT(j,1,22) is t_2
     &                          (Z_CONT(j,1,20)*inv_A(1,2)+
     &                          Z_CONT(j,1,22)*inv_A(2,2))*K_ct(2))    

                            ENDIF ! STATUS=3 frictional contact

                          ELSEIF (STATUS.EQ.2) THEN ! tied contact
C*** 29/09/08 JHC Removed Justin's tied contact
C                            K=(TIED_STIFF*normal_comp(nh_y,m)*
C     '                        normal_comp(nh_z,n)+TIED_STIFF*
C     '                        tang_comp(nh_y,m,1)*tang_comp(nh_z,n,1)+
C     '                        TIED_STIFF*tang_comp(nh_y,m,2)*
C     '                        tang_comp(nh_z,n,2))*
C     '                        basem*basen*scalem*scalen

C*** 29/09/08 JHC Modify tied contact to apply compression/tension 

C***                        **** FRICTIONLESS PART ****      
                            DO x=1,2
                              DO y=1,2
                                DOT=0.0d0
                                DO nh=1,3
                                  DOT=DOT+
     &                             tang_deri_comp(nh,x,y)*
     &                             normal_comp(nh,1)
                                ENDDO
                                kappa(x,y)=DOT
                              ENDDO
                            ENDDO

                            inc_t_n=CONT_STIFF*
     &                        normal_comp(nh_y,m)*scalem*basem*
     &                        normal_comp(nh_z,n)*scalen*basen

                            DO x=1,2
                              IF(m.EQ.1) THEN
                                N_m(x)=0.0d0
                              ELSE
                                N_m(x)=scalem*basem_deri(x)*
     &                            normal_comp(nh_y,m)
                              ENDIF
                              IF(n.EQ.1) THEN
                                N_n(x)=0.0d0
                              ELSE
                                N_n(x)=scalen*basen_deri(x)*
     &                            normal_comp(nh_z,n)
                              ENDIF
                              T_m(x)=scalem*basem*tang_comp(nh_y,m,x)
                              T_n(x)=scalen*basen*tang_comp(nh_z,n,x)
                            ENDDO
                            DO x=1,2 
                              D_m(x)=0.0d0
                              D_n(x)=0.0d0
                              DO y=1,2  
                                D_m(x)=D_m(x)+inv_A(x,y)
     &                            *(T_m(y)+norm_gap*N_m(y))
                                D_n(x)=D_n(x)+inv_A(x,y)
     &                            *(T_n(y)+norm_gap*N_n(y))
                              ENDDO
                            ENDDO

                            inc_var_gap=0.0d0

C*** 03/03/08 JHC Only add this term for non-host mesh problems
                            IF((KTYP5H(nr_cont(m)).EQ.0).AND.
     &                        (KTYP5H(nr_cont(n)).EQ.0).AND.
     &                        (CONV_IT.GE.ADD_GEOM)) THEN
                              DO x=1,2
                                DO y=1,2
                                  tmp_a=0.0d0
                                  tmp_b=0.0d0
                                  DO z=1,2
                                    tmp_a=tmp_a+kappa(x,z)*D_m(z)
                                    tmp_b=tmp_b+kappa(y,z)*D_n(z)
                                  ENDDO ! z

                                  inc_var_gap=inc_var_gap+
     &                              kappa(y,x)*D_m(y)*D_n(x)+
     &                              norm_gap*inv_m(x,y)*
     &                              (N_m(x)-tmp_a)*
     &                              (N_n(y)-tmp_b)  
                                ENDDO ! y

                                inc_var_gap=inc_var_gap-
     &                            D_m(x)*N_n(x)-N_m(x)*D_n(x)
                              ENDDO ! x
                            ENDIF ! non-host mesh

                            K=inc_t_n+inc_var_gap*Z_CONT(j,1,18)

C*** 29/09/08 JHC Apply elastic tangential resistance
                            DO x=1,2
                              DO y=1,2
                                DO z=1,2
                                  DOT=0.0d0
                                  DO nh=1,3
                                    DOT=DOT+
     &                               tang_deri_comp(nh,x,y)*
     &                               tang_comp(nh,1,z)
                                  ENDDO
                                  omega(x,y,z)=DOT
                                ENDDO
                              ENDDO
                            ENDDO
                              
                            DO x=1,2
                              DO y=1,2
                                IF (m.EQ.1) THEN
                                  Ttang_m(x,y)=0.0d0
                                  Ntang_m(x,y)=0.0d0
                                ELSE
                                  Ttang_m(x,y)=basem_deri(y)*
     &                              scalem*tang_comp(nh_y,m,x)
                                  Ntang_m(x,y)=basem_double_deri(x,y)*
     &                              scalem*normal_comp(nh_y,m)
                                ENDIF

                                IF (n.EQ.1) THEN
                                  Ttang_n(x,y)=0.0d0
                                  Ntang_n(x,y)=0.0d0
                                ELSE
                                  Ttang_n(x,y)=basen_deri(y)*
     &                              scalen*tang_comp(nh_z,n,x)
                                  Ntang_n(x,y)=basen_double_deri(x,y)*
     &                              scalen*normal_comp(nh_z,n)
                                ENDIF
                              ENDDO   
                            ENDDO

                            DO x=1,2
                              K_ct(x)=0.0d0
                              tmp_a=0.0d0
                              tmp_b=0.0d0
                              DO y=1,2
                                K_ct(x)=K_ct(x)+
     &                            Ttang_m(y,x)*D_n(y)+
     &                            D_m(y)*Ttang_n(y,x)+
     &                            Ttang_m(x,y)*D_n(y)+
     &                            D_m(y)*Ttang_n(x,y)+
     &                            norm_gap*
     &                            (Ntang_m(x,y)*D_n(y)+
     &                            D_m(y)*Ntang_n(x,y))
                                tmp_a=tmp_a+kappa(x,y)*D_m(y)
                                tmp_b=tmp_b+kappa(x,y)*D_n(y)
                                DO z=1,2
                                  K_ct(x)=K_ct(x)-
     &                              omega(y,z,x)*D_m(y)*D_n(z)-
     &                              omega(x,z,y)*
     &                              (D_m(y)*D_n(z)+D_m(z)*D_n(y))
                                  tmp_c=0.0d0
                                  tmp_d=0.0d0
                                  DO v=1,2
                                    tmp_c=tmp_c+omega(v,z,z)*D_m(v)
                                    tmp_d=tmp_d+omega(v,z,z)*D_n(v)
                                  ENDDO
                                  K_ct(x)=K_ct(x)-inv_m(y,z)*T_m(y)*
     &                              (Ttang_n(z,x)-tmp_d)-inv_m(y,z)*
     &                              (Ttang_m(z,x)-tmp_c)*T_n(y) 
                                ENDDO
                              ENDDO
                              K_ct(x)=K_ct(x)-normal_comp(nh_y,m)*
     &                          scalem*basem*(N_n(x)-tmp_b)-
     &                          (N_m(x)-tmp_a)*normal_comp(nh_z,n)*
     &                          scalen*basen
                            ENDDO

C*** 29/09/08 JHC Stick-only
                            K_direct=0.0d0
                            DO x=1,2
                              G_n(x)=0.0d0
                              DO y=1,2
                                G_n(x)=G_n(x)-slip(y)*(
     &                            Ttang_n(x,y)+Ttang_n(y,x))
                                DO z=1,2
                                  G_n(x)=G_n(x)+slip(y)*
     &                              (omega(x,z,y)+omega(y,z,x))*
     &                              D_n(z)
                                ENDDO ! z
                              ENDDO ! y
                            ENDDO ! x
                            DO x=1,2
                              inc_t_trial(x)=G_n(x)
                              DO y=1,2
                                inc_t_trial(x)=inc_t_trial(x)+
     &                            (cov_m(x,y)*D_n(y))
                              ENDDO
                              inc_t_trial(x)=TIED_STIFF*
     &                          inc_t_trial(x)
                            ENDDO
                            DO x=1,2
                              K_direct=K_direct+
     &                          D_m(x)*inc_t_trial(x)
                            ENDDO

                            K=K+K_direct+( 
     &                        (Z_CONT(j,1,20)*inv_A(1,1)+              ! (Z_CONT(j,1,20) is t_1
     &                        Z_CONT(j,1,22)*inv_A(2,1))*K_ct(1)+      ! (Z_CONT(j,1,22) is t_2
     &                        (Z_CONT(j,1,20)*inv_A(1,2)+
     &                        Z_CONT(j,1,22)*inv_A(2,2))*K_ct(2))

C*** 18/03/08 JHC Removed Justin's frictional contact and replaced with Coulomb's frictional contact
C                          ELSEIF (STATUS.EQ.3) THEN ! contact with friction   
C                            IF (norm_gap.LE.0.0d0) THEN !outside volume
C                              K=((CONT_STIFF*
C     '                          exp((CONT_STIFF/F_hat)*norm_gap))*
C     '                          normal_comp(nh_y,m)*normal_comp(nh_z,n)+
C     '                          (CONT_STIFF*
C     '                          exp((CONT_STIFF/F_hat)*norm_gap))*
C     '                          FRIC_COEFF*tang_comp(nh_y,m,1)*
C     '                          tang_comp(nh_z,n,1)+
C     '                          (CONT_STIFF*
C     '                          exp((CONT_STIFF/F_hat)*norm_gap))*
C     '                          FRIC_COEFF*tang_comp(nh_y,m,2)*
C     '                          tang_comp(nh_z,n,2))*
C     '                          basem*basen*scalem*scalen    
C                            ELSE !inside volume
C                              K=(CONT_STIFF*
C     '                          normal_comp(nh_y,m)*normal_comp(nh_z,n)+
C     '                          CONT_STIFF*FRIC_COEFF*
C     '                          tang_comp(nh_y,m,1)*tang_comp(nh_z,n,1)+
C     '                          CONT_STIFF*FRIC_COEFF*
C     '                          tang_comp(nh_y,m,2)*
C     '                          tang_comp(nh_z,n,2))*basem*basen*scalem*
C     '                          scalen
C                            ENDIF
                          ENDIF ! frictionless, tied or frictional

C*** 27/02/08 JHC Removed if-statement below. 
C                          nyrfix=NYNP(nkr,nvr,nh_y,npr,0,1,nr_cont(m))
C                          nycfix=NYNP(nkc,nvc,nh_z,npc,0,1,nr_cont(n))
C                          IF((.NOT.FIX(nyrfix,1)).AND.
C     '                      (.NOT.FIX(nycfix,1))) THEN                                        
                          CALL SPARSE(nyr,nyc,NYT(1,1,nx),nz,NZ_GK_M,
     '                      NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,
     '                      *9999) 
C*** 25/02/08 JHC removed variable stiffness_comp and use K instead
C                          stiffness_comp=K*basem*basen*scalem*scalen
C                            GK(nz)=GK(nz)+stiffness_comp*weight*jacobian
                          GK(nz)=GK(nz)+K*weight*jacobian
C                          ENDIF
                        ENDDO !nnkc
                      ENDDO !nnc
                    ENDDO !nnkr
                  ENDDO !nnr
                ENDDO !nh_z
              ENDDO !nh_y
            ENDDO !n
          ENDDO !m
        ENDIF ! Z_CONT_LIST(nd,1,3).EQ.1
      ENDDO !j - contact point
      
      CALL EXITS('CONTACT_STIFFNESS')
      RETURN
 9999 CALL ERRORS('CONTACT_STIFFNESS',ERROR)
      CALL EXITS('CONTACT_STIFFNESS')
      RETURN 1
      END

