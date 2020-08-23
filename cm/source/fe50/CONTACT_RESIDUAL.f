      SUBROUTINE CONTACT_RESIDUAL(IBT,IDO,INP,NBH,NBJ,NBJF,NFF,NHE,
     '  NKEF,NKHE,NKJE,NNF,NPF,NPNE,nr_loop,NRE,NVHE,NVJE,NW,nx,NYNP,
     '  Z_CONT_LIST,CURVCORRECT,SE,YP,XA,XP,Z_CONT,ZA,ZE,ZP,FIX,
     '  ERROR,*)


C#### Subroutine: CONTACT_RESIDUAL
C###  Description:
C###    Calculates contact residual component and modifies global residual
C###    vector YP(ny,4) with contact contribution.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nonl00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NFF(6,NEM),NHE(NEM),
     '  NKEF(0:4,16,6,NBFM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),nr_loop,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),Z_CONT_LIST(NDM,2,7)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  YP(NYM,NIYM),XA(NAM,NJM,NEM),XP(NKM,NVM,NJM,NPM),
     '  Z_CONT(NDM,2,67),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)

!     Local Variables
      REAL*8 base,base_deri_1,base_deri_2,
     &  contra_tang1(3),contra_tang2(3),dxj_dxi(3,3),inv_A(2,2),
     &  inv_m(2,2),jacobian,jac1,jac2,jac3,lamda_n,lamda_t(2),
     &  m(2,2),magnitude,normal_comp(3),norm_gap,
     &  Phi_trial,residual,residual_friction,slip(2),sum,
     &  t(2),tang1_comp(3),tang2_comp(3),tang1_gap,tang2_gap,
     &  t_trial(2),tb_trial(3),weight,XI(3),XI_new(2)
      INTEGER nk,nv,nh,np,nb,ny,nn,i,j,ne,nr,row,col,nu,nsf,nse,
     '  nf,nef,nnk,STATUS

!     Functions
      REAL*8 PSI1

!     Variables setup in CALC_FACE_INFORMATION_DEP
      REAL*8 SF(NSM,NBFM),XE(NSM,NJM)
      INTEGER NKJF(NKM,NNM,NJM),NPNF(NNM,NBFM),NVJF(NNM,NBFM,NJM)

      CALL ENTERS('CONTACT_RESIDUAL',*9999)

C*** 22/02/08 JHC CONT_WT removed from common block nonl00.cmn
C      weight=1.0d0/DBLE(CONT_WT**2) ! integration point weight

      DO i=1,2  ! Loop over Slave/Master
        DO j=1,NDT   ! Loop over contact points
C*** 22/02/08 JHC allow different number of contact points at different face
          weight=1.0d0/DBLE(Z_CONT_LIST(j,1,6)) ! integration point weight

          IF((Z_CONT_LIST(j,1,3).EQ.1).OR.  ! in contact
     '      (Z_CONT_LIST(j,1,4).EQ.2)) THEN  ! tied contact
            ne=Z_CONT_LIST(j,i,1) ! global elem # for contact pt k
            nr=NRE(ne) ! region # for elem ne
            IF(KTYP5H(nr).EQ.0) THEN ! do face projection
              nef=Z_CONT_LIST(j,i,2) ! local face # for contact pt k
              nf=NFF(nef,ne) ! global face #
            ENDIF
C*** 28/02/08 JHC Added for checking new contact points and initialise their contact forces            
            IF((CONT_IT.EQ.1).OR. ! 1st contact step
     &        (Z_CONT_LIST(j,i,5).EQ.1)) THEN ! New contact points added so need to initialise its Lagrang. multipliers
C*** 10/10/08 JHC DO NOT initialise contact forces when they are read in through ipdata file
              IF(.NOT.CALL_DATA_CONT) THEN
                Z_CONT(j,i,17)=0.0d0 ! Initialise normal lagrangian multiplier
                Z_CONT(j,i,19)=0.0d0 ! Initialise tangent1 lagrangian multiplier
                Z_CONT(j,i,21)=0.0d0 ! Initialise tangent2 lagrangian multiplier
              ENDIF

              IF(Z_CONT_LIST(j,i,7).NE.0) THEN ! New contact points coincide geometrically with other exsisting contact point
                Z_CONT(j,i,17)=Z_CONT(Z_CONT_LIST(j,i,7),i,17) ! copy normal lagrangian multiplier from the identical exising contact point
                Z_CONT(j,i,19)=Z_CONT(Z_CONT_LIST(j,i,7),i,19) ! copy tangent1 lagrangian multiplier
                Z_CONT(j,i,21)=Z_CONT(Z_CONT_LIST(j,i,7),i,21) ! copy tangent2 lagrangian multiplier
              ENDIF
              Z_CONT_LIST(j,i,5)=2 ! new contact point's frictional forces are updated, switch the flag off
            ENDIF

            lamda_n=Z_CONT(j,i,17) ! Current estimate of normal contact force.
            lamda_t(1)=Z_CONT(j,i,19) ! Current estimate of tangential force 1
            lamda_t(2)=Z_CONT(j,i,21) ! Current estimate of tangential force 2
                                                   
            STATUS=Z_CONT_LIST(j,i,4) ! frictionless, tied contact or friction

            IF(nr.EQ.nr_loop) THEN ! only do nr # from ZPRP loop

              IF(KTYP5H(nr).EQ.1) THEN ! adding residual to element host

                CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '            NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '            ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)

              ELSE ! ! adding residual to face 

C*** 03/03/08 JHC Integral over contact surface must be done in the reference configuration, XE
C                CALL CALC_FACE_INFORMATION_DEP(NBH(1,1,ne),
C     '            NBHF(1,1,nf),nef,NHE(ne),
C     '            NKHE(1,1,1,ne),NKEF,NKHF,NNF,NPNE(1,1,ne),
C     '            NPNF,NVHE(1,1,1,ne),NVHF,nx,SE(1,1,ne),SF,
C     '            ERROR,*9999)

C                CALL ZPZE(NBHF(1,1,nf),1,3,
C     '            NKHF,NPF(1,nf),NPNF,nr,NVHF,NW(ne,1),nx,
C     '            CURVCORRECT(1,1,1,ne),SF,ZA(1,1,1,ne),
C     '            ZE,ZP,ERROR,*9999)
                CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
     '            NBJF(1,nf),nef,NKJE(1,1,1,ne),
     '            NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
     '            NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)

                CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,NVJF,
     '            SF,XA(1,1,ne),XE,XP,ERROR,*9999)

              ENDIF !KTYP5H 

              DO nh=1,3   ! Loop over Z_CONT fields
                XI(nh)=Z_CONT(j,i,nh)
                normal_comp(nh)=Z_CONT(j,i,nh+3)
                tang1_comp(nh)=Z_CONT(j,i,nh+6)
                tang2_comp(nh)=Z_CONT(j,i,nh+9)
              ENDDO !nh
                                
              norm_gap=Z_CONT(j,i,13)
              tang1_gap=Z_CONT(j,i,14)
              tang2_gap=Z_CONT(j,i,15)

              IF((STATUS.EQ.1).OR.(STATUS.EQ.3)) THEN ! frictionless/frictional contact
C*** 28/02/08 JHC Added penalty method
                IF(PENALTY.EQ.1) THEN ! penalty method
C It should be noted that AUGMENT must be used here rather than
C AUG_IT as AUG_IT gets re-initialised to be zero on next LOAD step
C which then would use penalty method rather than augmented Lag. method
C Same for CONTACT_STIFFNESS.f
                  IF(AUGMENT.EQ.0) THEN ! standard penalty method
                    IF(norm_gap.LT.0.0d0) THEN !outside volume
                      Z_CONT(j,i,18)=0.0d0 
                    ELSE !inside volume
                      Z_CONT(j,i,18)=CONT_STIFF*norm_gap
                    ENDIF !norm_gap
                  ELSE ! Augmented Lagrangian method
                    Z_CONT(j,i,18)=lamda_n+CONT_STIFF*norm_gap
                    Z_CONT(j,i,46)=Z_CONT(j,i,18)
                    IF(Z_CONT(j,i,46).LT.0.0d0) THEN !outside volume
                      Z_CONT(j,i,18)=0.0d0 
                    ENDIF !norm_gap
                  ENDIF 
                ELSE ! X-constraint method
                  IF(norm_gap.LT.0.0d0) THEN !outside volume
C*** 26/02/08 JHC Added check for division by zero
                    IF(DABS(lamda_n).LE.ZERO_TOL) THEN
                      Z_CONT(j,i,18)=0.0d0
                    ELSE 
                      Z_CONT(j,i,18)=(lamda_n*exp((CONT_STIFF/lamda_n)*
     '                  norm_gap))!Store normal contact force 
                    ENDIF 
                  ELSE !inside volume
                    Z_CONT(j,i,18)=(lamda_n+(CONT_STIFF*norm_gap))!Store normal contact force 
                  ENDIF !norm_gap
                ENDIF ! Penalty or X-constraint

C*** 18/03/08 JHC Added Coulomb's frictional contact
                IF(STATUS.EQ.3) THEN ! frictional contact
C*** Not implemented (at least not tested) with host-mesh problems
                  IF(KTYP5H(nr).EQ.1) THEN ! host element basis
                    CALL ASSERT(.FALSE.,
     '                '>>Frictional contact not tested for host meshes',
     '                ERROR,*9999)
                  ENDIF
                  slip(1)=Z_CONT(j,i,43)
                  slip(2)=Z_CONT(j,i,44)
 
                  m(1,1)=Z_CONT(j,i,47)
                  m(1,2)=Z_CONT(j,i,48)
                  m(2,1)=Z_CONT(j,i,49)
                  m(2,2)=Z_CONT(j,i,50)

                  inv_m(1,1)=Z_CONT(j,i,51)
                  inv_m(1,2)=Z_CONT(j,i,52)
                  inv_m(2,1)=Z_CONT(j,i,53)
                  inv_m(2,2)=Z_CONT(j,i,54)

                  inv_A(1,1)=Z_CONT(j,i,55)
                  inv_A(1,2)=Z_CONT(j,i,56)
                  inv_A(2,1)=Z_CONT(j,i,57)
                  inv_A(2,2)=Z_CONT(j,i,58)

                  t_trial(1)=lamda_t(1)+FRIC_STIFF*
     &              (m(1,1)*slip(1)+m(1,2)*slip(2))
                  t_trial(2)=lamda_t(2)+FRIC_STIFF*
     &              (m(2,1)*slip(1)+m(2,2)*slip(2))

                  magnitude=DSQRT(t_trial(1)*inv_m(1,1)*t_trial(1)+
     &             t_trial(1)*inv_m(1,2)*t_trial(2)+
     &             t_trial(2)*inv_m(2,1)*t_trial(1)+ 
     &             t_trial(2)*inv_m(2,2)*t_trial(2))

                  Z_CONT(j,i,59)=magnitude

                  Phi_trial=magnitude-FRIC_COEFF*Z_CONT(j,i,18)
            
                  Z_CONT(j,i,45)=Phi_trial

                  IF(Phi_trial.LE.0.0d0) THEN ! stick

                    Z_CONT(j,i,20)=t_trial(1) ! store t_1 at n+1
                    Z_CONT(j,i,22)=t_trial(2) ! store t_2 at n+1

                    Z_CONT(j,i,60)=1.0d0

C***                Set Pt as zero    
                    Z_CONT(j,i,61)=0.0d0
                    Z_CONT(j,i,62)=0.0d0
                    Z_CONT(j,i,63)=0.0d0

C***                Set Pt_1,2 covariant contravariant as zero    
                    Z_CONT(j,i,64)=0.0d0
                    Z_CONT(j,i,65)=0.0d0

C***                Set Pt^1,2 contravariant as zero    
                    Z_CONT(j,i,66)=0.0d0
                    Z_CONT(j,i,67)=0.0d0

                 ELSE                         ! slip

                    IF((DABS(magnitude).LT.ZERO_TOL).AND.i.EQ.1) THEN
                      WRITE(OP_STRING,'('' >>Warning: magnitude of '
     &                  //'frictional forces at nd='',I5,'
     &                  //''' is less than zero => '
     &                  //'Division By Zero'')') j
                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                    ENDIF

                    Z_CONT(j,i,20)=FRIC_COEFF*Z_CONT(j,i,18)*t_trial(1) ! store t_1 at n+1
     &                /magnitude
                    Z_CONT(j,i,22)=FRIC_COEFF*Z_CONT(j,i,18)*t_trial(2) ! store t_2 at n+1
     &                /magnitude

                    Z_CONT(j,i,60)=-1.0d0
                 
                    DO nh=1,3
                      contra_tang1(nh)=inv_m(1,1)*Z_CONT(j,1,nh+6)+ ! tangent1
     &                  inv_m(1,2)*Z_CONT(j,1,nh+9)                 ! tangent2
                      contra_tang2(nh)=inv_m(2,1)*Z_CONT(j,1,nh+6)+
     &                  inv_m(2,2)*Z_CONT(j,1,nh+9)
                      tb_trial(nh)=t_trial(1)*contra_tang1(nh)+
     &                  t_trial(2)*contra_tang2(nh)
                    ENDDO

C***                store Pt
                    Z_CONT(j,i,61)=tb_trial(1)/magnitude
                    Z_CONT(j,i,62)=tb_trial(2)/magnitude
                    Z_CONT(j,i,63)=tb_trial(3)/magnitude
  
C***                compute and store Pt_1,2 covariant
                    Z_CONT(j,i,64)=t_trial(1)/magnitude
                    Z_CONT(j,i,65)=t_trial(2)/magnitude

C***                compute and store Pt^1,2 contravariant
                    Z_CONT(j,i,66)=(inv_m(1,1)*t_trial(1)+
     &                inv_m(1,2)*t_trial(2))/magnitude
                    Z_CONT(j,i,67)=(inv_m(2,1)*t_trial(1)+
     &                inv_m(2,2)*t_trial(2))/magnitude

                  ENDIF

                ENDIF

              ELSEIF(STATUS.EQ.2) THEN ! tied contact
C*** 29/09/08 JHC Replaced Justin's tied contact with modified frictional contact (stick only).
C                Z_CONT(j,i,18)=(lamda_n+(TIED_STIFF*norm_gap))!Store normal contact force
C                Z_CONT(j,i,20)=(lamda_t(1)+(TIED_STIFF*tang1_gap))!Store tan1 contact force
C                Z_CONT(j,i,22)=(lamda_t(2)+(TIED_STIFF*tang2_gap))!Store tan2 contact force  

C*** 29/09/08 JHC Modify tied contact to apply compression/tension
                IF(AUGMENT.EQ.0) THEN ! Regardless of X-constraint or penalty method.
                  Z_CONT(j,i,18)=CONT_STIFF*norm_gap
                ELSE ! Augmented Lagrangian method
                  Z_CONT(j,i,18)=lamda_n+CONT_STIFF*norm_gap
                ENDIF 

C*** 29/09/08 JHC Apply elasic resistance: only stick tangential force
                slip(1)=Z_CONT(j,i,43)
                slip(2)=Z_CONT(j,i,44)
 
                m(1,1)=Z_CONT(j,i,47)
                m(1,2)=Z_CONT(j,i,48)
                m(2,1)=Z_CONT(j,i,49)
                m(2,2)=Z_CONT(j,i,50)

                inv_m(1,1)=Z_CONT(j,i,51)
                inv_m(1,2)=Z_CONT(j,i,52)
                inv_m(2,1)=Z_CONT(j,i,53)
                inv_m(2,2)=Z_CONT(j,i,54)

                inv_A(1,1)=Z_CONT(j,i,55)
                inv_A(1,2)=Z_CONT(j,i,56)
                inv_A(2,1)=Z_CONT(j,i,57)
                inv_A(2,2)=Z_CONT(j,i,58)

                t_trial(1)=lamda_t(1)+TIED_STIFF*
     &            (m(1,1)*slip(1)+m(1,2)*slip(2))
                t_trial(2)=lamda_t(2)+TIED_STIFF*
     &            (m(2,1)*slip(1)+m(2,2)*slip(2))

                magnitude=DSQRT(t_trial(1)*inv_m(1,1)*t_trial(1)+
     &           t_trial(1)*inv_m(1,2)*t_trial(2)+
     &           t_trial(2)*inv_m(2,1)*t_trial(1)+ 
     &           t_trial(2)*inv_m(2,2)*t_trial(2))

                Z_CONT(j,i,59)=magnitude

                Z_CONT(j,i,20)=t_trial(1) ! store t_1 at n+1
                Z_CONT(j,i,22)=t_trial(2) ! store t_2 at n+1
                Z_CONT(j,i,60)=1.0d0

C***            Set Pt as zero    
                Z_CONT(j,i,61)=0.0d0
                Z_CONT(j,i,62)=0.0d0
                Z_CONT(j,i,63)=0.0d0

C***            Set Pt_1,2 covariant contravariant as zero    
                Z_CONT(j,i,64)=0.0d0
                Z_CONT(j,i,65)=0.0d0

C***            Set Pt^1,2 contravariant as zero    
                Z_CONT(j,i,66)=0.0d0
                Z_CONT(j,i,67)=0.0d0

C*** 18/03/08 JHC Removed Justin's frictional contact and 
C                 replaced with Coulomb's frictional contact
C              ELSEIF (STATUS.EQ.3) THEN !including friction   
C                IF(norm_gap.LE.0.0) THEN !outside volume 
C                  Z_CONT(j,i,18)=(F_hat*exp((CONT_STIFF/F_hat)*
C     '            norm_gap))!Store normal contact force  
C                  Z_CONT(j,i,20)=Z_CONT(j,i,18)*FRIC_COEFF*tang1_gap !Store tan1 contact force 
C                  Z_CONT(j,i,22)=Z_CONT(j,i,18)*FRIC_COEFF*tang2_gap !Store tan2 contact force        
C                ELSE !inside volume
C                  Z_CONT(j,i,18)=(F_hat+(CONT_STIFF*norm_gap))!Store normal contact force  
C                  Z_CONT(j,i,20)=Z_CONT(j,i,18)*FRIC_COEFF*tang1_gap !Store tan1 contact force 
C                  Z_CONT(j,i,22)=Z_CONT(j,i,18)*FRIC_COEFF*tang2_gap !Store tan2 contact force        
C                ENDIF !norm_gap  
              ENDIF !STATUS

              DO nh=1,3   ! loop over x,y,z
                IF(KTYP5H(nr).EQ.1) THEN ! host element basis
                  nb=NBH(nh,1,ne)  ! bases function #
                ELSE ! face basis
C                  nb=NBHF(nh,1,nf)  ! bases function #
                  nb=NBJF(nh,nf)  ! bases function #
                ENDIF

                IF(KTYP5H(nr).EQ.0) THEN ! 2D surface
                  IF(i.EQ.1) THEN ! slave surface
C ***               Calculate Jacobian for slave surface only
                    DO row=1,2 !
                      DO col=1,3 ! x,y,z
                        sum=0.0d0
                        DO nn=1,NNT(nb)
                          DO nk=1,NKT(nn,nb)
                            nsf=nk+(nn-1)*NKT(nn,nb)
                            IF(row.EQ.1) THEN
                              nu=2 ! du/dXi1
                            ELSEIF(row.EQ.2) THEN
                              nu=4 ! du/dXi2
                            ENDIF
                            XI_new(1)=XI(NPF(1,nf))
                            XI_new(2)=XI(NPF(3,nf))
C*** 03/03/08 JHC Integration must be done over undeformed surface
C                            sum=sum+ZE(nsf,col)*
C     '                        PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI_new)
C*** 01/09/2008 JHC should be using correct basis function number
C                            sum=sum+XE(nsf,col)*
C     '                        PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI_new)
                            sum=sum+XE(nsf,col)*PSI1(
     &                        IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     &                        nb,nu,nk,nn,XI_new)
                          ENDDO ! nk
                          dxj_dxi(row,col)=sum
                        ENDDO !nn
                      ENDDO !col
                    ENDDO !row
                    jac1=(dxj_dxi(2,3)*dxj_dxi(1,2))-
     '                (dxj_dxi(2,2)*dxj_dxi(1,3))
                    jac2=(dxj_dxi(2,3)*dxj_dxi(1,1))-
     '                (dxj_dxi(2,1)*dxj_dxi(1,3))
                    jac3=(dxj_dxi(2,2)*dxj_dxi(1,1))-
     '                (dxj_dxi(2,1)*dxj_dxi(1,2))
                    jacobian=dsqrt(jac1**2+jac2**2+jac3**2)

                    Z_CONT(j,1,16)=jacobian ! store jacobian in slave
                    Z_CONT(j,2,16)=jacobian ! store jacobian in master
                  ENDIF !i
                ELSE !KTYP5H(nr).EQ.1 (3D volume host)
                  Z_CONT(j,1,16)=1.0d0 ! store jacobian in slave
                  Z_CONT(j,2,16)=1.0d0 ! store jacobian in master
                ENDIF !KTYP5H
                                                                       
                jacobian=Z_CONT(j,i,16)

                DO nn=1,NNT(nb)  ! weight contact pt. out to global nodes
                  IF(KTYP5H(nr).EQ.1) THEN
                    nv=NVHE(nn,nb,nh,ne)  ! version #
                    np=NPNE(nn,nb,ne)   ! global node #
                  ELSE
C                    nv=NVHF(nn,nb,nh)  ! version #
                    nv=NVJF(nn,nb,nh)  ! version #
                    np=NPNF(nn,nb)   ! global node #
                  ENDIF
                  DO nnk=1,NKT(nn,nb)  ! loop over nk
                    IF(KTYP5H(nr).EQ.1) THEN
                      nk=NKHE(nnk,nn,nh,ne)
                    ELSE
C                      nk=NKHF(nnk,nn,nh)
                      nk=NKJF(nnk,nn,nh)
                    ENDIF
                    IF(KTYP5H(nr).eq.1) THEN  
C*** 01/09/2008 JHC should be using correct basis function number
C                      base=PSI1(IBT,IDO,INP,nb,1,nnk,nn,XI)
                      base=PSI1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     &                  nb,1,nnk,nn,XI)
                    ELSE
                      XI_new(1)=XI(NPF(1,nf))
                      XI_new(2)=XI(NPF(3,nf))
C*** 01/09/2008 JHC should be using correct basis function number
C                      base=PSI1(IBT,IDO,INP,nb,1,nnk,nn,XI_new)
                      base=PSI1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     &                  nb,1,nnk,nn,XI_new)
C*** 18/03/08 JHC for frictional contact
                      IF(i.EQ.1) THEN
                        base_deri_1=0.0d0
                        base_deri_2=0.0d0
                      ELSE
C*** 01/09/2008 JHC should be using correct basis function number
C                        base_deri_1=PSI1(IBT,IDO,INP,nb,2,nnk,nn,XI_new)
C                        base_deri_2=PSI1(IBT,IDO,INP,nb,4,nnk,nn,XI_new)
                        base_deri_1=PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
     &                    INP(1,1,nb),nb,2,nnk,nn,XI_new)
                        base_deri_2=PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
     &                    INP(1,1,nb),nb,4,nnk,nn,XI_new)
                      ENDIF
                    ENDIF

C*** 28/02/08 JHC Compute contact residual
C*** 29/09/08 JHC No longer need IF statement. Tied contact follows the same eq. as frictional contact.
C                    IF((STATUS.EQ.1).OR.(STATUS.EQ.3)) THEN ! frictionless/frictional contact

                    residual=Z_CONT(j,i,18)*normal_comp(nh)*base ! normal part
                    IF(STATUS.EQ.3.OR.STATUS.EQ.2) THEN ! tangential part for frictional/tied contact
                      residual_friction=
     &                  Z_CONT(j,i,20)*   
     &                  (inv_A(1,1)*(tang1_comp(nh)*base+
     &                   norm_gap*base_deri_1*normal_comp(nh))+
     &                   inv_A(1,2)*(tang2_comp(nh)*base+
     &                   norm_gap*base_deri_2*normal_comp(nh)))+
     &                  Z_CONT(j,i,22)*
     &                  (inv_A(2,1)*(tang1_comp(nh)*base+
     &                   norm_gap*base_deri_1*normal_comp(nh))+
     &                   inv_A(2,2)*(tang2_comp(nh)*base+
     &                   norm_gap*base_deri_2*normal_comp(nh)))
                      residual=residual-residual_friction
                    ENDIF

C*** 29/09/08 JHC Removed Justin's tied contact. Tied contact follows the same eq. as frictional contact.
C                    ELSEIF(STATUS.EQ.2) THEN ! tied contact 
C                      residual=(Z_CONT(j,i,18)*normal_comp(nh)+
C     &                  Z_CONT(j,i,20)*tang1_comp(nh)+
C     &                  Z_CONT(j,i,22)*tang2_comp(nh))*base
C
C*** 18/03/08 JHC removed Justin's frictional contact and replaced with Coulomb's frictional contact
C                    ELSE ! frictional
C                      residual=(Z_CONT(j,i,18)*normal_comp(nh)+
C     &                  Z_CONT(j,i,20)*tang1_comp(nh)+
C     &                  Z_CONT(j,i,22)*tang2_comp(nh))*base
C                    ENDIF

                    ny=NYNP(nk,nv,nh,np,1,1,nr)
                    IF(KTYP5H(nr).EQ.1) THEN
                      nse=nnk+(nn-1)*NKT(nn,nb)
                      YP(ny,4)=YP(ny,4)-residual*
     &                  SE(nse,nb,ne)*weight*jacobian
                    ELSE
                      nsf=nnk+(nn-1)*NKT(nn,nb)
C*** 27/02/08 JHC Removed if statement below. 
C                      nyfix=NYNP(nk,nv,nh,np,0,1,nr)
C                      IF(.NOT.FIX(nyfix,1)) THEN
                      YP(ny,4)=YP(ny,4)-residual*
     &                  SF(nsf,nb)*weight*jacobian
C                      ENDIF
                      YP(ny,6)=YP(ny,6)-residual*    ! store nodal contact residual
     &                  SF(nsf,nb)*weight*jacobian   ! can be accessed via "fem list node contact_force"
                    ENDIF
                  ENDDO !nk
                ENDDO !nn
              ENDDO !nh
            ENDIF !nr_loop
          ENDIF !Z_CONT_LIST(nd,1,3).EQ.1
        ENDDO !j
      ENDDO !i


      CALL EXITS('CONTACT_RESIDUAL')
      RETURN
 9999 CALL ERRORS('CONTACT_RESIDUAL',ERROR)
      CALL EXITS('CONTACT_RESIDUAL')
      RETURN 1
      END

