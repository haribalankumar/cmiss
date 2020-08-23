      SUBROUTINE PERFUSION_OUTPUTS(nb,NEELEM,NORD,NPNE,nx,NXI,NVJE,NYNE,
     &  NYNP,XAB,XP,YP,ERROR,*)

C#### Subroutine: PERFUSION OUTPUTS
C###  Description:
C###    PERFUSION_OUTPUTS calculates anything important that we want to 
C###    output from pulmonary perfusion models

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER nb,NEELEM(0:NE_R_M),NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),nx,
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM)
      REAL*8 XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     local variables
      INTEGER ne,ne0,nh,nj,noelem,np1,np2,nv1,nv2,sum_capillaries
      REAL*8  HEIGHT(NJT),L_in,L_out,LPM_FUNC,P1,P2,Ppl,Q01,R_in,R_out,
     &  stretch,x_cap,y_cap,z_cap,k_factor,mean_k_factor,
     &  max_shear_stress,P_inlet,vol_art,vol_vein,LENGTH,
     &  volume_art,volume_sum,avrad,
     &  volume_vein


      CALL ENTERS('PERFUSION_OUTPUTS',*9999)
      IF (LADDER.EQ.1)THEN
      ELSEIF(LADDER.EQ.2)THEN
C!!! LADDER MODEL: Open 2 files to record output. The first 
C!!! ('micro_flow_unit.out')records the properties of the whole unit.
C!!! The second ('micro_flow_ladder.out') records the propertes of each generation  
                CALL OPENF(IOFILE2,'DISK','micro_flow_unit.out','NEW',
     &             'SEQUEN','FORMATTED',132,ERROR,*9999)
                WRITE(OP_STRING,
     &           '('' ne |  x |  y |  z | Pin |'//
     &            ' Pout | Qtot |sum Qsheet | '//
     &            ' Rtot | Blood_vol |'//
     &            ' sheet_area | ave_TT |ave_H |'//
     &            ' Ppl'')')
                CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)
                CALL OPENF(IOFILE3,'DISK','micro_flow_ladder.out','NEW',
     &             'SEQUEN','FORMATTED',132,ERROR,*9999)
                WRITE(OP_STRING,
     &           '('' ne | x | y | z | gen | Pin |'//
     &            'Pout_art | Pin_vein | Pout |'//
     &            'Qtot |'//
     &            'Qgen |Qsheet | Hdiff | '//
     &            'Rsheet |Rtot |'//
     &            'RBC_TT | Hart | Hvein | zone |'//
     &            'cap_vol | cap_SA | recruited '')')
                CALL WRITES(IOFILE3,OP_STRING,ERROR,*9999)
      ELSEIF(LADDER.eq.3)THEN
       CALL OPENF(IOFILE2,'DISK','micro_flow_sheet.out','NEW',
     &                 'SEQUEN','FORMATTED',132,ERROR,*9999)
            WRITE(OP_STRING,
     &            '(''ne | Gen | S_ord | Pin | Pout |'//
     &             'Pdrop | Q |Res |'//
     &             'Cap_vol | SA=10 | Hart |'//
     &             'Hvein | zone |  RBC_tt |'//
     &             'H_diff  | x | y | z | recruited |'//
     &             'Ppl '')')
            CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)
      ENDIF
      volume_sum=0.d0           !Sum artery/vein volume
      vol_art=0.d0              !Sum extra-acinar arterial volume
      vol_vein=0.d0             !Sum extra-acinar venous volume
      max_shear_stress=0.d0     !Initialise maximum shear stress value
      P_inlet=YP(NYNP(1,1,1,np_in,0,1),1) !Pressure at inlet node
      mean_k_factor=0.d0
      sum_capillaries=0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nh=NH_LOC(1,nx)
        IF(XAB(nej_cap,ne).EQ.1.d0)THEN !check its a capillary
          sum_capillaries=sum_capillaries+1
          ne0=NXI(-1,1,ne)!upstream element number
          P1=YP(NYNP(1,1,nh,NPNE(1,nb,ne),0,1),1) !pressure at start node of capillary element
          P2=YP(NYNP(1,1,nh,NPNE(2,nb,ne),0,1),1)!pressure at end node of capillary element
          Q01=YP(NYNE(1,nh,1,1,ne0),1) !flow in element upstream of capillary element
          IF(COMPLIANCE.EQ.3) THEN
            Ppl=XP(1,1,7,NPNE(1,nb,ne))*98.06d0 !cm H2O -!Pa temporarily hard coded! 15.02.08
            stretch=1.d0 !No mechanics therefore stretch factor=1.d0
          ELSEIF(COMPLIANCE.EQ.6) THEN !Coupled to mechanics pleural & vessel stretch fields
            Ppl=XP(1,1,nj_pleural,NPNE(1,nb,ne))*98.07d0 !cmH2O->Pa
            IF(nj_stretch.EQ.0) THEN
              stretch=1.d0 !Set to 1 i.e. no effect if undefined
            ELSE
              stretch=XP(1,1,nj_stretch,NPNE(1,nb,ne))
            ENDIF
          ENDIF
          R_in=XP(1,1,nj_radius_R0,NPNE(1,nb,ne))*
     &      sqrt(1.d0/stretch) !Including mechanics stretch factor (as is CALC_PRESS_AREA)
          R_out=XP(1,1,nj_radius_R0,NPNE(2,nb,ne))*
     &      sqrt(1.d0/stretch)
C... For hypoxia simulations which include a vascular constriction factor:
           IF(nj_hypoxia.NE.0) THEN
            k_factor=XP(1,1,nj_hypoxia,NPNE(1,nb,ne)) !Arterial constriction factor
            IF(k_factor.GT.1.d0) THEN
             !Use WRITE(OP_STRING...) and CALL WRITES for this
              write(*,*) "CHANGE CODE k_factor",NPNE(1,nb,ne),k_factor
             !k_factor=1.d0 !If undefined (=2) then have equal to 1
            ENDIF
          ELSE
            k_factor=1.d0      !If undefined Then have equal to 1 (no constriction)
          ENDIF
          mean_k_factor=mean_k_factor+k_factor

          !Length of up and downstream elts   
          np1=NPNE(1,nb,ne0)
          np2=NPNE(2,nb,ne0)
          L_in=((XP(1,1,1,np2)-XP(1,1,1,np1))**2.d0+
     &     (XP(1,1,2,np2)-XP(1,1,2,np1))**2.d0+
     &     (XP(1,1,3,np2)-XP(1,1,3,np1))**2.d0)**0.5d0
          L_out=L_in !Trees are currently always the same at terminals
          x_cap=XP(1,1,1,NPNE(1,nb,ne))
          y_cap=XP(1,1,2,NPNE(1,nb,ne))
          z_cap=XP(1,1,3,NPNE(1,nb,ne))
          DO nj=1,NJT
            HEIGHT(nj)=XP(1,1,nj,NPNE(1,nb,ne))-
     &       XP(1,1,nj,np_in) !gravitational head - same for art & vein now!
          ENDDO
          CALL CAP_FLOW_PARAM(ne,L_in,L_out,Ppl,
     &      R_in,R_out,stretch,ERROR,*9999)
          IF(LADDER.EQ.1)THEN
          ELSEIF(LADDER.EQ.2)THEN
           CALL CAP_FLOW_LADDER(ne,HEIGHT,LPM_FUNC,P1,P2,Ppl,Q01,R_in,
     &            R_out,x_cap,y_cap,z_cap,k_factor,.TRUE.,ERROR,*9999)
          ELSEIF(LADDER.EQ.3)THEN
            CALL CAP_FLOW_MBA(ne,NORD,LPM_FUNC,P1,P2,Ppl,Q01,      
     &        x_cap,y_cap,z_cap,.TRUE.,ERROR,*9999)
          ENDIF!LADDER
      ELSEIF(XAB(nej_cap,ne).EQ.0.d0.OR.XAB(nej_cap,ne).EQ.2.d0)THEN !Artery or vein
         np1=NPNE(1,nb,ne)
         np2=NPNE(2,nb,ne)
         LENGTH=((XP(1,1,1,np2)-XP(1,1,1,np1))**2.d0+
     &        (XP(1,1,2,np2)-XP(1,1,2,np1))**2.d0+
     &        (XP(1,1,3,np2)-XP(1,1,3,np1))**2.d0)**0.5d0
         nv1=NVJE(1,nb,nj_radius,ne)
         nv2=NVJE(2,nb,nj_radius,ne)            
         avrad=0.5d0*(XP(1,nv1,nj_radius,np1)+
     &        XP(1,nv2,nj_radius,np2))
         volume_sum=volume_sum+(pi*avrad**2.d0*length)/
     &        1000.d0           !ml
         Q01=YP(NYNE(1,nh,1,1,ne),1) !flow in element upstream of capillary element
c         IF(XAB(nej_cap,ne).EQ.2.d0) THEN   !Venous, this allocated and -ve direction
         IF(Q01.LT.0.d0) THEN   !Venous, this allocated and -ve direction
            vol_vein=vol_vein+(pi*avrad**2.d0*length)/
     &           1000.d0        !ml
c         ELSEIF(XAB(nej_cap,ne).EQ.0.d0) THEN
         ELSE
            vol_art=vol_art+(pi*avrad**2.d0*length)/
     &           1000.d0        !ml
C         ELSE
C            WRITE(*,*) "----- ERROR ---- PERFUSION_OUTPUTS.f this 
C     &         element has no vessel type!",ne
         ENDIF                  !Not a capillary element
         IF(XAB(7,ne).GT.max_shear_stress) THEN
            max_shear_stress=XAB(7,ne)
         ENDIF         
         mean_k_factor=mean_k_factor/sum_capillaries

        ENDIF !XAB art cap ven
       ENDDO
      IF (LADDER.EQ.1)THEN
      ELSEIF(LADDER.EQ.2)THEN
       CALL CLOSEF(IOFILE2,ERROR,*9999)
       CALL CLOSEF(IOFILE3,ERROR,*9999)
      ELSEIF(LADDER.eq.3)THEN
       CALL CLOSEF(IOFILE2,ERROR,*9999)
      ENDIF

C... KSB: Writing out blood volumes etc 18/02/11
         !Use WRITE(OP_STRING...) and CALL WRITES for this
         write(*,*) "CHANGE CODE ARTERY/VEIN VOLUME (ml)=",volume_sum 
         volume_art=volume_art/1000.d0 !Intra-acinar arteries
         volume_vein=volume_vein/1000.d0 !Intra-acinar veins
         
C...  Write volume data into a file for output:
         CALL OPENF(IOFILE3,'DISK','blood_volumes.out','NEW',
     &        'SEQUEN','FORMATTED',132,ERROR,*9999)
         WRITE(OP_STRING,'('' TOTAL | Artery | Vein | '//
     &        'Arteriole | Venule | Sheet | Below Sat (ml) | '//
     &        'Max shear stress'')')
         CALL WRITES(IOFILE3,OP_STRING,ERROR,*9999)
         WRITE(IOFILE3,'(6(F9.3,X),2(F9.5,X))') volume_sum,
     &        vol_art,vol_vein,max_shear_stress
         CALL CLOSEF(IOFILE3,ERROR,*9999)
                  


      CALL EXITS('PERFUSION_OUTPUTS')
      RETURN
 9999 CALL ERRORS('PERFUSION_OUTPUTS',ERROR)
      CALL EXITS('PERFUSION_OUTPUTS')
      RETURN 1
      END




