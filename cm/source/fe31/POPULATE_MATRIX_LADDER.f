      SUBROUTINE POPULATE_MATRIX_LADDER(ne,NonZeros,MatrixSize,
     & submatrixsize,SparseCol,SparseRow,area,Pin,Pout,Ppl,Pressure,
     & Q01_mthrees,Q_sheet,RHS,SparseVal,k_factor,ERROR,*)
C###  Description:
C###  Sets up ladder matrix entries that are NOT independent of iteration
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'!
      INCLUDE 'lung_nej00.cmn'

!      !INPUT AND OUTPUT PARAMETER LIST
      INTEGER ne,NonZeros,MatrixSize,submatrixsize,SparseCol(NonZeros),
     &      SparseRow(MatrixSize+1)
      REAL*8 area,Pin,Pout,Ppl,Pressure(submatrixsize),Q01_mthrees,
     &       Q_sheet(num_symm_gen),RHS(MatrixSize),SparseVal(Nonzeros),
     &       k_factor
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER gen,count,j
      REAL*8 Q,Palv_Pa,P_exta,P_extv,radupdate,R_art1,R_art2,
     &   R_ven1,R_ven2,SHEET_RES,Q_c,
     &   R_sheet(num_symm_gen),Q_gen,Hart,
     &   Hven,RBC_TT,zone,Pin_sheet,Pout_sheet,area_new,test,
     &   recruited,volume_vessels

      CALL ENTERS('POPULATE_MATRIX_LADDER',*9999)
C###  RESISTANCE CALCULATIONS STEPPING THROUGH THE GENERATIONS 
C###  AND CONSERVATION OF FLOW!!
C...  initialising vessel volumes
      volume_vessels=0.d0

C...  Previous iterations estimate for total flow through the system
      !ALYS: at the moment we aren't iterating so use Q01
       Q=Q01_mthrees  
       radupdate=0.d0
      DO gen=1,num_symm_gen-1
C...    FIRST HALF OF ARTERIOLE
C...    Update radius of arteriole based on inlet pressure 
         IF(rad_a(gen).LT.100.d-6) THEN 
           P_exta=Palv*98.06d0 ! From Yen Alveolar pressure dominates vessels <200um diam
         ELSE
           P_exta=-Ppl
         ENDIF
         IF ((Pressure(4*gen-3)-P_exta).LE.Pub_a_v)THEN
c            radupdate=rad_a(gen)*(1+(Pressure(4*gen-3)-P_exta)*alpha_a)
            radupdate=rad_a(gen)+alpha_a*(Pressure(4*gen-3)-P_exta)*
     &       (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*
     &       (Pressure(4*gen-3)-P_exta)/(2.d0-2.d0*num_symm_gen)
         ELSE
c            radupdate=rad_a(gen)*(1+Pub_a_v*alpha_a)
            radupdate=rad_a(gen)+alpha_a*Pub_a_v*(gen-num_symm_gen)
     &       /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v
     &       /(2.d0-2.d0*num_symm_gen)
         ENDIF
         IF(nj_hypoxia.NE.0) THEN
           IF(rad_a(gen).LE.0.25d0)THEN
              radupdate=radupdate*k_factor   
           ENDIF
         ENDIF
         volume_vessels=volume_vessels+
     &     (2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0)

C...  Calculate Poiseuille resistance in first half of arteriole - (only 
C...  half total generation length)
           R_art1=(8.d0*mu_app(gen)*L_a(gen)/2.d0)/(pi*radupdate**4.d0)

C...    FIRST HALF OF VENULE        
C...    Update radius of venule based on inlet pressure NB: May need to
C...     update this later to give an average of inlet and outlet radii
           IF (Pressure(4*gen-1)-Palv_pa.LE.Pub_a_v)THEN
           radupdate=rad_v(gen)*(1+(Pressure(4*gen-1)-Palv_pa)*alpha_v)
           ELSE
              radupdate=rad_v(gen)*(1+Pub_a_v*alpha_v)
           ENDIF
         IF(rad_v(gen).LT.100.d-6) THEN 
           P_extv=Palv*98.06d0 ! From Yen Alveolar pressure dominates vessels <200um diam
         ELSE
           P_extv=-Ppl
         ENDIF
         IF((Pressure(4*gen-1)-P_extv).LE.Pub_a_v)THEN
c           radupdate=rad_v(gen)*(1+(Pressure(4*gen-1)-P_extv)*alpha_v)
            radupdate=rad_v(gen)+alpha_v*(Pressure(4*gen-1)-P_extv)*
     &       (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*
     &       (Pressure(4*gen-1)-P_extv)/(2.d0-2.d0*num_symm_gen)
         ELSE
c            rad_update=rad_v(gen)*(1+Pub_a_v*alpha_v)
            radupdate=rad_v(gen)+alpha_v*Pub_a_v*(gen-num_symm_gen)
     &       /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v
     &       /(2.d0-2.d0*num_symm_gen)
         ENDIF
         volume_vessels=volume_vessels+
     &     (2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0)

C...   Calculate Poiseuille resistance in first half of venule
           R_ven1=(8.d0*mu_app(gen)*L_v(gen)/2.d0)/(pi*radupdate**4.d0);

C...   CAPILLARY ELEMENT (arteriole + venule + capillary)  
C...    pressure into the capillaries
           Pin_sheet=Pressure(4*gen-3)
           Pout_sheet=Pressure(4*gen-1)
C...     calulate resistance and transit time through a single capillary
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, 
     &   zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited,
     &   ERROR,*9999)
            Q_sheet(gen)=(Pin_sheet-Pout_sheet)/SHEET_RES
            Q_gen=Q_sheet(gen)*2**gen
            R_sheet(gen)=SHEET_RES
C...    Update soln matrix (delta p=QR across a cap)
           SparseVal(NonZeros-num_symm_gen-1-3*(num_symm_gen-gen))=
     &       -SHEET_RES
 
C...   SECOND HALF OF ARTERIOLE
C...    Update radius of arteriole based on inlet pressure NB: May need to
C....   update this later to give an average of inlet and outlet radii
      IF (Pressure(4*gen-2)-P_exta.LE.Pub_a_v)THEN
c            radupdate=rad_a(gen)*(1+(Pressure(4*gen-2)-P_exta)*alpha_a)
            radupdate=rad_a(gen)+alpha_a*(Pressure(4*gen-2)-P_exta)*
     &       (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*
     &       (Pressure(4*gen-2)-P_exta)/(2.d0-2.d0*num_symm_gen)
         ELSE
c            radupdate=rad_a(gen)*(1+Pub_a_v*alpha_a)
            radupdate=rad_a(gen)+alpha_a*Pub_a_v*(gen-num_symm_gen)
     &       /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v
     &       /(2.d0-2.d0*num_symm_gen)
         ENDIF
         IF(nj_hypoxia.NE.0) THEN
           IF(rad_a(gen).LE.0.25d0)THEN
              radupdate=radupdate*k_factor   
           ENDIF
         ENDIF
         volume_vessels=volume_vessels+
     &     (2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_a(gen)/2.d0)*1000.d0**3.d0)

C...  Calculate Poiseuille resistance in second half of arteriole - (only
C...   half total generation length)
      R_art2=(8*mu_app(gen)*L_a(gen)/2.d0)/(pi*radupdate**4.d0);

C###   SECOND HALF OF VENULE
C...    Update radius - linear with pressure or constant at high pressure
      IF (Pressure(4*gen)-P_extv.LE.Pub_a_v)THEN
c           radupdate=rad_v(gen)*(1+(Pressure(4*gen)-P_extv)*alpha_v)
            radupdate=rad_v(gen)+alpha_v*(Pressure(4*gen)-P_extv)*
     &       (gen-num_symm_gen)/(1.d0-num_symm_gen)+alpha_c*(1-gen)*
     &       (Pressure(4*gen)-P_extv)/(2.d0-2.d0*num_symm_gen)
         ELSE
c            rad_update=rad_v(gen)*(1+Pub_a_v*alpha_v)
            radupdate=rad_v(gen)+alpha_v*Pub_a_v*(gen-num_symm_gen)
     &       /(1.d0-num_symm_gen)+alpha_c*(1-gen)*Pub_a_v
     &       /(2.d0-2.d0*num_symm_gen)
         ENDIF
         volume_vessels=volume_vessels+
     &     (2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0) !mm3
       test=(2.d0**gen*(pi*radupdate**2.d0*L_v(gen)/2.d0)*1000.d0**3.d0)
C...  Poiseuille resistance in second half of venule
      R_ven2=(8.d0*mu_app(gen)*L_v(gen)/2.d0)/(pi*radupdate**4.d0)

C...  First generation
       IF (gen.EQ.1) THEN
          SparseVal(2)=-R_art1/2.d0
          SparseVal(5)=R_ven1
          SparseVal(6)=-R_ven1/2.d0
          SparseVal(8)=-R_art2/2.d0
          SparseVal(11)=R_ven2
          SparseVal(12)=-R_ven2/2.d0
          count=12
       ELSE
        DO j=2,gen
          SparseVal(count+3+(j-2))=R_art1/(2.d0**(gen+1-j))
        ENDDO
      count=count+gen+2
      SparseVal(count)=-R_art1/2.d0**gen
        DO j=2,gen+1
          SparseVal(count+3+(j-2))=R_ven1/(2.d0**(gen+1-j))
        ENDDO
      count=count+gen+3
      SparseVal(count)=-R_ven1/2.d0**gen
         DO j=2,gen
          SparseVal(count+3+(j-2))=R_art2/(2.d0**(gen+1-j))
        ENDDO      
      count=count+gen+2
      SparseVal(count)=-R_art2/2.d0**gen
        DO j=2,gen+1
          SparseVal(count+3+(j-2))=R_ven2/(2.d0**(gen+1-j))
        ENDDO
      count=count+gen+3
      SparseVal(count)=-R_ven2/2.d0**gen
      ENDIF

      ENDDO
C...  ------------FINAL GENERATION----------------------------
C... --------The capillaries covering the alveolar sacs---------
C...These are just capillary beds without an associated arteriole/venule - 
        Pin_sheet=Pressure(4*num_symm_gen-6); !%pressure into final capillary sheets
        Pout_sheet=Pressure(4*num_symm_gen-4); !%pressure out of final capillary sheets
        CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, 
     &   zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited,
     &   ERROR,*9999)
        SparseVal(NonZeros-num_symm_gen-1)=-SHEET_RES
            Q_sheet(num_symm_gen)=(Pin_sheet-Pout_sheet)/SHEET_RES
            Q_gen=Q_sheet(num_symm_gen)*2**gen
            R_sheet(num_symm_gen)=SHEET_RES


C###  EXIT SUBROUTINE 'POPULATE_MATRIX_LADDER'
      CALL EXITS('POPULATE_MATRIX_LADDER')
      RETURN
 9999 CALL ERRORS('POPULATE_MATRIX_LADDER',ERROR)
      CALL EXITS('POPULATE_MATRIX_LADDER')
      RETURN 1
      END

