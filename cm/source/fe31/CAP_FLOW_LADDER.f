      SUBROUTINE CAP_FLOW_LADDER(ne,HEIGHT,LPM_FUNC,Pin,Pout,Ppl,
     &  Q01,R_in,R_out,x,y,z,k_factor,OUTPUT_PERFUSION,ERROR,*)  

C#### Subroutine CAP_FLOW_LADDER
C###  CREATED: JUNE 2009 by ARC

C###  Description:
C###  This function is uses the ladder model for perfusion in the acinus.
C###  It uses a symmetric, bifuracting arteriole and venule tree with N
C###  generations and calls a function to calculate resistance across
C###  a capillary sheet at each generation (CAP_FLOW_SHEET.f).  


C###  -----------------------INPUT-------------------------
C###  The input to this subroutine is: 
C###  ne=element number
C###  Pin= Pressure into the acinus
C###  Pout=Pressure out of the acinus
C###  Ppl= Pleural pressure

C###  ----------------------OUTPUT----------------------------
C###  Important output to large vessel models 
C###  LPM_FUNC= Resistance across the acinus:

C###  This needs to be fed back into the large vessel model. Then 
C###  Pin-Pout=LPM_FUNC*Q is solved for Pin, Pout and Q as part 
C###  that system. Pin and Pout can then be fed back iteratively 
C###  into this subroutine...  
C###  
C###  In addition this subroutine outputs:
C###  Pressure at each each vessel interection and flow, resistance
C###  and RBC transit times through each capillary element. 


C###  UNITS. The units that are used here are m, 

C!!  NB/ Palv is stored in common block lung00.cmn


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi10.cmn'

     
      !INPUT AND OUTPUT PARAMETER LIST
      INTEGER ne
      REAL*8 HEIGHT(3),LPM_FUNC,Pin,Pout,Ppl,Q01,R_in,R_out,
     &  x,y,z,k_factor
      LOGICAL OUTPUT_PERFUSION
      CHARACTER ERROR*(*)

      !Local variables
      INTEGER MatrixSize,NonZeros,submatrixsize
      INTEGER i,N
      REAL*8  area,
     &  Pressure(4*num_symm_gen-4),
     &  Q01_mthrees, 
     &  sheet_number

C###  ENTER SUBROUTINE 'CAP_FLOW_LADDER'
      CALL ENTERS('CAP_FLOW_LADDER',*9999)
      CALL ASSERT(num_symm_gen.LT.13,'>>max symmetrics [13] exceeded!',
     &   ERROR,*9999)

C###  Number of non-zero entries in solution matrix. 
      NonZeros=3
      DO i=2,num_symm_gen
         NonZeros=NonZeros+4*i+10
      ENDDO
C###  The size of the solution matrix (number of unknown pressures and flows)
      MatrixSize=5*num_symm_gen-3
C###  The number of unknown pressures
      submatrixsize=4*num_symm_gen-4
  
C...  ---INITIALISATION
C...  The input Q01 gives us an estimate for flow into the acinus from the large
C...  vessel model.
C...  This is in mm^3/s and needs to be converted to m^3/s to use in calculating
C...  arteriole and venule resistance
      Q01_mthrees=Q01/1.d9 !mm3/s->m3/s 
C###  Sheet area (unscaled): 
C...  We define a sheet area for input into the capillary model. 
C...  This area is at full inflation and will be scaled within CAP_FLOW_SHEET  
C...  Area of an individual sheet
      sheet_number=0
      DO i=1,num_symm_gen
         sheet_number=sheet_number+2.d0**i
      ENDDO
      area=total_cap_area/sheet_number!m^2
C...  Initial guess for pressure distribution lets say all arterial pressures are the same 
C...  and all the venous pressures are the same solution appears independent of this.
      DO i=1,num_symm_gen-1
          Pressure(4*i-3)=1000 ! Pa
          Pressure(4*i-2)=1000
          Pressure(4*i-1)=100
          Pressure(4*i)=100
      ENDDO
     
C###  ---CALL THE FUNCTIONS THAT CALCULATE THE FLOW ACROSS THE LADDER FOR A GIVEN PRESSURE DROP--
      CALL EVALUATE_LADDER(ne,NonZeros,MatrixSize,submatrixsize,
     & area,HEIGHT,Pin,Pout,Ppl,Pressure,Q01_mthrees,x,y,z,k_factor,
     & OUTPUT_PERFUSION,ERROR,*9999)

C###  ---FINAL FUNCTION OUTPUT (Resistance across ladder)---
C...  This takes difference between the inlet and outlet pressures 
C...  (Pin and Pout) and divides by an updated flow (Q01_mthrees)
C...  to give updated resistance across the ladder. This feeds back
C...  to the large vessel model. 
      LPM_FUNC=(Pin-Pout)/(Q01_mthrees*1000.d0**3) !Pa.s/m^3->pa.s/mm^3

C###  EXIT SUBROUTINE 'CAP_FLOW_LADDER'
      CALL EXITS('CAP_FLOW_LADDER')
      RETURN
 9999 CALL ERRORS('CAP_FLOW_LADDER',ERROR)
      CALL EXITS('CAP_FLOW_LADDER')
      RETURN 1
      END


