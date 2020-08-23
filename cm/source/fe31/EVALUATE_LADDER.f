      SUBROUTINE EVALUATE_LADDER(ne,NonZeros,MatrixSize,submatrixsize,
     & area,HEIGHT,Pin,Pout,Ppl,Pressure,Q01_mthrees,x,y,z,k_factor,
     & OUTPUT_PERFUSION,ERROR,*)
C#### Subroutine EVALUATE_LADDER
C###  CREATED: JAN 2010 by ARC

C###  Description:
C###  Sets up and solves matrix equations for ladder model

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'
      !INPUT AND OUTPUT PARAMETER LIST
      INTEGER ne,NonZeros,MatrixSize,submatrixsize
      REAL*8 area,HEIGHT(3),Pin,Pout,Ppl,Pressure(submatrixsize),
     & Q01_mthrees,x,y,z,k_factor
      LOGICAL OUTPUT_PERFUSION
      CHARACTER ERROR*(*)
      INTEGER i,iter,j,gen,SparseCol(NonZeros),
     & SparseRow(MatrixSize+1),zone,nj,num_sheet
      REAL*8 area_new,ErrorEstimate,Hart,Hven,Pin_SHEET,Pout_SHEET,
     &      Q_c,Qtot,Qgen,Q_sheet(num_symm_gen),RBC_TT,
     &      RHS(MatrixSize),Rtot,SHEET_RES,Solution(MatrixSize),
     &      SolutionLast(MatrixSize),SparseVal(Nonzeros),TOTAL_CAP_VOL,
     &      TOTAL_SHEET_H,TOTAL_SHEET_SA,P_inlet,P_outlet,
     &      R_upstream,R_downstream,GRAVITY,
     &      recruited,TT_TOTAL


      CALL ENTERS('EVALUATE_LADDER',*9999)
C###  INITIAL SOLUTION GUESS 
      DO i=1,submatrixsize
      solution(i)=Pressure(i)
      ENDDO
      DO i=submatrixsize+1,matrixSize-1
       solution(i)=Q01_mthrees/2**num_symm_gen
      ENDDO
      Solution(Matrixsize)=Q01_mthrees
C###  INITIALISE SOLUTIONLAST
      DO j=1,MatrixSize
        SolutionLast(j)=Solution(j)
      ENDDO
     
C### INPUT TO THE LADDER MODEL THAT IS INDEPENDENT OF ITERATION
      CALL LADDERSOL_MATRIX(NonZeros,MatrixSize,submatrixsize,
     &    SparseCol,SparseRow,SparseVal,RHS,Pin,Pout,ERROR,*9999)

C### ITERATIVE LOOP
      iter=0 
      ErrorEstimate=1.d10
      DO WHILE(ErrorEstimate.GT.1.0d-9.AND.iter.LT.100)
        iter=iter+1
C...  CALCULATE RESISTANCE GIVEN CURRENT PRESSURE AND FLOW - THEN UPDATE 
C...  SparseVal- THese are the only elements of the solution matrix that need
C.... iteratively updating
        CALL POPULATE_MATRIX_LADDER(ne,NonZeros,MatrixSize,
     &    submatrixsize,SparseCol,SparseRow,area,Pin,Pout,Ppl,Pressure,
     &    Q01_mthrees,Q_sheet,RHS,SparseVal,k_factor,ERROR,*9999)
 
        CALL MINI_LINEAR_SYSTEM(MatrixSize,NonZeros,SparseCol,SparseRow,
     &   Solution,SparseVal,RHS,ERROR,*9999)
         DO j=1,submatrixsize
            Pressure(j)=Solution(j)
         ENDDO
         Q01_mthrees=Solution(MatrixSize)   
C     Estimating Error in solution
      ErrorEstimate=0.d0
        DO i=1,MatrixSize
          ErrorEstimate=ErrorEstimate+
     &     DABS((Solution(i)-SolutionLast(i))**2.d0
     &     /Solution(i)**2.d0)
          SolutionLast(i)=Solution(i)
        ENDDO
      ErrorEstimate=ErrorEstimate/MatrixSize
      ENDDO

      Qtot=0
      Do i=1,num_symm_gen
        Qtot=Qtot+Q_sheet(i)*2.d0**i
      ENDDO
      Rtot=(Pin-Pout)/Q01_mthrees
       
       IF(OUTPUT_PERFUSION)THEN
C###  GET SOLUTIONS TO WRITE TO FILE 
        TOTAL_CAP_VOL=0.d0
        TOTAL_SHEET_SA=0.d0
        TT_TOTAL=0.d0
        TOTAL_SHEET_H=0.d0
        num_sheet=0
C... Ladder output
        DO i=1,num_symm_gen-1
           gen=i
           Pin_sheet=Pressure(4*i-3)
           Pout_sheet=Pressure(4*i-1)
C...     calulate resistance and transit time through a single capillary
           CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, 
     &      zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited,
     &      ERROR,*9999)
           num_sheet=num_sheet+2**gen
           Qgen=Q_c*2.d0**i
         WRITE(IOFILE3,
     &   '(I6,X,3(F9.2,X),I6,X,4(F8.2,X),4(F8.5,X),
     &    2(F10.2,X),3(F8.4,X),I6,X,2(F10.5,X),2(F8.4,X),
     &    (F10.2,X))')
     &    ne,x,y,z,gen,Pin,Pin_sheet,Pout_sheet,Pout,Qtot*1.d9,
     &    Qgen*1.d9,Q_c*1.d9,(Hart-Hven)*1.d6,SHEET_RES/1000.d0**3.d0,
     &    Rtot/1000.d0**3.d0,RBC_tt,Hart*1.d6,Hven*1.d6,zone,
     &    (Hart/2.d0+Hven/2.d0)*area_new*1.d9*(2.d0**gen),
     &    area_new*1.d6*(2.d0**gen),recruited


           TOTAL_CAP_VOL=TOTAL_CAP_VOL
     &       +(Hart/2.d0+Hven/2.d0)*area_new*1.d9*(2.d0**gen)
           TOTAL_SHEET_H=TOTAL_SHEET_H
     &       +(Hart/2.d0+Hven/2.d0)*(2.d0**gen)*1.d6
          TOTAL_SHEET_SA=TOTAL_SHEET_SA
     &       +area_new*1.d6*(2.d0**gen)
          TT_TOTAL=TT_TOTAL+RBC_tt*(2.d0**gen)
        ENDDO

        gen=num_symm_gen 
        Pin_sheet=Pressure(4*num_symm_gen-6); !%pressure into final capillary sheets
        Pout_sheet=Pressure(4*num_symm_gen-4); !%pressure out of final capillary sheets
         CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT, 
     &   zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited,
     &   ERROR,*9999)
         num_sheet=num_sheet+2**gen
           Qgen=Q_C*2.d0**num_symm_gen
         WRITE(IOFILE3,
     &   '(I6,X,3(F9.2,X),I6,X,4(F8.2,X),4(F8.5,X),
     &   2(F10.2,X),3(F8.4,X),I6,X,2(F10.5,X),3(F8.4,X))')
     &    ne,x,y,z,gen,Pin,Pin_sheet,Pout_sheet,Pout,Qtot*1.d9,
     &    Qgen*1.d9,
     &    Q_c*1.d9,(Hart-Hven)*1.d6,SHEET_RES/1000.d0**3.d0,
     &    Rtot/1000.d0**3.d0,RBC_tt,Hart*1.d6,Hven*1.d6,zone,
     &    R_upstream,R_downstream,
     &    (Hart/2.d0+Hven/2.d0)*area_new*1.d9*(2.d0**gen),
     &    area_new*1.d6*(2.d0**gen),recruited

          TOTAL_CAP_VOL=TOTAL_CAP_VOL
     &       +(Hart/2.d0+Hven/2.d0)*area_new*1.d9*(2.d0**gen)
           TOTAL_SHEET_H=TOTAL_SHEET_H
     &       +(Hart/2.d0+Hven/2.d0)*(2.d0**gen)*1.d6
          TOTAL_SHEET_H=TOTAL_SHEET_H/num_sheet
          TOTAL_SHEET_SA=TOTAL_SHEET_SA
     &       +area_new*1.d6*(2.d0**gen)
            TT_TOTAL=TT_TOTAL+RBC_tt*(2.d0**gen)
            TT_TOTAL=TT_TOTAL/num_sheet

C... General output
! ne=1  |  x=2  |  y=3  |  z=4  | Pin=5 Pa |Pout=6 Pa | Qtot=7 mm^3/s |sum Qsheet=8 mm^3 /s|
! Rtot=9 Pa/mm^3 | Blood_vol=10 mm^3| sheet_area= 11 mm^2 | ave_TT=12 s |ave_H=13 um |Ppl=14 Pa
          WRITE(IOFILE2,
     &   '(I6,X,5(F9.2,X),2(F8.5,X),F10.2,X,F8.4,X,F10.4,X,F10.3,X,F8.4,
     &    X,F9.4,X)') 
     &    ne,x,y,z,Pin,Pout,
     &    Q01_mthrees*1.d9,Qtot*1.d9,Rtot/1000.d0**3.d0,TOTAL_CAP_VOL,
     &    TOTAL_SHEET_SA,TT_TOTAL,TOTAL_SHEET_H,Ppl
        ENDIF
C###  EXIT SUBROUTINE 'EVALUATE_LADDER'
      CALL EXITS('EVALUATE_LADDER')
      RETURN
 9999 CALL ERRORS('EVALUATE_LADDER',ERROR)
      CALL EXITS('EVALUATE_LADDER')
      RETURN 1
      END

