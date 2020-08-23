      SUBROUTINE CAP_FLOW_MBA(ne,NORD,LPM_FUNC,Pin,Pout,Ppl,Q01,
     & x,y,z,OUTPUT_PERFUSION,ERROR,*)

C#### Subroutine: CAP_FLOW_SIMPLE
C###  Description:
C###    This subroutine uses the simplified sheet flow model developed 
C###    by Fung & Sobin (1969) to describe blood flow through the 
c###    densely packed pulmonary capillary blood vessels to develop a 
C###    lumped parameter model to coupled arterial-capillary-venous 
C###    flow in the lung.
C###
C###    Pressures are in Pa. Pin/Pout=blood pressues in and out, 
C###    respectively, Ppl=pleural pressure, PA=alveolar pressure.

C###   This model also includes x # symmetric generations of arteriole 
C###   and venules vessels assuming Poiseuille resistance within them.

C**** Created by AJS & KSB, April 2008
C!!  NB/ Palv is stored in common block lung00.cmn
C**** Edited by ARC, July 2009: The sheet flow element of small vessel flow
C**** has been saved in a separate subroutine (CAP_FLOW_SHEET.f) This means 
C***  that the same code can be used in both the symmetric branching acinus
C***  the ladder acinus.   

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'

      !Parameter list
      INTEGER ne,NORD(5,ne)
      REAL*8 LPM_FUNC,Pin,Pout,Ppl,Q01,
     & x,z,y
      LOGICAL OUTPUT_PERFUSION
      CHARACTER ERROR*(*)

      !Local variables
      INTEGER gener,order,zone
      REAL*8 area,area_new,Hart,Hven,Q01_mthrees,
     &   Q_c,RBC_tt,recruited,SHEET_RES
	
      CALL ENTERS('CAP_FLOW_MBA',*9999)
C...  The input Q01 gives us an estimate for flow into the acinus from the large vessel model.
C...  This is in mm^3/s and needs to be converted to m^3/s to use in calculating arteriole and 
C...  venule resistance
      Q01_mthrees=Q01/1.d9 !mm3/s->m3/s 
C###  Sheet area (unscaled): 
C...  We need to define a sheet area for input into the capillary model. T
C...  This area should be at full inflation and will be scaled within CAP_FLOW_SHEET  
C...  Area of an individual sheet
      gener=NORD(1,ne)
      order=NORD(3,ne)
      IF(gener.LE.1)THEN!1st generation is transitional bronchiole
        area=1.d-20
      ELSE !alveolarised
        area=total_cap_area
      ENDIF
C...  CALL SHEET FLOW MODEL
      CALL CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_tt,zone,
     &   Pin,Pout,Ppl,area,area_new,recruited,*9999)
C...  Total resistance across arterioles, capillaries and venules:
      LPM_FUNC=SHEET_RES/1000.d0**3.d0  !Change units to Pa.s/mm^3      
      IF(OUTPUT_PERFUSION)THEN
        IF(gener.GT.1)THEN
!        TOTAL_CAP_VOL=(Hart/2.d0+Hven/2.d0)*area_new*1.d9
!        volume_sheet=volume_sheet+TOTAL_CAP_VOL
          WRITE(IOFILE2,
     &             '(I6,X,I6,X,I6,X,3(F8.2,X),F8.5,X,F10.1,X,4(F8.5,X),
     &              I6,X,F12.2,X,3(F8.2,X),F8.2,X,2(F8.4,X))')
     &            ne,gener,order,Pin,Pout,Pin-Pout,
     &            Q_c*1.d12,LPM_FUNC,(Hart/2.d0+Hven/2.d0)*area_new*1.d9 
     &            ,area_new*1.d6,Hart*1.d6,
     &            Hven*1.d6,zone,RBC_tt,(Hart-Hven)*1.d6,x,y,z,
     &            recruited,Ppl
        ENDIF
      ENDIF

      CALL EXITS('CAP_FLOW_MBA')
      RETURN
 9999 CALL ERRORS('CAP_FLOW_MBA',ERROR)
      CALL EXITS('CAP_FLOW_MBA')
      RETURN 1
      END
