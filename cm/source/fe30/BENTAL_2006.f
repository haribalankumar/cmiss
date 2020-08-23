      SUBROUTINE BENTAL_2006(n,t,Y,F,RPAR,IPAR)
   
C#### Subroutine: BENTAL_2006
C###  Description:
C###    System of equations for the Ben-Tal 2006 gas exchange model. 
C###    This is called from RADAU5.f (ode solver). The code assumes
C###    a breath-hold (i.e. q=0) over the timestep.
C**** Created by AJS, May 2006

      IMPLICIT NONE    
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'odes00.cmn'  
        
!     Parameter List
      INTEGER n      
      REAL*8 F(n),IPAR(0:MAX_PAR),RPAR(0:MAX_PAR),t,Y(n)
!     Local Variables
      INTEGER i
      REAL*8 Dc,Doxy,dfdp,E,fom,fcm,Pm,Th,Va,Vc,q,Qa 
        
      CALL ENTERS('BENTAL_2006',*9999) 

C     Get parameters
      Pm=GASPAR(1)
      fom=GASPAR(2)
      fcm= GASPAR(3)
      E=GASPAR(4)
      Va=GASPAR(5)
      Vc=GASPAR(6) 
      Th=GASPAR(7)
      Doxy=GASPAR(8)
      Dc=GASPAR(9)

C     Calculate parameters
      q=0.0d0 !flow into unit
      foi=Y(2)
      fci=Y(3)
      ycd=Y(3)
      yod=Y(2)
      Qa=q + Dc*(Y(5)-(Y(3)*(Y(1)-p_water)))+Doxy*(Y(4)-
     &  (Y(2)*(Y(1)-p_water)))

C     Initialise
      DO i=1,n
        F(i)=0.0d0
      ENDDO

C     System of equations for derivative functions  Y = [Pa, fo, fc, po, pc, z]  
      F(1)=Pm*E*Qa/Y(1)

      IF(CONST_ALV_PRESSURES)THEN!hold fo, fc constant
        F(2)=0.0d0
        F(3)=0.0d0
      ELSE 
        F(2)=1/Va*(Doxy*(Y(4)-Y(2)*(Y(1)-p_water))+ 
     &    (foi - Y(2))*q-Y(2)*(Dc*(Y(5)-Y(3)*(Y(1)-
     &    p_water))+Doxy*(Y(4)-Y(2)*(Y(1)-p_water))))
        F(3)=1/Va*(Dc*(Y(5)-Y(3)*(Y(1)-p_water))+ 
     &    (fci-Y(3))*q-Y(3)*(Doxy*(Y(4)-Y(2)*(Y(1)-
     &    p_water))+Dc*(Y(5)-Y(3)*(Y(1)-p_water))))
      ENDIF

      dfdp=((L*(1+kt*sigma_o*Y(4))**4+(1+kr*sigma_o*Y(4))**4)*
     &  (3*L*kt**2*sigma_o**2*Y(4)*(1+kt*sigma_o*Y(4))**2+ 
     &  L*kt*sigma_o*(1+kt*sigma_o*Y(4))**3+3*kr**2*sigma_o**2*Y(4)
     &  *(1+kr*sigma_o*Y(4))**2+kr*sigma_o*(1+kr*sigma_o*Y(4))**3) 
     &  -(L*kt*sigma_o*Y(4)*(1+kt*sigma_o*Y(4))**3+kr*sigma_o*Y(4)
     &  *(1+kr*sigma_o*Y(4))**3)*(4*L*kt*sigma_o*(1+kt*sigma_o
     &  *Y(4))**3+4*kr*sigma_o*(1+kr*sigma_o*Y(4))**3))/ 
     &  ((L*(1+kt*sigma_o*Y(4))**4+(1+kr*sigma_o*Y(4))**4)**2)

      F(4)=Doxy/standard_molar_vol/(sigma_o*Vc)*(1/(1+(4*Th/sigma_o)*
     &  dfdp))*(Y(2)*(Y(1)-p_water)-Y(4))
    
      F(5)=Dc/standard_molar_vol/(sigma_c*Vc)*(Y(3)*(Y(1)-
     &   p_water)-Y(5))+delta*ltwo/sigma_c*hconc*Y(6)-
     &   delta*rtwo*Y(5)
            
      
      F(6)=delta*rtwo*sigma_c*Y(5)-delta*ltwo*hconc*Y(6)
      
 9999 CALL EXITS('BENTAL_2006')
      RETURN
      END 
