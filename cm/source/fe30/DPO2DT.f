      SUBROUTINE DPO2DT(n,t,Y,F,RPAR,IPAR)
   
C#### Subroutine: BENTAL_2006
C###  Description:
C###    ODE for rate of change of blood PO2, based on the Ben-Tal 2006 
C###    gas exchange model. This is called from RADAU5.f (ode solver). 
C###    Y = po (PO2 in capillary blood)
C**** Created by AJS, June 2009

      IMPLICIT NONE    
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'odes00.cmn'  
        
!     Parameter List
      INTEGER n
      REAL*8 F(n),IPAR(0:MAX_PAR),RPAR(0:MAX_PAR),t,Y(n)
!     Local Variables
      REAL*8 Doxy,dfdp,fO2,PA,Th,Vc
        
      CALL ENTERS('DPO2DT',*9999) 

C     Get parameters
      PA=RPAR(1)
      Vc=RPAR(2) 
      Doxy=RPAR(3)
      fO2=RPAR(4)
      Th=RPAR(5)

      IF(Vc.EQ.0.d0.OR.Doxy.EQ.0.d0)THEN
        F(1)=0.d0

      ELSE
C       Derivative function
        dfdp=((L*(1+kt*sigma_o*Y(1))**4+(1+kr*sigma_o*Y(1))**4)*
     &    (3*L*kt**2*sigma_o**2*Y(1)*(1+kt*sigma_o*Y(1))**2+ 
     &    L*kt*sigma_o*(1+kt*sigma_o*Y(1))**3+3*kr**2*sigma_o**2*Y(1)
     &    *(1+kr*sigma_o*Y(1))**2+kr*sigma_o*(1+kr*sigma_o*Y(1))**3) 
     &    -(L*kt*sigma_o*Y(1)*(1+kt*sigma_o*Y(1))**3+kr*sigma_o*Y(1)
     &    *(1+kr*sigma_o*Y(1))**3)*(4*L*kt*sigma_o*(1+kt*sigma_o
     &    *Y(1))**3+4*kr*sigma_o*(1+kr*sigma_o*Y(1))**3))/ 
     &    ((L*(1+kt*sigma_o*Y(1))**4+(1+kr*sigma_o*Y(1))**4)**2)

        F(1)=Doxy/standard_molar_vol/(sigma_o*Vc)*(1/(1+(4*Th/sigma_o)*
     &    dfdp))*(fO2*(PA-p_water)-Y(1))

      ENDIF ! Vc or Doxy or Dc = 0
!       write(*,*) 'Doxy=',Doxy,' standard_molar_vol=',standard_molar_vol
!       write(*,*) 'sigma_o=',sigma_o,' Vc=',Vc,' Th=',Th,' dfdp=',dfdp
!       write(*,*) 'fO2=',fO2,' PA=',PA,' pw=',p_water
!       write(*,*) 'PbO2(initial)=',Y(1),' F=',F(1)

 9999 CALL EXITS('DPO2DT')
      RETURN
      END 
