      SUBROUTINE SUM_RESIS_COMPL(BBM,CE,COMP_SUM,CW,
     & NEELEM,NELIST,NXI,NORD,NO_SUBTENDED,
     & RESIS_SUM,VOL_SUM,undef,ERROR,*)
    
      INCLUDE 'b00.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'
      
!     Parameter List
      INTEGER NEELEM(0:NE_R_M),NELIST(0:NEM),NORD(5,NE_R_M),
     & NO_SUBTENDED(NE_R_M),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),COMP_SUM(NE_R_M),CW,
     &  RESIS_SUM(NE_R_M),VOL_SUM(NE_R_M),RESIS_AV(NE_R_M),
     &  COMP_AV(NE_R_M),undef
      CHARACTER ERROR*(*)     
!     Local variables
      INTEGER noelem,ne,noelem2,ne2,term,start,subtended,
     & term_ne(NE_R_M),
     & asym_match,test_above,test1,test2,start_gen,kount,
     & end_gen,sym_resisG,sym_resisT,asym_resis,
     & SYM_SUBTENDED(NE_R_M)  
      REAL*8 sum_vol,sum_comp,sum_resis,inv_sum_resis,check,check2,
     & sum_av,ref_ratio,vol_equiv,ratio,lambda,lambda15,
     & exp_term,comp,cc,min_E,a,b,c,RFROMC   
     
      CALL ENTERS('SUM_RESIS_COMPL',*9999)
!    This subroutine calcualtes the summed resistance and compliance for 
!    each tree. Note on entry BBM(1,ne) and CW(nm_C,ne) are >0.0d0 only for t
!    terminal

C     indices for CE arrays
      nm_R=1
      nm_C=2
      nm_Ppl=3
      nm_Rt=4
      nm_lambda=5
      nm_dPl=6
      nm_volumes=7
      nm_vinit=8
      nm_vol_bel=9
      nm_R_trunc=10
      nm_C_trunc=11
      nm_ST_trunc=12
      nm_C_av=13
      nm_R_av=14

!    Have to intialise these variables, they seem to get messy 
   	  Do noelem=1,NEELEM(0)
      	  ne=NEELEM(noelem)
      	  VOL_SUM(ne)=0.0d0
      	  RESIS_SUM(ne)=0.0d0
      	  COMP_SUM(ne)=0.0d0
      	  NO_SUBTENDED(ne)=0
      	  SYM_SUBTENDED(ne)=0
      	  CE(nm_C_av,ne)=0.0d0
      	  CE(nm_R_av,ne)=0.0d0
      	  term_ne(ne)=0
    	 ENDDO   

!intialise the volumes/compliances at the terminals
       term=0
       check=0
       DO noelem=1,NEELEM(0)
	     	ne=NEELEM(noelem)
		IF(NXI(1,0,ne).EQ.0)THEN
			VOL_SUM(ne)=BBM(1,ne)
			check=check+VOL_SUM(ne)
			COMP_SUM(ne)=CE(nm_C,ne)
			term=term+1
			term_ne(term)=ne
			NO_SUBTENDED(ne)=1
		ENDIF	
       ENDDO
!    Volume and Compliance,no_subtended in asymmetric model   
 	       DO noelem=NEELEM(0),1,-1
     		ne=NEELEM(noelem)
    	 	sum_vol=0.0d0
    	 	subtended=0
     		  !note terminal elements will have no children
     		  IF(NXI(1,0,ne).GT.0)THEN
     		  	DO noelem2=1,NXI(1,0,ne) !for each child
             			ne2=NXI(1,noelem2,ne) !get the child element number
             			sum_vol=sum_vol+VOL_SUM(ne2) !add in parallel
             			subtended=subtended+NO_SUBTENDED(ne2)
          	 	ENDDO
          	 ENDIF
         	!add the volumes and compliances to the one above 
         	VOL_SUM(ne)=VOL_SUM(ne)+sum_vol
                CE(nm_vinit,ne)=CE(nm_vinit,ne)+smum_vol
         	NO_SUBTENDED(ne)=NO_SUBTENDED(ne)+subtended
         	!only add if sum_comp non zero and dont invert possible zeros
c         	IF(sum_comp.gt.0.0d0)THEN
c                        IF(COMP_SUM(ne).lt.1.0d-6)THEN
c         		COMP_SUM(ne)=COMP_SUM(ne)+sum_comp
c         		ELSE
c         		COMP_SUM(ne)=1.0d0/COMP_SUM(ne)+1.0d0/sum_comp
c         		COMP_SUM(ne)=1.0d0/COMP_SUM(ne)
c         		ENDIF
c         	ENDIF
     		ENDDO
                
       
!find the equivalent number of terminals in a symmetric model       
        Do noelem=1,neelem(0)
       		ne=NEELEM(noelem)
			IF(NELIST(ne).eq.1)THEN
			  ne0=NXI(1,1,ne)
			    IF(NELIST(ne0).eq.0)THEN	
				!find the best match for the number of terminals
				test=0
				asym_match=0
				Do while(test<=NO_SUBTENDED(ne))
					asym_match=asym_match+1
					test=2**asym_match
					test_above=2**(asym_match+1)
				ENDDO
					test1=abs(test-NO_SUBTENDED(ne))
					test2=abs(test2-NO_SUBTENDED(ne))
					!test which estimate of no generations is closer to nos of acinii and save
					IF(test1>test2)THEN
						asym_match=asym_match+1
					ENDIF 	
						SYM_SUBTENDED(ne)=2**asym_match
			ENDIF !terminal NELIST(ne0).eq.0
		 ENDIF !NELIST(ne).eq.1 (in tree)	
      	ENDDO !each element
      	
!recalculate the compliance for each unit
 !next step is for each unit I will need to scale the reference volume so that the compliance scales correctly
	   !store in XAB(5,ne)=reference volume for each terminal unit
	  		  a=SEDF_COEFFS(1)  !dimensionless
      			  b=SEDF_COEFFS(2) !dimensionless
      			  c=SEDF_COEFFS(3) !Pa
      			  cc=c
	   Do noelem=1,neelem(0)
       		ne=NEELEM(noelem)
			IF(NELIST(ne).eq.1)THEN
			  ne0=NXI(1,1,ne)
			  IF(NELIST(ne0).eq.0)THEN
			  check=(NO_SUBTENDED(ne))
			  check2=(SYM_SUBTENDED(ne))	
			  ref_ratio=check/check2 !ratio of number of units in asymmetric tree
			  vol_equiv=VOL_SUM(ne)/SYM_SUBTENDED(ne)
			  !CALCULATE the LUMPED TISSUE COMPLIANCE
          		  ratio=vol_equiv/(undef*ref_ratio)
         		  lambda = ratio**(1.d0/3.d0) !uniform extension ratio 
                 exp_term=DEXP(0.75d0*(3.d0*a+b)*(lambda**2-1.d0)**2)
		          IF(ratio.GE.1.5d0)THEN
		           comp=cc*exp_term/6.d0*(3.d0*(3.d0*a+b)**2
     &        *(lambda**2-1.d0)**2/lambda**2+(3.d0*a+b)
     &        *(lambda**2+1.d0)/lambda**4)
         		 ELSE
            		    lambda15=1.5d0**(1.d0/3.d0)
     		           min_E=cc*exp_term/6.d0*(3.d0*(3.d0*a+b)**2
     &        *(lambda15**2-1.d0)**2/lambda15**2+(3.d0*a+b)
     &        *(lambda15**2+1.d0)/lambda15**4)
            		   comp=cc*0.17d0+2.d0*(ratio-1.d0)*
     &        (min_E-cc*0.17d0)
          		ENDIF
	         comp=vol_equiv/comp !in units of volume/pressure
	         IF(CW.GT.0.0d0)THEN
	        	 comp=1.0d0/(1.0d0/comp+1.d0/CW)
	         ENDIF
	         COMP_SUM(ne)=comp*SYM_SUBTENDED(ne) !because all paths are equal
	         ENDIF
	        ENDIF 
	      ENDDO   
	      
	      
	       DO noelem=NEELEM(0),1,-1
     		ne=NEELEM(noelem)
     		IF(NELIST(ne).eq.1)THEN
                  ne0=NXI(1,1,ne)
     		  IF(NELIST(ne0).eq.0)THEN
     		  	DO noelem2=1,NXI(1,0,ne) !for each child
             			ne2=NXI(1,noelem2,ne) !get the child element number
             			sum_comp=sum_comp+COMP_SUM(ne2) !add in parallel
          	 	ENDDO
          	 ENDIF
         	!only add if sum_comp non zero and dont invert possible zeros
         	IF(sum_comp.gt.0.0d0)THEN
                        IF(COMP_SUM(ne).lt.1.0d-6)THEN
         		COMP_SUM(ne)=COMP_SUM(ne)+sum_comp
         		ELSE
         		COMP_SUM(ne)=1.0d0/COMP_SUM(ne)+1.0d0/sum_comp
         		COMP_SUM(ne)=1.0d0/COMP_SUM(ne)
         		ENDIF
         	ENDIF
         	ENDIF
     		ENDDO

     
     
!initialise the resistances  and update the compliances  
	     DO noelem=1,NEELEM(0)
     		ne=NEELEM(noelem)
     		RESIS_SUM(ne)=CE(nm_R,ne)
     	     ENDDO
     
! Resistance     
 	       DO noelem=NEELEM(0),1,-1
     		ne=NEELEM(noelem)
     		inv_sum_resis=0.0d0
    	 	sum_resis=0.0d0
    	 	IF(NXI(1,0,ne).GT.0)THEN
     		  !note terminal elements will have no children
     		  DO noelem2=1,NXI(1,0,ne) !for each child
             		ne2=NXI(1,noelem2,ne) !get the child element number
             		inv_sum_reis=inv_sum_reis+(1.0d0/RESIS_SUM(ne2))
          	 ENDDO
          	ENDIF 
          	sum_reisis=1.0d0/inv_sum_resis
         	!add the resistances in series to the one above
         	RESIS_SUM(ne)=RESIS_SUM(ne)+sum_resis
     		ENDDO
        
     
     
      CALL EXITS('SUM_RESIS_COMPL')
      RETURN
 9999 CALL ERRORS('SUM_RESIS_COMPL',ERROR)
      CALL EXITS('SUM_RESIS_COMPL')
      RETURN 1
      END
