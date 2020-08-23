      SUBROUTINE SUM_RESIS_COMPL2(BBM,CE,COMP_SUM,
     & NEELEM,NELIST,NXI,NORD,NO_SUBTENDED,
     & RESIS_SUM,undef,VOL_SUM,ERROR,*)


    
    
      INCLUDE 'b00.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'
      
!     Parameter List
      INTEGER NEELEM(0:NE_R_M),NELIST(0:NEM),NORD(5,NE_R_M),
     & NO_SUBTENDED(NE_R_M),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),COMP_SUM(NE_R_M),
     &  RESIS_SUM(NE_R_M),VOL_SUM(NE_R_M),RESIS_SUM2(NE_R_M)
      CHARACTER ERROR*(*)     
!     Local variables
      INTEGER noelem,ne,noelem2,ne2,term,start,
     & subtended,term_ne(NE_R_M),total_iterations,ne_save,
     & TRUNC_TREE(NE_R_M),TERMINALS(NE_R_M),nep,
     & PARENT_LIST(NE_R_M),parent,last_parent,neoelem0,ne0,subscript      
      REAL*8 sum_vol,sum_comp,sum_resis,inv_sum_resis,check,
     & sum_av,c1,c2,r1,r2,Ceff,Re,w,A,B,K,M,undef   
      LOGICAL DONE,GO_ON,IN_PARENT_LIST,MY_STOP
     
      CALL ENTERS('SUM_RESIS_COMPL2',*9999)
!    This subroutine calculates the summed resistance and compliance for 
!    each tree.
!    Note on entry BBM(1,ne) and CW(nm_C,ne) are >0.0d0 only for t
!    terminals
!    The subroutine at each call starts at the terminal branches and using the 
!    compliance and resistance determined for each terminal branch calculates the 
!    correct compliance and resistance sum -it is an iterative method and uses the 
!    algorithm used to truncate the tree by the automatic method. 


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
      
      ! value for omega
      w=2*pi/T_cycle 

!    Have to intialise these variables, they seem to get messy 
   	  Do noelem=1,NEELEM(0)
      	  ne=NEELEM(noelem)
      	  VOL_SUM(ne)=0.0d0
      	  RESIS_SUM(ne)=0.0d0
      	  RESIS_SUM2(ne)=0.0d0
      	  COMP_SUM(ne)=0.0d0
      	  NO_SUBTENDED(ne)=0
      	  term_ne(ne)=0
      	  TRUNC_TREE(ne)=1
      	  TERMINALS(ne)=0
    	 ENDDO   

!intialise the volume/compliance/resistance at the terminals
       term=0
       DO noelem=1,NEELEM(0)
	     	ne=NEELEM(noelem)
		IF(NXI(1,0,ne).EQ.0)THEN
			VOL_SUM(ne)=BBM(1,ne)
			COMP_SUM(ne)=CE(nm_C,ne) !from tissue complaince
			term=term+1
			term_ne(term)=ne
			NO_SUBTENDED(ne)=1
			TERMINALS(ne)=1
			RESIS_SUM(ne)=CE(nm_R,ne) !from branch resistance	
		ENDIF
		RESIS_SUM2(ne)=CE(nm_R,ne)	
       ENDDO

!   Determine the regional volumes and the number subtended by any given unit.
  	       DO noelem=NEELEM(0),1,-1
     		ne=NEELEM(noelem)
    	 	sum_vol=0.0d0
    	 	sum_comp=0.0d0
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
                CE(nm_vinit,ne)=CE(nm_vinit,ne)+sum_vol
         	NO_SUBTENDED(ne)=NO_SUBTENDED(ne)+subtended
         	!only add if sum_comp non zero and dont invert possible zeros
     		ENDDO
       

!Resitance and Complaince are calcualted using iterative grouping of paired terminal
!units and Equations (13) and (14) Otis 1956 JAP  Vol 8 January 427-443
   	DONE=.FALSE.
   	MY_STOP=.FALSE.
   	total_iterations=0	   
   	!use a done loop to calculate the whole tree resistance as opposed to 
   	!iterations which determine the truncation level.
   	DO while(.NOT.DONE)      
   	     total_iterations=total_iterations+1
c#      Modification to deal with truncation of multi-branching trees   
c#      ie if there are refined elements in the tree want the chosen INITIAL TERMINAL element to be the first one after the
c#      junction (later need to review calculation of truncted resistance)     
        DO noelem=1,term
            ne=TERM_NE(noelem)
            ne_save=ne
            GO_ON=.true.
            !look up the tree until the parent has two daughters, that are terminal
            DO while(GO_ON)
           	ne2=NXI(-1,1,ne)
            	IF(NXI(1,0,ne2).eq.2)THEN
            		GO_ON=.false.
		ELSE
			TRUNC_TREE(ne)=0 !remove the terminal elements in a multi-element tree
                 	TERMINALS(ne)=0
			ne=ne2
		ENDIF
		IF(ne_save.ne.ne)THEN
		   TERM_NE(noelem)=ne !update ne in the terminal list	
		   TERMINALS(ne)=1 !update the terminal list
		ENDIF
	    ENDDO
         ENDDO

!CRITERIA FOR BEING A NEW TERMINAL IS THAT IT SUBTENDS 2 TERMINAL BRANCHES
	 !identify branches which will become terminal
	 parent=0
	 last_parent=0
	 DO noelem=1,term
	      ne=TERM_NE(noelem)
	      nep=NXI(-1,1,ne)
	      ne0=NXI(-1,1,nep)
	      IF(nep.eq.1)THEN
	      	!nep is the trachea and this is the last iteration
	      	DONE=.true.
	      ENDIF
	      test=0
	      DO noelem2=1,NXI(1,0,nep) !for each child of the parent

             	 ne2=NXI(1,noelem2,nep) !get the child element number
             	 IF(TERMINALS(ne2).eq.1)THEN
             	 	test=test+1
             	 ENDIF	
             ENDDO

             !create a reduced parent list (works for ordered trees)
             IF(test.eq.2)THEN
                IF(nep.NE.last_parent)THEN
             		parent=parent+1
             		PARENT_LIST(parent)=nep 
             	ENDIF
             	last_parent=nep	
             ENDIF
         ENDDO
          
          !update the terminals and the truncated tree
          !assumption is that in all of the tree below the "truncation point" the volume
          !increases by the smae amount
          DO noelem=1,parent
	    ne=PARENT_LIST(noelem)
            subscript=0
	    DO noelem2=1,NXI(1,0,ne) !for each child of the parent
	         subscript=subscript+1
             	 ne2=NXI(1,noelem2,ne) !get the child element number
                 TRUNC_TREE(ne2)=0
                 TERMINALS(ne2)=0
                 IF(subscript.eq.1)Then
                 	  t1=resis_sum2(ne2)*comp_sum(ne2)
                 	  r1=resis_sum2(ne2)
                 	  c1=comp_sum(ne2)
                 ELSEIF(subscript.eq.2)THEN
                        t2=resis_sum2(ne2)*comp_sum(ne2)
                        r2=resis_sum2(ne2)
                        c2=comp_sum(ne2)
                 ENDIF
             ENDDO
             A=(w**2)*r1*c1*r2*c2-1.0d0
             B=w*(r1*c1+r2*c2)
             K=(w**2)*c2*c1*(r2+r1)
             M=w*(c1+c2)
             !calculate effective compliance(Otis(13))
             Ceff=(K**2+M**2)/(w*(B*K-A*M))
             COMP_SUM(ne)=Ceff
             !calculate effective resistance(Otis(14)) of subtended
             !add the effective resistances in parallel
             Re=(A*K+B*M)/(K**2+M**2)
             RESIS_SUM(ne)=Re+RESIS_SUM(ne)
             RESIS_SUM2(ne)=Re+RESIS_SUM2(ne)
             TERMINALS(ne)=1
          ENDDO 
	 
	  !update the terminal list
	  DO noelem=1,term
	      TERM_NE(noelem)=0
          ENDDO
          term=0
       	   DO noelem=1,NEELEM(0)
	     ne=NEELEM(noelem)
	     IF(TERMINALS(ne).eq.1)THEN
	     	term=term+1
	     	TERM_NE(term)=ne
	     ENDIF	
	  ENDDO
	 
	  IF(term.eq.2)THEN
	  	DONE=.true.
	  ENDIF	
	ENDDO !not done 
        
        open(98,file='comp_resis_Otis.txt')
         Do noelem=1,NEELEM(0)
      	  ne=NEELEM(noelem)
      	  write(98,*)'ne=',ne,'  R=',resis_sum(ne),'   C=',comp_sum(ne)
      	 ENDDO
      	close(98)  
     
     
      CALL EXITS('SUM_RESIS_COMPL2')
      RETURN
 9999 CALL ERRORS('SUM_RESIS_COMPL2',ERROR)
      CALL EXITS('SUM_RESIS_COMPL2')
      RETURN 1
      END
