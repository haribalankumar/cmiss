      SUBROUTINE TRUNCATE(NBJ,NDP,NEELEM,NELIST,NENP,NORD,NPNE,
     &  NPNODE,NRLIST,NVJE,NVJP,NXI,BBM,CE,XAB,XP,ZD,ERROR,*)
     

C#### Subroutine: TRUNCATE
C###  Description:
C###  Truncate will truncate a one dimensional airway mesh
C###  Calculate the lumped compliance and resistance for each level
C###  of truncation
C###  Created by Jennine Mitchell Feb 2011...Note same lists as EVFLOW.f

C###  fem truncate

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn' !MEM_INIT
      INCLUDE 'geom00.cmn' 
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc' !DP_TYPE
      INCLUDE 'valu00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NDP(NDM),NEELEM(0:NE_R_M,0:NRM),
     &  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM,NXM),XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER N3CO,noelem,IFROMC,nr,ne,MECHANICS_FILETYPE,
     & NJJ_FLOW,
     & NJJ_VRATIO,NO_SUBTENDED(NE_R_M),no_iterations,test,
     & njj_truncate,TERMINALS(NE_R_M),noelem2,ne2,nep,
     & terminal_daughters,my_it,nv1,nv2,nb,np,np2,term,
     & TERM_LIST(NE_R_M),parent,TRUNC_TREE(NE_R_M),
     & PARENT_LIST(NE_R_M),last_parent,noelem0,ne_save,
     & asym_match,test_above,test1,test2,start_gen,kount,
     & end_gen,sym_resisG,sym_resisT,asym_resis

      REAL*8 RFROMC,CW,sumvolume,undef,VINIT(NEM),
     & RESIS_SUM(NE_R_M),VOL_SUM(NE_R_M),COMP_SUM(NE_R_M),
     & volume,ref_ratio,vol_equiv,ratio,lambda,lambda15,
     & exp_term,comp,cc,min_E,a,b,c,dpl
      LOGICAL CBBREV,READ_VOLUMES,ITERATE,IN_PARENT_LIST,
     & GO_ON,FIRST 
     
      CALL ENTERS('TRUNCATE',*9999)
c      #intended options
c      (1)
c      # truncate iteratively ntimes 
c      #this is truncation so that if there are two daughters they get 'merged'
c      #this truncation continues iteratively
c      (2)
c      #truncate by field ie if the tree is not to be truncated, then use
c      # a signal in field X 0/1 ie 0 truncate 1 retain (allows to truncate
c      #back to MDCT tree (the truncation will still be based on the iterative)
c      #truncation routine
c      #Note have imported NORD but not considering symmetry at current time 
c      ##ALSO sets a logical TRUNCATED to tell the EVFLOW & subsidary 
c      #subroutines that the mesh has been truncated..this will be used
c      #to govern the recalculation of volumes/compliance and resitance
c      #that would otherwise occur
      
      !set a logical TRUNCATED which will in EVFLOW.f and the 
      !subroutines that it calls govern the calculations of volume
      !tissue compliance and airway resistance

      !indices for CE arrays
      nm_R=1
      nm_C=2
      nm_Ppl=3
      nm_Rt=4
      nm_lambda=5
      nm_dPl=6
      nm_volumes=7
      nm_vinit=8
      nm_vol_bel=9
      !was 6,7,8,9,10
      nm_R_trunc=10
      nm_C_trunc=11
      nm_ST_trunc=12
      nm_C_av=13
      nm_R_av=14

      !same defaults as EVFLOW.f
      NJJ_FLOW=3 !default
      nj_flow=NJ_LOC(NJL_FIEL,NJJ_FLOW,nr) !airway flow


      IF(CBBREV(CO,'region',3,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
      ELSE
          nr=1
      ENDIF

      CW=0.2d0
      
      IF(CBBREV(CO,'READ_VOLUMES',4,noco+1,NTCO,N3CO)) THEN
        READ_VOLUMES=.TRUE.
        IF(CBBREV(CO,'UNDEFORMED',3,noco+1,NTCO,N3CO)) THEN
          undef=RFROMC(CO(N3CO+1))
        ELSE
          undef=1.d0
        ENDIF

        MECHANICS_FILETYPE=1 !default
        ! nj_Vratio=5 !default
        IF(CBBREV(CO,'DATAPOINT',3,noco+1,NTCO,N3CO))THEN
          nj_Vratio=IFROMC(CO(N3CO+1))
          MECHANICS_FILETYPE=1
        ELSEIF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO))THEN
          NJJ_VRATIO=IFROMC(CO(N3CO+1))
          nj_Vratio=NJ_LOC(NJL_FIEL,NJJ_VRATIO,nr) !deformed to undeformed acinar volume ratio
          MECHANICS_FILETYPE=2
        ENDIF
      ENDIF

      IF(CBBREV(CO,'LUMPED_SYMMETRY',3,noco+1,NTCO,N3CO)) THEN
          LUMPED_SYMMETRY=.true.
      ELSE
          LUMPED_SYMMETRY=.false.
      ENDIF
      
      IF(CBBREV(CO,'TIME_DEP_SOLN',3,noco+1,NTCO,N3CO)) THEN
          TIME_DEP_SOLN=.true.
      ELSE
          TIME_DEP_SOLN=.false.
      ENDIF
      
        IF(CBBREV(CO,'ITERATE',3,noco+1,NTCO,N3CO))THEN
          no_iterations=IFROMC(CO(N3CO+1))
          njj_truncate=IFROMC(CO(N3CO+2))
          nj_truncate=NJ_LOC(NJL_FIEL,NJJ_TRUNCATE,nr)
          write(*,*)'nj_truncate=',nj_truncate !tick
          ITERATE=.TRUE.
        ELSEIF(CBBREV(CO,'TR_FIELD',3,noco+1,NTCO,N3CO))THEN
          ITERATE=.FALSE.
           njj_truncate=IFROMC(CO(N3CO+1))
           write(*,*)'njj_truncate=',njj_truncate !tick
          nj_truncate=NJ_LOC(NJL_FIEL,NJJ_TRUNCATE,nr)
        ENDIF



        IF(CBBREV(CO,'T_CYCLE',3,noco+1,NTCO,N3CO)) THEN
        T_CYCLE=RFROMC(CO(N3CO+1))
        ENDIF

c      !read in the volume of each terminal unit
c      !note BBM(1,ne)>0.0d0 only for terminals
          CALL READINITIALVOLUME(MECHANICS_FILETYPE,NBJ,NDP,
     &      NEELEM(0,nr),NENP(1,1,1),NPNE,NXI,BBM,
     &      CE,sumvolume,undef,XP,ZD,READ_VOLUMES,ERROR,*9999)

c     !calculate the resistance
c   	  !.......XP(nj=9) == Poiseulle/Pedley resistance
       CALL BRANCHRESISTANCE(NBJ,NEELEM(0,nr),NPNE,NVJE,
     &  CE(1,1,1),XP,ERROR,*9999)

c      !calculate the tissue compliance note CE(nmC,ne)>0.0d0 only for terminals
      CALL TISSUECOMPLIANCE(NEELEM(0,nr),NELIST,NPNE,NXI,BBM,
     & CE(1,1,1),CW,dPl,undef,FIRST,ERROR,*9999)

	
c     !sum resistances and complainces throughout the tree
c     !vol_sum, comp_sum refers to the sum of the compliances and volumes
c     !subtending a given branch
c     !resis_sum is the resistance in a given branch plus the resistance
c     !of any branches below it

      IF(LUMPED_SYMMETRY)THEN
        CALL SUM_RESIS_COMPL(BBM,CE,COMP_SUM,CW,
     & NEELEM(0,nr),NELIST,NXI,NORD,NO_SUBTENDED,
     & RESIS_SUM,VOL_SUM,undef,ERROR,*9999)
       ELSE
c     calculate the resitance and compliance before truncation
      	CALL SUM_RESIS_COMPL2(BBM,CE,COMP_SUM,NEELEM(0,nr),
     & NELIST,NXI,NORD,NO_SUBTENDED,RESIS_SUM,undef,
     & VOL_SUM,ERROR,*9999)
      ENDIF


c    !need a binary tree map (TRUNC_TREE(elem)=(1/0) 1 means remains in
c      truncated model
c    !to identify the terminal elements for each tree TERMINALS(ne)=1 
c      identies a terminal branch, TERM_LIST(term)=ne the elements that are 
c      terminal

c	!initialise
	DO noelem=1,NEELEM(0,nr)
	     	ne=NEELEM(noelem,nr)
	     	TERMINALS(ne)=0
	     	TRUNC_TREE(ne)=1
	ENDDO
c	!set the terminals up to be 1
	term=0
	 DO noelem=1,NEELEM(0,nr)
	     	ne=NEELEM(noelem,nr)
		IF(NXI(1,0,ne).EQ.0)THEN
        		TERMINALS(ne)=1
        		term=term+1
        		TERM_LIST(term)=ne
		ENDIF	
       	ENDDO 

     
	IF(ITERATE)THEN
        write(*,*)'iteratively defined truncation field'
	DO my_it=1,no_iterations
c#      Modification to deal with truncation of multi-branching trees   
c#      ie if there are refined elements in the tree want the chosen INITIAL TERMINAL element to be the first one after the
c#      junction (later need to review calculation of truncted resistance)     
        DO noelem=1,term
            ne=TERM_LIST(noelem)
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
		   TERM_LIST(noelem)=ne !update ne in the terminal list	
		   TERMINALS(ne)=1 !update the terminal list
		ENDIF
	    ENDDO
         ENDDO
         
	!CRITERIA FOR BEING A TERMINAL BRANCH IS THAT
	 !identify branches which will become terminal
	 parent=0
	 last_parent=0
	 DO noelem=1,term
	      ne=TERM_LIST(noelem)
	      nep=NXI(-1,1,ne)
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
          DO noelem=1,parent
	    ne=PARENT_LIST(noelem)
	    DO noelem2=1,NXI(1,0,ne) !for each child of the parent
             	 ne2=NXI(1,noelem2,ne) !get the child element number
                 TRUNC_TREE(ne2)=0
                 TERMINALS(ne2)=0
             ENDDO
             TERMINALS(ne)=1
          ENDDO 
	 
	  !update the terminal list
	  DO noelem=1,term
	      TERM_LIST(noelem)=0
          ENDDO
          term=0
          
       	   DO noelem=1,NEELEM(0,nr)
	     ne=NEELEM(noelem,nr)
	     IF(TERMINALS(ne).eq.1)THEN
	     	term=term+1
	     	TERM_LIST(term)=ne
	     ENDIF	
	  ENDDO
	       	
       		write(*,*)
       		write(*,*)'XXXXXXX'
       		write(*,*)'Iteration ',my_it
       		write(*,*)'Number of terminal units =',term

	ENDDO !my_it
	
	ELSE !user defined truncation field
	  write(*,*)'Using user define truncation field',nj_truncate
	 !update trunc tree from field
	 Do noelem=1,neelem(0,nr)
       		ne=NEELEM(noelem,nr)
       		nb=NBJ(1,ne)
       		nv1=NVJE(1,nb,nj_truncate,ne)
        	np2=NPNE(2,nb,ne)
        	TERMINALS(ne)=0
        	!if the terminal node in the truncation tree is zero
        	!the TRUNC_TREE array (elem based) is zero
        	IF(XP(1,1,nj_truncate,np2).lt.1.0d-3)THEN
        		TRUNC_TREE(ne)=0
        	ELSE
        		!write(*,*)'Truncated tree contains element', ne
        		TRUNC_TREE(ne)=1	
        	ENDIF	
         ENDDO
       
c	!set the terminals up to be 1
	term=0
	 DO noelem=1,NEELEM(0,nr)
	     	ne=NEELEM(noelem,nr)
	     	ne2=NXI(1,1,ne) !get the child element number
	     	IF(TRUNC_TREE(ne).eq.1)THEN
	     	    IF(TRUNC_TREE(ne2).eq.0)THEN
        		TERMINALS(ne)=1
        		term=term+1
        		TERM_LIST(term)=ne
        	     ENDIF	
		ENDIF	
       ENDDO
       ENDIF !iterate or field truncation
       
       write(*,*)'There are ',term,'terminals'
       !This writes out a list of the terminal elements in the truncated tree
       open(67,file='terminal_elem_lumped.txt')
       term=0
	 DO noelem=1,NEELEM(0,nr)
	     	ne=NEELEM(noelem,nr)
	     	IF(TERMINALS(ne).eq.1)THEN
        		term=term+1
        		write(67,*)ne
		ENDIF
	 TERM_GEN=term		
       ENDDO
       close(67)

       !update field versions in binary truncation field 
       Do noelem=1,neelem(0,nr)
       	ne=NEELEM(noelem,nr)
       	IF(TRUNC_TREE(ne).eq.0)THEN
       	nb=NBJ(1,ne)
        nv1=NVJE(1,nb,nj_truncate,ne)
        nv2=NVJE(1,nb,nj_truncate,ne)
        np=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
       	XP(1,1,nj_truncate,np2)=0.0d0
        XP(1,1,nj_truncate,np2)=0.0d0
        nep=NXI(-1,1,ne)
            DO noelem2=1,NXI(1,0,nep) !for each daughter branch
            ne2=NXI(1,noelem2,nep) !the daughter element number
            nb=NBJ(1,ne2)
            np2=NPNE(2,nb,ne2)
            XP(1,noelem2+1,nj_truncate,np)=XP(1,1,nj_truncate,np2)
          ENDDO
        ENDIF
       ENDDO	
       

       !intialise
        Do noelem=1,neelem(0,nr)
       		ne=NEELEM(noelem,nr)
       		NELIST(ne)=0
       		CE(nm_R_trunc,ne,1)=0.0d0
       		CE(nm_C_trunc,ne,1)=0.0d0
       	ENDDO
       	
       	
       	Do noelem=1,neelem(0,nr)
       	  ne=NEELEM(noelem,nr)
          IF(TRUNC_TREE(ne).eq.1)THEN
        	NELIST(ne)=1
	  ENDIF   
        ENDDO
       	
       	       
       IF(LUMPED_SYMMETRY)THEN
       !calculate the compliance and resistance after the truncation are identified
        CALL SUM_RESIS_COMPL(BBM,CE,COMP_SUM,CW,NEELEM(0,nr),
     & NELIST,NXI,NORD,NO_SUBTENDED,RESIS_SUM,
     & VOL_SUM,undef,ERROR,*9999)
       ENDIF

        !initialise      
        Do noelem=1,neelem(0,nr)
       	    ne=NEELEM(noelem,nr)
            BBM(1,ne)=0.0d0
        ENDDO

       noelem0=0
       term=0
       Do noelem=1,neelem(0,nr)
       	ne=NEELEM(noelem,nr)
       	!calculate standard Poiseulle Resistance
        IF(TRUNC_TREE(ne).eq.1)THEN
        	noelem0=noelem0+1
        	CE(nm_R_trunc,ne,1)=CE(nm_R,ne,1)
        	!NELIST(ne)=1
        	!overwrite with lumped resistance at terminals
        	!store the compliance
        	!set the daughter numbers to be zero
        	IF(TERMINALS(ne).eq.1)THEN
              		term=term+1
              		CE(nm_R_trunc,ne,1)=RESIS_SUM(ne) 
              		CE(nm_C_trunc,ne,1)=COMP_SUM(ne) 
              		CE(nm_ST_trunc,ne,1)=NO_SUBTENDED(ne) 
              		BBM(1,ne)=VOL_SUM(ne)
                ENDIF
	ENDIF   
       ENDDO
       
c     !set a logical that the truncation is done
c     !truncated tree thereafter.    
      
      TRUNCATED=.TRUE.
       volume=0.0d0
       Do noelem=1,neelem(0,nr)
       	 ne=NEELEM(noelem,nr)
             IF(TERMINALS(ne).eq.1)THEN
              	volume=volume+BBM(1,ne)
             ENDIF
       ENDDO
    

      CALL EXITS('TRUNCATE')
      RETURN
 9999 CALL ERRORS('TRUNCATE',ERROR)
      CALL EXITS('TRUNCATE')
      RETURN 1
      END



