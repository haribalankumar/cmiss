      SUBROUTINE CALC_TRANSIT_TIME(nb,NELIST,NENP,NPNE,NYNE,
     '  NYNP,CE,PATHWAY,PLOT_DATA,TIME,YP,ERROR,*)

C### Subroutine: CALC_TRANSIT_TIME
C### Description:
C###   CALC_TRANSIT_TIME calculates the red blood cell (RBC)
C###   transit time distribution through a pulmonary
C###   capillary network. Data for randomly selected pathways is
C###   calculated & the probability of flow in that path is used to
C###   weight the average values.

C**** Created by Kelly Burrowes, June 2002
C***  Method based on work by Huang et al, 2001.

C***  PATHWAY(path#,1) stores last element number of pathway
C***  PATHWAY(path#,2) stores sum of RBC transit time
C***  PATHWAY(path#,3) stores product of the RBC volume fraction in ne
C***  PATHWAY(path#,4) stores sum of the length of the pathway
C***  PATHWAY(path#,5) stores the number of segments in a pathway
C***  PATHWAY(path#,6) stores WBC transit time through path
C***  PATHWAY(path#,7) stores # time WBC stuck in pathway
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER nb,NELIST(0:NEM),NENP(NPM,0:NEPM),
     '  NPNE(NNM,NBFM,NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 CE(NMM,NEM),PATHWAY(0:MAX_PATH,FACTORS),
     '  PLOT_DATA(4500,FACTORS),TIME(NEM),YP(NYM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,ISEED,j,MAX_PATH_I,ne,ne2,ne3,NO_INT,noelem,noelem2,
     '  np1,np2,np3,ny,path,tint,TOTAL_PATHS,TOTAL_PATHS_I
      REAL*8 AVGE_WBC,CM_RANDOM_NUMBER,c_tt,flow,flow_in,
     '  flow_segment,high,length,length_avge,low,MIN_WBC,
     '  MAX_WBC,NO_BLOCK,NO_STUCK,NO_STUCK_PATH,num_segments,P1,P2,P_cr,
     '  R_vessel,R_WBC,RBC_fraction,RBC_time,SD,SD_WBC,SIZE_INT,stuck,t,
     '  tau,te,tol,total_flow,tp,volume,WBC_TE,WBC_time2,WBC_TP
      LOGICAL CONTINU,FOUND,OUT

      CALL ENTERS('CALC_TRANSIT_TIME',*9999)

      DO noelem=1,NEM
        TIME(noelem)=0.d0 !initialise
      ENDDO
C... finding total flow into network
      total_flow=0.d0
      DO i=1,INLETS(0)
        ne=INLETS(i)
        total_flow=total_flow+DABS(YP(NYNE(1,1,0,1,ne)))*CE(nm_Hd,ne)
      ENDDO
C... calculating RBC transit time in each capillary segment
      c_tt=1.4d0 !ratio of RBC/plasma transit times (ref:Presson)
      DO noelem=1,NELIST(0)
        ne=NELIST(noelem)
        volume=PI*CE(nm_a,ne)*CE(nm_b,ne)*CE(nm_length,ne)!*1000.d0**3.d0 !mm^3->um^3
        ny=NYNE(1,1,0,1,ne) !ny #
C        CE(nm_tt,ne)=DABS(volume/(YP(ny)*c_tt)) !RBC transit time in ne
        TIME(ne)=DABS(volume/(YP(ny)*c_tt)) !RBC transit time in ne
      ENDDO !noelem
C... always start with inlet element, randomly chose path from here
      rbc_time=0.d0
      TOTAL_PATHS=0
      DO i=0,MAX_PATH !initialise
        DO j=1,FACTORS !stores required factors for a given pathway
          PATHWAY(i,j)=0.d0
        ENDDO
      ENDDO
      DO i=1,INLETS(0)
        WRITE(OP_STRING,'('' Inlet element # = '',I5)')INLETS(i)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ISEED=2
        NO_BLOCK=0.d0
        AVGE_WBC=0.d0 !initialise
        DO j=1,FACTORS
          PATHWAY(0,j)=0.d0 !reseting totals
        ENDDO
        MAX_PATH_I=INT(MAX_PATH*(DABS(YP(NYNE(1,1,0,1,INLETS(i)))
     '    *CE(nm_Hd,INLETS(i))))/total_flow) !MAX_PATH*frac RBC into each inlet
        TOTAL_PATHS_I=0
        DO WHILE(TOTAL_PATHS_I.LT.MAX_PATH_I)
          TOTAL_PATHS=TOTAL_PATHS+1 !new pathway
          TOTAL_PATHS_I=TOTAL_PATHS_I+1 !pathway # for inlet(i)
          NO_STUCK=0.d0 !initialise # WBCs stuck
C... randomly give WBC diameter from range 6.8-8.3 um (Doerschuk:)
          MAX_WBC=6.8d-3/2.d0 !mm
          MIN_WBC=8.3d-3/2.d0 !mm
          R_WBC=CM_RANDOM_NUMBER(ISEED)*(MAX_WBC-MIN_WBC)+MIN_WBC
          AVGE_WBC=AVGE_WBC+R_WBC !sums WBC radius to find average
C... Rc=3.5um (Huang, et al, 2001), diameter=6.8-8.3um Doerschuk, 1999.
C... don't want time/length data from inlet or outlet vessels
          PATHWAY(TOTAL_PATHS,1)=INLETS(i) !multiple inlets
          PATHWAY(TOTAL_PATHS,3)=1.d0 !always =1 for inlet
          PATHWAY(TOTAL_PATHS,5)=0
          CONTINU=.TRUE.
          DO WHILE(CONTINU)
C... creating 1 complete pathway at a time, randomly select next flow ne
            ne=INT(PATHWAY(TOTAL_PATHS,1))
            np1=NPNE(1,nb,ne)
            np2=NPNE(2,nb,ne)
            P1=YP(NYNP(1,1,1,np1)) !find out direction of flow
            P2=YP(NYNP(1,1,1,np2))
            IF(P1.LT.P2) np2=np1 !else np2=np2
            FOUND=.FALSE.
            num=0
            DO WHILE(.NOT.FOUND.AND.num.LT.50)
              num=num+1 
              noelem=INT(CM_RANDOM_NUMBER(ISEED)*NENP(np2,0))+1 !random
              ne2=NENP(np2,noelem) !next element in pathway
              np1=NPNE(2,nb,ne2) !check if flow into ne2 possible
              IF(np1.EQ.np2) np1=NPNE(1,nb,ne2)
              P1=YP(NYNP(1,1,1,np1))
              IF(P1.LT.P2.AND.ne2.NE.ne.AND.CE(nm_Hd,ne2).GT.0.d0)
     '          FOUND=.TRUE. !ne2 in pathway
            ENDDO !WHILE
            IF(FOUND) THEN !(tau) changing dyn/cm ->N/m
              tau=0.035d0*1.d-3 !(N/m) avge tension in cell cortex
C... effective circular vessel radius
C             R_vessel=DSQRT(CE(nm_a,ne2)*CE(nm_b,ne2)) !mm
C             R_vessel=CE(nm_b,ne2) ! =h/2 (height of capillary sheet)
              R_vessel=(DSQRT(CE(nm_a,ne2)*CE(nm_b,ne2))+CE(nm_b,ne2))
     '          /2.d0 !average of min & max vessel approximation
C... critical pressure (Pa) for WBC entrance to vessel (units N/m/m=Pa)
              P_cr=(2.d0*tau/(R_WBC/1000.d0))*((R_WBC/R_vessel)-1.d0) !Pa
              tol=1.d-4
 !pressure must exceed P_cr+tolerance, otherwise get very long entrance time
              IF(P2-P1.GT.P_cr+tol) THEN !critical pressure exceeded
                length=CE(nm_length,ne2) !vessel length
                flow=DABS(YP(NYNE(1,1,0,1,ne2))) !flow in segment
                CALL WBC_TRANSIT(flow,length,P1,P2,P_cr,R_WBC,
     '            R_vessel,t,te,tp,ERROR,*9999)
                PATHWAY(TOTAL_PATHS,8)=PATHWAY(TOTAL_PATHS,8)+tp !passing time
                PATHWAY(TOTAL_PATHS,9)=PATHWAY(TOTAL_PATHS,9)+te !entrance time
                PATHWAY(TOTAL_PATHS,6)=PATHWAY(TOTAL_PATHS,6)+t
C... neutrophil is defined as "stopped" if entrance time(te) is >= 1 sec
                IF(te.GE.1.d0) NO_STUCK=NO_STUCK+1.d0 !will squeeze through
              ELSE
C... make resistance in this element (ne2) infinite, WBC blocks ne2
C... or more effectively - remove element from solution
                NO_STUCK=NO_STUCK+1.d0 !# times WBC is stuck in pathway
                NO_BLOCK=NO_BLOCK+1.d0 !# which stuck in network
              ENDIF
              flow_in=0.d0 !initialise
              DO noelem2=1,NENP(np2,0)
                ne3=NENP(np2,noelem2)
                IF(ne3.NE.ne2) THEN !ne3 isn't ne being evaluated
                  np3=NPNE(1,nb,ne3)
                  IF(np2.EQ.np3) np3=NPNE(2,nb,ne3) !other node
C... checks if flow is into junction
                  P1=YP(NYNP(1,1,1,np3))
                  P2=YP(NYNP(1,1,1,np2))
                  IF(P2.LT.P1)
     '              flow_in=flow_in+DABS(YP(NYNE(1,1,0,1,ne3)))
C                 '  * CE(nm_Hd,ne3) !(Hd=hematocrit) vol flow of RBC
C.. currently not sure whether RBC flow or total blood flow used
                ENDIF
              ENDDO !noelem2
              flow_segment=YP(NYNE(1,1,0,1,ne2)) !*CE(nm_Hd,ne2)
              RBC_fraction=DABS(flow_segment/flow_in) !frac in ne
              PATHWAY(TOTAL_PATHS,1)=ne2 !last ne # of pathway
              PATHWAY(TOTAL_PATHS,3)=PATHWAY(TOTAL_PATHS,3)*RBC_fraction
              OUT=.FALSE.
              DO j=1,OUTLETS(0)
                IF(ne2.EQ.OUTLETS(j)) THEN
                  OUT=.TRUE. !ne2=outlet vessel
                  CONTINU=.FALSE.
                ENDIF
              ENDDO !j
              IF(.NOT.OUT) THEN !don't add outlet vessel data
                PATHWAY(TOTAL_PATHS,2)=PATHWAY(TOTAL_PATHS,2)+
     &            TIME(ne2)
                PATHWAY(TOTAL_PATHS,4)=PATHWAY(TOTAL_PATHS,4)+
     '            CE(nm_length,ne2)
                PATHWAY(TOTAL_PATHS,5)=PATHWAY(TOTAL_PATHS,5)+1.d0
              ENDIF
            ELSE !IF(.NOT.FOUND)
              CONTINU=.FALSE.
              TOTAL_PATHS=TOTAL_PATHS-1 !pathway not completed, start again
            ENDIF !FOUND
          ENDDO !CONTINU
          PATHWAY(TOTAL_PATHS,7)=NO_STUCK
        ENDDO !WHILE
C... OPTCAP prints results to screen
        CALL OPTCAP(i,MAX_PATH_I,TOTAL_PATHS,AVGE_WBC,
     '    NO_BLOCK,PATHWAY,ERROR,*9999)
      ENDDO !INLETS
C... write transit times for each pathway to a file, so can plot results
      SIZE_INT=0.010d0
      low=0.d0
      high=45.d0
      NO_INT=INT((high-low)/SIZE_INT)
      IF(NO_INT.GT.4500) THEN
        WRITE(OP_STRING,'('' Increase # intervals to: '',I5)')
     '    NO_INT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      DO tint=1,NO_INT
        DO i=1,FACTORS
          PLOT_DATA(tint,i)=0.d0 !initialise
        ENDDO
      ENDDO
      PATHWAY(0,3)=0.d0
      PATHWAY(0,7)=0.d0
      NO_STUCK_PATH=0.d0
      DO path=1,TOTAL_PATHS
        PATHWAY(0,3)=PATHWAY(0,3)+PATHWAY(path,3) !RBC fraction
        PATHWAY(0,7)=PATHWAY(0,7)+PATHWAY(path,7) !# times WBC stuck in path
        IF(PATHWAY(path,7).GT.0.d0) NO_STUCK_PATH=NO_STUCK_PATH+1.d0
        DO tint=1,NO_INT
          IF(PATHWAY(path,2).GE.(low+(tint-1)*SIZE_INT).AND. !RBC
     '      PATHWAY(path,2).LT.(low+tint*SIZE_INT)) THEN !inside interval
            PLOT_DATA(tint,1)=PLOT_DATA(tint,1)+PATHWAY(path,3)!frac Q
            PLOT_DATA(tint,2)=PLOT_DATA(tint,2)+1.d0 !# paths in inter
          ENDIF
          IF(PATHWAY(path,6).GE.(low+(tint-1)*SIZE_INT).AND. !WBC
     '      PATHWAY(path,6).LT.(low+tint*SIZE_INT)) THEN
            PLOT_DATA(tint,4)=PLOT_DATA(tint,4)+PATHWAY(path,3) !frac Q
            PLOT_DATA(tint,5)=PLOT_DATA(tint,5)+1.d0 !# paths in inter
          ENDIF
          !wbc results
          IF(PATHWAY(path,8).GE.(low+(tint-1)*SIZE_INT).AND. !WBC tp
     '      PATHWAY(path,8).LT.(low+tint*SIZE_INT)) THEN
            PLOT_DATA(tint,8)=PLOT_DATA(tint,8)+1.d0 !# paths in inter
          ENDIF
          IF(PATHWAY(path,9).GE.(low+(tint-1)*SIZE_INT).AND. !WBC te
     '      PATHWAY(path,9).LT.(low+tint*SIZE_INT)) THEN
            PLOT_DATA(tint,9)=PLOT_DATA(tint,9)+1.d0 !# paths in inter
          ENDIF
        ENDDO !tint
        IF(PATHWAY(path,7).EQ.0) THEN !# times WBC have been stopped
          PLOT_DATA(1,3)=PLOT_DATA(1,3)+1 !in a pathway STUCK ZERO TIMES
        ELSE IF(PATHWAY(path,7).EQ.1) THEN
          PLOT_DATA(2,3)=PLOT_DATA(2,3)+1 !STUCK ONCE
        ELSE IF(PATHWAY(path,7).EQ.2) THEN
          PLOT_DATA(3,3)=PLOT_DATA(3,3)+1 !STUCK TWICE
        ELSE IF(PATHWAY(path,7).GE.3) THEN
          PLOT_DATA(4,3)=PLOT_DATA(4,3)+1 !STUCK THREE OR MORE TIMES
        ENDIF
      ENDDO !path
      SD=0.d0 !this calculates a sample standard deviation
      SD_WBC=0.d0
      RBC_time=0.d0
      WBC_time2=0.d0
      WBC_TE=0.d0
      WBC_TP=0.d0
      num_segments=0.d0
      length_avge=0.d0
      DO path=1,TOTAL_PATHS
        RBC_time=RBC_time+PATHWAY(path,2)*(PATHWAY(path,3)/PATHWAY(0,3))
        WBC_time2=WBC_time2+PATHWAY(path,6)*(PATHWAY(path,3)/
     '    PATHWAY(0,3))
        WBC_TE=WBC_TE+PATHWAY(path,9)*(PATHWAY(path,3)/
     '    PATHWAY(0,3))
        WBC_TP=WBC_TP+PATHWAY(path,8)*(PATHWAY(path,3)/
     '    PATHWAY(0,3))
        num_segments=num_segments+PATHWAY(path,5)
        length_avge=length_avge+PATHWAY(path,4)
      ENDDO
      DO path=1,TOTAL_PATHS
        SD=SD+(PATHWAY(path,2)-RBC_time)**2.d0
        SD_WBC=SD_WBC+(PATHWAY(path,6)-WBC_time2)**2.d0
      ENDDO
      SD=DSQRT(SD/(TOTAL_PATHS-1))
      SD_WBC=DSQRT(SD_WBC/(TOTAL_PATHS-1))
      WRITE(OP_STRING,'(''Standard Deviation (RBC),(WBC)= '',2(F12.4))')
     '  SD,SD_WBC
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      num_segments=num_segments/TOTAL_PATHS
      length_avge=length_avge/TOTAL_PATHS
      WRITE(OP_STRING,'('' Overall average RBC transit time '//
     '  '(s)= '',F12.4)') RBC_time
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Overall average # segments in a pathway '//
     '  '= '',F12.4)') num_segments
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Overall average WBC transit time (s) '//
     '  '= '',F12.4)') WBC_time2
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Overall average WBC entrance time (s) '//
     '  '= '',F12.4)') WBC_TE
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Overall average WBC passage time (s) '//
     '  '= '',F12.4)') WBC_TP      
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Overall average path length (mm) '//
     '  '= '',F12.4)') length_avge
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      PATHWAY(0,7)=PATHWAY(0,7)/TOTAL_PATHS
      WRITE(OP_STRING,'('' Average # times WBC stop in path= '',F12.4)')
     '  PATHWAY(0,7)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      stuck=(NO_STUCK_PATH/TOTAL_PATHS)*100.d0
      WRITE(OP_STRING,'('' % pathways WBCs get stopped in = '',F12.4)')
     '  stuck
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' % of total RBCs in all pathways: '',D12.4)')
     '  PATHWAY(0,3)*100.d0 !for test
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      PLOT_DATA(1,3)=(PLOT_DATA(1,3)/TOTAL_PATHS)*100.d0
      PLOT_DATA(2,3)=(PLOT_DATA(2,3)/TOTAL_PATHS)*100.d0
      PLOT_DATA(3,3)=(PLOT_DATA(3,3)/TOTAL_PATHS)*100.d0
      PLOT_DATA(4,3)=(PLOT_DATA(4,3)/TOTAL_PATHS)*100.d0      
      CALL OPENF(IOFILE2,'DISK','transit2.out','NEW',
     '  'SEQUEN','FORMATTED',132,ERROR,*9999)
C      WRITE(OP_STRING,
C     '  '('' TIME INTERVAL | % RBC FLOW | % RBC PATHWAYS |'//
C     '  ' % WBC FLOW | % WBC PATHWAYS '')')
C      CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)
      DO tint=1,NO_INT
        PLOT_DATA(tint,1)=(PLOT_DATA(tint,1)/PATHWAY(0,3))*100.d0
        !RBC: % of sampled flow in this time interval
        PLOT_DATA(tint,2)=(PLOT_DATA(tint,2)/TOTAL_PATHS)*100.d0
        !RBC: % of sampled pathways in this time interval
        PLOT_DATA(tint,4)=(PLOT_DATA(tint,4)/PATHWAY(0,3))*100.d0
        !WBC: % of sampled flow in this time interval
        PLOT_DATA(tint,5)=(PLOT_DATA(tint,5)/TOTAL_PATHS)*100.d0
        !WBC: % of sampled pathways in this time interval
        WRITE(IOFILE2,'(3X,5(F12.6))') (tint-1+tint)*SIZE_INT*0.5d0,
     '    PLOT_DATA(tint,1),PLOT_DATA(tint,2),PLOT_DATA(tint,4),
     '    PLOT_DATA(tint,5) !x @ center of interval
      ENDDO
      CALL CLOSEF(IOFILE2,ERROR,*9999)
      CALL OPENF(IOFILE3,'DISK','transit3.out','NEW',
     '  'SEQUEN','FORMATTED',132,ERROR,*9999)
C      WRITE(OP_STRING,
C     '  '('' TIME INTERVAL | WBC STUCK ONCE | TWICE | THRICE '')')
C      CALL WRITES(IOFILE3,OP_STRING,ERROR,*9999)
      !This is the # of times WBC gets stuck
      DO i=1,4
        WRITE(IOFILE3,'(3X,2(F12.6))') (i-1+i)*0.5d0,PLOT_DATA(i,3)
      ENDDO 
      CALL CLOSEF(IOFILE3,ERROR,*9999)

!testing WBC TRANSIT TIMES
      CALL OPENF(IOFILE2,'DISK','wbc_test.out','NEW','SEQUEN',
     '  'FORMATTED',132,ERROR,*9999)
C      DO path=1,TOTAL_PATHS
C        IF(PATHWAY(path,6).GT.100.d0) THEN
C          WRITE(IOFILE2,'(2X,I7,7(F12.4))') path,PATHWAY(path,6),
C     '      PATHWAY(PATH,1),PATHWAY(PATH,2),PATHWAY(PATH,3),
C     '      PATHWAY(PATH,4),PATHWAY(PATH,5),PATHWAY(PATH,7)
C        ENDIF
C     ENDDO

C      DO path=1,TOTAL_PATHS
      DO tint=1,NO_INT
        PLOT_DATA(tint,8)=(PLOT_DATA(tint,8)/TOTAL_PATHS)*100.d0
        PLOT_DATA(tint,9)=(PLOT_DATA(tint,9)/TOTAL_PATHS)*100.d0
        !WRITE(*,*),"TINT,TP,TE",TINT,PLOT_DATA(tint,8),PLOT_DATA(tint,8)
        WRITE(IOFILE2,'(3X,3(F14.6))') (tint-1+tint)*SIZE_INT*0.5d0,
     '    PLOT_DATA(tint,8),PLOT_DATA(tint,9)
      ENDDO
C     ENDDO
      
      CALL CLOSEF(IOFILE2,ERROR,*9999)

      CALL EXITS('CALC_TRANSIT_TIME')
      RETURN
 9999 CALL ERRORS('CALC_TRANSIT_TIME',ERROR)
      CALL EXITS('CALC_TRANSIT_TIME')
      RETURN 1
      END

