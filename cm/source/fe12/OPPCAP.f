      SUBROUTINE OPPCAP(nb,NEELEM,NENP,NPNE,NYNE,NYNP,CE,YP,ERROR,*)

C#### Subroutine: OPPCAP
C###  Description:
C###    OPPCAP outputs pulmonary capillary blood flow model results

C***  Created by Kelly Burrowes, May 2003.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter list
      INTEGER NB,NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),
     '  NPNE(NNM,NBFM,NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM)
      REAL*8 CE(NMM,NEM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,ne,noelem,no_perfused,np,NUM_BLOCKED
      REAL*8 area,diameter,diameter_feed,Dm,DP,flow,h,height,length,
     '  length_feed,max_diam,max_length,min_diam,min_length,
     '  per_perfused,TOL,TOTAL_NE,volume
      LOGICAL FEED
      
      CALL ENTERS('OPPCAP',*9999)

      TOL=1.d-8 !zero flow tolerance (if flow.LT.TOL, then no flow)
      length=0.d0
      length_feed=0.d0
      diameter=0.d0
      diameter_feed=0.d0
      flow=0.d0
      volume=0.d0
      height=0.d0
      no_perfused=0
      min_diam=1.d6
      max_diam=0.d0
      min_length=1.d6
      max_length=0.d0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        FEED=.FALSE.
        DO i=1,INLETS(0)
          IF(ne.EQ.INLETS(i)) FEED=.TRUE.
        ENDDO
        DO i=1,OUTLETS(0)
          IF(ne.EQ.OUTLETS(i)) FEED=.TRUE.
        ENDDO
        IF(.NOT.FEED) THEN
          length=length+CE(nm_length,ne) !length of segments
          diameter=diameter+CE(nm_Dh,ne) !hydraulic diameter
          area=PI*CE(nm_a,ne)*CE(nm_b,ne) !mm^2
          volume=volume+area*CE(nm_length,ne) !volume
C!!! NB this is total volume should only include perfused volume
          height=height+(2.0d0*CE(nm_b,ne)) !sheet thickness
          min_diam=MIN(CE(nm_Dh,ne),min_diam)
          max_diam=MAX(CE(nm_Dh,ne),max_diam)
          min_length=MIN(CE(nm_length,ne),min_length)
          max_length=MAX(CE(nm_length,ne),max_length)
          flow=flow+DABS(YP(NYNE(1,1,0,1,ne),1))
        ELSE !FEED
          length_feed=length_feed+CE(nm_length,ne)
          diameter_feed=diameter_feed+CE(nm_Dh,ne)
        ENDIF
C... overall perfusion fraction for entire capillary mesh
        h=2.d0*CE(nm_b,ne) !height of capillary section
        Dm=0.0027d0 !mm, min vessel height RBC can pass through
        IF(CE(nm_Hd,ne).GT.0.d0.AND.h.GT.Dm.AND.
     '    DABS(YP(NYNE(1,1,0,1, ne),1)).GT.TOL)
     '    no_perfused=no_perfused+1 !then segment perfused
      ENDDO !noelem
      length=length/(NEELEM(0)-INLETS(0)-OUTLETS(0)) !avge ne length
      WRITE(OP_STRING,'('' Average length of capillary segments in '//
     '  'whole mesh (um)='',F12.6)') length*1000.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Capillary length range (um)= '',F10.6,'//
     '  ''','',F12.6)') min_length*1000.d0,max_length*1000.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      length_feed=length_feed/(INLETS(0)+OUTLETS(0))
      WRITE(OP_STRING,'('' Average length of feed segments '//
     '  '(um)='',F12.6)') length_feed*1000.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      diameter=diameter/(NEELEM(0)-INLETS(0)-OUTLETS(0))
      WRITE(OP_STRING,'('' Average hydraulic diameter of capillary '//
     '  'segments in whole mesh (um):'',F12.6)') diameter*1000.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Hydraulic diameter range (um)= '',F12.6,'//
     '  ''','',F12.6)') min_diam*1000.d0,max_diam*1000.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      diameter_feed=diameter_feed/(INLETS(0)+OUTLETS(0))
      WRITE(OP_STRING,'('' Average hydraulic diameter of feed '//
     '  'vessels (um):'',F12.6)') diameter_feed*1000.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      height=height/(NEELEM(0)-INLETS(0)-OUTLETS(0))
      WRITE(OP_STRING,'('' Average capillary sheet height '//
     '  '(um):'',F12.6)') height*1000.d0
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      per_perfused=(DBLE(no_perfused)/DBLE(NEELEM(0)))*100.d0 !%
      WRITE(OP_STRING,'('' Percentage of perfused segments '//
     '  'in whole mesh:'',F12.6)') per_perfused
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Total blood flow in capillaries  '//
     '  '(mm**3/s):'',F10.6)') flow
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Average capillary flow per volume '//
     '  '(/s):'',F10.6)') flow/volume
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DP=YP(NYNP(1,1,1,NPNE(1,nb,INLETS(1))),1)-
     '  YP(NYNP(1,1,1,NPNE(2,nb,OUTLETS(1))),1)
      WRITE(OP_STRING,'('' Pressure drop over capillary bed '//
     '  '(Pa) ='',F12.6)') DP
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(INLETS(0).NE.1.OR.OUTLETS(0).NE.1) THEN
        WRITE(OP_STRING,
     '    '(''   Warning: Multiple inlets/outlets.'//
     '    '   Pressure drop only calculated between inlet '',I5,'//
     '    ''' and outlet '',I5)') INLETS(1),OUTLETS(1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(CAP_MESH.EQ.3) THEN
        volume=volume/N_ALVEOLI*300.d6*0.001d0 !mm**3->mL total vol
        WRITE(OP_STRING,'('' Total volume of capillaries  '//
     '    '(mL) ='',F12.6)') volume
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      NUM_BLOCKED=NE_BLOCK(0)
      TOTAL_NE=NEELEM(0)
      NUM_BLOCKED=100*NINT(NUM_BLOCKED/TOTAL_NE)
      WRITE(OP_STRING,'('' % of segments blocked '',F10.4)')
     '  NUM_BLOCKED
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO i=1,INLETS(0)
        ne=INLETS(i)
        IF(NENP(NPNE(1,nb,ne),0).EQ.1) THEN
          np=NPNE(1,nb,ne)
        ELSE
          np=NPNE(2,nb,ne)
        ENDIF
        WRITE(OP_STRING,
     '    '('' Inlet pressure (Pa)'',D10.4,''Inlet flow '//
     '    '(mm**3/s)='',D10.4)') YP(NYNP(1,1,1,np),1),
     '    YP(NYNE(1,1,0,1,ne),1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO
      DO i=1,OUTLETS(0)
        ne=OUTLETS(i)
        IF(NENP(NPNE(1,nb,ne),0).EQ.1) THEN
          np=NPNE(1,nb,ne)
        ELSE
          np=NPNE(2,nb,ne)
        ENDIF
        WRITE(OP_STRING,
     '    '('' Outlet pressure (Pa)'',D10.4,''Outlet flow '//
     '    '(mm**3/s)='',D10.4)') YP(NYNP(1,1,1,np),1),
     '    YP(NYNE(1,1,0,1,ne),1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO
      
      CALL EXITS('OPPCAP')
      RETURN
 9999 CALL ERRORS('OPPCAP',ERROR)
      CALL EXITS('OPPCAP')
      RETURN 1
      END    

      
