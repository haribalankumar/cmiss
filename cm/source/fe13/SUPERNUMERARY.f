      SUBROUTINE SUPERNUMERARY(NBJ,NEELEM,NELIST,NELIST2,NENP,NKJ,
     &  NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,CE,
     &  D_RATIO,DIAM_RATIO,D_MAX_ORDER,L_BRANCH,SE,SV_FREQ,XP,
     &  ZD,VEINS,ERROR,*)
      
C####  Subroutine: SUPERNUMERARY
C###   Description:
C###     SUPERNUMERARY generates the pulmonary arterial (and venous)
C###     supernumerary arteries (SA) which branch off the conventional
C###     arteries (CA) at right angles.
      
C***  Created by Kelly Burrowes, July 2003
C!!!  This subroutine is still under development. !!!

C!!!  NB/ DIAM_RATIO (Strahler order branching ratio) & D_MAX_ORDER
C!!!   (diameter of maximum order vessel) are temporarily passed in
C!!!   and used to allocate diameter values to vessels. 02/03/04
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
C      INCLUDE 'pulm00.cmn'
C      INCLUDE 'grou00.cmn'
!     Parameter values
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),
     &  NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),NP_INTERFACE(0:NPM,0:3),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),D_RATIO,L_BRANCH(*),
     '  SE(NSM,NBFM,NEM),SV_FREQ,XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
      LOGICAL VEINS
!     Local variables
      INTEGER CA_SA,i,ISEED,j,MAX_ORD,MIN_ORD,n_horsfield,nb,
     &  nd,ne,ne0,ne2,NELIST_TEMP(0:NE_R_M),ne_bifur,ne_nxi,ne2_old(10),
     &  ngen,nj,noelem,noelem2,noelem_nxi,no_it,nonode,nonode2,np,np1,
     &  np2,np_bifur,NPLIST(0:NP_R_M),NPLIST2(0:NP_R_M),NUM_CA,NUM_SA,
     &  nv1,nv2,ord,point,point1,remove,STRAHLER,STRAHLER_ADD,temp_nv
      REAL*8 A(NJT),ANG2V,angle,B(NJT),CA_SA_AVGE,
     &  CM_RANDOM_NUMBER,D1,D2,D_ORD,diameter,diameter1,D_MAX_ORDER,
     &  DIAM_RATIO,dist,FREQUENCY,L_D_RATIO,length_A,length_B,length_ne,
     &  length_ne2,max_angle,max_angle1,min_angle,min_angle1,min_dist,
     &  min_dist1,parent_diam,radius,RATIO,V(NJT)
      LOGICAL CONTINU
      CHARACTER STRING*255
      
      CALL ENTERS('SUPERNUMERARY',*9999)

      nb=NBJ(1,NEELEM(1,nr))
      ne=NET(0)
      noelem=NET(nr)
      np=NPT(0)
      nonode=NPT(nr)
      NPLIST(0)=0
      NUM_CA=NEELEM(0,nr) !number of conventional arteries
      NELIST(0)=0 !initialise element group - supernumerary vessels
      NELIST_TEMP(0)=0 !initialise element group - stem vessels
      NELIST2(0)=0
      WRITE(OP_STRING,'(''1ST SUPERNUMERARY ELEMENT'',I6)') ne+1
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''1ST SUPERNUMERARY NODE'',I6)') np+1
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)     
      
      MIN_ORD=1 !only modelling 11 Strahler orders (conducting vessels)
      MAX_ORD=11
      
      CALL ASSERT(NDT.GT.0,
     '  '>>No data points defined to generate supernumerary vessels'
     '  ,ERROR,*9999)
      ISEED=1
      CA_SA_AVGE=0.d0
C     DO noelem2=1,NEELEM(0,nr) !for each arterial element
      DO noelem2=NEELEM(0,nr),1,-1 !looping from terminals upwards to enable corrent Strahler ordering
        ne2=NEELEM(noelem2,nr)
        NELIST_TEMP(0)=NELIST_TEMP(0)+1
        NELIST_TEMP(NELIST_TEMP(0))=ne2 !stores stem vessels
        ngen=NORD(1,ne2) !generation #
        ord=NORD(3,ne2)
        IF(ngen.EQ.0) THEN
          CALL ASSERT(.FALSE.,'>>Branch orders must be evaluated first',
     &      ERROR,*9999)
        ELSEIF(ngen.GT.7) THEN !supernumeraries start from 8th generation
          np1=NPNE(1,nb,ne2)
          np2=NPNE(2,nb,ne2)
          length_ne=CE(1,ne2)
C         parent_diam=CE(4,ne2)*2.d0 !radius*2
          nv1=NVJE(1,nb,nj_radius,ne2) !version at nn=1 of element (upstream end)
          nv2=NVJE(2,nb,nj_radius,ne2)
          parent_diam=(XP(1,nv1,nj_radius,np1)+XP(1,nv2,nj_radius,np2))!average diam
          IF(ngen.GE.7.AND.ord.GE.MIN_ORD) THEN
            FREQUENCY=SV_FREQ !average supernumerary branch frequency
            CA_SA=0
            DO WHILE(FREQUENCY.GT.0.d0)
              IF(CM_RANDOM_NUMBER(ISEED).LT.FREQUENCY) CA_SA=CA_SA+1
              FREQUENCY=FREQUENCY-1.d0                              
            ENDDO !should give average frequency=SV_FREQ
          ELSE
            CA_SA=0
          ENDIF
          IF(CA_SA.GT.0) CA_SA_AVGE=CA_SA_AVGE+DBLE(CA_SA) !average frequency
          DO i=1,10
            ne2_old(i)=0 !initialise
          ENDDO
          DO i=1,CA_SA !for each supernumerary artery
            nonode=nonode+1
            np=np+1
            CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
            CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
            NPNODE(nonode,nr)=np
            DO nj=1,NJT
              V(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
              XP(1,1,nj,np)=XP(1,1,nj,np1)+DBLE(i)/DBLE((CA_SA+1))*V(nj)
            ENDDO
            CONTINU=.TRUE. !initialise
C...  Parent element is now split into 2 to accomodate supernumerary element
            ne=ne+1 !& create a new element for the second half
            noelem=noelem+1
            CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',
     '        ERROR,*9999)
            CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
            NELIST_TEMP(0)=NELIST_TEMP(0)+1
            NELIST_TEMP(NELIST_TEMP(0))=ne !stores stem vessels
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=np
            NPNE(2,nb,ne)=NPNE(2,nb,ne2)
            NPNE(2,nb,ne2)=np !changes end node of existing element
            NORD(1,ne)=NORD(1,ne2)+1 !generation number
            IF(NORD(1,ne).GT.GENM) THEN !maximum generation #
              WRITE(OP_STRING,'(''Generations exceeded maximim'',I6)')
     &          NORD(1,ne)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            NORD(2,ne2)=1 !reset for call to order_number, NORD(2,ne2) recalculated below
            CALL ORDER_NUMBER(n_horsfield,ne2,NORD,NXI,STRAHLER,
     &        STRAHLER_ADD,ERROR,*9999)
            !call ORDER_NUMBER with ne2 to use correct NXI values
            NORD(2,ne)=n_horsfield !Horsfield order
            NORD(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
            NXI(1,0,ne)=NXI(1,0,ne2)
            DO noelem_nxi=1,NXI(1,0,ne)
              NXI(1,noelem_nxi,ne)=NXI(1,noelem_nxi,ne2)
            ENDDO
 !fix NXI array for old element (ne2)
            NXI(1,0,ne2)=2
            NXI(1,1,ne2)=ne
            NXI(1,2,ne2)=ne+1 !new supernumerary element
C... ordering of element ne2 done below, after supernumerary elem is created
            temp_nv=NVJP(nj_radius,NPNE(2,nb,ne))
            CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(2,nb,ne),
     '        NPNE(1,nb,ne),nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
            !NB/ this resets all version # to be=1
            NVJP(nj_radius,NPNE(2,nb,ne))=temp_nv
           !ensures same number of versions as before call to GN1DNEJ
            NVJP(nj_radius,NPNE(1,nb,ne))=2 !2 versions
            length_ne=0.d0
            length_ne2=0.d0
            DO nj=1,NJT
              length_ne=length_ne+(XP(1,1,nj,NPNE(2,nb,ne))-
     '          XP(1,1,nj,NPNE(1,nb,ne)))**2.d0
              length_ne2=length_ne2+(XP(1,1,nj,NPNE(2,nb,ne2))-
     '          XP(1,1,nj,NPNE(1,nb,ne2)))**2.d0
            ENDDO
            CE(1,ne)=DSQRT(length_ne)
            CE(1,ne2)=DSQRT(length_ne2)
C           CE(4,ne)=CE(4,ne2) !radius
            nv1=NVJE(1,nb,nj_radius,ne2) !version at nn=1 of element
            NVJE(2,nb,nj_radius,ne)=NVJE(2,nb,nj_radius,ne2) !puts original end in new end
            NVJE(2,nb,nj_radius,ne2)=1 !use first version of new node
            NVJE(1,nb,nj_radius,ne)=2 !use second version of new node
            NVJP(nj_radius,np)=2 !2 versions at this node
            nv2=NVJE(2,nb,nj_radius,ne)
            radius=(XP(1,nv1,nj_radius,NPNE(1,nb,ne2))
     &        +XP(1,nv2,nj_radius,NPNE(2,nb,ne)))/2.d0
            XP(1,1,nj_radius,NPNE(1,nb,ne))=radius
            XP(1,2,nj_radius,NPNE(1,nb,ne))=radius
C           XP(1,1,nj_radius,NPNE(2,nb,ne))=radius !this node radius already defined NPNE(2,nb,ne)=NPNE(2,nb,ne2)
C...  Find closest terminal bronchiole node to the supernumerary node
            length_ne=CE(1,ne2)+CE(1,ne)
            min_angle=PI/2.d0-0.25d0
            max_angle=PI/2.d0+0.25d0
            min_dist=1.d6
            DO nj=1,NJT
              B(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
            ENDDO
            DO nd=1,NDT              
              dist=0.d0
              DO nj=1,NJT
                dist=dist+(ZD(nj,nd)-XP(1,1,nj,np))**2.d0
              ENDDO              
              dist=DSQRT(dist)
              IF(dist.LT.min_dist) THEN !check if angle to np is approx PI/2
                DO nj=1,NJT
                  A(nj)=ZD(nj,nd)-XP(1,1,nj,np) !using data points
C                  B(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
                ENDDO
                angle=ANG2V(A,B)
                IF(angle.GT.min_angle.AND.angle.LT.max_angle) THEN
                  min_dist=dist
                  point=nd !closest data point
                ENDIF
              ENDIF
            ENDDO !nd
C... Create node & element towards closest terminal airway node
            nonode=nonode+1
            np=np+1
            NENP(np,0,nr)=0 !initialise for new node
            CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
            CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
            NPNODE(nonode,nr)=np
            !testing different length/diameter ratios
            L_D_RATIO=5.0d0
            IF(VEINS) L_D_RATIO=6.d0
            diameter=D_RATIO*parent_diam!*length_ne/L_D_ratio !ratio of CA:SA diam
            length_A=0.d0
            DO nj=1,NJT
C              A(nj)=XP(1,1,nj,point)-XP(1,1,nj,np-1)
              A(nj)=ZD(nj,point)-XP(1,1,nj,np-1)
              length_A=length_A+A(nj)**2.d0
            ENDDO
            length_A=DSQRT(length_A) !length of vector A

C... To determine which Strahler order this vessel fits into
            ord=1
            !Diameter ratio
            D_ORD=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-ord))
            D1=(D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-(ord+1)))+D_ORD)/2.d0
            IF(diameter.LT.D1.AND.diameter.GT.0.d0) THEN
              CONTINU=.FALSE.
            ELSE
              CONTINU=.TRUE.
            ENDIF
            DO WHILE(ord.LT.MAX_ORD.AND.CONTINU)
              ord=ord+1
              D_ORD=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-ord))
              D1=(D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-(ord+1)))+D_ORD)
     &          /2.d0
              D2=(D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-(ord-1)))+D_ORD)
     &          /2.d0
              IF(diameter.LT.D1.AND.diameter.GT.D2) CONTINU=.FALSE.
            ENDDO !WHILE            
C... Create new element            
            noelem=noelem+1
            ne=ne+1
            CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',ERROR,
     &        *9999)
            CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
            NELIST(0)=NELIST(0)+1
            NELIST(NELIST(0))=ne !stores supernumerary elements
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=np-1
            NPNE(2,nb,ne)=np
            NORD(1,ne)=NORD(1,ne2)+1 !generation
            NORD(2,ne)=ord !Horsfield order # (same as Strahler)
            NORD(3,ne)=ord !Strahler order
            IF(DIAM_STRAHLER) THEN
              NORD(4,ne)=ord !currently just have diam-def Strahler=Strahler
            ENDIF
C... Now defining order for parent element ne2
            IF(NORD(3,ne).EQ.NORD(3,ne-1)) THEN
              NORD(3,ne2)=NORD(3,ne)+1
            ELSE
              NORD(3,ne2)=MAX(NORD(3,ne),NORD(3,ne-1))
            ENDIF
            IF(NORD(2,ne).EQ.NORD(2,ne-1)) THEN
              NORD(2,ne2)=NORD(2,ne)+1
            ELSE
              NORD(2,ne2)=MAX(NORD(2,ne),NORD(2,ne-1))+1
            ENDIF
            IF(i.GT.1) THEN !if # supernumeraries greater than 1 for this branch
              DO j=1,i-1
                NORD(2,ne2_old(j))=NORD(2,NXI(1,1,ne2_old(j)))+1
                IF(ord.GE.NORD(3,NXI(1,1,ne2_old(j)))) THEN
                  NORD(3,ne2_old(j))=
     &              MAX(NORD(3,NXI(1,1,ne2_old(j))),ord)+1
C                  WRITE(*,*),"STRAHLER INCREASED",ne2_old(j)
                ENDIF
              ENDDO
            ENDIF
            IF(DIAM_STRAHLER) THEN
              NORD(4,ne)=NORD(3,ne) !diameter-defined Strahler order
            ENDIF
            
            CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(2,nb,ne),NPNE(1,nb,ne),
     &        nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
            NVJP(nj_radius,NPNE(1,nb,ne))=3            
C           CE(4,ne)=diameter/2.d0 !radius
C           nv=NVJE(1,nb,nj_radius,ne) !version at nn=1 of element (upstream end)
C           NVJE(3,nb,nj_radius,ne)=3 !3rd version
            XP(1,3,nj_radius,NPNE(1,nb,ne))=diameter/2.d0 !already defined
            XP(1,1,nj_radius,NPNE(2,nb,ne))=diameter/2.d0
C           XP(1,nv,nj_radius,ne)=diameter/2.d0 !radius
C... arterial node only goes to terminal airway node (i.e. not past it)
            IF(L_D_RATIO*diameter.GT.length_A) THEN
              CE(1,ne)=length_A !length
C             diameter=length_A/L_D_RATIO !stop when reach terminal airway node
            ELSE IF(NORD(3,ne).LE.MIN_ORD) THEN
              CE(1,ne)=L_BRANCH(MIN_ORD) !OR = L_D_RATIO*diameter
              !won't bifurcate - already order 1 vessel
            ELSE
              CE(1,ne)=L_D_RATIO*diameter
              NPLIST(0)=NPLIST(0)+1 !this node will bifurcate
              NPLIST(NPLIST(0))=np
              NENP(np,0,nr)=NENP(np,0,nr)+1
              NENP(np,NENP(np,0,nr),nr)=ne !new element
            ENDIF
            DO nj=1,NJT
              XP(1,1,nj,np)=XP(1,1,nj,np-1)+(CE(1,ne)/length_A)*A(nj)
            ENDDO
            ne2_old(i)=ne2
            ne2=ne-1 !for next loop through CA_SA
          ENDDO !i
        ELSE !ngen.LE.7
          !apply correct order number to these elements
          NORD(2,ne2)=1
          CALL ORDER_NUMBER(n_horsfield,ne2,NORD,NXI,STRAHLER,
     &      STRAHLER_ADD,ERROR,*9999)
          NORD(2,ne2)=n_horsfield !Horsfield order
          NORD(3,ne2)=STRAHLER+STRAHLER_ADD !Strahler order
        ENDIF !ngen
      ENDDO
      min_angle=PI/4.d0-0.4d0 !angle limits for branching
      max_angle=PI/2.d0+0.4d0
      min_angle1=PI/4.d0-0.4d0
      max_angle1=PI/2.d0+0.4d0
      
C... For each node in NPLIST create supernumerary arterial bifurcation
      no_it=0
      DO WHILE(NPLIST(0).NE.0)
        no_it=no_it+1
        WRITE(OP_STRING,'(''NUMBER OF ITERATIONS'',2(I6))') no_it,
     &    NPLIST(0)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        NPLIST2(0)=0
        DO nonode2=1,NPLIST(0)
          np_bifur=NPLIST(nonode2)
          ne_bifur=NENP(np_bifur,NENP(np_bifur,0,nr),nr) !element number
          min_dist=1.d6
          min_dist1=1.d6
          point=0
          point1=0
C         diameter=CE(4,ne_bifur)*2.d0
C          nv=NVJE(1,nb,nj_radius,ne_bifur) !version at nn=1 of element
          diameter=XP(1,1,nj_radius,np_bifur)*2.d0
C... determine what order this size branch fits into
          IF(no_it.EQ.1) THEN !1st time through only
            ord=1
C           D1=D_BRANCH(ord+1) !Horsfield data
C           D2=D_BRANCH(ord-1)
            !Diameter ratio
            D_ORD=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-ord))
            D1=(D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-(ord+1)))+D_ORD)/2.d0
            IF(diameter.LT.D1.AND.diameter.GT.0.d0) THEN
              CONTINU=.FALSE.
            ELSE
              CONTINU=.TRUE.
            ENDIF
            DO WHILE(ord.LT.MAX_ORD.AND.CONTINU)
              ord=ord+1
              D_ORD=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-ord))
              D1=(D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-(ord+1)))+D_ORD)
     &          /2.d0
              D2=(D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-(ord-1)))+D_ORD)
     &          /2.d0
              IF(diameter.LT.D1.AND.diameter.GT.D2) CONTINU=.FALSE.
            ENDDO !WHILE
          ELSE
            ord=NORD(3,ne_bifur)
          ENDIF
          IF(ord.LE.MIN_ORD) THEN
            ord=MIN_ORD
          ELSE
            ord=ord-1 !new daughter branches order
          ENDIF
          DO nj=1,NJT
            B(nj)=XP(1,1,nj,NPNE(2,nb,ne_bifur))-
     '        XP(1,1,nj,NPNE(1,nb,ne_bifur))
          ENDDO
          DO nd=1,NDT
            dist=0.d0
            DO nj=1,NJT
              dist=dist+(ZD(nj,nd)-XP(1,1,nj,np_bifur))**2.d0
            ENDDO            
            dist=DSQRT(dist)
            IF(dist.LT.min_dist.OR.dist.LT.min_dist1) THEN
              DO nj=1,NJT
                A(nj)=ZD(nj,nd)-XP(1,1,nj,np_bifur)
              ENDDO
              angle=ANG2V(A,B)
              IF(dist.LT.min_dist) THEN
                IF(angle.GE.min_angle.AND.angle.LE.max_angle) THEN
                  min_dist=dist
                  point=nd
                ENDIF
              ELSE
                IF(angle.GE.min_angle1.AND.angle.LE.max_angle1) THEN
                  min_dist1=dist
                  point1=nd
                ENDIF
              ENDIF
            ENDIF
          ENDDO !nd
          IF(point.EQ.0.OR.point1.EQ.0) THEN
            IF(point.EQ.0) THEN
              WRITE(OP_STRING,
     '          '(''Error: point1 not found in SUPERNUMERARY'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ELSEIF(point1.EQ.0) THEN
              WRITE(OP_STRING,
     '          '(''Error: point2 not found in SUPERNUMERARY'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            length_A=0.d0
            length_B=0.d0
            DO nj=1,NJT
              A(nj)=ZD(nj,point)-XP(1,1,nj,np_bifur) !branch 1
              length_A=length_A+A(nj)**2.d0
              B(nj)=ZD(nj,point1)-XP(1,1,nj,np_bifur) !branch 2
              length_B=length_B+B(nj)**2.d0
            ENDDO
            length_A=DSQRT(length_A)
            length_B=DSQRT(length_B)
            nonode=nonode+1
            np=np+1
            NENP(np,0,nr)=0 !initialise for new node
            CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
            CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
            NPNODE(nonode,nr)=np
            diameter=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-ord)) !D_BRANCH(ord)
            L_D_ratio=6.d0
            IF(VEINS) L_D_RATIO=8.d0
            noelem=noelem+1
            ne=ne+1
            CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',
     '        ERROR,*9999)
            CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
            NELIST(0)=NELIST(0)+1
            NELIST(NELIST(0))=ne !stores supernumerary elements
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=np_bifur
            NPNE(2,nb,ne)=np
            NORD(1,ne)=NORD(1,ne_bifur)+1
            IF(NORD(1,ne).GT.GENM) THEN !maximum generation #
              WRITE(OP_STRING,'(''Generations exceeded maximum'',I6)')
     &          NORD(1,ne)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            NORD(2,ne)=ord !Horsfield order
            NORD(3,ne)=ord !Strahler order
            IF(DIAM_STRAHLER) THEN
              NORD(4,ne)=ord !diameter-defined Strahler (temporarily=ord)
            ENDIF
            CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(2,nb,ne),NPNE(1,nb,ne),
     '        nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
            NVJP(nj_radius,np_bifur)=3
C           CE(4,ne)=diameter/2.d0
            NVJE(2,nb,nj_radius,ne_bifur)=2
            NVJE(1,nb,nj_radius,ne)=1 !version at nn=1 of element
            XP(1,2,nj_radius,np_bifur)=diameter/2.d0 !already defined
            XP(1,1,nj_radius,NPNE(2,nb,ne))=diameter/2.d0
            D1=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-MIN_ORD))
            !D1=D_BRANCH(MIN_ORD)
            IF(diameter*L_D_ratio.GE.length_A) THEN
              CE(1,ne)=length_A !branch only goes to acinar unit
            ELSEIF(diameter.LE.D1.OR.ord.LE.MIN_ORD) THEN
              !dont put into NPLIST
              CE(1,ne)=L_BRANCH(ord) !OR=L_D_ratio*diameter
            ELSE
              CE(1,ne)=L_D_ratio*diameter
              NPLIST2(0)=NPLIST2(0)+1 !supernumerary branches again
              NPLIST2(NPLIST2(0))=np
              NENP(np,0,nr)=NENP(np,0,nr)+1
              NENP(np,NENP(np,0,nr),nr)=ne !new element
            ENDIF
            DO nj=1,NJT
              XP(1,1,nj,np)=XP(1,1,nj,np_bifur)+CE(1,ne)*A(nj)
     '          /length_A                            
            ENDDO

            nonode=nonode+1
            np=np+1
            NENP(np,0,nr)=0 !initialise
            CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
            CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
            NPNODE(nonode,nr)=np
            diameter1=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-ord)) !D_BRANCH(ord)
            noelem=noelem+1
            ne=ne+1
            CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',
     '        ERROR,*9999)
            CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
            NELIST(0)=NELIST(0)+1
            NELIST(NELIST(0))=ne !stores supernumerary elements
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=np_bifur
            NPNE(2,nb,ne)=np
            NORD(1,ne)=NORD(1,ne_bifur)+1
            IF(NORD(1,ne).GT.GENM) THEN !maximum generation #
              WRITE(OP_STRING,'(''Generations exceeded maximim'',I6)')
     &          NORD(1,ne)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            NORD(2,ne)=ord !Horsfield order
            NORD(3,ne)=ord !Strahler order
            IF(DIAM_STRAHLER) THEN
              NORD(4,ne)=ord !diameter-defined Strahler (temporarily=ord)
            ENDIF
            CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(2,nb,ne),NPNE(1,nb,ne),
     &        nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
            NVJP(nj_radius,np_bifur)=3 !3 versions
C           CE(4,ne)=diameter1/2.d0
C           nv=NVJE(1,nb,nj_radius,ne) !version at nn=1 of element
            NVJE(1,nb,nj_radius,ne)=1
C            NVJE(3,nb,nj_radius,ne_bifur)=3
            XP(1,3,nj_radius,np_bifur)=diameter1/2.d0 !already defined
            XP(1,1,nj_radius,NPNE(2,nb,ne))=diameter1/2.d0
            D1=D_MAX_ORDER/(DIAM_RATIO**(MAX_ORD-MIN_ORD))
            !D1=D_BRANCH(MIN_ORD)
            IF(diameter1*L_D_ratio.GE.length_B) THEN
              CE(1,ne)=length_B
            ELSEIF(diameter1.LE.D1.OR.ord.LE.MIN_ORD)
     '          THEN
              CE(1,ne)=L_D_ratio*diameter1
              !dont put into NPLIST
            ELSE
              CE(1,ne)=L_D_ratio*diameter1
              NPLIST2(0)=NPLIST2(0)+1 !supernumerary branches again
              NPLIST2(NPLIST2(0))=np
              NENP(np,0,nr)=NENP(np,0,nr)+1
              NENP(np,NENP(np,0,nr),nr)=ne !new element
            ENDIF
            DO nj=1,NJT
              XP(1,1,nj,np)=XP(1,1,nj,np_bifur)+CE(1,ne)*B(nj)/length_B
            ENDDO
          ENDIF !point.EQ.0
        ENDDO !nonode2
        NPLIST(0)=0
        DO nonode2=1,NPLIST2(0) !copy NPLIST2 into NPLIST
          np2=NPLIST2(nonode2)
          NPLIST(0)=NPLIST(0)+1
          NPLIST(NPLIST(0))=np2
        ENDDO
      ENDDO !WHILE
      NEELEM(0,nr)=noelem
      NPNODE(0,nr)=nonode
      NPNODE(0,0)=np
      NEELEM(0,0)=ne
      NPT(nr)=nonode
      NET(nr)=noelem
      NPT(0)=np
      NET(0)=ne
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*9999)
      STRING='supernumerary'
      CALL GRELEM_SUB(NELIST,STRING,.TRUE.,ERROR,*9999)
      NUM_SA=NELIST(0) !Number of supernumerary arteries
      WRITE(OP_STRING,'(''LAST SUPERNUMERARY ELEMENT'',I6)') ne
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      !elements in reverse order, put them into correct order
      DO noelem=1,NELIST_TEMP(0)
        NELIST2(0)=NELIST2(0)+1
        NELIST2(NELIST2(0))=NELIST_TEMP(NELIST_TEMP(0)-noelem+1)
      ENDDO
      
      !TEMPORARY
      IF(VEINS) THEN !remove 1st 4 elems from analysis (4 main veins feed into left atrium)
        remove=4 !# to remove
        NELIST2(0)=NELIST2(0)-remove
        DO noelem=1,NELIST2(0)
          ne=NELIST2(noelem)
          NELIST2(noelem)=NELIST2(noelem+remove)
        ENDDO
      ENDIF
        
      STRING='stem' !element group of stem vessels - used in OPMESH2
      CALL GRELEM_SUB(NELIST2,STRING,.TRUE.,ERROR,*9999)
      
      RATIO=DBLE(NUM_SA)/DBLE(NUM_CA)
      CA_SA_AVGE=CA_SA_AVGE/DBLE(NUM_CA) !average # SA emerging from CA
      WRITE(OP_STRING,'(''# conventional='',I6,''# supernumerary='',
     &  I6,''RATIO='',F8.3)') NUM_CA,NUM_SA,RATIO
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''Average # SVs emerging from CVs'',F8.3)')
     &  CA_SA_AVGE
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C... Recalculate generation number for each branch
      noelem=1
      ne=NEELEM(noelem,nr)
      ne0=NXI(-1,1,ne) !parent
      IF(ne0.ne.0) THEN
        WRITE(OP_STRING,
     &    '(''NEELEM(1,nr) not 1st element, generation '//
     &    ' numbers will be incorrect !! '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      NELIST2(0)=0
      NELIST(0)=1
      NELIST(NELIST(0))=ne
      DO WHILE(noelem.LT.NEELEM(0,nr).AND.NELIST(0).NE.0)
        noelem2=1
        DO WHILE(noelem2.LE.NELIST(0))
          ne=NELIST(noelem2)
          noelem=noelem+1 !counts number of elements done
          ne0=NXI(-1,1,ne) !parent
          IF(ne0.NE.0)THEN
            ngen=NORD(1,ne0) !parent generation
            IF(NXI(1,0,ne0).EQ.1)THEN !single daughter
              IF(NORD(5,ne).EQ.0)THEN !same generation
                NORD(1,ne)=ngen
              ELSE IF(NORD(5,ne).EQ.1)THEN !start of 'half' branch
                NORD(1,ne)=ngen+1
              ENDIF
            ELSE IF(NXI(1,0,ne0).GE.2)THEN
              NORD(1,ne)=ngen+1
            ENDIF
          ELSE
            NORD(1,ne)=1 !generation 1
          ENDIF
          DO noelem_nxi=1,NXI(1,0,ne) !put all forward elements into NELIST2
            ne_nxi=NXI(1,noelem_nxi,ne)
            NELIST2(0)=NELIST2(0)+1
            NELIST2(NELIST2(0))=ne_nxi
          ENDDO
          noelem2=noelem2+1 !increment loop
        ENDDO !noelem2
        NELIST(0)=0
        DO noelem2=1,NELIST2(0)
          NELIST(0)=NELIST(0)+1
          NELIST(NELIST(0))=NELIST2(noelem2)
        ENDDO
        NELIST2(0)=0
      ENDDO !noelem
      IF(NELIST(0).EQ.0.AND.noelem.LT.NEELEM(0,nr)) THEN
        WRITE(OP_STRING,
     &    '(''Generation loop exited before all elements numbered! '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      
      
      CALL EXITS('SUPERNUMERARY')
      RETURN
 9999 CALL ERRORS('SUPERNUMERARY',ERROR)
      CALL EXITS('SUPERNUMERARY')
      RETURN 1
      END

      
