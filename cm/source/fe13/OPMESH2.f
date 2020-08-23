      SUBROUTINE OPMESH2(NBJ,NBRANCHES,NEELEM,NELIST,NENP,NORD,NPLIST,
     &  NPNE,nr,NVJE,nx,NXI,NYNP,BRANCHES,CP,diameters,lengths,STATS,XP,
     &  YP,LISOLN,OP_COORD,ERROR,*)
      
C#### Subroutine: OPMESH2
C###  Description:updated
C###    OPMESH2 outputs lung mesh data. #Modified

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBRANCHES(5,NE_R_M),NEELEM(0:NE_R_M),
     &  NELIST(0:NEM),NENP(NPM,0:NEPM),NORD(5,NE_R_M),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),nr,NVJE(NNM,NBFM,NJM,NEM),nx,
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 BRANCHES(10,NE_R_M),CP(NMM,NPM),diameters(NEM),
     &  lengths(NEM),STATS(21,NE_R_M),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL LISOLN,OP_COORD
!     Local Variables
      INTEGER i,INDEX(4),j,N,nb,NBINS(5),N_BR,ne,
     &  ne0,ne1,ne2,ne_major,ne_minor,ne_next,ne_offset,ngen,nh,nj,
     &  NMAX_GEN(4),noelem,nonode,np,np0,np1,np2,N_SEGMENTS,
     &  NTALLY(4,6,GENM),N_TERMINAL(GENM),ntotal,NTOTALN(20),num_ddp,
     &  num_llp,NUM_SCHEMES,nv1,nv2,ny,ny1,ny2,OP_FORMAT,SUM_TERM
      REAL*8 angle,average_term_gen,BINS(5),
     &  length,MEAN_DIAM,MEANS(20),
     &  radius,RATIOS(4,3),R_SQ(4,3),SD(4,6,GENM),SDT(20),slope,
     &  SUM_MEAN(4,6,GENM),undefined,X(GENM),XP1(3),XP2(3),XP3(3),
     &  YREGRESS(GENM,3),sum_length,MAX_CONC
      REAL*8 CONC_EVAL,RADIUS_MEAN
      LOGICAL ADD
      CHARACTER STRING*255


C     INDEX(i) stores the order type (i=1,2,3,4 for generation, Horsfield
C     order, Strahler order, diameter-defined Strahler order
C     respectively) for the current element.

C     RATIOS(i,j) stores the branching (j=1), length (j=2), and diameter
C     (j=3) ratios for generations (i=1), Horsfield orders (i=2),
C     Strahler orders (i=3), and diameter-defined Strahler orders (i=4).

C     SD(i,j,n) stores the standard deviations for order type i,
C     parameter j, in order n.

C     MEANS(j) stores the mean value of parameter j for all branches.      
C     SDT(j) stores the SD for parameter j for all branches.
      
C     SUM_MEAN(i,j,n) stores the mean values for order type i, parameter
C     j, order n.
      
C     BRANCHES(j,n) stores the individual values of parameter j for
C     branch n.
      
      CALL ENTERS('OPMESH2',*9999)
      
      CALL ASSERT(nj_radius.GT.0,'>>Define radius field first',ERROR,
     &  *9999)

      nb=NBJ(1,1) !temp

      undefined=1.d4
      OP_FORMAT=2
      NUM_SCHEMES=3 !# of ordering schemes (generations, H order, S orders)
      IF(DIAM_STRAHLER) NUM_SCHEMES=4 !(+diameter-defined Strahler ordering)
      nh=NH_LOC(1,nx)

      IF(LISOLN)THEN
        DO nonode=1,NPLIST(0)
          np=NPLIST(nonode)
          ne=NENP(np,1) !first element
          IF(np.EQ.NPNE(2,nb,ne))THEN
            length=0.d0
            np1=NPNE(1,nb,ne)
            np2=NPNE(2,nb,ne)
            DO nj=1,NJT
              length=length+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2
            ENDDO !nj
            sum_length=DSQRT(length)
          ELSE
            sum_length=0.d0
          ENDIF
          ne0=NXI(-1,1,ne)
          DO WHILE(ne0.NE.0)
            length=0.d0
            np1=NPNE(1,nb,ne0)
            np2=NPNE(2,nb,ne0)
            DO nj=1,NJT
              length=length+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2
            ENDDO !nj
            sum_length=sum_length+DSQRT(length)
            ne0=NXI(-1,1,ne0)
          ENDDO
          
          ny1=NYNP(1,1,nh,np) !for temperature in K, or concentration
          IF(ITYP3(nr,nx).EQ.1)THEN
            WRITE(OP_STRING,'(2(F10.4))') sum_length,YP(ny1,1)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ITYP3(nr,nx).EQ.2)THEN
            ny2=NYNP(1,1,nh+1,np) !for humidity in g/mm^3
            MAX_CONC=CONC_EVAL(YP(ny1,1))
            WRITE(OP_STRING,'(6(F10.4))') sum_length,YP(ny1,1)-273.15d0,
     &        YP(ny2,1),YP(ny2,1)/MAX_CONC*100.d0,CP(4,np)-273.15d0,
     &        CP(8,np)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nonode

      ELSE IF(.NOT.OP_COORD.AND.(.NOT.LISOLN))THEN
        CALL ASSERT(NKM.GE.2,'>>Increase NKM to 2',ERROR,*9999)
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(1,ne)
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          diameters(ne)=2.d0*RADIUS_MEAN(NBJ,ne,NPNE,NVJE,XP)
          length=0.d0
          DO nj=1,NJT
            length=length+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2
          ENDDO !nj
          lengths(ne)=DSQRT(length)
        ENDDO !noelem
        
C...... Initialise arrays
        DO j=1,20
          MEANS(j)=0.d0
          SDT(j)=0.d0
          NTOTALN(j)=0
        ENDDO
        DO N=1,GENM
          N_TERMINAL(N)=0
          DO i=1,NUM_SCHEMES
            DO j=1,6
              SUM_MEAN(i,j,N)=0.d0
              SD(i,j,N)=0.d0
              NTALLY(i,j,N)=0
            ENDDO
          ENDDO
        ENDDO
        DO N=1,NE_R_M
          BRANCHES(1,n)=0.d0
        ENDDO
        DO N=1,5
          NBINS(N)=0
          BINS(N)=0.d0
        ENDDO !N
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          WRITE(STRING,'(''>>Increase NE_R_M'')')
          CALL ASSERT(ne.LE.NE_R_M,STRING,ERROR,*9999)
          STATS(2,ne)=undefined
          STATS(5,ne)=undefined
          STATS(6,ne)=undefined
          DO N=1,NXI(1,0,ne)
            ne2=NXI(1,N,ne)
            STATS(2,ne2)=undefined
            STATS(5,ne2)=undefined
            STATS(6,ne2)=undefined
          ENDDO !N
        ENDDO !noelem
        
        ntotal=0
        num_ddp=0
        num_llp=0
        N=0
        
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          ne0=NXI(-1,1,ne)
          DO i=1,NUM_SCHEMES
            INDEX(i)=NORD(i,ne)
            CALL ASSERT(INDEX(i).LE.GENM,'>>Increase GENM',ERROR,
     &        *9999)
          ENDDO
          ADD=.FALSE.
          IF(ne0.EQ.0)THEN
            ADD=.TRUE.
          ELSE IF(ne0.NE.0.AND.NORD(1,ne0).NE.INDEX(1))THEN
 !only count as extra branch if at start
            ADD=.TRUE.
          ENDIF
          IF(ADD)THEN
            N=N+1
            DO i=1,NUM_SCHEMES
              NBRANCHES(i,N)=INDEX(i) !generation, H order, S order, D-D S order
            ENDDO !i
            IF(ne0.NE.0)THEN
              NBRANCHES(5,N)=NORD(3,ne0) !Strahler order of parent
            ELSE
              NBRANCHES(5,N)=0
            ENDIF
C.......... Add length of all segments along branch, calculate mean
C.......... diameter
            
            N_SEGMENTS=1
            
            MEAN_DIAM=diameters(ne)
            BRANCHES(1,N)=lengths(ne)
            
            ne_next=ne
            DO WHILE(NXI(1,0,ne_next).EQ.1)
              ne_next=NXI(1,1,ne_next) !next segment
              BRANCHES(1,N)=BRANCHES(1,N)+lengths(ne_next) !sum lengths
              MEAN_DIAM=MEAN_DIAM+diameters(ne_next) !sum diameters
              N_SEGMENTS=N_SEGMENTS+1 !count number of segments in branch
            ENDDO
            STATS(6,ne)=MEAN_DIAM/DBLE(N_SEGMENTS)
            BRANCHES(2,N)=MEAN_DIAM/DBLE(N_SEGMENTS) !record mean diameter
            BRANCHES(5,N)=BRANCHES(1,N)/BRANCHES(2,N) !L/D
c            IF(BRANCHES(5,N).GT.8.d0)THEN
c              BRANCHES(5,N)=undefined
c            ENDIF
C.......... Calculate branching angle to parent
            IF(INDEX(1).GT.1)THEN
              nb=NBJ(1,ne)
              np0=NPNE(1,nb,ne0) !start of parent
              np1=NPNE(1,nb,ne) !start node
              np2=NPNE(2,nb,ne) !end node
              DO nj=1,3
                XP1(nj)=XP(1,1,nj,np0)
                XP2(nj)=XP(1,1,nj,np1)
                XP3(nj)=XP(1,1,nj,np2)
              ENDDO !nj
              CALL MESH_ANGLE(angle,XP1,XP2,XP3,ERROR,*9999)
              STATS(2,ne)=angle !temporary storage of angle
              BRANCHES(3,N)=angle*180.d0/PI !store the branching angle to parent
              
              ntotal=ntotal+1
              IF(diameters(ne0).GT.0.d0.AND.diameters(ne).GT.0.d0)THEN
                IF(diameters(ne)/diameters(ne0).LE.1.d0)THEN
                  num_ddp=num_ddp+1
                ENDIF
              ENDIF
              IF(diameters(ne0).GE.4.d0)THEN
                NBINS(1)=NBINS(1)+1
                BINS(1)=BINS(1)+angle
              ELSE IF(diameters(ne0).GE.2.d0)THEN
                NBINS(2)=NBINS(2)+1
                BINS(2)=BINS(2)+angle
              ELSE IF(diameters(ne0).GE.1.d0)THEN
                NBINS(3)=NBINS(3)+1
                BINS(3)=BINS(3)+angle
              ELSE IF(diameters(ne0).GE.0.7d0)THEN
                NBINS(4)=NBINS(4)+1
                BINS(4)=BINS(4)+angle
              ENDIF
            ELSE
              BRANCHES(3,N)=undefined
            ENDIF
            
            STATS(5,ne)=lengths(ne)
            ne_next=ne
            IF(NXI(1,0,ne_next).EQ.1)THEN
              ne_next=NXI(1,1,ne_next)
              STATS(5,ne)=STATS(5,ne)+lengths(ne_next)
            ENDIF
c            STATS(5,ne)=length_ne !temporary storage of summed length
c            L_D_mean=L_D_mean+length_ne/(diameters(ne))
c            L_D(1,ngen)=L_D(1,ngen)+length_ne/(diameters(ne))
c            L_D(2,norder)=L_D(2,norder)+length_ne/(diameters(ne))
c            L_D(3,nsorder)=L_D(3,nsorder)+length_ne/(diameters(ne))
          ENDIF
          IF(NXI(1,0,ne).EQ.0)THEN
            N_TERMINAL(INDEX(1))=N_TERMINAL(INDEX(1))+1
          ENDIF
C... Geometric properties of mesh
          BRANCHES(4,N)=undefined !initialise to no rotation angle
          IF(NXI(-1,0,ne).GT.0.AND.NXI(1,0,ne).GT.1)THEN
            ne0=NXI(-1,1,ne)
            IF(NXI(1,0,ne0).GT.1)THEN
              CALL MESH_PLANE_ANGLE(NBJ,ne,NPNE,NXI,angle,XP,ERROR,
     &          *9999)
              BRANCHES(4,N)=angle*180.d0/PI !rotation angle
            ENDIF
          ENDIF
        ENDDO !noelem
        
        N_BR=N
        CALL ASSERT(NMM.GE.21,'>>Increase NMM to 21',ERROR,*9999)
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          DO j=11,21 !initialise values for the summary statistics
            STATS(j,ne)=undefined
          ENDDO !j
          ne0=NXI(-1,1,ne)
          IF(ne0.NE.0)THEN
            IF(NORD(1,ne0).NE.NORD(1,ne))THEN
              IF(STATS(5,ne)/STATS(5,ne0).LE.1.d0)THEN
                num_llp=num_llp+1
              ENDIF
c              IF(STATS(5,ne)/(diameters(ne)).LT.8.d0)THEN
              STATS(19,ne)=STATS(5,ne)/STATS(5,ne0) !L/Lparent
              IF(diameters(ne0).GT.0.d0.AND.diameters(ne).GT.0.d0)THEN
                STATS(16,ne)=diameters(ne)/diameters(ne0) !D/Dparent
              ENDIF
c              ENDIF
            ENDIF
          ENDIF
          IF(NXI(1,0,ne).GE.2)THEN !'bi'furcations only
            ne1=NXI(1,1,ne) !first child
            ne2=NXI(1,2,ne) !second child
            
C.......... Summary statistics
            IF(STATS(6,ne1).LT.undefined.AND.STATS(6,
     &        ne2).LT.undefined)THEN
              IF(STATS(6,ne1).GE.STATS(6,ne2))THEN !diameter classification
                ne_major=ne1
                ne_minor=ne2
              ELSE
                ne_major=ne2
                ne_minor=ne1
              ENDIF
              IF(STATS(2,ne_minor).LT.undefined.AND.STATS(2,
     &          ne_major).LT.undefined)THEN
                STATS(11,ne)=STATS(2,ne_minor)*180.d0/PI
                STATS(12,ne)=STATS(2,ne_major)*180.d0/PI
              ENDIF
              
              IF(diameters(ne_minor).GT.0.d0.AND.
     &          diameters(ne_major).GT.0.d0)THEN
                STATS(13,ne)=STATS(5,ne_minor)/diameters(ne_minor) !L/D minor
                STATS(14,ne)=STATS(5,ne_major)/diameters(ne_major) !L/D major
                STATS(15,ne)=diameters(ne_minor)/diameters(ne_major) !minor D / major D
                STATS(17,ne)=diameters(ne_minor)/diameters(ne) !minor D / D parent
                STATS(18,ne)=diameters(ne_major)/diameters(ne) !major D / D parent
c               write(*,*) NORD(1,ne),STATS(17,ne),STATS(18,ne)
              ENDIF
              IF(STATS(5,ne1).LE.STATS(5,ne2))THEN !length classification
                ne_major=ne1
                ne_minor=ne2
              ELSE
                ne_major=ne2
                ne_minor=ne1
              ENDIF !length criteria
c              IF(STATS(5,ne_minor)/(diameters(ne_minor)).LT.8.d0.AND.STATS(5,
c     &          ne_major)/(diameters(ne_major)).LT.8.d0)THEN
              STATS(20,ne)=STATS(5,ne_major)/STATS(5,ne_minor)
c              ENDIF
            ENDIF
            
            IF(DOP)THEN
              WRITE(OP_STRING,'('' ne'',I5,''  angle'',F6.2,'
     &          //'''  L/Dmin'',F5.2,''  L/Dmaj'',F5.2)')
     &          ne,STATS(2,ne)*180.d0/PI,STATS(13,ne),STATS(14,ne)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF !DOP
          ENDIF !NXI
        ENDDO !noelem
        
C... Calculate mean branching statistics from values in BRANCHES
        DO N=1,N_BR
          DO i=1,NUM_SCHEMES !for generations, Horsfield orders, Strahler orders
            INDEX(i)=NBRANCHES(i,N)
C.......... length and diameter            
            DO j=1,2
              SUM_MEAN(i,j,INDEX(i))=SUM_MEAN(i,j,INDEX(i))+
     &          BRANCHES(j,N)
              IF(i.EQ.3.AND.j.EQ.1)THEN
                IF(INDEX(i).NE.NBRANCHES(5,N))THEN !not same as parent 
                  NTALLY(i,j,INDEX(i))=NTALLY(i,j,INDEX(i))+1
                ENDIF
              ELSE
                NTALLY(i,j,INDEX(i))=NTALLY(i,j,INDEX(i))+1
              ENDIF
            ENDDO !j
C.......... branching angle and rotation angle            
            DO j=3,4
              IF(BRANCHES(j,N).LT.undefined)THEN
                SUM_MEAN(i,j,INDEX(i))=SUM_MEAN(i,j,INDEX(i))
     &            +BRANCHES(j,N)
                NTALLY(i,j,INDEX(i))=NTALLY(i,j,INDEX(i))+1
              ENDIF
            ENDDO !j
C.......... ratio of L:D            
            j=5
            IF(BRANCHES(j,N).LT.undefined)THEN
              SUM_MEAN(i,j,INDEX(i))=SUM_MEAN(i,j,INDEX(i))
     &          +BRANCHES(j,N)
              NTALLY(i,j,INDEX(i))=NTALLY(i,j,INDEX(i))+1
            ENDIF
          ENDDO !i
          
C........ Summary statistics from BRANCHES
          DO j=3,5 !branching angle, rotation angle, L/D
            IF(BRANCHES(j,N).LT.undefined)THEN
              MEANS(j-2)=MEANS(j-2)+BRANCHES(j,N)
            ENDIF
          ENDDO !j
          
        ENDDO !N
        
        DO N=1,GENM
          DO j=3,5
            NTOTALN(j-2)=NTOTALN(j-2)+NTALLY(1,j,N)
          ENDDO !j
        ENDDO !N
        
        DO N=1,GENM
          DO i=1,NUM_SCHEMES
            DO j=1,5
              IF(NTALLY(i,j,N).GT.0)THEN
                SUM_MEAN(i,j,N)=SUM_MEAN(i,j,N)/DBLE(NTALLY(i,j,N))
                NMAX_GEN(i)=N
              ELSE
                SUM_MEAN(i,j,N)=0.d0
              ENDIF
            ENDDO !j
          ENDDO !i
        ENDDO !N
        
        DO N=1,5
          IF(NBINS(N).NE.0)THEN
            BINS(N)=BINS(N)/DBLE(NBINS(N))*180.d0/PI
          ENDIF
        ENDDO !N
        
C...... Summary statistics from BRANCHES
        DO j=3,5 !branching angle, rotation angle, L/D
          IF(NTOTALN(j-2).NE.0)THEN
            MEANS(j-2)=MEANS(j-2)/DBLE(NTOTALN(j-2))
          ELSE
            MEANS(j-2)=0.d0
          ENDIF
        ENDDO !j
        
        i=2 !Horsfield orders
        j=6 !Nw/Nw-1
        DO N=1,GENM-1
          IF(NTALLY(i,1,N).GT.0.AND.NTALLY(i,1,N+1).GT.0)THEN
            SUM_MEAN(i,j,N)=DBLE(NTALLY(i,1,N))/DBLE(NTALLY(i,1,N
     &        +1))
          ELSE
            SUM_MEAN(i,j,N)=0.d0
          ENDIF
        ENDDO !N
        
C...... Summary statistics from CE
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          DO j=11,21
            IF(STATS(j,ne).LT.undefined)THEN
              MEANS(j-7)=MEANS(j-7)+STATS(j,ne)
              NTOTALN(j-7)=NTOTALN(j-7)+1
            ENDIF
          ENDDO !j
        ENDDO !noelem
        
        DO j=11,21
          IF(NTOTALN(j-7).GT.0)THEN
            MEANS(j-7)=MEANS(j-7)/DBLE(NTOTALN(j-7))
          ENDIF
        ENDDO !j
C... End of mean calculation
        
C... Calculate the standard deviations
C...... sum of (value-mean)^2
        DO N=1,N_BR
          DO i=1,NUM_SCHEMES !for generations, Horsfield orders, Strahler orders
            INDEX(i)=NBRANCHES(i,N)
            DO j=1,5 !length, diameter, branching angle, rotation angle, L/D
              IF(BRANCHES(j,N).LT.undefined)THEN
                SD(i,j,INDEX(i))=SD(i,j,INDEX(i))
     &            +(BRANCHES(j,N)-SUM_MEAN(i,j,INDEX(i)))**2
              ENDIF
            ENDDO !j
          ENDDO !i
          DO j=3,5 !branching angle, rotation angle, L/D
            IF(BRANCHES(j,N).LT.undefined)THEN
              SDT(j-2)=SDT(j-2)
     &          +(BRANCHES(j,N)-MEANS(j-2))**2
            ENDIF
          ENDDO !j
        ENDDO !N
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          DO j=11,21
            IF(STATS(j,ne).LT.undefined)THEN
              SDT(j-7)=SDT(j-7)+(STATS(j,ne)-MEANS(j-7))**2
            ENDIF
          ENDDO !j
        ENDDO !noelem
        
C...... SD = sqrt(1/(n-1)*sum)
        DO N=1,GENM
          DO i=1,NUM_SCHEMES !for generations, Horsfield orders, Strahler orders
            DO j=1,5 !length, diameter, branching angle, rotation angle, L/D
              IF(NTALLY(i,j,N).GT.1)THEN
                SD(i,j,N)=DSQRT(SD(i,j,N)/DBLE(NTALLY(i,j,N)-1))
              ELSE
                SD(i,j,N)=0.d0
              ENDIF
            ENDDO !j
          ENDDO !i
        ENDDO !N
        DO j=1,13
          IF(NTOTALN(j).GT.1)THEN
            SDT(j)=DSQRT(SDT(j)/DBLE(NTOTALN(j)-1))
          ENDIF
        ENDDO !j
C.. End of standard deviation calculation        
        
C... Output statistics
        average_term_gen=0.d0
        SUM_TERM=0
        WRITE(OP_STRING,'(/'' Generation  #branches  #terminal'
     &    //'   Length          Diameter         Branching'
     &    //'        Rotation         ratio L:D'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''                        branches'
     &    //'     (mm)             (mm)           angle(deg)'
     &    //'      angle(deg)'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(115(''-''))')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        
        i=1
        DO N=1,NMAX_GEN(i)
          WRITE(OP_STRING,'(3(I10),5(F8.2,'' ('',F6.2,'')''))')
     '      N,NTALLY(i,1,N),N_TERMINAL(N),SUM_MEAN(i,1,N),SD(i,1,N),
     '      SUM_MEAN(i,2,N),SD(i,2,N),SUM_MEAN(i,3,N),SD(i,3,N),
     '      SUM_MEAN(i,4,N),SD(i,4,N),SUM_MEAN(i,5,N),SD(i,5,N)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          average_term_gen=average_term_gen+N_TERMINAL(N)*N
          SUM_TERM=SUM_TERM+N_TERMINAL(N)
        ENDDO
        IF(SUM_TERM.GT.0)THEN
          average_term_gen=average_term_gen/DBLE(SUM_TERM)
        ELSE
          average_term_gen=0.d0
        ENDIF
        
        WRITE(OP_STRING,'(/'' Horsfield   #branches     Length'
     '    //'           Diameter       Branching        Rotation'
     '    //'         ratio L:D      Nw/Nw-1'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    order                    (mm)'
     '    //'              (mm)         angle(deg)'
     &    //'     angle(deg)'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(115(''-''))')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        
        i=2
        DO N=1,NMAX_GEN(i)
          WRITE(OP_STRING,'(2(I10),5(F8.2,'' ('',F6.2,'')''),F8.2)')
     '      N,NTALLY(2,1,N),SUM_MEAN(i,1,N),SD(i,1,N),
     '      SUM_MEAN(i,2,N),SD(i,2,N),SUM_MEAN(i,3,N),SD(i,3,N),
     '      SUM_MEAN(i,4,N),SD(i,4,N),SUM_MEAN(i,5,N),SD(i,5,N),
     '      SUM_MEAN(i,6,N)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        
        WRITE(OP_STRING,'(/''   Strahler  #branches    Length'
     '    //'          Diameter        Branching'
     '    //'        Rotation          ratio L:D'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''    order                  (mm)'
     '    //'             (mm)          angle(deg)'
     '    //'      angle(deg)'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(115(''-''))')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        i=3
        DO N=1,NMAX_GEN(i)
          WRITE(OP_STRING,'(2(I10),5(F8.2,'' ('',F6.2,'')''))')
     '      N,NTALLY(3,1,N),SUM_MEAN(i,1,N),SD(i,1,N),
     '      SUM_MEAN(i,2,N),SD(i,2,N),SUM_MEAN(i,3,N),SD(i,3,N),
     '      SUM_MEAN(i,4,N),SD(i,4,N),SUM_MEAN(i,5,N),SD(i,5,N)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO  
        
        IF(DIAM_STRAHLER) THEN
          WRITE(OP_STRING,
     '      '(/'' Diam-Def Strahler  #branches    Length'
     '      //'          Diameter        Branching'
     '      //'        Rotation          ratio L:D'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''    order                  (mm)'
     '      //'             (mm)          angle(deg)'
     '      //'      angle(deg)'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(115(''-''))')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          
          i=4
          DO N=1,NMAX_GEN(i)
            WRITE(OP_STRING,'(2(I10),5(F8.2,'' ('',F6.2,'')''))')
     '        N,NTALLY(i,1,N),SUM_MEAN(i,1,N),SD(i,1,N),
     '        SUM_MEAN(i,2,N),SD(i,2,N),SUM_MEAN(i,3,N),SD(i,3,N),
     '        SUM_MEAN(i,4,N),SD(i,4,N),SUM_MEAN(i,5,N),SD(i,5,N)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO            
        ENDIF !DIAM_STRAHLER
        
        DO i=2,NUM_SCHEMES !Horsfield and Strahler orders
          DO N=1,NMAX_GEN(i)
            X(N)=N
            YREGRESS(N,1)=DLOG10(DBLE(NTALLY(i,1,N)))
            YREGRESS(N,2)=DLOG10(SUM_MEAN(i,1,N))
            YREGRESS(N,3)=DLOG10(SUM_MEAN(i,2,N))
          ENDDO !N
          DO j=1,3 !number of branches, length, diameter
            CALL LINREGRESS(NMAX_GEN(i),R_SQ(i,j),slope,X,
     '        YREGRESS(1,j),ERROR,*9999)
            RATIOS(i,j)=10.d0**DABS(slope)
          ENDDO !j
        ENDDO !i
        WRITE(OP_STRING,'(/''SUMMARY OF MEAN GEOMETRY STATISTICS''
     &    )')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(60(''-''))')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' terminal generation  = '',F7.3,'
     '    //'/'' branching angle      = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' rotation angle       = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' minor angle          = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' major angle          = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' L/D                  = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' L/D minor child      = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' L/D major child      = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' minor D/major D      = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' D/Dparent            = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' %D/Dparent  < 1      = '',F7.3,'
     '    //'/'' Dmin/Dparent         = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' Dmaj/Dparent         = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' L/Lp                 = '',F7.3,'' ('',F6.3,'')'','
     '    //'/'' %L/Lp < 1            = '',F7.3,'
     '    //'/'' L1/L2 (L1 < L2)      = '',F7.3,'' ('',F6.3,'')'')')
        
     '    average_term_gen,MEANS(1),SDT(1),MEANS(2),SDT(2),MEANS(4),
     '    SDT(4),MEANS(5),SDT(5),MEANS(3),SDT(3),MEANS(6),SDT(6),
     '    MEANS(7),SDT(7),MEANS(8),SDT(8),MEANS(9),SDT(9),
     &    DBLE(num_ddp)/DBLE(ntotal)*100.d0,MEANS(10),SDT(10),
     &    MEANS(11),SDT(11),MEANS(12),SDT(12),DBLE(num_llp)
     &    /DBLE(NELIST(0)-1)*100.d0,MEANS(13),SDT(13)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        
        WRITE(OP_STRING,'('' Rb Strahler          = '',F7.3,'
     '    //''' Rsq ='',F6.3,'//
     '    '/'' Rl Strahler          = '',F7.3,'' Rsq = '',F6.3,'//
     '    '/'' Rd Strahler          = '',F7.3,'' Rsq = '',F6.3,'//
     '    '/'' Rb Horsfield         = '',F7.3,'' Rsq = '',F6.3,'//
     '    '/'' Rl Horsfield         = '',F7.3,'' Rsq = '',F6.3,'//
     '    '/'' Rd Horsfield         = '',F7.3,'' Rsq = '',F6.3)')
     '    RATIOS(3,1),R_SQ(3,1),RATIOS(3,2),R_SQ(3,2),RATIOS(3,3),
     '    R_SQ(3,3),RATIOS(2,1),R_SQ(2,1),RATIOS(2,2),R_SQ(2,2),
     '    RATIOS(2,3),R_SQ(2,3)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(DIAM_STRAHLER) THEN
          WRITE(OP_STRING,
     &      '('' Rb diam-def Strahler  ='',F7.3,'' Rsq = '',F6.3,'//
     &      '/'' Rl diam-def Strahler  ='',F7.3,'' Rsq = '',F6.3,'//
     &      '/'' Rd diam-def Strahler  ='',F7.3,'' Rsq = '',F6.3)')
     &      RATIOS(4,1),R_SQ(4,1),RATIOS(4,2),R_SQ(4,2),RATIOS(4,3),
     &      R_SQ(4,3)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        
        WRITE(OP_STRING,'('' mean angle Dp 4.0+   = '',F7.3,'//
     &    '/''  mean angle Dp 3.0+   = '',F7.3,'//
     &    '/''  mean angle Dp 2.0+   = '',F7.3,'//
     &    '/''  mean angle Dp 1.0+   = '',F7.3,'//
     &    '/''  mean angle Dp 0.7+   = '',F7.3)')
     '    (BINS(j),j=1,5)

        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(OP_COORD)THEN!output to Ohio/Boston format
        ne_offset=NELIST(1)-1 !make sure it comes back to zero
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          ne0=NXI(-1,1,ne)
          nb=NBJ(1,ne)
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          ngen=NORD(1,ne)
          length=0.d0
          DO nj=1,NJT
            length=length+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2
          ENDDO !nj
          lengths(ne)=DSQRT(length)
          nv1=NVJE(1,nb,nj_radius,ne)
          nv2=NVJE(2,nb,nj_radius,ne)
          radius=0.5d0*(XP(1,nv1,nj_radius,np1)+XP(1,nv2,nj_radius,
     &      np2))
          diameters(ne)=2.0d0*radius
          IF(OP_FORMAT.EQ.1)THEN
 !original format needed by Ashish
c            WRITE(OP_STRING,'(I6,7(D14.4))') ngen-1,XP(1,1,1,np1),
c     &        XP(1,1,2,np1),XP(1,1,3,np1),XP(1,1,1,np2),XP(1,1,2,np2),
c     &        XP(1,1,3,np2),2.d0*radius
            WRITE(OP_STRING,'(I6,7(D14.4),2(I6))') ngen-1,XP(1,1,1,
     &        np1),XP(1,1,2,np1),XP(1,1,3,np1),XP(1,1,1,np2),XP(1,1,
     &        2,np2),XP(1,1,3,np2),2.d0*radius,ne,ne0
c            WRITE(OP_STRING,'(2(I6),6(D14.4))') noelem,ngen,XP(1,1,1,
c     '        np1),XP(1,1,2,np1),XP(1,1,3,np1),XP(1,1,1,np2),XP(1,1,2,
c     '        np2),XP(1,1,3,np2)
          ELSE IF(OP_FORMAT.EQ.2)THEN !Boston format
            !Jennines file in this case the Boston Format has
            ! been modified from the general case
	    
            WRITE(OP_STRING,'(2(I6),8(D14.4),3(I6),1(D14.4))')
     &     ne-ne_offset,ne0-ne_offset,lengths(ne),diameters(ne),
     &       XP(1,1,1,np1),XP(1,1,2,np1),XP(1,1,3,np1),XP(1,1,1,np2),
     &       XP(1,1,2,np2),XP(1,1,3,np2),ngen,NORD(5,ne),NORD(2,ne),
     &       XP(1,1,nj_alveoli,np2)
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !noelem(ne)
      ELSE IF(LISOLN)THEN
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          ne0=NXI(-1,1,ne)
          IF(ne0.NE.0) lengths(ne)=lengths(ne)+lengths(ne0)
        ENDDO !noelem
        
        nh=NH_LOC(1,nx)
        nb=NBJ(1,ne)
        ne=NELIST(1)
        np=NPNE(2,nb,ne)
        ny=NYNP(1,1,nh,np) 
        WRITE(OP_STRING,'(2(D14.4))') 0.d0,YP(ny,1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne)
          ny=NYNP(1,1,nh,np)
          WRITE(OP_STRING,'(2(D14.4))') lengths(ne),YP(ny,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !noelem
        
      ENDIF !OP_COORD
      
        
      CALL EXITS('OPMESH2')
      RETURN
 9999 CALL ERRORS('OPMESH2',ERROR)
      CALL EXITS('OPMESH2')
      RETURN 1
      END


