      SUBROUTINE DIAM_DEF_STRAHLER(MIN_ORDER,NEELEM,NORD,nr,NXI,CE,
     &  D_BRANCH,ERROR,*)
      
C#### Subroutine: DIAM_DEF_STRAHLER
C###  Description:
C###    DIAM_DEF_STRAHLER allocates diameter defined Strahler orders to
C###    branches in pulmonary trees.
      
C***  Created February, 2003.      
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
!     Parameter list
      INTEGER MIN_ORDER,NEELEM(0:NE_R_M,0:NRM),NORD(5,NE_R_M),nr,
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),D_BRANCH(*)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER gen_next,MAX_ORD,ne,ne0,ngen,noelem,ord
      REAL*8 D(GENM),D1,D2,diameters(NE_R_M),D_PREV(GENM),DIAM,ERR,
     &  L_D_RATIO,N(GENM),SD(GENM),SUM_LENGTH
      LOGICAL CONTINU,FIRST
      
      CALL ENTERS('DIAM_DEF_STRAHLER',*9999)
      
      DO ngen=1,GENM
        D_PREV(ngen)=0.d0 !for diameter-defined ordering
      ENDDO !ngen
      MAX_ORD=0 
      DO noelem=NEELEM(0,nr),1,-1
        ne=NEELEM(noelem,nr)
        ord=MAX(NORD(2,ne),1)  
        SUM_LENGTH=CE(1,ne)
        IF(NXI(1,0,ne).NE.0)THEN
          ne0=NXI(1,1,ne)
          gen_next=NORD(1,ne0) !generation
          DO WHILE(gen_next.EQ.NORD(1,ne))
            SUM_LENGTH=SUM_LENGTH+CE(1,ne0)
            IF(NXI(1,1,ne0).NE.0)THEN
              gen_next=NORD(1,NXI(1,1,ne0))
              ne0=NXI(1,1,ne0)
            ELSE
              gen_next=NORD(1,ne)+1
            ENDIF
          ENDDO
        ENDIF
        IF(NXI(-1,0,ne).NE.0)THEN
          ne0=NXI(-1,1,ne)
          gen_next=NORD(1,ne0)
          DO WHILE(gen_next.EQ.NORD(1,ne))
            SUM_LENGTH=SUM_LENGTH+CE(1,ne0)
            IF(NXI(-1,1,ne0).NE.0)THEN
              gen_next=NORD(1,NXI(-1,1,ne0))
              ne0=NXI(-1,1,ne0)
            ELSE
              gen_next=NORD(1,ne)-1
            ENDIF
          ENDDO
        ENDIF
        DIAM=SUM_LENGTH/L_D_RATIO
        
        IF(DIAM.LT.D_BRANCH(NORD(3,ne)+MIN_ORDER))THEN
          DIAM=D_BRANCH(NORD(3,ne)+MIN_ORDER)
          diameters(ne)=DIAM
        ENDIF
      ENDDO !noelem
      
      FIRST=.TRUE.
C... Iterate until converged
      ERR=1.d0
      CONTINU=.TRUE.
      DO WHILE(ERR.GT.0.1d0)
C... Calculate diam average (D) & standard deviation (SD) in each order
        DO ngen=1,GENM
          D(ngen)=0.d0 !initialise 
          SD(ngen)=0.d0
          N(ngen)=0.d0
        ENDDO !ngen
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(FIRST) NORD(4,ne)=NORD(3,ne) !original Strahler into NORD(4,ne)
          ord=NORD(4,ne) !diameter-defined Strahler order
          DIAM=diameters(ne)
          D(ord)=D(ord)+DIAM
          N(ord)=N(ord)+1.d0 !sums number of branches in order
          MAX_ORD=MAX(MAX_ORD,ord)
        ENDDO !noelem
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          ord=NORD(4,ne)
          DIAM=CE(4,ne)*2.d0
          SD(ord)=SD(ord)+(DIAM-(D(ord)/N(ord)))**2.d0
        ENDDO
        DO ord=1,MAX_ORD
          IF(N(ord).GT.0.d0) THEN
            SD(ord)=DSQRT(SD(ord)/N(ord))
            D(ord)=D(ord)/N(ord) !average diameter
          ELSE IF(N(ord).EQ.0.d0) THEN
            WRITE(OP_STRING,'('' No branches in order: '',I5)')
     '        ord
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            CONTINU=.FALSE.
          ENDIF
        ENDDO
C... Implement diameter criteria (from Jiang:1994)
        IF(CONTINU) THEN
          DO noelem=NEELEM(0,nr),1,-1
            ne=NEELEM(noelem,nr)
            ord=NORD(4,ne) !Strahler order
            IF(ord.GT.1.AND.ord.LT.MAX_ORD) THEN
              D1=((D(ord-1)+SD(ord-1))+(D(ord)-SD(ord)))/2.d0
 !lower diameter limit
              D2=((D(ord)+SD(ord))+(D(ord+1)-SD(ord+1)))/2.d0
 !upper diameter limit
              DIAM=diameters(ne)
C... Vessel considered to be order ord if it's diam is between D1 & D2
              IF(DIAM.GE.D1.AND.DIAM.LE.D2) THEN
                NORD(4,ne)=ord
              ELSE IF(DIAM.LT.D1) THEN
                NORD(4,ne)=ord-1
              ELSE IF(DIAM.GT.D2) THEN
                NORD(4,ne)=ord+1
              ENDIF
            ENDIF
          ENDDO !noelem
          ERR=0.d0
          DO ord=1,MAX_ORD
            ERR=ERR+DABS((D_PREV(ord)-D(ord)))
            D_PREV(ord)=D(ord) !average diameter for order
          ENDDO
        ELSE !if N(ord).EQ.0 exit
          ERR=0.d0
        ENDIF !CONTINU
        FIRST=.FALSE.
      ENDDO !WHILE
      
      CALL EXITS('DIAM_DEF_STRAHLER')
      RETURN
 9999 CALL ERRORS('DIAM_DEF_STRAHLER',ERROR)
      CALL EXITS('DIAM_DEF_STRAHLER')
      RETURN 1
      END  

      
