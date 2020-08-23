
      SUBROUTINE BRANCH_ORD(NELIST,NORD,NXI,ERROR,*)
      
C#### Subroutine: BRANCH_ORD
C###  Description:
C###    BRANCH_ORD allocates orders to branches in pulmonary trees.
C###    Defines generations, Horsfield orders, Strahler orders.
      
C***  Created February, 2003.      
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
!     Parameter list
      INTEGER NELIST(0:NEM),NORD(5,NE_R_M),NXI(-NIM:NIM,0:NEIM,0:NEM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ne,ne0,n_generation,noelem,n_horsfield,STRAHLER,
     &  STRAHLER_ADD
      
      CALL ENTERS('BRANCH_ORD',*9999)

C.....Find all elements with no adjacent element in -Xi1
c      NELIST_START(0)=0
c      DO noelem=1,NELIST(0)
c        ne=NELIST(noelem)
c        IF(NXI(-1,0,ne).EQ.0)THEN
cC         No adjacent element in -Xi1 direction
c          NELIST_START(0)=NELIST_START(0)+1
c          NELIST_START(NELIST_START(0))=ne
c          NORD(1,ne)=1 !generation 1
c        ENDIF !NXI
c      ENDDO !noelem

C.....Calculate branch generations
      maxgen=1
      DO noelem=1,NELIST(0)
        ne=NELIST(noelem)
        ne0=NXI(-1,1,ne) !parent
        IF(ne0.NE.0)THEN
          n_generation=NORD(1,ne0) !parent generation
          IF(NXI(1,0,ne0).EQ.1)THEN !single daughter
            IF(NORD(5,ne).EQ.0)THEN !same generation
              NORD(1,ne)=n_generation
            ELSE IF(NORD(5,ne).EQ.1)THEN !start of 'half' branch
              NORD(1,ne)=n_generation+1
            ENDIF
          ELSE IF(NXI(1,0,ne0).GE.2)THEN
            NORD(1,ne)=n_generation+1
          ENDIF
        ELSE
          NORD(1,ne)=1 !generation 1
        ENDIF
        maxgen=max(maxgen,NORD(1,ne))
      ENDDO

C.....Calculate the branch orders
      DO noelem=NELIST(0),1,-1
        ne=NELIST(noelem)
        !ORDER_NUMBER find the Horsfield and Strahler order number
        CALL ORDER_NUMBER(n_horsfield,ne,NORD,NXI,STRAHLER,STRAHLER_ADD,
     &    ERROR,*9999)
        
        NORD(2,ne)=n_horsfield !store the Horsfield order
        NORD(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
      ENDDO !noelem
      
      CALL EXITS('BRANCH_ORD')
      RETURN
 9999 CALL ERRORS('BRANCH_ORD',ERROR)
      CALL EXITS('BRANCH_ORD')
      RETURN 1
      END  

      
