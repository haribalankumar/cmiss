
      SUBROUTINE ORDER_NUMBER(n_horsfield,ne,NORD,NXI,STRAHLER,
     &  STRAHLER_ADD,ERROR,*)

C#### Subroutine: ORDER_NUMBER
C###  Description:
C###  ORDER_NUMBER defines the Strahler and Horsfield order number for a branch
C###  in a tree-like network.

C*** Created by K.Burrowes, August 2004.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter list
      INTEGER n_horsfield,ne,NORD(5,NE_R_M),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  STRAHLER,STRAHLER_ADD
      CHARACTER ERROR*(*)
!     Local variables      
      INTEGER n_children,ne2,noelem,temp1
      
      CALL ENTERS('ORDER_NUMBER',*9999)
      
      n_horsfield=MAX(NORD(2,ne),1)  
      n_children=NXI(1,0,ne) !number of child branches

      IF(n_children.EQ.1)THEN
        IF(NORD(1,NXI(1,1,ne)).EQ.0)  n_children=0
      ENDIF
      STRAHLER=0
      STRAHLER_ADD=1
      IF(n_children.GE.2)THEN !branch has two or more daughters
        STRAHLER=NORD(3,NXI(1,1,ne)) !first daughter
        DO noelem=1,n_children !for all daughters
          ne2=NXI(1,noelem,ne) !global element # of daughter
          temp1=NORD(2,ne2) !Horsfield order of daughter 
          IF(temp1.GT.n_horsfield) n_horsfield=temp1
          IF(NORD(3,ne2).LT.STRAHLER)THEN
            STRAHLER_ADD=0
          ELSE IF(NORD(3,ne2).GT.STRAHLER)THEN
            STRAHLER_ADD=0
            STRAHLER=NORD(3,ne2) !highest daughter
          ENDIF
        ENDDO !noelem (ne2)
        n_horsfield=n_horsfield+1 !Horsfield ordering
      ELSE IF(n_children.EQ.1)THEN
        ne2=NXI(1,1,ne) !global element # of daughter
        IF(NORD(5,ne).EQ.1.d0)THEN
          n_horsfield=NORD(2,ne2)+1
          STRAHLER_ADD=NORD(3,ne2)+1
        ELSE
          n_horsfield=NORD(2,ne2) !same order as daughter
          STRAHLER_ADD=NORD(3,ne2)
        ENDIF
      ENDIF !NXI
        
      
      CALL EXITS('ORDER_NUMBER')
      RETURN
 9999 CALL ERRORS('ORDER_NUMBER',ERROR)
      CALL EXITS('ORDER_NUMBER')
      RETURN 1
      END  
      
