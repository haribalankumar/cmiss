      SUBROUTINE DEFINETUMOUR(NPNODE,radiation_diam,tumour_coord,
     &  tumour_diam,XP,RADIATION,ERROR,*)

C#### Subroutine: DEFINETUMOUR
C###  Description:
C###    DEFINETUMOUR finds nodes in the tumour and irradiated tissue
C###    spaces. Tumour/irradiated tissue volumes (diameter) and centre
C###    location are specified in the demesh command.
C###    
C###    radiation_diam = diameter of (spherical) irradiated tissue
C###    tumour_coord = fractional coordinate (relative to total node space)
C###    tumour_diam = diameter of (spherical) tumour

C***   Created by Annalisa Swan, Nov 2010.


      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
!       INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'

!     Parameter list
      INTEGER NPNODE(0:NP_R_M)
      REAL*8 radiation_diam,tumour_coord(3),tumour_diam,
     &  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL RADIATION
!     Local variables
      INTEGER i,nonode,np,NPLIST3(0:NP_R_M),NPLIST4(0:NP_R_M)
      REAL*8 dist_node,length(3),maxcoord(3),mincoord(3),
     &  tumour_centre(3)
      CHARACTER STRING*255

      CALL ENTERS('DEFINETUMOUR',*9999)


C...  Find min/max coordinate values for all nodes
      DO i=1,3
        maxcoord(i)=XP(1,1,i,NPNODE(1))
        mincoord(i)=XP(1,1,i,NPNODE(1))
      ENDDO !i
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO i=1,3
          IF(XP(1,1,i,np).GE.maxcoord(i))THEN
            maxcoord(i)=XP(1,1,i,np)
          ENDIF
          IF(XP(1,1,i,np).LE.mincoord(i))THEN
            mincoord(i)=XP(1,1,i,np)
          ENDIF
        ENDDO !i
      ENDDO !nonode

C...  Find tumour centre coordinates 
      DO i=1,3
        length(i)=maxcoord(i)-mincoord(i) !length in ith direction for all nodes
        tumour_centre(i)=mincoord(i)+(tumour_coord(i)*length(i)) !actual coordinate of tumour centre
      ENDDO !i


C...  Find nodes that are (1) within tumour radius and (2) outside tumour and within radiation radius
      NPLIST3(0)=0
      NPLIST4(0)=0
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        dist_node=0.0d0 !distance to node from tumour centre
        DO i=1,3
          dist_node=dist_node+(tumour_centre(i)-XP(1,1,i,np))**2.0d0
        ENDDO !i
        dist_node=SQRT(dist_node)
C       If tumour radius is greater than distance to node then node is inside tumour
        IF(dist_node.LE.(0.5d0*tumour_diam))THEN 
          NPLIST3(0)=NPLIST3(0)+1
          NPLIST3(NPLIST3(0))=np
        ENDIF
        IF(RADIATION)THEN !check whether node is within radiation radius
          IF((dist_node.GT.(0.5d0*tumour_diam)).AND. !greater than tumour radius
     &      (dist_node.LE.(0.5d0*radiation_diam)))THEN !less than radiation radius
            NPLIST4(0)=NPLIST4(0)+1
            NPLIST4(NPLIST4(0))=np
          ENDIF
        ENDIF !RADIATION
      ENDDO
      CALL ASSERT(NPLIST3(0).NE.0,
     &  '>> No nodes occuring in tumour location. Redefine tumour.',
     &  ERROR,*9999)
      IF(RADIATION)THEN
        CALL ASSERT(NPLIST4(0).NE.0,
     &    '>> No nodes occuring in radiation volume. Redefine.',
     &    ERROR,*9999)
      ENDIF

C...  Create node groups
      STRING='TUMOUR'! tumour node group name
      CALL GRNODE_SUB(NPLIST3,STRING,.TRUE.,ERROR,*9999) !AS=.TRUE.=1
      STRING='IRRADIATED'! irradiated node group name
      write(*,*) 'RADIATION=',RADIATION
      IF(RADIATION) CALL GRNODE_SUB(NPLIST4,STRING,.TRUE.,ERROR,*9999) !AS=.TRUE.=1


      CALL EXITS('DEFINETUMOUR')
      RETURN
 9999 CALL ERRORS('DEFINETUMOUR',ERROR)
      CALL EXITS('DEFINETUMOUR')
      RETURN 1
      END

