      SUBROUTINE FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,NQXI,NQS,nqsc,
     '  ERROR,*)

C#### Subroutine: FIND_LATT_NEIJK
C###  Description:
C###    FIND_LATT_NEIJK searches through NLATNE for the element
C###    which holds the lattice point nlat and returns this
C###    global element number as ne, along with the lattice
C###    points ijk lattice indice numbers (the position of the
C###    lattice number within the element lattice cube), and nqsc
C###    the number of this elements grid scheme.      
      
      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
!     Parameter List        
      INTEGER ne,nlat,NLATNE(NEQM+1),i,j,k,NQS(NEQM),NQXI(0:NIM,NQSCM),
     &  nqsc
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ne_next
      LOGICAL NOTFOUND
      
      CALL ENTERS('FIND_LATT_NEIJK',*9999)

      NOTFOUND=.TRUE.

      ne=1
      
C     Find the NE which nlat is within
      DO WHILE(NOTFOUND)
        CALL ASSERT(ne.LE.NEQM,'>> Lattice point not in valid element',
     &    ERROR,*9999)
        IF(NLATNE(ne).EQ.0) THEN
          ne=ne+1
        ELSEIF (NLATNE(ne).LE.nlat) THEN
          ne_next=ne+1
          DO WHILE(NLATNE(ne_next).EQ.0)
            ne_next=ne_next+1
          ENDDO
          IF (NLATNE(ne_next).GT.nlat) THEN
            NOTFOUND=.FALSE.
          ELSE
            ne=ne_next
          ENDIF
        ENDIF
      ENDDO
      
      nqsc=NQS(ne)
      k=INT(((nlat-(NLATNE(ne)-1))-1)/(NQXI(1,nqsc)*NQXI(2,nqsc)))+1
      j=INT(((nlat-(NLATNE(ne)-1)-(k-1)*(NQXI(1,nqsc)*NQXI(2,nqsc)))-1)/
     &  NQXI(1,nqsc))+1
      i=nlat-(NLATNE(ne)-1)-((k-1)*NQXI(2,nqsc)+j-1)*NQXI(1,nqsc)      

      CALL EXITS('FIND_LATT_NEIJK')
      RETURN
 9999 CALL ERRORS('FIND_LATT_NEIJK',ERROR)
      CALL EXITS('FIND_LATT_NEIJK')
      RETURN 1      
      END
