      SUBROUTINE VERT_LATT_INDICES(IJK,NEVS,maxi,maxj,maxk,ERROR,*)

C#### Subroutine: VERT_LATT_INDICES
C###  Description:
C###    VERT_LATT_INDICES calculates the lattice indices for
C###    the vertices of an element. These lattice indices are
C###    then used by the lattice jump functions to correctly
C###    position grid points.      
      
      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn' 
!     Parameter List
      INTEGER IJK(8,3),NEVS(0:12),maxi,maxj,maxk
      CHARACTER ERROR*(*) 
!     Local Variables
      INTEGER IMASK(2,8),JMASK(2,8),KMASK(2,8),nn

      DATA IMASK /1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1/
      DATA JMASK /1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1/
      DATA KMASK /1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1/      

      CALL ENTERS('VERT_LATT_INDICES',*9999)

      DO nn=1,8
        IJK(nn,1) = IMASK(1,NEVS(nn+1))+IMASK(2,NEVS(nn+1))*maxi
        IJK(nn,2) = JMASK(1,NEVS(nn+1))+JMASK(2,NEVS(nn+1))*maxj
        IJK(nn,3) = KMASK(1,NEVS(nn+1))+KMASK(2,NEVS(nn+1))*maxk
      ENDDO
 
      CALL EXITS('VERT_LATT_INDICES')
      RETURN
 9999 CALL ERRORS('VERT_LATT_INDICES',ERROR)
      CALL EXITS('VERT_LATT_INDICES')
      RETURN 1
      END
      
