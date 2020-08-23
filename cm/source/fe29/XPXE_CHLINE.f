      SUBROUTINE XPXE_CHLINE(CHXE,njf,nl,NPL,nr,XP,ERROR,*)

C#### Subroutine: XPXE_CHLINE
C###  Description:
C###    XPXE_CHLINE populates XP into CHXE, the element nodal values,
C###    where all basis functions are converted to
C###    cubic-Hermite basis functions.
C**** Written by Duane Malcolm, 06 November 2003

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER njf,nk,nl,NPL(5,0:3,NLM),nr
      REAL*8 CHXE(NSM,4),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER nj
C      REAL*8 
C      LOGICAL 

      CALL ENTERS('XPXE_CHLINE',*9999)
      
      IF(NPL(1,1,nl).EQ.1)THEN ! Linear Lagrange basis
        ! Here the linear basis function is converted to a cubic-Hermite
        ! basis function and the nodal values are stored in CHXE
        CALL ASSERT(.FALSE.,'>> Line basis function not implemented',
     &    ERROR,*9999)
      ELSEIF(NPL(1,1,nl).EQ.2)THEN ! Quadratic Lagrange basis
        ! Here the quadratic basis function is converted to a cubic-Hermite
        ! basis function and the nodal values are stored in CHXE
        CALL ASSERT(.FALSE.,'>> Line basis function not implemented',
     &    ERROR,*9999)
      ELSEIF(NPL(1,1,nl).EQ.3)THEN ! Cubic Lagrange basis
        ! Here the cubic basis function is converted to a cubic-Hermite
        ! basis function and the nodal values are stored in CHXE
        CALL ASSERT(.FALSE.,'>> Line basis function not implemented',
     &    ERROR,*9999)
      ELSEIF(NPL(1,1,nl).EQ.4)THEN ! Cubic-Hermite basis
        ! Checking Xi direction of line and setting derivative required
        ! for the cubic-Hermite description of line
        IF(NPL(1,0,nl).EQ.1) nk=2
        IF(NPL(1,0,nl).EQ.2) nk=3
        IF(NPL(1,0,nl).EQ.3) nk=4
        
        ! Looping over number of dimension and getting nodal values 
        ! for line
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          CHXE(1,nj)=XP(1,1,NJ_LOC(NJL_GEOM,nj,nr),NPL(2,1,nl)) ! x
          CHXE(2,nj)=XP(nk,1,NJ_LOC(NJL_GEOM,nj,nr),NPL(2,1,nl)) ! dx/dxi
          CHXE(3,nj)=XP(1,1,NJ_LOC(NJL_GEOM,nj,nr),NPL(3,1,nl)) ! x
          CHXE(4,nj)=XP(nk,1,NJ_LOC(NJL_GEOM,nj,nr),NPL(3,1,nl)) ! dx/dxi
        ENDDO
        
        ! Looping over unset dimentions and setting nodal values to zero
        IF(NJ_LOC(NJL_GEOM,0,nr).LT.3)THEN
          DO nj=NJ_LOC(NJL_GEOM,0,nr)+1,3
            CHXE(1,nj)=0.0D0 ! x
            CHXE(2,nj)=0.0D0 ! dx/dxi
            CHXE(3,nj)=0.0D0 ! x
            CHXE(4,nj)=0.0D0 ! dx/dxi
          ENDDO
        ENDIF
        
        ! Copying field njf into CHXE
        CHXE(1,4)=XP(1,1,NJ_LOC(NJL_FIEL,njf,nr),NPL(2,1,nl)) ! f
        CHXE(2,4)=XP(nk,1,NJ_LOC(NJL_FIEL,njf,nr),NPL(2,1,nl)) ! df/dxi
        CHXE(3,4)=XP(1,1,NJ_LOC(NJL_FIEL,njf,nr),NPL(3,1,nl)) ! f
        CHXE(4,4)=XP(nk,1,NJ_LOC(NJL_FIEL,njf,nr),NPL(3,1,nl)) ! df/dxi
        
      ELSE
        CALL ASSERT(.FALSE.,'>> Line basis function not implemented',
     &    ERROR,*9999)
      ENDIF                                                             
      
      CALL EXITS('XPXE_CHLINE')
      RETURN
 9999 CALL ERRORS('XPXE_CHLINE',ERROR)
      CALL EXITS('XPXE_CHLINE')
      RETURN 1
      END


