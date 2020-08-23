      SUBROUTINE CALC_NH_LOC(nhx_MAX,nx,ERROR,*)

C#### Subroutine: CALC_NH_LOC
C###  Description:
C###    CALC_NH_LOC sets up the NH_LOC and NH_TYPE arrays
C###    for the current nx.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER nhx_MAX,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER FREELIST(NH_LOC_MX),nh,nhx,NUMFREE,nx1,nxx

      CALL ENTERS('CALC_NH_LOC',*9999)

C     Set up NH_LOC(nhx,nx)
C     Clear out any existing nh's
      ! For nx>1, NH_LOC(0,nx) is always 0 so nothing gets cleared out.
      ! This breaks field fitting problems that use problem class 2
      ! if explicit mappings are defined because it means that NH_TYPE
      ! doesn't get cleared and causes it to use the next available
      ! nh value from problem class 1 instead of starting again at 1.
      DO nhx=1,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        NH_TYPE(nh,1)=0 !nhx for this nh
        NH_TYPE(nh,2)=0 !nx for this nh
        NH_LOC(nhx,nx)=0
      ENDDO !nhx
      NH_LOC(0,nx)=0
C     Reset max# nh's
      NH_LOC(0,0)=0
      DO nxx=1,NX_LIST(0)
        nx1=NX_LIST(nxx)
        DO nhx=1,NH_LOC(0,nx1)
          nh=NH_LOC(nhx,nx1)
          IF(nh.GT.NH_LOC(0,0)) NH_LOC(0,0)=nh
        ENDDO !nhx
      ENDDO !nxx (nx1)

C     Store free nh's in FREELIST
      ! If the nh's were cleared out above then why is this needed?
      nh=0
      NUMFREE=0
      DO WHILE((NUMFREE.NE.nhx_MAX).AND.(nh.LE.NHM))
        nh=nh+1
        IF(NH_TYPE(nh,1).EQ.0) THEN !this nh is not in use
          NUMFREE=NUMFREE+1
          CALL ASSERT(NUMFREE.LE.NH_LOC_MX,'>>Increase NH_LOC_MX in '
     '      //'loc00.cmn',ERROR,*9999)
          FREELIST(NUMFREE)=nh
        ENDIF
      ENDDO !while...
      CALL ASSERT(nh.LE.NHM,' >>Increase NHM',ERROR,*9999)
C     Store field in free space
      DO nhx=1,NUMFREE
        nh=FREELIST(nhx)
        NH_LOC(nhx,nx)=nh
        NH_TYPE(nh,1)=nhx
        NH_TYPE(nh,2)=nx
        IF(nh.GT.NH_LOC(0,0)) NH_LOC(0,0)=nh
      ENDDO !nhx
      NH_LOC(0,nx)=nhx_MAX !is max# nhx variables

      CALL EXITS('CALC_NH_LOC')
      RETURN
 9999 CALL ERRORS('CALC_NH_LOC',ERROR)
      CALL EXITS('CALC_NH_LOC')
      RETURN 1
      END


