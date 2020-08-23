      INTEGER FUNCTION INT_LATT_JUMP(i,j,k,MINS,MAXS,VALS,COL)

C#### Function: INT_LATT_JUMP
C###  Type: INTEGER
C###  Description:
C###    INT_LATT_JUMP calls the three child functions
C###    LATT_JUMP1, LATT_JUMP2 and LATT_JUMP3 and interpolates
C###    the results. INT_LATT_JUMP is called by
C###    the lattice grid point generation subroutine
C###    to discretely interpolate values in i,j and k.      
C###    LATT_JUMP1 is used of collapsing occurs along
C###    a particular direction. LATT_JUMP2 is used
C###    if there is no collapsing. LATT_JUMP3 is used
C###    if three points on a face are collapsed.
C###    The variable COL is 1 if the interpolation is across
C###    vertices, 0 if it is across elements.
      
      IMPLICIT NONE      
!     Parameter List
      INTEGER i,j,k,MINS(3),MAXS(3),VALS(8),COL
!     Local Variables
      INTEGER L1,L2,L3,L4,L12,L34,FTOV(4,6),count,valmax,nf3,f,
     '  numatvalmax
      
!     Functions
      INTEGER LATT_JUMP1,LATT_JUMP2,LATT_JUMP3

      DATA FTOV /1,3,5,7,
     '  2,4,6,8,
     '  1,2,5,6,
     '  3,4,7,8,
     '  1,2,3,4,
     '  5,6,7,8/
      
C!!!  DAH 18/7/02 Assert calls need to be added to ensure
C!!!              that no invalid types of collapse are made.    

      valmax=0
      nf3=0
      IF(COL.EQ.1) THEN
        DO count=1,8
          IF(valmax.LT.VALS(count)) THEN
            valmax=VALS(count)
          ENDIF
        ENDDO !count
        DO f=1,6
          numatvalmax=0
          IF(VALS(FTOV(1,f)).EQ.valmax) numatvalmax=numatvalmax+1
          IF(VALS(FTOV(2,f)).EQ.valmax) numatvalmax=numatvalmax+1
          IF(VALS(FTOV(3,f)).EQ.valmax) numatvalmax=numatvalmax+1
          IF(VALS(FTOV(4,f)).EQ.valmax) numatvalmax=numatvalmax+1
          IF (numatvalmax.EQ.3) THEN
            nf3=nf3+1
          ENDIF
        ENDDO !f
      ENDIF
C     Interpolate in the i direction
      IF(VALS(1).EQ.VALS(2)) THEN
        L1=VALS(1)
      ELSE
        L1=LATT_JUMP2(i,MINS(1),MAXS(1),VALS(1),VALS(2))
      ENDIF
      IF(VALS(3).EQ.VALS(4)) THEN
        L2=VALS(3)
      ELSE
        L2=LATT_JUMP2(i,MINS(1),MAXS(1),VALS(3),VALS(4))
      ENDIF
      IF(VALS(5).EQ.VALS(6)) THEN
        L3=VALS(5)
      ELSE
        L3=LATT_JUMP2(i,MINS(1),MAXS(1),VALS(5),VALS(6))
      ENDIF
      IF(VALS(7).EQ.VALS(8)) THEN
        L4=VALS(7)
      ELSE
        L4=LATT_JUMP2(i,MINS(1),MAXS(1),VALS(7),VALS(8))
      ENDIF

C     Interpolate in the j direction
      IF(VALS(1).EQ.VALS(3).OR.VALS(2).EQ.VALS(4)) THEN
        IF(nf3.GE.2) THEN
          L12=LATT_JUMP3(j,MINS(2),MAXS(2),L1,L2)
        ELSE
          L12=LATT_JUMP1(j,MINS(2),MAXS(2),L1,L2)
        ENDIF
      ELSE
        L12=LATT_JUMP2(j,MINS(2),MAXS(2),L1,L2)
      ENDIF
      IF(VALS(5).EQ.VALS(7).OR.VALS(6).EQ.VALS(8)) THEN
        IF(nf3.GE.2) THEN
          L34=LATT_JUMP3(j,MINS(2),MAXS(2),L3,L4)
        ELSE        
          L34=LATT_JUMP1(j,MINS(2),MAXS(2),L3,L4)
        ENDIF
      ELSE
        L34=LATT_JUMP2(j,MINS(2),MAXS(2),L3,L4)
      ENDIF
      IF(VALS(1).EQ.VALS(5).OR.VALS(2).EQ.VALS(6).OR.
     '  VALS(3).EQ.VALS(7).OR.VALS(4).EQ.VALS(8)) THEN
        IF(nf3.GE.2) THEN
          INT_LATT_JUMP=LATT_JUMP3(k,MINS(3),MAXS(3),L12,L34)
        ELSE         
          INT_LATT_JUMP=LATT_JUMP1(k,MINS(3),MAXS(3),L12,L34)
        ENDIF
      ELSE
        INT_LATT_JUMP=LATT_JUMP2(k,MINS(3),MAXS(3),L12,L34)
      ENDIF
      RETURN
      END
      
