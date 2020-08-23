      SUBROUTINE READ_GEOMETRY(BDRY,INODE,IREGION,ITYPE,N_BDRY,N_IBDRY,
     '  N_INTNL,NPLIST,BC,NODE,XP,LGAPSQ,ERROR,*)

C#### Subroutine: READ_GEOMETRY
C###  Description:
C###    Reads boundary geometry from XP_B,XP_IB,XP_IN.
CC JMB 15-NOV-2001

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER BDRY(NBOUNDM),INODE(N_GM),IREGION(0:N_GM),ITYPE(N_GM),
     '  N_BDRY,N_IBDRY,N_INTNL,NPLIST(0:NPM)
      REAL*8 BC(LDBC,3),NODE(LDNODE,3),XP(NKM,NVM,NJM,NPM)
      LOGICAL LGAPSQ
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,no_b,np

      CALL ENTERS('READ_GEOMETRY',*9999)

      ! Initialization
      I=0
      LGAPSQ=.FALSE.
      DO no_b=1,N_BDRY
        i=i+1
        CALL ASSERT(i.LE.N_GM,'>>Increase N_GM in genmesh.cmn',
     '    ERROR,*9999)
        INODE(i)=BDRY(no_b) !location can have more than 1 B node
        ITYPE(i)=0
        np=NPLIST(no_b) !defines the B nodes
        DO j=1,3
          BC(i,j)=0.d0
c         NODE(i,j)=XP_B(j,no_b)
          NODE(i,j)=XP(1,1,j,np)
        ENDDO
        IREGION(i)=1
      ENDDO !no_b
      DO no_b=1,N_IBDRY
        i=i+1
        np=NPLIST(N_BDRY+no_b)
        DO j=1,3
c         NODE(i,j)=XP_IB(j,no_b)
          NODE(i,j)=XP(1,1,j,np)
        ENDDO
        IREGION(i)=2
      ENDDO !no_b
      DO no_b=1,N_INTNL
        i=i+1
        np=NPLIST(N_BDRY+N_IBDRY+no_b)
        DO j=1,3
c         NODE(i,j)=XP_IN(j,no_b)
          NODE(i,j)=XP(1,1,j,np)
        ENDDO !j
        IREGION(i)=3
      ENDDO !no_b
c      GAPSQ=GAPSQ !should be controlled from command line
      LGAPSQ=.TRUE.
      IREGION(0)=I

      CALL EXITS('READ_GEOMETRY')
      RETURN
 9999 CALL ERRORS('READ_GEOMETRY',ERROR)
      CALL EXITS('READ_GEOMETRY')
      RETURN 1
      END


