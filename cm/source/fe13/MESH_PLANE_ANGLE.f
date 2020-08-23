      SUBROUTINE MESH_PLANE_ANGLE(NBJ,ne,NPNE,NXI,ANGLE,XP,ERROR,*)

C#### Subroutine: MESH_PLANE_ANGLE
C###  Description:
C###  MESH_PLANE_ANGLE calculates the rotation angle between the
C###  branching planes of a parent (ne) and its child branches.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER NBJ(NJM,NEM),ne,NPNE(NNM,NBFM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 ANGLE,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb,ne0,ne1,nj,np1,np2,np3,np4,np5
      REAL*8 norm_1(4),norm_2(4),SCALAR,XPOINT(3,5),temp
      REAL*8 CALC_ANGLE

      CALL ENTERS('MESH_PLANE_ANGLE',*9999)

      CALL ASSERT(NXI(-1,0,ne).GT.0,
     '  '>>Branch has no parent, coding error',ERROR,*9999)

      nb=NBJ(1,ne)
      ne0=NXI(-1,1,ne) !parent element
      np1=NPNE(1,nb,ne) !start node
      np2=NPNE(2,nb,ne) !end node
      ne1=NXI(1,1,ne0)
      IF(DOP)THEN
        WRITE(OP_STRING,'(5(I6))') nb,ne0,np1,np2,ne1
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(ne1.EQ.ne) ne1=NXI(1,2,ne0) !sibling
      np3=NPNE(2,nb,ne1) !end node of sibling
      np4=NPNE(2,nb,NXI(1,1,ne)) !end node of first child
      np5=NPNE(2,nb,NXI(1,2,ne)) !end node of second child
      DO nj=1,3
        XPOINT(nj,1)=XP(1,1,nj,np1)
        XPOINT(nj,2)=XP(1,1,nj,np2)
        XPOINT(nj,3)=XP(1,1,nj,np3)
        XPOINT(nj,4)=XP(1,1,nj,np4)
        XPOINT(nj,5)=XP(1,1,nj,np5)
      ENDDO !nj
      CALL PLANE_FROM_3_PTS(norm_1,2,XPOINT(1,1),XPOINT(1,2),
     '  XPOINT(1,3),ERROR,*9999)
      CALL NORMALISE2(3,norm_1,temp,ERROR,*9999) !unit vector
      CALL PLANE_FROM_3_PTS(norm_2,2,XPOINT(1,2),XPOINT(1,4),
     '  XPOINT(1,5),ERROR,*9999)
      CALL NORMALISE2(3,norm_2,temp,ERROR,*9999) !unit vector
      ANGLE=CALC_ANGLE(norm_1,norm_2)
c      ANGLE=SCALAR(3,norm_1,norm_2)
c      ANGLE=MAX(-1.d0,ANGLE)
c      ANGLE=MIN(1.d0,ANGLE)
c      ANGLE=DACOS(DABS(ANGLE))

      CALL EXITS('MESH_PLANE_ANGLE')
      RETURN
 9999 CALL ERRORS('MESH_PLANE_ANGLE',ERROR)
      CALL EXITS('MESH_PLANE_ANGLE')
      RETURN 1
      END


