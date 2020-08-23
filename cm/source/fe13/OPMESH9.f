      SUBROUTINE OPMESH9(ERROR,*)

C#### Subroutine: OPMESH9
C###  Description:
C###    OPMESH! outputs regular mesh data.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mesh01.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,m,n,ni,nj
      CHARACTER OPT2(4)*20

      DATA OPT2(1)/'Even spacing'/,
     '     OPT2(2)/'Spacing specified by blocks'/,
     '     OPT2(3)/'Spacing specified by position'/,
     '     OPT2(4)/'Unused'/

      CALL ENTERS('OPMESH9',*9999)
      WRITE(OP_STRING,'(/'' Mesh type is '',A)') OPT2(MESH1_TYPE)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'( '' Basis function# for mesh is '',I1)')
     '  MESH1_NB
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      WRITE(OP_STRING,'(/'' Coordinates of corner nodes:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO N=1,2**NIT(MESH1_NB) !loops over corner nodes
        WRITE(OP_STRING,'( ''  Corner node '',I1,'':'',3D12.4)')
     '    n,(MESH1_COORD(n,nj),nj=1,NJT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      WRITE(OP_STRING,'(/'' Element numbers:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO ni=1,NIT(MESH1_NB)
        IF(MESH1_TYPE.EQ.1) THEN !even spacing
          WRITE(OP_STRING,'( '' #elements in s('',I1,'
     '      //''') direction = '',I3)') ni,MESH1_S(ni,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(MESH1_TYPE.EQ.2) THEN !spacing specified by blocks
          WRITE(OP_STRING,'(/'' #elements in 3 blocks in s('',I1,'
     '      //''') direction = '',3I3)') ni,(MESH1_S(ni,i),i=1,3)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' ratio of coarse to fine spacing ='''
     '      //',I3)') MESH1_R(ni)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(MESH1_TYPE.EQ.3) THEN !spacing specified by position
          WRITE(OP_STRING,'( '' #elements in s('',I1,'
     '      //''') direction = '',I3)') ni,MESH1_S(ni,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Relative positions along s('','
     '      //'I1,''): '',25D12.4)')
     '      ni,(MESH1_X(ni,m),m=0,MESH1_S(ni,0))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ENDIF
      ENDDO

      CALL EXITS('OPMESH9')
      RETURN
 9999 CALL ERRORS('OPMESH9',ERROR)
      CALL EXITS('OPMESH9')
      RETURN 1
      END



