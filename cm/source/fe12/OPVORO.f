      SUBROUTINE OPVORO(NFVC,NODENVC,NPNODE,
     '  nr,VC,XNFV,XP,ERROR,*)

C#### Subroutine: OPVORO
C###  Description:
C###    OPVORO lists Voronoi cell information

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),NPNODE(0:NP_R_M,0:NRM),
     '  nr
      REAL*8 VC(0:NVCM),XNFV(-(NJM+1):NJM,NFVM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nfvl,nj,nonode,np,nvc

      CALL ENTERS('OPVORO',*9999)

      WRITE(OP_STRING,'('' Listing Voronoi cells:'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)
        CALL ASSERT( np.NE.0, '>> Global node not found. Check ' //
     '    ' that the correct region is specified',ERROR,*9999)
        WRITE(OP_STRING,'(/$,''Voronoi cell:'',I6)') nvc
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''Node:        '',
     '    I6,'' (nonode = '',I6,'')'')') np,nonode
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''Coordinates:    '',3D12.4)')
     '    (XP(1,1,nj,np),nj=1,NJ_LOC(NJL_GEOM,0,nr))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/$,''Adjacency Molecule:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''==================='')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        
C TVK 06/04/99 Changed output to accomdate 2D and 3D
        
        IF(NJT.EQ.2) THEN
          WRITE(OP_STRING,'(''  Node  Face        Area'//
     '      '        Dist  Centroid         '//
     '      '       Norm'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nfvl=1,NFVC(1,0,nvc)
            WRITE(OP_STRING,
     '        '(I6,I6,D12.4,D12.4,D12.4,D12.4,2D12.4)')
     '        NFVC(1,nfvl,nvc),NFVC(2,nfvl,nvc),
     '        XNFV(FAREA,NFVC(2,nfvl,nvc)),
     '        1.0D0/XNFV(IDIST,NFVC(2,nfvl,nvc)),
     '        XNFV(-2,NFVC(2,nfvl,nvc)),
     '        XNFV(-3,NFVC(2,nfvl,nvc)),
     '        (XNFV(nj,NFVC(2,nfvl,nvc)),nj=1,NJ_LOC(NJL_GEOM,0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ELSE
          WRITE(OP_STRING,'(''  Node  Face        Area'//
     '      '        Dist  Centroid                     '//
     '      '       Norm'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nfvl=1,NFVC(1,0,nvc)
            WRITE(OP_STRING,
     '        '(I6,I6,D12.4,D12.4,D12.4,D12.4,D12.4,3D12.4)')
     '        NFVC(1,nfvl,nvc),NFVC(2,nfvl,nvc),
     '        XNFV(FAREA,NFVC(2,nfvl,nvc)),
     '        1.0D0/XNFV(IDIST,NFVC(2,nfvl,nvc)),
     '        XNFV(-2,NFVC(2,nfvl,nvc)),
     '        XNFV(-3,NFVC(2,nfvl,nvc)),
     '        XNFV(-4,NFVC(2,nfvl,nvc)),
     '        (XNFV(nj,NFVC(2,nfvl,nvc)),nj=1,NJ_LOC(NJL_GEOM,0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF !NJT
        WRITE(OP_STRING,'(/$,''Volume:'',D12.4)') VC(nvc)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO
      WRITE(OP_STRING,'(/$,''Total Volume:'',D12.4)') VC(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL EXITS('OPVORO')
      RETURN
 9999 CALL ERRORS('OPVORO',ERROR)
      CALL EXITS('OPVORO')
      RETURN 1
      END
      
      
