      SUBROUTINE NENXI_VORO(NBJ,NEELEM,NENP,NPNE,NPNODE,nr,NXI,
     '  ERROR,*)

C#### Subroutine: NENXI_VORO
C###  Description:
C###    NENXI_VORO finds the element surrounding element ne. When
C###    USE_VORONOI is switched on and simplex elements are being used,
C###    this subroutine calculates the surrounding elements and loads
C###    them into NXI(0,nei,ne).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NXI(-NIM:NIM,0:NEIM,0:NEM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER OTHER_NODES(3),NE2_NODES(4),noelem,ne,j,nonode,np,
     '  nb,nn,nn_np,OTHER_NODE_TOT,np2,noelem2,ne2
      LOGICAL WORK(4),ENTRIES_MATCH

      CALL ENTERS('NENXI_VORO',*9999)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO j=0,NEIM
          NXI(0,j,ne)=0
        ENDDO
      ENDDO
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        DO noelem=1,NENP(np,0,nr)
          ne=NENP(np,noelem,nr)
          nb=NBJ(1,ne)
          OTHER_NODE_TOT=0
          DO nn=1,NNT(nb)
            IF(NPNE(nn,nb,ne).NE.np) THEN
              OTHER_NODE_TOT=OTHER_NODE_TOT+1
              OTHER_NODES(OTHER_NODE_TOT)=NPNE(nn,nb,ne)
            ELSE
              nn_np=nn
            ENDIF
          ENDDO
          np2=OTHER_NODES(1)
          DO noelem2=1,NENP(np2,0,nr)
            ne2=NENP(np2,noelem2,nr)
            IF(ne2.NE.ne) THEN
              DO nn=1,NNT(nb)
                NE2_NODES(nn)=NPNE(nn,nb,ne2)
              ENDDO
              IF(ENTRIES_MATCH(NE2_NODES,NNT(nb),OTHER_NODES,
     '          OTHER_NODE_TOT,WORK)) THEN
                NXI(0,0,ne)=NXI(0,0,ne)+1
                NXI(0,nn_np,ne)=ne2
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('NENXI_VORO')
      RETURN
 9999 CALL ERRORS('NENXI_VORO',ERROR)
      CALL EXITS('NENXI_VORO')
      RETURN 1
      END


