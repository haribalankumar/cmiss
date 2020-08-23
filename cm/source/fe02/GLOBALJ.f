      SUBROUTINE GLOBALJ(NBJ,NEELEM,NKJ,NPLIST,NPNE,NPNODE,ERROR,*)

C#### Subroutine: GLOBALJ
C###  Description:
C###    GLOBALJ calculates nodal arrays
C###    NBJ(nj,ne) and NPNE(nn,nb,ne)

C**** 23-6-92 AJP.  NKJ was originally set up when defining nodes and
C**** its value could only be increased here.  Now the values are reset
C**** to 1 and all are recalculated in this subroutine.
C**** 16/2/93 CPB NKJ no longer reset.  AJP overruled.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJ(NJM,NPM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nb1,ne,nj,nn,noelem,nonode,np,nr

      CALL ENTERS('GLOBALJ',*9999)

      DO nr=1,NRT
        ! We will use NPLIST as a means of storing
        ! whether a node is defined in this region
        ! Initialise NPLIST
        DO np=0,NPM
          NPLIST(np)=0
        ENDDO
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          NPLIST(np)=1
        ENDDO
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            IF(NPLIST(np).NE.0) THEN

              DO nj=1,NJ_LOC(0,0,nr)
                nb1=NBJ(nj,ne)
                IF(nb1.NE.0) THEN
                  IF(NKT(nn,nb1).GT.NKJ(nj,np)) THEN
                    NKJ(nj,np)=NKT(nn,nb1)
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('GLOBALJ')
      RETURN
 9999 CALL ERRORS('GLOBALJ',ERROR)
      CALL EXITS('GLOBALJ')
      RETURN 1
      END


