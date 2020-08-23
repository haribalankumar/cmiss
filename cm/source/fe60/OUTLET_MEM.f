      SUBROUTINE OUTLET_MEM(FACE_TOT,MAX_LOC_FACES,NBJ,NENP,
     '  NPLIST,NPNE,NPNODE,nr,NVCB,NVCNODE,N_OUTLET,ERROR,*)

C#### Subroutine: OUTLET_MEM
C###  Description:
C###    OUTLET_MEM computes the amount of memory required for the outlet
C###    Poisson equation arrays

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'voro00.inc'
!     Parameter list
      INTEGER N_OUTLET,FACE_TOT,MAX_LOC_FACES,NBJ(NJM,NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER bnvc,nonode,np,elemcount,noelem,ne,nb,nn,ADJNODE,
     '  cnonode,INTCOUNT,SIMPLEXLIST(60),ADJNODELIST(0:30),FTOT,nfvl,cnp
      LOGICAL ADD,FOUND

      CALL ENTERS('OUTLET_MEM',*9999)

      N_OUTLET=0
      FACE_TOT=0
      MAX_LOC_FACES=0

C     ..Compute the nonode numbers for each node np
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NPLIST(np)=nonode
      ENDDO

      DO bnvc=1,NVBT

C       ..Outlet nodes only
        IF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN
          nonode=NVCB(1,bnvc)
          N_OUTLET=N_OUTLET+1
          np=NPNODE(nonode,nr)

C         ..Compute the boundary simplices only
          elemcount=0
          DO noelem=1,NENP(np,0,nr)
            ne=NENP(np,noelem,nr)
            nb=NBJ(1,ne)
            INTCOUNT=0
            DO nn=1,NNT(nb)
              ADJNODE=NPNE(nn,nb,ne)
              cnonode=NPLIST(ADJNODE)

C             ..Split up into if..elseif statement in case we get an
C               array out of bounds on NVCB
              IF(NVCNODE(TYPE,cnonode).EQ.INTERNAL.OR.
     '          NVCNODE(TYPE,cnonode).EQ.INTBOUN) THEN
                INTCOUNT=INTCOUNT+1
              ELSEIF(NVCNODE(TYPE,cnonode).EQ.BOUNDARY.AND.
     '            NVCB(BCTYPE,NVCNODE(MAP,cnonode)).EQ.WALL) THEN
                INTCOUNT=INTCOUNT+1
              ENDIF
            ENDDO

            IF(INTCOUNT.LT.NNT(nb)) THEN
              elemcount=elemcount+1
              SIMPLEXLIST(elemcount)=ne
            ENDIF
          ENDDO

          ADJNODELIST(0)=0
          DO noelem=1,elemcount
            ne=SIMPLEXLIST(noelem)

C           ..Loop over the element's local nodes to get the surrounding
C           nodes
            nb=NBJ(1,ne)
            DO nn=1,NNT(nb)
              ADJNODE=NPNE(nn,nb,ne)

              IF(ADJNODE.NE.np) THEN

                cnonode=NPLIST(ADJNODE)
C               ..Only IB-IB connections
                IF(NVCNODE(TYPE,cnonode).EQ.INTBOUN) THEN

C                 ..Make sure that we haven't done this node already
                  ADD=.TRUE.
                  FOUND=.FALSE.
                  nfvl=1
                  DO WHILE(.NOT.FOUND.AND.nfvl.LE.ADJNODELIST(0))
                    cnp=NPNODE(ADJNODELIST(nfvl),nr)
                    IF(ADJNODE.EQ.cnp) THEN
                      ADD=.FALSE.
                      FOUND=.TRUE.
                    ELSE
                      nfvl=nfvl+1
                    ENDIF
                  ENDDO
                ELSE
                  ADD=.FALSE.
                ENDIF

C               ..Increment arrays
                IF(ADD) THEN
                  FTOT=ADJNODELIST(0)+1
                  ADJNODELIST(FTOT)=cnonode
                  ADJNODELIST(0)=FTOT
                  IF(FTOT.GT.MAX_LOC_FACES) MAX_LOC_FACES=FTOT
                  FACE_TOT=FACE_TOT+1
                ENDIF
              ENDIF !intboun
            ENDDO !nn
          ENDDO !noelem
        ENDIF !outlet
      ENDDO !bnvc

      CALL EXITS('OUTLET_MEM')
      RETURN
 9999 CALL ERRORS('OUTLET_MEM',ERROR)
      CALL EXITS('OUTLET_MEM')
      RETURN 1
      END


