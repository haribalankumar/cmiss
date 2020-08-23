      SUBROUTINE NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*)

C#### Subroutine: NENXI
C###  Description:
C###    NENXI finds elements surrounding element ne.

C#### Variable: NXI(-ni:ni,0:nei,0:ne)
C###  Type: INTEGER
C###  Set_up: NENXI
C###  Description:
C###    <HTML> <PRE>
C###    NXI(-3,1,ne) is element in -Xi(3) direction
C###         3      "    "    "   Xi(3)     "
C###        -2      "    "    "  -Xi(2)     "
C###         2      "    "    "   Xi(2)     "
C###        -1      "    "    "  -Xi(1)     "
C###         1      "    "    "   Xi(1)     "
C###    NXI(.,1,0) is 0.
C###    <BR>
C###    For Element Branching:
C###
C###    NXI(.,0,ne) is the number of elements adjacent to ne
C###    NXI(.,nei,ne) are the adjacent element numbers
C###    <BR>
C###    <B>Note:</B> Element neighbours are listed for first the current region,
C###    then all other regions starting from number 1 to the highest region.
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NNB(4,4,4,NBFM),NPNE(NNM,NBFM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,idrn,IEND,ini,j,INPI(3),INPIT(3),nb,ne,
     '  NE_LIST(0:NEPM),ne2,nei,nei_found,nei0,NEIMAX,NEIMI,nep,ni,ni1,
     '  ni2,ni_check,ni_search,NIS(2),NITB,nn,nn1,nn2,noelem,np,np2,npi,
     '  NP_LIST(0:NNM),nr,nri,nr2,NXI0NR(0:NEPM,NRT)
      LOGICAL NOTFOUND,FACE_COLLAPSED(-NIM:NIM),KEEP,XI_COLLAPSED

      CALL ENTERS('NENXI',*9999)

C     This is probably in case ne is zero when using NXI
      DO ini=-NIM,NIM
        DO nei=0,1
          NXI(ini,nei,0)=0
        ENDDO
      ENDDO

      NEIMAX=0 !max nei

C     Loop over all defined elements and determine in which region they are.
      DO nr=1,NRT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          NITB=NIT(nb) !assume <= 3

C         Find the number of nodes in each direction.
C         Assumes that the largest number of nodes in each
C         direction can be found by searching the first line in NNB.
          DO ni=1,3
            INPI(ni)=1
          ENDDO !ni
          DO ni=1,3
            INPI(ni)=2
            NOTFOUND=.TRUE.
            DO WHILE(INPI(ni).LE.4.AND.NOTFOUND)
              IF(NNB(INPI(1),INPI(2),INPI(3),nb).EQ.0) THEN
                NOTFOUND=.FALSE.
              ELSE
                INPI(ni)=INPI(ni)+1
              ENDIF
            ENDDO !INPI(ni)
            INPIT(ni)=INPI(ni)-1
            INPI(ni)=1
          ENDDO !ni

C         Initialize.
C         The first element is initialized because many routines do not
C         check the number of adjacent elements.
          DO ini=-NIM,NIM
            DO nei=0,1
              NXI(ini,nei,ne)=0
            ENDDO
          ENDDO
C         Place the current element in NXI.
          NXI(0,0,ne)=1
          NXI(0,1,ne)=ne

C         Determine which faces are collapsed.
C         Loop over face normals.
          DO ni=1,NITB
C           Determine xi directions in surface
            NIS(1)=1+MOD(ni,3)
            NIS(2)=6-ni-NIS(1)
C           Loop over two faces with this face normal.
            INPI(ni)=1
            DO idrn=-1,1,2
              ini=idrn*ni
              FACE_COLLAPSED(ini)=.FALSE.
              DO j=1,2
                ni_check=NIS(j)
                IF(ni_check.LE.NITB) THEN
                  ni_search=NIS(3-j)
                  XI_COLLAPSED=.TRUE.
                  INPI(ni_search)=0
                  DO WHILE(INPI(ni_search).LT.INPIT(ni_search)
     '              .AND.XI_COLLAPSED)
                    INPI(ni_search)=INPI(ni_search)+1
                    INPI(ni_check)=1
                    nn1=NNB(INPI(1),INPI(2),INPI(3),nb)
                    INPI(ni_check)=2
                    nn2=NNB(INPI(1),INPI(2),INPI(3),nb)
                    IF(nn2.NE.0.AND.nn1.NE.0) THEN
                      IF(NPNE(nn2,nb,ne).NE.NPNE(nn1,nb,ne)) THEN
                        XI_COLLAPSED=.FALSE.
                      ENDIF !NPNE
                    ENDIF !nn
                  ENDDO !inpi(ni_search)
                  IF(XI_COLLAPSED) FACE_COLLAPSED(ini)=.TRUE.
                ENDIF !ni_check<=NITB
              ENDDO !j
              INPI(ni)=INPIT(ni) !for second loop
            ENDDO !idrn
          ENDDO !ni

C         Loop over the elements of NXI(*,,ne).
C         The zeroth ini is done first so that equivalent elements can
C         be removed from the adjacent element inis
          DO i=1,2*NITB+1

C           Get list of nodes that must match in this direction.
            IF(i.EQ.1) THEN
              ini=0
C             Look for equivalent element.  All nodes must match.
              NP_LIST(0)=NNT(nb)
              DO nn=1,NP_LIST(0)
                NP_LIST(nn)=NPNE(nn,nb,ne)
              ENDDO
            ELSE !i.ne.1
C             Look for adjacent element.  Get nodes on a face
              ni=i/2
              IF(i.EQ.2*ni) THEN
                ini=ni
                INPI(ni)=INPIT(ni)
              ELSE
                ini=-ni
                INPI(ni)=1
              ENDIF
C             Determine xi directions in surface
              NIS(1)=1+MOD(ni,3)
              NIS(2)=6-ni-NIS(1)
C             If the face is collapsed then don't look in this
C             direction.  The exception is if the opposite face is also
C             collapsed.  This may indicate that we have one of those
C             funny elements used in non-rc coords that goes
C             around the central axis back to itself.
              IF(FACE_COLLAPSED(ini).AND..NOT.FACE_COLLAPSED(-ini)) THEN
C               no nodes to search will be interpreted no search
                NP_LIST(0)=0
                NE_LIST(0)=0
              ELSE !not COLLAPSED
C               Collect nodes
                npi=0
                DO ni1=1,INPIT(NIS(1))
                  INPI(NIS(1))=ni1
                  DO ni2=1,INPIT(NIS(2))
                    INPI(NIS(2))=ni2
                    nn=NNB(INPI(1),INPI(2),INPI(3),nb)
                    IF(nn.NE.0) THEN
                      npi=npi+1
                      NP_LIST(npi)=NPNE(nn,nb,ne)
                    ENDIF !nn!=0
                  ENDDO !ni2
                ENDDO !ni1
                NP_LIST(0)=npi
              ENDIF !COLLAPSED
            ENDIF !i==1

            IF(NP_LIST(0).NE.0) THEN

C             Remove duplicate nodes to reduce number of searches
C             required.
              CALL ILISTRMDUP(NP_LIST(0),NP_LIST(1),ERROR,*9999)

C             size required for current list if greater than NEIM
              NEIMI=0

CC LKC 8-NOV-1999
CC   Solve the problem of neighbouring elements by searching first
CC   in the current region then search the other regions starting from 1.
CC   The order of the region
CC   numbers may/will effect the neighbours which are found.
CC
CC
CCC LKC 17-JUL-1999 We really only want to loop over adjacent elements in
CCC  the same region otherwise interfaces elements may be considered to
CCC  be adjacent.
CCC
CCC news mpn 14-Sep-95
CCC         loop over all regions for adjacent elements
CCC          DO nrr=1,NRT
CCC            DO noelem1=1,NEELEM(0,nrr)
CCC              nee=NEELEM(noelem1,nrr)
CCCC old          DO nee=1,NET(nr)
CCC          DO noelem1=1,NEELEM(0,nr)
CCC            nee=NEELEM(noelem1,nr)
CC
CC CS 18-OCT-1999 LKC comments not true. We do want to know about
CC  adjacent elements in different regions eg pressure coupling.
CC  Reverting back for now till a solution is agreed upon.
CC rgb lets solve the problem of element looping by instead
CC     going over all elements surrounding an element

              DO nri=1,NRT

C               Loop over current region first then other regions so
C               that a routine that only expects one adjacent element
C               will find one in the current region (if one exists).
                IF(nri.EQ.1) THEN
                  nr2=nr
                  IF(ini.EQ.0) THEN
                    NXI0NR(1,nr2)=ne !current list of equivalent elements
                    NXI0NR(0,nr2)=1
                  ENDIF
                ELSE
                  IF(nri.LE.nr) THEN
                    nr2=nri-1
                  ELSE
                    nr2=nri
                  ENDIF
                  IF(ini.EQ.0) NXI0NR(0,nr2)=0 !no known equivalent elements
                ENDIF

C               Set up a list of possible neighbouring elements from a
C               node.  Only elements of the same dimension are considered.
                np=NP_LIST(1)
C!!   We assume that each element is only in one region
                nei=1
                DO nep=1,NENP(np,0,nr2)
                  ne2=NENP(np,nep,nr2)
                  IF(NIT(NBJ(1,ne2)).EQ.NITB) THEN
                    NE_LIST(nei)=ne2
                    nei=nei+1
                  ENDIF
                ENDDO !nep
                NE_LIST(0)=nei-1

                IF(NXI0NR(0,nr2).NE.0) THEN
C                 Remove known equivalent elements from possibles.
C                 NXI0NR and NE_LIST are in increasing order.
                  nei=1
                  nep=1
                  nei0=1
                  DO WHILE(nep.LE.NE_LIST(0))
                    IF(nei0.GT.NXI0NR(0,nr2)) THEN
                      KEEP=.TRUE.
                    ELSE
C                     NE_LIST(nep) should never be > NXI0NR(nei0,nr2)
C                     due to equivalent elements always being in NE_LIST.
                      KEEP=NE_LIST(nep).NE.NXI0NR(nei0,nr2)
                    ENDIF
                    IF(KEEP) THEN
                      NE_LIST(nei)=NE_LIST(nep)
                      nep=nep+1
                      nei=nei+1
                    ELSE
                      nep=nep+1
                      nei0=nei0+1
                    ENDIF
                  ENDDO !nei0,nep
                  NE_LIST(0)=nei-1
                ENDIF !(NXI0NR(0,nr2).NE.0)

C               Check each other node
                DO npi=2,NP_LIST(0)
                  np2=NP_LIST(npi)
C                 Eliminate elements that are not connected to this
C                 node.
C                 NENP and NE_LIST are in increasing order.
                  nei=1
                  nei_found=0
                  nep=1
                  DO WHILE(nei.LE.NE_LIST(0).AND.nep.LE.NENP(np2,0,nr2))
                    IF(NE_LIST(nei).GT.NENP(np2,nep,nr2)) THEN
                      nep=nep+1
                    ELSEIF(NE_LIST(nei).LT.NENP(np2,nep,nr2)) THEN
                      nei=nei+1
                    ELSE !NE_LIST(nei).EQ.NENP(np2,nep,nr2)
                      nei_found=nei_found+1
                      NE_LIST(nei_found)=NE_LIST(nei)
                      nei=nei+1
                      nep=nep+1
                    ENDIF
                  ENDDO !nei,nep
                  NE_LIST(0)=nei_found
                ENDDO !np2

                IF(ini.EQ.0) THEN
C                 (Re)setup NXI0NR
                  nei=1
                  IF(nr2.EQ.nr) THEN
C                   need to include current element
                    NOTFOUND=.TRUE.
                    DO WHILE(nei.LE.NE_LIST(0).AND.NOTFOUND)
                      IF(NE_LIST(nei).LT.ne) THEN
                        NXI0NR(nei,nr2)=NE_LIST(nei)
                        nei=nei+1
                      ELSE
                        NOTFOUND=.FALSE.
                      ENDIF
                    ENDDO !nei
                    NXI0NR(nei,nr2)=ne
                    nei0=nei+1
                  ELSE
                    nei0=nei
                  ENDIF !nr2==nr
                  DO WHILE(nei.LE.NE_LIST(0))
                    NXI0NR(nei0,nr2)=NE_LIST(nei)
                    nei=nei+1
                    nei0=nei0+1
                  ENDDO !nei
                  NXI0NR(0,nr2)=nei0-1
                ENDIF !ini==0
C               Add the elements found to NXI
                IF(NXI(ini,0,ne)+NE_LIST(0).GT.NEIM) THEN
                  IF(NEIMI.EQ.0) THEN
                    NEIMI=NXI(ini,0,ne)+NE_LIST(0)
                    NE_LIST(0)=NEIM-NXI(ini,0,ne) !add only elements that fit
                  ELSE
                    NEIMI=NEIMI+NE_LIST(0)
                    NE_LIST(0)=0 !don't add any elements
                    IF(NEIMI.GT.NEIMAX) NEIMAX=NEIMI
                  ENDIF
                ENDIF !NEIM
                nei=NXI(ini,0,ne)
                DO nei_found=1,NE_LIST(0)
                  nei=nei+1
                  NXI(ini,nei,ne)=NE_LIST(nei_found)
                ENDDO
                NXI(ini,0,ne)=nei

              ENDDO !nri
            ENDIF !NP_LIST(0).NE.0

          ENDDO !i
        ENDDO !noelem
      ENDDO !nr
      IF(NEIMAX.GT.NEIM) THEN
        IEND=0
        CALL APPENDC(IEND,'Increase NEIM to ',ERROR)
        CALL APPENDI(IEND,NEIMAX,ERROR)
        CALL FLAG_ERROR(1,ERROR(:IEND))
        GOTO 9999
      ENDIF

      IF(DOP) THEN
        DO nr=1,NRT
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ini=-NIM,NIM
              WRITE(OP_STRING,
     '          '('' NXI('',I2,'',...,'',I5,''):'',7(X,I5))')
     '          ini,ne,(NXI(ini,nei,ne),nei=1,NXI(ini,0,ne))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !ini
          ENDDO !noelem
        ENDDO !nr
      ENDIF

      CALL EXITS('NENXI')
      RETURN
 9999 CALL ERRORS('NENXI',ERROR)
      CALL EXITS('NENXI')
      RETURN 1
      END


