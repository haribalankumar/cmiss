      SUBROUTINE GRGAUS(GRNGLIST,NBJ,NEELEM,NELIST,
     '  NGLIST,NKJE,NPF,NPNE,NRLIST,NVJE,PG,SE,XA,
     '  XE,XG,XP,STRING,ERROR,*)

C#### Subroutine: GRGAUS
C###  Description:
C###    GRGAUS groups Gauss points.
C**** NTGRGA is number of Gauss groups currently defined.
C**** LAGRGA(nogrga) is label given to group number NOGRGA.
C**** LIGRGA(0,nogrga) is number of elements in the group number NOGRGA.
C**** LIGRGA(1...,nogrga) contains the following ordered list for each element
C**** LIGRGA(1,nogrga) is the element number.
C**** LIGRGA(2,nogrga) is the number of points grouped in the element
C**** LIGRGA(3..,nogrga) is the list of local Gauss point numbers in the element
C**** Created by Carey Stevens May 1999

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER GRNGLIST(0:NEGM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NGLIST(0:NGM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IBEG2,IEND,IEND2,ILISTMBR,
     '  index,index2,N3CO,max_index,nb,NBJ_temp(12),
     '  ne,ne_grouped,ne_num_points,ne_num_points_grouped,
     '  ng_grouped,ng,njj1,njj2,noelem,noelem2,nogrga,
     '  N1GRGA,nj,nr,nrr,NTR,num_points
      REAL*8 distance,origin(3),radius(1),Z(3)
      CHARACTER CHAR*30,CHAR2*2,LABEL*30
      LOGICAL ALL_REGIONS,CBBREV,FOUND,NEW_ELEM

      CALL ENTERS('GRGAUS',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR2,'(I2)') NTGRGR+1
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM group gauss
C###  Parameter:     <element (ELEMENT#s/NAME/all)[all]>
C###    Specify a list or name of a group of elements which contain
C###    the Gauss points to group
C###  Parameter:     <gauss (GAUSS_PTS#s/all)[all]>
C###    Specify a list of the Gauss points to group
C###  Parameter:     <as LABEL[gauss_1]>
C###    Specifies the name of the gauss group
C###  Parameter:     <region #[1]>
C###    Specify the region in which to group the Gauss points
C###  Description:
C###    Groups Gauss points into a group.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<element (ELEMENT#s/NAME/all)'
     '    //'[all]>'
        OP_STRING(3)=BLANK(1:15)//'<gauss (GAUSS_PTS#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//
     '    '<as LABEL[gauss_1'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(5)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group gauss function sphere
C###  Parameter:     <origin X#[0.0],Y#[0.0],Z#[0.0]>
C###   Specify the origin of the sphere
C###  Parameter:     <radius #[1.0]>
C###   Specify the radius of the sphere
C###  Parameter:     <as LABEL[gauss_1]>
C###    Specifies the name of the gauss group
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the regions in which to group the Gauss points
C###    The all value specifies all currently defined regions.
C###  Description:
C###    Group Gauss point numbers within a user specified sphere
C###    into a group.

        OP_STRING(1)=STRING(1:IEND)//' function sphere'
        OP_STRING(2)=BLANK(1:15)//'<origin X#[0.0],Y#[0.0],Z#[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<radius #[1.0]>'
        OP_STRING(4)=BLANK(1:15)//
     '    '<as LABEL[gauss_1'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group gauss ungrouped
C###  Parameter:     <element (ELEMENT#s/NAME/all)[all]>
C###    Specify a list or name of a group of elements which contain
C###    the Gauss points to group
C###  Parameter:     <as LABEL[gauss_1]>
C###    Specifies the name of the gauss group
C###  Parameter:     <region #[1]>
C###    Specify the region in which to group the Gauss points
C###  Description:
C###    Groups Gauss points that are not already in a group
C###    into a group.

        OP_STRING(1)=STRING(1:IEND)//' ungrouped'
        OP_STRING(2)=BLANK(1:15)//'<element (ELEMENT#s/NAME/all)'
     '    //'[all]>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<as LABEL[gauss_1'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','GRGAUS',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GAUSS_PT_MATERIALS.NE.0,
     '    '>>Set USE_GAUSS_PT_MATERIALS to 1',ERROR,*9999)
        IF(CBBREV(CO,'FUNCTION',2,noco,NTCO,N3CO)) THEN
C         Group Gauss points by function
          noco=N3CO
          IF(CBBREV(CO,'SPHERE',1,noco+1,NTCO,N3CO)) THEN
            IF(CBBREV(CO,'RADIUS',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),1,NTR,radius,ERROR,*9999)
            ELSE
              radius(1)=1.0d0
            ENDIF
            IF(CBBREV(CO,'ORIGIN',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTR,origin,ERROR,*9999)
            ELSE
              origin(1)=0.d0
              origin(2)=0.d0
              origin(3)=0.d0
            ENDIF
            CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,
     '        ALL_REGIONS,ERROR,*9999)
            DO nrr=1,NRLIST(0)
              nr=NRLIST(nrr)
              CALL PARSE_ELEMENTS(NEELEM,NELIST,noco-1,NRLIST,NTCO,CO,
     '          ERROR,*9999)
              CALL ASSERT(NELIST(0).NE.0,
     '          '>>No elements were found.',ERROR,*9999)
C             Set up GRNGLIST
              num_points=0
              index=1
              GRNGLIST(0)=0
              DO noelem=1,NELIST(0)
                ne=NELIST(noelem)
                FOUND=.FALSE.
                ne_num_points=0
                DO njj1=1,3 !geom/fibres/field
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    CALL ASSERT(nj.LE.12,
     '                '>>ERROR: increase size of NBJ_temp '
     '                //'to NJM',ERROR,*9999)
                    NBJ_temp(nj)=NBJ(nj,ne)
                  ENDDO !njj2
                ENDDO !njj1
                CALL XPXE(NBJ_temp,NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                DO ng=1,NGT(NBJ(1,ne))
                  CALL XEXG(NBJ_temp,ng,nr,PG,XE,XG,ERROR,*9999)
                  CALL XZ(ITYP10(nr),XG(1,1),Z)
                  distance=0.0d0
                  njj1=1 !geom
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    distance=distance+ABS(origin(nj)-Z(nj))**2
                  ENDDO !njj2
                  distance=SQRT(distance)
                  IF(distance.LE.radius(1)) THEN
                    num_points=num_points+1
                    ne_num_points=ne_num_points+1
                    FOUND=.TRUE.
                    GRNGLIST(index)=ne
                    GRNGLIST(index+1)=ne_num_points
                    GRNGLIST(index+1+ne_num_points)=ng
                  ENDIF
                ENDDO !ng
                IF(FOUND) THEN
                  index=index+2+ne_num_points
                  GRNGLIST(0)=GRNGLIST(0)+1
                ENDIF
              ENDDO !noelem
            ENDDO !nrr
            max_index=index
          ENDIF

        ELSE IF(CBBREV(CO,'UNGROUPED',2,noco,NTCO,N3CO)) THEN
C         Group the Gauss points that are not already in a group
          noco=N3CO
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,
     '      ALL_REGIONS,ERROR,*9999)
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco-1,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            CALL ASSERT(NELIST(0).NE.0,
     '        '>>No elements were found.',ERROR,*9999)
            num_points=0
            GRNGLIST(0)=0
            index2=1
            DO noelem=1,NELIST(0)
              ne=NELIST(noelem)
              nb=NBJ(1,ne)
              ne_num_points_grouped=0
              NEW_ELEM=.FALSE.
              DO ng=1,NGT(nb)
C               Check if ng in ne is in any other group, if not add
C               it to this group
                FOUND=.FALSE.
                DO nogrga=1,NTGRGA
                  index=2
                  DO noelem2=1,ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),1)
                    ne_grouped=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),index)
                    index=index+1
                    ne_num_points=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),
     '                  index)
                    IF(ne_grouped.EQ.ne) THEN
                      DO i=index+1,index+ne_num_points
                        ng_grouped=ILISTMBR(%VAL(LIGRGA_PTR(nogrga)),i)
                        IF(ng_grouped.EQ.ng) then
                          FOUND=.TRUE.
                        ENDIF !ng_grouped.EQ.ng
                      ENDDO !i
                    ENDIF !ne_grouped.EQ.ne
                    index=index+ne_num_points+1
                  ENDDO !noelem2
                ENDDO !nogrga
                IF(.NOT.FOUND) THEN
                  ne_num_points_grouped=ne_num_points_grouped+1
                  GRNGLIST(index2)=ne
                  GRNGLIST(index2+1)=ne_num_points_grouped
                  GRNGLIST(index2+1+ne_num_points_grouped)=ng
                  num_points=num_points+1
                  NEW_ELEM=.TRUE.
                ENDIF !.NOT.FOUND
              ENDDO !ng
              IF(NEW_ELEM) THEN
                index2=index2+ne_num_points_grouped+2
                max_index=index2
                GRNGLIST(0)=GRNGLIST(0)+1
              ENDIF !NEW_ELEM
            ENDDO !noelem
          ENDDO !nrr

        ELSE
C         Group Gauss points listed on the command line
          noco=noco+1
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,
     '      ALL_REGIONS,ERROR,*9999)
          DO nrr=1,NRLIST(0)
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco-1,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            CALL ASSERT(NELIST(0).NE.0,
     '        '>>No elements were found.',ERROR,*9999)
C           Create Gauss point list
            NGLIST(0)=0
            IF(CBBREV(CO,'GAUSS',1,noco,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),
     '          NGLIST(1),ERROR,*9999)
            ENDIF
            IF(NGLIST(0).EQ.0) THEN
              nb=NBJ(1,NELIST(1))
              NGLIST(0)=NGT(nb)
              DO ng=1,NGLIST(0)
                NGLIST(ng)=ng
              ENDDO
            ENDIF
C           Set up GRNGLIST
            index=1
            DO noelem=1,NELIST(0)
              ne=NELIST(noelem)
              GRNGLIST(index)=ne
              index=index+1
              GRNGLIST(index)=NGLIST(0)
              index=index+1
              DO ng=1,NGLIST(0)
                GRNGLIST(index)=ng
                index=index+1
              ENDDO
            ENDDO
            max_index=index
            num_points=NELIST(0)*NGLIST(0)
            GRNGLIST(0)=NELIST(0)
          ENDDO
        ENDIF

        IF(CBBREV(CO,'AS',1,noco,NTCO,N3CO)) THEN
          !Check whether group name already exists
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),CHAR)
          N1GRGA=0
          DO nogrga=1,NTGRGA
            CALL CUPPER(LAGRGA(nogrga),LABEL)
            CALL STRING_TRIM(LABEL,IBEG2,IEND2)
            IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG2:IEND2)) THEN
              N1GRGA=nogrga !is existing group label ID
              GO TO 100
            ENDIF
          ENDDO
 100      IF(N1GRGA.EQ.0) THEN !need new group label
            NTGRGA=NTGRGA+1 !is new total #groups
            CALL ASSERT(NTGRGA.LE.GRGA_MAXGRP,
     '        '>>Can''t create any more groups. '//
     '        'Consider reusing an existing group.',ERROR,*9999)
            N1GRGA=NTGRGA
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            LAGRGA(N1GRGA)=CO(N3CO+1)(IBEG:IEND) !is new group label
          ENDIF

        ELSE
          NTGRGA=NTGRGA+1 !is new total
          CALL ASSERT(NTGRGA.LE.GRGA_MAXGRP,
     '      '>>Can''t create any more groups. '//
     '      'Consider reusing an existing group.',ERROR,*9999)
          N1GRGA=NTGRGA
          WRITE(CHAR2,'(I2)') N1GRGA
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          LAGRGA(N1GRGA)='gaus_'//CHAR2(IBEG2:IEND2) !new group label
        ENDIF
        NLIGRGA(N1GRGA)=num_points
        CALL ALLOCATE_MEMORY(max_index,0,
     '    INTTYPE,LIGRGA_PTR(N1GRGA),MEM_INIT,ERROR,*9999)
        CALL ILIST_COPY(max_index,GRNGLIST(0),%VAL(LIGRGA_PTR(N1GRGA)))
      ENDIF

      CALL EXITS('GRGAUS')
      RETURN
 9999 CALL ERRORS('GRGAUS',ERROR)
      CALL EXITS('GRGAUS')
      RETURN 1
      END


