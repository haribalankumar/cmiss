      SUBROUTINE EVORDR(NBJ,NEELEM,NELIST,NENP,NORD,NPNE,NPNODE,NRLIST,
     ' NXI,ERROR,*)

C#### Subroutine: EVORDR
C###  Description:
C###    EVORDR calculates generations, Horsfield orders, and Strahler
C###    orders for 1D trees.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IFROMC,N3CO,nb,noelem,nonode,np,np2,inlets,outlets,ne,nr,
     '  num_elems,i,nn
      LOGICAL ALL_REGIONS,CBBREV,DIAMETER_DEFINED,DISCONNECT,
     &  DUPLICATE
      CHARACTER STRING*255
      
      CALL ENTERS('EVORDR',*9999)

      IF(CO(noco+1).EQ.'?') THEN

C---------------------------------------------------------------------

C#### Command: FEM evaluate order
C###  Parameter:      <diameter_defined>
C###    Specifies that diameter-defined Strahler ordering will be
C###    calculated. 
C###  Parameter:      <elements (#s/all) [all]>
C###    Specifies the elements for which ordering will be done. 
C###    All elements with no adjacent element in the -xi1 direction are
C###    assumed to be generation 1.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Description:
C###    EVORDR evaluates the mesh ordering (generations, Horsfield
C###    ordering, Strahler ordering, diameter-defined Strahler ordering
C###    if specified).  

        OP_STRING(1)=BLANK(1:15)//'<diameter_defined>'
        OP_STRING(2)=BLANK(1:15)//'<elements (#s/all) [all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all) [all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVORDR',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        
        DIAMETER_DEFINED=.FALSE.
        IF(CBBREV(CO,'DIAMETER_DEFINED',4,noco+1,NTCO,N3CO)) THEN
          DIAMETER_DEFINED=.TRUE.
        ENDIF

C       Calculate generations, Horsfield orders, Strahler orders
        CALL BRANCH_ORD(NELIST,NORD,NXI,ERROR,*9999)
        IF(DIAMETER_DEFINED)THEN
C         Calculate diameter-based Strahler orders           
c          CALL GNDIAM()
        ENDIF

c       Check for disconnected nodes and number of inlets and outlets
        IF(CBBREV(CO,'REGION',4,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        DUPLICATE=.FALSE.
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(1,ne)
          np=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          IF(np.EQ.np2)THEN
            DUPLICATE=.TRUE.
            WRITE(OP_STRING,'('' Element '',I5,'
     &        //''' repeats nodes'')')ne
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
         ENDIF
        ENDDO
        
        CALL ASSERT(.NOT.DUPLICATE,'Remove duplicate nodes',
     &    ERROR,*9999)

        DISCONNECT=.FALSE.
        INLETS=0
        OUTLETS=0
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(1,ne)
          DO nn=1,2
          np=NPNE(nn,nb,ne)
          num_elems=NENP(np,0,nr)
          IF(num_elems.EQ.0)THEN
            DISCONNECT=.TRUE.
            WRITE(OP_STRING,'('' Node '',I5,'
     &        //''' not connected'')')np
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSEIF(num_elems.EQ.1)THEN
            ne=NENP(np,1,nr)
            IF(NXI(1,0,ne).EQ.0)OUTLETS=OUTLETS+1
            IF(NXI(-1,0,ne).EQ.0)then
                INLETS=INLETS+1
                write(*,*) 'inlet order',ne,NORD(3,ne)
           endif
          ELSEIF(num_elems.GT.3)THEN
            WRITE(OP_STRING,'('' Node '',I5,'
     &        //''' is attached to '',i5,'
     &        //''' elements'')')np,num_elems
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          ENDDO
        ENDDO

        WRITE(OP_STRING,'('' Number of inlets '',i5)')inlets
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Number of outlets '',i5)')outlets
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        CALL ASSERT(.NOT.DISCONNECT,'>>Disconnected nodes',ERROR,*9999)
        IF(OUTLETS.GT.0)THEN 
                NTERMINAL=OUTLETS
                CALL SET_USER_INTEGER('NTERMINAL',NTERMINAL,ERROR)
        ENDIF

      ENDIF
      
      CALL EXITS('EVORDR')
      RETURN
 9999 CALL ERRORS('EVORDR',ERROR)
      CALL EXITS('EVORDR')
      RETURN 1
      END



