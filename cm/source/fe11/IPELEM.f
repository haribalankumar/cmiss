      SUBROUTINE IPELEM(IBT,nbasis,NBJ,NEELEM,NELIST,NENP,NE_OFFSET,
     &  NKJE,NPNE,NPNODE,NP_OFFSET,nr,NRLIST,NRE,NVJE,NVJP,
     &  BASIS_OVERRIDE,TYPE,TREE,ERROR,*)

C#### Subroutine: IPELEM
C###  Description:
C###    IPELEM inputs element topology.

C**** 17-4-93  Element definitons for basis functions in the same family
C****          are now set up automatically - only the definitions for
C****          one basis function in the family needs to be entered.

C#### Variable: NKJE(nk,nn,nj,ne)
C###  Type: INTEGER
C###  Set_up: IPELEM,IPMESH1,IPMESH2,IPMESH3,IPMESH5
C###  Description:
C###    NKJE(nk,nn,nj,ne) is the global derivative number of local
C###    derivative nk of node nn of element ne for variable nj.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'fibr01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),nbasis,NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),NE_OFFSET,
     '  NKJE(NKM,NNM,NJM,NEM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NP_OFFSET,nr,NRE(NEM),NRLIST(0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      LOGICAL BASIS_OVERRIDE,TREE
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,ICHAR,IEND,IEND1,IEND2,INFO,mm,MP,
     '  n1elem,N2ELEM,nb,nb2,nb_actual,nb_extended,ne,ne2,nj,nj2,njj,
     &  njmax,njtype,njj2,nk,nn,noelem,noelem2,nonode,NOQUES,np,nr2,
     '  NTELEM,OC_COUNT,N3CO,NELIST(0:NEM)
      CHARACTER CHAR1*5,CHAR2*3,CHAR3*2,CHAR4*2,CHAR5*5,CHAR6*6
      LOGICAL CONTINU,DOBASIS,EXISTS,FIBREFIELD,FILEIP,FIN,
     '  FIRST,FIRST_VERSION,
     '  FOUND,NBLOG,NE_EXISTS,NP_EXISTS,SAMETYPE,SAMENUMXI,
     '  CBBREV

      CALL ENTERS('IPELEM',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C GDR 17/1/05 Make the "...;w elements #"option work  
      IF(IOTYPE.EQ.3) THEN
        IF(CBBREV(CO,'ELEMENTS',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,
     '      NTCO,CO,ERROR,*9999)
        ELSE
          DO noelem=1,NEELEM(0,nr)
            NELIST(noelem)=NEELEM(noelem,nr)
          ENDDO !noelem
          NELIST(0)=NEELEM(0,nr)
        ENDIF          
        NTELEM=NELIST(0)
        IDATA(1)=NTELEM
      ENDIF
      
      IF(TYPE(1:8).EQ.'GEOMETRY') THEN
        FORMAT='($,'' The number of elements is [1]: '',I5)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NEM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NTELEM=IDATA(1)
      ELSE
        IF(CALL_ELEM) THEN
          IDEFLT(1)=NEELEM(0,nr)
        ELSE
          IDEFLT(1)=1
        ENDIF
        WRITE(CHAR1,'(I5)') IDEFLT(1)
        FORMAT='($,'' The number of elements is ['//CHAR1(1:5)
     '    //']: '',I5)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NEM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C LKC 28-JUN-1998 new assert
C        IF(IOTYPE.NE.3) NTELEM=IDATA(1)
        IF(IOTYPE.NE.3) THEN
          NTELEM=IDATA(1)
          IF(TYPE(1:5).EQ.'FIELD') THEN
            CALL ASSERT(NTELEM.GE.NEELEM(0,nr),
     '        '>> Field elements <> geometry elements',ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF
      CALL ASSERT(NTELEM.GT.0,'>>No elements being defined',ERROR,*9999)
      CALL ASSERT(NTELEM.LE.NE_R_M,'>>Increase NE_R_M',ERROR,*9999)
      N2ELEM=1

      IF(IOTYPE.NE.3.AND.TYPE(1:8).EQ.'GEOMETRY'.AND..NOT.ADD) THEN
        NEELEM(0,nr)=0
      ENDIF

      DO noelem=1,NTELEM

        IF(TYPE(1:8).EQ.'GEOMETRY'.OR.NEELEM(0,nr).EQ.0) THEN
          IF(noelem.EQ.1) THEN
            IF(ADD) THEN
              IDEFLT(1)=NET(nr)+1
            ELSE
              IF(nr.EQ.1) THEN
                IDEFLT(1)=1
              ELSE !Choose 1+maximum element number to date
                IDEFLT(1)=NET(0)+1
              ENDIF
            ENDIF
          ELSE
            IDEFLT(1)=NEELEM(N2ELEM,nr)+1
          ENDIF
        ELSE
          IDEFLT(1)=NEELEM(noelem,nr)
        ENDIF
        WRITE(CHAR1,'(I5)') IDEFLT(1)
        IF(IOTYPE.EQ.3) THEN
C GDR 17/1/05 Make the "...;w elements #"option work  
C          ne=NEELEM(noelem,nr)
          ne=NELIST(noelem)
          IDATA(1)=ne
        ENDIF
        FORMAT='(/$,'' Element number ['//CHAR1(1:5)//']: '',I5)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NEM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ne=IDATA(1)+NE_OFFSET

        IF(TYPE(1:8).EQ.'GEOMETRY') THEN
          NRE(ne)=nr
          IF(noelem.EQ.1) THEN
            IDEFLT(1)=NJT
          ELSE
            IDEFLT(1)=NJ_LOC(NJL_GEOM,0,nr)
          ENDIF
          WRITE(CHAR2,'(I1)') IDEFLT(1)
          FORMAT='($,'' The number of geometric Xj-coordinates is ['//
     '      CHAR2(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NJ_LOC(NJL_GEOM,0,nr)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NJ_LOC(NJL_GEOM,0,nr)=IDATA(1)
          IF(JTYP6.EQ.2) THEN
            IDEFLT(1)=ITYP10(nr)
            WRITE(CHAR1,'(I1)') IDEFLT(1)
            FORMAT='($,'' The coordinate system is ['//CHAR1(1:1)//
     '        ']: '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=ITYP10(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ITYP10(nr)=IDATA(1)
          ENDIF
          njtype=NJL_GEOM
          njmax=NJ_LOC(NJL_GEOM,0,nr)
        ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
          njtype=NJL_FIBR
          njmax=NJ_LOC(NJL_FIBR,0,nr)
        ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
          njtype=NJL_FIEL
          njmax=NJ_LOC(NJL_FIEL,0,nr)
        ENDIF
        CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)

        DO njj=1,njmax
          nj=NJ_LOC(njtype,njj,nr)
          IF(TYPE(1:8).EQ.'GEOMETRY'.OR.TYPE(1:5).EQ.'FIBRE') THEN
            IF(noelem.EQ.1) THEN
              IDEFLT(1)=1
            ELSE IF(noelem.GT.1) THEN
              IDEFLT(1)=NBJ(nj,NEELEM(N2ELEM,nr))
            ENDIF
C DMAL/MPN 2003Feb17: fixing default field element basis functions
          ELSEIF(TYPE(1:5).EQ.'FIELD') THEN
            IDEFLT(1)=NBJ(NJ_LOC(NJL_GEOM,1,nr),ne)
          ELSE   ! 'FIELD' or other
            IF(CALL_ELEM) THEN
              IF(njj.LE.NJ_LOC(NJL_GEOM,0,nr)) THEN
C               Default to basis function of associated geom var in ne
                IDEFLT(1)=NBJ(NJ_LOC(NJL_GEOM,njj,nr),ne)
              ELSE
C               Default to basis function of first geom var in ne
                IDEFLT(1)=NBJ(NJ_LOC(NJL_GEOM,1,nr),ne)                
              ENDIF
C old              IDEFLT(1)=NBJ(njj,ne)
            ELSE
              IF(noelem.EQ.1) THEN
                IDEFLT(1)=1
              ELSE IF(noelem.GT.1) THEN
C               Default to basis function of previous element
                IDEFLT(1)=NBJ(nj,NEELEM(N2ELEM,nr))
              ENDIF
            ENDIF              
          ENDIF
          WRITE(CHAR1,'(I2)') njj
          WRITE(CHAR2,'(I3)') IDEFLT(1)
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          IF(TYPE(1:8).EQ.'GEOMETRY') THEN
            FORMAT='($,'' The basis function type for geometric'//
     '        ' variable '//CHAR1(IBEG1:IEND1)//' is ['//
     '        CHAR2(IBEG2:IEND2)//']: '',I2)'
          ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
            IF(njj.EQ.1) THEN
              FORMAT='($,'' The basis function type for the fibre '
     '          //'angle is ['//CHAR2(IBEG2:IEND2)//']: '',I2)'
            ELSE IF(njj.EQ.2) THEN
              FORMAT='($,'' The basis function type for the '
     '          //'imbrication angle is ['//CHAR2(IBEG2:IEND2)
     '          //']: '',I2)'
            ELSE IF(njj.EQ.3) THEN
              FORMAT='($,'' The basis function type for the sheet '
     '          //'angle is ['//CHAR2(IBEG2:IEND2)//']: '',I2)'
            ENDIF
          ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
            FORMAT='($,'' The basis function type for field'//
     '        ' variable '//CHAR1(IBEG1:IEND1)//' is ['//
     '        CHAR2(IBEG2:IEND2)//']: '',I2)'
          ENDIF
          IF(IOTYPE.EQ.3) IDATA(1)=NBJ(nj,ne)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(BASIS_OVERRIDE)THEN
              nb=nbasis
              nb_actual=IDATA(1)
            ELSE
              nb=IDATA(1)
            ENDIF
            NBJ(nj,ne)=nb
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,ne)=nk
              ENDDO !nk
            ENDDO !nn
          ENDIF ! IOTYPE.NE.3
        ENDDO !njj

        FIRST=.TRUE.
        FIRST_VERSION=.TRUE.
        DO nb=1,NBFT
          IF(IOTYPE.NE.3) THEN
C           initialise NVJE with default value
C            DO njj1=1,3 !geometry/fibres/field
            DO njj2=1,NJ_LOC(njtype,0,nr)
              nj=NJ_LOC(njtype,njj2,nr)
              DO nn=1,NNT(nb)
                IF(TREE)THEN
                  !default version for node in N'th element around node is N
                  np=NPNE(nn,nb,ne)
                  noelem2=1
                  ne2=NENP(np,noelem2,nr)
                  DO WHILE(ne2.NE.ne)
                    noelem2=noelem2+1
                    ne2=NENP(np,noelem2,nr)
                  ENDDO
                  NVJE(nn,nb,nj,ne)=noelem2
                ELSE
                  NVJE(nn,nb,nj,ne)=1
                ENDIF
              ENDDO !nn
            ENDDO !njj2
C            ENDDO !njj1
          ENDIF
          SAMETYPE=.FALSE.
          SAMENUMXI=.FALSE.
          NBLOG=.FALSE.
          IF(TYPE(1:8).EQ.'GEOMETRY') THEN
C new MPN 29Apr97:
            DO njj2=1,NJ_LOC(njtype,0,nr)
              nj=NJ_LOC(njtype,njj2,nr)
C old            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              IF(nb.EQ.NBJ(nj,ne)) NBLOG=.TRUE.
              IF(NBC(nb).EQ.NBC(NBJ(nj,ne)).OR.NBC(nb).EQ.7)
     '          SAMETYPE=.TRUE. !Same basis type or extended basis
              IF(NIT(nb).EQ.NIT(NBJ(nj,ne))) SAMENUMXI=.TRUE.
            ENDDO !nj
            DOBASIS=NNT(nb).GT.0.AND.SAMETYPE.AND.(NBLOG.OR.
     '        SAMENUMXI.OR.IBT(1,NIT(nb),nb).EQ.9)
            FIBREFIELD=.FALSE.
          ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
            DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
              nj=NJ_LOC(NJL_FIBR,njj,nr)
              IF(nb.EQ.NBJ(nj,ne)) NBLOG=.TRUE.
            ENDDO !njj
            DOBASIS=NNT(nb).GT.0.AND.NBLOG
            FIBREFIELD=.TRUE.
          ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
            DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
              nj=NJ_LOC(NJL_FIEL,njj,nr)
              IF(nb.EQ.NBJ(nj,ne)) NBLOG=.TRUE.
            ENDDO !njj
            DOBASIS=NNT(nb).GT.0.AND.NBLOG
            FIBREFIELD=.TRUE.
          ENDIF
          IF(DOBASIS) THEN
            WRITE(CHAR1,'(I2)') NNT(nb)
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            IF(BASIS_OVERRIDE)THEN
              WRITE(CHAR2,'(I3)') nb_actual
            ELSE
              WRITE(CHAR2,'(I3)') nb
            ENDIF
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            IF(FIBREFIELD) THEN
              FOUND=.FALSE.
              nb2=1
C cpb 30/4/97 Changing back
CC new MPN 29Apr97:
C              DO WHILE(.NOT.FOUND.AND.nb2.LT.nb)
CC old              DO WHILE(.NOT.FOUND.AND.nb2.LE.NBT)
              DO WHILE(.NOT.FOUND.AND.nb2.LE.NBT)
                IF(NBC(nb2).EQ.NBC(nb)) THEN
C new MPN 29Apr97:
                  DO njj2=1,NJ_LOC(njtype,0,nr)
                    nj=NJ_LOC(njtype,njj2,nr)
C old                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    IF(NBJ(nj,ne).EQ.nb2) FOUND=.TRUE.
                  ENDDO !nj
                  IF(.NOT.FOUND) nb2=nb2+1
                ELSE
                  nb2=nb2+1
                ENDIF
              ENDDO !nb2
              IF(FOUND) THEN
                DO nn=1,NNT(nb)
                  IDEFLT(nn)=NPNE(nn,nb2,ne)
                ENDDO !nn
                FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//
     '            ' numbers for basis '//CHAR2(IBEG2:IEND2)//
     '            ' [prev]: '',I5,14(1X,I5))'
              ELSE
                DO nn=1,NNT(nb)
                  IDEFLT(nn)=1
                ENDDO !nn
                FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//
     '            ' global numbers for basis '//CHAR2(IBEG2:IEND2)//
     '            ': '',I5,14(1X,I5))'
              ENDIF
            ELSE
              IF(FIRST) THEN
C               note: FIRST is set to .FALSE. below
                DO nn=1,NNT(nb)
                  IDEFLT(nn)=1
                ENDDO !nn
                FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//
     '            ' global numbers for basis '//CHAR2(IBEG2:IEND2)//
     '            ': '',I5,14(1X,I5))'
              ELSE
                nb2=nb-1
C If nb2=0 this condition set may crash
C                DO WHILE(nb2.GT.0.AND.NNT(nb2).NE.NNT(nb))
C                  nb2=nb2-1
C                ENDDO !nb2
                CONTINU=.TRUE.
                DO WHILE(CONTINU)
                  IF(nb2.GT.0) THEN
                    IF(NNT(nb2).NE.NNT(nb)) THEN
                      nb2=nb2-1
                    ELSE
                      CONTINU=.FALSE.
                    ENDIF
                  ELSE
                    CONTINU=.FALSE.
                  ENDIF
                ENDDO !continu
                IF(nb2.EQ.0) nb2=nb-1
                DO nn=1,NNT(nb)
                  IDEFLT(nn)=NPNE(nn,nb2,ne)
                ENDDO !nn
                FORMAT='($,'' Enter the '//CHAR1(IBEG:IEND)//
     '            ' numbers for basis '//CHAR2(IBEG2:IEND2)//
     '            ' [prev]: '',I5,14(1X,I5))'
              ENDIF
            ENDIF

 10         IF(IOTYPE.EQ.3) THEN
              DO nn=1,NNT(nb)
                IDATA(nn)=NPNE(nn,nb,ne)
              ENDDO !nn
            ENDIF
            CDATA(1)='NODES' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NNT(nb),
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
C news MPN 17-Jul-95: checking for existence of element nodes in global
C                     node list for current region
              DO nn=1,NNT(nb)
                np=IDATA(nn)+NP_OFFSET
                NP_EXISTS=.FALSE.
                nonode=0
                DO WHILE(nonode.LT.NPNODE(0,nr).AND..NOT.NP_EXISTS)
                  nonode=nonode+1
                  IF(np.EQ.NPNODE(nonode,nr)) NP_EXISTS=.TRUE.
                ENDDO !nonode

                IF(.NOT.NP_EXISTS) THEN
C                 node np is not defined in current region
                  WRITE(CHAR5,'(I5)') np
                  WRITE(OP_STRING,'(''>>Node '//CHAR5(1:5)
     '              //' is not defined - reenter'')')
                  CALL WRITES(IOER,OP_STRING,ERROR,*9999)
                  GO TO 10
                ENDIF
              ENDDO !nn
C newe

C news MPN 13-Dec-94
              DO nn=1,NNT(nb)
                NPNE(nn,nb,ne)=IDATA(nn)+NP_OFFSET
C KAT 23Feb01: now handled by NKJE above
C                DO nk=1,NKT(nn,nb)
C                  NKE(nk,nn,nb,ne)=nk
C                ENDDO !nk
              ENDDO !nn
c cpb 27/1/97 Redundant now
CC             if fibres/imbric/sheets already set up the store basis
CC             function number for ne in NBJ
C              IF(NJ_LOC(NJL_FIBR,1,nr).GT.0) THEN  !fibres defined
C                NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)=NB_FIBRE
C              ENDIF
C              IF(NJ_LOC(NJL_FIBR,2,nr).GT.0) THEN  !imbrics defined
C                NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)=NB_IMBRIC
C              ENDIF
C              IF(NJ_LOC(NJL_FIBR,3,nr).GT.0) THEN  !sheets defined
C                NBJ(NJ_LOC(NJL_FIBR,3,nr),ne)=NB_SHEET
C              ENDIF
            ENDIF
C           Prompt for different versions for repeated element nodes
C           Check for repeated nodes for the current element/basis
            DO nn=1,NNT(nb)
              np=NPNE(nn,nb,ne)-NP_OFFSET
C             prompt for version numbers if an nj of the current
C             element vertex has more than one version
C cpb 27/1/97 adding fibre/field elements
C              DO njj1=1,3  !geometry/fibres/field
C                DO njj2=1,NJ_LOC(njj1,0,nr)
              DO njj2=1,NJ_LOC(njtype,0,nr)
                nj=NJ_LOC(njtype,njj2,nr)
                IF(nb.EQ.NBJ(nj,ne).OR.NBC(nb).EQ.7) THEN
C                 only prompt for versions of current nj if nb is
C                 the basis fn number used to interpolate nj
C                 -- or if collocation basis (GBS 12/6/96)
                  IF(NVJP(nj,np+NP_OFFSET).GT.1) THEN !>1 global version
                    OC_COUNT=1
C                   count the number of previous occurances of the
C                   current element vertex in the element list
                    DO mm=1,nn-1
                      mp=NPNE(mm,nb,ne)-NP_OFFSET
                      IF(np.EQ.mp) OC_COUNT=OC_COUNT+1
                    ENDDO !mm
C                   prompt for the version number of nj for the
C                   OC_COUNT'th occurance of np in the element list
                    WRITE(CHAR1,'(I2)') OC_COUNT
                    WRITE(CHAR5,'(I5)') np
                    WRITE(CHAR6,'(I6)') np
                    IF(FIBREFIELD) THEN
                      IF(njj2.LT.10) THEN
                       WRITE(CHAR3,'(I1)') njj2
                      ELSE
                       WRITE(CHAR3,'(I2)') njj2
                      ENDIF
                    ELSE
                      WRITE(CHAR3,'(I1)') nj
                    ENDIF
                    IF(FIRST_VERSION) THEN
                      IF(.NOT.TREE)THEN
                        IDEFLT(1)=1
                      ELSE
                        IDEFLT(1)=NVJE(nn,nb,nj,ne)
                      ENDIF
                    ELSE IF(.NOT.FIRST_VERSION) THEN
                      nb2=nb-1
C                     look for previous basis with prompts for this
C                     version
C poor coding,replaced DO WHILE(nb2.GT.0.AND.
C     '                  (nb2.NE.NBJ(nj,ne).AND.NBC(nb2).NE.7))
C                        nb2=nb2-1
C                      ENDDO !nb2
                      FIN=.FALSE.
                      DO WHILE(.NOT.FIN)
                        IF(nb2.LE.0) THEN
                          FIN=.TRUE.
                        ELSE
                          IF(NBC(nb2).EQ.7) THEN
                            FIN=.TRUE.
                          ELSE
                            IF(nb2.EQ.NBJ(nj,ne)) THEN
                              FIN=.TRUE.
                            ELSE
                              nb2=nb2-1
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO
                      IF(nb2.EQ.0) THEN
                        IF(TREE)THEN
                          IDEFLT(1)=NVJE(nn,nb,nj,ne)
                        ELSE
C MPN 21Jan2005:        if nj is a field variable and the current np has the
C                       same number of versions as the corresponding
C                       geometric variable, then almost certainly going
C                       to do a  geometric fitting problem, so best to use the
C                       version number from the geometric variable as
C                       the most sensible default
C                       (this fixes fem def elem;d field)
                          IF(njtype.EQ.NJL_FIEL) THEN !field variable
                            nj2=NJ_LOC(NJL_GEOM,njj2,nr) !equivalent geom nj
                            IF(nj2.NE.0)THEN
                              IF(NVJP(nj,np).EQ.NVJP(nj2,np)) THEN
                                IDEFLT(1)=NVJE(nn,nb,nj2,ne)
                              ENDIF !NVJP(nj,np).EQ.NVJP(nj2,np)
                            ELSE
                              IDEFLT(1)=1
                            ENDIF !nj2.NE.0
                          ELSE
                            IDEFLT(1)=1
c                              IDEFLT(1)=NVJE(nn,nb2,nj,ne)
C                              IDEFLT(1)=NVJE(mm,nb2,nj,ne)
                          ENDIF !njtype.EQ.NJL_FIEL
                        ENDIF

C MPN 15-Apr97: bug fix
                      ENDIF
                    ENDIF !FIRST_VERSION
                    IF(TYPE(1:5).EQ."FIBRE") THEN
                      IDEFLT(1)=NVJE(nn,nb,1,ne)
                    ENDIF
                    WRITE(CHAR4,'(I2)') IDEFLT(1)
C MPN 15-Apr97: new
                    IF(njj2.LT.10) THEN
                     FORMAT='($,'' The version number for'
     '                //' occurrence '//CHAR1(1:2)
     '                //' of node '//CHAR5(1:5)
     '                //', njj='//CHAR3(1:1)//' is ['
     '                //CHAR4(1:2)//']: '',I2)'
                    ELSE
                     FORMAT='($,'' The version number for'
     '                //' occurrence '//CHAR1(1:2)
     '                //' of node '//CHAR5(1:5)
     '                //', njj='//CHAR3(1:2)//' is ['
     '                //CHAR4(1:2)//']: '',I2)'
                    ENDIF

C MPN 15-Apr97: redundant
C                    IF(FIBREFIELD) THEN
C                      FORMAT='($,'' The version number for'
C     '                  //' occurrence '//CHAR1(1:2)
C     '                  //' of node '//CHAR5(1:5)
C     '                  //', njj='//CHAR3(1:1)//' is ['
C     '                  //CHAR4(1:2)//']: '',I2)'
C                    ELSE
C                      FORMAT='($,'' The version number for'
C     '                  //' occurrence '//CHAR1(1:2)
C     '                  //' of node '//CHAR5(1:5)
C     '                  //', njj='//CHAR3(1:1)//' is ['
C     '                  //CHAR4(1:2)//']: '',I2)'
C                    ENDIF
C MPN 15-Apr97: bug fix
                    IF(IOTYPE.EQ.3) IDATA(1)=NVJE(nn,nb,nj,ne)
C                    IF(IOTYPE.EQ.3) IDATA(1)=NVJE(mm,nb,nj,ne)
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,
     '                NOQUES,FILEIP,FORMAT,1,
     '                ADATA,ADEFLT,CDATA,CDEFLT,
     '                ICHAR,IDATA,IDEFLT,1,NVJP(nj,np+NP_OFFSET),
     '                LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '                INFO,ERROR,*9999)
C MPN 15-Apr97: bug fix
                    IF(IOTYPE.NE.3) NVJE(nn,nb,nj,ne)=IDATA(1)
C                    IF(IOTYPE.NE.3) NVJE(mm,nb,nj,ne)=IDATA(1)
                    IF(FIRST_VERSION) FIRST_VERSION=.FALSE.
                  ENDIF !NVJP(nj,np).GT.1
                ELSE
C MPN 15-Apr97: new
C                 Initialise NVJE to 1
                  NVJE(nn,nb,nj,ne)=1
                ENDIF !nb.EQ.NBJ(nj,ne).OR.NBC(nb).EQ.7
              ENDDO !njj2
C              ENDDO !njj1
            ENDDO !nn

            IF(FIRST) FIRST=.FALSE.
C newe

          ENDIF
        ENDDO !nb

        IF(IOTYPE.NE.3) THEN
          IF(TYPE(1:8).EQ.'GEOMETRY') THEN
            EXISTS=.FALSE.
            DO n1elem=1,NEELEM(0,nr)
              IF(NEELEM(n1elem,nr).EQ.ne) THEN
C             ne exists in current region
                EXISTS=.TRUE.
                N2ELEM=n1elem
              ENDIF
            ENDDO !n1elem
            IF(.NOT.EXISTS) THEN !new element NE
              NEELEM(0,nr)=NEELEM(0,nr)+1
              CALL ASSERT(NEELEM(0,nr).LE.NE_R_M,'>>Increase NE_R_M',
     '          ERROR,*9999)
              NEELEM(NEELEM(0,nr),nr)=ne
              N2ELEM=NEELEM(0,nr)
            ENDIF

C news MPN 22-Jun-94
            NE_EXISTS=.FALSE.
            DO nr2=1,NRT
              IF((nr2.NE.nr).AND.(.NOT.NE_EXISTS))THEN
                IF(NEELEM(0,nr2).GT.0)THEN !Some elements defined in nr2
                  DO noelem2=1,NEELEM(0,nr2)
                    IF(ne.EQ.NEELEM(noelem2,nr2)) NE_EXISTS=.TRUE.
                  ENDDO !noelem2
                ENDIF
              ENDIF
            ENDDO !nr2
            IF(.NOT.(EXISTS.OR.NE_EXISTS))THEN
C             NE does not exists in any region so increment global #
              NEELEM(0,0)=NEELEM(0,0)+1
            ENDIF
C newe
          ELSE
            EXISTS=.FALSE.
            IF(CALL_ELEM) THEN
              DO n1elem=1,NEELEM(0,nr)
                IF(NEELEM(n1elem,nr).EQ.ne) THEN
C               ne exists in current region
                  EXISTS=.TRUE.
                  N2ELEM=n1elem
                ENDIF
              ENDDO !n1elem
            ENDIF
            IF(.NOT.EXISTS) THEN !new element NE
              N2ELEM=NEELEM(0,nr)+1
            ENDIF
          ENDIF
        ELSE
C         Default is last node for backward compatibility
C         (but don't need a default at all).
          N2ELEM=noelem
        ENDIF !IOTYPE eq/ne 3
      ENDDO !noelem

      IF(TYPE(1:8).EQ.'GEOMETRY') THEN
        IF(IOTYPE.NE.3) THEN
          CALL ISORT(NEELEM(0,nr),NEELEM(1,nr))
          NET(nr)=0
          DO noelem=1,NEELEM(0,nr)
            IF(NEELEM(noelem,nr).GT.NET(nr))
     '        NET(nr)=NEELEM(noelem,nr)
          ENDDO
C Code below added 22-10-92 AJP & 31Mar93 PJH
          IF(NET(nr).GT.NET(0)) NET(0)=NET(nr)
C old MPN 22-Jun-94: moved up into noelem loop
c        DO noelem=1,NEELEM(0,nr) !Set up NEELEM(np,0)
c          ne=NEELEM(noelem,nr)
c          EXISTS=.FALSE.
c          DO nr2=1,NRT
c            IF((nr2.NE.nr).AND.(.NOT.EXISTS))THEN
c              IF(NEELEM(0,nr2).GT.0)THEN !Some elements defined in nr2
c                DO NOELEM2=1,NEELEM(0,nr2)
c                  IF(ne.EQ.NEELEM(NOELEM2,nr2))EXISTS=.TRUE.
c                ENDDO
c              ENDIF
c            ENDIF
c          ENDDO
c          IF(.NOT.EXISTS)THEN
c            NEELEM(0,0)=NEELEM(0,0)+1
c          ENDIF
c        ENDDO
        ENDIF
      ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
C  Setup version info for extended basis
        nb_extended=1
        DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
          nb_extended=nb_extended+1
        ENDDO
        IF(NBC(nb_extended).EQ.7) THEN
          WRITE(OP_STRING,
     '      '(/'' >>Version numbers updated for extended basis.'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(IOTYPE.NE.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NB_FIBRE=NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
            DO nn=1,NNT(NB_FIBRE)
              IF(NBC(nb_extended).EQ.7) THEN
                NVJE(nn,nb_extended,NJ_LOC(NJL_FIBR,1,nr),ne)=1
              ENDIF
            ENDDO !nn
            IF(NJ_LOC(NJL_FIBR,0,nr).GE.2) THEN
              NB_IMBRIC=NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)
              DO nn=1,NNT(NB_IMBRIC)
                IF(NBC(nb_extended).EQ.7) THEN
                  NVJE(nn,nb_extended,NJ_LOC(NJL_FIBR,2,nr),ne)=1
                ENDIF
              ENDDO !nn
            ENDIF
          ENDDO !noelem
        ENDIF
      ENDIF

      CALL EXITS('IPELEM')
      RETURN
 9999 CALL ERRORS('IPELEM',ERROR)
      CALL EXITS('IPELEM')
      RETURN 1
      END


