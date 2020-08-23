      SUBROUTINE IPMESH2_1(NBJ,NBJ_TEMP,NEELEM,NEELEM_TEMP,NELIST2,NENP,
     &  NKJ,NKJ_TEMP,NKJE,NKJE_TEMP,NORD,NORD_TEMP,NP_ATTACH,
     &  NP_INTERFACE,NPLIST,NPNE,NPNE_TEMP,NPNODE,NPNODE_TEMP,
     &  nr,NRE,NRE_TEMP,NVJE,NVJE_TEMP,NVJP,NVJP_TEMP,NXI,SE,XP,XP_TEMP,
     &  ERROR,*)

C#### Subroutine: IPMESH2_1
C###  Description:
C###    IPMESH2_1 in/out-puts global coordinates, fields, elements for
C###    lung tree meshes.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'parameters.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBJ_TEMP(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NEELEM_TEMP(NE_R_M),NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     &  NKJ(NJM,NPM),NKJ_TEMP(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NKJE_TEMP(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),NORD_TEMP(NE_R_M),
     &  NP_ATTACH,NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNE_TEMP(NNM,NBFM,NEM),NPLIST(0:NPM),
     &  NPNODE(0:NP_R_M,0:NRM),NPNODE_TEMP(NP_R_M),nr,NRE(NEM),
     &  NRE_TEMP(NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJE_TEMP(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NVJP_TEMP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),
     &  XP_TEMP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,nb,N2ELEM,ne,ne2,ne_offset,
     &  ne_parent,nj,nj2,nj_geometry,njj,njj1,njj2,NJJ_END,
     &  NJJ_START,NJTOT,njtype,nk,NKJP(NJ_LOC_MX),nn,noelem,
     &  noelem2,noelem_parent,nonode,nonode2,NOQUES,np,
     &  np_current,np_default,np_offset,np_parent,np1,NP2,
     &  ns,NTELEM,NTNODE,ntype,nu,NUM_FIELD,nv,NV_MAX,SYMM_INDEX
      REAL*8 length,XP_OFFSET(3)
      CHARACTER CHAR1*5,CHAR2*1,CHAR4*12,CHAR5*2
      LOGICAL FILEIP,FIRST_VERSION,PROMPT_NV(NJ_LOC_MX),SYMMETRY_PROMPT

      CALL ENTERS('IPMESH2_1',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      NJJ_START=1 !for mesh definition the field starts anew

      DO np=1,NPM
        NPLIST(np)=0
      ENDDO
      NPLIST(0)=0

      IF(IOTYPE.NE.3) THEN
        IF(ADD)THEN
          IF(NELIST2(0).EQ.0)THEN !list of parent elements has not been defined
C***    Get the list of parent elements        
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(NXI(1,0,ne).EQ.0)THEN !APPEND MODEL HERE
                NELIST2(0)=NELIST2(0)+1
                NELIST2(NELIST2(0))=ne
              ENDIF
            ENDDO !noelem
          ENDIF        
          np_current=NPT(0) !current highest global node #
        ELSE
          NELIST2(0)=1
          NELIST2(1)=0
        
          np_current=0
          NPT(nr)=0
          IF(NEELEM(0,nr).GT.0)THEN
            WRITE(OP_STRING,'('' Overwriting existing elements''
     &        /''Use FEM DEFINE;ADD MESH to append to existing mesh'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !ADD
        np_offset=0
        ne_offset=0
      ELSE
        np_current=0
        np_offset=0
        ne_offset=0
      ENDIF!IOTYPE

      DO nj=1,NJT
        XP_OFFSET(nj)=0.d0
      ENDDO !nj
      
      IDEFLT(1)=1
      WRITE(CHAR1,'(I5)') IDEFLT(1)
      IF(IOTYPE.EQ.3) THEN
        NTNODE=NPNODE(0,nr)
        IDATA(1)=NTNODE
      ENDIF
      FORMAT='($,'' The number of nodes is ['//CHAR1(1:5)//']: '',I5)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,LDATA,LDEFLT,
     &  RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NTNODE=IDATA(1)
      CALL ASSERT(NTNODE.LE.NP_R_M,'>>Increase NP_R_M',ERROR,
     &  *9999)

      IDEFLT(1)=NJT
      WRITE(CHAR1,'(I1)') IDEFLT(1)
      FORMAT='($,'' Number of coordinates ['//CHAR1(1:1)//']: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NJT
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,LDATA,LDEFLT,
     &  RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADD)THEN
          CALL ASSERT(IDATA(1).EQ.NJT,'>>Inconsistent coordinates',
     '      ERROR,*9999)
        ELSE
          CALL ASSERT(IDATA(1).EQ.NJT,'>>Define coordinates first',
     '      ERROR,*9999)
        ENDIF
      ENDIF

C*** Setting NJ_LOC(NJL_GEOM) here rather than at elements      
      IF(IOTYPE.NE.3) NJ_LOC(NJL_GEOM,0,nr)=NJT !used assert, so ok for add

C*** Get field information
      IDEFLT(1)=1 !default to a single field
      WRITE(CHAR2,'(I1)') IDEFLT(1)
      IF(IOTYPE.EQ.3) THEN
        IDATA(1)=NJ_LOC(NJL_FIEL,0,nr)
      ENDIF
      FORMAT='($,'' Number of field variables ['//CHAR2(1:1)//
     &  ']: '',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=NJ_LOC(NJL_FIEL,0,nr)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,
     &  RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NUM_FIELD=IDATA(1)
        IF(ADD)THEN
          CALL ASSERT(IDATA(1).EQ.NJ_LOC(NJL_FIEL,0,nr),
     &      '>>Inconsistent number of fields',ERROR,*9999)
        ELSE
          CALL IPFIEL_NJLOC(NJJ_START,NKJ,NPNODE,nr,NUM_FIELD,NVJP,
     &      ERROR,*9999)
          CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).LE.NJ_LOC_MX,
     '      '>>Increase NJ_LOC_MX',ERROR,*9999)
        ENDIF
      ENDIF

      IF(NTNODE.GT.0) THEN
C***  Prompt for versions of coordinates and field variables
C***    Versions for coordinates
        DO ntype=1,2
          IF(ntype.EQ.1)THEN
            NJJ_END=NJT
          ELSE
            NJJ_END=NJ_LOC(NJL_FIEL,0,nr)
          ENDIF
          DO njj=1,NJJ_END
            WRITE(CHAR5,'(I2)') njj
            IF(ntype.EQ.1)THEN !geometry
              nj=njj
              FORMAT='($,'' Do you want prompting for different '//
     &          'versions of nj='//CHAR5(1:2)//' [N]? '',A)'
            ELSE !field variables
              nj=NJ_LOC(NJL_FIEL,njj,nr)
              FORMAT='($,'' Do you want prompting for different '//
     &          'versions of field variable'//CHAR5(1:2)//' [N]? '',A)'
            ENDIF
            IF(IOTYPE.EQ.3) THEN
              PROMPT_NV(nj)=.FALSE.
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(NVJP(nj,np).GT.1) PROMPT_NV(nj)=.TRUE.
              ENDDO
              IF(PROMPT_NV(nj)) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                PROMPT_NV(nj)=.TRUE.
              ELSE
                PROMPT_NV(nj)=.FALSE.
              ENDIF
            ENDIF
          ENDDO !njj
        ENDDO !ntype

C***  Prompt for derivatives of coordinates and field variables
C***  (There is a maximum of 1 derivative for this mesh type)
        DO ntype=1,2
          IF(ntype.EQ.1)THEN
            NJJ_END=NJT
          ELSE
            NJJ_END=NJ_LOC(NJL_FIEL,0,nr)
          ENDIF
          DO njj=1,NJJ_END
            WRITE(CHAR5,'(I2)') njj
            IF(ntype.EQ.1)THEN !geometry
              nj=njj
              FORMAT='($,'' The number of derivatives for coordinate '
     '          //CHAR5(1:2)//' is [0]: '',I1)'
            ELSE !field variables
              nj=NJ_LOC(NJL_FIEL,njj,nr)
              FORMAT=
     &          '($,'' The number of derivatives for field '//
     &          'variable '//CHAR5(1:2)//' is [0]: '',I2)'
            ENDIF
            IF(IOTYPE.EQ.3) THEN !get max # of derivs
              NKJP(nj)=0
              DO nonode=1,NTNODE
                IF(NKJ(nj,NPNODE(nonode,nr)).GT.NKJP(nj)) THEN
                  NKJP(nj)=NKJ(nj,NPNODE(nonode,nr))
                ENDIF
              ENDDO
              IDATA(1)=NKJP(nj)-1
            ENDIF
            IDEFLT(1)=0
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,1,
     &        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NKJP(nj)=IDATA(1)+1
          ENDDO !njj
        ENDDO !ntype

C        N2NODE=0
        IF(IOTYPE.NE.3.AND..NOT.ADD) THEN
          NPNODE(0,nr)=0
        ENDIF

        np_default=np_current
        IF(ADD)THEN
          ne_parent=NELIST2(1)
          np_parent=NPNE(2,NBJ(1,ne_parent),ne_parent)
        ENDIF

        DO nonode=1,NTNODE
          IDEFLT(1)=np_default+1 !previous np + 1
          WRITE(CHAR1,'(I5)') IDEFLT(1)
          IF(IOTYPE.EQ.3) THEN
            np=NPNODE(nonode,nr)
            IDATA(1)=np
          ENDIF
          FORMAT='(/$,'' Node number ['//CHAR1(1:5)//']: '',I5)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &      NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) np=IDATA(1)

          NJTOT=NJT
          DO nj=1,NJTOT
            IF(PROMPT_NV(nj)) THEN !prompt for different versions of nj
              WRITE(CHAR1,'(I1)') nj
              FORMAT='($,'' The number of versions for nj='
     '          //CHAR1(1:1)//' is [1]: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=NVJP(nj,np)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &          NVM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                CALL ASSERT(IDATA(1).LE.NVM,'>>Increase NVM',
     '            ERROR,*9999)
                IF(ADD)THEN
                  NVJP_TEMP(nj,np)=IDATA(1)
                ELSE
                  NVJP(nj,np)=IDATA(1)
                ENDIF
              ENDIF
            ELSE
              IF(ADD)THEN
                NVJP_TEMP(nj,np)=1
              ELSE
                NVJP(nj,np)=1
              ENDIF
            ENDIF
            
            IF(ADD)THEN
              NV_MAX=NVJP_TEMP(nj,np)
            ELSE
              NV_MAX=NVJP(nj,np)
            ENDIF
            
            DO nv=1,NV_MAX
              IF(NV_MAX.GT.1) THEN !ask for diff nj versions
                WRITE(CHAR1,'(I2)') nv
                FORMAT='('' For version number'//CHAR1(1:2)//':'')'
                CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,
     &            FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,
     &            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     &            INFO,ERROR,*9999)
              ENDIF
              WRITE(CHAR1,'(I1)') nj
              IF(IOTYPE.NE.3)THEN
                IF(ADD)THEN
                  NKJ_TEMP(nj,np)=NKJP(nj)
                ELSE
                  NKJ(nj,np)=NKJP(nj)
                ENDIF
              ENDIF

              DO nk=1,NKJP(nj)
                nu=nk !note: in IPNODE would be nu=NUK(nk)
                IF(IOTYPE.NE.3)THEN
                  RDEFLT(1)=0.d0
                ELSE
                  RDEFLT(1)=XP(nk,nv,nj,np)
                ENDIF
                WRITE(CHAR4,'(E12.5)') RDEFLT(1)
                IF(nu.EQ.1) THEN
                  FORMAT='($,'' The Xj('//CHAR1(1:1)//') coordinate'//
     '              ' is ['//CHAR4(1:12)//']: '',G25.17)'
                ELSE IF(nu.EQ.2) THEN !for 1D elements only
                  CHAR2='1'
                  FORMAT='($,'' The derivative wrt direction '//
     '              CHAR2(1:1)//' is ['//CHAR4(1:12)//']: '',G25.17)'
                ENDIF
                IF(IOTYPE.NE.3) THEN
                  IF(.NOT.ADD.AND.nk.EQ.1.AND.nv.EQ.1.AND.nj.EQ.1)THEN
                    NPNODE(0,0)=NPNODE(0,0)+1
                  ENDIF
                ELSE !IF(IOTYPE.EQ.3) THEN
                  RDATA(1)=XP(nk,nv,nj,np)
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,
     &            FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '            ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '            RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  IF(ADD)THEN
                    XP_TEMP(nk,nv,nj,np)=RDATA(1) !+XP_OFFSET(nj)
                  ELSE
c                  XP(nk,nv,nj,np+np_offset)=RDATA(1) !+XP_OFFSET(nj)
                    XP(nk,nv,nj,np)=RDATA(1) !+XP_OFFSET(nj)
                  ENDIF
                ENDIF
              ENDDO !nk
            ENDDO !nv
          ENDDO !nj
          
          DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            WRITE(CHAR5,'(I2)') njj
            IF(ADD)THEN
              NKJ_TEMP(nj,np)=NKJP(nj)
            ELSE
              NKJ(nj,np)=NKJP(nj)
            ENDIF
            IF(PROMPT_NV(nj)) THEN !prompt for diff versions of field
              FORMAT='($,'' The number of versions for'
     '          //' field variable'//CHAR5(1:2)//' is [1]: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=NVJP(nj,np)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &          NVM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3)THEN
                IF(ADD)THEN
                  NVJP_TEMP(nj,np)=IDATA(1)
                ELSE
c                    NVJP(nj,np+np_offset)=IDATA(1)
                  NVJP(nj,np)=IDATA(1)
                ENDIF
              ENDIF
            ELSE
              IF(ADD)THEN
                NVJP_TEMP(nj,np)=1
              ELSE
c              NVJP(nj,np+np_offset)=1
                NVJP(nj,np)=1
              ENDIF
            ENDIF
            IF(ADD)THEN
              DO nv=1,NVJP_TEMP(nj,np)
                IF(NVJP_TEMP(nj,np).GT.1) THEN !ask for diff nj versions
                  WRITE(CHAR5,'(I2)') nv
                  FORMAT='('' For version number'//CHAR5(1:2)//':'')'
                  CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              RMIN,RMAX,INFO,ERROR,*9999)
                ENDIF
                WRITE(CHAR5,'(I2)') njj
                DO nk=1,NKJP(nj)
                  nu=nk !note: in IPFIEL nu=NUK(nk), but here limit is 1 deriv
                  IF(np.EQ.1) THEN
                    RDEFLT(1)=0.0d0
                  ELSE IF(np.GT.1) THEN
                    CALL ASSERT(NKM.GE.nk,'>>Increase NKM',ERROR,
     &                *9999)
                    RDEFLT(1)=XP_TEMP(nk,nv,nj,np-1)
                  ENDIF
                  WRITE(CHAR4,'(D12.5)') RDEFLT(1)
                  IF(nu.EQ.1) THEN
                    FORMAT='($,'' The field variable'//CHAR5(1:2)
     '                //' value is ['//CHAR4//']: '',D12.5)'
                  ELSE IF(nu.EQ.2) THEN
                    CHAR1='1'
                    FORMAT='($,'' The derivative wrt direction '//
     '                CHAR1//' is ['//CHAR4//']: '',D12.5)'
                  ENDIF
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              -RMAX,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    XP_TEMP(nk,nv,nj,np)=RDATA(1)
                  ENDIF
                ENDDO !nk
              ENDDO !nv
            ELSE !not add
              DO nv=1,NVJP(nj,np)
                IF(NVJP(nj,np).GT.1) THEN !ask for diff nj versions
                  WRITE(CHAR5,'(I2)') nv
                  FORMAT='('' For version number'//CHAR5(1:2)//':'')'
                  CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              RMIN,RMAX,INFO,ERROR,*9999)
                ENDIF
                WRITE(CHAR5,'(I2)') njj
                DO nk=1,NKJP(nj)
                  nu=nk !note: in IPFIEL nu=NUK(nk), but here limit is 1 deriv
c                      IF(np+np_offset.EQ.1) THEN
                  IF(np.EQ.1) THEN
                    RDEFLT(1)=0.0d0
c                      ELSE IF(np+np_offset.GT.1) THEN
                  ELSE IF(np.GT.1) THEN
                    CALL ASSERT(NKM.GE.nk,'>>Increase NKM',ERROR,
     &                *9999)
c                      RDEFLT(1)=XP(nk,nv,nj,np+np_offset-1)
                    RDEFLT(1)=XP(nk,nv,nj,np-1)
                  ENDIF
                  WRITE(CHAR4,'(D12.5)') RDEFLT(1)
                  IF(nu.EQ.1) THEN
                    FORMAT='($,'' The field variable'//CHAR5(1:2)
     '                //' value is ['//CHAR4//']: '',D12.5)'
                  ELSE IF(nu.EQ.2) THEN
                    CHAR1='1'
                    FORMAT='($,'' The derivative wrt direction '//
     '                CHAR1//' is ['//CHAR4//']: '',D12.5)'
                  ENDIF
                  IF(IOTYPE.EQ.3) RDATA(1)=XP(nk,nv,nj,np)
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              -RMAX,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
c                        XP(nk,nv,nj,np+np_offset)=RDATA(1)
                    XP(nk,nv,nj,np)=RDATA(1)
                  ENDIF
                ENDDO !nk
              ENDDO !nv
            ENDIF !ADD
          ENDDO !njj
          
C*** If adding on to an existing mesh, then the first node in the added
C*** mesh is placed at the same location as a node in the exisitng
C*** mesh.i.e. the node does not exist in the new mesh. To do this an
C*** offset for the coordinates of the added mesh is calculated, the
C*** appended element in the existing mesh is adapted, and the appended
C*** node is updated with appropriate versions of coordinates and
C*** fields.
        IF(IOTYPE.NE.3.AND..NOT.ADD) THEN
            NPNODE(0,nr)=NPNODE(0,nr)+1
            CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase
     &        NP_R_M',ERROR,*9999)
            NPNODE(NPNODE(0,nr),nr)=np
          ELSE
            NPNODE_TEMP(nonode)=np
            NPLIST(nonode)=np
          ENDIF
        ENDDO !End of nonode loop

        NPLIST(0)=NTNODE
        
        IF(IOTYPE.NE.3) THEN
          CALL ISORT(NPNODE(0,nr),NPNODE(1,nr))
          NPT(nr)=0
          DO nonode=1,NPNODE(0,nr)
            IF(NPNODE(nonode,nr).GT.NPT(nr)) NPT(nr)=NPNODE(nonode,nr)
          ENDDO
          IF(NPT(nr).GT.NPT(0)) NPT(0)=NPT(nr)
C         Setup regions that interface nodes belong to
          CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
        ENDIF
        
      ENDIF !if ntnode.gt.0
      
C***  Read/write the element information, and if reading then set up
C***  the element field arrays at the same time
      IF(IOTYPE.EQ.3) THEN
        NTELEM=NEELEM(0,nr)
        IDATA(1)=NTELEM
      ENDIF
      FORMAT='(/$,'' The number of elements is [1]: '',I5)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NEM,LDATA,LDEFLT,
     &  RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NTELEM=IDATA(1)
      CALL ASSERT(NTELEM.GT.0,'>>No elements being defined',ERROR,*9999)
      CALL ASSERT(NTELEM.LE.NE_R_M,'>>Increase NE_R_M',ERROR,*9999)
      N2ELEM=1

C***  Check for symmetric elements      
      IF(IOTYPE.EQ.3) THEN
        IF(SYMMETRIC_NE)THEN
          ADATA(1)='Y'
          SYMMETRY_PROMPT=.TRUE.
        ELSE
          ADATA(1)='N'
          SYMMETRY_PROMPT=.FALSE.
        ENDIF
      ENDIF
      FORMAT='($,'' Are any elements symmetric [N]? '',A)'
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     &  RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        IF(ADATA(1).EQ.'Y')THEN
          SYMMETRIC_NE=.TRUE.
          SYMMETRY_PROMPT=.TRUE.
        ELSE
          IF(ADD.AND.SYMMETRIC_NE)THEN
            SYMMETRIC_NE=.TRUE.
          ELSE
            SYMMETRIC_NE=.FALSE.
          ENDIF
          SYMMETRY_PROMPT=.FALSE.
        ENDIF
      ENDIF
      
      IF(IOTYPE.NE.3.AND..NOT.ADD) THEN
        NEELEM(0,nr)=0
      ENDIF

      ne_parent=NELIST2(1)
      DO noelem=1,NTELEM
        IF(noelem.EQ.1) THEN
          IDEFLT(1)=NET(0)+1
        ELSE
          IDEFLT(1)=NEELEM(N2ELEM,nr)+1
        ENDIF
        WRITE(CHAR1,'(I5)') IDEFLT(1)
        IF(IOTYPE.EQ.3) THEN
          ne=NEELEM(noelem,nr)
          IDATA(1)=ne
        ENDIF
        FORMAT='(/$,'' Element number ['//CHAR1(1:5)//']: '',I5)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NEM,LDATA,
     &    LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ne=IDATA(1)
        IF(ADD)THEN
          NRE_TEMP(ne)=nr
        ELSE
          NRE(ne)=nr
        ENDIF
        
C*** The tree mesh has basis functions the same for all nj_geom
        nj_geometry=NJ_LOC(NJL_GEOM,1,nr)
        IF(noelem.EQ.1) THEN !set default basis function
          IDEFLT(1)=1
        ELSE IF(noelem.GT.1) THEN
          IDEFLT(1)=NBJ(nj_geometry,NEELEM(1,nr))
        ENDIF
        WRITE(CHAR5,'(I2)') IDEFLT(1)
        FORMAT='($,'' The basis function type for geometry'//
     '    ' is ['//CHAR5(1:2)//']: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=NBJ(nj_geometry,ne)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,LDATA,
     &    LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          nb=IDATA(1)
          DO nj=1,NJT
            IF(ADD)THEN
              NBJ_TEMP(nj,ne)=nb
            ELSE
              NBJ(nj,ne)=nb
            ENDIF
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                IF(ADD)THEN
                  NKJE_TEMP(nk,nn,nj,ne)=nk
                ELSE
                  NKJE(nk,nn,nj,ne)=nk
                ENDIF
              ENDDO !nk
            ENDDO !nn
          ENDDO !nj
        ENDIF ! IOTYPE.NE.3
        
        DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
          nj=NJ_LOC(NJL_FIEL,njj,nr)
          IF(noelem.EQ.1) THEN !set default basis function
            IDEFLT(1)=1
          ELSE IF(noelem.GT.1) THEN
            IDEFLT(1)=NBJ(nj,NEELEM(1,nr))
          ENDIF
          WRITE(CHAR1,'(I1)') njj
          WRITE(CHAR5,'(I2)') IDEFLT(1)
          FORMAT='($,'' The basis function type for field'//
     '      ' variable '//CHAR1(1:1)//' is ['//CHAR5(1:2)//']: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NBJ(nj,ne)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &      NBT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            nb=IDATA(1)
            IF(ADD)THEN
              NBJ_TEMP(nj,ne)=nb
            ELSE
              NBJ(nj,ne)=nb
            ENDIF
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                IF(ADD)THEN
                  NKJE_TEMP(nk,nn,nj,ne)=nk
                ELSE
                  NKJE(nk,nn,nj,ne)=nk
                ENDIF
              ENDDO !nk
            ENDDO !nn
          ENDIF ! IOTYPE.NE.3
        ENDDO !njj
        
        FIRST_VERSION=.TRUE.
        IF(IOTYPE.NE.3) THEN
C           initialise NVJE with default value
          DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj2,nr)
            IF(ADD)THEN
              nb=NBJ_TEMP(1,ne)
              DO nn=1,NNT(nb)
                NVJE_TEMP(nn,nb,nj,ne)=1
              ENDDO !nn
            ELSE
              nb=NBJ(1,ne) !nb for geometry
              DO nn=1,NNT(nb)
                NVJE(nn,nb,nj,ne)=1
              ENDDO !nn
            ENDIF
          ENDDO !njj2
          DO njj2=1,NJ_LOC(NJL_FIEL,0,nr)
            nj=NJ_LOC(NJL_FIEL,njj2,nr)
            IF(ADD)THEN
              nb=NBJ_TEMP(nj,ne) !nb for field
              DO nn=1,NNT(nb)
                NVJE_TEMP(nn,nb,nj,ne)=1
              ENDDO !nn
            ELSE
              nb=NBJ(nj,ne) !nb for field
              DO nn=1,NNT(nb)
                NVJE(nn,nb,nj,ne)=1
              ENDDO !nn
            ENDIF
          ENDDO !njj2
        ENDIF

        IF(ADD)THEN
          nb=NBJ_TEMP(1,ne)
        ELSE
          nb=NBJ(1,ne) !basis function for geometry
        ENDIF
        
        DO nn=1,NNT(nb)
          IDEFLT(nn)=1
        ENDDO !nn
        WRITE(CHAR2,'(I1)') NNT(nb)
        WRITE(CHAR5,'(I2)') nb
        
        FORMAT='($,'' Enter the '//CHAR2(1:1)//' global numbers for'//
     '    ' basis'//CHAR5(1:2)//': '',I5,14(1X,I5))'
        IF(IOTYPE.EQ.3) THEN
          DO nn=1,NNT(nb)
            IDATA(nn)=NPNE(nn,nb,ne)
          ENDDO !nn
        ENDIF
        CDATA(1)='NODES' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &    NNT(nb),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO nn=1,NNT(nb)
            IF(ADD)THEN
              NPNE_TEMP(nn,nb,ne)=IDATA(nn)
            ELSE
              NPNE(nn,nb,ne)=IDATA(nn)
            ENDIF
          ENDDO !nn
        ENDIF

        DO ntype=1,2 !geometry, fields 
          IF(ntype.EQ.1)THEN
            njtype=1 !geometry
          ELSE
            njtype=3 !field
          ENDIF
          IF(ADD)THEN
            nb=NBJ_TEMP(njtype,ne)
          ELSE
            nb=NBJ(njtype,ne)
          ENDIF
          DO nn=1,NNT(nb)
            IF(ADD)THEN
              np=NPNE_TEMP(nn,nb,ne)
            ELSE
              np=NPNE(nn,nb,ne)
            ENDIF
            np2=np
            DO njj2=1,NJ_LOC(njtype,0,nr)
              nj=NJ_LOC(njtype,njj2,nr)
              IF(ADD)THEN
                IF(NVJP_TEMP(nj,np).GT.1) THEN !>1 global version
                  WRITE(CHAR1,'(I5)') np2 !original format
                  WRITE(CHAR2,'(I1)') njj2 !using njj2'th field 
                  IDEFLT(1)=1
                  WRITE(CHAR4,'(I2)') IDEFLT(1)
                  IF(ntype.EQ.1)THEN
                    FORMAT='($,'' The version number for node '
     '                //CHAR1(1:5)//', nj='//CHAR2(1:1)//' is ['
     &                //CHAR4(1:2)//']: '',I2)'
                  ELSE
                    FORMAT='($,'' The version number for node '
     '                //CHAR1(1:5)//', field '//CHAR2(1:1)//' is ['
     &                //CHAR4(1:2)//']: '',I2)'
                  ENDIF
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,1,NVJP_TEMP(nj,np),LDATA,LDEFLT,RDATA,
     &              RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) NVJE_TEMP(nn,nb,nj,ne)=IDATA(1)
                  IF(FIRST_VERSION) FIRST_VERSION=.FALSE.
                ELSE
                  NVJE_TEMP(nn,nb,nj,ne)=1
                ENDIF !NVJP(nj,np).GT.1
              ELSE
                IF(NVJP(nj,np).GT.1) THEN !>1 global version
                  WRITE(CHAR1,'(I5)') np2 !original format
                  WRITE(CHAR2,'(I1)') njj2 !using njj2'th field 
                  IDEFLT(1)=1
                  WRITE(CHAR4,'(I2)') IDEFLT(1)
                  IF(ntype.EQ.1)THEN
                    FORMAT='($,'' The version number for node '
     '                //CHAR1(1:5)//', nj='//CHAR2(1:1)//' is ['
     &                //CHAR4(1:2)//']: '',I2)'
                  ELSE
                    FORMAT='($,'' The version number for node '
     '                //CHAR1(1:5)//', field '//CHAR2(1:1)//' is ['
     &                //CHAR4(1:2)//']: '',I2)'
                  ENDIF
                  IF(IOTYPE.EQ.3) IDATA(1)=NVJE(nn,nb,nj,ne)
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,1,NVJP(nj,np),LDATA,LDEFLT,RDATA,
     &              RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) NVJE(nn,nb,nj,ne)=IDATA(1)
                  IF(FIRST_VERSION) FIRST_VERSION=.FALSE.
                ELSE
                  NVJE(nn,nb,nj,ne)=1
                ENDIF !NVJP(nj,np).GT.1
              ENDIF
            ENDDO !njj2
          ENDDO !nn
        ENDDO !ntype
        
C***  Prompt for symmetry information if the mesh contains symmetric
C***  elements
        IF(SYMMETRY_PROMPT)THEN
          IF(IOTYPE.EQ.3) THEN
            IF(NORD(5,ne).EQ.1.OR.NORD(5,ne).EQ.-1)THEN !symmetric element
              ADATA(1)='Y'
              IDATA(1)=NORD(5,ne)                      
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          WRITE(CHAR1,'(I5)') ne
          FORMAT='($,'' Is element '//CHAR1(1:5)//
     &      ' symmetric [N]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &      FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
            FORMAT='($,'' Is symmetry diverging [1] or converging '//
     &        '(-1) for element '//CHAR1(1:5)//'? '',I2)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        -1,1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            SYMM_INDEX=IDATA(1)
          ENDIF
          IF(IOTYPE.NE.3)THEN
            SYMM_INDEX=IDATA(1) !symmetry index
            IF(ADD)THEN
              IF(ADATA(1).EQ.'Y')THEN
                NORD_TEMP(ne)=SYMM_INDEX
              ELSE
                NORD_TEMP(ne)=0
              ENDIF
            ELSE
              IF(ADATA(1).EQ.'Y')THEN
                NORD(5,ne+ne_offset)=SYMM_INDEX
              ELSE
                NORD(5,ne+ne_offset)=0
              ENDIF
            ENDIF
          ENDIF
        ELSE
          IF(ADD)THEN
            NORD_TEMP(ne)=0
          ELSE
            NORD(5,ne)=0
          ENDIF
        ENDIF !SYMMETRY
        
        IF(IOTYPE.NE.3.AND..NOT.ADD) THEN
          NEELEM(0,nr)=NEELEM(0,nr)+1
          CALL ASSERT(NEELEM(0,nr).LE.NE_R_M,'>>Increase NE_R_M',
     '      ERROR,*9999)
          NEELEM(NEELEM(0,nr),nr)=ne
        ELSE
          NEELEM_TEMP(noelem)=ne
        ENDIF
        IF(.NOT.ADD.AND.IOTYPE.NE.3)THEN
          NEELEM(0,0)=NEELEM(0,0)+1
        ENDIF
      ENDDO !noelem
      
      IF(IOTYPE.NE.3) THEN
        CALL ISORT(NEELEM(0,nr),NEELEM(1,nr))
        NET(nr)=0
        DO noelem=1,NEELEM(0,nr)
          IF(NEELEM(noelem,nr).GT.NET(nr))
     '      NET(nr)=NEELEM(noelem,nr)
          ne=NEELEM(noelem,nr)
        ENDDO
        IF(NET(nr).GT.NET(0)) NET(0)=NET(nr)
      ENDIF

      IF(ADD)THEN
        CALL ISORT(NPLIST(0),NPLIST(1)) !non-decreasing order, to get offset
        ne=NET(0)
        np=NPT(0)
        noelem=NEELEM(0,nr)
        nonode=NPNODE(0,nr)
        
        DO noelem_parent=1,NELIST2(0) !for each of the specified parents
          ne_parent=NELIST2(noelem_parent) !attachment element
          nb=NBJ(1,ne_parent) !basis function
          np_parent=NPNE(2,nb,ne_parent) !node for attachment
          np_offset=NPT(0)-NPLIST(1)+1
          IF(NPLIST(1).EQ.NP_ATTACH) np_offset=NPT(0)-NPLIST(2)+1
          DO nj=1,NJT
            XP_OFFSET(nj)=XP(1,1,nj,np_parent)-XP_TEMP(1,1,nj,
     &        NP_ATTACH)
          ENDDO !nj

C***      Set up values for attachment to parent element
          DO nonode2=1,NTNODE !for non-attachment nodes
            np2=NPNODE_TEMP(nonode2)
            IF(np2.EQ.NP_ATTACH)THEN
C***        Attach NP_ATTACH at np_parent
              DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
                nj=NJ_LOC(NJL_FIEL,njj,nr)
                DO nk=1,NKJP(nj)
                  DO nv=1,NVJP_TEMP(nj,np_parent)
                    XP(nk,nv+NVJP(nj,np_parent),nj,np_parent)=
     &                XP_TEMP(nk,nv,nj,np2)
                  ENDDO !nv
                  NVJP(nj,np_parent)=NVJP(nj,np_parent)+NVJP_TEMP(nj,
     &              np2)
                ENDDO !nk
                NKJ(nj,np_parent)=NKJP(nj)
              ENDDO !njj
            ELSE
              np=NPNODE_TEMP(nonode2)+np_offset
              nonode=nonode+1
              CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',ERROR,
     &          *9999)
              DO ntype=1,2
                IF(ntype.EQ.1)THEN
                  njtype=1 !geometry
                ELSE
                  njtype=3 !field
                ENDIF
                DO njj=1,NJ_LOC(njtype,0,nr)
                  nj=NJ_LOC(njtype,njj,nr)
                  NKJ(nj,np)=NKJ_TEMP(nj,np2)
                  NVJP(nj,np)=NVJP_TEMP(nj,np2)
                ENDDO !njj
              ENDDO !ntype
              
              DO nj=1,NJT !coordinates
                XP(1,1,nj,np)=XP_TEMP(1,1,nj,np2)+XP_OFFSET(nj)
                IF(NKJ(nj,np).GT.1)THEN
                  XP(2,1,nj,np)=XP_TEMP(2,1,nj,np2) !in case derivative is set
                ENDIF
                XP(1,2,nj,np)=XP_TEMP(1,2,nj,np2) !direction, temporary
              ENDDO !nj
              DO njj=1,NJ_LOC(NJL_FIEL,0,nr) !fields
                nj=NJ_LOC(NJL_FIEL,njj,nr)
                DO nk=1,NKJ(nj,np)
                  DO nv=1,NVJP(nj,np)
                    XP(nk,nv,nj,np)=XP_TEMP(nk,nv,nj,np2)
                  ENDDO !nv
                ENDDO !nk
              ENDDO !njj

              NP_INTERFACE(np,0)=1
              NP_INTERFACE(np,1)=nr
              NPNODE(nonode,nr)=np
              NPNODE(0,nr)=NPNODE(0,nr)+1
              NPNODE(0,0)=NPNODE(0,0)+1
              NPT(0)=NPT(0)+1
            ENDIF
          ENDDO !nonode2

          DO noelem2=1,NTELEM
            ne2=NEELEM_TEMP(noelem2)
            ne=ne+1
            noelem=noelem+1
            CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',ERROR,
     &        *9999)

            NRE(ne)=NRE_TEMP(ne2)
            DO ntype=1,2
              IF(ntype.EQ.1)THEN
                njtype=1 !geometry
              ELSE
                njtype=3 !field
              ENDIF
              DO njj=1,NJ_LOC(njtype,0,nr)
                nj=NJ_LOC(njtype,njj,nr)
                nb=NBJ_TEMP(nj,ne2)
                NBJ(nj,ne)=nb
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    NKJE(nk,nn,nj,ne)=NKJE_TEMP(nk,nn,nj,ne2)
                  ENDDO !nk
                  IF(NPNE_TEMP(nn,nb,ne2).EQ.NP_ATTACH)THEN
                    DO njj2=1,NJ_LOC(njtype,0,nr)
                      nj2=NJ_LOC(njtype,njj2,nr)
                      NVJE(nn,nb,nj2,ne)=1
                    ENDDO !njj2
                  ELSE
                    NVJE(nn,nb,nj,ne)=NVJE_TEMP(nn,nb,nj,ne2)
                  ENDIF
                ENDDO !nn
              ENDDO !njj
            ENDDO !ntype
            nb=NBJ(1,ne)
            DO nn=1,NNT(nb)
              IF(NPNE_TEMP(nn,nb,ne2).EQ.NP_ATTACH)THEN
                NPNE(nn,nb,ne)=NPNE(2,nb,ne_parent)
              ELSE
                NPNE(nn,nb,ne)=NPNE_TEMP(nn,nb,ne2)+np_offset
              ENDIF
            ENDDO!nn
            NORD(5,ne)=NORD_TEMP(ne2)
            
            NEELEM(noelem,nr)=ne
            NEELEM(0,nr)=NEELEM(0,nr)+1
            NEELEM(0,0)=NEELEM(0,0)+1
            NET(0)=NET(0)+1
          ENDDO !noelem2
          
        ENDDO !noelem_parent
      ENDIF !ADD
      nb=NBJ(1,NEELEM(1,nr)) !basis function for geometry
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      
      IF(IOTYPE.NE.3) THEN
        IF(NPT(0).GT.0) THEN
          CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO nb=1,NBFT
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,ne)=1.d0
              ENDDO !ns
            ENDDO !nb
            nb=NBJ(1,ne)
            np1=NPNE(1,nb,ne) !start node
            np2=NPNE(2,nb,ne) !end node
            length=0.d0
            DO nj=1,NJT
              length=length+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2
            ENDDO !nj
            length=DSQRT(length)
            DO nj=1,NJT
              XP(1,2,nj,np2)=(XP(1,1,nj,np2)-XP(1,1,nj,np1))/length
            ENDDO !nj
          ENDDO !noelem
        ENDIF
      ENDIF

      CALL ISORT(NPNODE(0,nr),NPNODE(1,nr))
      NPT(nr)=NPNODE(NPNODE(0,nr),nr) !highest node number

      CALL EXITS('IPMESH2_1')
      RETURN
 9999 CALL ERRORS('IPMESH2_1',ERROR)
      DO njj1=NJJ_START,NJ_LOC(NJL_FIEL,0,nr)
        NJ_LOC(NJL_FIEL,njj1,nr)=0
      ENDDO
      NJ_LOC(NJL_FIEL,0,nr)=NJJ_START-1
      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          IF(NJ_LOC(njj1,njj2,nr).GT.NJ_LOC(0,0,nr))
     '      NJ_LOC(0,0,nr)=NJ_LOC(njj1,njj2,nr)
        ENDDO !njj2
      ENDDO !njj1
      NJ_LOC(0,0,0)=0
      NJ_LOC(NJL_FIEL,0,0)=0
      DO nr=1,NRT
        IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0)) NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
        IF(NJ_LOC(NJL_FIEL,0,nr).GT.NJ_LOC(NJL_FIEL,0,0))
     '    NJ_LOC(NJL_FIEL,0,0)=NJ_LOC(NJL_FIEL,0,nr)
      ENDDO !nr
      CALL EXITS('IPMESH2_1')
      RETURN 1
      END

