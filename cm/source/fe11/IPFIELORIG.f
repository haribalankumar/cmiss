      SUBROUTINE IPFIELORIG(NEELEM,NENP,NJJ_END,
     &  NJJ_START,NKJ,NPNODE,nr,NUM_FIELD_DEFAULT,NVJP,XAB,XP,DERIVS,
     &  ELEM_FIELD,FLOW_NODAL,
     &  MECH,HYPOXIA,STRETCH,TREE,VALUES,VERSIONS,ERROR,*)

C#### Subroutine: IPFIEL
C###  Description:
C###    IPFIEL inputs additional field variables
C###    1..NJ_LOC(NJL_FIEL,0) for region nr.
C###    VALUES is a logical variable to indicate that only a subset of 
C###    fields are read/written. The field indices are specified by
C###    NJJ_START and NJJ_END.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mesh00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     &  NJJ_END,NJJ_START,NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),nr,
     &  NUM_FIELD_DEFAULT,
     &  NVJP(NJM,NPM)
      REAL*8 XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL DERIVS,ELEM_FIELD,FLOW_NODAL,MECH,HYPOXIA,STRETCH,TREE,
     &  VALUES,VERSIONS
!     Local Variables
      INTEGER FREELIST(NEJ_LOC_MX),IBEG,IEND,ICHAR,index,INFO,ne,nj,njj,
     '  njj1,nk,NKJP(3*NJ_LOC_MX),noelem,nonode,NOQUES,
     '  np,nrr,NTELEM,NTNODE,nu,NUK(8),NUM_FIELD,numf,NUMFREE,nv,NVJP_IT
C unused in move to IPFIEL_NJLOC: FREELIST(NJ_LOC_MX)      
      CHARACTER CHAR2*2,CHAR1*1,CHAR3*1,CHAR4*6,
     '  CHAR12*12
      LOGICAL FILEIP,PROMPT_NV(NJ_LOC_MX)

      DATA NUK/1,2,4,6,7,9,10,11/
c      DATA NJJ_START/1/ !Initialised for 9999 error handling

      CALL ENTERS('IPFIEL',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C      IF(NET(nr).GT.0) THEN
C        ELEM=.TRUE.
C      ELSE
C        ELEM=.FALSE.
C      ENDIF

C     Determine the first new field
      IF(.NOT.ELEM_FIELD)THEN
        IF(.NOT.VALUES)THEN !AJS fix for option to set start and end field number
          IF(ADD) THEN
            NJJ_START=NJ_LOC(NJL_FIEL,0,nr)+1
          ELSE
            NJJ_START=1
          ENDIF
        ENDIF

        WRITE(CHAR2,'(I2)') NJT
        CALL STRING_TRIM(CHAR2,IBEG,IEND)
        IDEFLT(1)=NUM_FIELD_DEFAULT
        FORMAT='($,'' The number of field variables is ['
     '    //CHAR2(IBEG:IEND)//']: '',I2)'
        IF(IOTYPE.EQ.3)THEN
          IDATA(1)=NJ_LOC(NJL_FIEL,0,nr)
          IF(VALUES) IDATA(1)=NJJ_END-NJJ_START+1 !AJS number of fields determined by start and end index
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NJM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NUM_FIELD=IDATA(1)
          IF(.NOT.VALUES)THEN
            CALL IPFIEL_NJLOC(NJJ_START,NKJ,NPNODE,nr,NUM_FIELD,NVJP,
     '        ERROR,*9999)
            NJJ_END=NJ_LOC(NJL_FIEL,0,nr)
          ENDIF
        ENDIF

        CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).LE.NJ_LOC_MX,
     '  '>>Increase dimension of PROMPT_NV',ERROR,*9999)
        IDEFLT(1)=NPNODE(0,nr)
        WRITE(CHAR4,'(I6)') IDEFLT(1)
        IF(IOTYPE.EQ.3) THEN
          NTNODE=NPNODE(0,nr)
          IDATA(1)=NTNODE
        ENDIF
        FORMAT='($,'' The number of nodes is ['//CHAR4(1:6)//']: '',I6)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NTNODE=IDATA(1)

C       ask for version prompting
        DO njj=NJJ_START,NJJ_END
c          DO njj=NJJ_START,NJ_LOC(NJL_FIEL,0,nr)
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            WRITE(CHAR2,'(I2)') njj
            CALL STRING_TRIM(CHAR2,IBEG,IEND)
            FORMAT='($,'' Do you want prompting for different versions '
     '      //'of field variable '//CHAR2(IBEG:IEND)//' [N]? '',A)'
            IF(IOTYPE.EQ.3) THEN
              PROMPT_NV(njj)=.FALSE.
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(NVJP(nj,np).GT.1)
     '            PROMPT_NV(njj)=.TRUE.
              ENDDO
              IF(PROMPT_NV(njj)) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ELSE
              PROMPT_NV(njj)=.TRUE.
              ADATA(1)='Y'
            ENDIF

            IF(VERSIONS.AND.IOTYPE.NE.3)THEN
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
            ELSE
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
            ENDIF
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                PROMPT_NV(njj)=.TRUE.
              ELSE
                PROMPT_NV(njj)=.FALSE.
              ENDIF
            ENDIF !iotype

          ENDDO !njj

!           DO njj=NJJ_START,NJ_LOC(NJL_FIEL,0,nr) !AJS fix for option to set start and end field number
          DO njj=NJJ_START,NJJ_END
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            WRITE(CHAR2,'(I2)') njj
            CALL STRING_TRIM(CHAR2,IBEG,IEND)
            IF(IOTYPE.EQ.3) THEN
C MPN 19Feb98: bug
              NKJP(nj)=0
C old            NKJP(NJ_LOC(NJL_FIEL,1,nr))=0
              DO nonode=1,NTNODE
                IF(NKJ(nj,NPNODE(nonode,nr)).GT.NKJP(nj)) THEN
                  NKJP(nj)=NKJ(nj,NPNODE(nonode,nr))
                ENDIF
              ENDDO
              IDATA(1)=NKJP(nj)-1
            ENDIF
            IF(DERIVS)THEN
              IDEFLT(1)=NKJ(1,NPNODE(1,nr))-1
            ELSE
             IDEFLT(1)=0
            ENDIF
            FORMAT='($,'' The number of derivatives for field variable '
     '        //CHAR2(IBEG:IEND)//' is [0]: '',I1)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     '        NKM-1,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '        *9999)
            IF(IOTYPE.NE.3) NKJP(nj)=IDATA(1)+1
          ENDDO !njj

C      ELSE
C        NTNODE=0
C      ENDIF

        IF(NTNODE.GT.0) THEN
          DO nonode=1,NTNODE
            IDEFLT(1)=NPNODE(nonode,nr)
            WRITE(CHAR4,'(I6)') IDEFLT(1)
            IF(IOTYPE.EQ.3) THEN
              np=NPNODE(nonode,nr)
              IDATA(1)=NP
            ENDIF
            FORMAT='(/$,'' Node number ['//CHAR4(1:6)//']: '',I6)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) np=IDATA(1)

c            DO njj=NJJ_START,NJ_LOC(NJL_FIEL,0,nr)
            DO njj=NJJ_START,NJJ_END
              nj=NJ_LOC(NJL_FIEL,njj,nr)
              WRITE(CHAR1,'(I1)') njj
              NKJ(nj,np)=NKJP(nj)
              IF(PROMPT_NV(njj)) THEN  !prompt for diff versions of field
                IF(IOTYPE.EQ.3) IDATA(1)=NVJP(nj,np)
                IF(VERSIONS)THEN !option for defining, not o/p
                  WRITE(CHAR3,'(I1)') NENP(np,0,nr) !number of elements at node
                  FORMAT='($,'' The number of versions for'
     '              //' field variable '//CHAR1//' is ['//CHAR3//
     &              ']: '',I2)'
                  IF(IOTYPE.NE.3) IDATA(1)=NENP(np,0,nr)
                  IF(IOTYPE.NE.3) IDEFLT(1)=NENP(np,0,nr)
                ELSE
                  FORMAT='($,'' The number of versions for'
     '              //' field variable '//CHAR1//' is [1]: '',I2)'
                  IF(VERSIONS)THEN
                    IF(IOTYPE.NE.3)  IDEFLT(1)=NVJP(1,np)
                  ELSE
                    IF(IOTYPE.NE.3) IDEFLT(1)=1
                  ENDIF
                ENDIF
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &            1,NVM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &            *9999)
                IF(IOTYPE.NE.3)THEN
                  IF(.NOT.VALUES)THEN
                    NVJP(nj,np)=IDATA(1)
                  ENDIF
                  NVJP_IT=IDATA(1)
                ELSE
                  NVJP_IT=IDATA(1)
                ENDIF
              ELSE
                NVJP(nj,np)=1
                NVJP_IT=1
              ENDIF
c              DO nv=1,NVJP(nj,np)
              DO nv=1,NVJP_IT
c                IF(NVJP(nj,np).GT.1) THEN !ask for diff nj versions
                IF(NVJP_IT.GT.1) THEN !ask for diff nj versions
                  WRITE(CHAR2,'(I2)') nv
                  CALL STRING_TRIM(CHAR2,IBEG,IEND)
                  FORMAT='('' For version number '//CHAR2(IBEG:IEND)
     '              //':'')'
                  CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,
     '              FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,
     '              IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '              INFO,ERROR,*9999)
c               CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
                ENDIF
                WRITE(CHAR2,'(I2)') njj
                CALL STRING_TRIM(CHAR2,IBEG,IEND)
                DO nk=1,NKJP(nj)
                  IF(nonode.EQ.1) THEN
                    RDEFLT(1)=0.0d0
                  ELSE IF(np.GT.1) THEN
C LKC 10-SEP-97 add assert line below - for accessing XP array
C LKC 13-OCT-97 Moved from outside DO loop
                    CALL ASSERT(NKM.GE.nk,'>>Increase NKM',ERROR,*9999)
                    IF(IOTYPE.NE.2) THEN
                      IF(nv.GT.1) THEN
                        ! Use previous version of same node
                        RDEFLT(1)=XP(nk,nv-1,nj,np)
                      ELSE
                        RDEFLT(1)=XP(nk,nv,nj,NPNODE(nonode-1,nr))
                      ENDIF
                    ELSE
                      RDEFLT(1)=0.d0 !if reading this value not used anyway
                    ENDIF
                  ENDIF
                  WRITE(CHAR12,'(D12.5)') RDEFLT(1)
                  nu=NUK(nk)
                  IF(nu.EQ.1) THEN
                    FORMAT='($,'' The field variable '//CHAR2(IBEG:IEND)
     '                //' value is ['//CHAR12//']: '',D12.5)'
                  ELSE IF(nu.EQ.2.OR.nu.EQ.4.OR.nu.EQ.7) THEN
                    IF(nu.EQ.2) THEN
                      CHAR1='1'
                    ELSE IF(nu.EQ.4) THEN
                      CHAR1='2'
                    ELSE IF(nu.EQ.7) THEN
                      CHAR1='3'
                    ENDIF
                    FORMAT='($,'' The derivative wrt direction '//
     '                CHAR1//' is ['//CHAR12//']: '',D12.5)'
                  ELSE IF(nu.EQ.6.OR.nu.EQ.9.OR.nu.EQ.10) THEN
                    IF(nu.EQ.6) THEN
                      CHAR1='1'
                      CHAR3='2'
                    ELSE IF(nu.EQ.9) THEN
                      CHAR1='1'
                      CHAR3='3'
                    ELSE IF(nu.EQ.10) THEN
                      CHAR1='2'
                      CHAR3='3'
                    ENDIF
                    FORMAT='($,'' The derivative wrt directions '//
     '                CHAR1//' & '//CHAR3//
     '                ' is ['//CHAR12//']: '',D12.5)'
                  ELSE IF(nu.EQ.11) THEN
                    FORMAT='($,'' The derivative wrt directions '//
     '                  '1, 2 & 3 is ['//CHAR12//']: '',D12.5)'
                  ENDIF
                  IF(IOTYPE.EQ.3) RDATA(1)=XP(nk,nv,nj,np)
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '              IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
C GMH 8/1/97 Update cmgui link
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    XP(nk,nv,nj,np)=RDATA(1)
                  ENDIF
                ENDDO !nk
              ENDDO !nv
            ENDDO !njj
          ENDDO !nonode

C        IF(IOTYPE.NE.3) THEN
C          IF(ELEM) THEN
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
C                nj=NJ_LOC(NJL_FIEL,njj,nr)
C                nb=NBJ(nj,ne)
C                DO nn=1,NNT(nb)
C                  NVJE(nn,nb,nj,ne)=1 !default version
C                ENDDO !nn
C              ENDDO !nj
C            ENDDO !noelem
C            REDEFINE=.FALSE.
C            DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
C              IF(PROMPT_NV(njj)) REDEFINE=.TRUE.
C            ENDDO
C            IF(REDEFINE) THEN
C              WRITE(OP_STRING,'('' >> Redefine elements to use correct '
C     '          //'versions for the field'')')
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDIF
C        ENDIF
        ENDIF !ntnode>0

        IF(MECH) nj_pleural=NJ_LOC(NJL_FIEL,NJJ_START,nr)
        IF(STRETCH) nj_stretch=NJ_LOC(NJL_FIEL,NJJ_START,nr)
        IF(FLOW_NODAL) THEN
          nj_flow=NJ_LOC(NJL_FIEL,NJJ_START,nr)
c          NJ_LOC(NJL_FIEL,0,nr)=NJJ_START-1
          write(*,*) "CHANGE CODE nj_flow=",nj_flow 
          !Use WRITE(OP_STRING...) and CALL WRITES for this - see above (in comments)
        ENDIF
        IF(HYPOXIA) THEN
          nj_hypoxia=NJ_LOC(NJL_FIEL,NJJ_START,nr)
c          NJ_LOC(NJL_FIEL,0,nr)=NJJ_START-1
          !write(*,*) "nj_hypoxia=",nj_hypoxia,NJL_FIEL,NJJ_START
          !Use WRITE(OP_STRING...) and CALL WRITES for this - see above (in comments)
        ENDIF  

      ELSE   !ELEM_FIELD
        IF(ADD) THEN
          NJJ_START=NEJ_LOC(0,nr)+1
        ELSE
          NJJ_START=1
        ENDIF

        WRITE(CHAR2,'(I2)') NJT
        CALL STRING_TRIM(CHAR2,IBEG,IEND)
        IDEFLT(1)=NUM_FIELD_DEFAULT
        FORMAT='($,'' The number of field variables is ['
     '    //CHAR2(IBEG:IEND)//']: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=NEJ_LOC(0,nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NJM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NUM_FIELD=IDATA(1)
          nj=NEJ_LOC(0,nr)
          NUMFREE=0
          DO WHILE((NUMFREE.NE.NUM_FIELD).AND.(nj.LE.NJM))
            nj=nj+1
            CALL ASSERT(nj.LE.NEJ_LOC_MX,'>>Increase NEJ_LOC_MX in '
     '        //'loc00.cmn',ERROR,*9999)
            NUMFREE=NUMFREE+1
            CALL ASSERT(NUMFREE.LE.NEJ_LOC_MX,'>>Increase NJ_LOC_MX in '
     '        //'loc00.cmn',ERROR,*9999)
            FREELIST(NUMFREE)=nj
          ENDDO
          CALL ASSERT(nj.LE.NJM,' >>Increase NJM',ERROR,*9999)
C     Clear any existing fields unless doing an add
          IF(.NOT.ADD) THEN
            DO njj1=1,NEJ_LOC(0,nr)
              NEJ_LOC(njj1,nr)=0
            ENDDO
            NEJ_LOC(0,nr)=0
          ENDIF !NOT ADD
          index=NJJ_START
          DO numf=1,NUMFREE
            nj=FREELIST(numf)
            NEJ_LOC(index,nr)=nj
            IF(nj.GT.NEJ_LOC(0,nr)) NEJ_LOC(0,nr)=nj
            index=index+1
          ENDDO
          NEJ_LOC(0,nr)=NUM_FIELD
         DO nrr=1,NRT
            IF(NEJ_LOC(0,nrr).GT.NEJ_LOC(0,0))
     '        NEJ_LOC(0,0)=NEJ_LOC(0,nrr)
          ENDDO !nrr
        ENDIF

        IDEFLT(1)=NEELEM(0,nr)
        WRITE(CHAR4,'(I6)') IDEFLT(1)
        IF(IOTYPE.EQ.3) THEN
          NTELEM=NEELEM(0,nr)
          IDATA(1)=NTELEM
        ENDIF
        FORMAT='($,'' The number of elements is ['//CHAR4(1:6)//
     '     ']: '',I6)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NTELEM=IDATA(1)
        IF(NTELEM.GT.0) THEN
          DO noelem=1,NTELEM
            IDEFLT(1)=NEELEM(noelem,nr)
            WRITE(CHAR4,'(I6)') IDEFLT(1)
            IF(IOTYPE.EQ.3) THEN
              ne=NEELEM(noelem,nr)
              IDATA(1)=ne
            ENDIF
            FORMAT='(/$,'' Element number ['//CHAR4(1:6)//']: '',I6)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ne=IDATA(1)

            DO njj=NJJ_START,NEJ_LOC(0,nr)
              nj=NEJ_LOC(njj,nr)
              WRITE(CHAR1,'(I1)') njj
              CALL STRING_TRIM(CHAR1,IBEG,IEND)
              CALL STRING_TRIM(CHAR2,IBEG,IEND)
              IF(noelem.EQ.1)THEN
                RDEFLT(1)=0.d0
              ELSE
                IF(IOTYPE.NE.2) THEN
                  RDEFLT(1)=XAB(nj,NEELEM(noelem-1,nr))
                ELSE
                  RDEFLT(1)=0.d0 !if reading this value not used anyway
                ENDIF
              ENDIF
              WRITE(CHAR12,'(D12.5)') RDEFLT(1)
              FORMAT='($,'' The field variable '//CHAR2(IBEG:IEND)
     '          //' value is ['//CHAR12//']: '',D12.5)'
              IF (IOTYPE.EQ.3) RDATA(1)=XAB(nj,ne)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '          NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '         ERROR,*9999)
              IF(IOTYPE.NE.3)then
                XAB(nj,ne)=RDATA(1)
              endif
            ENDDO !njj
          ENDDO !nonode

C        IF(IOTYPE.NE.3) THEN
C          IF(ELEM) THEN
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
C                nj=NJ_LOC(NJL_FIEL,njj,nr)
C                nb=NBJ(nj,ne)
C                DO nn=1,NNT(nb)
C                  NVJE(nn,nb,nj,ne)=1 !default version
C                ENDDO !nn
C              ENDDO !nj
C            ENDDO !noelem
C            REDEFINE=.FALSE.
C            DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
C              IF(PROMPT_NV(njj)) REDEFINE=.TRUE.
C            ENDDO
C            IF(REDEFINE) THEN
C              WRITE(OP_STRING,'('' >> Redefine elements to use correct '
C     '          //'versions for the field'')')
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDIF
C        ENDIF
        ENDIF !ntnode>0
      ENDIF !elem_field

      CALL EXITS('IPFIEL')
      RETURN
 9999 CALL ERRORS('IPFIEL',ERROR)
      CALL EXITS('IPFIEL')
      RETURN 1
      END


