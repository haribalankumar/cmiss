      SUBROUTINE IPFIBR(NKJ,NPNODE,nr,NVJP,XP,
     '  FROM_TYPE,ERROR,*)

C#### Subroutine: IPFIBR
C###  Description:
C###    IPFIBR inputs fibre direction field for region nr.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVJP(NJM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER FROM_TYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER FREELIST(NJ_LOC_MX),IBEG1,ICHAR,IEND1,INFO,
     '  nj,njj,njj1,njj2,nk,NKJP(NJM),nonode,
     '  NOQUES,np,np2,nrr,NTNODE,nu,NUK(8),numf,NUM_FIBRES,NUMFREE,nv
      CHARACTER CHAR1*5,CHAR2*1,CHAR3*1,CHAR4*12,CHAR5*2,
     '  TYPE*11
      LOGICAL FILEIP,PROMPT_NV(3)

      DATA NUK/1,2,4,6,7,9,10,11/

      CALL ENTERS('IPFIBR',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C      IF(NET(nr).GT.0) THEN
C        ELEM=.TRUE.
C      ELSE
C        ELEM=.FALSE.
C      ENDIF

C cpb 10/3/97 Merging sheets with fibres

C      FORMAT='('' Specify whether [1]: '''//
C     '  '/''   (1) fibre angle in Xi_1&2 plane'''//
C     '  '/''   (2) fibre and imbrication angles'''//
C     '  '/''   (3) '''//
C     '  '/$,''    '',I1)'
C      IF(IOTYPE.EQ.3) THEN
C        IF(NJ_LOC(NJL_FIBR,0,nr).EQ.1) THEN
C          IDATA(1)=1
C        ELSE
C          IDATA(1)=2
C        ENDIF
C      ENDIF
      IDEFLT(1)=1
      FORMAT='($,'' The number of fibre variables is [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NJ_LOC(NJL_FIBR,0,nr)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NUM_FIBRES=IDATA(1)

C ***   Check to see if there is enough room to store the fibres
        nj=0
        NUMFREE=0
        DO WHILE((NUMFREE.NE.NUM_FIBRES).AND.(nj.LE.NJM))
          nj=nj+1
          CALL ASSERT(nj.LE.3*NJ_LOC_MX,'>>Increase NJ_LOC_MX in '
     '      //'geom00.cmn',ERROR,*9999)
          IF(NJ_TYPE(nj,1).EQ.0.OR.NJ_TYPE(nj,1).EQ.NJL_FIBR) THEN
            !Empty space or fibre info in that nj location
            NUMFREE=NUMFREE+1
            CALL ASSERT(NUMFREE.LE.NJ_LOC_MX,'>>Increase NJ_LOC_MX in '
     '        //'geom00.cmn',ERROR,*9999)
            FREELIST(NUMFREE)=NJ
          ENDIF
        ENDDO
        CALL ASSERT(nj.LE.NJM,' >>Increase NJM',ERROR,*9999)

C ***   Clear any existing fibres etc
        DO njj1=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,njj1,nr)
          NJ_TYPE(nj,1)=0
          NJ_TYPE(nj,2)=0
          NJ_LOC(NJL_FIBR,njj1,nr)=0
        ENDDO
        NJ_LOC(NJL_FIBR,0,nr)=0
        NJ_LOC(0,0,nr)=0
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=NJ
          ENDDO
        ENDDO
C ***   Store the fibres in free space
        DO numf=1,NUMFREE
          nj=FREELIST(numf)
          NJ_LOC(NJL_FIBR,numf,nr)=nj
          NJ_TYPE(nj,1)=NJL_FIBR
          NJ_TYPE(nj,2)=numf
          IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
          IF(nj.GT.NJ_LOC(NJL_FIBR,0,0)) NJ_LOC(NJL_FIBR,0,0)=nj
        ENDDO
        NJ_LOC(NJL_FIBR,0,nr)=NUM_FIBRES
        DO nrr=1,NRT
          IF(NJ_LOC(0,0,nrr).GT.NJ_LOC(0,0,0))
     '      NJ_LOC(0,0,0)=NJ_LOC(0,0,nrr)
        ENDDO !NRR
C KAT 21May99: Initialize NKJ in case nodes are not set and CALC_NUNK
C       crashes.  NVJP is initialized in FEMINI.  Check this in case the
C       node has been defined in another region.
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
            nj=NJ_LOC(NJL_FIBR,njj,nr)
            IF(NVJP(nj,np).EQ.0) NKJ(nj,np)=0
          ENDDO !nj
        ENDDO !np
      ENDIF

      FORMAT='('' Specify how fibre angle is defined [1]: '''//
     '  '/''   (1) wrt Xi(1) coordinate'''//
     '  '/''   (2) wrt Xi(2) coordinate'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=JTYP12
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) JTYP12=IDATA(1)
      CALL ASSERT(JTYP12.EQ.1,'>>>ERROR: Sorry this option is not '
     '  //'yet implemented',ERROR,*9999)

      FORMAT='('' Specify whether angles entered in [1]: '''//
     '  '/''   (1) degrees'''//
     '  '/''   (2) radians'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=JTYP13
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) JTYP13=IDATA(1)

C cpb 10/3/97 fibre element information now set up witht a define
C element fibre

Cc cpb 16/9/95 Adding element dependent fibre basis types
C      FORMAT='($,'' Is the basis function for the fibre angle '
C     '    //'element dependent [N]? '',A)'
C      CDEFLT(1)='N'
C      IF(IOTYPE.EQ.3) THEN
C        ADATA(1)='N'
C        nb1=NBJ(NJ_LOC(NJL_FIBR,1,nr),NEELEM(1,nr))
C        DO noelem=2,NEELEM(0,nr)
C          ne=NEELEM(noelem,nr)
C          nb2=NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
C          IF(nb2.NE.nb1) ADATA(1)='Y'
C        ENDDO !ne
C      ENDIF
C      CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C      IF(ADATA(1).EQ.'N') THEN
C        BYELEMFIB=.FALSE.
C        FORMAT='($,'' The basis function type number for the '
C     '    //'fibre angle is [1]: '',I1)'
C        IF(IOTYPE.EQ.3) THEN
C          NB_FIBRE=NBJ(NJ_LOC(NJL_FIBR,1,nr),NEELEM(1,nr))
C          IDATA(1)=NB_FIBRE
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) THEN
C          NB_FIBRE=IDATA(1)
C          DO noelem=1,NEELEM(0,nr)
C            ne=NEELEM(noelem,nr)
C            NB_GEOM=NBJ(1,ne)
C            IF(NIT(NB_GEOM).NE.NIT(NB_FIBRE)) THEN
C              NB_FIBRE=NB_GEOM
C            ENDIF
C            NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)=NB_FIBRE
C          ENDDO !noelem (ne)
C          DO nonode=1,NPNODE(0,nr)
C            np=NPNODE(nonode,nr)
C            NKJ(NJ_LOC(NJL_FIBR,1,nr),np)=NKT(0,NB_FIBRE)
C          ENDDO !nonode (np)
C        ENDIF
C      ELSE
C        CALL ASSERT(CALL_ELEM,'>>Define elements first',ERROR,*9999)
C        BYELEMFIB=.TRUE.
C        NKTMAX=1
C        DO noelem=1,NEELEM(0,nr)
C          ne=NEELEM(noelem,nr)
C          IDEFLT(1)=NBJ(1,ne)
C          WRITE(CHAR2,'(I1)') IDEFLT(1)
C          WRITE(CHAR6,'(I5)') ne
C          FORMAT='($,'' The basis function type number for the fibre '
C     '      //'angle in element '//CHAR6(1:5)//' is ['//CHAR2//']: '','
C     '      //'I1)'
C          IF(IOTYPE.EQ.3) THEN
C            NB_FIBRE=NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
C            IDATA(1)=NB_FIBRE
C          ENDIF
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            NB_FIBRE=IDATA(1)
C            NB_GEOM=NBJ(1,ne)
C            IF(NIT(NB_GEOM).NE.NIT(NB_FIBRE)) THEN
C              NB_FIBRE=NB_GEOM
C            ENDIF
C            NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)=NB_FIBRE
C          ENDIF !iotype
C          IF(NKT(0,NB_FIBRE).GT.NKTMAX) NKTMAX=NKT(0,NB_FIBRE)
C        ENDDO !noelem
C        DO nonode=1,NPNODE(0,nr)
C          np=NPNODE(nonode,nr)
C          NKJ(NJ_LOC(NJL_FIBR,1,nr),np)=NKTMAX
C        ENDDO !nonode (np)
C      ENDIF
C
C      IF(NJL_LOC(NJL_FIBR,0,nr).GE.2) THEN
Cc cpb 16/9/95 Adding element dependent fibre basis types
C        FORMAT='($,'' Is the basis function for the imbrication angle '
C     '      //'element dependent [N]? '',A)'
C        CDEFLT(1)='N'
C        IF(IOTYPE.EQ.3) THEN
C          ADATA(1)='N'
C          nb1=NBJ(NJ_LOC(NJL_FIBR,2,nr),NEELEM(1,nr))
C          DO noelem=2,NEELEM(0,nr)
C            ne=NEELEM(noelem,nr)
C            nb2=NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)
C            IF(nb2.NE.nb1) ADATA(1)='Y'
C          ENDDO !ne
C        ENDIF
C        CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(ADATA(1).EQ.'N') THEN
C          BYELEMIMB=.FALSE.
C          FORMAT='($,'' The basis function type number for the '
C     '      //'imbrication angle is [1]: '',I1)'
C          IF(IOTYPE.EQ.3) THEN
C            NB_IMBRIC=NBJ(NJ_LOC(NJL_FIBR,2,nr),NEELEM(1,nr))
C            IDATA(1)=NB_IMBRIC
C          ENDIF
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            NB_IMBRIC=IDATA(1)
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)=NB_IMBRIC
C            ENDDO
C            DO nonode=1,NPNODE(0,nr)
C              np=NPNODE(nonode,nr)
C              NKJ(NJ_LOC(NJL_FIBR,2,nr),np)=NKT(0,NB_IMBRIC)
C            ENDDO
C          ENDIF
C        ELSE
C          CALL ASSERT(CALL_ELEM,'>>Define elements first',ERROR,*9999)
C          BYELEMIMB=.TRUE.
C          NKTMAX=1
C          DO noelem=1,NEELEM(0,nr)
C            ne=NEELEM(noelem,nr)
C            IDEFLT(1)=NBJ(1,ne)
C            WRITE(CHAR2,'(I1)') IDEFLT(1)
C            WRITE(CHAR6,'(I5)') ne
C            FORMAT='($,'' The basis function type number for the '
C     '        //'imbrication angle in element '//CHAR6(1:5)
C     '        //' is ['//CHAR2//']: '',I1)'
C            IF(IOTYPE.EQ.3) THEN
C              NB_IMBRIC=NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)
C              IDATA(1)=NB_IMBRIC
C            ENDIF
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(IOTYPE.NE.3) THEN
C              NB_IMBRIC=IDATA(1)
C              NB_GEOM=NBJ(1,ne)
C              IF(NIT(NB_GEOM).NE.NIT(NB_IMBRIC)) THEN
C                NB_IMBRIC=NB_GEOM
C              ENDIF
C              NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)=NB_IMBRIC
C            ENDIF !iotype
C            IF(NKT(0,NB_IMBRIC).GT.NKTMAX) NKTMAX=NKT(0,NB_IMBRIC)
C          ENDDO !noelem
C          DO nonode=1,NPNODE(0,nr)
C            np=NPNODE(nonode,nr)
C            NKJ(NJ_LOC(NJL_FIBR,2,nr),np)=NKTMAX
C          ENDDO !nonode (np)
C        ENDIF
C      ENDIF

C      IF(NNT(NB_FIBRE).GT.0) THEN


        IDEFLT(1)=NPNODE(0,nr)
        WRITE(CHAR1,'(I5)') IDEFLT(1)
        IF(IOTYPE.EQ.3) THEN
          NTNODE=NPNODE(0,nr)
          IDATA(1)=NTNODE
        ENDIF
        FORMAT='($,'' The number of nodes is ['//CHAR1(1:5)//']: '',I5)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NTNODE=IDATA(1)

!news MPN 16-Nov-94
        !ask for version prompting

        DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,njj,nr)
          IF(njj.EQ.1) THEN
            FORMAT='($,'' Do you want prompting for different versions '
     '        //'of the fibre angle [N]? '',A)'
          ELSE IF(njj.EQ.2) THEN
            FORMAT='($,'' Do you want prompting for different versions '
     '        //'of the imbrication angle [N]? '',A)'
          ELSE IF(njj.EQ.3) THEN
            FORMAT='($,'' Do you want prompting for different versions '
     '        //'of the sheet angle [N]? '',A)'
          ENDIF
          IF(IOTYPE.EQ.3) THEN
            PROMPT_NV(njj)=.FALSE.
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(NVJP(nj,np).GT.1) PROMPT_NV(njj)=.TRUE.
            ENDDO
            IF(PROMPT_NV(njj)) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              PROMPT_NV(njj)=.TRUE.
            ELSE
              PROMPT_NV(njj)=.FALSE.
            ENDIF
          ENDIF
        ENDDO !njj


        DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,njj,nr)
          IF(njj.EQ.1) THEN
            FORMAT='($,'' The number of derivatives for the fibre angle'
     '        //' is [0]: '',I1)'
          ELSE IF(njj.EQ.2) THEN
            FORMAT='($,'' The number of derivatives for the imbrication'
     '        //' angle is [0]: '',I1)'
          ELSE IF(njj.EQ.3) THEN
            FORMAT='($,'' The number of derivatives for the sheet angle'
     '        //' is [0]: '',I1)'
          ENDIF
          IF(IOTYPE.EQ.3) THEN
C           Find the maximum number of derivatives in the node list
            NKJP(nj)=0
            DO nonode=1,NTNODE
              IF(NKJ(nj,NPNODE(nonode,nr)).GT.NKJP(nj)) THEN
                NKJP(nj)=NKJ(nj,NPNODE(nonode,nr))
              ENDIF
            ENDDO
            IDATA(1)=NKJP(nj)-1
          ENDIF
          IDEFLT(1)=0
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NKM-1,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NKJP(nj)=IDATA(1)+1
        ENDDO !njj

C      ELSE
C        NTNODE=0
C      ENDIF !nnt

      IF(NTNODE.GT.0) THEN
        np2=0
        DO nonode=1,NTNODE
          IDEFLT(1)=NPNODE(nonode,nr)
          WRITE(CHAR1,'(I5)') IDEFLT(1)
          IF(IOTYPE.EQ.3) THEN
            np=NPNODE(nonode,nr)
            IDATA(1)=NP
          ENDIF
          FORMAT='(/$,'' Node number ['//CHAR1(1:5)//']: '',I5)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) np=IDATA(1)

          DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
            nj=NJ_LOC(NJL_FIBR,njj,nr)
            WRITE(CHAR1,'(I1)') nj
            NKJ(nj,np)=NKJP(nj)
            IF(njj.EQ.1) THEN
              TYPE='fibre'
            ELSE IF(njj.EQ.2) THEN
              TYPE='imbrication'
            ELSE IF(njj.EQ.3) THEN
              TYPE='sheet'
            ENDIF
            CALL STRING_TRIM(TYPE,IBEG1,IEND1)
            IF(PROMPT_NV(njj)) THEN !prompt for diff versions of fibres
              IF(IOTYPE.EQ.3)  THEN
                IDATA(1)=NVJP(nj,np)
              ELSE
                IF (FROM_TYPE(1:8).EQ."GEOMETRY") THEN
                  IDATA(1)=NVJP(1,np)
                ENDIF
              ENDIF
              IDEFLT(1)=IDATA(1)
              WRITE(CHAR5,'(I2)') IDEFLT(1)
              IF(njj.EQ.1) THEN
                FORMAT='($,'' The number of versions for the fibre '
     '            //'angle is ['//CHAR5(1:2)//']: '',I2)'
              ELSE IF(njj.EQ.2) THEN
                FORMAT='($,'' The number of versions for the '
     '            //'imbrication angle is ['//CHAR5(1:2)//']: '',I2)'
              ELSE IF(njj.EQ.3) THEN
                FORMAT='($,'' The number of versions for the sheet '
     '            //'angle is ['//CHAR5(1:2)//']: '',I2)'
              ENDIF
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          NVM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) NVJP(nj,np)=IDATA(1)
            ELSE
              NVJP(nj,np)=1
            ENDIF

            DO nv=1,NVJP(nj,np)
              IF(NVJP(nj,np).GT.1) THEN
C               ask for different nj versions
                WRITE(CHAR5,'(I2)') nv
                FORMAT='('' For version number'//CHAR5(1:2)//':'')'
                CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
              ENDIF

              DO nk=1,NKJP(nj)
                IF(np2.NE.0) THEN
C                 Only useful if it has enough versions and derivatives
                  IF(nv.GT.NVJP(nj,np2).OR.nk.GT.NKJ(nj,np2)) np2=0
                ENDIF
                IF(np2.EQ.0) THEN
C CS 29/1/2001 removing PI/2 rotation see MAT_VEC_ROTATE
C                  IF(njj.EQ.3) THEN ! defult for sheet angle
C                    IF(JTYP13.EQ.1) THEN ! in degrees
C                      RDEFLT(1)=90.0d0
C                    ELSE ! radians
C                      RDEFLT(1)=PI/2.0d0
C                    ENDIF
C                  ELSE
C                    RDEFLT(1)=0.0d0
C                  ENDIF
                  RDEFLT(1)=0.0d0
                ELSE
                  IF(JTYP13.EQ.1) THEN !in degrees
                    RDEFLT(1)=XP(nk,nv,nj,np2)*180.0d0/PI
                  ELSE IF(JTYP13.EQ.2) THEN !in radians
                    RDEFLT(1)=XP(nk,nv,nj,np2)
                  ENDIF
                ENDIF
                WRITE(CHAR4,'(D12.5)') RDEFLT(1)
                nu=NUK(nk)
                IF(nu.EQ.1) THEN
                  FORMAT='($,'' The '//TYPE(IBEG1:IEND1)//' angle '
     '              //'is ['//CHAR4(1:12)//']: '',D12.5)'
                ELSE IF(nu.EQ.2.OR.nu.EQ.4.OR.nu.EQ.7) THEN
                  IF(nu.EQ.2) THEN
                    CHAR2='1'
                  ELSE IF(nu.EQ.4) THEN
                    CHAR2='2'
                  ELSE IF(nu.EQ.7) THEN
                    CHAR2='3'
                  ENDIF
                  FORMAT='($,'' The derivative wrt direction '//
     '              CHAR2(1:1)//' is ['//CHAR4(1:12)//']: '',D12.5)'
                ELSE IF(nu.EQ.6.OR.nu.EQ.9.OR.nu.EQ.10) THEN
                  IF(nu.EQ.6) THEN
                    CHAR2='1'
                    CHAR3='2'
                  ELSE IF(nu.EQ.9) THEN
                    CHAR2='1'
                    CHAR3='3'
                  ELSE IF(nu.EQ.10) THEN
                    CHAR2='2'
                    CHAR3='3'
                  ENDIF
                  FORMAT='($,'' The derivative wrt directions '//
     '              CHAR2(1:1)//' & '//CHAR3(1:1)//
     '              ' is ['//CHAR4(1:12)//']: '',D12.5)'
                ELSE IF(nu.EQ.11) THEN
                  FORMAT='($,'' The derivative wrt directions '//
     '                '1, 2 & 3 is ['//CHAR4(1:12)//']: '',D12.5)'
                ENDIF
                IF(IOTYPE.EQ.3) THEN
                  IF(JTYP13.EQ.1) THEN
                    RDATA(1)=XP(nk,nv,nj,np)*180.0D0/PI
                  ELSE
                    RDATA(1)=XP(nk,nv,nj,np)
                  ENDIF
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  IF(JTYP13.EQ.1) THEN !in degrees
                    XP(nk,nv,nj,np)=RDATA(1)*PI/180.0D0
                  ELSE IF(JTYP13.EQ.2) THEN !in radians
                    XP(nk,nv,nj,np)=RDATA(1)
                  ENDIF
                ENDIF
              ENDDO !nk
            ENDDO !nv
          ENDDO !njj
          np2=np
        ENDDO !np


C cpb 10/3/97 fibre element information now set up with a define
C element fibre


C!news MPN 16-Nov-94
CC  Setup version info for extended basis
C        nb_extended=1
C        DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
C          nb_extended=nb_extended+1
C        ENDDO
C        IF(NBC(nb_extended).EQ.7) THEN
C          WRITE(OP_STRING,
C     '      '('' Version numbers updated for extended basis.'')')
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        ENDIF
C        IF(IOTYPE.NE.3) THEN
C          IF(ELEM) THEN
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              NB_FIBRE=NBJ(NJ_LOC(NJL_FIBR,1,nr),ne)
C              DO nn=1,NNT(NB_FIBRE)
C                NVJE(nn,NB_FIBRE,NJ_LOC(NJL_FIBR,1,nr),ne)=1 !deflt version
C                IF(NBC(nb_extended).EQ.7) THEN
C                  NVJE(nn,nb_extended,NJ_LOC(NJL_FIBR,1,nr),ne)=1
C                ENDIF
C              ENDDO !nn
C              IF(NJ_LOC(NJL_FIBR,0,nr).GE.2) THEN
C                NB_IMBRIC=NBJ(NJ_LOC(NJL_FIBR,2,nr),ne)
C                DO nn=1,NNT(NB_IMBRIC)
C                  NVJE(nn,NB_IMBRIC,NJ_LOC(NJL_FIBR,2,nr),ne)=1
C                  IF(NBC(nb_extended).EQ.7) THEN
C                    NVJE(nn,nb_extended,NJ_LOC(NJL_FIBR,2,nr),ne)=1
C                  ENDIF
C                ENDDO !nn
C              ENDIF
C            ENDDO !noelem
C            IF(PROMPT_NV(1).OR.PROMPT_NV(2)) THEN
C              WRITE(OP_STRING,'(''>>Redefine elements to use correct '
C     '          //'versions for fibre (& imbric) field(s)'')')
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDIF
C        ENDIF
C!newe
C      ELSE IF(NAT(NB_FIBRE).GT.0) THEN
C        DO noelem=1,NEELEM(0,nr)
C          ne=NEELEM(noelem,nr)
C          IF(NIT(NBJ(1,ne)).GT.1) THEN
C            WRITE(CHAR1,'(I4)') ne
C            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
C            FORMAT='('' element number '//CHAR1(IBEG1:IEND1)//': '')'
C            CALL GINOUT(IOTYPE,0,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '        ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            nj=NJ_LOC(NJL_FIBR,1,nr)
C            WRITE(CHAR1,'(I1)') nj
C            IF(ne.EQ.1) THEN
C              RDEFLT(1)=0.0D0
C            ELSE IF(ne.GT.1) THEN
C              RDEFLT(1)=XA(1,nj,ne-1)
C            ENDIF
C            WRITE(CHAR4,'(E12.5)') RDEFLT(1)
C            FORMAT='($,'' the Xj('//CHAR1(1:1)//
C     '        ') coordinate is ['//CHAR4(1:12)//']: '',E12.5)'
C            IF(IOTYPE.EQ.3) THEN
C              IF(JTYP13.EQ.1) THEN
C                RDATA(1)=XA(1,nj,ne)*180.0D0/PI
C              ELSE IF(JTYP13.EQ.2) THEN
C                RDATA(1)=XA(1,nj,ne)
C              ENDIF
C            ENDIF
C            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C            IF(IOTYPE.NE.3) THEN
C              IF(JTYP13.EQ.1) THEN
C                XA(1,nj,ne)=RDATA(1)*PI/180.0D0
C              ELSE IF(JTYP13.EQ.2) THEN
C                XA(1,nj,ne)=RDATA(1)
C              ENDIF
C            ENDIF
C          ENDIF
C        ENDDO

      ENDIF

      CALL EXITS('IPFIBR')
      RETURN
 9999 CALL ERRORS('IPFIBR',ERROR)
      DO njj1=1,NJ_LOC(NJL_FIBR,0,nr)
        NJ_LOC(NJL_FIBR,njj1,nr)=0
      ENDDO
      NJ_LOC(NJL_FIBR,0,nr)=0
      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          IF(NJ_LOC(njj1,njj2,nr).GT.NJ_LOC(0,0,nr))
     '      NJ_LOC(0,0,nr)=NJ_LOC(njj1,njj2,nr)
        ENDDO !njj2
      ENDDO !njj1
      NJ_LOC(0,0,0)=0
      NJ_LOC(NJL_FIBR,0,0)=0
      DO nr=1,NRT
        IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0)) NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
        IF(NJ_LOC(NJL_FIBR,0,nr).GT.NJ_LOC(NJL_FIBR,0,0))
     '    NJ_LOC(NJL_FIBR,0,0)=NJ_LOC(NJL_FIBR,0,nr)
      ENDDO !nr
      CALL EXITS('IPFIBR')
      RETURN 1
      END


