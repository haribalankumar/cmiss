      SUBROUTINE IPBAS5(NBH,NBJ,NEELEM,NHE,NHP,NPNE,NPNODE,nr,nx,
     '  ERROR,*)

C#### Subroutine: IPBAS5
C###  Description:
C###    IPBAS5 inputs basis functions for dependent variables of finite
C###    elasticity problems.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM),NHP(NPM,0:NRM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ib,IBEG,ICHAR,ie,IEND,INFO,nb_offset,
     '  ne,ne1,nh,nhx,nhx_MAX,nj,noelem,NOQUES
      LOGICAL FILEIP
      CHARACTER CHAR1*3,CHAR2*3
C this data statement is now in DENODS, FE24. When reading deformed
C nodal coords without defining equation this array needs to be set up
C      DATA NHT50/2,2,3, 2,2,3, 3,4,4, 3,3,4, 3,3,4, 3,3,4/
C MLB 21/3/03 Moved initialization to BLK50

      CALL ENTERS('IPBAS5',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=0 !temporary (needs deleting later)

      FORMAT='('' Specify whether the material is [2]: '''//
     '  '/''   (1) Compressible  '''//
     '  '/''   (2) Incompressible'''//
     '  '/''   (3) Incomp + fluid'''//
     '  '/''   (4) Compressible + fluid'''//
     '  '/''   (5) Incomp + inextensible'''//
     '  '/''   (6) Compressible + fluid for lung'''//
     '  '/''   (7) Compressible + face pressure bcs'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=2
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP52(nr)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP52(nr)=IDATA(1)

      IF(IOTYPE.NE.3) THEN
C       Fill in NHT50: Other NHT50 values are const
C       (defined in data statement in DENODS, FE24 see comment above)
        NHT50(1,1)=NJ_LOC(NJL_GEOM,0,nr)
        NHT50(2,1)=NJ_LOC(NJL_GEOM,0,nr)
        NHT50(1,4)=NJ_LOC(NJL_GEOM,0,nr)+1
        NHT50(2,4)=NJ_LOC(NJL_GEOM,0,nr)
        NHT50(1,5)=NJ_LOC(NJL_GEOM,0,nr)
        NHT50(2,5)=NJ_LOC(NJL_GEOM,0,nr)
      ENDIF

      IF(IOTYPE.NE.3) THEN
C       set up # dep vars per element (NHE)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          NHE(ne)=NHT50(KTYP52(nr),KTYP51(nr))
        ENDDO !noelem (ne)
      ENDIF !iotype.ne.3

C     Determine max# nhx's and set up NH_LOC for the current nx
      nhx_MAX=0
      DO noelem=1,NEELEM(0,nr) !to set #dependent variables
        ne=NEELEM(noelem,nr)
        IF(NHE(ne).GT.nhx_MAX) nhx_MAX=NHE(ne)
      ENDDO !noelem
      CALL CALC_NH_LOC(nhx_MAX,nx,ERROR,*9999)

C     Set NBH for the problem.
      DO nhx=1,NHT50(KTYP52(nr),KTYP51(nr))
        nh=NH_LOC(nhx,nx)
        WRITE(CHAR1,'(I1)') nh
C cpb 12/11/97 Adding in the ability to specify the dependent
C variable basis to be the same as the independent variable basis
C        FORMAT='('' Do you want to specify the basis function type'
C     '    //' for dependent variable '//CHAR1(1:1)//''''
C     '    //'/$,'' separately for each element [N]? '',A)'
        FORMAT='(/'' Is the basis function type'//
     '    ' for dependent variable '//CHAR1(1:1)//' '''//
     '    '/$,'' different in each element [N]? '',A)'
        IF(IOTYPE.EQ.3) ADATA(1)='Y'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          FORMAT='(/'' Is the basis function type'//
     '      ' for dependent variable '//CHAR1(1:1)//' '''//
     '      '/$,'' the same as the geometric basis type [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'N') THEN

C LKC 17-MAR-1999 Adding the ability to offset basis by a certain amount.
C*** Is the basis function an offset of the geometric basis?
            FORMAT='(/'' Is the basis function type'//
     '        ' for dependent variable '//CHAR1(1:1)//' '''//
     '        '/$,'' an offset of the geometric basis type [N]? '',A)'
            IF(IOTYPE.EQ.3) ADATA(1)='N'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ANO,CDATA,CDEFLT
     '        ,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(ADATA(1).EQ.'N') THEN

              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(noelem.EQ.1) THEN
                  IDEFLT(1)=0
                ELSE
                  ne1=NEELEM(noelem-1,nr)
                  IDEFLT(1)=NBH(nh,1,ne1)
                ENDIF
                WRITE(CHAR2,'(I3)') IDEFLT(1)
                CALL STRING_TRIM(CHAR2,ib,ie)
                WRITE(CHAR1,'(I3)') ne
                CALL STRING_TRIM(CHAR1,IBEG,IEND)
                FORMAT='($,'' Element '//CHAR1(IBEG:IEND)//' ['
     '            //CHAR2(ib:ie)//']: '',I1)'
                IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,ne)
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,
     '            ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,20,
     '            LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) NBH(nh,1,ne)=IDATA(1)
              ENDDO !noelem (ne)

            ELSE

C LKC 8-APR-1999 initialise
                IF(IOTYPE.EQ.3) nb_offset=-1

C*** Asking question for the offset
C*** Check that basis numbers of defined for these offsets
                FORMAT='($,'' The basis offset number is [1]: '',I2)'
                IF(IOTYPE.EQ.3) IDATA(1)=nb_offset
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,-NBM,NBM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  nb_offset=IDATA(1)
                ENDIF

                CALL ASSERT(ABS(nb_offset).LT.NBT,
     '            '>> Not enough basis functions for this offset',
     '            ERROR,*9999)

C LKC 19-APR-1999 added a new assert
C Do a quick check for the first nj of the first element of the nr
                CALL ASSERT(
     '            NBASEF(NBJ(1,NEELEM(1,nr))+nb_offset,0).GT.0,
     '            '>> No basis function exists for this offset',
     '            ERROR,*9999)

                nj=NJ_LOC(NJL_GEOM,nhx,nr)
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  NBH(nh,1,ne)=NBJ(nj,ne)+nb_offset
                ENDDO !noelem (ne)
              ENDIF ! basis function is an offset

          ELSE
C*** Setting basis to be the same as geometric basis
            nj=NJ_LOC(NJL_GEOM,nhx,nr)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NBH(nh,1,ne)=NBJ(nj,ne)
            ENDDO !noelem (ne)
          ENDIF ! basis same as geometric

        ELSE !AJP 25/1/96
C old       ELSE IF(ADATA(1).EQ.'N') THEN
          FORMAT='($,'' The basis type number is [1]: '',I1)'
          IF(IOTYPE.EQ.3) THEN
            ne=NEELEM(1,nr) !const so use first elem
            IDATA(1)=NBH(nh,1,ne)
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBFM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
C news HS,MPN 21-April-2004
C Just checking if enough basis function types are defined
            CALL ASSERT(IDATA(1).LE.NBT,
     &         '>> Basis funtion type not defined',
     &         ERROR,*9999)
C newe
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NBH(nh,1,ne)=IDATA(1)
            ENDDO !noelem (ne)
          ENDIF
        ENDIF

        IF(IOTYPE.NE.3) THEN
C         copy NBH for nc=2 from nc=1. Force dep var basis is assumed
C         to be the same as displacement dep var basis.
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NBH(nh,2,ne)=NBH(nh,1,ne)
          ENDDO !noelem (ne)
        ENDIF
      ENDDO !nh

C     Set NHP for the problem.
      IF(IOTYPE.NE.3) THEN
        CALL CALC_NHP(NBH,NEELEM,NHE,NHP,NPNE,NPNODE,nr,nx,
     '    ERROR,*9999)
      ENDIF !iotype.ne.3

      CALL EXITS('IPBAS5')
      RETURN
 9999 CALL ERRORS('IPBAS5',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPBAS5')
      RETURN 1
      END


