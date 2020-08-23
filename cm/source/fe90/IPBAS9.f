      SUBROUTINE IPBAS9(IBT,IDO,INP,NBH,NBJ,NEELEM,NELIST,NHE,NHP,
     '  NKJE,NPF,NP_INTERFACE,NPNE,NPNODE,nr,NVJE,nx,CE,
     '  CURVCORRECT,SE,XA,XE,XP,ERROR,*)

C#### Subroutine: IPBAS9
C###  Description:
C###    IPBAS9 inputs basis functions for dependent variables.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NHE(NEM),
     '  NHP(NPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVJE(NNM,NBFM,NJM,NEM),nx
      REAL*8 CE(NMM,NEM),CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO, !SMAR009 22/12/98 INITIAL,
     '  n,n1,nb,nb_offset,nc,ne,nh,nhx,nhx_MAX,nj,noelem,NOQUES,np
      REAL*8 DXNDS(3,2)
      CHARACTER CHAR1*1,CHAR2*3
      LOGICAL FILEIP,HERMITE_3D,HERMITE_Q,INLIST,INTERFACE,NORMB,NOCROSS

      CALL ENTERS('IPBAS9',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C***  Set NHE for the problem.

      nhx_MAX=0
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
C*** cpb 28/11/94 NHE(ne)=1 temporary
        NHE(ne)=1
        IF(NHE(ne).GT.nhx_MAX) nhx_MAX=NHE(ne)
      ENDDO
      CALL CALC_NH_LOC(nhx_MAX,nx,ERROR,*9999)

      DO nhx=1,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        NORMB=.FALSE.
        IF(call_mesh.and.iotype.NE.3) THEN
          !Mesh has automatically been calculated so set up the
          !initial conditions automatically as well (mesh is for
          !eccentric spheres model or one of the dipole solutions).
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NBH(nh,1,ne)=NBJ(1,ne)
            NBH(nh,2,ne)=NBH(nh,1,ne)
          ENDDO
        ELSE
          WRITE(CHAR1,'(I1)') nhx
          FORMAT='(/'' Do you want to specify basis function type'//
     '      ' for dependent variable '//CHAR1(1:1)//' '''//
     '      '/$,'' and its normal derivative separately [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='Y'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
C news AJP 21-6-94
            FORMAT='('' Normal derivative basis function:'')'
            CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            NORMB=.TRUE.
            WRITE(CHAR1,'(I1)') nhx
C cpb 12/11/97 Adding in the ability to specify the dependent
C variable basis to be the same as the independent variable basis
C        FORMAT='(/'' Do you want to specify basis function type'//
C     '    ' for normal derivative '//CHAR1(1:1)//' '''//
C     '    '/$,'' separately for each element [N]? '',A)'
            FORMAT='(/'' Is the basis function type'//
     '        ' for normal derivative '//CHAR1(1:1)//' '''//
     '        '/$,'' different in each element [N]? '',A)'
            IF(IOTYPE.EQ.3) ADATA(1)='Y'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(ADATA(1).EQ.'Y') THEN
C cpb 19/4/97 Adding list input.
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
CC cpb 3/6/95 Changing the default to be the geometric basis
CC                IF(ne.EQ.1) THEN
CC                  IDEFLT(1)=1
CC                ELSE
CC                  IDEFLT(1)=NBH(nh,2,ne-1)
CC                ENDIF
C                IDEFLT(1)=NBJ(nh,ne)
C                WRITE(CHAR2,'(I3)') ne
C                CALL STRING_TRIM(CHAR2,IBEG,IEND)
C                WRITE(CHAR3,'(I3)') IDEFLT(1)
C                CALL STRING_TRIM(CHAR3,IB,IE)
C                FORMAT='($,'' Element '//CHAR2(IBEG:IEND)//' ['
C     '            //CHAR3(IB:IE)//']:'',I2)'
C                IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,2,ne)
C                CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '            IDEFLT,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
C     '            INFO,ERROR,*9999)
C                IF(iotype.NE.3) THEN
C                  NBH(nh,2,ne)=IDATA(1)
C                  CALL ASSERT(NBC(IDATA(1)).EQ.5.OR.NBC(IDATA(1)).EQ.6,
C     '              '>>Not a boundary element basis',ERROR,*9999)
C                ENDIF
C              ENDDO
              FORMAT='(/'' Is the basis function type'//
     '          ' for normal derivative '//CHAR1(1:1)//' '''//
     '          '/$,'' the same as the geometric basis type [N]? '',A)'
              IF(IOTYPE.EQ.3) ADATA(1)='N'
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(ADATA(1).EQ.'N') THEN
                IF(IOTYPE.NE.3) THEN
                  DO noelem=1,NEELEM(0,nr)
                    NBH(nh,2,ne)=0
                  ENDDO !noelem (ne)
                ENDIF
                noelem=0
 6100           FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
                IF(IOTYPE.EQ.3) THEN
                  noelem=noelem+1
                  IF(noelem.LE.NEELEM(0,nr)) THEN
                    ne=NEELEM(noelem,nr)
                    IDATA(1)=ne
                    NELIST(0)=1
                    NELIST(1)=ne
                  ELSE
                    IDATA(0)=0
                    IDATA(1)=0
                  ENDIF
                ENDIF
 6200           CDATA(1)='ELEMENTS' !for use with group input
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '            0,NET(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
                IF(IDATA(1).NE.0) THEN !not default exit
                  IF(IOTYPE.NE.3) THEN
                    NELIST(0)=IDATA(0)
                    DO n=1,IDATA(0)
                      NELIST(n)=IDATA(n)
                      ne=IDATA(n)
                      IF(.NOT.INLIST(ne,NEELEM(1,nr),NEELEM(0,nr),N1))
     '                  THEN
                        WRITE(OP_STRING,'('' >>Element '',I5,'' is not '
     '                    //'in the current region'')') ne
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        GOTO 6200
                      ENDIF
                    ENDDO !n
                  ENDIF !IOTYPE.NE.3
C                 Define basis for first element in group. Rest of
C                 group filled in at the end.
                  ne=NELIST(1)
                  IDEFLT(1)=NBJ(nh,ne)
                  WRITE(CHAR2,'(I2)') IDEFLT(1)
                  FORMAT='($,'' The basis type number is ['
     '              //CHAR2(1:2)//']: '',I2)'
                  IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,2,ne)
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    NBH(nh,2,ne)=IDATA(1)
                    CALL ASSERT(NBC(IDATA(1)).EQ.5.OR.
     '                NBC(IDATA(1)).EQ.6,
     '                '>>Not a boundary element basis',ERROR,*9999)
C                   Apply to all elements in the group
                    DO n=2,NELIST(0)
                      ne=NELIST(n)
                      NBH(nh,2,ne)=IDATA(1)
                    ENDDO !n
                  ENDIF
                  GOTO 6100 !For more elements
                ENDIF
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(NBH(nh,2,ne).EQ.0) THEN
                    WRITE(ERROR,'('' >>Element '',I5,'' does not have '
     '                //'a basis type set'')') ne
                    GOTO 9999
                  ENDIF
                ENDDO !noelem
              ELSE
                nj=NJ_LOC(NJL_GEOM,nhx,nr)
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  NBH(nh,2,ne)=NBJ(nj,ne)
                ENDDO !noelem (ne)
              ENDIF
            ELSE
              FORMAT='($,'' The basis type number is [1]: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,NEELEM(1,nr))
              IDEFLT(1)=1
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(iotype.NE.3) THEN
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  NBH(nh,2,ne)=IDATA(1)
                ENDDO
                CALL ASSERT(NBC(IDATA(1)).EQ.5.OR.NBC(IDATA(1)).EQ.6,
     '            '>>Not a boundary element basis',ERROR,*9999)
              ENDIF
            ENDIF
          ENDIF
          WRITE(CHAR1,'(I1)') nhx
C cpb 12/11/97 Adding in the ability to specify the dependent
C variable basis to be the same as the independent variable basis
C        FORMAT='(/'' Do you want to specify basis function type'//
C     '    ' for dependent variable '//CHAR1(1:1)//' '''//
C     '    '/$,'' separately for each element [N]? '',A)'
          FORMAT='(/'' Is the basis function type'//
     '      ' for dependent variable '//CHAR1(1:1)//' '''//
     '      '/$,'' different in each element [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='Y'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
C CPB 19/4/97 Adding list input
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
CC cpb 3/6/95 Changing the default to be the geometric basis
CC              IF(ne.EQ.1) THEN
CC                IDEFLT(1)=1
CC              ELSE
CC                IDEFLT(1)=NBH(nh,1,ne-1)
CC              ENDIF
C              IDEFLT(1)=NBJ(nh,ne)
C              WRITE(CHAR2,'(I3)') ne
C              CALL STRING_TRIM(CHAR2,IBEG,IEND)
C              WRITE(CHAR3,'(I3)') IDEFLT(1)
C              CALL STRING_TRIM(CHAR3,IB,IE)
C              FORMAT='($,'' Element '//CHAR2(IBEG:IEND)//' ['
C     '          //CHAR3(IB:IE)//']:'',I2)'
C              IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,ne)
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
C     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C              IF(iotype.NE.3) THEN
C                NBH(nh,1,ne)=IDATA(1)
C                CALL ASSERT(NBC(IDATA(1)).EQ.5.OR.NBC(IDATA(1)).EQ.6,
C     '            '>>Not a boundary element basis',ERROR,*9999)
C              ENDIF
C            ENDDO
            FORMAT='(/'' Is the basis function type'//
     '        ' for dependent variable '//CHAR1(1:1)//' '''//
     '        '/$,'' the same as the geometric basis type [N]? '',A)'
            IF(IOTYPE.EQ.3) ADATA(1)='N'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(ADATA(1).EQ.'N') THEN

C LKC 17-MAR-1999 Adding the ability to offset basis by a certain amount.
C*** Is the basis function an offset of the geometric basis?
              FORMAT='(/'' Is the basis function type'//
     '          ' for dependent variable '//CHAR1(1:1)//' '''//
     '          '/$,'' an offset of the geometric basis type [N]? '',A)'
              IF(IOTYPE.EQ.3) ADATA(1)='N'
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ANO,CDATA,CDEFLT
     '          ,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(ADATA(1).EQ.'N') THEN
                IF(IOTYPE.NE.3) THEN
                  DO noelem=1,NEELEM(0,nr)
                    NBH(nh,1,ne)=0
                  ENDDO !noelem (ne)
                ENDIF
                noelem=0
 6300           FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
                IF(IOTYPE.EQ.3) THEN
                  noelem=noelem+1
                  IF(noelem.LE.NEELEM(0,nr)) THEN
                    ne=NEELEM(noelem,nr)
                    IDATA(1)=ne
                    NELIST(0)=1
                    NELIST(1)=ne
                  ELSE
                    IDATA(0)=0
                    IDATA(1)=0
                  ENDIF
                ENDIF
 6400           CDATA(1)='ELEMENTS' !for use with group input
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '            0,NET(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
                IF(IDATA(1).NE.0) THEN !not default exit
                  IF(IOTYPE.NE.3) THEN
                    NELIST(0)=IDATA(0)
                    DO n=1,IDATA(0)
                      NELIST(n)=IDATA(n)
                      ne=IDATA(n)
                      IF(.NOT.INLIST(ne,NEELEM(1,nr),NEELEM(0,nr),N1))
     '                  THEN
                        WRITE(OP_STRING,'('' >>Element '',I5,'' is not '
     '                    //'in the current region'')') ne
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        GOTO 6400
                      ENDIF
                    ENDDO !n
                  ENDIF !IOTYPE.NE.3
C               Define basis for first element in group. Rest of
C               group filled in at the end.
                  ne=NELIST(1)
                  IDEFLT(1)=NBJ(nh,ne)
                  WRITE(CHAR2,'(I2)') IDEFLT(1)
                  FORMAT='($,'' The basis type number is ['
     '              //CHAR2(1:2)//']: '',I2)'
                  IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,ne)
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    CALL ASSERT(
     '                NBC(IDATA(1)).EQ.5.OR.NBC(IDATA(1)).EQ.6,
     '                '>>Not a boundary element basis',ERROR,*9999)
                    NBH(nh,1,ne)=IDATA(1)
C                 Apply to all elements in the group
                    DO n=2,NELIST(0)
                      ne=NELIST(n)
                      NBH(nh,1,ne)=IDATA(1)
                    ENDDO !n
                  ENDIF
                  GOTO 6300 !For more elements
                ENDIF

                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(NBH(nh,1,ne).EQ.0) THEN
                    WRITE(ERROR,'('' >>Element '',I5,'' does not have '
     '                //'a basis type set'')') ne
                    GOTO 9999
                  ENDIF
                ENDDO !noelem
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
            ENDIF
          ELSE
            FORMAT='($,'' The basis type number is [1]: '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,NEELEM(1,nr))
            IDEFLT(1)=1
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(iotype.NE.3) THEN
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                NBH(nh,1,ne)=IDATA(1)
              ENDDO
              CALL ASSERT(NBC(IDATA(1)).EQ.5.OR.NBC(IDATA(1)).EQ.6,
     '          '>>Not a boundary element basis',ERROR,*9999)
            ENDIF
          ENDIF !basis different for different elements
        ENDIF !end of mesh_call
        IF(iotype.NE.3.AND..NOT.NORMB) THEN

C LKC 2-MAR-2000 new assert
          CALL ASSERT(NCM.GE.2,'>> Increase NCM to 2',
     '      ERROR,*9999)

          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NBH(nh,2,ne)=NBH(nh,1,ne)
          ENDDO
        ENDIF
      ENDDO !nhx

! Check if there are any bicubic hermite elements.  If so determine
! whether cross derivative coefficients are to be included in the
! interpolation of the dependent variable (if so BE matrix system
! likely to be poorly conditioned).

      HERMITE_3D = .FALSE.
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr) !This identifies an elem in the BE region.
        nb=NBASEF(NBH(NH_LOC(1,nx),1,ne),1) ! First dependent variable
        IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
          IF(NIM.EQ.1) THEN
            IF(IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4) THEN
            !Hermite in each direction or hermite simplex or hermite sec
              HERMITE_3D=.TRUE.
              GOTO 100
            ENDIF
          ELSE
            IF((IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2).OR.
     '        (IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4).OR.
     '        (((IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6).AND.
     '        IBT(2,1,nb).EQ.4).AND.IBT(1,2,nb).EQ.2).OR.
     '        (((IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6).AND.
     '        IBT(2,2,nb).EQ.4).AND.IBT(1,1,nb).EQ.2)) THEN
            !Hermite in each direction or hermite simplex or hermite sec
              HERMITE_3D=.TRUE.
              GOTO 100
            ENDIF
          ENDIF
        ENDIF
      ENDDO !noelem
 100  CONTINUE
      noelem=1
      HERMITE_Q=.FALSE.
      DO WHILE(.NOT.HERMITE_Q.AND.noelem.LE.NEELEM(0,nr))
        ne=NEELEM(noelem,nr)
        nb=NBASEF(NBH(NH_LOC(1,nx),2,ne),1)
        IF(NIM.EQ.1) THEN
          IF((IBT(1,1,nb).EQ.2).OR.
     '      (IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4).OR.
     '      ((IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6).AND.
     '      IBT(2,1,nb).EQ.2))
     '      THEN
            HERMITE_Q=.TRUE.
          ELSE
            noelem=noelem+1
          ENDIF
        ELSE
          IF((IBT(1,1,nb).EQ.2).OR.
     '      (IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4).OR.
     '      ((IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6).AND.
     '      IBT(2,1,nb).EQ.2).OR.((IBT(1,2,nb).EQ.5.OR.
     '      IBT(1,2,nb).EQ.6).AND.IBT(2,2,nb).EQ.2))
     '      THEN
            HERMITE_Q=.TRUE.
          ELSE
            noelem=noelem+1
          ENDIF
        ENDIF
      ENDDO
C cpb 20/11/96 Adding curvature correction for BEM.
      IF(HERMITE_Q) THEN
        FORMAT='(/$,'' Do you want to use the curvature corrections '
     '    //'for GQ [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(BEMCURVATURECORRECTION) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            BEMCURVATURECORRECTION=.TRUE.
          ELSE
            BEMCURVATURECORRECTION=.FALSE.
          ENDIF
        ENDIF
      ELSE
        BEMCURVATURECORRECTION=.FALSE.
      ENDIF
      IF(HERMITE_3D) THEN
        DO nc=1,2
          KTYP93(nc,nr)=0
        ENDDO !nc
C KAT 13Jan00: Cross derivatives may already be eliminated in ipbase.
        NOCROSS=.FALSE.
        DO nb=1,NBFT
          IF(NBCD(nb).EQ.1) NOCROSS=.TRUE.
        ENDDO
        IF(.NOT.NOCROSS) THEN
          FORMAT='(/$,'' Do you want to set cross derivatives of'//
     '      ' the dependent variables to zero [Y]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='Y'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') KTYP93(1,nr)=1
          IF(KTYP93(1,nr).EQ.1) THEN !Check normal deriv interpolation
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              ! This identifies an element in the BE region.
              nb=NBASEF(NBH(NH_LOC(1,nx),2,ne),1) ! First dependent variable
              IF(NIM.EQ.1) THEN
                IF(IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4) THEN
                  !Hermite in each direction or hermite simplex or herm sec
                  KTYP93(2,nr)=1
                  GOTO 110
                ENDIF
              ELSE
                IF((IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2).OR.
     '            (IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4).OR.
     '            (((IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6).AND.
     '            IBT(2,1,nb).EQ.4).AND.IBT(1,2,nb).EQ.2).OR.
     '            (((IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6).AND.
     '            IBT(2,2,nb).EQ.4).AND.IBT(1,1,nb).EQ.2)) THEN
                  !Hermite in each direction or hermite simplex or herm sec
                  KTYP93(2,nr)=1
                  GOTO 110
                ENDIF
              ENDIF
            ENDDO
 110        CONTINUE
          ENDIF
        ENDIF !.NOT.NOCROSS
      ELSE
        KTYP93(1,nr)=0
        KTYP93(2,nr)=0
      ENDIF

      IF(BEMCURVATURECORRECTION) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          np=NPNE(1,nb,ne)
          INTERFACE=(NP_INTERFACE(np,0).GT.1).AND.
     '      (NP_INTERFACE(np,1).NE.nr)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '      nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '      ERROR,*9999)
          CALL CALC_CURVCORRECT(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      nb,ne,NPNE(1,1,ne),nr,CE(1,ne),CURVCORRECT(1,1,1,ne),DXNDS,
     '      SE(1,1,ne),XE,XP,INTERFACE,ERROR,*9999)
        ENDDO !noelem (ne)
      ENDIF

C     Set NHP for the problem.
      IF(IOTYPE.NE.3) THEN
        CALL CALC_NHP(NBH,NEELEM,NHE,NHP,NPNE,NPNODE,nr,nx,
     '    ERROR,*9999)
      ENDIF !iotype.ne.3

      CALL EXITS('IPBAS9')
      RETURN
 9999 CALL ERRORS('IPBAS9',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPBAS9')
      RETURN 1
      END


