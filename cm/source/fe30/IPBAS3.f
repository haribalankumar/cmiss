      SUBROUTINE IPBAS3(NBJ,NBH,NEELEM,NELIST,NHE,NHP,NPNE,
     '  NPNODE,nr,NW,nx,ERROR,*)

C#### Subroutine: IPBAS3
C###  Description:
C###    IPBAS3 inputs basis functions for dependent variables for
C###    FE30 problems.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NHE(NEM),NHP(NPM,0:NRM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NW(NEM,3),nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,
     '  INITIAL(1),n,nb_offset,nc,n1,ne,nh,nhx,nhx_MAX,nj,noelem,NOQUES,
     '  nr_loop
      LOGICAL BDRY_ELEMENT,FILEIP,INLIST,NORMB
      CHARACTER CHAR1*3,CHAR2*3

      CALL ENTERS('IPBAS3',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
!news AJP 21/3/96
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        NW(ne,1)=1 !Initialise for elliptic problems
      ENDDO
!newe

! Set NHE for the problem.
      nhx_MAX=0
      DO noelem=1,NEELEM(0,nr) !to set #dependent variables
        ne=NEELEM(noelem,nr)
        IF(ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4) THEN
C         Static or quasi-static analysis
          IF(ITYP2(nr,nx).EQ.7) THEN !Maxwells
            IF(ITYP3(nr,nx).EQ.1) THEN !Electrostatic
              NHE(ne)=1
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN  !Magnetostatic
              NHE(ne)=NJ_LOC(NJL_GEOM,0,nr)
            ELSE IF(ITYP3(nr,nx).EQ.3) THEN  !Magnetostatic
              NHE(ne)=NJ_LOC(NJL_GEOM,0,nr)
            ENDIF            
          ELSE IF(ITYP2(nr,nx).EQ.8) THEN       !Fluid mechanics
            IF(ITYP3(nr,nx).EQ.1) THEN       !Prandtl bdry layer eqtns
              NHE(ne)=NJ_LOC(NJL_GEOM,0,nr)
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN  !Fluid interface stability
              NHE(ne)=1
            ENDIF
          ELSE IF(ITYP2(nr,nx).EQ.9) THEN  !Oxygen transport
            IF(ITYP3(nr,nx).EQ.1) THEN       !Multi-field oxy transport
              NHE(ne)=KTYP15
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN  !Glucose-oxygen transport
              NHE(ne)=1
            ENDIF
          ELSE
            NHE(ne)=1
          ENDIF

        ELSE IF(ITYP5(nr,nx).EQ.2) THEN  !Time integration
          IF(ITYP2(nr,nx).EQ.3) THEN !Advection-diffusion equations
            NHE(ne)=KTYP3A(nx)
          ELSE IF(ITYP2(nr,nx).EQ.5) THEN       !Navier-Stokes equations
            IF(ITYP3(nr,nx).EQ.1) THEN       !Fluid in elastic tube
              NHE(ne)=4
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN  !Lung gas transport
              NHE(ne)=1
            ELSE IF(ITYP3(nr,nx).EQ.3) THEN !General Navier Stokes
              NHE(ne)=NJ_LOC(NJL_GEOM,0,nr) + 1
            ELSE IF(ITYP3(nr,nx).EQ.4) THEN !Stokes flow (no advection)
              NHE(ne)= NJ_LOC(NJL_GEOM,0,nr) + 1
            ELSE IF(ITYP3(nr,nx).EQ.5) THEN !Bidomain Stokes flow
              NHE(ne)= 2*NJ_LOC(NJL_GEOM,0,nr) + 2
            ENDIF
          ELSE IF(ITYP2(nr,nx).EQ.10) THEN !Oxygen transport
            IF(ITYP3(nr,nx).EQ.1) THEN       !Multi-field oxy transport
              NHE(ne)=KTYP15
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN  !Glucose-oxygen transport
              NHE(ne)=1
            ENDIF
          ELSE IF(ITYP2(nr,nx).EQ.11) THEN !pulmonary transport
            IF(ITYP3(nr,nx).EQ.1)THEN
              NHE(ne)=1 !gas concentration
            ELSE IF(ITYP3(nr,nx).EQ.2)THEN
              NHE(ne)=2 !nh=1 for temp; nh=2 for water vapour
            ELSE IF(ITYP3(nr,nx).GE.3)THEN
              NHE(ne)=1 !blood or simple flow
            ELSE
              NHE(ne)=1 !change as necessary
            ENDIF
          ELSE
            NHE(ne)=1
          ENDIF

        ELSE IF(ITYP5(nr,nx).EQ.5) THEN !Wavefront paths
          NHE(ne)=1 !incidence times

        ENDIF
        IF(NHE(ne).GT.nhx_MAX) nhx_MAX=NHE(ne)
      ENDDO !noelem

      CALL CALC_NH_LOC(nhx_MAX,nx,ERROR,*9999)

! Set NBH for the problem.
C LC 25/2/97 section archived : news AJP 25/1/96

!     Check whether any bdry elements
      BDRY_ELEMENT=.FALSE.
      DO nr_loop=1,NRT
        IF(ITYP4(nr_loop,nx).EQ.2) THEN
          BDRY_ELEMENT=.TRUE.
        ENDIF
      ENDDO

C KAT 2001-04-12: moved to IPEQUA
C      CALL ASSERT(NCM.GE.2,'>>Increase NCM to >= 2',ERROR,*9999)

      DO nhx=1,NH_LOC(0,nx)
        nh=nh_loc(nhx,nx)
        WRITE(CHAR1,'(I1)') nhx

        NORMB=.FALSE.
        IF(BDRY_ELEMENT) THEN !BEM in some region
C THE FOLLOWING QUESTON SHOULD BE ASKED IN EVERY CASE ??? AJP 3/2/95
          FORMAT='(/'' Do you want to specify basis function type'//
     '      ' for dependent variable '//CHAR1(1:1)//' '''//
     '      '/$,'' and its normal derivative separately [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='Y'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
            NORMB=.TRUE.
C cpb 19/4/97 Adding list input
C            FORMAT='($,'' The basis type number for the normal'//
C     '        ' derivative is [1]: '',I1)'
C            IF(IOTYPE.EQ.3) THEN
C              IDATA(1)=NBH(nh,2,NEELEM(1,nr))
C              IF(IDATA(1).EQ.0) IDATA(1)=NBH(nh,1,NEELEM(1,nr))
C            ENDIF
C            INITIAL=1
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,INITIAL,1,99,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(iotype.ne.3) THEN
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                NBH(nh,2,ne)=IDATA(1)
C              ENDDO
C            ENDIF
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
                  IDEFLT(1)=NBJ(1,ne)
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
              INITIAL(1)=1
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,INITIAL,
     &          1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(iotype.NE.3) THEN
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  NBH(nh,2,ne)=IDATA(1)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDIF !End of coupled loop question

        WRITE(CHAR1,'(I1)') nhx
C cpb 12/11/97 Adding in the ability to specify the dependent
C variable basis to be the same as the independent variable basis
C        FORMAT='(/'' Do you want to specify basis function type'//
C     '    ' for dependent variable '//CHAR1(1:1)//' '''//
C     '    '/$,'' separately for each element [N]? '',A)'
        FORMAT='(/'' Is the basis function type'//
     '    ' for dependent variable '//CHAR1(1:1)//' '''//
     '    '/$,'' different in each element [N]? '',A)'
        IF(IOTYPE.EQ.3) ADATA(1)='Y'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
C cpb 19/4/97 Adding list input
C          DO noelem=1,NEELEM(0,nr)
C            ne=NEELEM(noelem,nr)
Cc cpb 27/8/95 Changing the default to be the geometric basis
CC              IF(ne.EQ.1) IDEFLT(1)=1
CC              ne1=NEELEM(noelem-1,nr)
CC              IF(ne.ne.1) IDEFLT(1)=NBH(nh,1,ne1)
C            IDEFLT(1)=NBJ(1,ne)
C            WRITE(CHAR1,'(I3)') ne
C            CALL STRING_TRIM(CHAR1,IBEG,IEND)
C            WRITE(CHAR2,'(I3)') IDEFLT(1)
C            CALL STRING_TRIM(CHAR2,IB,IE)
C            FORMAT='($,'' Element '//CHAR1(IBEG:IEND)//' ['
C     '        //CHAR2(IB:IE)//']: '',I1)'
C            IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,ne)
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(iotype.ne.3) THEN
C              NBH(nh,1,ne)=IDATA(1) !variables
C              NBH(nh,2,ne)=IDATA(1) !equations
C              IF((KTYP90.EQ.2).AND.(.NOT.NORMB))THEN !Coupled FE/BE
C                NBH(nh,2,ne)=NBH(nh,1,ne)
C              ENDIF
C            ENDIF
C          ENDDO !noelem
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
C                 Define basis for first element in group. Rest of
C                 group filled in at the end.
                  ne=NELIST(1)
                  IDEFLT(1)=NBJ(1,ne)
                  WRITE(CHAR2,'(I2)') IDEFLT(1)
                  FORMAT='($,'' The basis type number is ['
     '              //CHAR2(1:2)//']: '',I2)'
                  IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,ne)
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    NBH(nh,1,ne)=IDATA(1)
C                   Apply to all elements in the group
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

              ENDIF !basis is an offset


C CHANGES END

          ELSE
C*** Setting basis to be the same as geometric
            nj=NJ_LOC(NJL_GEOM,nhx,nr)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NBH(nh,1,ne)=NBJ(nj,ne)
            ENDDO !noelem (ne)
          ENDIF !geomtric basis

          IF(IOTYPE.NE.3.AND..NOT.NORMB.AND.NCT(nr,nx).GT.1) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NBH(nh,2,ne)=NBH(nh,1,ne)
            ENDDO
          ENDIF
        ELSE !AJP 25/1/96
          FORMAT='($,'' The basis type number is [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,NEELEM(1,nr))
          INITIAL(1)=NBJ(1,NEELEM(1,nr))
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,INITIAL,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
C??? KAT 2007-01-30: Should there be some sort of consistency check
C??? here, for number of xi directions?  Number of Guass points?
C??? Anything else?
              DO nc=1,NCT(nr,nx)
                NBH(nh,nc,ne)=IDATA(1)
              ENDDO
C KAT 2001-04-12: already done
C              IF((KTYP90.EQ.2).AND.(.NOT.NORMB))THEN !Coupled FE/BE
C                NBH(nh,2,ne)=NBH(nh,1,ne)
C              ENDIF
            ENDDO !noelem
          ENDIF !iotype
        ENDIF !adata(1)
      ENDDO !nhx

C 25/2/97 LC removed section from :
C     cpb 11/9/95 Adding zero cross derivative basis functions.

C     Set NHP for the problem.
      IF(IOTYPE.NE.3) THEN
        IF(ITYP2(nr,nx).EQ.11)THEN
          !special case for large 1D lung trees & 1D capillaries
          CALL CALC_NHP_1D(NHE(NEELEM(1,nr)),NHP,NPNODE,nr,ERROR,*9999)
        ELSE
          CALL CALC_NHP(NBH,NEELEM,NHE,NHP,NPNE,NPNODE,nr,nx,
     '      ERROR,*9999)
        ENDIF !ITYP3
      ENDIF !iotype.ne.3

      CALL EXITS('IPBAS3')
      RETURN
 9999 CALL ERRORS('IPBAS3',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPBAS3')
      RETURN 1
      END



