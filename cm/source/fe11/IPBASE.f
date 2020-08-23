      SUBROUTINE IPBASE(IBT,IDO,INP,NAN,NBB,NDET,NGAP,NKB,NNB,NSB,
     '  DET,PG,WG,XIG,ERROR,*)

C#### Subroutine: IPBASE
C###  Description:
C###    IPBASE inputs basis function parameters.

C**** For most boundary element analysis the user specifies a low order
C**** and high order quadrature scheme. These are stored in
C**** NGLIMITS(ni,nb,1) and NGLIMITS(ni,nb,2).
C**** Several other quadrature schemes are
C**** defined  with orders varying between the low and high order
C**** schemes.
C**** ALIM(ni,nb) and BLIM(ni,nb) give limits on integrals for which the
C**** quadrature schemes are reqd. e.g. ALIM(1,nb)=0 and BLIM(1,nb)=1
C**** for normal quadrature; ALIM(1,nb)=0, BLIM(1,nb)=1/2 for quadrature
C**** from 0 to 1/2 (required when e.g. a 1D quad element needs to be
C**** split in two).

C#### Variable: IDO(nk,nn,0:ni,nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    IDO(nk,nn,ni,nbf) is an index for derivative order:
C###    1-zeroth order, 2-first  order. Thus IDO(nk,1..,nbf)=1,2,1
C###    implies that derivative  number nk is a first derivative wrt
C###    Xi(2). IDO(nk, nn,0,nbf) is nu (partial derivative) number
C###    corresponding to nk.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBB,NDET(NBFM,0:NNM),NGAP(NIM,NBM),
     '  NKB(2,2,2,NNM,NBFM),NNB(4,4,4,NBFM),NSB(NKM,NNM,NBFM)
      REAL*8 DET(NBFM,0:NNM,NGM,6),PG(NSM,NUM,NGM,NBM),
     '  WG(NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ICHAR,IDOI(3),ido1,ido2,ido3,IEND,INFO,INPI(3),inp1,
     '  inp2,inp3,INPMC,nb,nbf,nbf1,nbf2,add_nbf,ng,ni,nif,nix,NITB,nk,
     '  nn,NOQUES,ns
      CHARACTER CHAR1*10
      LOGICAL CONTINUE,FILEIP

      CALL ENTERS('IPBASE',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      add_nbf=0

      IF(NBB.EQ.0) THEN
        FORMAT='($,'' Enter the number of types of basis'//
     '    ' function [1]: '',I2)'
        IF(IOTYPE.NE.3) THEN
          IF(.NOT.ADD) THEN
            NBFT=0
            NBT=0
          ENDIF
          nbf1=NBFT+1
          nb=NBT
        ELSE !IF(IOTYPE.EQ.3) THEN
          nbf1=1
          nb=0
          IDATA(1)=NBFT
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
           NBT=NBT+IDATA(1)
           NBFT=NBFT+IDATA(1)
        ENDIF
        CALL ASSERT(NBFT.LE.NBFM,'>>Increase NBFM',ERROR,*9999)
        CALL ASSERT(NBT.LE.NBM,'>>Increase NBM',ERROR,*9999)
        nbf2=NBFT
      ELSE
        nbf1=NBB
        nbf2=NBB
        nb=NBASEF(nbf1,1)-1
      ENDIF

      DO nbf=nbf1,nbf2
        IF(IOTYPE.EQ.3) THEN
          nb=NBASEF(nbf1,1)
        ELSE !IF(IOTYPE.NE.3) THEN
          nb=nb+1
          NAT(nbf)=0
          NBCD(nbf)=0
          NKT(0,nbf)=0
          NNT(nbf)=0
          NST(nbf)=0
          NBASEF(nbf,0)=1
          NBASEF(nbf,1)=nb
          DO ni=1,NIM
            IBT(1,ni,nbf)=0
            IBT(2,ni,nbf)=0
            DO nn=1,NNM
              INP(nn,ni,nbf)=1
              DO nk=1,NKM
                IDO(nk,nn,ni,nbf)=1
              ENDDO
            ENDDO
          ENDDO
          DO inp3=1,4
            DO inp2=1,4
              DO inp1=1,4
                NNB(inp1,inp2,inp3,nbf)=0
              ENDDO
            ENDDO
          ENDDO
          DO nn=1,NNM
            NKT(nn,nbf)=0
            DO nk=1,NKM
              IDO(nk,nn,0,nbf)=0
            ENDDO
            DO ido3=1,2
              DO ido2=1,2
                DO ido1=1,2
                  NKB(ido1,ido2,ido3,nn,nbf)=0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        IF(.NOT.ADD) THEN
          WRITE(CHAR1,'(I2)') nbf
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
        ELSE
          add_nbf=add_nbf+1
          WRITE(CHAR1,'(I2)') add_nbf
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
        ENDIF

        FORMAT='(/'' For basis function type '//CHAR1(IBEG:IEND)//
     '         ' the type of nodal interpolation is [1]:'''//
     '      '/''   (0) Auxiliary basis only'''//
     '      '/''   (1) Lagrange/Hermite tensor prod'''//
     '      '/''   (2) Simplex/Serendipity/Sector'''//
     '      '/''   (3) B-spline tensor product'''//
     '      '/''   (4) Fourier Series/Lagrange/Hermite tensor prod'''//
     '      '/''   (5) Boundary Element Lagrange/Hermite tensor pr.'''//
     '      '/''   (6) Boundary Element Simplex/Serendipity/Sector'''//
     '      '/''   (7) Extended Lagrange (multigrid collocation)'''//
     '      '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NBC(nbf)
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,7,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NBC(nbf)=IDATA(1)

        IF(NBC(nbf).EQ.1) THEN      !Lagrange/Hermite tensor prod basis
          CALL BASIS1(IBT(1,1,nbf),IDO(1,1,0,nbf),INP(1,1,nbf),nbf,
     '      NGAP,PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)

        ELSE IF(NBC(nbf).EQ.2) THEN !Simplex/Serendipity/Sector basis
          CALL BASIS2(IBT(1,1,nbf),IDO(1,1,0,nbf),INP(1,1,nbf),nbf,
     '      NGAP,PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)

        ELSE IF(NBC(nbf).EQ.3) THEN !B-spline tensor product basis
          CALL BASIS3(IBT(1,1,nbf),INP(1,1,nbf),nbf,NGAP,ERROR,*9999)

        ELSE IF(NBC(nbf).EQ.4) THEN !Fourier Series tensor prod basis
          CALL BASIS4(IBT(1,1,nbf),IDO(1,1,0,nbf),INP(1,1,nbf),nbf,
     '      NGAP,PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)

        ELSE IF(NBC(nbf).EQ.5) THEN !Boundary element tensor prod basis
          CALL ASSERT(USE_BEM.EQ.1,
     '      '>> Must have USE_BEM set',ERROR,*9999)
          CALL BASIS5(IBT,IDO,INP,nbf,NDET,NGAP,
     '      DET,PG,WG,XIG,ERROR,*9999)

        ELSE IF(NBC(nbf).EQ.6) THEN !Boundary simplex/sector basis
          CALL ASSERT(USE_BEM.EQ.1,
     '      '>> Must have USE_BEM set',ERROR,*9999)
          CALL BASIS6(IBT,IDO,INP,nbf,NDET,NGAP,
     '      DET,PG,WG,XIG,ERROR,*9999)

        ELSE IF(NBC(nbf).EQ.7) THEN !Extended Lagrange basis
          CALL BASIS7(IBT(1,1,nbf),IDO(1,1,0,nbf),INP(1,1,nbf),nbf,
     '      NGAP,PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)

        ENDIF

C KAT   Invert INP,IDO to obtain NNB,NKB and calculate NSB
        ns=1
        DO ni=1,3
          INPI(ni)=1
          IDOI(ni)=1
        ENDDO
        NITB=MIN(NIT(nbf),3)
        DO nn=1,NNT(nbf)
          DO ni=1,NITB
            INPI(ni)=INP(nn,ni,nbf)
          ENDDO
          NNB(INPI(1),INPI(2),INPI(3),nbf)=nn
          DO nk=1,NKT(nn,nbf)
            DO ni=1,NITB
              IDOI(ni)=IDO(nk,nn,ni,nbf)
            ENDDO
            NKB(IDOI(1),IDOI(2),IDOI(3),nn,nbf)=nk
            NSB(nk,nn,nbf)=ns
            ns=ns+1
          ENDDO
        ENDDO
        IF(NBC(nbf).EQ.2.OR.NBC(nbf).EQ.6) THEN !could be simplex/sector
C         Set up NNB for collapsed nodes
          IF(NIM.GT.1) THEN
C           This budget code was essentially moved from NENXI fe02 v1.424
            IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.3) THEN
C             Hermite simplex
              IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                NNB(2,1,1,nbf)=NNB(1,1,1,nbf)
              ELSE !Apex at node 3
                NNB(2,2,1,nbf)=NNB(1,2,1,nbf)
              ENDIF
            ENDIF
          ENDIF !NIM
          DO ni=1,NITB !dirn of collapse
            IF(IBT(1,ni,nbf).EQ.5.OR.IBT(1,ni,nbf).EQ.6) THEN !collapsed sector
              nix=IBT(3,ni,nbf) !out of face
C             Find node position index of collapsed face
              IF(IBT(1,ni,nbf).EQ.5) THEN !collapsed at xi=0
                INPI(nix)=1
              ELSE !IBT(1,ni,nb).EQ.6 collapsed at xi=1
                IF(IBT(1,nix,nbf).EQ.1) THEN !Lagrange
                  INPI(nix)=IBT(2,nix,nbf)+1
                ELSE !IBT(1,nix,nb).EQ.2 Hermite
                  INPI(nix)=2
                ENDIF !IBT(1,nix,nb)
              ENDIF !IBT(1,ni,nb)
C             Find highest position index in direction of collapse
              IF(IBT(2,ni,nbf).EQ.4) THEN !Hermite
                INPMC=2
              ELSE !Lagrange
                INPMC=IBT(2,ni,nbf)+1
              ENDIF !IBT(2,ni,nb)
C             Copy value of NNB for INPI(ni)=1 to other vertices in
C             collapsed face.
              nif=6-ni-nix !direction in face perpendicular to collapse
              INPI(nif)=1
              CONTINUE=.TRUE.
              DO WHILE(CONTINUE) !perpen to collapse
                INPI(ni)=1 !in collapsed dirn
                nn=NNB(INPI(1),INPI(2),INPI(3),nbf)
                IF(nn.EQ.0) THEN
                  CONTINUE=.FALSE.
                ELSE
                  DO WHILE(CONTINUE.AND.INPI(ni).LT.INPMC) !collapse dirn
                    INPI(ni)=INPI(ni)+1
                    NNB(INPI(1),INPI(2),INPI(3),nbf)=nn
                  ENDDO !INPI(ni)
                  INPI(nif)=INPI(nif)+1
                  CONTINUE=INPI(nif).LE.4
                ENDIF !nn
              ENDDO !INPI(nif)
            ENDIF !collapsed sector
          ENDDO !ni
        ENDIF !NBC(nb).EQ.2

        IF(NBC(nbf).LE.2) THEN !Add auxillary basis
          FORMAT='($,'' Enter the number of auxiliary element'//
     '      ' parameters [0]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NAT(nbf)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,11,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NAT(nbf)=IDATA(1)

          CALL ASSERT(NAM.GE.NAT(nbf),'>>NAM too small',ERROR,*9999)
          CALL ASSERT(NSM.GE.NAT(nb)+NST(nb),'>>NSM too small',ERROR,
     '      *9999)
          IF(NAT(nbf).GT.0) THEN
            FORMAT='(/'' Auxiliary basis is [1]:'''//
     '      '/''   (1) Legendre'''//
     '      '/''   (2) Fourier'''//
     '      '/''   (3) Pressure'''//
     '      '/''   (4) Unused'''//
     '      '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=NABTYP(nbf)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     &        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NABTYP(nbf)=IDATA(1)
            CALL BASIS8(NBC(nbf),NAN(1,1,nbf),nbf,NGAP,PG(1,1,1,nb),
     '        WG(1,nb),XIG(1,1,nb),ERROR,*9999)
          ELSE
            NABTYP(nbf)=0
          ENDIF

C MLB 18 August 1997
C Multigrid levels now set up at solve time under new grid scheme
C
C        ELSE IF(NBC(nb).EQ.7) THEN !Extended Lagrange (mg collocation)
C          FORMAT='($,'' Enter the number of multigrid levels [1]: '','
C     '      //'I2)'
C          IF(IOTYPE.EQ.3) IDATA(1)=NMGT
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,11,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) NMGT=IDATA(1)
C          CALL ASSERT(NAM.GE.NMGT,'>>NAM too small',ERROR,*9999)
C          NABTYP(nb)=0
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,'('' IBT(1,ni,'',I2,''): '',3I3)')
     '      nbf,(IBT(1,ni,nbf),ni=1,NIT(nbf))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' IBT(2,ni,'',I2,''): '',3I3)')
     '      nbf,(IBT(2,ni,nbf),ni=1,NIT(nbf))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO ng=1,NGT(nb)
            WRITE(OP_STRING,'('' XIG(ni,'',I3,'','',I2,'
     '         //'''): '',3E11.3)') ng,nb,(XIG(ni,ng,nb),NI=1,NIT(nbf))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

        IF(NKT(0,nbf).GT.1.AND.NNT(nbf).GT.0) THEN
        IF(.NOT.ADD) THEN
          WRITE(CHAR1,'(I2)') nbf
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
        ELSE
          WRITE(CHAR1,'(I2)') add_nbf
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
        ENDIF
C cpb 5/12/96 Re-adding arc length scale factors
c cpb 8/10/94 Changing around the options to specify scale factors for
C either C1 or G1 continuous interpolation
C          FORMAT='(/'' For basis function type '//
C     '      CHAR1(IBEG:IEND)//' specify option [4]: '''//
C     '      '/''   (1) Unit scale factors'''//
C     '      '/''   (2) Element scale factors read in'''//
C     '      '/''   (3) Global scale factors read in'''//
C     '      '/''   (4) Scale factors calc.d from arc length'''//
C     '      '/''   (5) Scale factors calc.d from angle change'''//
C     '      '/$,''    '',I1)'
C          IDEFLT(1)=4
C KAT 28Aug98: Adding harmonic mean and changing format
C          FORMAT='(/'' For basis function type '//
C     '      CHAR1(IBEG:IEND)//' specify option [5]: '''//
C     '      '/''   (1) Unit scale factors'''//
C     '      '/''   (2) Element scale factors read in'''//
C     '      '/''   (3) Global scale factors read in'''//
C     '      '/''   (4) Scale factors calc.d from angle change'''//
C     '      '/''   (5) Scale factors calc.d from average arc length'''//
C     '      '/''   (6) Scale factors calc.d from arc length'''//
C     '      '/$,''    '',I1)'
          FORMAT='(/'' For basis function type '//
     '      CHAR1(IBEG:IEND)//' scale factors are [6]: '''//
     '      '/''   (1) Unit'''//
     '      '/''   (2) Read in - Element based'''//
     '      '/''   (3) Read in - Node based'''//
     '      '/''   (4) Calculated from angle change'''//
     '      '/''   (5) Calculated from arc length'''//
     '      '/''   (6) Calculated from arithmetic mean arc length'''//
     '      '/''   (7) Calculated from harmonic mean arc length'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=6
          IF(IOTYPE.EQ.3) IDATA(1)=NBI(nbf)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3)NBI(nbf)=IDATA(1)
        ELSE
          NBI(nbf)=1
        ENDIF

        IF(NBC(nbf).EQ.5.OR.NBC(nbf).EQ.6) THEN !bem or singular basis
          IF(IOTYPE.NE.3)THEN
            NBT=NBT+NBASEF(nbf,0)-1
          ENDIF
        ENDIF

      ENDDO

      DO nbf=1,NBFT
        CALL ASSERT(NNT(nbf).LE.64,'>>NNT exceeds dimension of'
     '   //' NKT(nn,nbf) in GEOM00.CMN',ERROR,*9999)
      ENDDO !nb

      CALL EXITS('IPBASE')
      RETURN
 9999 CALL ERRORS('IPBASE',ERROR)
      CALL EXITS('IPBASE')
      RETURN 1
      END


