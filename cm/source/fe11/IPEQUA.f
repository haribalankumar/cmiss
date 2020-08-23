      SUBROUTINE IPEQUA(IBT,IDO,INP,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,LGE,
     '  NBH,NBHF,NBJ,NEELEM,NELIST,NENP,NFF,NHE,NHP,NHQ,
     '  NKH,NKHE,NKJE,NLL,NNF,NPF,NP_INTERFACE,NPL,NPNE,NPNODE,NPNY,
     '  NQNY,nr,NRLIST,NVHE,NVHP,NVJE,NVJP,NW,nx,NXI,NYNE,NYNP,NYNQ,
     '  NYNR,NYQNR,CE,CURVCORRECT,SE,XA,XAB,XE,XP,FIX,ERROR,*)

C#### Subroutine: IPEQUA
C###  Description:
C###    IPEQUA defines equation for region nr and problem nx.

C#### Variable: ITYP1(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP1(nr,nx) is 3,4,5,6 or 9 for use of FE30,FE40,FE50,FE60
C###    or FE90.

C#### Variable: ITYP2(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP2(nr,nx) is 1..15 for equation type.

C#### Variable: ITYP3(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP3(nr,nx) is equation type qualifier (=ktyp51(nr)).

C#### Variable: ITYP4(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP4(nr,nx) is 1..7 for Galerkin FEM/Direct BEM/Indirect BEM/
C###    collocation-multigrid/Finite volume/Grid-based Finite Element/
C###    Grid-based Finite Volume.

C#### Variable: ITYP5(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP5(nr,nx) is 1..6 for static/time/modal/Quasi/wavefront/
C###    buckling.

C#### Variable: ITYP6(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP6(nr,nx) is 1,2 for linear/nonlinear problem.

C#### Variable: ITYP7(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP7(nr,nx) is nonlinear equation/bcs.

C#### Variable: ITYP8(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP8(nr,nx) is radiation bc.

C#### Variable: ITYP12(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###  ITYP12(nr,nx) specifies type of Navier-Stokes flow in elastic tube
C###  solution; = 1 for coronary (or other) types, = 2 for pulmonary
C###  blood flow, =3 for pulmonary flow using external pressure field

C#### Variable: NKHE(nk,nn,nh,ne)
C###  Type: INTEGER
C###  Set_up: IPEQUA,IPFIT
C###  Description:
C###    NKHE(nk,nn,nh,ne) is the global derivative number of local
C###    derivative nk of node nn of element ne for dependent variable nh.

C#### Variable: KTYP3A
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP3A is the number of dependent variables in the
C###    time-dependent advection-diffusion problem type.

C#### Variable: KTYP3C
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP3C defines whether the equations for the time-dependent
C###    advection-diffusion equations are coupled.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'nonl00.cmn'      
C      INCLUDE 'cmiss$reference:mesh00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'stab00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NFF(6,NEM),NHE(NEM),
     '  NHP(NPM,0:NRM),NHQ(NRM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  NQNY(2,NYQM,0:NRCM),nr,NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NW(NEM,3),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNQ(NHM,NQM,0:NRCM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      REAL*8 CE(NMM,NEM),CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,ie,IEND,IFROMC,INFO,LAST_TYPE,MAXOPT,
     '  n,n1,n3co,nb,nc,ne,nh,nhx,njj,nk,nn,noelem,NOQUES,nonode,nonrr,
     '  np,nrr,NT_versions(9),elemtype,IBEG2,IEND2
      INTEGER*4 WORK_PTR,WORK_ISC_PTR,WORK_ISR_PTR
      REAL*8 RFROMC
      CHARACTER CHAR*1,CHAR1*4,CHAR2*20,CHAR_2*2
      LOGICAL BEMREGION,CALC_SPARSITY,CBBREV,FILEIP,INLIST
!     External functions
      INTEGER IDIGITS
      
      CALL ENTERS('IPEQUA',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      KTYP3C(nx)=0

      IF(CBBREV(CO,'COMPLIANCE',4,noco+1,NTCO,n3co)) THEN
        IF(CBBREV(CO,'LAMBERT',4,noco+1,NTCO,n3co)) THEN
          COMPLIANCE=1
        ELSEIF(CBBREV(CO,'NONE',4,noco+1,NTCO,n3co)) THEN
          COMPLIANCE=2
        ELSEIF(CBBREV(CO,'VESSEL',6,noco+1,NTCO,n3co)) THEN
          COMPLIANCE=3
        ELSEIF(CBBREV(CO,'MECHANICS',4,noco+1,NTCO,n3co)) THEN
          COMPLIANCE=4
        ELSEIF(CBBREV(CO,'TRANSPULMONARY',4,noco+1,NTCO,n3co)) THEN
          COMPLIANCE=5
        ELSEIF(CBBREV(CO,'VESSEL_MECHANICS',8,noco+1,NTCO,n3co)) THEN
          COMPLIANCE=6
        ELSE
          WRITE(OP_STRING,'('' Not implemented '')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE
        COMPLIANCE=1   !Default to lambert model
      ENDIF !Compliance

C     Scaling factor for turbulence correction factor - used in CALC_RESIS_FLOW
      IF(CBBREV(CO,'SCALE_TURBULENCE',4,noco+1,NTCO,n3co)) THEN
        SCALE_TURBULENCE=RFROMC(CO(N3CO+1)) 
      ELSE
        SCALE_TURBULENCE=1.0d0  !No scaling
      ENDIF !SCALE_TURBULENCE
      IF(CBBREV(CO,'LPM_FIELD',5,noco+1,NTCO,N3CO))THEN
        nej_cap=IFROMC(CO(N3CO+1))
      ENDIF 
      IF(CBBREV(CO,'CAP_VALUE',5,noco+1,NTCO,N3CO))THEN
        CAP_FIELD=RFROMC(CO(N3CO+1))
      ELSE
        CAP_FIELD=1.d0
      ENDIF            
   
C*** 23/04/08 JHC Removed writing out scale_turbulence to screen
C      write(*,*) 'SCALE_TURBULENCE=',SCALE_TURBULENCE

      IDEFLT(1)=0 !AJP 27-3-95
      IF(nr.GT.1) THEN
        IDEFLT(1)=ITYP5(nr-1,nx)
      ENDIF
      IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
      WRITE(CHAR,'(I1)') IDEFLT(1)
      FORMAT='('' Specify whether ['//CHAR//']:'''//
     '  '/''   (1) Static analysis  '''//
     '  '/''   (2) Time integration '''//
     '  '/''   (3) Modal analysis   '''//
     '  '/''   (4) Quasi-static analysis'''//
     '  '/''   (5) Wavefront path analysis'''//
     '  '/''   (6) Buckling analysis'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=ITYP5(nr,nx)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,6,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) ITYP5(nr,nx)=IDATA(1)

      NOQUES=0
      IF(nr.GT.1) THEN
        IDEFLT(1)=ITYP2(nr-1,nx)
      ENDIF
      IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
      WRITE(CHAR,'(I1)') IDEFLT(1)

      IF(ITYP5(nr,nx).EQ.1) THEN                    !Static analysis
        FORMAT='('' Specify equation ['//CHAR//']:'''//
     '    '/''   (1) Linear elasticity           '''//
     '    '/''   (2) Finite elasticity           '''//
     '    '/''   (3) Laplace equation            '''//
     '    '/''   (4) Helmholtz equation          '''//
     '    '/''   (5) div(k.grad(u))=f            '''//
     '    '/''   (6) Linear 2nd order elliptic   '''//
     '    '/''   (7) Maxwells equations          '''//
     '    '/''   (8) Fluid mechanics             '''//
     '    '/''   (9) Oxygen transport            '''//
     '    '/$,''    '',I1)'
        MAXOPT=9

      ELSE IF(ITYP5(nr,nx).EQ.2) THEN               !Time integration
        FORMAT='('' Specify equation ['//CHAR//']:'''//
     '    '/''   (1) Linear elasticity           '''//
     '    '/''   (2) Finite elasticity           '''//
     '    '/''   (3) Advection-diffusion         '''//
     '    '/''   (4) Wave equation               '''//
     '    '/''   (5) Navier-Stokes equations     '''//
     '    '/''   (6) Bio-heat transfer           '''//
     '    '/''  *(7) Maxwells equations          '''//
     '    '/''   (8) Huygens activation          '''//
     '    '/''   (9) Cellular based modelling    '''//
     '    '/''  (10) Oxygen transport            '''//
     '    '/''  (11) Pulmonary transport         '''//
     '    '/''  (12)                             '''//
     '    '/$,''   '',I2)'
        MAXOPT=12
        
      ELSE IF(ITYP5(nr,nx).EQ.3) THEN               !Modal analysis
        FORMAT='('' Specify equation ['//CHAR//']:'''//
     '    '/''   (1) Linear elasticity           '''//
     '    '/''   (2)                             '''//
     '    '/''   (3) Laplace equation            '''//
     '    '/''   (4) Helmholtz equation          '''//
     '    '/''   (5) div(k.grad(u))=f            '''//
     '    '/''   (6) Linear 2nd order elliptic   '''//
     '    '/''   (7) Biharmonic equation         '''//
     '    '/''   (8) Vocal tract equations       '''//
     '    '/''   (9)                             '''//
     '    '/$,''    '',I1)'
        MAXOPT=9

c cpb 1/5/95 Replacing Fourier analysis with Static analysis
c      ELSE IF(ITYP5(nr,nx).EQ.4) THEN               !Fourier Transform
c        FORMAT='('' Specify equation ['//CHAR//']:'''//
c     '    '/''   (1) Linear elasticity           '''//
c     '    '/''   (2) Finite elasticity           '''//
c     '    '/$,''    '',I1)'
c        MAXOPT=2

      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static analysis
        FORMAT='('' Specify equation ['//CHAR//']:'''//
     '    '/''   (1) Linear elasticity           '''//
     '    '/''   (2) Unused                      '''//
     '    '/''   (3) Laplace equation            '''//
     '    '/''   (4) Unused                      '''//
     '    '/''   (5) div(k.grad(u))=f            '''//
     '    '/''   (6) Unused                      '''//
     '    '/''   (7) Unused                      '''//
     '    '/''   (8) Unused                      '''//
     '    '/''   (9) Unused                      '''//
     '    '/$,''    '',I1)'
        MAXOPT=5

      ELSE IF(ITYP5(nr,nx).EQ.5) THEN        !Wavefront path analysis
        IDEFLT(1)=3
        IF(nr.GT.1) THEN
          IDEFLT(1)=ITYP2(nr-1,nx)
        ENDIF
        WRITE(CHAR,'(I1)') IDEFLT(1)
        FORMAT='('' Specify equation ['//CHAR//']:'''//
     '    '/''   (1)           '''//
     '    '/''   (2)           '''//
     '    '/''   (3) Elliptic 1st order eikonal  '''//
     '    '/''   (4)           '''//
     '    '/''   (5)           '''//
     '    '/''   (6)           '''//
     '    '/''   (7)           '''//
     '    '/''   (8) Huygen''''s activation        '''//
     '    '/''   (9)           '''//
     '    '/$,''    '',I1)'
        MAXOPT=9

      ELSE IF(ITYP5(nr,nx).EQ.6) THEN               !Buckling analysis
        FORMAT='('' Specify equation ['//CHAR//']:'''//
     '    '/''   (1) Linear elasticity           '''//
     '    '/''   (2) Finite elasticity           '''//
     '    '/$,''    '',I1)'
        MAXOPT=2
      ENDIF

      IF(IOTYPE.EQ.3) IDATA(1)=ITYP2(nr,nx)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,MAXOPT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) ITYP2(nr,nx)=IDATA(1)

      IF(ITYP2(nr,nx).EQ.1) THEN !linear elasticity
        CALL ASSERT(CALL_FIBR,
     '    ' Fibres must be defined for linear elasticity problems',
     '    ERROR,*9999)
        CALL ASSERT(CALL_ELFB,
     '    ' Element fibres are required for linear elasticity problems',
     '    ERROR,*9999)
      ELSEIF(ITYP2(nr,nx).EQ.2) THEN !finite elasticity
        CALL ASSERT(CALL_FIBR,
     '    ' Fibres must be defined for finite elasticity problems',
     '    ERROR,*9999)
        CALL ASSERT(CALL_ELFB,
     '    ' Element fibres are required for finite elasticity problems',
     '    ERROR,*9999)
      ENDIF

      NOQUES=0
      IF(ITYP2(nr,nx).EQ.1) THEN !linear elasticity
        NCT(nr,nx)=2 !variable and reaction

C 24/2/97 LC removed section : GMH 5/12/95 Change method of asking for element type

        !Initialise first, so we ensure that a type is given to each
        !element.
        IF(IOTYPE.NE.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NW(ne,1)=0 !invalid element type
          ENDDO
        ENDIF
        FORMAT='('' Enter element type:'''//
     '    '/''   (1) Truss or cable      (5) Membrane           '//
     '    ' ( 9) 3-dimensional'''//
     '    '/''   (2) Batten              (6) Plate (Kirchhoff)  '//
     '    ' (10) Tank Bottom '''//
     '    '/''   (3) Beam  (Kirchhoff)   (7) Shell              '//
     '    ' (11) Plane stress '''//
     '    '/''   (4) Link                (8) Shell/Fluid I-face '//
     '    ' (12) Plane strain '''//
     '    '/'' '')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C GMH 5/12/95 Copied from fe50.f IPINI5
        IF(NJT.EQ.3) THEN
          LAST_TYPE=9 !3D element
        ELSE
          LAST_TYPE=1 !plane stress
        ENDIF
 6720   FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
        IF(IOTYPE.EQ.3) THEN
          noelem=noelem+1
          IF(noelem.LE.NEELEM(0,nr)) THEN
            ne=NEELEM(noelem,nr)
            IDATA(1)=ne
            NELIST(0)=1
            NELIST(1)=ne
          ELSE
            IDATA(0)=0
          ENDIF
        ENDIF
 6760   CDATA(1)='ELEMENTS' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '    0,NET(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '    ERROR,*9999)
        IF(IDATA(1).NE.0) THEN !not default exit
          NELIST(0)=IDATA(0)
          DO n=1,IDATA(0)
            NELIST(n)=IDATA(n)
            ne=IDATA(n)
            IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '        NEELEM(0,nr),N1)) THEN
              WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '          //'in the current region'')') ne
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              GOTO 6760
            ENDIF
          ENDDO !n

C         Get type for first element in group
          ne=NELIST(1) !rest of group filled at end of loop
          IDEFLT(1)=LAST_TYPE
          WRITE(CHAR2,'(I3)') IDEFLT(1)
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          FORMAT='($,'' Element type is ['//CHAR2(IBEG2:IEND2)//
     '      ']: '',I3)'
          IF(IOTYPE.EQ.3) IDATA(1)=NW(ne,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,12,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) LAST_TYPE=IDATA(1)

C         Apply type to all of elements group
          DO n=1,NELIST(0)
            ne=NELIST(n)
            NW(ne,1)=LAST_TYPE
            IF(NW(ne,1).EQ.4) THEN
C ***       Link elements: length approx zero so set SE=1 for y-deriv
              nb=NBJ(2,ne)
              SE(2,nb,ne)=1.d0
              SE(4,nb,ne)=1.d0
            ELSE IF(NW(ne,1).EQ.6) THEN !Plate element s.b. bilin geom
              nb=NBJ(1,ne)
              CALL ASSERT(NNT(nb).EQ.4,
     '          '>>Plate must have bilinear geometry',ERROR,*9999)
              IF(NJT.EQ.3) THEN
                CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).GT.0,
     '            '>>3D Plate must have fibre field',ERROR,*9999)
                nc=1 !Temporary AJP 17-12-91
                CALL ASSERT(NHM*NSM.GE.12*NKT(0,NBH(3,nc,ne)),
     '            '>>NVM is too small',ERROR,*9999)
              ENDIF
            ELSE IF(NW(ne,1).EQ.9) THEN
              CALL ASSERT(NJT.EQ.3,'>>Need 3D coordinates',ERROR,*9999)
              nb=NBJ(1,ne)
              CALL ASSERT(NIT(nb).EQ.3,'>>Need 3D basis',ERROR,*9999)
            ENDIF
          ENDDO !n

          GO TO 6720 !for more elements
        ENDIF !idata(1).ne.0

C       Check that all elements have been defined
        DO ie=1,12
          ETYP(ie)=.FALSE.
        ENDDO
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          CALL ASSERT(NW(ne,1).NE.0,
     '      '>>Element type not defined for all elements',ERROR,*9999)
          ETYP(NW(ne,1))=.TRUE.
        ENDDO

        FORMAT=
     '   '(/'' Specify whether the thermal strain effects are [0]: '''//
     '    '/''   (0) Not included'''//
     '    '/''   (1) Included as fixed initial strain'''//
     '    '/''   (2) Included but constrained by displacement b.c.s'''//
     '    '/''   (3) Included with full coupling to mechanics'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=0
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP43
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP43=IDATA(1)

        FORMAT='(''  '')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

      ELSE IF(ITYP2(nr,nx).EQ.2) THEN !finite elasticity
        NCT(nr,nx)=2 !variable and reaction
        FORMAT='('' Specify problem type [1]: '''//
     '    '/''   (1) Plane stress      '''//
     '    '/''   (2) Plane strain      '''//
     '    '/''   (3) 3-dimensional     '''//
     '    '/''   (4) Membrane theory   '''//
     '    '/''   (5) String theory     '''//
     '    '/''  *(6) Shell theory      '''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP51(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,6,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          ITYP3(nr,nx)=IDATA(1)
          KTYP51(nr)=ITYP3(nr,nx)
        ENDIF

        FORMAT=
     '    '('' Specify whether the dependent variables are [1]: '''//
     '    '/''   (1) Geometric coordinates'''//
     '    '/''   (2) Displacements'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP58(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP58(nr)=IDATA(1)

      ENDIF

C MLB 15-Oct-1998 use THETA in IPSOLV file
C      IF(ITYP5(nr,nx).EQ.2) THEN !time integration
C        IF(ITYP2(nr,nx).EQ.9) THEN !Ionic current activation
C          FORMAT='('' Specify equation type [1]:'''//
C     '      '/''   (1) Explicit   '''//
C     '      '/''   (2) Implicit   '''//
C     '      '/''   (3) Crank Nicholson   '''//
C     '      '/$,''    '',I1)'
C          IDEFLT(1)=1
C          IF(IOTYPE.EQ.3) IDATA(1)=ITYP16(nr,nx)
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) ITYP16(nr,nx)=IDATA(1)
C        ENDIF
C      ENDIF

      MAXOPT=0
      NOQUES=0
      IF(ITYP5(nr,nx).EQ.1) THEN !static analysis
        IF(ITYP2(nr,nx).EQ.3) THEN !Laplace equation
          IF(nr.GT.1) THEN
            IDEFLT(1)=ITYP3(nr-1,nx)
          ENDIF
          IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
          WRITE(CHAR,'(I1)') IDEFLT(1)
          NCT(nr,nx)=2 !variable and flux
! Added AJP 24-2-94 to include conductivity in Laplace equation
          FORMAT='('' Specify equation ['//CHAR//']:'''//
     '      '/''   (1) Standard Laplace               '''//
     '      '/''   (2) Generalised Laplace            '''//
     '      '/''   (3) Potential flow around aerofoil '''//
     '      '/''   (4) Unused                         '''//
     '      '/''   (5) Unused                         '''//
     '      '/$,''    '',I2)'
          MAXOPT=3
        ELSE IF(ITYP2(nr,nx).EQ.4) THEN !Helmholtz equation
          NCT(nr,nx)=2 !variable and flux
          IDEFLT(1)=1
          FORMAT='('' Specify equation [1]:'''//
     '      '/''   (1) Standard Helmholtz          '''//
     '      '/''   (2) Modified Helmholtz (Yukawa) '''//
     '      '/$,''    '',I2)'
          MAXOPT=2
! news AJP 3/9/97
        ELSE IF(ITYP2(nr,nx).EQ.5) THEN !Poisson equation
          IF(nr.GT.1) THEN
            IDEFLT(1)=ITYP3(nr-1,nx)
          ENDIF
          IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
          WRITE(CHAR,'(I1)') IDEFLT(1)
          NCT(nr,nx)=2 !variable and flux
          FORMAT='('' Specify form of f ['//CHAR//']:'''//
     '      '/''   (1) f=General                      '''//
     '      '/''   (2) f=-div(s.grad(g))              '''//
     '      '/''   (3) Unused                         '''//
     '      '/$,''    '',I2)'
          MAXOPT=3
! newe AJP
        ELSE IF(ITYP2(nr,nx).EQ.7) THEN !Maxwells
          IDEFLT(1)=1
          WRITE(CHAR,'(I1)') IDEFLT(1)
          NCT(nr,nx)=2 !variable and flux
          FORMAT='('' Specify equation ['//CHAR//']:'''//
     '      '/''  *(1) Electrostatic                  '''//
     '      '/''   (2) Magnetostatic - no gauge       '''//
     '      '/''   (3) Magnetostatic - Coulomb gauge  '''//
     '      '/''   (4) Unused                         '''//
     '      '/''   (5) Unused                         '''//
     '      '/$,''    '',I2)'
          MAXOPT=3

        ELSE IF(ITYP2(nr,nx).EQ.8) THEN !Fluid mechanics
          NCT(nr,nx)=2 !variable and reaction
          KTYP51(nr)=3 !3D analysis
          IDEFLT(1)=1
          FORMAT='('' Specify equation [1]:'''//
     '      '/''   (1) Prandtl boundary layer eqtns'''//
     '      '/''   (2) Fluid interface stability   '''//
     '      '/''   (3) Constant volume constraint  '''//
     '      '/''   (4)                             '''//
     '      '/''   (5)                             '''//
     '      '/$,''    '',I1)'
          MAXOPT=5
        ELSE IF(ITYP2(nr,nx).EQ.9) THEN !Oxygen transport
          NCT(nr,nx)=2 !variable and ???
          IDEFLT(1)=1
          FORMAT='('' Specify equation [1]:'''//
     '      '/''   (1) Multi-field oxygen transport'''//
     '      '/''   (2) Glucose-Oxygen transport    '''//
     '      '/$,''    '',I1)'
          MAXOPT=2
        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.2) THEN !time integration
        IF(ITYP2(nr,nx).EQ.5) THEN !Navier-Stokes equations
          IF(nr.GT.1) THEN
            IDEFLT(1)=ITYP3(nr-1,nx)
          ENDIF
          IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
          WRITE(CHAR,'(I1)') IDEFLT(1)
          NCT(nr,nx)=2 !variable and ???
          FORMAT='('' Specify equation ['//CHAR//']:'''//
     '      '/''   (1) Fluid in elastic tube       '''//
     '      '/''   (2) Lung gas transport          '''//
     '      '/''   (3) General Navier-Stokes eqns  '''//
     '      '/''   (4) Stokes flow (no advection)  '''//
     '      '/''   (5) Bidomain Stokes flow        '''//
     '      '/$,''    '',I1)'
          MAXOPT=5
        ELSE IF(ITYP2(nr,nx).EQ.9) THEN !cellular based modelling
          CALL ASSERT(NQT.GT.0,'Need to define grid first',ERROR,*9999)
          NCT(nr,nx)=1 !variable
          IDEFLT(1)=1
          FORMAT='('' Specify cellular model type [1]:'''//
     '      '/''   (1) Electrical                    '''//
     '      '/''   (2) Mechanical                    '''//
     '      '/''   (3) Metabolism                    '''//
     '      '/''   (4) Signalling Pathways           '''//
     '      '/''   (5) Drug Interaction              '''//
     '      '/''   (6)                               '''//
     '      '/''   (7) Coupled                       '''//
     '      '/$,''    '',I1)'
          MAXOPT=7
          IF(IOTYPE.EQ.3) IDATA(1)=ITYP19(nr,nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,MAXOPT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ITYP19(nr,nx)=IDATA(1)

          IF(ITYP19(nr,nx).EQ.1) THEN !Electrical cellular modelling
            IDEFLT(1)=1
            FORMAT='('' Specify electrical model [1]:'''//
     '        '/''   (1) Cubic - no recovery           '''//
     '        '/''   (2) FitzHugh-Nagumo               '''//
     '        '/''   (3) van Capelle-Durrer            '''//
     '        '/''   (4) Beeler-Reuter                 '''//
     '        '/''   (5) Jafri-Rice-Winslow            '''//
     '        '/''   (6) Luo-Rudy                      '''//
     '        '/''   (7) diFrancesco-Noble             '''//
     '        '/''   (8) Noble-98                      '''//
     '        '/''   (9) Hodgkin-Huxley                '''//
     '        '/''  (10) User defined                  '''//
     '        '/$,''    '',I1)'
            MAXOPT=10
          ELSEIF(ITYP19(nr,nx).EQ.2) THEN !Mechanical cellular modelling
            IDEFLT(1)=1
            FORMAT='('' Specify mechanical model [1]:'''//
     '        '/''   (1) Hunter-McCulloch-ter Keurs    '''//
     '        '/''   (2) Fading Memory                 '''//
     '        '/''   (3) Distribution Moment           '''//
     '        '/''   (4) Infarct                       '''//
     '        '/''   (5) User defined                  '''//
     '        '/$,''    '',I1)'
            MAXOPT=5
          ELSEIF(ITYP19(nr,nx).EQ.3) THEN !Metabolic cellular modelling
            IDEFLT(1)=1
            FORMAT='('' Specify mechanical model [1]:'''//
     '        '/''   (1) myocite    '''//
     '        '/$,''    '',I1)'
            MAXOPT=1
          ELSEIF(ITYP19(nr,nx).EQ.4) THEN !Signalling pathway modelling
            CALL ASSERT(.FALSE.,'Cell signalling not yet implemented',
     '        ERROR,*9999)
          ELSEIF(ITYP19(nr,nx).EQ.5) THEN !Drug interaction modelling
            CALL ASSERT(.FALSE.,'Drug interaction not yet implemented',
     '        ERROR,*9999)
          ELSEIF(ITYP19(nr,nx).EQ.6) THEN !??
            CALL ASSERT(.FALSE.,'Not implemented',
     '        ERROR,*9999)
          ELSEIF(ITYP19(nr,nx).EQ.7) THEN !Coupled cellular modelling
            IDEFLT(1)=1
            FORMAT='('' Specify coupled model [1]:'''//
     '        '/''   (1) LR - HMT                      '''//
     '        '/''   (2) Noble 98 - HMT                '''//
     '        '/''   (3) HMT - SPM                     '''//
     '        '/''   (4) User Defined                  '''//
     '        '/$,''    '',I1)'
            MAXOPT=4
          ENDIF!ITYP19(nr,nx)
        ELSE IF(ITYP2(nr,nx).EQ.10) THEN !Oxygen transport
          NCT(nr,nx)=2 !variable and ???
          IDEFLT(1)=1
          FORMAT='('' Specify equation [1]:'''//
     '      '/''   (1) Multi-field oxygen transport'''//
     '      '/''   (2) Glucose-Oxygen transport    '''//
     '      '/$,''    '',I1)'
          MAXOPT=2
        ELSE IF(ITYP2(nr,nx).EQ.11) THEN !pulmonary transport
          IF(ITYP3(nr,nx).LE.2) THEN
            NCT(nr,nx)=2
          ELSE IF(ITYP3(nr,nx).GE.3) THEN !capillary, poiseuelle, MCT
            NCT(nr,nx)=1
          ENDIF
          IDEFLT(1)=1
          FORMAT='('' Specific equation [1]:'''//
     &      '/''   (1) Airway advection-diffusion           '''//
     &      '/''   (2) Coupled airway water/heat transfer   '''//
     &      '/''   (3) Capillary blood flow                 '''//
     &      '/''   (4) Simple pressure-resistance-flow      '''//
     &      '/''   (5) Mucociliary transport                '''//
     &      '/''   (6) Arteriole-cap-venule blood flow '''//
     &      '/''   (7) Gas exchange with pre-defined trasport'''//
     &      '/$,''    '',I1)'
          MAXOPT=7
        ENDIF
      ELSE IF(ITYP5(nr,nx).EQ.3) THEN !Modal analysis
        IF(ITYP2(nr,nx).EQ.4) THEN !Helmholtz equation
          NCT(nr,nx)=2 !variable and ???
          IDEFLT(1)=1
          FORMAT='('' Specify equation [1]:'''//
     '      '/''   (1) Standard Helmholtz          '''//
     '      '/''   (2) Modified Helmholtz (Yukawa) '''//
     '      '/$,''    '',I2)'
          MAXOPT=2
        ENDIF

c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier Transform
c        NCT(nr,nx)=2 !variable and ???

      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static analysis
        IF(ITYP2(nr,nx).EQ.3) THEN !Laplace equation
          IF(nr.GT.1) THEN
            IDEFLT(1)=ITYP3(nr-1,nx)
          ENDIF
          IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
          WRITE(CHAR,'(I1)') IDEFLT(1)
          NCT(nr,nx)=2 !variable and flux
          FORMAT='('' Specify equation ['//CHAR//']:'''//
     '      '/''   (1) Standard Laplace          '''//
     '      '/''   (2) Generalised Laplace       '''//
     '      '/''   (3) Unused                    '''//
     '      '/''   (4) Unused                    '''//
     '      '/''   (5) Unused                    '''//
     '      '/$,''    '',I2)'
          MAXOPT=2
! news AJP 3/9/97
        ELSE IF(ITYP2(nr,nx).EQ.5) THEN !Poisson equation
          IF(nr.GT.1) THEN
            IDEFLT(1)=ITYP3(nr-1,nx)
          ENDIF
          IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
          WRITE(CHAR,'(I1)') IDEFLT(1)
          NCT(nr,nx)=2 !variable and flux
          FORMAT='('' Specify form of f ['//CHAR//']:'''//
     '      '/''   (1) f=General                      '''//
     '      '/''   (2) f=-div(s.grad(g))              '''//
     '      '/''   (3) Unused                         '''//
     '      '/$,''    '',I2)'
          MAXOPT=3
! newe AJP
        ENDIF
      ELSE IF(ITYP5(nr,nx).EQ.5) THEN !Wavefront path analysis
        IF(ITYP2(nr,nx).EQ.3) THEN !Elliptic eikonal
          IDEFLT(1)=1
          IF(nr.GT.1) THEN
            IDEFLT(1)=ITYP3(nr-1,nx)
          ENDIF
          WRITE(CHAR,'(I1)') IDEFLT(1)
          FORMAT='('' Specify nature of material ['//CHAR//']:'''//
     '      '/''   (1) Isotropic                 '''//
     '      '/''   (2) Orthotropic monodomain    '''//
     '      '/''  *(3) Orthotropic bidomain      '''//
     '      '/''   (4) Unused                    '''//
     '      '/''   (5) Unused                    '''//
     '      '/$,''    '',I1)'
          MAXOPT=2
C         NCT is set later
        ENDIF
      ENDIF

      IF(MAXOPT.EQ.0) THEN
        ITYP3(nr,nx)=1
      ELSE IF(MAXOPT.GT.0) THEN
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP3(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,MAXOPT,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP3(nr,nx)=IDATA(1)
      ENDIF

C    AJS Sep 2007 Pressure-flow-resistance equations for blood
      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.11.
     &  AND.ITYP3(nr,nx).EQ.4)THEN
        COUPLE_VIA_LPM='N' !Initialise 
        IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6)THEN !vessels only
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=NJ_TYPE(nj_radius_R0,2) !nj_radius_R0 defined in lung00.cmn
          ENDIF
          IDEFLT(1)=1
          FORMAT='(/$,'' The field number for the unstressed vessel '
     &      //'radius (R0) is [2]: '',I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &      0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &      *9999)
          CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &    '>>Define field first',ERROR,*9999)
          IF(IOTYPE.NE.3)THEN
            njj=IDATA(1)
            nj_radius_R0=NJ_LOC(NJL_FIEL,njj,nr)
          ENDIF
C.... KSB Option to allow coupling arterial-capillary-venous flow
C.... circuit
          FORMAT='($,'' Is circulatory network coupled via micro-'//
     &      'circulation LPM [Y]? '',A)'
          ADEFLT(1)='Y'
          IF(IOTYPE.EQ.3) ADATA(1)=COUPLE_VIA_LPM
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &       IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     &      ERROR,*9999)
          IF(IOTYPE.NE.3) COUPLE_VIA_LPM=ADATA(1)        
          FORMAT='($,'' Do you want to include sheet model [Y]? '',A)'
          ADEFLT(1)='Y'
          IF(IOTYPE.EQ.3) ADATA(1)=SHEET_FLOW
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &      IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     &      ERROR,*9999)
          IF(IOTYPE.NE.3) SHEET_FLOW=ADATA(1)       
          IF(SHEET_FLOW.EQ.'Y') THEN
             FORMAT='('' Ladder model? [1]: '''//
     &            '/''   (1) No ladder, sheets at terminal     '''//
     &            '/''   (2) Ladder model included             '''//
     &            '/''   (3) Multibranching acinar ladder      '''//
     &            '/''   (4) Include supernumerary vessels     '''//
     &            '/$,''    '',I1)'
             IDEFLT(1)=1
             IF(IOTYPE.EQ.3) IDATA(1)=LADDER
             CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &            1,4,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &            *9999)
             IF(IOTYPE.NE.3) LADDER=IDATA(1)             
          ENDIF
       ENDIF
      ENDIF

C    AJS Sep 2007 Pulmonary gas exchange with pre-defined transport
      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.11.
     &  AND.ITYP3(nr,nx).EQ.7)THEN
        ADEFLT(1)='N'
        FORMAT='(/$,'' Hold alveolar air partial pressures constant '
     &    //'[N]? '',A)'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &    0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &    *9999)
        IF(IOTYPE.EQ.3) ADATA(1)='N'
        IF(ADATA(1).EQ.'Y') CONST_ALV_PRESSURES=.TRUE.
        IF(ADATA(1).EQ.'N') CONST_ALV_PRESSURES=.FALSE.
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NJ_TYPE(nj_Vdot,2)
        ENDIF
        IDEFLT(1)=1
        FORMAT='(/$,'' The field number for the volume change '
     &    //'(Vdot) is [1]: '',I2)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &    0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &    *9999)
        CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &  '>>Define field first',ERROR,*9999)
        IF(IOTYPE.NE.3)THEN
          njj=IDATA(1)
          nj_Vdot=NJ_LOC(NJL_FIEL,njj,nr)
        ENDIF
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NJ_TYPE(nj_PA,2)
        ENDIF
        IDEFLT(1)=1
        FORMAT='(/$,'' The field number for the alveolar pressure '
     &    //'(PA) is [1]: '',I2)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &    0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &    *9999)
        CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &  '>>Define field first',ERROR,*9999)
        IF(IOTYPE.NE.3)THEN
          njj=IDATA(1)
          nj_PA=NJ_LOC(NJL_FIEL,njj,nr)
        ENDIF
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NJ_TYPE(nj_Qdot,2)
        ENDIF
        IDEFLT(1)=1
        FORMAT='(/$,'' The field number for the blood flow '
     &    //'(Qdot) is [1]: '',I2)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &    0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &    *9999)
        CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &  '>>Define field first',ERROR,*9999)
        IF(IOTYPE.NE.3)THEN
          njj=IDATA(1)
          nj_Qdot=NJ_LOC(NJL_FIEL,njj,nr)
        ENDIF
      ENDIF

      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.1)
     &  THEN !Navier-stokes; elastic tubes:
        !pulmonary or coronary (or other) flow problems
        FORMAT=
     '    '('' Specify whether flow solution setting [1]: '''//
     '    '/''   (1) Coronary'''//
     '    '/''   (2) Pulmonary'''//
     '    '/''   (3) Pulmonary + external pressure field'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP12(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP12(nr,nx)=IDATA(1)

      ELSE IF(ITYP2(nr,nx).EQ.10.AND.ITYP3(nr,nx).EQ.1) THEN !Multi-field
        NCT(nr,nx)=2 !variable and ???                    !oxygen
        FORMAT='($,'' Specify # of fields [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP15
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP15=IDATA(1)

      ELSE IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.8.
     '  AND.ITYP3(nr,nx).EQ.2) THEN             !Fluid interface
        NCT(nr,nx)=2 !variable and ???            !stability
        RDEFLT(1)=1.0d0
        FORMAT='($,'' Specify the time increment [1.0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=TINCR
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TINCR=RDATA(1)
        RDEFLT(1)=9.81d0
        FORMAT='($,'' Specify the gravitational acceleration '
     '    //'[9.81]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=G_ACCEL
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) G_ACCEL=RDATA(1)

c AJS 11/2010 Add gas exchange option for advection-diffusion problem
      ELSE IF(ITYP2(nr,nx).EQ.11.AND.(ITYP3(nr,nx).EQ.7.OR.
     &  ITYP3(nr,nx).EQ.1))THEN !lung gas exchange 
        IDEFLT(1)=1
        FORMAT='('' Specify the gas exchange model [1]:'''//
     &    '/''   (1) No gas exchange                      '''//
     &    '/''   (2) Ben-Tal 2006 model        '''//
     &    '/$,''    '',I1)'
          MAXOPT=2
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP7(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,MAXOPT,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP7(nr,nx)=IDATA(1)
!         CALL ASSERT((ITYP3(nr,nx).EQ.7.AND.ITYP7(nr,nx).GT.1),
!      '    'Cannot choose option 1 for gas exchange with '//
!      '    'pre-defined transport',ERROR,*9999) !error not relevant anymore
      ENDIF

      IF((ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.5).
     '  AND.ITYP2(nr,nx).EQ.8) !Huygen's activation
     '  THEN
        ITYP5(nr,nx)=2 !time integration
        NCT(nr,nx)=1 !variable
        ITYP4(nr,nx)=4 !collocation
        IF(ITYP3(nr,nx).EQ.1) THEN
          FORMAT='($,'' Specify whether model is (1)forwards,'
     '      //' or (2)backwards in time [1]: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP31
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP31=IDATA(1)
        ELSE
          KTYP31=1
        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN !cellular based modelling
        NCT(nr,nx)=1 !variable
C SGM 10-10-00  ITYP4 now user defined (below)
C       ITYP4(nr,nx)=4 !collocation
        ITYP1(nr,nx)=3 !fe30 problem
        ITYP6(nr,nx)=1 !linear

        IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7) THEN
!electrical or coupled modelling
          FORMAT=
     '      '($,'' Enter (1) monodomain, or (2) bidomain [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP32
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP32=IDATA(1)
C MLB why is this necessary?
C        IF(KTYP32.EQ.2) THEN  !Bidomain model
C          WRITE(OP_STRING,'('' All dimensions must be in mm. '')')
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        ENDIF

          IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.1) THEN !Cubic
            FORMAT='($,'' Enter (1) cubic, (2) quintic or '
     '        //'(3) heptic <order 7> [1]:'',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP33
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP33=IDATA(1)

          ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.2) THEN !FHN
            FORMAT='($,'' Enter (1) standard, (2) Rogers or '
     '        //'(3) Panfilov [1]:'',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP33
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP33=IDATA(1)

          ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.3) THEN !van Capelle-Durer
            FORMAT=
     '        '($,'' Enter (1) standard, or (2) Calif.mod [1]:'',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP33
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP33=IDATA(1)

            IF(KTYP33.EQ.2) THEN !VCD-California
              FORMAT=
     '          '($,'' Enter (1) normal, or (2) ischemic [1]:'',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=KTYP34
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) KTYP34=IDATA(1)
            ENDIF

          ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.4) THEN !Beeler-Reuter
            IDEFLT(1)=1
            WRITE(CHAR,'(I1)') IDEFLT(1)
            FORMAT='('' Specify sodium kinetics ['//CHAR//']:'''//
     '        '/''   (1) Standard              '''//
     '        '/''   (2) Ebihara-Johnson       '''//
     '        '/''   (3) Drouhard-Roberge      '''//
     '        '/''   (4) Defibrillation DR     '''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP33
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP33=IDATA(1)

            IF(KTYP33.EQ.4) THEN !defib
              FORMAT='($,'' Include electroporation (Y/N) [N]? '',A)'
              IF(IOTYPE.EQ.3) THEN
                IF(KTYP39.EQ.1) THEN
                  ADATA(1)='Y'
                ELSE
                  ADATA(1)='N'
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     '          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                IF(ADATA(1).EQ.'Y') THEN
                  KTYP39=1
                ELSE
                  KTYP39=0
                ENDIF
              ENDIF

              RDEFLT(1)=8.0d0
              WRITE(CHAR2,'(F8.2)') RDEFLT(1)
              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
              FORMAT='($,'' Enter the scale factor for reducing the '
     '          //'APD ['//CHAR2(IBEG2:IEND2)//']: '',F8.2)'
              IF(IOTYPE.EQ.3) RDATA(1)=KTYP39R
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,99.0d0,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) KTYP39R=RDATA(1)
            ELSE
              !value used for all non-defib BR models
              KTYP39R=1.0d0
            ENDIF

          ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.5) THEN !Jaffri-Rice-Winslow
            IDEFLT(1)=1
            WRITE(CHAR,'(I1)') IDEFLT(1)
            FORMAT='('' Specify JRW model type ['//CHAR//']:'''//
     '        '/''   (1) Standard              '''//
     '        '/''   (2) Princeton             '''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP33
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP33=IDATA(1)
          ELSE IF((ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.6).OR.
     '       (ITYP19(nr,nx).EQ.7.AND.ITYP3(nr,nx).EQ.1)) THEN
!Luo-Rudy or LR-HMT
          ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.7) THEN !diFrancesco-Noble
          ELSE IF((ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.8).OR.
     '        (ITYP19(nr,nx).EQ.7.AND.ITYP3(nr,nx).EQ.2)) THEN !Noble98 or Noble98 - HMT
          ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.9)  THEN !Hodgkin-Huxley
          ELSE IF((ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.10).OR.
     '        (ITYP19(nr,nx).EQ.7.AND.ITYP3(nr,nx).EQ.4)) THEN
!User defined
            IDEFLT(1)=1
            WRITE(CHAR,'(I1)') IDEFLT(1)
            FORMAT=
     '        '('' Specify user defined model type ['//CHAR//']:'''//
     '        '/''   (1) Model 1    '''//
     '        '/''   (2) Model 2    '''//
     '        '/''   (3) Model 3    '''//
     '        '/''   (4) Model 4    '''//
     '        '/''   (5) Model 5    '''//
     '        '/''   (6) cellML     '''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP33
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,6,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP33=IDATA(1)
          ENDIF !ityp3
        ELSEIF(ITYP19(nr,nx).EQ.2) THEN !mechanical cellular modelling
          IF(ITYP3(nr,nx).EQ.1) THEN !HMT mechanics model
            !temp hack for now!!
            CALL_MATE=.TRUE.
          ELSEIF(ITYP3(nr,nx).EQ.2) THEN !Fading Memory model
            CALL ASSERT(.FALSE.,
     '        'Fading memory model not yet implemented',ERROR,*9999)
          ELSEIF(ITYP3(nr,nx).EQ.3) THEN !DM model
            CALL ASSERT(.FALSE.,
     '        'Distribution moment model not yet implemented',
     '        ERROR,*9999)
          ELSEIF(ITYP3(nr,nx).EQ.4) THEN !Infarct
            CALL_MATE=.TRUE.
          ELSEIF(ITYP3(nr,nx).EQ.5) THEN !User defined
            CALL_MATE=.TRUE.
            IDEFLT(1)=1
            WRITE(CHAR,'(I1)') IDEFLT(1)
            FORMAT=
     '        '('' Specify user defined model type ['//CHAR//']:'''//
     '        '/''   (1) Model 1    '''//
     '        '/''   (2) Model 2    '''//
     '        '/''   (3) Model 3    '''//
     '        '/''   (4) Model 4    '''//
     '        '/''   (5) Model 5    '''//
     '        '/''   (6) cellML     '''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP33
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,6,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP33=IDATA(1)
          ENDIF !ityp3
        ENDIF !ityp19

C SGM 10-10-00 Added menu below to define ITYP4
C MLT 29Nov02 Added Finite volume technique to menu
        IDEFLT(1)=4
        WRITE(CHAR,'(I1)') IDEFLT(1)
        FORMAT='('' Specify whether solution is by ['//CHAR//']:'''//
     '    '/''   (1)  '''//
     '    '/''   (2)  '''//
     '    '/''   (3)  '''//
     '    '/''   (4) Collocation               '''//
     '    '/''   (5)  '''//
     '    '/''   (6) Grid-based Finite Element '''//
     '    '/''   (7) Grid-based Finite Volume   '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP4(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,4,7,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP4(nr,nx)=IDATA(1)

      ELSE IF(ITYP5(nr,nx).EQ.5) THEN !non-Huygen wavefront path analysis
        NCT(nr,nx)=2 !variable and fluxy thing

        IDEFLT(1)=1
        IF(nr.GT.1) THEN
          IDEFLT(1)=ITYP4(nr-1,nx)
        ENDIF
        WRITE(CHAR,'(I1)') IDEFLT(1)
        FORMAT='('' Specify whether solution is by ['//CHAR//']:'''//
     '    '/''   (1) Finite elements  '''//
     '    '/''   (2)  '''//
     '    '/''   (3)  '''//
     '    '/''   (4)  '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP4(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP4(nr,nx)=IDATA(1)

        IDEFLT(1)=3
        IF(nr.GT.1) THEN
          IDEFLT(1)=ITYP15(nr-1,nx)
        ENDIF
        WRITE(CHAR,'(I1)') IDEFLT(1)
        FORMAT='('' Specify finite element method ['//CHAR//']:'''//
     '    '/''   (1) Galerkin                                 '''//
     '    '/''   (2) Petrov-Galerkin with material derivatives'''//
     '    '/''   (3) Petrov-Galerkin with xi derivatives      '''//
     '    '/''   (4) Stabilized Galerkin                      '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP15(nr,nx)+1
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP15(nr,nx)=IDATA(1)-1

        ITYP1(nr,nx)=3 !FE30
        ITYP6(nr,nx)=2 !nonlinear
        ITYP7(nr,nx)=1 !nonlinear equation

      ELSE !everything bar activation problems
        NCT(nr,nx)=2 !variable and flux
        IF(nr.GT.1) THEN
          IDEFLT(1)=ITYP4(nr-1,nx)
        ENDIF
        IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
        WRITE(CHAR,'(I1)') IDEFLT(1)
C MLT 29Nov02 Added Finite volume technique to menu
        FORMAT='('' Specify whether solution is by ['//CHAR//']:'''//
     '    '/''   (1) Galerkin finite elements  '''//
     '    '/''   (2) Direct boundary elements  '''//
     '    '/''   (3) Finite Difference         '''//
     '    '/''   (4) Collocation               '''//
     '    '/''   (5) Finite volume technique   '''//
     '    '/''   (6) Grid-based Finite Element '''//
     '    '/''   (7) Grid-based Finite Volume   '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP4(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP4(nr,nx)=IDATA(1)

        IF (ITYP4(nr,nx).EQ.3) THEN
          IF(nr.GT.1) THEN
            IDEFLT(1)=ITYP16(nr-1,nx)
          ENDIF
          IF(IDEFLT(1).EQ.0) IDEFLT(1)=4
C the only option implimented 31/10/96 NPS
          WRITE(CHAR,'(I1)') IDEFLT(1)
          FORMAT='('' Specify whether solution is by ['//CHAR//']:'''//
     '      '/''   (1) Fully explicit*            '''//
     '      '/''   (2) Fully implicit*           '''//
     '      '/''   (3) Crank Nicholson*           '''//
     '      '/''   (4) Lax Wendroff              '''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=ITYP16(nr,nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ITYP16(nr,nx)=IDATA(1)
        ENDIF

        IF(ITYP2(nr,nx).EQ.1) THEN !Linear elasticity
          ITYP1(nr,nx)=4 !FE40
        ELSE IF(ITYP2(nr,nx).EQ.2) THEN !Finite elasticity
          ITYP1(nr,nx)=5 !FE50
        ELSE IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
          ITYP1(nr,nx)=5 !Constant volume constraint
        ELSE IF(ITYP4(nr,nx).EQ.2) THEN !Direct BEM
          ITYP1(nr,nx)=9 !FE90
        ELSE IF(ITYP4(nr,nx).EQ.5) THEN !Finite volume
          ITYP1(nr,nx)=6 !FE60
        ELSE
          ITYP1(nr,nx)=3 !FE30
        ENDIF

C PM 26-JUL-01
C       IF(ITYP4(nr,nx).NE.2.AND.ITYP4(nr,nx).NE.3) THEN !not BEM
        IF(ITYP4(nr,nx).NE.2) THEN !not BEM
          IF(ITYP5(nr,nx).LE.4) THEN
C           static, time, quasistatic or modal analysis
            IF(ITYP5(nr,nx).EQ.1) THEN !static analysis
              IF(ITYP2(nr,nx).EQ.1.OR.ITYP2(nr,nx).EQ.3.OR
     '          .ITYP2(nr,nx).EQ.4.OR.ITYP2(nr,nx).EQ.5.OR
     '          .ITYP2(nr,nx).EQ.6.OR.ITYP2(nr,nx).EQ.7) THEN
                IDEFLT(1)=1
              ELSE IF(ITYP2(nr,nx).EQ.2) THEN !Finite elasticity
                IDEFLT(1)=2
              ELSE IF(ITYP2(nr,nx).EQ.8) THEN !Fluid mechanics
                IF(ITYP3(nr,nx).EQ.1) THEN      !Boundary layer eqns
                  IDEFLT(1)=2
                ELSE IF(ITYP3(nr,nx).EQ.2) THEN !Fluid interface
                  IDEFLT(1)=1                   !stability
                ELSE IF(ITYP3(nr,nx).EQ.3) THEN !Constant volume
                  IDEFLT(1)=2                   !constraint
                ENDIF
              ENDIF
            ELSE IF(ITYP5(nr,nx).EQ.2) THEN !time integration
              IF(ITYP2(nr,nx).EQ.1.OR.ITYP2(nr,nx).EQ.3.OR
     '          .ITYP2(nr,nx).EQ.4.OR.ITYP2(nr,nx).EQ.6.OR
     '          .ITYP2(nr,nx).EQ.7.OR
     '          .ITYP2(nr,nx).EQ.8.OR.ITYP2(nr,nx).EQ.9) THEN
                IDEFLT(1)=1
              ELSE IF(ITYP2(nr,nx).EQ.2.OR.ITYP2(nr,nx).EQ.10) THEN
                IDEFLT(1)=2
              ELSE IF(ITYP2(nr,nx).EQ.5) THEN
                IF(ITYP3(nr,nx).EQ.1) THEN      !Fluid in elastic tube
                  IDEFLT(1)=1
                ELSE IF(ITYP3(nr,nx).EQ.2) THEN !Lung gas transport
                  IDEFLT(1)=1
                ENDIF
              ELSE IF(ITYP2(nr,nx).EQ.11) THEN !pulmonary transport
                IDEFLT(1)=1
              ENDIF
            ELSE IF(ITYP5(nr,nx).EQ.3) THEN !modal analysis
              IDEFLT(1)=1
            ELSE IF(ITYP5(nr,nx).EQ.4) THEN !quasistatic analysis
              IDEFLT(1)=1
            ENDIF
            WRITE(CHAR,'(I1)') IDEFLT(1)
            FORMAT='('' Specify whether ['//CHAR(1:1)//']:'''//
     '        '/''   (1) Linear analysis    '''//
     '        '/''   (2) Nonlinear analysis '''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=ITYP6(nr,nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ITYP6(nr,nx)=IDATA(1)
          ELSE
            ITYP6(nr,nx)=1
          ENDIF
          IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
            FORMAT='('' Specify whether [1]:'''//
     '        '/''   (1) Nonlinear equation '''//
     '        '/''   (2) Nonlinear bdry conditions'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=ITYP7(nr,nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ITYP7(nr,nx)=IDATA(1)
            IF(ITYP7(nr,nx).EQ.2) THEN !Nonlinear bdry conditions
              FORMAT='('' Specify whether [1]:'''//
     '          '/''   (1) Radiation bdry conditions'''//
     '          '/''   (2) Unused'''//
     '          '/$,''    '',I1)'
              IDEFLT(1)=1
              IF(IOTYPE.EQ.3) IDATA(1)=ITYP8(nr,nx)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ITYP8(nr,nx)=IDATA(1)
            ENDIF
          ENDIF !ityp6

        ELSE !BEM
          ITYP6(nr,nx)=1
          !GMH 5/4/95.  If it is linear elasticity, then make sure that
          !all elements are the same, and either plane stress or strain
          IF(ITYP2(nr,nx).EQ.1) THEN !linear elasticity
            ne = NEELEM(1,nr)
            elemtype = NW(ne,1)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CALL ASSERT((NW(ne,1).EQ.11).OR.(NW(ne,1).EQ.12),'>>All'
     '          //' elements must be either plane stress or strain',
     '          ERROR,*9999)
              CALL ASSERT(NW(ne,1).EQ.elemtype,
     '          '>>All elements must be of the same type',ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF

        IF(ITYP5(nr,nx).EQ.2) THEN !time integration
          FORMAT='($,'' Are any equation coeffs time-varying [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
            FORMAT='(/'' Specify whether these coeffs are [1]:'''//
     '        '/''   (1) Constant with respect to time'''//
     '        '/''   (2) Defined in subroutine USER_IPEQUA'''//
     '        '/''   (3) Read from File.IPEQUA_time at each step'''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP3_equa(nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP3_equa(nx)=IDATA(1)
          ELSE
            KTYP3_equa(nx)=1
          ENDIF
          IF(ITYP2(nr,nx).EQ.3) THEN !advection-diffusion
            FORMAT='($,'' Specify the number of dependent variables '//
     '             '[1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP3A(nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP3A(nx)=IDATA(1)
          ELSE
            KTYP3A(nx)=1
          ENDIF

        ELSE IF(ITYP5(nr,nx).NE.2) THEN !not time integration
          KTYP3_equa(nx)=1
          KTYP3A(nx)=1
        ENDIF !ityp5

      ENDIF

      IF(NCT(nr,nx).GT.NCM) THEN
        IEND=0
        CALL APPENDC(IEND,'NCM needs to be at least ',ERROR)
        CALL APPENDI(IEND,NCT(nr,nx),ERROR)
        GOTO 9999
      ENDIF ! NCT > NCM

C new CS 15 Feb 2000 infarct modelling
      IF(.NOT.(ITYP19(nr,nx).EQ.2.AND.ITYP3(nr,nx).EQ.4)) THEN !Infarct
C old MPN 13-Aug-95: Should not be needed anymore
C      KTYP1 =ITYP2(nr,nx) !temporary only
C       KTYP11=ITYP3(nr,nx) !temporary only

C     Find out dependent variable basis and calculate NHE and NHP
        IF(ITYP1(nr,nx).EQ.3) THEN !FE30 problems
          IF(ITYP4(nr,nx).EQ.3) THEN !finite difference
            IF(ITYP3(nr,nx).EQ.1) THEN !flow in elastic tubes
              NHQ(nr)=6 ! radius,velocity,pressure
            ELSE
              NHQ(nr)=4 !defult value
            ENDIF
            IF(NHM.LT.NHQ(nr)) THEN
              WRITE(CHAR,'(I1)') IDIGITS(NHQ(nr))
              WRITE(ERROR,'(''>>Increase NHM to '',I'//CHAR//')')
     &          NHQ(nr)
              GO TO 9999
            ENDIF
          ELSE
C??? KAT 2007-01-30: Should this call be skipped for grid-based FEM/FV?
C???                 These elements don't really have dependent variable
C???                 basis functions.
            CALL IPBAS3(NBJ,NBH,NEELEM,NELIST,NHE,NHP,NPNE,
     '        NPNODE,nr,NW,nx,ERROR,*9999)
C KAT 28Apr99:  NHQ should not be used
C          NHQ(nr)=1 !not finite difference default
          ENDIF
        ELSE IF(ITYP1(nr,nx).EQ.4) THEN !FE40 problems
          CALL IPBAS4(NBJ,NBH,NEELEM,NHE,NHP,NPNE,NPNODE,nr,NW,nx,
     '      ERROR,*9999)
        ELSE IF(ITYP1(nr,nx).EQ.5) THEN !FE50 problems
          CALL IPBAS5(NBH,NBJ,NEELEM,NHE,NHP,NPNE,NPNODE,nr,nx,
     '      ERROR,*9999)
        ELSE IF(ITYP1(nr,nx).EQ.6) THEN !FE60 problems
          CALL IPBAS6(NBJ,NBH,NEELEM,NHE,NHP,
     '      NPNODE,nr,NW,nx,ERROR,*9999)
          KTYP24=1 ! Always store sparsely in compressed row format
        ELSE IF(ITYP1(nr,nx).EQ.9) THEN !FE90 problems
          IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.1) THEN
C         Static linear elasticity
            CALL IPBAS4(NBJ,NBH,NEELEM,NHE,NHP,NPNE,NPNODE,
     '        nr,NW,nx,ERROR,*9999)
          ELSE !rest of BEM problems
            CALL IPBAS9(IBT,IDO,INP,NBH,NBJ,NEELEM,NELIST,NHE,NHP,
     '        NKJE,NPF,NP_INTERFACE,NPNE,NPNODE,nr,NVJE,
     '        nx,CE,CURVCORRECT,SE,XA,XE,XP,ERROR,*9999)
          ENDIF
        ELSE
          ERROR='>>No ipbas# routine for this setup'
          GOTO 9999
        ENDIF !ityp1

C     Setup dependent version mapping arrays
      IF(ITYP1(nr,nx).EQ.3) THEN
C FE30 problems
        IF(ITYP4(nr,nx).NE.3) THEN
C excluding finite difference
          IF(ITYP5(nr,nx).EQ.1.AND !static analysis
     '      .ITYP2(nr,nx).EQ.3.AND !Laplace equation
     '      .ITYP3(nr,nx).EQ.3) THEN !aerofoil analysis
            CALL IPAERO(NEELEM,NLL,NPL,nr,ERROR,*9999)
            CALL CALC_VERSIONS_DEP(NBH,NBJ,NEELEM,NHE,NHP(1,nr),
     '        NPNE,NPNODE,nr,NVHE,NVHP,NVJE,NVJP,nx,ERROR,*9999)
          ELSE !non-aerofoil cases are prompted
C         Set up defaults
            CALL CALC_VERSIONS_DEP(NBH,NBJ,NEELEM,NHE,NHP(1,nr),
     '        NPNE,NPNODE,nr,NVHE,NVHP,NVJE,NVJP,nx,ERROR,*9999)
C         Prompt for special cases.
C         NVHE,NVHP are already set up for version mappings
            IF(JTYP2A.NE.1) THEN !
              DO nhx=1,NHE(NEELEM(1,nr)) !Note use of first ne of region
                nh=NH_LOC(nhx,nx)
                WRITE(CHAR1,'(I1)') nhx
                FORMAT='($,'' The max # of versions of variable '
     '            //CHAR1(1:1)//' is [1]: '',I1)'
                IF(IOTYPE.EQ.3) THEN
c news PJH 18Sep95
                  NT_versions(nh)=NVHP(nh,NPNODE(1,nr),1,nr)
c newe
                  IDATA(1)=NT_versions(nh)
                ENDIF
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     '            NVM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '            *9999)
                IF(IOTYPE.NE.3) THEN
                  NT_versions(nh)=IDATA(1)
                ENDIF

                IF(NT_versions(nh).GT.1) THEN
                  FORMAT='($'' Enter nodes with >1 version: '',20I5)'
                  IF(IOTYPE.EQ.3) THEN
                    n=0
                    DO nonode=1,NPNODE(0,nr)
                      np=NPNODE(nonode,nr)
                      IF(NVHP(nh,np,1,nr).GT.1) THEN
                        n=n+1
                        IDATA(n)=np
                      ENDIF
                      IDATA(0)=n
                    ENDDO
                  ENDIF
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,NPT(nr),LDATA,LDEFLT,RDATA,RDEFLT,
     &              RMIN,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    DO n=1,IDATA(0)
                      NVHP(nh,IDATA(n),1,nr)=NT_versions(nh)
                      NVHP(nh,IDATA(n),2,nr)=NT_versions(nh)
                    ENDDO
                  ENDIF
                ENDIF !nt_versions(nh)>1
              ENDDO !nh
C             Find max NVHP for zero'th nr location
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nc=1,NCT(nr,nx)
                  DO nhx=1,NH_LOC(0,nx)
                    nh=NH_LOC(nhx,nx)
                    IF(NVHP(nh,np,nc,nr).GT.NVHP(nh,np,nc,0))
     '            NVHP(nh,np,nc,0)=NVHP(nh,np,nc,nr)
                  ENDDO !nh
                ENDDO !nc
              ENDDO !nonode (np)
            ENDIF !JTYP2.NE.1
          ENDIF !aerofoil/others
        ENDIF !ityp4 not finite difference
      ELSE !all other problem types
        CALL CALC_VERSIONS_DEP(NBH,NBJ,NEELEM,NHE,NHP(1,nr),
     '    NPNE,NPNODE,nr,NVHE,NVHP,NVJE,NVJP,nx,ERROR,*9999)
      ENDIF !ityp1

C     Calculate NKH
C      IF(IOTYPE.NE.3) THEN
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov2002 grid finite volume also
      IF((ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9)
     '  .OR.(ITYP4(nr,nx).EQ.3).OR.(ITYP4(nr,nx).EQ.4)
     '  .OR.(ITYP4(nr,nx).EQ.6).OR.(ITYP4(nr,nx).EQ.7)) THEN
C Don't do for cardiac activation or finite difference or collocation or grid-based FE
      ELSE
C CPB 9/2/96 Initialise all NKH in this region
        CALL CALC_NKH(NBH,NEELEM,NHP(1,nr),NKH,NPNE,NPNODE,
     '    nr,NW,nx,ERROR,*9999)
      ENDIF
C      ENDIF

      IF(IOTYPE.NE.3) THEN

        IF(ITYP4(nr,nx).EQ.1.OR.ITYP4(nr,nx).EQ.2) THEN! FEM or BEM
C         Calculate NKHE
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
C             This loop ensures that everything is set up for all nc
              DO nc=1,NCT(nr,nx)
                nb=NBH(nh,nc,ne)
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    NKHE(nk,nn,nh,ne)=nk
                  ENDDO !nk
                ENDDO !nn
              ENDDO !nc
            ENDDO !nhx
          ENDDO !noelem
        ENDIF! FEM or BEM

C       Calculate ny's etc
        IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN
          NHQ(nr)=1
          CALL CALC_NY_GRID_DEP(NHQ,NQNY,nr,nx,NYNQ,
     '      NYQNR,ERROR,*9999)
        ELSE IF(ITYP4(nr,nx).EQ.3) THEN ! finite difference
          CALL CALC_NY_GRID_DEP(NHQ,NQNY,nr,nx,NYNQ,
     '      NYQNR,ERROR,*9999)
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov2002 grid finite volume also
        ELSE IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6
     '          .OR.ITYP4(nr,nx).EQ.7) THEN ! collocation or grid-based FE/FV
          NHQ(nr)=1
          CALL CALC_NY_GRID_DEP(NHQ,NQNY,nr,nx,NYNQ,
     '      NYQNR,ERROR,*9999)
        ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.11.AND.
     '      (ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.OR.
     &      ITYP3(nr,nx).EQ.6)) THEN
C KSB 24Oct2001 Pulmonary capillary blood flow
          CALL CALC_NY_MAPS_DEP2(NBH,NEELEM,NHP,NKH,NPNE,
     '      NPNODE,NPNY,nr,NVHP,nx,NYNE,NYNP,NYNR,ERROR,*9999)
        ELSE
          CALL CALC_NY_MAPS_DEP(NBH,NEELEM,NHP,NKH,NP_INTERFACE,
     '      NPNODE,NPNY,nr,NVHP,nx,NYNE,NYNP,NYNR,ERROR,*9999)
        ENDIF
      ENDIF

C Calculate dependent variable face basis information
      IF(IOTYPE.NE.3.AND.(ITYP4(nr,nx).EQ.1.OR.ITYP4(nr,nx).EQ.2)
     '  .AND.ITYP4(nr,nx).NE.5) THEN
        CALL CALC_FACE_BASIS_DEP(IBT,NBH,NBHF,NEELEM,NFF,NHE,NNF,
     '    nr,nx,ERROR,*9999)
      ENDIF


      IF(nr.EQ.NRLIST(NRLIST(0)).AND.(ITYP4(nr,nx).NE.3)
     '  .AND.(ITYP4(nr,nx).NE.5)) THEN
C       Finite volume sparsity etc not set up yet
C       only calc sparsity after ny's for last region have been set up

C cpb 28/10/98 Ask about assembling global matrices for ASSEMBLE1,
C ASSEMBLE2, SOVE1, SOLVE2, AND SOLVE9 type problems at the moment.

        IF((ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4).AND.
     '    ITYP6(nr,nx).EQ.1.AND.
     '    (ITYP4(nr,nx).EQ.1.OR.ITYP4(nr,nx).EQ.2)) THEN
C         Problem is static or quasi-static, linear and uses finite
C         elements or boundary elements.
          ADEFLT(1)='N'
          FORMAT='(/$,'' Do you want to assemble directly to the '
     '      //'solution matrices [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(CALC_GLOBAL(nr,nx)) THEN
              ADATA(1)='N'
            ELSE
              ADATA(1)='Y'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              DO nonrr=1,NRLIST(0)
                nrr=NRLIST(nonrr)
                CALC_GLOBAL(nrr,nx)=.FALSE.
              ENDDO
            ELSE
              DO nonrr=1,NRLIST(0)
                nrr=NRLIST(nonrr)
                CALC_GLOBAL(nrr,nx)=.TRUE.
              ENDDO
            ENDIF
          ENDIF
        ELSE
          DO nonrr=1,NRLIST(0)
            nrr=NRLIST(nonrr)
            CALC_GLOBAL(nrr,nx)=.TRUE.
          ENDDO
        ENDIF

        IF(ITYP4(nr,nx).EQ.2) THEN !BEM
          ADEFLT(1)='N'
        ELSE !all other cases
          ADEFLT(1)='Y'
        ENDIF
        FORMAT='(/$,'' Do you want the global matrices stored as '
     '    //'sparse matrices ['//ADEFLT(1)(1:1)//']? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(KTYP24.EQ.0) THEN
            ADATA(1)='N'
          ELSE
            ADATA(1)='Y'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          KTYP24=1
          IF(IOTYPE.NE.3) CALL ASSERT(USE_SPARSE.NE.0,
     '      '>>Set USE_SPARSE to 1 to use sparse matrices',ERROR,*9999)
        ELSE
          KTYP24=0
        ENDIF

C SGM 26 Oct 2000 grid-based Finite element also
C MLT 29Nov02 grid Finite Volume also
        IF(ITYP4(nr,nx).NE.4.AND.ITYP4(nr,nx).NE.6
     '     .AND.ITYP4(nr,nx).NE.7) THEN !not collocation or grid-based FE or FV
          CALC_SPARSITY=.FALSE.
          IF(KTYP24.EQ.1) THEN !using sparse matrix storage
            FORMAT='($,'' Do you want to calculate the sparsity '
     '        //'pattern for the global matrices [Y]? '',A)'
            ADEFLT(1)='Y'
            IF(IOTYPE.EQ.3) THEN
              ADATA(1)='Y'
              CALC_SPARSITY=.TRUE.
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(ADATA(1).EQ.'Y') THEN
              IF(IOTYPE.NE.3) THEN
                CALC_SPARSITY=.TRUE.
              ENDIF
            ELSE
              WRITE(OP_STRING,'('' >>Remember to "read matrix;<FILE> '
     '          //'sparsity <ascii/binary>"'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            IF(ITYP2(nr,nx).EQ.11)THEN
               !special case for large 1D lung trees
              LUNG_SPARSITY=.FALSE.
              FORMAT='($,'' Is the sparsity pattern calculated '
     '          //'directly [Y]? '',A)'
              ADEFLT(1)='Y'
              IF(IOTYPE.EQ.3) THEN
                ADATA(1)='Y'
                LUNG_SPARSITY=.TRUE.
              ENDIF 
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     &          ERROR,*9999)
              IF(ADATA(1).EQ.'Y') THEN
                IF(IOTYPE.NE.3) THEN
                  LUNG_SPARSITY=.TRUE.
                ENDIF
              ENDIF

              IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5)THEN
                IF(IOTYPE.EQ.3) THEN
                  IDATA(1)=NJ_TYPE(nj_radius,2)
                ENDIF
                IDEFLT(1)=1
                FORMAT='(/$,'' The field number for the non-alveolar '
     &            //'radius is [1]: '',I2)'
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &            0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &            *9999)
                CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &            '>>Define field first',ERROR,*9999)
                IF(IOTYPE.NE.3)THEN
                  njj=IDATA(1)
                  nj_radius=NJ_LOC(NJL_FIEL,njj,nr)
                ENDIF
                IF(ITYP3(nr,nx).EQ.1)THEN ! advection diffusion
                  IF(IOTYPE.EQ.3) THEN
                    IDATA(1)=NJ_TYPE(nj_alveoli,2)
                  ENDIF
                  IDEFLT(1)=2
                  FORMAT='($,'' The field number for the duct:total '
     &              //'radius is [2]: '',I2)'
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                  CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &              '>>Define field first',ERROR,*9999)
                  IF(IOTYPE.NE.3)THEN
                    njj=IDATA(1)
                    nj_alveoli=NJ_LOC(NJL_FIEL,njj,nr)
                  ENDIF

                  IF(ITYP7(nr,nx).GT.1)THEN !with gas exchange
                    IDEFLT(1)=njj+1
                    WRITE(CHAR_2,'(I2)') IDEFLT(1)
                    FORMAT='($,'' The field number for the acinar '
     &                //'ventilation is ['//CHAR//']: '',I2)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                    CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &                '>>Define field first',ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      njj=IDATA(1)
                      nj_flow=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF

                    IDEFLT(1)=njj+1
                    WRITE(CHAR_2,'(I2)') IDEFLT(1)
                    FORMAT='($,'' The field number for the acinar '
     &                //'perfusion is ['//CHAR//']: '',I2)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                    CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &                '>>Define field first',ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      njj=IDATA(1)
                      nj_Qdot=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF

                    IDEFLT(1)=njj+1
                    WRITE(CHAR_2,'(I2)') IDEFLT(1)
                    FORMAT='($,'' The field number for the acinar '
     &                //'blood volume is ['//CHAR//']: '',I2)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                    CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &                '>>Define field first',ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      njj=IDATA(1)
                      nj_Vc=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF 

                    IDEFLT(1)=njj+1
                    WRITE(CHAR_2,'(I2)') IDEFLT(1)
                    FORMAT='($,'' The field number for the acinar '
     &                //'RBC transit time is ['//CHAR//']: '',I2)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                    CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &                '>>Define field first',ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      njj=IDATA(1)
                      nj_tt=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF

                    IDEFLT(1)=njj+1
                    WRITE(CHAR_2,'(I2)') IDEFLT(1)
                    FORMAT='($,'' The field number for the blood '
     &                //'partial pressure is ['//CHAR//']: '',I2)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                    CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &                '>>Define field first',ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      njj=IDATA(1)
                      nj_pob=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF

                    IDEFLT(1)=njj+1
                    WRITE(CHAR_2,'(I2)') IDEFLT(1)
                    FORMAT='($,'' The field number for the air '
     &                //'partial pressure is ['//CHAR//']: '',I2)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                    CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &                '>>Define field first',ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      njj=IDATA(1)
                      nj_poa=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF

                    IDEFLT(1)=njj+1
                    WRITE(CHAR_2,'(I2)') IDEFLT(1)
                    FORMAT='($,'' The field number for the source '
     &                //'term is ['//CHAR//']: '',I2)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                    CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &                '>>Define field first',ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      njj=IDATA(1)
                      nj_source=NJ_LOC(NJL_FIEL,njj,nr)
                    ENDIF
                  ENDIF !gas exchange

                ELSE IF(ITYP3(nr,nx).EQ.5)THEN !water and heat transfer
                  IF(IOTYPE.EQ.3) THEN
                    IDATA(1)=NJ_TYPE(nj_flow,2)
                  ENDIF
                  IDEFLT(1)=2
                  WRITE(CHAR,'(I1)') IDEFLT(1)
                  FORMAT='($,'' The field number for the mucus '
     &              //'velocity is ['//CHAR//']: '',I2)'
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                  CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &              '>>Define field first',ERROR,*9999)
                  IF(IOTYPE.NE.3)THEN
                    njj=IDATA(1)
                    nj_flow=NJ_LOC(NJL_FIEL,njj,nr)
                  ENDIF
                  IF(IOTYPE.EQ.3) THEN
                    IDATA(1)=NJ_TYPE(nj_source,2)
                  ENDIF
                  IDEFLT(1)=3
                  WRITE(CHAR,'(I1)') IDEFLT(1)
                  FORMAT='($,'' The field number for the mucus '
     &              //'source is ['//CHAR//']: '',I2)'
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,NJM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                  CALL ASSERT(NJ_LOC(NJL_FIEL,IDATA(1),nr).GT.0,
     &              '>>Define field first',ERROR,*9999)
                  IF(IOTYPE.NE.3)THEN
                    njj=IDATA(1)
                    nj_source=NJ_LOC(NJL_FIEL,njj,nr)
                  ENDIF
                ENDIF !ITYP3.EQ.1
              ENDIF !ITYP3.LE.2
            ENDIF !LUNG
          ENDIF

C CPB 6/11/95 Temporary work array allocation
          IF(ITYP2(nr,nx).EQ.11.AND.ITYP3(nr,nx).LE.6)THEN
            !special case for large 1D lung trees & capillary network
            IF(KTYP24.NE.0)THEN
                WORK_ISC_PTR=0
                WORK_ISR_PTR=0
              IF(.NOT.LUNG_SPARSITY)THEN
                CALL ALLOCATE_MEMORY(NYT(1,1,nx)*8,1,INTTYPE,
     '            WORK_ISC_PTR,MEM_INIT,ERROR,*9999)
                CALL ALLOCATE_MEMORY(NYT(1,1,nx)*8,1,INTTYPE,
     '            WORK_ISR_PTR,MEM_INIT,ERROR,*9999)
              ENDIF
              CALL CALC_SPARSE_SOLVE_1DTREE(NISC_GKM,NISR_GKM,ISC_GK,
     &          ISR_GK,NYT(1,1,nx),NYT(2,1,nx),NBJ,NEELEM,NENP,NHE(1),
     &          NPNE,NPNODE,nr,nx,NXI,NYNE,NYNP,NZ_GK_M,NZT(1,nx),
     &          %VAL(WORK_ISC_PTR),%VAL(WORK_ISR_PTR),ERROR,*9999) 
              IF(.NOT.LUNG_SPARSITY)THEN
                CALL FREE_MEMORY(WORK_ISC_PTR,ERROR,*9999)
                CALL FREE_MEMORY(WORK_ISR_PTR,ERROR,*9999)
              ENDIF
            ENDIF !pulmonary

          ELSEIF(ITYP4(nr,nx).NE.5) THEN !Sparsity set up for finite
            IF(KTYP24.NE.0) THEN         !volumes in fe60
              WORK_PTR=0
              CALL ALLOCATE_MEMORY(NYT(1,1,nx)*NYT(2,1,nx),1,CHARTYPE,
     '          WORK_PTR,MEM_INIT,ERROR,*9999)
            ENDIF
            IF(CALC_GLOBAL(nr,nx)) THEN
              CALL CALC_SPARSE_SOLVE(NISC_GKM,NISR_GKM,ISC_GK,ISR_GK,
     '          LGE,NYT(1,1,nx),NYT(2,1,nx),NBH,1,NEELEM,NHE,NPNE,
     '          NPNY,NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZ_GK_M,KTYP24,
     '          %VAL(WORK_PTR),FIX,CALC_SPARSITY,ERROR,*9999)
            ENDIF
          ENDIF
          IF(ITYP2(nr,nx).EQ.11.AND.ITYP3(nr,nx).LE.6)THEN
            !special case for large 1D lung trees, don't free memory
          ELSE
            IF(KTYP24.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
          ENDIF
          BEMREGION=.FALSE.
          DO nonrr=1,NRLIST(0)
            nrr=NRLIST(nonrr)
            IF(ITYP4(nrr,nx).EQ.2) BEMREGION=.TRUE.
          ENDDO
          IF(BEMREGION) THEN !BEM
C CPB 6/11/95 Temporary work array allocation
            IF(CALC_GLOBAL(nr,nx)) THEN
              IF(KTYP24.NE.0) THEN
                WORK_PTR=0
                CALL ALLOCATE_MEMORY(NYT(1,2,nx)*NYT(2,2,nx),1,CHARTYPE,
     '            WORK_PTR,MEM_INIT,ERROR,*9999)
              ENDIF
              CALL CALC_SPARSE_SOLVE(NISC_GQM,NISR_GQM,ISC_GQ,ISR_GQ,
     '          LGE,NYT(1,2,nx),NYT(2,2,nx),NBH,2,NEELEM,NHE,NPNE,
     '          NPNY,NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZ_GQ_M,KTYP24,
     '          %VAL(WORK_PTR),FIX,CALC_SPARSITY,ERROR,*9999)
              IF(KTYP24.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      ENDIF !IF(ITYP19(nr,nx).EQ.2.AND.ITYP3(nr,nx).NE.4)

C     JWF 6.3.2002 Contact Question
C     Currently setup for 3D Static finite or linear elasticity with Galerkin FEM.

C     Initialise KTYP and viscous term penalty (VIS_WT)
      KTYP5H(nr)=0 
      KTYP5J(nr)=0                  
      KTYP5K(nr)=0
      VIS_WT=0.0d0 

      IF(ITYP5(nr,nx).EQ.1.AND. !static analysis
C     &  ! KTYP51(nr).EQ.3 for 3D in finite elasticity
     &  ((ITYP2(nr,nx).EQ.2.AND.KTYP51(nr).EQ.3).OR.
C     &  ! and temporarily using LAST_TYPE.EQ.9 for 3D in linear elasticity
     &  (ITYP2(nr,nx).EQ.1.AND.LAST_TYPE.EQ.9)).AND.
     &  ITYP4(nr,nx).EQ.1) THEN !Galerkin finite elements
     
C       Set up NCT for the coupled region in contact problem
        NCT(0,nx)=2
        FORMAT='($,'' Does this problem involve contact'
     &    //' mechanics [N]? '',A)'
        ADEFLT(1)='N'
        IF(IOTYPE.EQ.3) THEN
          IF (KTYP5G(nr).GT.0) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &    IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &    *9999)
        IF(ADATA(1).EQ.'Y') THEN
          KTYP5G(nr)=1
C         Setup contact specific variables defined in nonl00.cmn
          CONT_IT=1
          AUG_IT=0
          LOAD_IT=1
          CONV_IT=1  
        ELSE
          KTYP5G(nr)=0
          CONT_IT=0
          AUG_IT=0
          LOAD_IT=0
          CONV_IT=0  
        ENDIF                   
C*** 28/02/08 JHC Initialise Augment
        AUGMENT=0 

        IF(KTYP5G(nr).EQ.1) THEN ! contact problem questions

C*** 30/05/08 JHC Add warning not to use sparse matrix for contact problems
          WRITE(OP_STRING,'('' >>Warning: Do not use sparse matrix '
     '      //'option in ipequa and ipsolv/irsolv files.'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          FORMAT='($,'' Is this body a host mesh [N]? '',A)'
          ADEFLT(1)='N'
          IF(IOTYPE.EQ.3) THEN
            IF (KTYP5H(nr).GT.0) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &      IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &      *9999)
          IF(ADATA(1).EQ.'Y') THEN
            KTYP5H(nr)=1
          ELSE
            KTYP5H(nr)=0
          ENDIF
                     
          FORMAT='($,'' Do you want the contact problem'
     &    //' to start with active constraints [N]? '',A)'
          ADEFLT(1)='N'
          IF(IOTYPE.EQ.3) THEN
            IF (KTYP5J(nr).GT.0) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &      IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &      *9999)
          IF(ADATA(1).EQ.'Y') THEN
            KTYP5J(nr)=1
          ELSE
            KTYP5J(nr)=0
          ENDIF
                       
          FORMAT='($,'' Do you want the contact problem'
     &    //' to include a psuedo-viscous term [N]? '',A)'
          ADEFLT(1)='N'
          IF(IOTYPE.EQ.3) THEN
            IF (KTYP5K(nr).GT.0) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &      IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &      *9999)
          IF(ADATA(1).EQ.'Y') THEN
            KTYP5K(nr)=1
          ELSE
            KTYP5K(nr)=0
          ENDIF
                  
          IF(KTYP5K(nr).EQ.1) THEN ! get penalty value
            RDEFLT(1)=0.0d0
            FORMAT='($,'' Enter the penalty weight [0.0]: '',D5.2)'
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=0.0d0
              VIS_WT=RDATA(1)
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              VIS_WT=RDATA(1)
            ENDIF
          ENDIF ! set penalty
                          
        ENDIF ! contact problem questions                                                             
                        
      ENDIF !static 3D elasticity, galerkin FEM
      CALL EXITS('IPEQUA')
      RETURN
 9999 CALL ERRORS('IPEQUA',ERROR)
      CALL EXITS('IPEQUA')
      RETURN 1
      END

