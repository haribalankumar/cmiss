      SUBROUTINE IPMESH11(IBT,NBJ,nb_face,NEELEM,NENFVC,NENP,
     '  NFVC,NKJ,NKJE,NODENVC,NODENVCB,NP_INTERFACE,NPLIST,
     '  NPNE,NPNODE,nr,NRE,NVCNODE,NVJE,NVJP,NXI,SE,VC,VC_INIT,
     '  XNFV,XP,ZA,ERROR,*)

C#### Subroutine: IPMESH11
C###  Description:
C###    IPMESH11 inputs and calculates the voronoi mesh

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),nb_face,NEELEM(0:NE_R_M,0:NRM),
     '  NENFVC(0:NFVCM,NFVM),NENP(NPM,0:NEPM,0:NRM),
     '  NFVC(2,0:NFVCM,NVCM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NODENVC(NVCM),NODENVCB(NVCBM),NP_INTERFACE(0:NPM,0:3),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NRE(NEM),NVCNODE(2,NP_R_M),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM),IBT(3,NIM,NBFM)
      REAL*8 SE(NSM,NBFM,NEM),VC(0:NVCM),VC_INIT(2,NVCM),
     '  XNFV(-(NJM+1):NJM,NFVM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ICHAR,INFO,n,NGROUP,NOQUES,NVIT,N_BDRY,N_IBDRY,N_INTNL
      LOGICAL FILEIP

      CALL ENTERS('IPMESH11',*9999)

C     ..Early exit if possible
      CALL ASSERT(USE_VORONOI.EQ.1,'>>Set USE_VORONOI = 1',ERROR,*9999)

C     ..Initialise Variables
      ICHAR=999
      NOQUES=0
      FILEIP=.FALSE.

C***  Voronoi cells are created from Delaunay triangulation
C     ..Initialisation
      NGROUP=1
      IDEFLT(1)=NGROUP
      IF(IOTYPE.EQ.3) IDATA(1)=NGROUP
      FORMAT='($,'' The number of node groups for the '//
     '  'Voronoi boundary (B) type is [1]: '',I2)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '  IDATA,IDEFLT,1,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '  RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NGROUP=IDATA(1)
      NVBT=0
      DO i=1,NGROUP
        FORMAT='($,'' Enter the node numbers of the '//
     '    'Voronoi boundary (B) type: '',I5)'
        CDATA(1)='NODES' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '    IDATA,IONE,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '    RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NVBT=NVBT+IDATA(0)
          N_BDRY=NVBT
          DO n=1,IDATA(0)
            NPLIST(n)=IDATA(n)
          ENDDO !n
        ENDIF !iotype.ne.3
      ENDDO !ngroup

      NGROUP=1
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=NGROUP
      FORMAT='($,'' The number of node groups for the '//
     '  'Voronoi internal boundary (IB) type is [1]: '',I2)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '  IDATA,IDEFLT,1,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '  RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NGROUP=IDATA(1)
      NVIBT=0
      DO i=1,NGROUP
        FORMAT='($,'' Enter the node numbers of the '//
     '    'Voronoi internal boundary (IB) type: '',I5)'
        CDATA(1)='NODES' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '    IDATA,IONE,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '    RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NVIBT=NVIBT+IDATA(0)
          N_IBDRY=NVIBT
          DO n=1,IDATA(0)
            NPLIST(NVBT+n)=IDATA(n)
          ENDDO !n
        ENDIF !iotype.ne.3
      ENDDO

      NGROUP=1
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=NGROUP
      FORMAT='($,'' The number of node groups for the '//
     '  'Voronoi internal (I) type is [1]: '',I2)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '  FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '  IDATA,IDEFLT,0,100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '  RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NGROUP=IDATA(1)
      NVIT=0
      DO i=1,NGROUP
        FORMAT='($,'' Enter the node numbers of the '//
     '    'Voronoi internal (I) type: '',I5)'
        CDATA(1)='NODES' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '    IDATA,IONE,1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '    RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NVIT=NVIT+IDATA(0)
          N_INTNL=NVIT
          DO n=1,IDATA(0)
            NPLIST(NVBT+NVIBT+n)=IDATA(n)
          ENDDO !n
        ENDIF !iotype.ne.3
      ENDDO

      CALL VORONOI(IBT,nb_face,N_BDRY,N_IBDRY,N_INTNL,NBJ,NEELEM,
     '  NENFVC,NENP,NFVC,NKJ,NKJE,1,NODENVC,NODENVCB,
     '  NP_INTERFACE,NPLIST,NPNE,NPNODE,nr,NRE,1,2,NVCNODE,NVJE,
     '  NVJP,NXI,SE,VC,VC_INIT,XNFV,XP,ZA,.FALSE.,ERROR,*9999)

      CALL EXITS('IPMESH11')
      RETURN
 9999 CALL ERRORS('IPMESH11',ERROR)
      CALL EXITS('IPMESH11')
      RETURN 1
      END


