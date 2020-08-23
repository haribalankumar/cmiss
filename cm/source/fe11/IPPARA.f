      SUBROUTINE IPPARA(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,ISIZE_MFI,
     '  ISIZE_TBH,MINIMAL,NDIPOLES,NDLT,NEELEM,NEL,NENP,
     '  NFFACE,NLLINE,NONY,NPNODE,NQET,NQXI,NVHP,NVJP,NYNO,YP,ERROR,*)

C#### Subroutine: IPPARA
C###  Description:
C###    IPPARA defines array dimension parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'maqloc00.cmn'
      INCLUDE 'nqloc00.cmn'
C      INCLUDE 'cmiss$reference:mesh00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'opti00.cmn'
      INCLUDE 'ptr01.cmn'
      INCLUDE 'quas00.cmn'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),ISIZE_MFI(3,NSSM),
     '  ISIZE_TBH(2),NDLT(NEM),NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),NFFACE(0:NF_R_M,NRM),
     '  NLLINE(0:NL_R_M,0:NRM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),NQXI(0:NIM,NQSCM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM)
      REAL*8 YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NR_MX
      PARAMETER (NR_MX=9)
      INTEGER i,ICHAR,INFO,N,NOQUES,nb,nc,ne,nh,nhx,NITMAX,nj,njj1,njj2,
     '  nl,NMAX,NMAX2,no,nonode,np,nq,nqsc,nr,nrc,nx,nxx,nw,ny,TOTAL
      LOGICAL FILEIP,MINIMAL,FOUND1,FOUND2

      CALL ENTERS('IPPARA',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C#### Variable: NAM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NAM is the maximum number of auxiliary parameters.

      FORMAT='($,'' Max# auxiliary parameters          (NAM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nb=1,NBT
            IF(NAT(nb).GT.NMAX) NMAX=NAT(nb)
          ENDDO !nb
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
          IF(CALL_GRID.AND.(NMAX.LT.3)) IDATA(1)=3
            !temp until YQ is sorted better
        ELSE
          IDATA(1)=NAM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NAM=IDATA(1)

C#### Variable: NBM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NBM is the maximum number of basis functions.

      FORMAT='($,'' Max# basis functions               (NBM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NBT.GT.0) THEN
            IDATA(1)=NBT
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NBM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NBM=IDATA(1)

C#### Variable: NCM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NCM is the maximum number of variable types for a given
C###    dependent variable.  Each variable type typically corresponds to
C###    a different derivative of the dependent variable, e.g. velocity
C###    and acceleration are different types of the dependent variable
C###    displacement.

      FORMAT='($,'' Max# var. types for a  dep. var.   (NCM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            TOTAL=0
            DO nr=1,NRT
              IF(NCT(nr,nx).GT.TOTAL) TOTAL=NCT(nr,nx)
C              TOTAL=TOTAL+NCT(nr,nx)
            ENDDO !nr
            IF(TOTAL.GT.NMAX) NMAX=TOTAL
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NCM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NCM=IDATA(1)

C#### Variable: NDM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NDM is the maximum number of data points.

      FORMAT='($,'' Max# data points                   (NDM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NDT.GT.0) THEN
            IDATA(1)=NDT
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NDM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NDM=IDATA(1)

C#### Variable: NEM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NEM is the maximum number of elements.

      FORMAT='($,'' Max# elements                      (NEM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NET(0).GT.0) THEN
            IDATA(1)=NET(0)
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NEM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NEM=IDATA(1)

C#### Variable: NE_R_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NE_R_M is the maximum number of elements in a region.  This is
C###    less than or equal to the maximum number of elements (NEM).

      FORMAT='($,'' Max# elements in a region       (NE_R_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            IF(NEELEM(0,nr).GT.NMAX) NMAX=NEELEM(0,nr)
          ENDDO !nr
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NE_R_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NE_R_M=IDATA(1)

C#### Variable: NEIM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NEIM is the maximum number of adjacent elements in an
C###    Xi direction.  This will usually be 1.

      FORMAT='($,'' Max# adjacent elements in Xi      (NEIM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IDATA(1)=1
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NEIM=IDATA(1)

C#### Variable: NFM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NFM is the maximum number of global face segments.

      FORMAT='($,'' Max# global face segments          (NFM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NFT.GT.0) THEN
            IDATA(1)=NFT
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NFM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NFM=IDATA(1)

C#### Variable: NF_R_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NF_R_M is the maximum number of faces in a region.  This is less
C###    than or equal to the maximum number of global face segments
C###    (NFM).

      FORMAT='($,'' Max# faces in a region          (NF_R_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            IF(NFFACE(0,nr).GT.NMAX) NMAX=NFFACE(0,nr)
          ENDDO !nr
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NF_R_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NF_R_M=IDATA(1)

C#### Variable: NFVCM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NFVCM is the maximum number of local Voronoi faces.

      FORMAT='($,'' Max# local Voronoi faces         (NFVCM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NFVCM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NFVCM=IDATA(1)

C#### Variable: NGM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NGM is the maximum number of Gauss points per element.

      FORMAT='($,'' Max# Gauss points per element      (NGM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nb=1,NBT
            IF(NGT(nb).GT.NMAX) NMAX=NGT(nb)
          ENDDO !nb
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NGM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NGM=IDATA(1)
        NGFM=IDNINT(NGM**(2.0d0/3))
      ENDIF

C#### Variable: NHM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NHM is the maximum number of dependent variables.  NHM must be
C###    equal to NJM for problems where the dependent variable array
C###    carries deformed coordinates (e.g. finite elasticity problems).

      FORMAT='($,'' Max# dependent variables           (NHM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NH_LOC(0,0).GT.0) THEN
            IDATA(1)=NH_LOC(0,0)
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NHM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NHM=IDATA(1)

C#### Variable: NIM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIM is the maximum number of local Xi-coordinates.

      FORMAT='($,'' Max# local Xi coordinates          (NIM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nb=1,NBT
            IF(NIT(nb).GT.NMAX) NMAX=NIT(nb)
          ENDDO !nb
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NIM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIM=IDATA(1)

C#### Variable: NJM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NJM is the maximum number of independent variables (geometry,
C###    field, and fibre).

      FORMAT='($,'' Max# global reference coordinates  (NJM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NJ_LOC(0,0,0).GT.0) THEN
            IDATA(1)=NJ_LOC(0,0,0)
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NJM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NJM=IDATA(1)

C#### Variable: NKM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NKM is the maximum number of derivatives per variable.

      FORMAT='($,'' Max# derivatives per variable      (NKM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nb=1,NBT
            IF(NKT(0,nb).GT.NMAX) NMAX=NKT(0,nb)
          ENDDO !nb
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NKM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NKM=IDATA(1)

C#### Variable: NLM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NLM is the maximum number of global line segments.

      FORMAT='($,'' Max# global line segments          (NLM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NLT.GT.0) THEN
            IDATA(1)=NLT
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NLM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NLM=IDATA(1)

C#### Variable: NL_R_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NL_R_M is the maximum number of lines in a region.  This is less
C###    than or equal to the maximum number of global line segments
C###    (NLM).

      FORMAT='($,'' Max# lines in a region          (NL_R_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            IF(NLLINE(0,nr).GT.NMAX) NMAX=NLLINE(0,nr)
          ENDDO !nr
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NL_R_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NL_R_M=IDATA(1)

C#### Variable: NMM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NMM is the maximum number of material parameters.

      FORMAT='($,'' Max# material parameters           (NMM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NMM.GT.0) THEN
            IDATA(1)=NMM
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NMM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NMM=IDATA(1)

C#### Variable: NNM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NNM is the maximum number of element nodes.

      FORMAT='($,'' Max# element nodes                 (NNM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nb=1,NBT
            IF(NNT(nb).GT.NMAX) NMAX=NNT(nb)
          ENDDO !nb
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NNM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NNM=IDATA(1)

C#### Variable: NOM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NOM is the maximum number of solution degrees of freedom.  For
C###    FEM and BEM problems the number of solution degrees of freedom
C###    is less than or equal to the number of mesh degrees of freedom.

      FORMAT='($,'' Max# degrees of freedom            (NOM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nx=1,NXM
            DO nr=0,NRT
              DO nrc=1,NRCM
                DO nc=1,4
C                  IF(NOT(nrc,1,nr,nx).GT.NMAX) NMAX=NOT(nrc,1,nr,nx)
C                  IF(NOT(nrc,4,nr,nx).GT.NMAX) NMAX=NOT(nrc,4,nr,nx)
                  IF(NOT(nrc,nc,nr,nx).GT.NMAX) NMAX=NOT(nrc,nc,nr,nx)
                ENDDO
              ENDDO !nrc
            ENDDO !nr
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NOM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NOM=IDATA(1)

C#### Variable: NPM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NPM is the maximum number of global nodes.

      FORMAT='($,'' Max# global nodes                  (NPM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NPT(0).GT.0) THEN
            IDATA(1)=NPT(0)
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NPM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NPM=IDATA(1)

C#### Variable: NP_R_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NP_R_M is the maximum number of global nodes in a region.  This
C###    is less than or equal to the maximum number of global nodes.

      FORMAT='($,'' Max# global nodes in a region   (NP_R_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            IF(NPNODE(0,nr).GT.NMAX) NMAX=NPNODE(0,nr)
          ENDDO !nr
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NP_R_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NP_R_M=IDATA(1)

C#### Variable: NQM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQM is the maximum number of global grid points.

      FORMAT='($,'' Max# global grid points            (NQM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NQT.GT.0) THEN
            IDATA(1)=NQT
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NQM
        ENDIF
      ENDIF

      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQM=IDATA(1)

C#### Variable: NSSM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NSSM is the maximum number of signal sets.
C###  See-Also: ISIZE_NSSM, ISIZE_PHI, MFI, PHI

      FORMAT='($,'' Max# signal sets                  (NSSM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IDATA(1)=1
          DO i=1,NSSM
            IF(ISIZE_MFI(1,i).GT.0) THEN !has electrodes
              IDATA(1)=i
            ENDIF
          ENDDO
        ELSE
          IDATA(1)=NSSM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NSSM=IDATA(1)

C#### Variable: NYQM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NYQM is the maximum number of grid dofs.

      FORMAT='($,'' Max# grid degrees of freedom      (NYQM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nx=1,NXM
            DO nr=0,NRT
              DO nrc=1,NRCM
                DO nc=1,4
                  IF(NOQT(nrc,nc,nr,nx).GT.NMAX) NMAX=NOQT(nrc,nc,nr,nx)
                ENDDO
              ENDDO !nrc
            ENDDO !nr
          ENDDO !nx
          IF(NMAX.LT.NQT) NMAX=NQT
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NYQM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NYQM=IDATA(1)

C#### Variable: NRM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NRM is the maximum number of regions and contours per element
C###    and variable segments.

      FORMAT='($,'' Max# regions                       (NRM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NRT.GT.0) THEN
            IDATA(1)=NRT
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NRM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NR_MX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
       IF(IOTYPE.NE.3) NRM=IDATA(1)

C#### Variable: NSM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NSM is the maximum number of element degrees of freedom per
C###    variable. This is less than or equal to the maximum number of
C###    element nodes multiplied by the maximum number of derivatives
C###    per variable (NNM*NKM).

      FORMAT='($,'' Max# element dofs per variable     (NSM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nb=1,NBT
            IF(NST(nb).GT.NMAX) NMAX=NST(nb)
          ENDDO !nb
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NSM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NSM=IDATA(1)
      ENDIF

C#### Variable: NSFM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NSFM is the maximum number of face dofs per variable.

      FORMAT='($,'' Max# face dofs per variable       (NSFM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=MIN(16,NSM-NSM/2)
          IDATA(1)=NMAX
        ELSE
          IDATA(1)=NSFM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NSFM=IDATA(1)
      ENDIF

C#### Variable: NTM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NTM is the maximum number of eigenvalues.

      FORMAT='($,'' Max# eigenvalues                   (NTM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NTM=KTYP17
          IF(NTM.GT.0) THEN
            IDATA(1)=NTM
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NTM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NTM=IDATA(1)

C#### Variable: NTSM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NTSM is the maximum number of time samples.

      FORMAT='($,'' Max# time samples                 (NTSM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NTSM.GT.0) THEN
            IDATA(1)=NTSM
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NTSM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NTSM=IDATA(1)

C#### Variable: NUM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NUM is the maximum number of derivative terms up to the 2nd
C###    order.
C###  See-Also: nu

      FORMAT='($,'' Max# derivatives up to 2nd order   (NUM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nb=1,NBT
            IF(NUT(nb).GT.NMAX) NMAX=NUT(nb)
          ENDDO !nb
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NUM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NUM=IDATA(1)

C#### Variable: NVCBM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NVCBM is the maximum number of Voronoi boundary nodes.

      FORMAT='($,'' Max# Voronoi boundary nodes      (NVCBM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NVCBM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NVCBM=IDATA(1)

C#### Variable: NVCM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NVCM is the maximum number of Voronoi cells.

      FORMAT='($,'' Max# Voronoi cells                (NVCM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NVCM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NVCM=IDATA(1)

C#### Variable: NVM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NVM is the maximum number of versions of a particular value.
C###    For example a derivative at a node as seen by one element may be
C###    different to that seen by another element. The different values
C###    seen by the different elements are called different versions in
C###    CMISS.

      FORMAT='($,'' Max# versions of a variable        (NVM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=0,NRT
            DO np=1,NPT(nr)
              DO njj1=1,3 !geometry,fibres,fields
                DO njj2=1,NJ_LOC(njj1,0,nr)
                  nj=NJ_LOC(njj1,njj2,nr)
                  IF(nj.GT.0) THEN
                    IF(NVJP(nj,np).GT.NMAX) NMAX=NVJP(nj,np)
                  ENDIF
                ENDDO !njj2
              ENDDO !njj1
              DO nxx=1,NX_LIST(0)
                nx=NX_LIST(nxx)
                DO nc=1,NCT(nr,nx)
                  DO nhx=1,NH_LOC(0,nx)
                    nh=NH_LOC(nhx,nx)
                    IF(NVHP(nh,np,nc,nr).GT.NMAX) NMAX=NVHP(nh,np,nc,nr)
                  ENDDO !nhx
                ENDDO !nc
              ENDDO !nx
            ENDDO !np
          ENDDO !nr
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NVM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NVM=IDATA(1)

C#### Variable: NWM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NWM is the maximum number of workstations for segment arrays.
C###    The use of the word workstation is a hangover from the VAX days,
C###    it may be more appropriately thought of as WINDOWS. For example
C###    to view a three dimensional object using GX, three workstations
C###    (windows) are opened, one viewing the XY-plane, one YZ-plane and
C###    another the ZX-plane.

      FORMAT='($,'' Max# workstations                  (NWM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          TOTAL=0
          DO nw=1,NWM
            IF(IWKS(nw).GT.0) TOTAL=TOTAL+1
          ENDDO
          IF(TOTAL.GT.0) THEN
            IDATA(1)=TOTAL
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NWM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NWM=IDATA(1)

C#### Variable: NXM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NXM is the maximum number of problem types.

      FORMAT='($,'' Max# problem types                 (NXM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NX_LIST(0).GT.0) THEN
            IDATA(1)=NX_LIST(0)
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NXM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NXM=IDATA(1)

C#### Variable: NYM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NYM is the maximum number of mesh degrees of freedom.

      FORMAT='($,'' Max# mesh dofs                     (NYM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
C          DO nxx=1,NX_LIST(0)
C            nx=NX_LIST(nxx)
C            TOTAL=0
C            DO nr=1,NRT
C              DO nc=1,NCT(nr,nx)
C                DO nrc=1,2
C                  TOTAL=TOTAL+NYT(nrc,nc,nx)
C                ENDDO
C              ENDDO !nc
C            ENDDO !nr
C            IF(TOTAL.GT.NMAX) NMAX=TOTAL
C          ENDDO !nx
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nrc=1,2
              DO nr=1,NRT
                TOTAL=0
                DO nc=1,NCT(nr,nx)
                  TOTAL=TOTAL+NYT(nrc,nc,nx)
                ENDDO
                IF(TOTAL.GT.NMAX) NMAX=TOTAL
              ENDDO
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NYM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NYM=IDATA(1)

C#### Variable: NY_R_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NY_R_M is the maximum number of mesh dofs in a
C###    region. This is less than or equal to the maximum number of mesh
C###    degrees of freedom (NYM).

      FORMAT='($,'' Max# mesh dofs in a region      (NY_R_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
C          DO nxx=1,NX_LIST(0)
C            nx=NX_LIST(nxx)
C            DO nr=1,NRT
C              TOTAL=0
C              DO nc=1,NCT(nr,nx)
C                DO nrc=1,2
C                  TOTAL=TOTAL+NYT(nrc,nc,nx)
C                ENDDO
C              ENDDO !nc
C              IF(TOTAL.GT.NMAX) NMAX=TOTAL
C            ENDDO !nr
C          ENDDO !nx
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              DO nc=1,NCT(nr,nx)
                DO nrc=1,2
                  IF(NYT(nrc,nc,nx).GT.NMAX) NMAX=NYT(nrc,nc,nx)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NY_R_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NY_R_M=IDATA(1)

C#### Variable: NZ_GD_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NZ_GD_M is the maximum number of nonzero elements in GD.
C###  See-Also: GD,SPARSITY STRUCTURES

      FORMAT='($,'' Max# dimension of GD           (NZ_GD_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            IF(NZT(3,nx).GT.NMAX) NMAX=NZT(3,nx)
            IF(NYT(1,3,nx).GT.NMAX) NMAX=NYT(1,3,nx)
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NZ_GD_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NZ_GD_M=IDATA(1)

C#### Variable: NZ_GK_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NZ_GK_M is the maximum number of nonzero elements in GK.
C###  See-Also: GK,SPARSITY STRUCTURES

      FORMAT='($,'' Max# dimension of GK           (NZ_GK_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            IF(NZT(1,nx).GT.NMAX) NMAX=NZT(1,nx)
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NZ_GK_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NZ_GK_M=IDATA(1)

C#### Variable: NZ_GKK_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NZ_GKK_M is the maximum number of nonzero elements in GKK.
C###  See-Also: GKK,SPARSITY STRUCTURES

      FORMAT='($,'' Max# dimension of GKK         (NZ_GKK_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=0,NRT
              IF(NZZT(1,nr,nx).GT.NMAX) NMAX=NZZT(1,nr,nx)
            ENDDO
            IF(SOLVELU_NZA(nx).GT.NMAX) NMAX=SOLVELU_NZA(nx)
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NZ_GKK_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NZ_GKK_M=IDATA(1)

C#### Variable: NZ_GM_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NZ_GM_M is the maximum number of nonzero elements in GM.
C###  See-Also: GM,SPARSITY STRUCTURES

      FORMAT='($,'' Max# dimension of GM           (NZ_GM_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            IF(NZT(4,nx).GT.NMAX) NMAX=NZT(4,nx)
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NZ_GM_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NZ_GM_M=IDATA(1)

C#### Variable: NZ_GMM_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NZ_GMM_M is the maximum number of nonzero elements in GMM.
C###  See-Also: GMM,SPARSITY STRUCTURES

      FORMAT='($,'' Max# dimension of GMM         (NZ_GMM_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=0,NRT
              DO nrc=1,NRCM
                IF(NZZT(4,nr,nx).GT.NMAX) NMAX=NZZT(4,nr,nx)
              ENDDO !nrc
            ENDDO !nr
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NZ_GMM_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NZ_GMM_M=IDATA(1)

C#### Variable: NZ_GQ_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NZ_GQ_M is the maximum number of nonzero elements in GQ.
C###  See-Also: GQ,SPARSITY STRUCTURES

      FORMAT='($,'' Max# dimension of GQ           (NZ_GQ_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            IF(NZT(2,nx).GT.NMAX) NMAX=NZT(2,nx)
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NZ_GQ_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NZ_GQ_M=IDATA(1)

C#### Variable: ISC_GDM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISC_GDM is the maximum dimension of ISC_GD.
C###  See-Also: ISC_GD

      FORMAT='($,'' Max# dimension of ISC_GD      (NISC_GDM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(NZT(3,nx).GT.NMAX) NMAX=NZT(3,nx)
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISC_GDM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISC_GDM=IDATA(1)

C#### Variable: ISR_GDM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISR_GDM is the maximum dimension of ISR_GD.
C###  See-Also: ISR_GD

      FORMAT='($,'' Max# dimension of ISR_GD      (NISR_GDM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((NYT(1,3,nx)+1).GT.NMAX) NMAX=NYT(1,3,nx)+1
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISR_GDM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISR_GDM=IDATA(1)

C#### Variable: ISC_GKM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISC_GKM is the maximum dimension of ISC_GK.
C###  See-Also: ISC_GK

      FORMAT='($,'' Max# dimension of ISC_GK      (NISC_GKM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(NZT(1,nx).GT.NMAX) NMAX=NZT(1,nx)
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISC_GKM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISC_GKM=IDATA(1)

C#### Variable: ISR_GKM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISR_GKM is the maximum dimension of ISR_GK.
C###  See-Also: ISR_GK

      FORMAT='($,'' Max# dimension of ISR_GK      (NISR_GKM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((NYT(1,1,nx)+1).GT.NMAX) NMAX=NYT(1,1,nx)+1
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISR_GKM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISR_GKM=IDATA(1)

C#### Variable: ISC_GKKM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISC_GKKM is the maximum dimension of ISC_GKK.
C###  See-Also: ISC_GKK

      FORMAT='($,'' Max# dimension of ISC_GKK    (NISC_GKKM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((SOLVEROPTION(nx).EQ.1).AND.(SPARSEGKK(nx).GT.0)) THEN
                IF(NZZT(1,nr,nx).GT.NMAX) NMAX=NZZT(1,nr,nx)
              ELSEIF(((SOLVEROPTION(nx).EQ.3).OR.
     '          (SOLVEROPTION(nx).EQ.4)).AND.(SPARSEGKK(nx).GT.0)) THEN
                IF(NZZT(1,nr,nx).GT.NMAX) NMAX=NZZT(1,nr,nx)
              ENDIF
              IF(SOLVELU_NZA(nx).GT.NMAX) NMAX=SOLVELU_NZA(nx)
              IF(SPARSEGKK(nx).EQ.5) THEN
                IF(2*SOLVELU_NZA(nx).GT.NMAX) NMAX=2*SOLVELU_NZA(nx)
              ENDIF
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISC_GKKM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISC_GKKM=IDATA(1)

C#### Variable: ISR_GKKM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISR_GKKM is the maximum dimension of ISR_GKK.
C###  See-Also: ISR_GKK

      FORMAT='($,'' Max# dimension of ISR_GKK    (NISR_GKKM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((SOLVEROPTION(nx).EQ.1).AND.(SPARSEGKK(nx).GT.0).AND.
     '          (SPARSEGKK(nx).NE.5)) THEN
                IF(NZZT(1,nr,nx).GT.NMAX) NMAX=NZZT(1,nr,nx)
              ELSEIF(((SOLVEROPTION(nx).EQ.3).OR.
     '          (SOLVEROPTION(nx).EQ.4)).AND.(SPARSEGKK(nx).GT.0)) THEN
                IF((NYT(1,1,nx)+1).GT.NMAX) NMAX=NYT(1,1,nx)+1
                IF((NOT(1,1,nr,nx)+1).GT.NMAX) NMAX=NOT(1,1,nr,nx)+1
              ENDIF
              IF(SPARSEGKK(nx).NE.5) THEN
                IF(SOLVELU_NZA(nx).GT.NMAX) NMAX=SOLVELU_NZA(nx)
              ENDIF
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISR_GKKM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISR_GKKM=IDATA(1)

C#### Variable: ISC_GMM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISC_GMM is the maximum dimension of ISC_GM.
C###  See-Also: ISC_GM

      FORMAT='($,'' Max# dimension of ISC_GM      (NISC_GMM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(NZT(4,nx).GT.NMAX) NMAX=NZT(4,nx)
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISC_GMM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISC_GMM=IDATA(1)

C#### Variable: ISR_GMM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISR_GMM is the maximum dimension of ISR_GM.
C###  See-Also: ISR_GM

      FORMAT='($,'' Max# dimension of ISR_GM      (NISR_GMM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((NYT(1,4,nx)+1).GT.NMAX) NMAX=NYT(1,4,nx)+1
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISR_GMM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISR_GMM=IDATA(1)

C#### Variable: ISC_GMMM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISC_GMMM is the maximum dimension of ISC_GMM.
C###  See-Also: ISC_GMM

      FORMAT='($,'' Max# dimension of ISC_GMM    (NISC_GMMM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((SOLVEROPTION(nx).EQ.1).AND.(SPARSEGKK(nx).GT.0)) THEN
                IF(NZZT(4,nr,nx).GT.NMAX) NMAX=NZZT(4,nr,nx)
              ELSEIF(((SOLVEROPTION(nx).EQ.3).OR.
     '          (SOLVEROPTION(nx).EQ.4)).AND.(SPARSEGKK(nx).GT.0)) THEN
                IF(NZZT(4,nr,nx).GT.NMAX) NMAX=NZZT(4,nr,nx)
              ENDIF
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISC_GMMM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISC_GMMM=IDATA(1)

C#### Variable: ISR_GMMM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISR_GMMM is the maximum dimension of ISR_GMM.
C###  See-Also: ISR_GMM

      FORMAT='($,'' Max# dimension of ISR_GMM    (NISR_GMMM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=NISR_GMMM
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((SOLVEROPTION(nx).EQ.1).AND.(SPARSEGKK(nx).GT.0)) THEN
                IF(NZZT(4,nr,nx).GT.NMAX) NMAX=NZZT(4,nr,nx)
              ELSEIF(((SOLVEROPTION(nx).EQ.3).OR.
     '          (SOLVEROPTION(nx).EQ.4)).AND.(SPARSEGKK(nx).GT.0)) THEN
                IF((NYT(1,4,nx)+1).GT.NMAX) NMAX=NYT(1,4,nx)+1
                IF((NOT(1,4,nr,nx)+1).GT.NMAX) NMAX=NOT(1,4,nr,nx)+1
              ENDIF
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISR_GMMM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISR_GMMM=IDATA(1)

C#### Variable: ISC_GQM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISC_GQM is the maximum dimension of ISC_GQ.
C###  See-Also: ISC_GQ

      FORMAT='($,'' Max# dimension of ISC_GQ      (NISC_GQM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(NZT(2,nx).GT.NMAX) NMAX=NZT(2,nx)
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISC_GQM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISC_GQM=IDATA(1)

C#### Variable: ISR_GQM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    ISR_GQM is the maximum dimension of ISR_GQ.
C###  See-Also: ISR_GQ

      FORMAT='($,'' Max# dimension of ISR_GQ      (NISR_GQM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          IF(KTYP24.EQ.1) THEN
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((NYT(1,2,nx)+1).GT.NMAX) NMAX=NYT(1,2,nx)+1
            ENDDO
          ENDIF
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NISR_GQM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NISR_GQM=IDATA(1)

C#### Variable: NZ_MINOSM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NZ_MINOSM is the maximum size of the minos sparse arrays.

      FORMAT='($,'' Max# size of Minos arrays    (NZ_MINOSM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          nl=0
          IF(KTYP29.EQ.2) nl=1
          NMAX=NTOPTI*(NTCNTR+nl)
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NZ_MINOSM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NZ_MINOSM=IDATA(1)

C#### Variable: NBFM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NBFM is the maximum number of basis function families.  A basis
C###    function family consists of a parent basis function and its
C###    children.  A family is defined by the parent basis function; the
C###    children are variations upon this; e.g.,in a BEM application,
C###    the children may differ from the parent by the number of Gauss
C###    points used.

      FORMAT='($,'' Max# basis function families      (NBFM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NBFT.GT.0) THEN
            IDATA(1)=NBFT
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NBFM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NBFM=IDATA(1)

C#### Variable: NCOM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NCOM is the maximum number of non-linear constraints
C###    for optimisation problems.

      FORMAT='($,'' Max# nonlin. optim.n constraints  (NCOM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NTCNTR.GT.0) THEN
            IDATA(1)=NTCNTR
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NCOM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NCOM=IDATA(1)


C#### Variable: NDEM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NDEM is the maximum number of data points in one element.

      FORMAT='($,'' Max# data points in one element   (NDEM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(USE_DATA.EQ.0) THEN
            IDATA(1)=1
          ELSE
            NMAX=-IMAX
            DO ne=1,NET(0)
              IF(NDLT(ne).GT.NMAX) NMAX=NDLT(ne)
            ENDDO
            IF(.NOT.CALL_DATA) NMAX=0
            IF(NMAX.GT.0) THEN
              IDATA(1)=NMAX
            ELSE
              IDATA(1)=1
            ENDIF
          ENDIF
        ELSE
          IDATA(1)=NDEM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NDEM=IDATA(1)

C#### Variable: NDIPOLEM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NDIPOLEM is the maximum number of dipoles.

      FORMAT='($,'' Max# dipoles in a region      (NDIPOLEM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(USE_DIPOLE.EQ.0) THEN
            IDATA(1)=1
          ELSE
            NMAX=-IMAX
            DO nr=1,NRM
              DO nx=1,NXM
                IF(NDIPOLES(nr,nx).GT.NMAX) NMAX=NDIPOLES(nr,nx)
              ENDDO
            ENDDO
            IF(NMAX.LE.0) NMAX=1
            NDIPOLEM=NMAX
          ENDIF
        ELSE
          IDATA(1)=NDIPOLEM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NDIPOLEM=IDATA(1)

C#### Variable: NDIPTIMM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NDIPTIMM is the maximum number of dipole time points.

      FORMAT='($,'' Max# time points for a dipole (NDIPTIMM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(USE_DIPOLE.EQ.1) THEN
            NMAX=-IMAX
            DO nr=1,NRM
              DO nx=1,NXM
                DO i=1,NDIPOLES(nr,nx)
                  IF(DIPOLE_CEN_NTIME(i,nr,nx).GT.NMAX)
     '              NMAX=DIPOLE_CEN_NTIME(i,nr,nx)
                  IF(DIPOLE_DIR_NTIME(i,nr,nx).GT.NMAX)
     '              NMAX=DIPOLE_DIR_NTIME(i,nr,nx)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            IDATA(1)=1
          ENDIF
          IF(NMAX.LE.0) NMAX=1
          NDIPTIMM=NMAX
        ELSE
          IDATA(1)=NDIPTIMM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NDIPTIMM=IDATA(1)

C#### Variable: NELM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NELM is the maximum of elements adjacent to a line.

      FORMAT='($,'' Max# elements along a line        (NELM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nl=1,NLT
            IF(NEL(0,nl).GT.NMAX) NMAX=NEL(0,nl)
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NELM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NELM=IDATA(1)

C#### Variable: NEPM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NEPM is the maximum of elements a node can belong to.

      FORMAT='($,'' Max# elements a node can be in    (NEPM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nr=1,NRT
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(NENP(np,0,0).GT.NMAX) NMAX=NENP(np,0,0)
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NEPM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NEPM=IDATA(1)

C#### Variable: NGRSEGM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NGRSEGM is the maximum of graphics segments.

      FORMAT='($,'' Max# segments                  (NGRSEGM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
           NGRSEGM=2*NTLINE
          IF(NGRSEGM.GT.0) THEN
            IDATA(1)=NGRSEGM
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NGRSEGM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NGRSEGM=IDATA(1)

C#### Variable: NIQM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIQM is the maximum number of variables per grid point.

      FORMAT='($,'' Max# variables per grid point     (NIQM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NIQ_LIST(0).EQ.0) THEN
            IDATA(1)=1*USE_GRID
          ELSE
            IDATA(1)=NIQ_LIST(0)
          ENDIF
        ELSE
          IDATA(1)=NIQM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIQM=IDATA(1)

C#### Variable: NIQSM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIQSM is the maximum number of cell variables (state+derived).

      FORMAT='($,'' Max# cell state variables        (NIQSM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(CALL_CELL) THEN
            IDATA(1)=NIQST
          ELSE
            IDATA(1)=1*USE_CELL
          ENDIF
        ELSE
          IDATA(1)=NIQSM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIQSM=IDATA(1)

C#### Variable: NIFEXTM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIFEXTM is the maximum number fibre extension variables.
C###  See-Also: FEXT

      FORMAT='($,'' Max# variables for fibre extens(NIFEXTM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          FOUND1=.FALSE.
          FOUND2=.FALSE.
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(ITYP1(nr,nx).EQ.5) FOUND1=.TRUE.
            ENDDO
            IF(KTYP53(nr).EQ.3) FOUND2=.TRUE.
          ENDDO
          IF(FOUND1.AND.FOUND2) THEN
            NMAX=7
          ELSE
            NMAX=0
          ENDIF
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF((ITYP6(nr,nx).EQ.2).AND.(NMAX.LT.7)) NMAX=7
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NIFEXTM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIFEXTM=IDATA(1)

C#### Variable: NIYM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIYM is the maximum number of variables per mesh dof.

      FORMAT='($,'' Max# variables per mesh dof       (NIYM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=0
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO ny=1,NYM
              TOTAL=0
              DO nw=1,NIYM
                IF(YP(ny,nw,nx).GT.0.0d0) TOTAL=nw
              ENDDO
              IF(TOTAL.GT.NMAX) NMAX=TOTAL
            ENDDO
          ENDDO
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              IF((ITYP1(nr,nx).EQ.9).AND.(NMAX.LT.4)) NMAX=4
              IF((ITYP6(nr,nx).EQ.2).AND.(NMAX.LT.7)) NMAX=7
            ENDDO
          ENDDO
          IF(CALL_ANAL.AND.(NMAX.LT.7)) NMAX=7
          IF((KTYP22.EQ.1).AND.(NMAX.LT.10)) NMAX=10
          IF((KTYP22.EQ.2).AND.(NMAX.LT.13)) NMAX=13
          IF((KTYP22.EQ.3).AND.(NMAX.LT.16)) NMAX=16
          IF(COUPLED_QUASI.AND.(NMAX.LT.16)) NMAX=16
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NIYM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIYM=IDATA(1)

C#### Variable: NIYFIXM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIYFIXM is the maximum number of boundary condition status
C###    variables per mesh dof.
C###  See-Also: FIX

      FORMAT='($,'' Max# variables / mesh dof(fix) (NIYFIXM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NIYFIXM.GT.0) THEN
            IDATA(1)=NIYFIXM
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NIYFIXM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIYFIXM=IDATA(1)

C#### Variable: NIYGM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIYGM is the maximum number of variables in the element Gauss
C###    point array YG.
C###  See-Also: YG

      FORMAT='($,'' Max# vars. at each gauss point   (NIYGM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=0
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3) THEN
                NMAX2=-IMAX
                DO nb=1,NBT
                  IF(NIT(nb).GT.NMAX2) NMAX2=NIT(nb)
                ENDDO !nb
                IF(NMAX2.GT.0) NMAX2=NMAX2*(NMAX2+1)/22+NMAX
                IF(NMAX2.GT.NMAX) NMAX=NMAX2
              ENDIF
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NIYGM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIYGM=IDATA(1)

C#### Variable: NIYGFM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIYGFM is the maximum number of variables in the face Gauss
C###    point array YGF.
C###  See-Also: YGF

      FORMAT='($,'' Max# vars. at face gauss points (NIYGFM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=0
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3) THEN
                NMAX2=-IMAX
                DO nb=1,NBT
                  IF(NIT(nb).GT.NMAX2) NMAX2=NIT(nb)
                ENDDO !nb
                IF(NMAX2.GT.0) NMAX2=NMAX2*(NMAX2-1)/2
                IF(NMAX2.GT.NMAX) NMAX=NMAX2
              ENDIF
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NIYGFM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIYGFM=IDATA(1)

C#### Variable:NLCM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NLCM is the maximum number of linear constraints
C###    for optimisation problems.

      FORMAT='($,'' Max# linear optimis.n constraints (NLCM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NLCM=0
          IF(KTYP29.EQ.2) NLCM=1
          IF(NLCM.GT.0) THEN
            IDATA(1)=NLCM
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NLCM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NLCM=IDATA(1)

C#### Variable: NMAQM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NMAQM is the maximum number of auxiliary grid parameters.

      FORMAT='($,'' Max# auxiliary grid parameters   (NMAQM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(MAQ_LIST(0).EQ.0) THEN
            IDATA(1)=1*USE_GRID
          ELSE
            IDATA(1)=MAQ_LIST(0)
          ENDIF
        ELSE
          IDATA(1)=NMAQM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NMAQM=IDATA(1)

C#### Variable: NMQM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NMQM is the maximum number of cell material parameters.

      FORMAT='($,'' Max# cell material parameters     (NMQM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(CALL_CELL) THEN
            IDATA(1)=NMQT
          ELSE
            IDATA(1)=1*USE_CELL
          ENDIF
        ELSE
          IDATA(1)=NMQM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NMQM=IDATA(1)

C#### Variable: NOPM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NOPM is the maximum number of optimisation variables.

      FORMAT='($,'' Max# optimisation variables       (NOPM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NTOPTI.GT.0) THEN
            IDATA(1)=NTOPTI
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NOPM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NOPM=IDATA(1)

C#### Variable: NORM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NORM is the maximum size of fractal tree order arrays.

      FORMAT='($,'' Max size fractal tree order array (NORM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NORM=0
          IF(NT_GEN.GT.0) NORM=NE_R_M
          IF(NORM.GT.0) THEN
            IDATA(1)=NORM
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NORM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NORM=IDATA(1)

C#### Variable: NOYM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NOYM is the maximum number of solution dofs attached to one
C###    mesh dof.

      FORMAT='($,'' Max# soln dofs for mesh dof       (NOYM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              IF(ITYP4(nr,nx).NE.4.AND.ITYP4(nr,nx).NE.6) THEN
                !not collocation or gird-base FE
                DO nrc=1,NRCM
                  DO ny=1,NYT(nrc,1,nx) ! nc=1 ?
                    IF(NONY(0,ny,nrc,nr,nx).GT.NMAX) THEN
                      NMAX=NONY(0,ny,nrc,nr,nx)
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF !ITYP4
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NOYM
         ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NOYM=IDATA(1)

C#### Variable: NPDM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NPDM is the maximum number of domains for Boundary Element
C###    problems.

      FORMAT='($,'' Max# domain nodes for BE problems (NPDM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IDATA(1)=1
        ELSE
          IDATA(1)=NPDM
         ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NPDM=IDATA(1)

C#### Variable: NQEM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQEM is the maximum number of grid points in an element

      FORMAT='($,'' Max# grid points per element      (NQEM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nq=1,NQSCM
            IF(NQET(nq).GT.NMAX) NMAX=NQET(nq)
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NQEM
         ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQEM=IDATA(1)

C#### Variable: NQGM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQGM is the maximum number of non-zeros in grid matrix rows

      FORMAT='($,'' Max# non-zeros in grid matrix row (NQGM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NITMAX=0
          DO nqsc=1,NQSCT
            IF(NQXI(0,nqsc).GT.NITMAX) NITMAX=NQXI(0,nqsc)
          ENDDO
          NMAX=0
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              IF(ITYP4(nr,nx).EQ.4) THEN !collocation or grid-based FE
                IF(NITMAX.EQ.3) THEN
                  N=22
                ELSE
                  N=3**NITMAX
                ENDIF
              ELSEIF(ITYP4(nr,nx).EQ.6.OR.ITYP4(nr,nx).EQ.7) THEN !grid-based FE
                N=3**NITMAX ! 3x3x3 support for a given grid point in GBFE
C MLT 29Nov02 added grid finite volume method
              ELSEIF(ITYP4(nr,nx).EQ.7) THEN !grid-based FV
                N=2*NITMAX+1 ! 3(1D), 5(2D), 7(3D)
              ELSE
                N=0
              ENDIF !ITYP4
              IF(N.GT.NMAX) NMAX=N
            ENDDO !nr
          ENDDO !nx
          IDATA(1)=NMAX
        ELSE
          IDATA(1)=NQGM
        ENDIF !MINIMAL
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQGM=IDATA(1)

C#### Variable: NQIM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQIM is the maximum number of cell integer parameters.

      FORMAT='($,'' Max# cell integer variables       (NQIM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(CALL_CELL) THEN
            IDATA(1)=NQIT
          ELSE
            IDATA(1)=1*USE_CELL
          ENDIF
        ELSE
          IDATA(1)=NQIM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQIM=IDATA(1)


C#### Variable: NQRM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQRM is the maximum number of cell real parameters.

      FORMAT='($,'' Max# cell real variables          (NQRM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(CALL_CELL) THEN
            IDATA(1)=NQRT
          ELSE
            IDATA(1)=1*USE_CELL
          ENDIF
        ELSE
          IDATA(1)=NQRM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQRM=IDATA(1)

C#### Variable: NQISVM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQISVM is the maximum number of spatially varying cell
C###    integer parameters.

      FORMAT='($,'' Max# spatial var cell int vars  (NQISVM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(CALL_CELL) THEN
            IDATA(1)=NQISVT
          ELSE
            IDATA(1)=1*USE_CELL
          ENDIF
        ELSE
          IDATA(1)=NQISVM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQISVM=IDATA(1)

C#### Variable: NQRSVM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQRSVM is the maximum number of spatially varying cell
C###    real parameters.

      FORMAT='($,'' Max# spatial var cell real vars (NQRSVM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(CALL_CELL) THEN
            IDATA(1)=NQRSVT
          ELSE
            IDATA(1)=1*USE_CELL
          ENDIF
        ELSE
          IDATA(1)=NQRSVM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQRSVM=IDATA(1)

C#### Variable: NQSCM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQSCM is the maximum number of grid schemes.

      FORMAT='($,'' Max# number of grid schemes      (NQSCM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=NQSCT
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NQSCM
         ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQSCM=IDATA(1)

C#### Variable: NQVM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NQVM is the maximum number of cell variants within a
C###    particular cell model.

      FORMAT='($,'' Max# cell variants                (NQVM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(CALL_CELL) THEN
            IDATA(1)=NQVT
          ELSE
            IDATA(1)=1*USE_CELL
          ENDIF
        ELSE
          IDATA(1)=NQVM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NQVM=IDATA(1)

C#### Variable: NRCM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NRCM is the maximum number of rows and columns.  This should
C###    always be 2.

      FORMAT='($,'' Max# rows and columns (sb 2)      (NRCM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=NRCM
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NRCM=IDATA(1)

C#### Variable: NREM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NREM is the maximum number of residuals for optimisation
C###    problems.

      FORMAT='($,'' Max# optimisation residuals       (NREM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(NT_RES.GT.0) THEN
            IDATA(1)=NT_RES
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NREM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NREM=IDATA(1)

C#### Variable: NTIMEPOINTSM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NTIMEPOINTSM is the maximum number of time points which
C###    are allowed to be set for a time variable.

C#### Variable: NTIMEPOINTST
C###  Type: INTEGER
C###  Set_up: FEMINI,IPTIME
C###  Description:
C###    NTIMEPOINTST is the maximum number of time points that
C###    have been set through the define time command.

      FORMAT='($,'' Max# time points          (NTIMEPOINTSM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IDATA(1)=NTIMEPOINTST
          IF(NTIMEPOINTST.LE.0) IDATA(1)=1
        ELSE
          IDATA(1)=NTIMEPOINTSM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NTIMEPOINTSM=IDATA(1)

C#### Variable: NTIMEVARSM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NTIMEVARSM is the maximum number of time variables which
C###    can be set up.

C#### Variable: NTIMEVARST
C###  Type: INTEGER
C###  Set_up: FEMINI,IPTIME
C###  Description:
C###    NTIMEVARSM is the maximum number of time variables that
C###    have been set through the define time command.

      FORMAT='($,'' Max# time variables         (NTIMEVARSM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IDATA(1)=NTIMEVARST
          IF(NTIMEVARST.LE.0) IDATA(1)=1
        ELSE
          IDATA(1)=NTIMEVARSM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NTIMEVARSM=IDATA(1)

C#### Variable: NYOM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NYOM is the maximum number of mesh dofs attached to one
C###    solution dof.

      FORMAT='($,'' Max# mesh dofs for soln dof       (NYOM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO no=1,NOM
            DO nr=0,NRT
              DO nrc=1,NRCM
                DO nxx=1,NX_LIST(0)
                  nx=NX_LIST(nxx)
                  IF(NYNO(0,no,nrc,nr,nx).GT.NMAX)
     '              NMAX=NYNO(0,no,nrc,nr,nx)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NYOM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,500,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NYOM=IDATA(1)

C#### Variable: NYROWM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NYROWM is the maximum number of rows in the unreduced right hand
C###    side vector GR.

      FORMAT='($,'' Max# rows in a problem          (NYROWM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              DO nc=1,NCT(nr,nx)
                IF(NYT(1,nc,nx).GT.NMAX) NMAX=NYT(1,nc,nx)
              ENDDO !nc
            ENDDO !nr
          ENDDO !nx
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NYROWM
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NYROWM=IDATA(1)

C#### Variable: NIMAGEM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NIMAGEM is the maximum image cell array dimension.

      FORMAT='($,'' Max image cell array dimension (NIMAGEM)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        NIMAGEM=1
        IDATA(1)=NIMAGEM
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIMAGEM=IDATA(1)

C#### Variable: NY_TRANSFER_M
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NY_TRANSFER_M is the maximum size of the transfer matrix.

      FORMAT='($,'' Size of transfer matrix  (NY_TRANSFER_M)[1]: '',I9)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          NMAX=-IMAX
          DO nrc=1,2
            IF(ISIZE_TBH(nrc).GT.NMAX) NMAX=ISIZE_TBH(nrc)
          ENDDO
          IF(NMAX.GT.0) THEN
            IDATA(1)=NMAX
          ELSE
            IDATA(1)=1
          ENDIF
        ELSE
          IDATA(1)=NY_TRANSFER_M
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NY_TRANSFER_M=IDATA(1)

C#### Variable: NYYM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NYYM is the maximum number of mesh dofs mapped to one
C###    mesh dof.

      FORMAT='($,'' Max# mesh dofs map to 1 mesh dof  (NYYM)[1]: '',I9)'
      IDEFLT(1)=0
      IF(IOTYPE.EQ.3) THEN
        IDATA(1)=NYYM
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NYYM=IDATA(1)

C#### Variable: USE_BEM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_BEM determines whether to allocate memory required for
C###    boundary elements (1) or not (0).

      FORMAT='($,'' USE_BEM       (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_BEM=0
          DO nb=1,NBT
            IF((NBC(nb).EQ.5).OR.(NBC(nb).EQ.6)) THEN
              USE_BEM=1
            ENDIF
          ENDDO
        ENDIF
        IDATA(1)=USE_BEM
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_BEM=IDATA(1)

C#### Variable: USE_CELL
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_CELL determines whether to allocate memory required for cell
C###    models (1) or not (0).

      FORMAT='($,'' USE_CELL      (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_CELL=0
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            DO nr=1,NRT
              IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) USE_CELL=1
            ENDDO
          ENDDO
        ENDIF
        IDATA(1)=USE_CELL
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_CELL=IDATA(1)

C#### Variable: USE_DATA
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_DATA determines whether to allocate memory required for data
C###    points (1) or not (0).

      FORMAT='($,'' USE_DATA      (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_DATA=0
          IF(CALL_DATA) USE_DATA=1
        ENDIF
        IDATA(1)=USE_DATA
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_DATA=IDATA(1)

C#### Variable: USE_DIPOLE
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_DIPOLE determines whether to allocate memory required for
C###    dipoles (1) or not (0).

      FORMAT='($,'' USE_DIPOLE    (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(USE_DIPOLE.EQ.1) THEN
            USE_DIPOLE=0
            DO nr=1,NRM
              DO nx=1,NXM
                IF(NDIPOLES(nr,nx).GT.0) USE_DIPOLE=1
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        IDATA(1)=USE_DIPOLE
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_DIPOLE=IDATA(1)

C#### Variable: USE_GAUSS_PT_MATERIALS
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_GAUSS_PT_MATERIALS determines whether to allocate memory
C###    required for definition of any material propterties by Gauss
C###    points (1) or not (0).

      FORMAT='($,'' USE_GAUSS_PT_MATERIALS  (0 or 1)[0]: '',I1)'
      IDEFLT(1)=0
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_GAUSS_PT_MATERIALS=0
        ENDIF
        IDATA(1)=USE_GAUSS_PT_MATERIALS
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_GAUSS_PT_MATERIALS=IDATA(1)

C#### Variable: USE_GRAPHICS
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_GRAPHICS determines whether to allocate memory required for
C###    graphics commands prefixed by 'FEM' (1) or not (0).

      FORMAT='($,'' USE_GRAPHICS  (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_GRAPHICS=0
          DO nw=1,NWM
            IF(IWKS(nw).GT.0) USE_GRAPHICS=1
          ENDDO
        ENDIF
        IDATA(1)=USE_GRAPHICS
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_GRAPHICS=IDATA(1)

C#### Variable: USE_GRID
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_GRID determines whether to allocate memory required for
C###    collocation grids (1) or not (0).

      FORMAT='($,'' USE_GRID      (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_GRID=0
          IF(CALL_GRID) USE_GRID=1
        ENDIF
        IDATA(1)=USE_GRID
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_GRID=IDATA(1)

C#### Variable: USE_LUNG
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_LUNG determines whether to allocate memory required for lung
C###    problems (1) or not (0).

      FORMAT='($,'' USE_LUNG      (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_LUNG=0
          IF(NT_GEN.GT.0) USE_LUNG=1
        ENDIF
        IDATA(1)=USE_LUNG
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_LUNG=IDATA(1)

C#### Variable: USE_MAGNETIC
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_MAGNETIC determines whether to allocate memory
C###    required for magnetic problems (1) or not (0).

      FORMAT='($,'' USE_MAGNETIC  (0 or 1)[0]: '',I1)'
      IDEFLT(1)=0
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF((USE_MAGNETIC.EQ.1).AND.(ISIZE_MFI(1,1).GT.0)) THEN
            USE_MAGNETIC=1
          ELSE
            USE_MAGNETIC=0
          ENDIF
        ENDIF
        IDATA(1)=USE_MAGNETIC
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_MAGNETIC=IDATA(1)



C#### Variable: USE_MAPS
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_MAPS determines whether to allocate memory required for the
C###    special version mapping (1) or not (0).

      FORMAT='($,'' USE_MAPS      (0 or 1)[0]: '',I1)'
      IDEFLT(1)=0
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          IF(JTYP2A.EQ.1.OR.JTYP2B.EQ.1.OR.JTYP2C.EQ.1) THEN
            USE_MAPS=1
          ELSE
            USE_MAPS=0
          ENDIF
        ENDIF
        IDATA(1)=USE_MAPS
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_MAPS=IDATA(1)

C#### Variable: USE_MINOS
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_MINOS determines whether to allocate memory required for the
C###    nonlinear solver package MINOS (1) or not (0).

      FORMAT='($,'' USE_MINOS     (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_MINOS=0
          IF(KTYP29.EQ.2) USE_MINOS=1
        ENDIF
        IDATA(1)=USE_MINOS
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_MINOS=IDATA(1)

C#### Variable: USE_NLSTRIPE
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_NLSTRIPE determines whether to allocate memory required for
C###    non-linear strip fitting (1) or not (0).

      FORMAT='($,'' USE_NLSTRIPE  (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) USE_NLSTRIPE=0
        IDATA(1)=USE_NLSTRIPE
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_NLSTRIPE=IDATA(1)

C#### Variable: USE_NONLIN
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_NONLIN determines whether to allocate memory required for
C###    non-linear problems involving elasticity or pressure boundary
C###    conditions (1) or not (0).

      FORMAT='($,'' USE_NONLIN    (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_NONLIN=0
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(ITYP6(nr,nx).EQ.2) USE_NONLIN=1
            ENDDO
          ENDDO
        ENDIF
        IDATA(1)=USE_NONLIN
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_NONLIN=IDATA(1)

C#### Variable: USE_NPSOL
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_NPSOL determines whether to allocate memory required for the
C###    non-linear solver package NPSOL (1) or not (0).

      FORMAT='($,'' USE_NPSOL     (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_NPSOL=0
          IF(KTYP29.EQ.1) USE_NPSOL=1
        ENDIF
        IDATA(1)=USE_NPSOL
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_NPSOL=IDATA(1)

C#### Variable: USE_SPARSE
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_SPARSE determines whether to allocate memory required for
C###    sparse matrix structures (1) or not (0).

      FORMAT='($,'' USE_SPARSE    (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_SPARSE=0
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            IF(SPARSEGKK(nx).GT.0) USE_SPARSE=1
          ENDDO
          IF(KTYP24.EQ.1) USE_SPARSE=1
        ENDIF
        IDATA(1)=USE_SPARSE
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_SPARSE=IDATA(1)

C#### Variable: USE_TRANSFER
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_TRANSFER determines whether to allocate memory required for
C###    a transfer coefficient matrix (1) or not (0).

      FORMAT='($,'' USE_TRANSFER  (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_TRANSFER=0
          IF(EVALUATE_TRANSFER) USE_TRANSFER=1
        ENDIF
        IDATA(1)=USE_TRANSFER
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_TRANSFER=IDATA(1)

C#### Variable: USE_TRIANGLE
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_TRIANGLE does not exist.

      FORMAT='($,'' USE_TRIANGLE  (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IDATA(1)=0
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

C#### Variable: USE_VORONOI
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_VORONOI determines whether to allocate and set up memory
C###    required for Voronoi triangulations (currently only for fluid
C###    flow) (1) or not (0).  If set to 1, element-node connectivity
C###    array NENP, lines, and faces will not be calculated.

      FORMAT='($,'' USE_VORONOI   (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_VORONOI=0
        ENDIF
        IDATA(1)=USE_VORONOI
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,7,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_VORONOI=IDATA(1)

C AJP 13-11-97
C#### Variable: USE_TIME
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    USE_TIME determines whether to allocate memory required for
C###    time-dependent arrays (transfer matrix stuff only at the moment)
C###    (1) or not (0).

      FORMAT='($,'' USE_TIME      (0 or 1)[1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) THEN
        IF(MINIMAL) THEN
          USE_TIME=USE_TRANSFER
        ENDIF
        IDATA(1)=USE_TIME
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,1,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) USE_TIME=IDATA(1)

      NNEPM=MAX(NEM,NPM)
C      NEQM=MAX(NEM,NQM)
      NLISTM=MAX(NEM,NPM,NLM,NDM,NQM)

C#### Variable: NOOPM
C###  Type: INTEGER
C###  Set_up: IPPARA,SYNTAX
C###  Description:
C###    NOOPM is the maximum number of optimisation or solution
C###    variables.
C###    It is set as NOOPM=MAX(NOM,NOPM)

      NOOPM=MAX(NOM,NOPM)
      NFVM=(NFVCM/2)*NVCM

      USE_OPTI=MAX(USE_NPSOL,USE_MINOS)

C     Making sure that memory is freed for NPE_PTR
      IF(NPE_PTR.NE.0) THEN
         CALL FREE_MEMORY(NPE_PTR,ERROR,*9999)
      ENDIF

C removing warning as arrays are now reallocated by default in deparra
C NPS 15/6/99

      CALL EXITS('IPPARA')
      RETURN
 9999 CALL ERRORS('IPPARA',ERROR)
      CALL EXITS('IPPARA')
      RETURN 1
      END



c cpb 26/10/94 Called from fe28! should be changed.
