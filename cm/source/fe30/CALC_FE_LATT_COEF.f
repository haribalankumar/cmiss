      SUBROUTINE CALC_FE_LATT_COEF(nb,NITB,NLATNE,NQGP,NQGP_PIVOT,
     &  NQNLAT,NQS,NQXI,nr,nx_ext,COEFFSEXT,CQ,GM,NQGW,PG,WG,XQ,
     &  BIDOMAIN,FIXQ,SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: CALC_FE_LATT_COEF
C###  Description:
C###    CALC_FE_LATT_COEF is used for the Grid-based Finite Element
C###    scheme for lattice-based elements. It finds the coeffs of
C###    the nonzeros in stiffness matrix (NQGW) and adds coefficients
C###    to global mass matrix (GM) for all elements. The number of
C###    nonzeros depends on the number of neighbouring elements.
C***  Created by Greg Sands 28 August 2003

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc' !NJL_FIBR
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER nb,NITB,NLATNE(NEQM+1),NQGP(0:NQGM,NQM),
     &  NQGP_PIVOT(NQGM,NQM),NQNLAT(NEQM*NQEM),NQS(NEQM),
     &  NQXI(0:NIM,NQSCM),nr,nx_ext
      REAL*8 COEFFSEXT(NQGM,NQM),CQ(NMM,NQM),GM(NZ_GM_M),NQGW(NQGM,NQM),
     &  PG(NSM,NUM,NGM,NBM),WG(NGM,NBM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,FIXQ(NYQM,NIYFIXM,NXM),SOLVEEIGHTPROBLEM

!     Local variables
      INTEGER GB_ELEM(8),i,m,ne,NGTB,ni,nj,nk,NNTB,nq,nqsc,offset,
     &  NSTELEM,GELEMVMIN,GELEMVMAX,IVMIN(3),IVMAX(3),GGPMIN,GGPMAX,
     &  NPOELEM,NNEELEM,NZEELEM
      REAL*8 SFE_VOLUME,ELEM_VOL,NEGLIGIBLE_VOLUME,MIN_VOL,MAX_VOL,
     &  MEAN_VOL
      LOGICAL INT_ELEM,EXT_ELEM
      DATA NEGLIGIBLE_VOLUME /1.0D-7/

      CALL ENTERS('CALC_FE_LATT_COEF',*9999)

      NNTB = NNT(nb)
      NGTB = NGT(nb)

      IF(SOLVEEIGHTPROBLEM) THEN
        EXT_ELEM=.TRUE.
        INT_ELEM=.FALSE.
      ELSE
        IF(BIDOMAIN) THEN
          EXT_ELEM=.TRUE.
        ELSE
          EXT_ELEM=.FALSE.
        ENDIF
        INT_ELEM=.TRUE.
      ENDIF

C!!! LKC 12-AUG-2005 Note that there are effectively 3 nq loops in this routine,
C      maybe this can be streamlined somehow ...

C initialise NQGP, NQGW, GM
      DO nq=NQR(1,nr),NQR(2,nr)
        NQGP(0,nq)=0
        DO i=1,NQGM ! initialise variables
          NQGP(i,nq)=0
          COEFFSEXT(i,nq)=0.0d0
          NQGW(i,nq)=0.0d0
          GM(i+NQGM*(nq-1))=0.0d0
        ENDDO !i
      ENDDO !nq

C check that there are no nodal derivatives stored at nodes as only
C     linear basis functions can be interpolated at the moment
      CALL ASSERT(NKT(0,nb).EQ.1.AND.NNTB.EQ.2**NITB,
     &  '>>Must use a linear basis to interpolate'
     &  //' dependent variable for grid-based FE',ERROR,*9999)

C check that the number of Xi directions for the element is equal to
C     the number of global coordinates
      CALL ASSERT(NITB.EQ.NJ_LOC(NJL_GEOM,0,nr),
     &  '>>#Xi directions for element must equal #global'
     &  //' coordinates (for grid-based FE method)',ERROR,*9999)

C MLT2Sep05 Initialise structured element mesh statistics
      MIN_VOL = 1.0D10
      MAX_VOL = 0.0D0
      MEAN_VOL = 0.0D0
      NSTELEM = 0
C Counters for positive, negative and zero element volume respectively
      NPOELEM = 0
      NNEELEM = 0
      NZEELEM = 0

C  MLT25Aug05 possible parallelisation paths include:
C   1. construct structured element contributions internal to 
C      host element in parallel per host element, then do a 
C      serial loop and just add in the boundary layer of 
C      structured elements per host element.
C   2. construct GM, NQGW and COEFFSEXT host element-wise
C      in parallel and then simply add together (adding appropriate
C      overlap) when GKK is assembled - this would be quite different 
C      to non-lattice GFEM or GFVM, although perhaps those two other
C      methods could also be done this way? 
      DO ne=1,NEQM

        nqsc=NQS(ne)
C***    Loop over subelements (if any)
        DO nk=0,MAX(NQXI(3,nqsc)-2,0)
          DO nj=0,MAX(NQXI(2,nqsc)-2,0)
            DO ni=0,MAX(NQXI(1,nqsc)-2,0)
C***          Offset for first node in subelement
              offset=ni+nj*NQXI(1,nqsc)+nk*NQXI(1,nqsc)*NQXI(2,nqsc)
C***          Construct local subelement
              DO m=0,NNTB-1
                GB_ELEM(m+1)=NQNLAT(NLATNE(ne)+offset+MOD(m,2)+
     &            NQXI(1,nqsc)*MOD(INT(m/2),2)+
     &            NQXI(1,nqsc)*NQXI(2,nqsc)*MOD(INT(m/4),2))
              ENDDO

C MLT25Aug05  Compute volume of GB_ELEM - if it is non-zero then compute
C             the contribution to mass and stiffness matrices.
C             The volume is computed numerically. Perhaps a more efficient
C             way to do this would be to use the formula derived in
C             SIAM J. Sci. Comput. 23(4), pp. 1274-1290, 2001.
              ELEM_VOL = SFE_VOLUME(GB_ELEM,NGTB,NITB,NNTB,
     &                              nr,PG(1,1,1,nb),XQ)

C JCZ24OCT07  Compute the numbers of positive element volume, zero
C             element volume, and negative volume. 

              IF (ELEM_VOL .gt. NEGLIGIBLE_VOLUME) THEN
                NPOELEM = NPOELEM + 1
              ELSE IF (DABS(ELEM_VOL) .le. NEGLIGIBLE_VOLUME) THEN
                NZEELEM = NZEELEM + 1
              ELSE
                NNEELEM = NNEELEM + 1 
              END IF

C MLT25Aug05  Check that element volume is positive and significant
              IF (ELEM_VOL.GT.NEGLIGIBLE_VOLUME) THEN

C MLT25Aug05  Compute some structured element statistics
                NSTELEM = NSTELEM + 1
                IF (ELEM_VOL.LT.MIN_VOL) THEN
                  MIN_VOL = ELEM_VOL
                  GELEMVMIN = ne
                  IVMIN(1) = ni+1
                  IVMIN(2) = nj+1
                  IVMIN(3) = nk+1
                  GGPMIN = NQNLAT(NLATNE(ne)+offset+MOD(0,2)+
     &              NQXI(1,nqsc)*MOD(0,2)+
     &              NQXI(1,nqsc)*NQXI(2,nqsc)*MOD(0,2))
                ENDIF
                IF (ELEM_VOL.GT.MAX_VOL) THEN
                  MAX_VOL = ELEM_VOL
                  GELEMVMAX = ne
                  IVMAX(1) = ni+1
                  IVMAX(2) = nj+1
                  IVMAX(3) = nk+1
                  GGPMAX = NQNLAT(NLATNE(ne)+offset+MOD(0,2)+
     &              NQXI(1,nqsc)*MOD(0,2)+
     &              NQXI(1,nqsc)*NQXI(2,nqsc)*MOD(0,2))
                ENDIF
                MEAN_VOL = MEAN_VOL + ELEM_VOL

C Compute structured element contributions for each node and add contributions
C to mass matrix and striffness matrices
                DO m=1,NNTB
                  nq=GB_ELEM(m)
                  CALL CALC_FE_COEF(GB_ELEM,m,NGTB,NITB,NNTB,nq,
     &              NQGP(0,nq),nr,nx_ext,COEFFSEXT(1,nq),CQ,GM,
     &              NQGW(1,nq),PG(1,1,1,nb),WG(1,nb),XQ,BIDOMAIN,
     &              EXT_ELEM,FIXQ,INT_ELEM,SOLVEEIGHTPROBLEM,
     &              ERROR,*9999)
                ENDDO !nq/m
              ENDIF ! non-zero volume
            ENDDO !ni
          ENDDO !nj
        ENDDO !nk
      ENDDO !ne

      CALL ASSERT(MEAN_VOL.GT.NEGLIGIBLE_VOLUME,
     & 'No volume in FE mesh', ERROR, *9999)

C MLT2Sep05 Report structured element statistics
      WRITE(*,*) ' #### Structured tics #### '
      WRITE(*,*) '      Total number of structured elements: ',NSTELEM
      WRITE(*,*) '      Total mesh volume (mm^3): ',MEAN_VOL
      MEAN_VOL = MEAN_VOL/DBLE(NSTELEM)
      WRITE(*,*) '      The number of postive element volume (mm^3): ',
     &             NPOELEM 
      WRITE(*,*) '      The number of negative element volume (mm^3): ',
     &             NNEELEM
      WRITE(*,*) '      The number of zero element volume (mm^3): ',
     &             NZEELEM
      WRITE(*,*) '      Mean structured element volume (mm^3): ',
     &             MEAN_VOL
      WRITE(*,*) '      Minimum non-negative/non-zero structured ',
     &             'element volume (mm^3): ',MIN_VOL
      WRITE(*,*) '        Located in global element: ',GELEMVMIN
      WRITE(*,*) '        Local lattice indices: ',IVMIN(1),IVMIN(2),
     &             IVMIN(3)
      WRITE(*,*) '        Global grid point number of origin: ',GGPMIN
      WRITE(*,*) '      Maximum element volume (mm^3): ',MAX_VOL
      WRITE(*,*) '        Located in global element: ',GELEMVMAX
      WRITE(*,*) '        Local lattice indices: ',IVMAX(1),IVMAX(2),
     &             IVMAX(3)
      WRITE(*,*) '        Global grid point number of origin: ',GGPMAX
     

C essential bc for solve8 problem or no surrounding elements (intracellular)
      DO nq=NQR(1,nr),NQR(2,nr)
C!!! KAT 21Feb01: Need to set diagonal of GM if there are no surrounding
C!!! intra elements even if there are surrounding extra elements.
        IF (NQGP(0,nq).LE.1) THEN
          NQGP(0,nq) = 1
          NQGP(1,nq) = nq
          IF(SOLVEEIGHTPROBLEM) THEN !essential bc
            COEFFSEXT(1,nq)= 1.0d0 !fix grid point to rhs value
          ELSE !no surrounding elements
            NQGW(1,nq)= 0.0d0
            COEFFSEXT(1,nq)= 0.0d0 !??? IS THIS RIGHT ???
            GM(1+NQGM*(nq-1)) = 1.0d0
          ENDIF
        ENDIF
C KAT 21Feb01: Should move outside loop after check no surrounding
C       elements set COEFFSEXT diagonal to 1 for analytic or
C       essential bc
        IF (BIDOMAIN) THEN
          IF (FIXQ(nq,1,nx_ext).OR.FIXQ(nq,3,nx_ext))THEN
            DO i=1,NQGP(0,nq)
              IF (NQGP(i,nq).EQ.nq) COEFFSEXT(i,nq) = 1.0d0
            ENDDO !i
          ENDIF !analytic or essential bc
        ENDIF !bidomain
        
C sort NQGP
        CALL ISORTP(NQGP(0,nq),NQGP(1,nq),NQGP_PIVOT(1,nq))
         
      ENDDO !nq
      
      write(*,*) ' exiting lattice '

      CALL EXITS('CALC_FE_LATT_COEF')
      RETURN
 9999 CALL ERRORS('CALC_FE_LATT_COEF',ERROR)
      CALL EXITS('CALC_FE_LATT_COEF')
      RETURN 1
      END

      REAL*8 FUNCTION SFE_VOLUME(E,NGTB,NITB,NNTB,nr,PG,XQ)

C#### Function: SFE_VOLUME
C###  Type: REAL*8
C###  Description:
C###    SFE_VOLUME determines volume of structured finite elements.
C***  Created by Mark Trew 25 August 2005

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc' !NJL_FIBR

!     Parameter List
      INTEGER E(8),NGTB,NITB,NNTB,nr
      REAL*8 PG(NSM,NUM,NGM),XQ(NJM,NQM)

!     Local variables
      INTEGER ng,nj1,nj,ni,nn,nu(3)
      REAL*8 DET,dxdXi(3,3)
      DATA nu /2,4,7/

C initialise element volume
      SFE_VOLUME = 0.0d0

      DO ng=1,NGTB !gauss point loop
C create dXidx matrix and store Jacobian
        DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,nj1,nr)
          DO ni=1,NITB
            dxdXi(nj,ni) = 0.0d0
            DO nn=1,NNTB
              dxdXi(nj,ni) = dxdXi(nj,ni) + XQ(nj,E(nn)) *
     &          PG(nn,nu(ni),ng)
            ENDDO !nn
          ENDDO !ni
        ENDDO !nj1
        SFE_VOLUME = SFE_VOLUME + DET(dxdXi) 
      ENDDO
      RETURN
      END
