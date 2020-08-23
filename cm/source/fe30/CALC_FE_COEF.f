      SUBROUTINE CALC_FE_COEF(GB_ELEM,m,NGTB,NITB,NNTB,nq,NQGP,
     &  nr,nx_ext,COEFFSEXT,CQ,GM,NQGW,PG,WG,XQ,
     &  BIDOMAIN,EXT_ELEM,FIXQ,INT_ELEM,SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: CALC_FE_COEF
C###  Description:
C###    CALC_FE_COEF is used for the Grid-based Finite Element
C###    scheme. It finds the coeffs of the nonzeros in stiffness
C###    matrix (NQGW) and adds coefficients to global mass matrix (GM).
C###    This code is shared between CALC_FE_LATT_COEF and CALC_FE_GRID_COEF.
C***  Created by Greg Sands 2 September 2003

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc' !NJL_FIBR
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER GB_ELEM(8),m,NGTB,NITB,NNTB,nq,NQGP(0:NQGM),nr,nx_ext
      REAL*8 COEFFSEXT(NQGM),CQ(NMM,NQM),GM(NZ_GM_M),NQGW(NQGM),
     &  PG(NSM,NUM,NGM),WG(NGM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,EXT_ELEM,FIXQ(NYQM,NIYFIXM,NXM),INT_ELEM,
     &  SOLVEEIGHTPROBLEM
!     Local variables
      INTEGER EXT_OFFSET,i,index,INT_OFFSET,ng,ni,nj,nj1,nm,NMT,nu(3),
     &  nn,pj
      REAL*8 CG(NMM),CMAM,kdPsidx(3),kdPsidx_ext(3),dxdNu(3,3),
     &  kdPsidNu(3),kdPsidNu_ext(3),dPsidx(3),dxdXi(3,3),dXidx(3,3),
     &  dXidxdPsidx(3),dXidxdPsidx_ext(3),JACOBIAN,K_INTEGRAND,
     &  K_INTEGRAND_ext,M_INTEGRAND,XG(NJM,NUM)
      LOGICAL INLIST
      DATA nu /2,4,7/

      CALL ENTERS('CALC_FE_COEF',*9999)

C calculate the max number of material parameters (used for this subroutine)
C CM and AM are in first two positions for non solve8 problems
      IF(SOLVEEIGHTPROBLEM) THEN
        EXT_OFFSET = 0
        NMT = 3
      ELSE
        INT_OFFSET = 2
        IF(BIDOMAIN) THEN
          EXT_OFFSET = INT_OFFSET + 3
          NMT = EXT_OFFSET + 3
        ELSE
          NMT = INT_OFFSET + 3
        ENDIF
      ENDIF

C evaluate integrands using gaussian quadrature to get terms for the
C stiffness matrices (NQGW and COEFFSEXT) and global mass matrix (GM)
      DO ng=1,NGTB !gauss point loop

C evaluate material parameters at gauss point
        DO nm=1,NMT
          CG(nm) = 0.0d0
          DO nn=1,NNTB
            CG(nm) = CG(nm) + CQ(nm,GB_ELEM(nn)) * PG(nn,1,ng)
          ENDDO !nn
        ENDDO !nm

        IF(.NOT.SOLVEEIGHTPROBLEM) THEN
          CMAM = CG(1) * CG(2)
        ENDIF

C create dXidx matrix and store Jacobian
        DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,nj1,nr)
          DO ni=1,NITB
            dxdXi(nj,ni) = 0.0d0
            DO nn=1,NNTB
              dxdXi(nj,ni) = dxdXi(nj,ni) + XQ(nj,GB_ELEM(nn)) *
     &          PG(nn,nu(ni),ng)
            ENDDO !nn
          ENDDO !ni
C MLT25Aug05 There are some combinations of (x,y,z) and (xi1,xi2,xi3)
C            where this test is not desirable, e.g. xi2 is aligned with
C            x and xi1 is aligned with y - gives 0 diagonal, but 
C            Jacobian array is invertible. For finite elements built
C            using lattice mappings, better to trap unviable elements 
C            with negative or negligibly small volumes before
C            entering CALC_FE_COEF.
C          IF(DABS(dxdXi(nj,nj)).LT.ZERO_TOL) dxdXi(nj,nj)=1.0d0
        ENDDO !nj1

        CALL INVERT(NITB,dxdXi,dXidx,JACOBIAN)

C create dxdNu (fibre orientation)
        DO nj1=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,nj1,nr)
          XG(nj,1) = 0.0d0
          DO nn=1,NNTB
            XG(nj,1) = XG(nj,1) + XQ(nj,GB_ELEM(nn)) * PG(nn,1,ng)
          ENDDO !nn
        ENDDO !nj1

        CALL MAT_VEC(NITB,nr,dxdNu(1,1),dxdNu(1,2),dxdNu(1,3),dxdXi,XG,
     &    ERROR,*9999)

C calculate dPsidx = dXidx * dPsidXi
        DO nj=1,NITB
          dPsidx(nj) = 0.0d0
          DO ni=1,NITB
            dPsidx(nj) = dPsidx(nj) + dXidx(ni,nj) * PG(m,nu(ni),ng)
          ENDDO !ni
        ENDDO !nj

C!      KAT 21Feb01: can calculate dPsidNu here.
C!      dXidNu would be good for bidomain.

C calculate dXidxdPsidx for the intracellular domain
        IF(INT_ELEM) THEN
C calculate kdPsidNu = d * dxdNu * dPsidx
          DO pj=1,NITB
            kdPsidNu(pj) = 0.0d0
            DO nj=1,NITB
              kdPsidNu(pj) = kdPsidNu(pj) + dxdNu(nj,pj) * dPsidx(nj)
            ENDDO !nj
            kdPsidNu(pj) = kdPsidNu(pj) * CG(INT_OFFSET+pj)
          ENDDO !pj

C calculate kdPsidx = dxdNu * kdPsidNu
          DO nj=1,NITB
            kdPsidx(nj) = 0.0d0
            DO pj=1,NITB
              kdPsidx(nj) = kdPsidx(nj) + dxdNu(nj,pj) * kdPsidNu(pj)
            ENDDO !pj
          ENDDO !nj

C calculate dXidxdPsidx(ni) = dXidx * kdPsidx
          DO ni=1,NITB
            dXidxdPsidx(ni) = 0.0d0
            DO nj=1,NITB
              dXidxdPsidx(ni)=dXidxdPsidx(ni) + dXidx(ni,nj)*kdPsidx(nj)
            ENDDO !nj
          ENDDO !ni
        ENDIF !int_elem

C calculate dXidxdPsidx for the extracellular domain (dXidxdPsidx_ext)
        IF(EXT_ELEM) THEN
C calculate kdPsidNu_ext = d * dxdNu * dPsidx
          DO pj=1,NITB
            kdPsidNu_ext(pj) = 0.0d0
            DO nj=1,NITB
              kdPsidNu_ext(pj)=kdPsidNu_ext(pj)+dxdNu(nj,pj)*dPsidx(nj)
            ENDDO !nj
            kdPsidNu_ext(pj) = kdPsidNu_ext(pj) * CG(EXT_OFFSET+pj)
          ENDDO !pj

C calculate kdPsidx_ext = dxdNu * kdPsidNu_ext
          DO nj=1,NITB
            kdPsidx_ext(nj) = 0.0d0
            DO pj=1,NITB
              kdPsidx_ext(nj) = kdPsidx_ext(nj) + dxdNu(nj,pj) *
     &          kdPsidNu_ext(pj)
            ENDDO !pj
          ENDDO !nj

C calculate dXidxdPsidx_ext(ni) = dXidx * kdPsidx_ext
          DO ni=1,NITB
            dXidxdPsidx_ext(ni) = 0.0d0
            DO nj=1,NITB
              dXidxdPsidx_ext(ni) = dXidxdPsidx_ext(ni) + dXidx(ni,nj) *
     &          kdPsidx_ext(nj)
            ENDDO !nj
          ENDDO !ni
        ENDIF !ext_elem

C evaluate integrands for each nodal position (grid point)
        DO nn=1,NNTB
          IF(INT_ELEM) THEN
            K_INTEGRAND = 0.0d0 !stiffness
            DO ni=1,NITB
              K_INTEGRAND = K_INTEGRAND + PG(nn,nu(ni),ng) *
     &          dXidxdPsidx(ni)
            ENDDO !ni
C integrands are multiplied by gauss weights
            K_INTEGRAND = K_INTEGRAND * JACOBIAN * WG(ng)
            M_INTEGRAND = PG(nn,1,ng) * PG(m,1,ng) * JACOBIAN * WG(ng) !mass
          ENDIF !int_elem

          IF(EXT_ELEM) THEN
            K_INTEGRAND_ext = 0.0d0 !stiffness
            DO ni=1,NITB
              K_INTEGRAND_ext = K_INTEGRAND_ext + PG(nn,nu(ni),ng) *
     &          dXidxdPsidx_ext(ni)
            ENDDO !ni
            K_INTEGRAND_ext = K_INTEGRAND_ext * JACOBIAN * WG(ng)
          ENDIF !ext_elem

C add grid point to NQGP if it isn't there already
C then add weighted integrand (at gauss point) to NQGW and GM
C note: NQGP is used to access COEFFSEXT, NQGW and GM
          INLIST = .FALSE.
          DO i=1,NQGP(0)
            IF (NQGP(i).EQ.GB_ELEM(nn)) THEN
              INLIST = .TRUE.
              index = i
            ENDIF
          ENDDO !i
          IF (.NOT.INLIST) THEN
            NQGP(0) = NQGP(0) + 1
            NQGP(NQGP(0)) = GB_ELEM(nn)
            index = NQGP(0)
          END IF

C add value of integrand to matrix terms
C NOTE: Negative sign is here because there is also a negative sign in
C         ASSEMBLE10_FE due to the sign convention of NQGW for
C         collocation
          IF(INT_ELEM) THEN
            NQGW(index) = NQGW(index) - K_INTEGRAND
C mass matrix only calculated for time dependence (not solve8)
            IF(.NOT.SOLVEEIGHTPROBLEM) THEN
              GM(index+NQGM*(nq-1)) = GM(index+NQGM*(nq-1)) + CMAM *
     &          M_INTEGRAND
            ENDIF
          ENDIF

          IF (BIDOMAIN) THEN
            IF (.NOT.(FIXQ(nq,1,nx_ext).OR.FIXQ(nq,3,nx_ext))) THEN
C COEFFSEXT is the sum of K(intracellular) and K(extracellular)
              IF(INT_ELEM)  COEFFSEXT(index) = 
     &          COEFFSEXT(index) - K_INTEGRAND
              IF (EXT_ELEM) COEFFSEXT(index) =
     &          COEFFSEXT(index) - K_INTEGRAND_ext
            ENDIF !not analytic or essential bc
          ENDIF !bidomain

        ENDDO !nn
      ENDDO !ng

      CALL EXITS('CALC_FE_COEF')
      RETURN
 9999 CALL ERRORS('CALC_FE_COEF',ERROR)
      CALL EXITS('CALC_FE_COEF')
      RETURN 1
      END
