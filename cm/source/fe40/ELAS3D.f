      SUBROUTINE ELAS3D(NBH,NBJ,nr,nw,nx,CG,
     '  ED,EM,ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*)

C#### Subroutine: ELAS3D
C###  Description:
C###    ELAS3D calculates the element load vector ER and the stiffness
C###    matrix ES for 3D elasticity stress (nw=9).
C GMH 3/1/96 Major changes to allow for fibre directions and therefore
C            anisotropic materials.
C            For working, see CMISS linear elasticity fibre working
C            sheet of PJH.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b15.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),nr,nw,nx
      REAL*8 CG(NMM,NGM),ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),
     '  ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM),WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k,k1,k2,ms,nb,ng,NGTB,ns,ns1,ns2,
     '  ns_k1,ns_k2,NSTB
      REAL*8 B(3,3,3,64),CC(6,6),C1,C2,C3,DENS,DPSINU(3),
     '  DXIX(3,3),DXINU(3,3),
     '  GL(3,3),GLNU(3,3),GU(3,3),GUNU(3,3),
     '  RT(3,3),RWG,S_thermal(3)
      LOGICAL ROTATE_FIBRE

      CALL ENTERS('ELAS3D',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' 3D elasticity'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      !Ensure 3D as all arrays are dimensioned to 3
      nb=NBH(NH_LOC(1,nx)) !is dependent variable basis
      NGTB=NGT(nb)
      NSTB=NST(nb)
C GMH 3/1/96 Should B be dimensioned in FE20 and passed thru?
      CALL ASSERT(NSTB.LE.64,'Last index of B insufficient',ERROR,*9999)

      ROTATE_FIBRE=.FALSE.
      IF(IMT(nw).EQ.5) THEN !Special case anisotropic
        IF(NJ_LOC(NJL_FIBR,0,nr).GT.0) THEN
          ROTATE_FIBRE=.TRUE.
        ENDIF
      ENDIF
      !if fibre angles, then set RT to the identity
      IF(.NOT.ROTATE_FIBRE) THEN
        RT(1,1)=1.0D0
        RT(1,2)=0.0D0
        RT(1,3)=0.0D0
        RT(2,1)=0.0D0
        RT(2,2)=1.0D0
        RT(2,3)=0.0D0
        RT(3,1)=0.0D0
        RT(3,2)=0.0D0
        RT(3,3)=1.0D0
      ENDIF

C GMH 3/1/96 The following is calculated according to material
C            coordinates.  If no material coordinates, then the
C            material coordinates drop down to being the coordinate
C            directions.  A bit more elegant this way.
C            For working, see CMISS linear elasticity fibre working
C            sheet of PJH.
      DO ng=1,NGTB
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        !get material coord info in DXIXN etc
        !if not fibre, then DXIXN==DXIX
        IF(ROTATE_FIBRE) THEN !Calculate wrt nu coords
          CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXINU,GLNU,GUNU,
     '      RG(ng),XG,ERROR,*9999)
        ELSE
          CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXINU,GLNU,GUNU,
     '      RG(ng),XG,ERROR,*9999)
        ENDIF
        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
        !Calculate B from e=Bu, where e is in material coordinates,
        !  and u is the geometric displacements.
        DO ns=1,NSTB
          IF(ROTATE_FIBRE) THEN !Calculate rotation matrix at gauss pt
            CALL MAT_VEC_NG(3,nr,RT(1,1),RT(1,2),RT(1,3),XG,
     '        ERROR,*9999)
          ENDIF
          DPSINU(1)=PG(ns,2,ng,nb)*DXINU(1,1)
     '      +PG(ns,4,ng,nb)*DXINU(2,1)
     '      +PG(ns,7,ng,nb)*DXINU(3,1)
          DPSINU(2)=PG(ns,2,ng,nb)*DXINU(1,2)
     '      +PG(ns,4,ng,nb)*DXINU(2,2)
     '      +PG(ns,7,ng,nb)*DXINU(3,2)
          DPSINU(3)=PG(ns,2,ng,nb)*DXINU(1,3)
     '      +PG(ns,4,ng,nb)*DXINU(2,3)
     '      +PG(ns,7,ng,nb)*DXINU(3,3)
          !Inefficient - should test for ROTATE_FIBRE (when RT=I)
          DO i=1,3
            DO j=1,3
              DO k=1,3
                B(i,j,k,ns)=0.5D0*(RT(k,i)*DPSINU(j)+RT(k,j)*DPSINU(i))
              ENDDO
            ENDDO !j
          ENDDO !i
        ENDDO !ms
        RWG=RG(ng)*WG(ng,nb)
        CALL CGC9(nr,nw,nx,CC,C1,C2,C3,CG(1,ng),DENS,S_thermal,ERROR,
     '    *9999)
! Isotropic or Anisotropic (special case)
        IF(IMT(nw).EQ.1.OR.IMT(nw).EQ.5) THEN
          DO k1=1,3
            ns_k1=(k1-1)*NSTB
            DO ns1=1,NSTB
              !Thermal is only implemented for isotropic, so it should
              !be wrt element coords
              IF(KTYP43.GT.0) THEN !include thermal stress effects
                ER(ns_k1+ns1) = ER(ns_k1+ns1) - S_thermal(k1)*
     '            (PG(ns1,2,ng,nb)*DXIX(1,k1)
     '            +PG(ns1,4,ng,nb)*DXIX(2,k1)
     '            +PG(ns1,7,ng,nb)*DXIX(3,k1))*RWG
              ENDIF !ktyp43
              DO k2=1,3
                ns_k2=(k2-1)*NSTB
                DO ns2=1,NSTB
C GMH 26/12/95 Should this be wrt element or material coords???
                  IF(k1.EQ.k2) THEN !mass lumping, so on diagonals
                    EM(ns_k1+ns1,ns_k2+ns2) =EM(ns_k1+ns1,ns_k2+ns2)+
     '                DENS*PG(ns1,1,ng,nb)*PG(ns2,1,ng,nb)*RWG
                  ENDIF
                  ES(ns_k1+ns1,ns_k2+ns2)=ES(ns_k1+ns1,ns_k2+ns2)+RWG*(
     '              B(1,1,k1,ns1)*(C1*B(1,1,k2,ns2)+C2*B(2,2,k2,ns2)+
     '              C2*B(3,3,k2,ns2))+
     '              B(2,2,k1,ns1)*(C2*B(1,1,k2,ns2)+C1*B(2,2,k2,ns2)+
     '              C2*B(3,3,k2,ns2))+
     '              B(3,3,k1,ns1)*(C2*B(1,1,k2,ns2)+C2*B(2,2,k2,ns2)+
     '              C1*B(3,3,k2,ns2))+
     '              2.0D0*C3*(B(1,2,k1,ns1)*B(1,2,k2,ns2)+
     '              B(1,3,k1,ns1)*B(1,3,k2,ns2)+
     '              B(2,3,k1,ns1)*B(2,3,k2,ns2)))
                ENDDO !ns2
              ENDDO !k2
            ENDDO !ns1
          ENDDO !k1
! Transversly isotropic (fibre in 1 dir.n)
        ELSE IF(IMT(nw).EQ.2) THEN
          ERROR='>>3D Transversly isotropic not implemented'
          GOTO 9999
! Transversly isotropic (fibre in 2 dir.n)
        ELSE IF(IMT(nw).EQ.3) THEN
          ERROR='>>3D Transversly isotropic not implemented'
          GOTO 9999
! Orthotropic (fibre in 1 dir.n)
        ELSE IF(IMT(nw).EQ.4) THEN
          ERROR='>>3D orthotropic not implemented'
          GOTO 9999
        ELSE
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF !anisotropy
      ENDDO !ng

      DO k1=1,3
        ns_k1=(k1-1)*NSTB
        DO ns1=1,NSTB
          ER(ns_k1+ns1) = ER(ns_k1+ns1) *SE(ns1,nb)
          DO k2=1,3
            ns_k2=(k2-1)*NSTB
            DO ns2=1,NSTB
              ES(ns_k1+ns1,ns_k2+ns2)=ES(ns_k1+ns1,ns_k2+ns2)*
     '          SE(ns1,nb)*SE(ns2,nb)
            ENDDO !ns2
          ENDDO !k2
        ENDDO !ns1
      ENDDO !k1

      IF(DAMPING_FACTOR1.GT.0.d0.OR.DAMPING_FACTOR2.GT.0.d0) THEN
C GMH 3/1/96 Why is this loop to 2*NSTB
        DO ms=1,2*NSTB
          DO ns=1,2*NSTB
            ED(ms,ns)=DAMPING_FACTOR1*DAMPING_FREQUENCY*EM(ms,ns)
     '               +DAMPING_FACTOR2*DAMPING_FREQUENCY*ES(ms,ns)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('ELAS3D')
      RETURN
 9999 CALL ERRORS('ELAS3D',ERROR)
      CALL EXITS('ELAS3D')
      RETURN 1
      END


