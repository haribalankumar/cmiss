      SUBROUTINE TOFFEL(ICOOR,nb,nr,CHTOFF,DBM,GU,XG,X3G,CURVE,ERROR,*)

C#### Subroutine: TOFFEL
C###  Description:
C###    TOFFEL evaluates Christoffel symbol of 2nd kind CHTOFF(ia,ib,ic)
C###    and, if CURVE is .true., the curvature tensor DBM(ia,ib,ic)
C###    in Xi-coordinate system at current Gauss point.
C**** X3G(nj,nd) are 3rd derivs of Xj (see subroutine X3XG)
C**** G(nj,1) is tangent to surface in Xi(1) direction
C**** G(nj,2) is    "     "    "     " Xi(2)     "
C**** G(nj,3) is    "     "    "      "Xi(3)     "
C**** DG(nj,ib,ic) are cpts of base vector G(nj,ib) derivs wrt Xi(ic)
C**** ICOOR is the global coordinate system type

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER ICOOR,nb,nr
      REAL*8 CHTOFF(3,3,3),DBM(3,3,3),GU(3,3),XG(NJM,NUM),X3G(4,3)
      CHARACTER ERROR*(*)
      LOGICAL CURVE
!     Local Variables
      INTEGER i,ia,ib,ic,id,ig,ij,il,j,jj,k,kj,lj,m,
     '  ND3(2,2,2),NITB,ni,nj,nk,nl,NU1(0:3),NU2(3,3)
      REAL*8 COSHL,COSM,COSP,COST,DDG(3,3,3,3),DG(3,3,3),dXrc_dXref,DZX,
     '  DZXX,DZXXX(3,3,3,3),G(3,3),X(3),
     '  LAMDA,MU,PHI,RAD,SINHL,SINM,SINP,SINT,SUM,SUM1,SUM2,THETA

      DATA NU1/1,2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/
      DATA ND3/1,2,2,3,2,3,3,4/

      CALL ENTERS('TOFFEL',*9999)

C     Initialise all real arrays to zero.
      DO ni=1,3
        DO nj=1,3
          G(nj,ni)=0.0d0
          DO nk=1,3
            DG(nk,nj,ni)=0.0d0
            DBM(ni,nj,nk)=0.0d0
            DO nl=1,3
              DDG(nl,nk,nj,ni)=0.0d0
              DZXXX(nl,nk,nj,ni)=0.0d0
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      NITB=NIT(nb)
!      IF(ITYP10(nr).EQ.1) THEN !rectangular cartesian coords
      IF(ICOOR.EQ.1) THEN !rectangular cartesian coords
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          DO j=1,NITB
            G(nj,j)=XG(nj,NU1(j))
          ENDDO
          DO j=1,NITB
            DO k=1,NITB
              DG(nj,j,k)=XG(nj,NU2(j,k))
            ENDDO
          ENDDO
        ENDDO
        DO i=1,NITB
          DO j=1,NITB
            DO k=1,NITB
              SUM=0.0D0
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                DO m=1,NITB
                  SUM=SUM+DG(nj,j,k)*G(nj,m)*GU(i,m)
                ENDDO
              ENDDO
              IF(DABS(SUM).GT.1.0D-8) THEN
                CHTOFF(i,j,k)=SUM
              ELSE
                CHTOFF(i,j,k)=0.0D0
              ENDIF
            ENDDO
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' CHTOFF('',I1,'','',I1,'',k'''
     '          //''') = '',3E13.5)') i,j,(CHTOFF(i,j,k),k=1,NITB)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDDO
        ENDDO

!      ELSE IF(ITYP10(nr).GT.1) THEN !curvilinear coordinates
      ELSE IF(ICOOR.GT.1) THEN !curvilinear coordinates

C 21/2/97 LC archived section :
C old MPN 24-Apr-96: DZX and DZXX functions in FE01 now used

        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          X(nj)=XG(nj,1)
        ENDDO!nj

        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          DO j=1,NITB
            SUM=0.0D0
            DO ij=1,NJ_LOC(NJL_GEOM,0,nr)
!              dXrc_dXref=DZX(ITYP10(nr),nj,ij,X)
              dXrc_dXref=DZX(ICOOR,nj,ij,X)
              SUM=SUM+dXrc_dXref*XG(ij,NU1(j))
            ENDDO
            G(nj,j)=SUM
            DO k=1,NITB
              SUM=0.0D0
              DO ij=1,NJ_LOC(NJL_GEOM,0,nr)
!                dXrc_dXref=DZX(ITYP10(nr),nj,ij,X)
                dXrc_dXref=DZX(ICOOR,nj,ij,X)
                SUM=SUM+dXrc_dXref*XG(ij,NU2(j,k))
                DO lj=1,NJ_LOC(NJL_GEOM,0,nr)
!                  SUM=SUM+DZXX(ITYP10(nr),nj,ij,lj,X)*
                  SUM=SUM+DZXX(ICOOR,nj,ij,lj,X)*
     '              XG(ij,NU1(j))*XG(lj,NU1(k))
                ENDDO
              ENDDO
              IF(DABS(SUM).GT.1.0D-8) THEN
                DG(nj,j,k)=SUM
              ELSE
                DG(nj,j,k)=0.0D0
              ENDIF
            ENDDO !k
          ENDDO !j
        ENDDO !nj

        DO i=1,NITB
          DO j=1,NITB
            DO k=1,NITB
              SUM=0.0D0
              DO m=1,NITB
                SUM1=0.0D0
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  SUM1=SUM1+DG(nj,j,k)*G(nj,m)
                ENDDO
                SUM=SUM+SUM1*GU(i,m)
              ENDDO !m
              IF(DABS(SUM).GT.1.0D-8) THEN
                CHTOFF(i,j,k)=SUM
              ELSE
                CHTOFF(i,j,k)=0.0D0
              ENDIF
            ENDDO !k
          ENDDO !j
        ENDDO !i
      ENDIF
C
C  If curvature tensors and curvature change tensors are required,
C  find with CURVE=.TRUE.
C
      IF(CURVE) THEN
        IF(NITB.EQ.2) THEN !surface geometry
          G(1,3)=G(2,1)*G(3,2)-G(3,1)*G(2,2)
          G(2,3)=G(3,1)*G(1,2)-G(1,1)*G(3,2)
          G(3,3)=G(1,1)*G(2,2)-G(2,1)*G(1,2)
          SUM=0.0D0
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            SUM=SUM+G(nj,3)**2
          ENDDO
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            G(nj,3)=G(nj,3)/DSQRT(SUM) !sqrt(sum)=rg
          ENDDO
          DO j=1,NITB
            DO k=1,NITB
              SUM=0.0D0
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                SUM=SUM+DG(nj,j,k)*G(nj,3)
              ENDDO
              IF(DABS(SUM).GT.1.0D-8) THEN
                CHTOFF(3,j,k)=SUM
              ELSE
                CHTOFF(3,j,k)=0.0D0
              ENDIF
            ENDDO !k
          ENDDO !j

!          IF(ITYP10(nr).EQ.1) THEN !rectangular cartesian coords
          IF(ICOOR.EQ.1) THEN !rectangular cartesian coords
            DO nj=1,NITB
              DO jj=1,NITB
                DO lj=1,NITB
                  DO kj=1,NITB
                    DDG(nj,jj,lj,kj)=X3G(ND3(jj,lj,kj),nj)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!          ELSE IF(ITYP10(nr).GT.1) THEN !curvilinear coordinates
          ELSE IF(ICOOR.GT.1) THEN !curvilinear coordinates

C!!!        DZXXX should be put into a function like DZXX in FE01

!            IF(ITYP10(nr).EQ.2) THEN !cylindrical polar coords
            IF(ICOOR.EQ.2) THEN !cylindrical polar coords
              RAD=XG(1,1)
              THETA=XG(2,1)
              SINT=DSIN(THETA)
              COST=DCOS(THETA)
              DZXXX(1,2,1,2)=-COST
              DZXXX(2,2,1,2)=-SINT
              DZXXX(1,2,2,1)=-COST
              DZXXX(2,2,2,1)=-SINT
              DZXXX(1,1,2,2)=-COST
              DZXXX(1,2,2,2)=RAD*SINT
              DZXXX(2,1,2,2)=-SINT
              DZXXX(2,2,2,2)=-RAD*COST

!            ELSE IF(ITYP10(nr).EQ.3) THEN !spherical polar coords
            ELSE IF(ICOOR.EQ.3) THEN !spherical polar coords
              RAD  =XG(1,1)
              THETA=XG(2,1)
              PHI  =XG(3,1)
              SINT=DSIN(THETA)
              COST=DCOS(THETA)
              SINP=DSIN(PHI)
              COSP=DCOS(PHI)
              DZXXX(1,2,1,2)=-COST*COSP
              DZXXX(1,3,1,2)=SINT*SINP
              DZXXX(2,2,1,2)=-SINT*COSP
              DZXXX(2,3,1,2)=-COST*SINP
              DZXXX(1,2,1,3)=SINT*SINP
              DZXXX(1,3,1,3)=-COST*COSP
              DZXXX(2,2,1,3)=-COST*SINP
              DZXXX(2,3,1,3)=-SINT*COSP
              DZXXX(3,3,1,3)=-SINT
              DZXXX(1,2,2,1)=-COST*COSP
              DZXXX(1,3,2,1)=SINT*SINP
              DZXXX(2,2,2,1)=-SINT*COSP
              DZXXX(2,3,2,1)=-COST*SINP
              DZXXX(1,1,2,2)=-COST*COSP
              DZXXX(1,2,2,2)=RAD*SINT*COSP
              DZXXX(1,3,2,2)=RAD*COST*SINP
              DZXXX(2,1,2,2)=-SINT*COSP
              DZXXX(2,2,2,2)=-RAD*COST*COSP
              DZXXX(2,3,2,2)=RAD*SINT*SINP
              DZXXX(1,1,2,3)=SINT*SINP
              DZXXX(1,2,2,3)=RAD*COST*SINP
              DZXXX(1,3,2,3)=RAD*SINT*COSP
              DZXXX(2,1,2,3)=-COST*SINP
              DZXXX(2,2,2,3)=RAD*SINT*SINP
              DZXXX(2,3,2,3)=-RAD*COST*COSP
              DZXXX(1,2,3,1)=SINT*SINP
              DZXXX(1,3,3,1)=-COST*COSP
              DZXXX(2,2,3,1)=-COST*SINP
              DZXXX(2,3,3,1)=-SINT*COSP
              DZXXX(3,3,3,1)=-SINP
              DZXXX(1,1,3,2)=SINT*SINP
              DZXXX(1,2,3,2)=RAD*COST*SINP
              DZXXX(1,3,3,2)=RAD*SINT*COSP
              DZXXX(2,1,3,2)=-COST*SINP
              DZXXX(2,2,3,2)=RAD*SINT*SINP
              DZXXX(2,3,3,2)=RAD*COST*COSP
              DZXXX(1,1,3,3)=-COST*COSP
              DZXXX(1,2,3,3)=RAD*SINT*COSP
              DZXXX(1,3,3,3)=RAD*COST*SINP
              DZXXX(2,1,3,3)=-SINT*COSP
              DZXXX(2,2,3,3)=-RAD*COST*COSP
              DZXXX(2,3,3,3)=RAD*SINT*SINP
              DZXXX(3,1,3,3)=-SINP
              DZXXX(3,3,3,3)=-RAD*COSP

!            ELSE IF(ITYP10(nr).EQ.4) THEN !prolate spheroidal coords
            ELSE IF(ICOOR.EQ.4) THEN !prolate spheroidal coords
              LAMDA=XG(1,1)
              MU=XG(2,1)
              THETA=XG(3,1)
              SINHL=DSINH(LAMDA)
              COSHL=DCOSH(LAMDA)
              SINM=DSIN(MU)
              COSM=DCOS(MU)
              SINT=DSIN(THETA)
              COST=DCOS(THETA)
              DZXXX(1,1,1,1)=FOCUS*SINHL*COSM
              DZXXX(1,2,1,1)=-FOCUS*COSHL*SINM
              DZXXX(2,1,1,1)=FOCUS*COSHL*SINM*COST
              DZXXX(2,2,1,1)=FOCUS*SINHL*COSM*COST
              DZXXX(2,3,1,1)=-FOCUS*SINHL*SINM*SINT
              DZXXX(3,1,1,1)=FOCUS*COSHL*SINM*SINT
              DZXXX(3,2,1,1)=FOCUS*SINHL*COSM*SINT
              DZXXX(3,3,1,1)=FOCUS*SINHL*SINM*COST
              DZXXX(1,1,1,2)=-FOCUS*COSHL*SINM
              DZXXX(1,2,1,2)=-FOCUS*SINHL*COSM
              DZXXX(2,1,1,2)=FOCUS*SINHL*COSM*COST
              DZXXX(2,2,1,2)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,3,1,2)=-FOCUS*COSHL*COSM*SINT
              DZXXX(3,1,1,2)=FOCUS*SINHL*COSM*SINT
              DZXXX(3,2,1,2)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,3,1,2)=FOCUS*COSHL*COSM*COST
              DZXXX(2,1,1,3)=-FOCUS*SINHL*SINM*SINT
              DZXXX(2,2,1,3)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,3,1,3)=-FOCUS*COSHL*SINM*COST
              DZXXX(3,1,1,3)=FOCUS*SINHL*SINM*COST
              DZXXX(3,2,1,3)=FOCUS*COSHL*COSM*COST
              DZXXX(3,3,1,3)=-FOCUS*COSHL*SINM*SINT
              DZXXX(1,1,2,1)=-FOCUS*COSHL*SINM
              DZXXX(1,2,2,1)=-FOCUS*SINHL*COSHL
              DZXXX(2,1,2,1)=FOCUS*SINHL*COSM*COST
              DZXXX(2,2,2,1)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,3,2,1)=-FOCUS*COSHL*COSM*SINT
              DZXXX(3,1,2,1)=FOCUS*SINHL*COSM*SINT
              DZXXX(3,2,2,1)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,3,2,1)=FOCUS*COSHL*COSM*COST
              DZXXX(1,1,2,2)=-FOCUS*SINHL*COSM
              DZXXX(1,2,2,2)=FOCUS*COSHL*SINM
              DZXXX(2,1,2,2)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,2,2,2)=-FOCUS*SINHL*COSM*COST
              DZXXX(2,3,2,2)=FOCUS*SINHL*SINM*SINT
              DZXXX(3,1,2,2)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,2,2,2)=-FOCUS*SINHL*COSM*SINT
              DZXXX(3,3,2,2)=-FOCUS*SINHL*SINM*COST
              DZXXX(2,1,2,3)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,2,2,3)=FOCUS*SINHL*SINM*SINT
              DZXXX(2,3,2,3)=-FOCUS*SINHL*COSM*COST
              DZXXX(3,1,2,3)=FOCUS*COSHL*COSM*COST
              DZXXX(3,2,2,3)=-FOCUS*SINHL*SINM*COST
              DZXXX(3,3,2,3)=-FOCUS*SINHL*COSM*SINT
              DZXXX(2,1,3,1)=-FOCUS*SINHL*SINM*SINT
              DZXXX(2,2,3,1)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,3,3,1)=-FOCUS*COSHL*SINM*COST
              DZXXX(3,1,3,1)=FOCUS*SINHL*SINM*COST
              DZXXX(3,2,3,1)=FOCUS*COSHL*COSM*COST
              DZXXX(3,3,3,1)=-FOCUS*COSHL*SINM*SINT
              DZXXX(2,1,3,2)=-FOCUS*COSHL*COSM*SINT
              DZXXX(2,2,3,2)=FOCUS*SINHL*SINM*SINT
              DZXXX(2,3,3,2)=-FOCUS*SINHL*COSM*COST
              DZXXX(3,1,3,2)=FOCUS*COSHL*COSM*COST
              DZXXX(3,2,3,2)=-FOCUS*SINHL*SINM*COST
              DZXXX(3,3,3,2)=-FOCUS*SINHL*COSM*SINT
              DZXXX(2,1,3,3)=-FOCUS*COSHL*SINM*COST
              DZXXX(2,2,3,3)=-FOCUS*SINHL*COSM*COST
              DZXXX(2,3,3,3)=FOCUS*SINHL*SINM*SINT
              DZXXX(3,1,3,3)=-FOCUS*COSHL*SINM*SINT
              DZXXX(3,2,3,3)=-FOCUS*SINHL*COSM*SINT
              DZXXX(3,3,3,3)=-FOCUS*SINHL*SINM*COST
            ENDIF !icoor

            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              DO ia=1,2
                DO ib=1,2
                  DO ic=1,2
                    SUM=0.0D0
                    DO ij=1,NJ_LOC(NJL_GEOM,0,nr)
!                      dXrc_dXref=DZX(ITYP10(nr),nj,ij,X)
                      dXrc_dXref=DZX(ICOOR,nj,ij,X)
                      SUM=SUM+dXrc_dXref*X3G(ND3(ia,ib,ic),ij)
                      DO kj=1,NJ_LOC(NJL_GEOM,0,nr)
!                        SUM=SUM+DZXX(ITYP10(nr),nj,ij,kj,X)*
                        SUM=SUM+DZXX(ICOOR,nj,ij,kj,X)*
     '                    (XG(ij,NU2(ia,ic))*
     '                    XG(kj,NU1(ib))+XG(ij,NU1(ia))*
     '                    XG(kj,NU2(ib,ic))+XG(kj,NU1(ic))*
     '                    XG(ij,NU2(ia,ib)))
                        DO jj=1,NJ_LOC(NJL_GEOM,0,nr)
                          SUM=SUM+DZXXX(nj,ij,jj,kj)*(XG(kj,NU1(ic))
     '                      *XG(ij,NU1(ia))*XG(jj,NU1(ib)))
                        ENDDO !jj
                      ENDDO !kj
                    ENDDO !ij
                    IF(DABS(SUM).LT.1.0D-8) THEN
                      DDG(nj,ia,ib,ic)=0.0D0
                    ELSE
                      DDG(nj,ia,ib,ic)=SUM
                    ENDIF
                  ENDDO !ic
                ENDDO !ib
              ENDDO !ia
            ENDDO !nj
          ENDIF !rc/curvilinear coords

          DO ia=1,NITB
            DO ib=1,NITB
              DO ig=1,NITB
                SUM=0.0d0
                SUM2=0.0d0
                DO il=1,2
                  SUM1=0.0d0
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    SUM1=SUM1+G(nj,3)*DDG(nj,il,ib,ig)
                  ENDDO
                  SUM2=SUM2+GU(ia,il)*SUM1
                  DO id=1,2
                    SUM=SUM-GU(ia,il)*(CHTOFF(id,il,ig)*CHTOFF(3,id,ib)
     '                +CHTOFF(id,ig,ib)*CHTOFF(3,id,il)
     '                +CHTOFF(id,ib,il)*CHTOFF(3,id,ig))
                  ENDDO
                ENDDO
                IF(DABS(SUM+SUM2).GT.1.0d-8) THEN
                  DBM(ia,ib,ig)=SUM+SUM2
                ELSE
                  DBM(ia,ib,ig)=0.0d0
                ENDIF
              ENDDO !ig
            ENDDO !ib
          ENDDO !ia
        ENDIF !nit=2
      ENDIF !curve

      CALL EXITS('TOFFEL')
      RETURN
 9999 CALL ERRORS('TOFFEL',ERROR)
      CALL EXITS('TOFFEL')
      RETURN 1
      END


