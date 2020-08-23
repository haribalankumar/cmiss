      SUBROUTINE CALC_GRID_COEF(NENQ,NITB,nq,NQGP,NQS,
     '  NQXI,NWQ,NXQ,nx_ext,nx_trans,
     '  COEFFSEXT,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,NQGW,PROPQ,
     '  BIDOMAIN,FIXQ,IMPLICIT,SOLVEEIGHTPROBLEM,ERROR,*)

C#### Subroutine: CALC_GRID_COEF
C###  Description:
C###    CALC_GRID_COEF finds the coeffs of the nonzeros in a matrix NQGW
C###    formed for the explicit or implicit soln of a grid pt scheme.
C###    The number of nonzeros depends on the number of local xi
C###    coords. It is calculated for a single nq.
C###    If adaptive (NMGT>1) the coeffs depend on the grid level na.
C***  Created by Martin Buist, May 1997

      IMPLICIT NONE

C      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'

!     Parameter list
      INTEGER NENQ(0:8,NQM),NITB,nq,NQGP(0:NQGM,NQM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),NWQ(8),
     '  NXQ(-NIM:NIM,0:4,0:NQM),nx_ext,nx_trans
      REAL*8 COEFFSEXT(NQGM),CQ(NMM,NQM),DNUDXQ(3,3,NQM),
     '  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),GCHQ(3),GUQ(3,3),NQGW(NQGM),
     '  PROPQ(3,3,4,2)
      CHARACTER ERROR*(*)
      LOGICAL BIDOMAIN,FIXQ(NYQM,NIYFIXM,NXM),IMPLICIT,SOLVEEIGHTPROBLEM
!     Local variables
      INTEGER i,n,n2,ni,nii,nij
      REAL*8 SUM(15)

      CALL ENTERS('CALC_GRID_COEF',*9999)

C Now added in assemble10
C      DTCMAM=DT/((CQ(1,nq)*1.0D-6)*CQ(2,nq)) !10^-6 correction for Cm

      IF(BIDOMAIN) THEN
        DO i=1,NQGM
          COEFFSEXT(i)=0.0d0
        ENDDO
      ENDIF !bidomain

      n=0
      n2=0
      IF(NWQ(1).EQ.0) THEN !internal gp

        IF(NITB.EQ.1) THEN
          DO i=1,NXQ(-1,0,nq)
            n=n+1
            NQGW(n)=-(GUQ(1,1)*PROPQ(1,1,2,1))+
     '        (GCHQ(1)*PROPQ(1,1,1,1))+
     '        (4.0d0*GUQ(1,1)*PROPQ(1,1,1,1))
            NQGW(n)=NQGW(n)/DBLE(NXQ(-1,0,nq))
            IF(BIDOMAIN) THEN
              COEFFSEXT(n)=-(GUQ(1,1)*(PROPQ(1,1,2,1)+
     '          PROPQ(1,1,2,2)))+
     '          (GCHQ(1)*(PROPQ(1,1,1,1)+PROPQ(1,1,1,2)))+
     '          (4.0d0*GUQ(1,1)*(PROPQ(1,1,1,2)+
     '          PROPQ(1,1,1,1)))
              COEFFSEXT(n)=COEFFSEXT(n)/DBLE(NXQ(-1,0,nq))
            ENDIF
          ENDDO !i

          n=n+1
          NQGW(n)=-8.0d0*GUQ(1,1)*PROPQ(1,1,1,1)
          IF(BIDOMAIN) THEN
            COEFFSEXT(n)=-8.0d0*GUQ(1,1)*(PROPQ(1,1,1,2)+
     '        PROPQ(1,1,1,1))
          ENDIF !bidomain

          DO i=1,NXQ(1,0,nq)
            n=n+1
            NQGW(n)=-(GCHQ(1)*PROPQ(1,1,1,1))+
     '        (GUQ(1,1)*PROPQ(1,1,2,1))+
     '        (4.0d0*GUQ(1,1)*PROPQ(1,1,1,1))
            NQGW(n)=NQGW(n)/DBLE(NXQ(1,0,nq))
            IF(BIDOMAIN) THEN
              COEFFSEXT(n)=-(GCHQ(1)*(PROPQ(1,1,1,1)+
     '          PROPQ(1,1,1,2)))+
     '          (GUQ(1,1)*(PROPQ(1,1,2,1)+PROPQ(1,1,2,2)))+
     '          (4.0d0*GUQ(1,1)*(PROPQ(1,1,1,2)+
     '          PROPQ(1,1,1,1)))
              COEFFSEXT(n)=COEFFSEXT(n)/DBLE(NXQ(1,0,nq))
            ENDIF
          ENDDO !i

        ELSE IF(NITB.EQ.2) THEN
          DO ni=1,8
            SUM(ni)=0.0d0
          ENDDO !ni
          DO nii=1,NITB
!         C(1j) * G1j (upper)
            SUM(1)=SUM(1)+GUQ(nii,1)*PROPQ(1,nii,1,1)
!         C(2j) * G2j (upper)
            SUM(2)=SUM(2)+GUQ(nii,2)*PROPQ(2,nii,1,1)
!         C(1l) * [Christoffel*Gij (upper)](l)
            SUM(3)=SUM(3)+GCHQ(nii)*PROPQ(1,nii,1,1)
!         C(2l) * [Christoffel*Gij (upper)](l)
            SUM(4)=SUM(4)+GCHQ(nii)*PROPQ(2,nii,1,1)
!         C(1j) * G2j (upper)
            SUM(5)=SUM(5)+GUQ(nii,1)*PROPQ(2,nii,1,1)
!         C(2j) * G1j (upper)
            SUM(6)=SUM(6)+GUQ(nii,2)*PROPQ(1,nii,1,1)
            DO nij=1,NITB
!           C(1k),j * Gij (upper)
              SUM(7)=SUM(7)+GUQ(nii,nij)*PROPQ(1,nii,1+nij,1)
!           C(2k),j * Gij (upper)
              SUM(8)=SUM(8)+GUQ(nii,nij)*PROPQ(2,nii,1+nij,1)
            ENDDO !nij
          ENDDO !nii

          n=n+1
          NQGW(n)=SUM(5)+SUM(6)
          DO i=1,NXQ(-2,0,nq)
            n=n+1
            NQGW(n)=(-SUM(8)+SUM(4)+4.0d0*SUM(2))/DBLE(NXQ(-2,0,nq))
          ENDDO !i
          n=n+1
          NQGW(n)=-SUM(5)-SUM(6)
          DO i=1,NXQ(-1,0,nq)
            n=n+1
            NQGW(n)=(-SUM(7)+SUM(3)+4.0d0*SUM(1))/DBLE(NXQ(-1,0,nq))
          ENDDO
          n=n+1
          NQGW(n)=-8.0d0*(SUM(1)+SUM(2))
          DO i=1,NXQ(1,0,nq)
            n=n+1
            NQGW(n)=(SUM(7)-SUM(3)+4.0d0*SUM(1))/DBLE(NXQ(1,0,nq))
          ENDDO
          n=n+1
          NQGW(n)=-SUM(5)-SUM(6)
          DO i=1,NXQ(2,0,nq)
            n=n+1
            NQGW(n)=(SUM(8)-SUM(4)+4.0d0*SUM(2))/DBLE(NXQ(2,0,nq))
          ENDDO
          n=n+1
          NQGW(n)=SUM(5)+SUM(6)

          IF(BIDOMAIN) THEN
            DO ni=1,8
              SUM(ni)=0.0d0
            ENDDO !ni
            DO nii=1,NITB
!           (C(1j)int + C(1j)ext) * G1j(upper)
              SUM(1)=SUM(1)+GUQ(nii,1)*(PROPQ(1,nii,1,1)+
     '          PROPQ(1,nii,1,2))
!           (C(2j)int + C(2j)ext) * G2j(upper)
              SUM(2)=SUM(2)+GUQ(nii,2)*(PROPQ(2,nii,1,1)+
     '          PROPQ(2,nii,1,2))
!           (C(1l)int + C(1l)ext) * [Christoffel*Gij(upper)](l)
              SUM(3)=SUM(3)+GCHQ(nii)*(PROPQ(1,nii,1,1)+
     '          PROPQ(1,nii,1,2))
!           (C(2l)int + C(2l)ext) * [Christoffel*Gij (upper)](l)
              SUM(4)=SUM(4)+GCHQ(nii)*(PROPQ(2,nii,1,1)+
     '          PROPQ(2,nii,1,2))
!           (C(1j)int + C(1j)ext) * G2j (upper)
              SUM(5)=SUM(5)+GUQ(nii,1)*(PROPQ(2,nii,1,1)+
     '          PROPQ(2,nii,1,2))
!           (C(2j)int + C(2j)ext) * G1j (upper)
              SUM(6)=SUM(6)+GUQ(nii,2)*(PROPQ(1,nii,1,1)+
     '          PROPQ(1,nii,1,2))
              DO nij=1,NITB
!             (C(1k),j int + C(1k),j ext) * Gij (upper)
                SUM(7)=SUM(7)+GUQ(nii,nij)*(PROPQ(1,nii,1+nij,1)+
     '            PROPQ(1,nii,1+nij,2))
!             (C(2k),j int + C(2k),j ext) * Gij (upper)
                SUM(8)=SUM(8)+GUQ(nii,nij)*(PROPQ(2,nii,1+nij,1)+
     '            PROPQ(2,nii,1+nij,2))
              ENDDO !nij
            ENDDO !nii

            n2=n2+1
            COEFFSEXT(n2)=SUM(5)+SUM(6)
            DO i=1,NXQ(-2,0,nq)
              n2=n2+1
              COEFFSEXT(n2)=(-SUM(8)+SUM(4)+4.0d0*SUM(2))/
     '          DBLE(NXQ(-2,0,nq))
            ENDDO !i
            n2=n2+1
            COEFFSEXT(n2)=-SUM(5)-SUM(6)
            DO i=1,NXQ(-1,0,nq)
              n2=n2+1
              COEFFSEXT(n2)=(-SUM(7)+SUM(3)+4.0d0*SUM(1))/
     '          DBLE(NXQ(-1,0,nq))
            ENDDO !i
            n2=n2+1
            COEFFSEXT(n2)=-8.0d0*(SUM(1)+SUM(2))
            DO i=1,NXQ(1,0,nq)
              n2=n2+1
              COEFFSEXT(n2)=(SUM(7)-SUM(3)+4.0d0*SUM(1))/
     '          DBLE(NXQ(1,0,nq))
            ENDDO !i
            n2=n2+1
            COEFFSEXT(n2)=-SUM(5)-SUM(6)
            DO i=1,NXQ(2,0,nq)
              n2=n2+1
              COEFFSEXT(n2)=(SUM(8)-SUM(4)+4.0d0*SUM(2))/
     '          DBLE(NXQ(2,0,nq))
            ENDDO !i
            n2=n2+1
            COEFFSEXT(n2)=SUM(5)+SUM(6)
          ENDIF !bidomain

        ELSE IF(NITB.EQ.3) THEN
          DO ni=1,15
            SUM(ni)=0.0d0
          ENDDO !ni
          DO nii=1,NITB
!         C(1j) * G1j(upper)
            SUM(1)=SUM(1)+GUQ(nii,1)*PROPQ(1,nii,1,1)
!         C(2j) * G2j(upper)
            SUM(2)=SUM(2)+GUQ(nii,2)*PROPQ(2,nii,1,1)
!         C(3j) * G3j(upper)
            SUM(3)=SUM(3)+GUQ(nii,3)*PROPQ(3,nii,1,1)
!         C(1l) * [Christoffel*Gij(upper)](l)
            SUM(4)=SUM(4)+GCHQ(nii)*PROPQ(1,nii,1,1)
!         C(2l) * [Christoffel*Gij(upper)](l)
            SUM(5)=SUM(5)+GCHQ(nii)*PROPQ(2,nii,1,1)
!         C(3l) * [Christoffel*Gij(upper)](l)
            SUM(6)=SUM(6)+GCHQ(nii)*PROPQ(3,nii,1,1)
            DO nij=1,NITB
!           C(1k),j * Gij(upper)
              SUM(7)=SUM(7)+GUQ(nii,nij)*PROPQ(1,nii,1+nij,1)
!           C(2k),j * Gij(upper)
              SUM(8)=SUM(8)+GUQ(nii,nij)*PROPQ(2,nii,1+nij,1)
!           C(3k),j * Gij(upper)
              SUM(9)=SUM(9)+GUQ(nii,nij)*PROPQ(3,nii,1+nij,1)
            ENDDO !nij
!         C(1j) * G2j(upper)
            SUM(10)=SUM(10)+GUQ(nii,1)*PROPQ(2,nii,1,1)
!         C(2j) * G1j(upper)
            SUM(11)=SUM(11)+GUQ(nii,2)*PROPQ(1,nii,1,1)
!         C(1j) * G3j(upper)
            SUM(12)=SUM(12)+GUQ(nii,1)*PROPQ(3,nii,1,1)
!         C(3j) * G1j(upper)
            SUM(13)=SUM(13)+GUQ(nii,3)*PROPQ(1,nii,1,1)
!         C(2j) * G3j(upper)
            SUM(14)=SUM(14)+GUQ(nii,2)*PROPQ(3,nii,1,1)
!         C(3j) * G2j(upper)
            SUM(15)=SUM(15)+GUQ(nii,3)*PROPQ(2,nii,1,1)
          ENDDO !nii

          IF((NXQ(-3,0,nq).LE.1).AND.(NXQ(-2,0,nq).LE.1).AND.
     '      (NXQ(-1,0,nq).LE.1).AND.(NXQ(1,0,nq).LE.1).AND.
     '      (NXQ(2,0,nq).LE.1).AND.(NXQ(3,0,nq).LE.1)) THEN

            !Usual case, no element splitting
            NQGW(1)=SUM(14)+SUM(15)                !xi -2 -3
            NQGW(2)=SUM(12)+SUM(13)                !xi -1 -3
            NQGW(3)=-SUM(9)+SUM(6)+4.0d0*SUM(3)    !xi  0 -3
            NQGW(4)=-SUM(12)-SUM(13)               !xi  1 -3
            NQGW(5)=-SUM(14)-SUM(15)               !xi  2 -3
            NQGW(6)=SUM(10)+SUM(11)                !xi -1 -2
            NQGW(7)=-SUM(8)+SUM(5)+4.0d0*SUM(2)    !xi  0 -2
            NQGW(8)=-SUM(10)-SUM(11)               !xi  1 -2
            NQGW(9)=-SUM(7)+SUM(4)+4.0d0*SUM(1)    !xi  0 -1
            NQGW(10)=-8.0d0*(SUM(1)+SUM(2)+SUM(3)) !xi  0  0
            NQGW(11)=SUM(7)-SUM(4)+4.0d0*SUM(1)    !xi  0  1
            NQGW(12)=-SUM(10)-SUM(11)              !xi -1  2
            NQGW(13)=SUM(8)-SUM(5)+4.0d0*SUM(2)    !xi  0  2
            NQGW(14)=SUM(10)+SUM(11)               !xi  1  2
            NQGW(15)=-SUM(14)-SUM(15)              !xi -2  3
            NQGW(16)=-SUM(12)-SUM(13)              !xi -1  3
            NQGW(17)=SUM(9)-SUM(6)+4.0d0*SUM(3)    !xi  0  3
            NQGW(18)=SUM(12)+SUM(13)               !xi  1  3
            NQGW(19)=SUM(14)+SUM(15)               !xi  2  3

          ELSE
            !xi -2 -3
            DO i=1,NXQ(-2,0,NXQ(-3,1,nq))
              n=n+1
              NQGW(n)=(SUM(14)+SUM(15))/
     '          DBLE(NXQ(-2,0,NXQ(-3,1,nq)))
            ENDDO

            !xi -1 -3
            DO i=1,NXQ(-1,0,NXQ(-3,1,nq))
              n=n+1
              NQGW(n)=(SUM(12)+SUM(13))/
     '          DBLE(NXQ(-1,0,NXQ(-3,1,nq)))
            ENDDO

            !xi  0 -3
            DO i=1,NXQ(-3,0,nq)
              n=n+1
              NQGW(n)=(-SUM(9)+SUM(6)+4.0d0*SUM(3))/
     '          DBLE(NXQ(-3,0,nq))
            ENDDO

            !xi  1 -3
            DO i=1,NXQ(1,0,NXQ(-3,1,nq))
              n=n+1
              NQGW(n)=(-SUM(12)-SUM(13))/
     '          DBLE(NXQ(1,0,NXQ(-3,1,nq)))
            ENDDO

            !xi  2 -3
            DO i=1,NXQ(2,0,NXQ(-3,1,nq))
              n=n+1
              NQGW(n)=(-SUM(14)-SUM(15))/
     '          DBLE(NXQ(2,0,NXQ(-3,1,nq)))
            ENDDO

            !xi -1 -2
            DO i=1,NXQ(-1,0,NXQ(-2,1,nq))
              n=n+1
              NQGW(n)=(SUM(10)+SUM(11))/
     '          DBLE(NXQ(-1,0,NXQ(-2,1,nq)))
            ENDDO

            !xi  0 -2
            DO i=1,NXQ(-2,0,nq)
              n=n+1
              NQGW(n)=(-SUM(8)+SUM(5)+4.0d0*SUM(2))/
     '          DBLE(NXQ(-2,0,nq))
            ENDDO

            !xi  1 -2
            DO i=1,NXQ(1,0,NXQ(-2,1,nq))
              n=n+1
              NQGW(n)=(-SUM(10)-SUM(11))/
     '          DBLE(NXQ(1,0,NXQ(-2,1,nq)))
            ENDDO

            !xi  0 -1
            DO i=1,NXQ(-1,0,nq)
              n=n+1
              NQGW(n)=(-SUM(7)+SUM(4)+4.0d0*SUM(1))/
     '          DBLE(NXQ(-1,0,nq))
            ENDDO

            !xi  0  0
            n=n+1
            NQGW(n)=-8.0d0*(SUM(1)+SUM(2)+SUM(3))

            !xi  0  1
            DO i=1,NXQ(1,0,nq)
              n=n+1
              NQGW(n)=(SUM(7)-SUM(4)+4.0d0*SUM(1))/
     '          DBLE(NXQ(1,0,nq))
            ENDDO

            !xi -1  2
            DO i=1,NXQ(-1,0,NXQ(2,1,nq))
              n=n+1
              NQGW(n)=(-SUM(10)-SUM(11))/
     '          DBLE(NXQ(-1,0,NXQ(2,1,nq)))
            ENDDO

            !xi  0  2
            DO i=1,NXQ(2,0,nq)
              n=n+1
              NQGW(n)=(SUM(8)-SUM(5)+4.0d0*SUM(2))/
     '          DBLE(NXQ(2,0,nq))
            ENDDO

            !xi  1  2
            DO i=1,NXQ(1,0,NXQ(2,1,nq))
              n=n+1
              NQGW(n)=(SUM(10)+SUM(11))/
     '          DBLE(NXQ(1,0,NXQ(2,1,nq)))
            ENDDO

            !xi -2  3
            DO i=1,NXQ(-2,0,NXQ(3,1,nq))
              n=n+1
              NQGW(n)=(-SUM(14)-SUM(15))/
     '          DBLE(NXQ(-2,0,NXQ(3,1,nq)))
            ENDDO

            !xi -1  3
            DO i=1,NXQ(-1,0,NXQ(3,1,nq))
              n=n+1
              NQGW(n)=(-SUM(12)-SUM(13))/
     '          DBLE(NXQ(-1,0,NXQ(3,1,nq)))
            ENDDO

            !xi  0  3
            DO i=1,NXQ(3,0,nq)
              n=n+1
              NQGW(n)=(SUM(9)-SUM(6)+4.0d0*SUM(3))/
     '          DBLE(NXQ(3,0,nq))
            ENDDO

            !xi  1  3
            DO i=1,NXQ(1,0,NXQ(3,1,nq))
              n=n+1
              NQGW(n)=(SUM(12)+SUM(13))/
     '          DBLE(NXQ(1,0,NXQ(3,1,nq)))
            ENDDO

            !xi  2  3
            DO i=1,NXQ(2,0,NXQ(3,1,nq))
              n=n+1
              NQGW(n)=(SUM(14)+SUM(15))/
     '          DBLE(NXQ(2,0,NXQ(3,1,nq)))
            ENDDO
          ENDIF

          IF(BIDOMAIN) THEN
            DO ni=1,15
              SUM(ni)=0.0d0
            ENDDO !ni
            DO nii=1,NITB
!           (C(1j)int + C(1j)ext) * G1j (upper)
              SUM(1)=SUM(1)+GUQ(nii,1)*(PROPQ(1,nii,1,1)+
     '          PROPQ(1,nii,1,2))
!           (C(2j)int + C(2j)ext) * G2j (upper)
              SUM(2)=SUM(2)+GUQ(nii,2)*(PROPQ(2,nii,1,1)+
     '          PROPQ(2,nii,1,2))
!           (C(3j)int + C(3j)ext) * G3j (upper)
              SUM(3)=SUM(3)+GUQ(nii,3)*(PROPQ(3,nii,1,1)+
     '          PROPQ(3,nii,1,2))
!           (C(1l)int + C(1l)ext) * [Christoffel*Gij (upper)](l)
              SUM(4)=SUM(4)+GCHQ(nii)*(PROPQ(1,nii,1,1)+
     '          PROPQ(1,nii,1,2))
!           (C(2l)int + C(2l)ext) * [Christoffel*Gij (upper)](l)
              SUM(5)=SUM(5)+GCHQ(nii)*(PROPQ(2,nii,1,1)+
     '          PROPQ(2,nii,1,2))
!           (C(3l)int + C(3l)ext) * [Christoffel*Gij (upper)](l)
              SUM(6)=SUM(6)+GCHQ(nii)*(PROPQ(3,nii,1,1)+
     '          PROPQ(3,nii,1,2))
              DO nij=1,NITB
!             (C(1k),j int + C(1k),j ext) * Gij (upper)
                SUM(7)=SUM(7)+GUQ(nii,nij)*(PROPQ(1,nii,1+nij,1)+
     '            PROPQ(1,nii,1+nij,2))
!             (C(2k),j int + C(2k),j ext) * Gij (upper)
                SUM(8)=SUM(8)+GUQ(nii,nij)*(PROPQ(2,nii,1+nij,1)+
     '            PROPQ(2,nii,1+nij,2))
!             (C(3k),j int + C(3k),j ext) * Gij (upper)
                SUM(9)=SUM(9)+GUQ(nii,nij)*(PROPQ(3,nii,1+nij,1)+
     '            PROPQ(3,nii,1+nij,2))
              ENDDO !nij
!           (C(1j)int + C(1j)ext) * G2j (upper)
              SUM(10)=SUM(10)+GUQ(nii,1)*(PROPQ(2,nii,1,1)+
     '          PROPQ(2,nii,1,2))
!           (C(2j)int + C(2j)ext) * G1j (upper)
              SUM(11)=SUM(11)+GUQ(nii,2)*(PROPQ(1,nii,1,1)+
     '          PROPQ(1,nii,1,2))
!           (C(1j)int + C(1j)ext) * G3j (upper)
              SUM(12)=SUM(12)+GUQ(nii,1)*(PROPQ(3,nii,1,1)+
     '          PROPQ(3,nii,1,2))
!           (C(3j)int + C(3j)ext) * G1j (upper)
              SUM(13)=SUM(13)+GUQ(nii,3)*(PROPQ(1,nii,1,1)+
     '          PROPQ(1,nii,1,2))
!           (C(2j)int + C(2j)ext) * G3j (upper)
              SUM(14)=SUM(14)+GUQ(nii,2)*(PROPQ(3,nii,1,1)+
     '          PROPQ(3,nii,1,2))
!           (C(3j)int + C(3j)ext) * G2j (upper)
              SUM(15)=SUM(15)+GUQ(nii,3)*(PROPQ(2,nii,1,1)+
     '          PROPQ(2,nii,1,2))
            ENDDO !nii

            IF((NXQ(-3,0,nq).LE.1).AND.(NXQ(-2,0,nq).LE.1).AND.
     '        (NXQ(-1,0,nq).LE.1).AND.(NXQ(1,0,nq).LE.1).AND.
     '        (NXQ(2,0,nq).LE.1).AND.(NXQ(3,0,nq).LE.1)) THEN

              !Usual case, no element splitting
              COEFFSEXT(1)=SUM(14)+SUM(15)                !xi -2 -3
              COEFFSEXT(2)=SUM(12)+SUM(13)                !xi -1 -3
              COEFFSEXT(3)=-SUM(9)+SUM(6)+4.0d0*SUM(3)    !xi  0 -3
              COEFFSEXT(4)=-SUM(12)-SUM(13)               !xi  1 -3
              COEFFSEXT(5)=-SUM(14)-SUM(15)               !xi  2 -3
              COEFFSEXT(6)=SUM(10)+SUM(11)                !xi -1 -2
              COEFFSEXT(7)=-SUM(8)+SUM(5)+4.0d0*SUM(2)    !xi  0 -2
              COEFFSEXT(8)=-SUM(10)-SUM(11)               !xi  1 -2
              COEFFSEXT(9)=-SUM(7)+SUM(4)+4.0d0*SUM(1)    !xi  0 -1
              COEFFSEXT(10)=-8.0d0*(SUM(1)+SUM(2)+SUM(3)) !xi  0  0
              COEFFSEXT(11)=SUM(7)-SUM(4)+4.0d0*SUM(1)    !xi  0  1
              COEFFSEXT(12)=-SUM(10)-SUM(11)              !xi -1  2
              COEFFSEXT(13)=SUM(8)-SUM(5)+4.0d0*SUM(2)    !xi  0  2
              COEFFSEXT(14)=SUM(10)+SUM(11)               !xi  1  2
              COEFFSEXT(15)=-SUM(14)-SUM(15)              !xi -2  3
              COEFFSEXT(16)=-SUM(12)-SUM(13)              !xi -1  3
              COEFFSEXT(17)=SUM(9)-SUM(6)+4.0d0*SUM(3)    !xi  0  3
              COEFFSEXT(18)=SUM(12)+SUM(13)               !xi  1  3
              COEFFSEXT(19)=SUM(14)+SUM(15)               !xi  2  3

            ELSE
              !xi -2 -3
              DO i=1,NXQ(-2,0,NXQ(-3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(SUM(14)+SUM(15))/
     '            DBLE(NXQ(-2,0,NXQ(-3,1,nq)))
              ENDDO

              !xi -1 -3
              DO i=1,NXQ(-1,0,NXQ(-3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(SUM(12)+SUM(13))/
     '            DBLE(NXQ(-1,0,NXQ(-3,1,nq)))
              ENDDO

              !xi  0 -3
              DO i=1,NXQ(-3,0,nq)
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(9)+SUM(6)+4.0d0*SUM(3))/
     '            DBLE(NXQ(-3,0,nq))
              ENDDO

              !xi  1 -3
              DO i=1,NXQ(1,0,NXQ(-3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(12)-SUM(13))/
     '            DBLE(NXQ(1,0,NXQ(-3,1,nq)))
              ENDDO

              !xi  2 -3
              DO i=1,NXQ(2,0,NXQ(-3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(14)-SUM(15))/
     '            DBLE(NXQ(2,0,NXQ(-3,1,nq)))
              ENDDO

              !xi -1 -2
              DO i=1,NXQ(-1,0,NXQ(-2,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(SUM(10)+SUM(11))/
     '            DBLE(NXQ(-1,0,NXQ(-2,1,nq)))
              ENDDO

              !xi  0 -2
              DO i=1,NXQ(-2,0,nq)
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(8)+SUM(5)+4.0d0*SUM(2))/
     '            DBLE(NXQ(-2,0,nq))
              ENDDO

              !xi  1 -2
              DO i=1,NXQ(1,0,NXQ(-2,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(10)-SUM(11))/
     '            DBLE(NXQ(1,0,NXQ(-2,1,nq)))
              ENDDO

              !xi  0 -1
              DO i=1,NXQ(-1,0,nq)
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(7)+SUM(4)+4.0d0*SUM(1))/
     '            DBLE(NXQ(-1,0,nq))
              ENDDO

              !xi  0  0
              n2=n2+1
              COEFFSEXT(n2)=-8.0d0*(SUM(1)+SUM(2)+SUM(3))

              !xi  0  1
              DO i=1,NXQ(1,0,nq)
                n2=n2+1
                COEFFSEXT(n2)=(SUM(7)-SUM(4)+4.0d0*SUM(1))/
     '            DBLE(NXQ(1,0,nq))
              ENDDO

              !xi -1  2
              DO i=1,NXQ(-1,0,NXQ(2,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(10)-SUM(11))/
     '          DBLE(NXQ(-1,0,NXQ(2,1,nq)))
              ENDDO

              !xi  0  2
              DO i=1,NXQ(2,0,nq)
                n2=n2+1
                COEFFSEXT(n2)=(SUM(8)-SUM(5)+4.0d0*SUM(2))/
     '            DBLE(NXQ(2,0,nq))
              ENDDO

              !xi  1  2
              DO i=1,NXQ(1,0,NXQ(2,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(SUM(10)+SUM(11))/
     '            DBLE(NXQ(1,0,NXQ(2,1,nq)))
              ENDDO

              !xi -2  3
              DO i=1,NXQ(-2,0,NXQ(3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(14)-SUM(15))/
     '            DBLE(NXQ(-2,0,NXQ(3,1,nq)))
              ENDDO

              !xi -1  3
              DO i=1,NXQ(-1,0,NXQ(3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(-SUM(12)-SUM(13))/
     '            DBLE(NXQ(-1,0,NXQ(3,1,nq)))
              ENDDO

              !xi  0  3
              DO i=1,NXQ(3,0,nq)
                n2=n2+1
                COEFFSEXT(n2)=(SUM(9)-SUM(6)+4.0d0*SUM(3))/
     '            DBLE(NXQ(3,0,nq))
              ENDDO

              !xi  1  3
              DO i=1,NXQ(1,0,NXQ(3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(SUM(12)+SUM(13))/
     '            DBLE(NXQ(1,0,NXQ(3,1,nq)))
              ENDDO

              !xi  2  3
              DO i=1,NXQ(2,0,NXQ(3,1,nq))
                n2=n2+1
                COEFFSEXT(n2)=(SUM(14)+SUM(15))/
     '            DBLE(NXQ(2,0,NXQ(3,1,nq)))
              ENDDO
            ENDIF
          ENDIF
        ENDIF !NITB

      ELSE !external grid point
        IF(SOLVEEIGHTPROBLEM.AND.(NX_LIST(0).EQ.1)) THEN
CMHT 24-03-00 NQGP_PIVOT(22,NQM) removed from parameter list, not used
C          NQGP(0,nq)=1
C          NQGP(1,nq)=nq
C          NQGP_PIVOT(1,nq)=1
C          COEFFSEXT(1)=1.0d0
          IF(SOLVE8_FLUXBC) THEN
            DO n=1,NQGP(0,nq)
              COEFFSEXT(n)=0.0d0
            ENDDO !n
            IF(nq.EQ.1) THEN
              COEFFSEXT(1)=1.0d0
            ELSE
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,
     '          NQXI,NXQ,COEFFSEXT,CQ(1,nq),DNUDXQ,DXDXIQ,DXDXIQ2,
     '          ERROR,*9999)
            ENDIF
          ELSE
            DO n=1,NQGP(0,nq)
              COEFFSEXT(n)=0.0d0
            ENDDO !n
            COEFFSEXT(1)=1.0d0
          ENDIF
        ELSE
          IF(IMPLICIT) THEN
            IF(FIXQ(nq,1,nx_trans)) THEN
              DO n=1,NQGP(0,nq)
                NQGW(n)=0.0d0
              ENDDO !n
              NQGW(1)=1.0d0
            ELSEIF(FIXQ(nq,2,nx_trans)) THEN
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,
     '          NQXI,NXQ,NQGW,CQ(3,nq),DNUDXQ,DXDXIQ,DXDXIQ2,
     '          ERROR,*9999)
            ELSE
              WRITE(ERROR,'(''>>No boundary condition set at point'
     '          //' (trans.) '',I8)') nq
              GOTO 9999
            ENDIF
          ELSE !explicit - always no flux!
C            CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,
C     '        NQXI,NXQ,NQGW,CQ(3,nq),DNUDXQ,DXDXIQ,DXDXIQ2,ERROR,*9999)
          ENDIF
          IF(BIDOMAIN) THEN
            IF(FIXQ(nq,1,nx_ext)) THEN !potential boundary condition
              COEFFSEXT(1)=1.0d0
              COEFFSEXT(2)=0.0d0
              COEFFSEXT(3)=0.0d0
            ELSEIF(FIXQ(nq,2,nx_ext)) THEN !flux boundary condition
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,
     '          NQXI,NXQ,COEFFSEXT,CQ(6,nq),DNUDXQ,DXDXIQ,DXDXIQ2,
     '          ERROR,*9999)
            ELSEIF(FIXQ(nq,3,nx_ext)) THEN
              COEFFSEXT(1)=1.0d0
              COEFFSEXT(2)=0.0d0
              COEFFSEXT(3)=0.0d0
            ELSE
              WRITE(ERROR,'(''>>No boundary condition set at point'
     '          //' (ext.) '',I8)') nq
              GOTO 9999
            ENDIF
          ENDIF !bidomain
        ENDIF
      ENDIF !internal/external grid pt

      CALL EXITS('CALC_GRID_COEF')
      RETURN
 9999 CALL ERRORS('CALC_GRID_COEF',ERROR)
      CALL EXITS('CALC_GRID_COEF')
      RETURN 1
      END



